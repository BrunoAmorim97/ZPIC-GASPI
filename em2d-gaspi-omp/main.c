/*
Copyright (C) 2017 Instituto Superior Tecnico

This file is part of the ZPIC Educational code suite

The ZPIC Educational code suite is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The ZPIC Educational code suite is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with the ZPIC Educational code suite. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include <mpi.h>
#include <GASPI.h>
#include <omp.h>

#include "zpic.h"
#include "simulation.h"
#include "emf.h"
#include "current.h"
#include "particles.h"
#include "timer.h"

// Include Simulation parameters here
// #include "input/(simulation name)-(num cells x)*(num cells y)-(particles per cell)-(num iterations).c"

#include "input/weibel-128x128-256-500.c"
// #include "input/weibel-9x9-1-500.c"
// #include "input/weibel-512x512-64-500.c"
// #include "input/weibel-512x512-256-500.c"
// #include "input/lwfa-2000x256-8-1450.c"

gaspi_rank_t proc_rank;
gaspi_rank_t num_procs;

// Number of proc divisions on each axis of the global space
int dims[NUM_DIMS];

// Process coordinates
int proc_coords[NUM_DIMS];

// index 0 for simulation space left edge, 1 for right edge
bool is_on_edge[2];

// Number of guard cells for linear interpolation
const int gc[2][2] = {	{1,2},
						{1,2} };

// Process cell block low coords
int proc_block_low[NUM_DIMS];

// Process cell block high coords
int proc_block_high[NUM_DIMS];

// Process cell block sizes
int proc_block_size[NUM_DIMS];

// Rank on every direction
gaspi_rank_t neighbour_rank[NUM_ADJ];

// Simulation space size of neighbour procs
int neighbour_nx[NUM_ADJ][NUM_DIMS];

// Pointers to particle segments
t_part* particle_segments[NUM_ADJ];

// Size of the receive, and send, part of the particle segments (in number of particles)
unsigned int part_send_seg_size[NUM_ADJ];

// Pointers to current segments
t_vfld* current_segments[NUM_ADJ];

// Pointers to current kernel smoothing segments
t_vfld* current_kernel_smoothing_segments[NUM_ADJ];

// Current segment sizes
int curr_send_size[NUM_ADJ][NUM_DIMS];

// Top left coord of the curr cells to send
int curr_cell_to_send_starting_coord[NUM_ADJ][NUM_DIMS];

// Current kernel transmission sizes
int curr_kernel_size[NUM_ADJ][NUM_DIMS];

// Top left coord of the curr kernel cells to send
int curr_kernel_send_coord[NUM_ADJ][NUM_DIMS];

// Top left coord of the curr cells this proc will override
int curr_kernel_write_coord[NUM_ADJ][NUM_DIMS];

// EMF segment sizes (normal iteration), on the perspective of the sender
int emf_seg_size[NUM_ADJ][NUM_DIMS];

// Top left coord of the EMF cells to send (normal iteration), on the perspective of the sender
int emf_cell_to_send_starting_coord[NUM_ADJ][NUM_DIMS];

// Top left coord of the EMF cells this proc will receive (normal iteration), on the perspective of the receiver
int emf_cell_to_write_starting_coord[NUM_ADJ][NUM_DIMS];

// Pointers to EMF segments
t_vfld* emf_segments[NUM_ADJ];

// Segments used to receive EMF and current data for reporting
t_vfld* reporting_data;

t_vfld* emf_e_report_array;
t_vfld* emf_b_report_array;
t_vfld* current_report_array;

// Thread private struct to save thread current
t_current_reduce current_reduce_priv;
#pragma omp threadprivate(current_reduce_priv)

int main(int argc, char* argv[])
{
	if(argc >= 2)
		omp_set_num_threads(atoi(argv[1]));

	MPI_Init(&argc, &argv);
	SUCCESS_OR_DIE(gaspi_proc_init(GASPI_BLOCK));

	SUCCESS_OR_DIE(gaspi_proc_rank(&proc_rank));
	SUCCESS_OR_DIE(gaspi_proc_num(&num_procs));

	// printf("Hello from rank %d of %d\n", proc_rank, num_procs);

	create_dims(num_procs);
	cart_coords(proc_rank, proc_coords);

	// printf("I have proc coords x:%d y:%d\n", proc_coords[0], proc_coords[1]);

	if (proc_rank == ROOT)
	{
		printf("Dims:[%d,%d] with %d procs\n", dims[0], dims[1], num_procs);
	}

	// Initialize simulation
	t_simulation sim;
	sim_init(&sim);

	#pragma omp parallel
	alloc_private_current_reduce(&sim.current);

	// Run simulation
	int n;
	float t;
	uint64_t t0, t1;

	// printf("%d %d\n", sim.species[0].np, sim.species[1].np);

	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
	t0 = timer_ticks();

	if (proc_rank == ROOT)
	{
		printf("Running with max %d threads...\n", omp_get_max_threads());
		printf("Starting simulation ...\n\n"); fflush(stdout);
	}

	for (n = 0, t = 0.0; t <= sim.tmax; n++, t = n * sim.dt)
	{
		if (proc_rank == ROOT)
		{
			printf("n = %i, t = %f\n", n, t); fflush(stdout);
		}

		if (report(n, sim.ndump)) gaspi_report(&sim);

		sim_iter(&sim);

		// printf("proc %d has %7d %7d particles\n", proc_rank, sim.species[0].np, sim.species[1].np); fflush(stdout);
		// printf("proc %d has %7d particles\n", proc_rank, sim.species[0].np); fflush(stdout);
	}

	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));

	if (proc_rank == ROOT)
	{
		t1 = timer_ticks();

		printf("\nSimulation ended.\n\n"); fflush(stdout);

		// Simulation times
		sim_timings(&sim, t0, t1);
	}

	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
	
	printf("Proc %2d finished with ", proc_rank);
	for (int i = 0; i < sim.n_species; i++)
	{
		printf("%7d ", sim.species[i].np);
	}
	printf("particles\n"); fflush(stdout);

	SUCCESS_OR_DIE(gaspi_proc_term(GASPI_BLOCK));
	MPI_Finalize();

	return 0;
}