#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "simulation.h"
#include "timer.h"
#include "current.h"
#include "emf.h"

extern int proc_coords[NUM_DIMS];
extern int proc_block_low[NUM_DIMS];
extern int proc_block_high[NUM_DIMS];
extern int dims[NUM_DIMS];

extern const int gc[NUM_DIMS][NUM_DIMS];

extern int neighbour_nx[NUM_ADJ][NUM_DIMS];
extern gaspi_rank_t neighbour_rank[NUM_ADJ];
extern gaspi_rank_t proc_rank;
extern gaspi_rank_t num_procs;

extern t_part *particle_segments[NUM_ADJ];

extern int part_send_seg_size[NUM_ADJ];

extern t_vfld *emf_e_report_array;
extern t_vfld *emf_b_report_array;
extern t_vfld *current_report_array;

extern t_vfld *reporting_data;

int report(int n, int ndump)
{
	if (ndump > 0)
	{
		return !(n % ndump);
	}
	else
	{
		return 0;
	}
}

void sim_timings(t_simulation *sim, uint64_t t0, uint64_t t1)
{

	int npart = 0;
	int i;

	for (i = 0; i < sim->n_species; i++)
	{
		npart += sim->species[i].np;
	}

	// fprintf(stdout, "Time for spec. advance = %f s\n", spec_time());
	// fprintf(stdout, "Time for emf   advance = %f s\n", emf_time());
	fprintf(stdout, "Total simulation time on %d procs = %f s\n", num_procs, timer_interval_seconds(t0, t1));
	fprintf(stdout, "\n");

	/* 	if (spec_time() > 0)
	{
		double perf = spec_perf();
		fprintf(stdout, "Particle advance [nsec/part] = %f \n", 1.e9*perf);
		fprintf(stdout, "Particle advance [Mpart/sec] = %f \n", 1.e-6/perf);
	} */
}

// create the segments that will be used to transfer particles to and from other procs
void create_particle_segments(const int nx_local[NUM_DIMS], t_simulation *sim)
{
	// maximum number of particles that will move between procs on a single iteration, per cell, per species. Used to compute the segment size
	unsigned int max_particles_cell = 0;

	for (int i = 0; i < sim->n_species; i++)
	{
		//allocate enough segment size so that, on a single iteration, ppc * MAX_PPC_MULTIPLIER particles can move between procs, per border cell, per species
		max_particles_cell += sim->species[i].ppc[0] * sim->species[i].ppc[1] * MAX_PPC_MULTIPLIER;
	}

	const gaspi_size_t sizes[NUM_ADJ] =
		{
			//LEFT								CENTER								RIGHT
			max_particles_cell * 1,				max_particles_cell * nx_local[0],	max_particles_cell * 1,				// DOWN !!!
			max_particles_cell * nx_local[1], 										max_particles_cell * nx_local[1],	// CENTER
			max_particles_cell * 1,				max_particles_cell * nx_local[0],	max_particles_cell * 1				// UP !!!
		};

	gaspi_size_t size;
	gaspi_pointer_t array;

	// Create the particle segments, segments are created with size * 2 because they are both send and receive segments,
	// * 2 so that we have different zones depending on the iteration number, to prevent overwriting
	// Segment will be used as follows:
	// Segment => 	[Particle even Send zone, Particle even Receive zone, Particle odd Send zone, Particle odd Receive zone]
	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		size = sizes[dir];
		part_send_seg_size[dir] = size;

		SUCCESS_OR_DIE(gaspi_segment_alloc(dir, size * 2 * 2 * sizeof(t_part), GASPI_MEM_UNINITIALIZED));
		SUCCESS_OR_DIE(gaspi_segment_register(dir, neighbour_rank[dir], GASPI_BLOCK));

		SUCCESS_OR_DIE(gaspi_segment_ptr(dir, &array));
		particle_segments[dir] = (t_part *)array;
	}
}

void sim_iter(t_simulation *sim)
{
	// Number of particles this process will send to each direction, per species
	int num_part_to_send[sim->n_species][NUM_ADJ]; memset(num_part_to_send, 0, sim->n_species * NUM_ADJ * sizeof(int));

	// Segment offset multiplier, 0 if iteration is even, 2 if odd
	const int seg_offset_mult = (sim->species->iter % 2) == 0 ? 0 : 2;

	// First available position to write a particle to, for each particle segment. Depends on iteration number
	int part_seg_write_index[NUM_ADJ];
	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		part_seg_write_index[dir] = seg_offset_mult * part_send_seg_size[dir];
	}

	// Index of the fake particle, for each species for each particle segment
	int fake_part_index[sim->n_species][NUM_ADJ];



	t_current* const restrict current = &sim->current;
	t_emf* const restrict emf = &sim->emf;

	current_zero(current);

	// Advance species
	for (int spec_i = 0; spec_i < sim->n_species; spec_i++)
	{
		spec_advance(&sim->species[spec_i], emf, current);
	}

	send_current(current);

	// Advance EM field using Yee algorithm 
	yee_b(emf);

	// Wait and update current data
	wait_save_update_current(current);
	current_smooth(current);
	current->iter++;

	// Advance EM fields using Yee algorithm modified for having E and B time centered
	// yee_b(emf); // already done
	yee_e(emf, current);
	yee_b(emf);

	// true if window will be moved this iteration, false otherwise
	const bool moving_window_iter = emf->moving_window && ( ((emf->iter + 1) * emf->dt) > (emf->dx[0] * (emf->n_move + 1)) );
	send_emf_gc(emf, moving_window_iter);

	// Move emf window if needed
	if(moving_window_iter) emf_move_window(emf);
	
	for (int spec_i = 0; spec_i < sim->n_species; spec_i++)
	{
		// Check for particles leaving this proc, copy them to particle segments if needed
		check_leaving_particles(&sim->species[spec_i], num_part_to_send, fake_part_index, part_seg_write_index);

		// Send particles on the particle segments
		send_spec(&sim->species[spec_i], sim->n_species, num_part_to_send, fake_part_index);

		// Inject new particles, if needed
		inject_particles(&sim->species[spec_i]);
	}

	wait_save_emf_gc(emf, moving_window_iter);
	emf->iter++;

	// wait for particle writes and copy them from segments to species array
	wait_save_particles(sim->species, sim->n_species);
}

void sim_set_spec_moving_window(t_simulation *sim)
{
	for (int i = 0; i < sim->n_species; i++)
		sim->species[i].moving_window = 1;
}

void sim_new(t_simulation *sim, int nx[NUM_DIMS], float box[NUM_DIMS], float dt, float tmax, int ndump, t_species *species, int n_species,
			 const bool moving_window)
{
	// Probably not necessary, just to be sure
	assign_proc_blocks(nx);

	discover_neighbours(proc_coords, dims, nx);

	sim->n_species = n_species;
	sim->species = species;

	if (moving_window)
		sim_set_spec_moving_window(sim);

	sim->dt = dt;
	sim->tmax = tmax;
	sim->ndump = ndump;

	const int nx_local[NUM_DIMS] = {BLOCK_SIZE(proc_coords[0], dims[0], nx[0]), BLOCK_SIZE(proc_coords[1], dims[1], nx[1])};

	// printf("size at proc %2d is %d %d\n", proc_rank, nx_local[0], nx_local[1]);

	// If this check fails change number of nodes or increase simulation nx size, number of procs should not be a prime number
	// Proc simulation regions need to have at least 2 cells on each dimention
	// On moving window simulations, procs need at least 3 cells on the x axis
	if (!moving_window)
		assert(nx_local[0] >= 2 && nx_local[1] >= 2);
	else
		assert(nx_local[0] >= 3 && nx_local[1] >= 2);

	// If this check fails increase number of procs so that no proc has the same proc as a neighbour on 2 different directions.
	// Or increase simulation nx so that those 2 procs have at least 3 cells on the axis they meet on.
	// The method used to update guard cells current will lead to desynchronization between edge guard cells
	// that represent the same cell on specific configurations with few procs and small numbers of cells per proc.
	// To prevent those configurations:
	assert(!(neighbour_rank[DOWN] == neighbour_rank[UP] && neighbour_nx[DOWN][1] <= gc[1][1] && num_procs > 1));
	assert(!(neighbour_rank[LEFT] == neighbour_rank[RIGHT] && neighbour_nx[LEFT][0] <= gc[0][1] && num_procs > 1));

	create_particle_segments(nx_local, sim);

	emf_new(&sim->emf, nx, nx_local, box, dt, moving_window);
	current_new(&sim->current, nx, nx_local, box, dt, moving_window);

	// Check time step
	float cour = sqrtf(1.0f / (1.0f / (sim->emf.dx[0] * sim->emf.dx[0]) + 1.0f / (sim->emf.dx[1] * sim->emf.dx[1])));
	if (dt >= cour)
	{
		fprintf(stdout, "Invalid timestep, courant condition violation, dtmax = %f \n", cour);
		exit(-1);
	}
}

void sim_add_laser(t_simulation *sim, t_emf_laser *laser)
{
	emf_add_laser(&sim->emf, laser);
}

void sim_set_smooth(t_simulation *sim, t_smooth *smooth)
{
	curr_set_smooth(&sim->current, smooth);
}

void sim_report_energy(t_simulation *sim)
{
	int i;

	double emf_energy[6];
	double part_energy[sim->n_species];

	emf_get_energy(&sim->emf, emf_energy);
	double tot_emf = emf_energy[0];
	for (i = 0; i < 6; i++)
	{
		tot_emf += emf_energy[i];
	}

	double tot_part = 0;
	for (i = 0; i < sim->n_species; i++)
	{
		part_energy[i] = sim->species[i].energy;
		tot_part += part_energy[i];
	}

	printf("Energy (fields | particles | total) = %e %e %e\n",
		   tot_emf, tot_part, tot_emf + tot_part);
}

void create_root_reporting_segments(int nx[NUM_DIMS])
{
	gaspi_size_t size = (gc[0][0] + nx[0] + gc[0][1]) * (gc[1][0] + nx[1] + gc[1][1]) * sizeof(t_vfld);
	gaspi_pointer_t pointer;

	SUCCESS_OR_DIE(gaspi_segment_alloc(CURRENT_REPORT, size, GASPI_MEM_INITIALIZED));
	SUCCESS_OR_DIE(gaspi_segment_alloc(EMF_B_REPORT, size, GASPI_MEM_INITIALIZED));
	SUCCESS_OR_DIE(gaspi_segment_alloc(EMF_E_REPORT, size, GASPI_MEM_INITIALIZED));

	// Resgister segments on all procs
	for (int proc_i = 0; proc_i < num_procs; proc_i++)
	{
		SUCCESS_OR_DIE(gaspi_segment_register(CURRENT_REPORT, proc_i, GASPI_BLOCK));
		SUCCESS_OR_DIE(gaspi_segment_register(EMF_B_REPORT, proc_i, GASPI_BLOCK));
		SUCCESS_OR_DIE(gaspi_segment_register(EMF_E_REPORT, proc_i, GASPI_BLOCK));
	}

	SUCCESS_OR_DIE(gaspi_segment_ptr(CURRENT_REPORT, &pointer));
	current_report_array = (t_vfld *)pointer;

	SUCCESS_OR_DIE(gaspi_segment_ptr(EMF_B_REPORT, &pointer));
	emf_b_report_array = (t_vfld *)pointer;

	SUCCESS_OR_DIE(gaspi_segment_ptr(EMF_E_REPORT, &pointer));
	emf_e_report_array = (t_vfld *)pointer;
}

// Wait for report writes
void wait_report_writes()
{
	// For each proc
	for (int proc_i = 0; proc_i < num_procs; proc_i++)
	{
		gaspi_notification_id_t id;
		gaspi_notification_t value;

		// Wait for Current
		SUCCESS_OR_DIE(gaspi_notify_waitsome(
			CURRENT_REPORT,
			proc_i,
			1,
			&id,
			GASPI_BLOCK));
		SUCCESS_OR_DIE(gaspi_notify_reset(CURRENT_REPORT, id, &value));

		// Wait for EMF B
		SUCCESS_OR_DIE(gaspi_notify_waitsome(
			EMF_B_REPORT,
			proc_i,
			1,
			&id,
			GASPI_BLOCK));
		SUCCESS_OR_DIE(gaspi_notify_reset(EMF_B_REPORT, id, &value));

		// Wait for EMF E
		SUCCESS_OR_DIE(gaspi_notify_waitsome(
			EMF_E_REPORT,
			proc_i,
			1,
			&id,
			GASPI_BLOCK));
		SUCCESS_OR_DIE(gaspi_notify_reset(EMF_E_REPORT, id, &value));
	}
}

void gaspi_report(t_simulation *sim)
{
	#define NUM_REPORTING_DATA_TYPES 3

	// printf("GASPI REPORT\n"); fflush(stdout);
	static bool created_segments = false;

	t_current *current = &sim->current;
	const int size_x = current->nx_local[0];
	const int size_y = current->nx_local[1];
	const int nrow = current->nrow_local; // Local nrow

	const int global_nrow = gc[0][0] + current->nx[0] + gc[0][1];
	const int global_buff_offset = gc[0][0] + (gc[1][0] * global_nrow); //offset not in bytes

	if (!created_segments)
	{
		created_segments = true;
		if (proc_rank == ROOT)
		{
			create_root_reporting_segments(current->nx);
		}

		gaspi_size_t seg_size = size_x * size_y * NUM_REPORTING_DATA_TYPES * sizeof(t_vfld);
		gaspi_pointer_t array_pointer;

		// Reporting segment will be used as follows:
		// Reporting => [Curr, EMF B, EMF E]
		SUCCESS_OR_DIE(gaspi_segment_alloc(REPORTING, seg_size, GASPI_MEM_UNINITIALIZED));

		SUCCESS_OR_DIE(gaspi_segment_ptr(REPORTING, &array_pointer));
		reporting_data = (t_vfld *)array_pointer;

		SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
	}

	int copy_index = 0;
	const size_t size_x_bytes = size_x * sizeof(t_vfld);

	const t_vfld *const restrict data_pointers[NUM_REPORTING_DATA_TYPES] = {sim->current.J, sim->emf.B, sim->emf.E};
	const gaspi_segment_id_t remote_segments[NUM_REPORTING_DATA_TYPES] = {CURRENT_REPORT, EMF_B_REPORT, EMF_E_REPORT};

	for (int i = 0; i < NUM_REPORTING_DATA_TYPES; i++)
	{

		const t_vfld *const restrict data_array = data_pointers[i];
		const gaspi_segment_id_t remote_segment = remote_segments[i];

		// Copy data from simulation arrays to the segment
		for (int y = 0; y < size_y; y++)
		{
			// Copy each line
			memcpy(&reporting_data[copy_index], &data_array[y * nrow], size_x_bytes);
			copy_index += size_x;
		}

		gaspi_offset_t local_offset = i * size_x * size_y * sizeof(t_vfld);															  // in bytes
		gaspi_offset_t remote_offset = (global_buff_offset + proc_block_low[0] + (proc_block_low[1] * global_nrow)) * sizeof(t_vfld); // in bytes

		SUCCESS_OR_DIE(gaspi_wait(Q_REPORTING, GASPI_BLOCK));

		// Send data on segment
		for (int y = 0; y < size_y; y++)
		{
			// Send each line
			SUCCESS_OR_DIE(gaspi_write(
				REPORTING,		// The segment id where data is located.
				local_offset,	// The offset where the data is located.
				ROOT,			// The rank where to write and notify.
				remote_segment, // The remote segment id to write the data to.
				remote_offset,	// The remote offset where to write to.
				size_x_bytes,	// The size of the data to write.
				Q_REPORTING,	// The queue where to post the request.
				GASPI_BLOCK		// Timeout in milliseconds.
				));

			local_offset += size_x_bytes;
			remote_offset += global_nrow * sizeof(t_vfld);
		}

		// Send notification
		SUCCESS_OR_DIE(gaspi_notify(
			remote_segment, // The remote segment id.
			ROOT,			// The rank to notify.
			proc_rank,		// The notification id.
			1,				// The notification value.
			Q_REPORTING,	// The queue to post the notification request.
			GASPI_BLOCK		// Timeout in milliseconds
			));
	}

	if (proc_rank == ROOT)
	{
		// Save original array pointers
		t_vfld *B_real = sim->emf.B;
		t_vfld *B_buff_real = sim->emf.B_buff;

		t_vfld *E_real = sim->emf.E;
		t_vfld *E_buff_real = sim->emf.E_buff;

		t_vfld *J_real = sim->current.J;
		t_vfld *J_buff_real = sim->current.J_buff;

		// Change simulation struct emf and current array pointers to report segments
		sim->emf.B = emf_b_report_array + global_buff_offset;
		sim->emf.E = emf_e_report_array + global_buff_offset;
		sim->current.J = current_report_array + global_buff_offset;

		sim->emf.B_buff = emf_b_report_array;
		sim->emf.E_buff = emf_e_report_array;
		sim->current.J_buff = current_report_array;

		// Wait for writes
		wait_report_writes();

		// Call the original ZPIC report function
		sim_report(sim);

		// Switch simulation struct emf and current array pointers back
		sim->emf.B = B_real;
		sim->emf.B_buff = B_buff_real;
		sim->emf.E = E_real;
		sim->emf.E_buff = E_buff_real;
		sim->current.J = J_real;
		sim->current.J_buff = J_buff_real;
	}

	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
}

/* 		t_current* current = &sim->current;
		t_vfld* J = sim->current.J;
		for (int y = -gc[1][0]; y < current->nx[1] + gc[1][1]; y++)
		{
			if (y == 0 || y == current->nx[1])
			{
				printf("\n");
			}
			
			for (int x = -gc[0][0]; x < current->nx[0] + gc[0][1]; x++)
			{
				if (x == 0 || x == current->nx[0])
				{
					printf("    ");
				}
				
				printf("%f ", J[y * current->nrow + x].z);
			}
			printf("\n");
		}
		printf("\n===EMF===\n");

		t_emf* emf = &sim->emf;
		t_vfld* E = sim->emf.E;
		for (int y = -gc[1][0]; y < emf->nx[1] + gc[1][1]; y++)
		{
			if (y == 0 || y == emf->nx[1])
			{
				printf("\n");
			}
			
			for (int x = -gc[0][0]; x < emf->nx[0] + gc[0][1]; x++)
			{
				if (x == 0 || x == emf->nx[0])
				{
					printf("    ");
				}
				
				printf("%f ", E[y * emf->nrow + x].z);
			}
			printf("\n");
		}
		printf("\n"); fflush(stdout); */

/* printf("CURRENT X\n");
		for (int y = -gc[1][0]; y < current->nx[1] + gc[1][1]; y++)
		{
			if (y == 0 || y == current->nx[1])
			{
				printf("\n");
			}

			for (int x = -gc[0][0]; x < current->nx[0] + gc[0][1]; x++)
			{
				if (x == 0 || x == current->nx[0])
				{
					printf("    ");
				}
				
				printf("%f ", current->J[y * current->nrow + x].x);
			}
			printf("\n");
		}
		printf("\n");

		printf("CURRENT Y\n");
		for (int y = -gc[1][0]; y < current->nx[1] + gc[1][1]; y++)
		{
			if (y == 0 || y == current->nx[1])
			{
				printf("\n");
			}

			for (int x = -gc[0][0]; x < current->nx[0] + gc[0][1]; x++)
			{
				if (x == 0 || x == current->nx[0])
				{
					printf("    ");
				}
				
				printf("%f ", current->J[y * current->nrow + x].y);
			}
			printf("\n");
		}
		printf("\n");

		printf("CURRENT Z\n");
		for (int y = -gc[1][0]; y < current->nx[1] + gc[1][1]; y++)
		{
			if (y == 0 || y == current->nx[1])
			{
				printf("\n");
			}

			for (int x = -gc[0][0]; x < current->nx[0] + gc[0][1]; x++)
			{
				if (x == 0 || x == current->nx[0])
				{
					printf("    ");
				}
				
				printf("%f ", current->J[y * current->nrow + x].z);
			}
			printf("\n");
		}
		printf("\n"); */