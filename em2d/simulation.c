#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "simulation.h"
#include "timer.h"


int report( int n, int ndump )
{
	if (ndump > 0)
	{
		return ! (n % ndump);
	}
	else
	{
		return 0;
	}
}

void sim_iter( t_simulation* sim )
{
	// Advance particles and deposit current
	current_zero( &sim -> current );
	for (int i = 0; i < sim -> n_species; i++)
	{
		// printf("Advancing particles from species %d\n", i); fflush(stdout);
		spec_advance(&sim -> species[i], &sim -> emf, &sim -> current );
	}

	// Update current boundary conditions and advance iteration, CHANGES CURRENT MATRIX
	current_update( &sim -> current );

	// // Print particles
	// for (int i = 0; i < sim -> n_species; i++)
	// {
	// 	printf("part from species %d: %d\n", i, sim->species[i].np);
	// 	for (int j = 0; j < sim->species[i].np; j++)
	// 	{
	// 		t_part part = sim->species[i].part[j];
	// 		printf("Part: ix:%d, iy:%d, x:%f, y:%f, ux:%f, uy:%f, uz:%f\n", part.ix, part.iy, part.x, part.y, part.ux, part.uy, part.uz); fflush(stdout);
	// 	}	
	// }
	
	// Advance EM fields
	emf_advance( &sim -> emf, &sim -> current );
}

void sim_timings( t_simulation* sim, uint64_t t0, uint64_t t1 )
{

	int npart = 0;
	int i;

	for(i=0; i<sim -> n_species; i++)
	{
		npart += sim -> species[i].np;
	}

	fprintf(stdout, "Time for spec. advance = %f s\n", spec_time());
	fprintf(stdout, "Time for emf   advance = %f s\n", emf_time());
	fprintf(stdout, "Total simulation time  = %f s\n", timer_interval_seconds(t0, t1));
	fprintf(stdout, "\n");

	if (spec_time()>0)
	{
		double perf = spec_perf();
		fprintf(stdout, "Particle advance [nsec/part] = %f \n", 1.e9*perf);
		fprintf(stdout, "Particle advance [Mpart/sec] = %f \n", 1.e-6/perf);
	}
}

void sim_new( t_simulation* sim, int nx[], float box[], float dt, float tmax, int ndump, t_species* species, int n_species )
{
	sim -> dt = dt;
	sim -> tmax = tmax;
	sim -> ndump = ndump;

	emf_new(&sim -> emf, nx, box, dt);
	current_new(&sim -> current, nx, box, dt);

	sim -> n_species = n_species;
	sim -> species = species;

	// Check time step
	float cour = sqrtf( 1.0f/( 1.0f/(sim->emf.dx[0]*sim->emf.dx[0]) + 1.0f/(sim->emf.dx[1]*sim->emf.dx[1]) ) );
	if ( dt >= cour )
	{
		fprintf(stdout, "Invalid timestep, courant condition violation, dtmax = %f \n", cour );
		exit(-1);
	}

}

void sim_add_laser( t_simulation* sim, t_emf_laser* laser )
{
	emf_add_laser( &sim->emf, laser );
}

void sim_set_smooth( t_simulation* sim, t_smooth* smooth )
{
    if ( (smooth -> xtype != NONE) && (smooth -> xlevel <= 0) )
	{
    	fprintf(stdout, "Invalid smooth level along x direction\n");
    	exit(-1);
    }

    if ( (smooth -> ytype != NONE) && (smooth -> ylevel <= 0) )
	{
    	fprintf(stdout, "Invalid smooth level along y direction\n");
    	exit(-1);
    }

	sim -> current.smooth = *smooth;
}

void sim_set_moving_window( t_simulation* sim )
{

	sim -> emf.moving_window = 1;
	sim -> current.moving_window = 1;

    int i;
	for(i=0; i<sim -> n_species; i++)
		sim -> species[i].moving_window = 1;

}


void sim_report_energy( t_simulation* sim )
{
	int i;

	double emf_energy[6];
	double part_energy[ sim -> n_species ];

	emf_get_energy( &sim -> emf, emf_energy );
	double tot_emf = emf_energy[0];
	for( i = 0; i < 6; i++ )
	{
		tot_emf += emf_energy[i];
	}

	double tot_part = 0;
	for( i = 0; i < sim -> n_species; i++ )
	{
		part_energy[i] = sim -> species[i].energy;
		tot_part += part_energy[i];
	}

	printf("Energy (fields | particles | total) = %e %e %e\n",
		tot_emf, tot_part, tot_emf+tot_part);

}

void sim_delete( t_simulation* sim )
{

	int i;
	for (i = 0; i<sim->n_species; i++) spec_delete( &sim->species[i] );

	free( sim->species );

	current_delete( &sim->current );
	emf_delete( &sim->emf );

}

