/**
 * ZPIC - em2d
 *
 * Weibel instability

Time for spec. advance = 419.327684 s
Time for emf   advance = 1.399192 s
Total simulation time  = 420.847216 s

Particle advance [nsec/part] = 99.775959 
Particle advance [Mpart/sec] = 10.022454

 */

#include <stdlib.h>
#include "../simulation.h"

void sim_init( t_simulation* sim ){

	// Time step
	float dt = 0.07;
	float tmax = 35.0;

	// Simulation box
	int   nx[2]  = { 512, 512 };
	float box[2] = { 51.2, 51.2 };

	// Diagnostic frequency
	int ndump = 10;

    // Initialize particles
	const int n_species = 2;
	t_species* species = (t_species *) malloc( n_species * sizeof( t_species ));

	// Use 4x4 particles per cell
	int ppc[] = {4,4};

	// Initial fluid and thermal velocities
	t_part_data ufl[] = { 0.0, 0.0, 0.6 };
	t_part_data uth[] = { 0.1, 0.1, 0.1 };

	spec_new( &species[0], "electrons", -1.0, ppc, ufl, uth, nx, box, dt, NULL );

	ufl[2] = -ufl[2];
	spec_new( &species[1], "positrons", +1.0, ppc, ufl, uth, nx, box, dt, NULL );

	// Initialize Simulation data
	sim_new( sim, nx, box, dt, tmax, ndump, species, n_species);

}

void sim_report( t_simulation* sim )
{
	// Jx, Jy, Jz
	current_report( &sim->current, 0 );
	current_report( &sim->current, 1 );
	current_report( &sim->current, 2 );

	// Bx, By, Bz
	emf_report( &sim->emf, BFLD, 0 );
	emf_report( &sim->emf, BFLD, 1 );
	emf_report( &sim->emf, BFLD, 2 );

	// Ex, Ey, Ez
	emf_report( &sim->emf, EFLD, 0 );
	emf_report( &sim->emf, EFLD, 1 );
	emf_report( &sim->emf, EFLD, 2 );
}
