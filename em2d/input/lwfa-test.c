/**
 * ZPIC - em2d
 *
 * Laser Wakefield Acceleration
 */

#include <stdlib.h>
#include <math.h>

#include "../simulation.h"

void sim_init( t_simulation* sim ){

	// Time step
	float dt = 0.014;
	// float tmax = 20.314;
	float tmax = 0.7;

	// Simulation box
	int   nx[2]  = { 9, 9 };
	float box[2] = { 1.8, 1.8 };

	// Diagnostic frequency
	int ndump = 0;

    // Initialize particles
	const int n_species = 1;

	// Use 4x2 particles per cell
	int ppc[] = {1,1};

	// Density profile
	t_density density = { .type = STEP, .start = 0.1 };

	t_species* species = (t_species *) malloc( n_species * sizeof( t_species ));
	spec_new( &species[0], "electrons", -1.0, ppc, NULL, NULL, nx, box, dt, NULL);//&density );

	// Initialize Simulation data
	sim_new(sim, nx, box, dt, tmax, ndump, species, n_species);

	sim_set_moving_window(sim);
	
	// Add laser pulse (this must come after sim_new)
	t_emf_laser laser = {
		.type = GAUSSIAN,
		.start = 0.2,
		.fwhm  = 0.2,
		.a0 = 0.5,
		.omega0 = 1.0,
		.W0 = 0.2,
		.focus = 0.5,
		.axis = 0.4,
		.polarization = M_PI_2
    };
	sim_add_laser( sim, &laser );


	// Set current smoothing (this must come after sim_new)
	t_smooth smooth = {
		.xtype = COMPENSATED,
		.xlevel = 4
	};

	sim_set_smooth( sim, &smooth );
}


void sim_report( t_simulation* sim ){

	// Bx, By, Bz
	emf_report( &sim->emf, BFLD, 0 );
	emf_report( &sim->emf, BFLD, 1 );
	emf_report( &sim->emf, BFLD, 2 );

	/*
	// All electric field components
	emf_report( &sim->emf, EFLD, 0 );
	emf_report( &sim->emf, EFLD, 1 );
	emf_report( &sim->emf, EFLD, 2 );

	// Charge density
	spec_report( &sim->species[0], CHARGE, NULL, NULL );

    // x1u1 phasespace
	const int pha_nx[] = {1024,512};
	const float pha_range[][2] = {{0.0,20.0}, {-2.0,+2.0}};
	spec_report(&sim->species[0], PHASESPACE(X1,U1), pha_nx, pha_range);
	*/

}
