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

#include "zpic.h"
#include "simulation.h"
#include "emf.h"
#include "current.h"
#include "particles.h"
#include "timer.h"

// Include Simulation parameters here
// #include "input/lwfa-2000x512-8-1450.c"
// #include "input/weibel-9x9-1-500.c"
#include "input/weibel-128x128-256-500.c"


int main()
{
	// Initialize simulation
	t_simulation sim;
	sim_init( &sim );

    // Run simulation
	int n;
	float t;


	fprintf(stdout, "Starting simulation ...\n\n");

	uint64_t t0,t1;
	t0 = timer_ticks();

	for (n = 0, t = 0.0; t <= sim.tmax; n++, t = n * sim.dt)
	{
		fprintf(stdout,"n = %i, t = %f\n", n, t);

		if ( report ( n , sim.ndump ) )
		{
			sim_report( &sim );
		}

        sim_iter( &sim );
	}

	t1 = timer_ticks();
	fprintf(stdout, "\nSimulation ended.\n\n");

	// sim_report( &sim );

	// Simulation times
    sim_timings( &sim, t0, t1 );

	printf("Finished with ");
	for (int i = 0; i < sim.n_species; i++)
	{
		printf("%7d ", sim.species[i].np);
	}
	printf("particles\n"); fflush(stdout);

    // Cleanup data
    sim_delete( &sim );

	return 0;
}
