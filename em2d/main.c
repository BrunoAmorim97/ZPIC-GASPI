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
// #include "input/weibel.c"
// #include "input/weibel-small.c"
// #include "input/weibel-test.c"
// #include "input/weibel-test-large.c"
// #include "input/larger_weibel.c"
// #include "input/lwfa-test.c"
#include "input/lwfa.c"

int main()
{

	// Initialize simulation
	t_simulation sim;
	sim_init( &sim );

    // Run simulation
	int n;
	float t;

	// printf("%d %d\n", sim.species[0].np, sim.species[1].np); fflush(stdout);

	fprintf(stdout, "Starting simulation ...\n\n");

	uint64_t t0,t1;
	t0 = timer_ticks();

	for (n = 0, t = 0.0; t <= sim.tmax; n++, t = n * sim.dt)
	{
		fprintf(stdout,"n = %i, t = %f\n", n, t);

		if ( report ( n , sim.ndump ) )
		{
			// print_current(&sim.current);
			// print_emf_b(&sim.emf);
			// print_emf_e(&sim.emf);

			sim_report( &sim );
		}

        sim_iter( &sim );
	}

	t1 = timer_ticks();
	fprintf(stdout, "\nSimulation ended.\n\n");

	// sim_report( &sim );

	// Simulation times
    sim_timings( &sim, t0, t1 );

    // Cleanup data
    sim_delete( &sim );

	return 0;
}
