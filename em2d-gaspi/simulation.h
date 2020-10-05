#ifndef __SIMULATION__
#define __SIMULATION__

#include <stdint.h>
#include "particles.h"
#include "emf.h"
#include "current.h"

#define STATIC_WINDOW 0
#define MOVING_WINDOW 1

typedef struct {

	// Time step
	float dt;
	float tmax;

	// Diagnostic frequency
	int ndump;

	// Simulation data
	int n_species;
	t_species* species;
	t_emf emf;
	t_current current;

	int nx_local[NUM_DIMS];

} t_simulation;


void sim_init( t_simulation* sim );
void sim_report( t_simulation* sim );

void sim_iter( t_simulation* sim );
void sim_report_energy( t_simulation* sim );

void sim_new(t_simulation* sim, int nx[NUM_DIMS], float box[NUM_DIMS], float dt, float tmax, int ndump, t_species* species, int n_species, const char moving_window);
int report( int n, int ndump );
void gaspi_report(t_simulation* sim);
void sim_timings( t_simulation* sim, uint64_t t0, uint64_t t1 );
void sim_add_laser( t_simulation* sim,  t_emf_laser* laser );
void sim_delete( t_simulation* sim );
void create_reporting_segments(int nx[NUM_DIMS]);

void sim_set_smooth( t_simulation* sim,  t_smooth* smooth );

#endif
