/*
 *  particles.h
 *  zpic
 *
 *  Created by Ricardo Fonseca on 11/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#ifndef __PARTICLES__
#define __PARTICLES__

#include "zpic.h"
#include "emf.h"
#include "current.h"

#define MAX_SPNAME_LEN 32

typedef struct
{
	int ix, iy; 				// cell coordinates, local
	t_part_data x, y;			// particle coordenates relative to the cell (where the particle is inside the cell) ???
	t_part_data ux, uy, uz;		// velocity
} t_part;

enum density_type { UNIFORM, STEP, SLAB };

typedef struct {

	float n;					// reference density (defaults to 1.0, multiplies density profile)

	enum density_type type;		// Density profile type
	float start, end;			// Position of the plasma start/end, in simulation units

} t_density;


typedef struct {

	char name[MAX_SPNAME_LEN];
	int id;

	// Particle data buffer
	t_part* part;
	int np;
	int np_max;

	// Number of rows on each line of the global simulation space
	int nrow;
	// Number of rows on each line on this proc
	int nrow_local;

	// mass over charge ratio
	t_part_data m_q;

	// total kinetic energy
	double energy;

	// charge of individual particle
	t_part_data q;

	// Number of particles per cell
	int ppc[NUM_DIMS];

	// Density profile to inject
	t_density density;

	// Initial momentum of particles
	t_part_data ufl[3];
	t_part_data uth[3];

	// Simulation box info
	int nx[NUM_DIMS];
	int nx_local[NUM_DIMS];

	t_part_data dx[NUM_DIMS];
	t_part_data box[NUM_DIMS];

	// Time step
	float dt;

	// Iteration number
	int iter;

	// Moving window
	bool moving_window;
	int n_move;

} t_species;

void spec_new(t_species* spec, char name[], const t_part_data m_q, const int ppc[],
	const t_part_data ufl[], const t_part_data uth[],
	const int nx[], t_part_data box[], const float dt, t_density* density);

void spec_delete(t_species* spec);

void send_spec(t_species* spec, const int num_spec, int num_part_to_send[][NUM_ADJ], int fake_part_index[][NUM_ADJ]);

void spec_advance(t_species* spec, t_emf* emf, t_current* current, int part_seg_write_index[NUM_ADJ], int num_part_to_send[][NUM_ADJ]);

void wait_save_particles(t_species* species_array, const int n_spec);

void inject_particles(t_species* spec_array, const int num_spec);

void add_fake_particles(int fake_part_index[][NUM_ADJ], int part_seg_write_index[NUM_ADJ], int num_part_to_send[][NUM_ADJ],
						const bool moving_window, const int spec_id);

void spec_sort(t_species* spec);

double spec_time(void);
double spec_perf(void);

/*********************************************************************************************

 Diagnostics

 *********************************************************************************************/

#define CHARGE 		0x1000
#define PHA    		0x2000
#define PARTICLES   0x3000
#define X1     		0x0001
#define X2     		0x0002
#define U1     		0x0004
#define U2     		0x0005
#define U3     		0x0006

#define PHASESPACE(a,b) ((a) + (b)*16 + PHA)

void spec_deposit_pha(const t_species* spec, const int rep_type,
	const int pha_nx[], const float pha_range[][2], float* buf);

void spec_report(const t_species* spec, const int rep_type,
	const int pha_nx[], const float pha_range[][2]);

void spec_deposit_charge(const t_species* spec, float* charge);


#endif
