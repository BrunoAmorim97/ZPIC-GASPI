/*
 *  particles.c
 *  zpic
 *
 *  Created by Ricardo Fonseca on 11/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "particles.h"

#include "random.h"
#include "emf.h"
#include "current.h"

#include "zdf.h"
#include "timer.h"

static double _spec_time = 0.0;
static double _spec_npush = 0.0;

extern int proc_rank;

extern int proc_coords[NUM_DIMS];
extern char is_on_edge[2];
extern int proc_block_low[NUM_DIMS];
extern int proc_block_high[NUM_DIMS];
extern int dims[NUM_DIMS];
extern const int gc[NUM_DIMS][NUM_DIMS];

extern int part_send_seg_size[NUM_ADJ];

extern gaspi_rank_t neighbour_rank[NUM_ADJ];
extern int neighbour_nx[NUM_ADJ][NUM_DIMS];

extern t_part* particle_segments[NUM_ADJ];


/**
 * Returns the total time spent pushing particles (includes boundaries and moving window)
 * @return  Total time in seconds
 */
double spec_time( void )
{
	return _spec_time;
}

/**
 * Returns the performance achieved by the code (push time)
 * @return  Performance in seconds per particle
 */
double spec_perf( void )
{
	return (_spec_npush > 0 )? _spec_time / _spec_npush: 0.0;
}

/*********************************************************************************************
 Initialization
 *********************************************************************************************/

void spec_set_u( t_species* spec, const int start, const int range[NUM_DIMS][NUM_DIMS])
{
	const int npc = spec->ppc[0] * spec->ppc[1];
	int i = start;

	// Possible ranges:

	// const int range[NUM_DIMS][NUM_DIMS] =  {{0, nx[0]-1},
	// 										   {0, nx[1]-1}};

	// const int range[NUM_DIMS][NUM_DIMS] =  {{spec->nx[0]-1,spec->nx[0]-1},
	// 										   {            0,spec->nx[1]-1}};
	// ^^^ This one is only called on the procs that will generate particles ! ^^^

	// generate random numbers up to the first proc cell
	// starting from y = 0, up to the y of the first cell on this proc, excluding
	for (int y = range[1][0]; y < proc_block_low[1]; y++)
	{
		// for each row
		for (int x = range[0][0]; x <= range[0][1]; x++)
		{
			// for particles in this cell
			for (int part = 0; part < npc; part++)
			{
				rand_norm(); rand_norm(); rand_norm();
			}
		}
	}
	
	// for each line with cells belonging to this proc
	for (int y = proc_block_low[1]; y <= proc_block_high[1]; y++)
	{
		int x;
		// for each cell on this line before the first proc cell
		for (x = range[0][0]; x < proc_block_low[0]; x++)
		{
			// for particles in this cell
			for (int part = 0; part < npc; part++)
			{
				rand_norm(); rand_norm(); rand_norm();
			}
		}

		// for cells with x contained in this procs simulation zone
		for (x = proc_block_low[0]; x <= proc_block_high[0]; x++)
		{
			// if x is not inside the injection range
			if( x < range[0][0] || x > range[0][1] )
				continue;

			// for particles in this cell
			for (int part = 0; part < npc; part++, i++)
			{
				// save the values
				spec->part[i].ux = spec->ufl[0] + spec->uth[0] * rand_norm();
				spec->part[i].uy = spec->ufl[1] + spec->uth[1] * rand_norm(); 
				spec->part[i].uz = spec->ufl[2] + spec->uth[2] * rand_norm();
			}
		}

		// for each cell on this line, after the last proc cell
		for (x = proc_block_high[0] + 1; x <= range[0][1]; x++)
		{
			// for particles in this cell
			for (int part = 0; part < npc; part++)
			{
				rand_norm(); rand_norm(); rand_norm();
			}
		}
	}

	// for each line after the last line with cells belonging to this proc
	for (int y = proc_block_high[1] + 1; y <= range[1][1]; y++)
	{
		for (int x = range[0][0]; x <= range[0][1]; x++)
		{
			// for particle in this cell
			for (int part = 0; part < npc; part++)
			{
				rand_norm(); rand_norm(); rand_norm();
			}
		}
	}

	/* int index = 0;
	for (int y = 0; y < spec->nx_local[1]; y++)
	{
		for (int x = 0; x < spec->nx_local[0]; x++)
		{
			for (int i = 0; i < npc; i++)
			{
				printf("%f " , spec->part[index++].ux); fflush(stdout);
			}
		}
		printf("\n"); fflush(stdout);
	}
	printf("\n"); fflush(stdout);

	index = 0;
	for (int y = 0; y < spec->nx_local[1]; y++)
	{
		for (int x = 0; x < spec->nx_local[0]; x++)
		{
			for (int i = 0; i < npc; i++)
			{
				printf("%f " , spec->part[index++].uy); fflush(stdout);
			}
		}
		printf("\n"); fflush(stdout);
	}
	printf("\n"); fflush(stdout);

	index = 0;
	for (int y = 0; y < spec->nx_local[1]; y++)
	{
		for (int x = 0; x < spec->nx_local[0]; x++)
		{
			for (int i = 0; i < npc; i++)
			{
				printf("%f " , spec->part[index++].uz); fflush(stdout);
			}
		}
		printf("\n"); fflush(stdout);
	}
	printf("\n"); fflush(stdout); */

	// old implementation
	/* 
	for (i = start; i <= end; i++)
	{
		spec->part[i].ux = spec->ufl[0] + spec->uth[0] * rand_norm();
		spec->part[i].uy = spec->ufl[1] + spec->uth[1] * rand_norm(); 
		spec->part[i].uz = spec->ufl[2] + spec->uth[2] * rand_norm();
	}
	*/
}

void spec_set_x( t_species* spec, const int range[NUM_DIMS][NUM_DIMS] )
{
	int i, j, k;
	
	float* poscell;
	float start, end;
	
	// Calculate particle positions inside the cell
	const int npc = spec->ppc[0] * spec->ppc[1];
	const t_part_data dpcx = 1.0f / spec->ppc[0];
	const t_part_data dpcy = 1.0f / spec->ppc[1];
	
	poscell = malloc( 2 * npc * sizeof( t_part_data ) );
	int ip = 0;
	for (j = 0; j < spec->ppc[1]; j++)
	{
		for (i = 0; i < spec->ppc[0]; i++)
		{
			poscell[ip]   = dpcx * ( i + 0.5 );
			poscell[ip+1] = dpcy * ( j + 0.5 );
			ip+=2;
		}
	}

	ip = spec->np;
	
	// Set position of particles in the specified grid range according to the density profile
	switch ( spec->density.type )
	{
	case STEP: // Step like density profile
		
		// Get edge position normalized to cell size;
		start = spec->density.start / spec->dx[0] - spec->n_move;

		for (j = range[1][0]; j <= range[1][1]; j++)
		{
			for (i = range[0][0]; i <= range[0][1]; i++)
			{
				for (k = 0; k < npc; k++)
				{
					if ( proc_block_low[0] + i + poscell[2*k] > start )
					{
						spec->part[ip].ix = i;
						spec->part[ip].iy = j;
						spec->part[ip].x = poscell[2*k];
						spec->part[ip].y = poscell[2*k+1];
						ip++;
					}
				}
			}
		}
		break;

	case SLAB: // Slab like density profile
		
		// Get edge position normalized to cell size;
		start = spec->density.start / spec->dx[0] - spec->n_move;
		end   = spec->density.end / spec->dx[0] - spec->n_move;

		for (j = range[1][0]; j <= range[1][1]; j++)
		{
			for (i = range[0][0]; i <= range[0][1]; i++)
			{
				for (k = 0; k < npc; k++)
				{
					if ( proc_block_low[0] + i + poscell[2*k] > start && proc_block_high[0] + i + poscell[2*k] < end )
					{
						spec->part[ip].ix = i;
						spec->part[ip].iy = j;
						spec->part[ip].x = poscell[2*k];
						spec->part[ip].y = poscell[2*k+1];
						ip++;
					}
				}
			}
		}
		break;

	default: // Uniform density
		for (j = range[1][0]; j <= range[1][1]; j++)
		{
			for (i = range[0][0]; i <= range[0][1]; i++)
			{

				for (k=0; k<npc; k++)
				{
					spec->part[ip].ix = i;
					spec->part[ip].iy = j;
					spec->part[ip].x = poscell[2*k];
					spec->part[ip].y = poscell[2*k+1];
					ip++;
				}
			}
		}
	}
	
	spec->np = ip;
	
	free(poscell);	
}

// Check if buffer is large enough and if not reallocate
void check_part_buffer_size(t_species* spec, int np_inj)
{
	if ( spec->np + np_inj > spec->np_max )
	{
		spec->np_max = (( spec->np_max + np_inj )/1024 + 1) * 1024;
		spec->part = realloc( (void*) spec->part, spec->np_max * sizeof(t_part) );

		// check if realloc was successful
		assert(spec->part);
	}
}

void spec_inject_particles( t_species* spec, const int range[NUM_DIMS][NUM_DIMS], const int range_local[NUM_DIMS][NUM_DIMS])
{
	int start = spec->np;

	// Get maximum number of particles to inject, on this proc
	int np_inj = ( range_local[0][1] - range_local[0][0] + 1 ) * ( range_local[1][1] - range_local[1][0] + 1 ) * 
				spec->ppc[0] * spec->ppc[1];

	// Check if buffer is large enough and if not reallocate
	check_part_buffer_size(spec, np_inj);

	// Set particle positions
	spec_set_x( spec, range_local );
	
	// Set momentum of injected particles
	spec_set_u( spec, start, range );
}

void spec_new( t_species* spec, char name[], const t_part_data m_q, const int ppc[], 
				const t_part_data *ufl, const t_part_data * uth,
				const int nx[], t_part_data box[], const float dt, t_density* density )
{
	static int spec_id = 0;
	spec->id = spec_id++;

	assign_proc_blocks(nx);

	int i, npc;
	
	// Species name
	strncpy( spec->name, name, MAX_SPNAME_LEN );

	const int nx_local[NUM_DIMS] = {BLOCK_SIZE(proc_coords[0], dims[0], nx[0]), BLOCK_SIZE(proc_coords[1], dims[1], nx[1])};

	spec->nrow_local = gc[0][0] + nx_local[0] + gc[0][1];
	spec->nrow = gc[0][0] + nx[0] + gc[0][1];
	
	npc = 1;
	// Store species data
	for (i = 0; i < NUM_DIMS; i++)
	{
		spec->nx[i] = nx[i];
		spec->nx_local[i] = nx_local[i];

		spec->ppc[i] = ppc[i];
		npc *= ppc[i];
		
		spec->box[i] = box[i];
		spec->dx[i] = box[i] / nx[i];

		spec->box_local[i] = spec->dx[i] * nx_local[i];
	}
	
	spec->m_q = m_q;
	spec->q = copysign( 1.0f, m_q ) / npc;

	spec->dt = dt;
	
	// Initialize particle buffer
	spec->np_max = 0;
	spec->part = NULL;
	
	
	// Initialize density profile
	if ( density )
	{
		spec->density = *density;
		if ( spec->density.n == 0. ) spec->density.n = 1.0;
	}
	else
	{
		// Default values
		spec->density = (t_density) { .type = UNIFORM, .n = 1.0 };
	}

	// Initialize temperature profile
	if ( ufl )
	{
		for(i = 0; i < 3; i++)
			spec->ufl[i] = ufl[i];
	}
	else
	{
		for(i = 0; i < 3; i++)
			spec->ufl[i] = 0;
	}

	// Density multiplier
	spec ->q *= fabsf( spec->density.n );

	if ( uth )
	{
		for(i = 0; i < 3; i++)
			spec->uth[i] = uth[i];
	}
	else
	{
		for(i = 0; i < 3; i++)
			spec->uth[i] = 0;
	}

	// Reset iteration number
	spec->iter = 0;

	// Reset moving window information
	spec->moving_window = 0;
	spec->n_move = 0;

	// Inject initial particle distribution
	spec->np = 0;
	
	const int range[NUM_DIMS][NUM_DIMS] =  {{0, nx[0]-1},
											{0, nx[1]-1}};

	const int range_local[NUM_DIMS][NUM_DIMS] =  {{0, nx_local[0]-1},
												  {0, nx_local[1]-1}};

	spec_inject_particles( spec, range, range_local);
}

// Removes particles that have left the simulation space
void remove_outer_particles(t_species *spec)
{
	// If proc is on the left edge of the simulation space
	if (is_on_edge[0])
	{
		// Use absorbing boundaries along the left edge
		int i = 0;
		while ( i < spec->np )
		{
			if (spec->part[i].ix < 0)
			{
				spec->part[i] = spec->part[ --spec->np ];
				continue;
			}
			i++;
		}
	}

	// If proc is on the right edge of the simulation space
	if (is_on_edge[1])
	{
		const int max_ix = spec->nx_local[0] - 1;
		
		// Use absorbing boundaries along the right edge
		int i = 0;
		while ( i < spec->np )
		{
			if (spec->part[i].ix > max_ix)
			{
				spec->part[i] = spec->part[ --spec->np ];
				continue;
			}
			i++;
		}
	}
}

// Also removes particles that are outside the simulation space
void spec_move_window( t_species *spec )
{
	// If the window needs to be moved this iteration
	if ( (spec->iter * spec->dt ) > (spec->dx[0] * (spec->n_move + 1)) )
	{
		// Shift all particles 1 cell to the left
		for(int i = 0; i < spec->np; i++ )
		{
			spec->part[i].ix--;
		}

		// Remove particles that left the simulation space
		remove_outer_particles(spec);

		// Increase moving window counter
		spec->n_move++;

		// if proc is on the right edge of the simulation space
		if(is_on_edge[1])
		{
			// Inject particles on the right edge of the simulation box
			const int range[NUM_DIMS][NUM_DIMS] =	{{spec->nx[0]-1,spec->nx[0]-1},
													{			  0,spec->nx[1]-1}};

			const int range_local[NUM_DIMS][NUM_DIMS] =  {{spec->nx_local[0]-1, spec->nx_local[0]-1},
												  		  {					 0, spec->nx_local[1]-1}};

			spec_inject_particles(spec, range, range_local);
		}
	}
	else
	{
		// Remove particles that left the simulation space
		remove_outer_particles(spec);
	}
	
}

/*********************************************************************************************
 
 Current deposition
 
 *********************************************************************************************/

void dep_current_zamb(int ix, int iy, int di, int dj, 
					  float x0, float y0, float dx, float dy,
					  float qnx, float qny, float qvz,
					  t_current *current)
{
	// Split the particle trajectory
	
	typedef struct {
		float x0, x1, y0, y1, dx, dy, qvz;
		int ix, iy;
	} t_vp;
	
	t_vp vp[3];
	int vnp = 1;
	
	// split

	// old position
	vp[0].x0 = x0;
	vp[0].y0 = y0;

	// position delta
	vp[0].dx = dx;
	vp[0].dy = dy;
	
	// new position
	vp[0].x1 = x0+dx;
	vp[0].y1 = y0+dy;
	
	vp[0].qvz = qvz/2.0;

	vp[0].ix = ix;
	vp[0].iy = iy;
		
	// x split
	if ( di != 0 )
	{
		
		//int ib = ( di+1 )>>1;
		int ib = ( di == 1 );
		
		float delta = (x0+dx-ib)/dx;
		
		// Add new particle
		vp[1].x0 = 1-ib;
		vp[1].x1 = (x0 + dx) - di;
		vp[1].dx = dx*delta;
		vp[1].ix = ix + di;
		
		float ycross = y0 + dy*(1.0f-delta);

		vp[1].y0 = ycross; 
		vp[1].y1 = vp[0].y1;
		vp[1].dy = dy*delta;	
		vp[1].iy = iy;
		
		vp[1].qvz = vp[0].qvz*delta;
		
		// Correct previous particle
		vp[0].x1 = ib;
		vp[0].dx *= (1.0f-delta);
		
		vp[0].dy *= (1.0f-delta);
		vp[0].y1  = ycross;
		
		vp[0].qvz *= (1.0f-delta);
		
		vnp++;		
	}
	
	// ysplit
	if ( dj != 0 )
	{
		int isy = 1 - ( vp[0].y1<0.0f || vp[0].y1>=1.0f );
		
		// int jb = ( dj+1 )>>1; 
		int jb = (dj == 1);
		
		// The static analyser gets confused by this but it is correct
		float delta = (vp[isy].y1-jb)/vp[isy].dy;
		
		// Add new particle
		vp[vnp].y0 = 1-jb;
		vp[vnp].y1 = vp[isy].y1 - dj;
		vp[vnp].dy = vp[isy].dy*delta;
		vp[vnp].iy = vp[isy].iy + dj;
		
		float xcross = vp[isy].x0 + vp[isy].dx*(1.0f-delta); 
		
		vp[vnp].x0 = xcross;
		vp[vnp].x1 = vp[isy].x1;
		vp[vnp].dx = vp[isy].dx*delta;
		vp[vnp].ix = vp[isy].ix;

		vp[vnp].qvz = vp[isy].qvz*delta;
		
		// Correct previous particle
		vp[isy].y1  = jb;
		vp[isy].dy *= (1.0f-delta);
		
		vp[isy].dx *= (1.0f-delta);
		vp[isy].x1  = xcross;
		
		vp[isy].qvz *= (1.0f-delta);
		
		// Correct extra vp if needed
 		if ( isy < vnp -1) {
			vp[1].y0 -= dj;
			vp[1].y1 -= dj;
			vp[1].iy += dj;
		} 
		vnp++;		
	}

	// Deposit virtual particle currents
	int k;
	const int nrow = current->nrow_local; // Local nrow
	t_vfld* restrict const J = current->J;

	for (k = 0; k < vnp; k++)
	{
		float S0x[2], S1x[2], S0y[2], S1y[2];
		float wl1, wl2;
		float wp1[2],wp2[2];
		
		S0x[0] = 1.0f - vp[k].x0;
		S0x[1] = vp[k].x0;

		S1x[0] = 1.0f - vp[k].x1;
		S1x[1] = vp[k].x1;

		S0y[0] = 1.0f - vp[k].y0;
		S0y[1] = vp[k].y0;

		S1y[0] = 1.0f - vp[k].y1;
		S1y[1] = vp[k].y1;

		wl1 = qnx * vp[k].dx;
		wl2 = qny * vp[k].dy;
		
		wp1[0] = 0.5f*(S0y[0] + S1y[0]);
		wp1[1] = 0.5f*(S0y[1] + S1y[1]);
		
		wp2[0] = 0.5f*(S0x[0] + S1x[0]);
		wp2[1] = 0.5f*(S0x[1] + S1x[1]);
		
		J[ vp[k].ix + nrow*vp[k].iy     ].x += wl1 * wp1[0];
		J[ vp[k].ix + nrow*(vp[k].iy+1) ].x += wl1 * wp1[1];

		J[ vp[k].ix   + nrow*vp[k].iy ].y += wl2 * wp2[0];
		J[ vp[k].ix+1 + nrow*vp[k].iy ].y += wl2 * wp2[1];

		J[ vp[k].ix   + nrow*vp[k].iy    ].z += vp[k].qvz * (S0x[0]*S0y[0]+S1x[0]*S1y[0]+(S0x[0]*S1y[0]-S1x[0]*S0y[0])/2.0f);
		J[ vp[k].ix+1 + nrow*vp[k].iy    ].z += vp[k].qvz * (S0x[1]*S0y[0]+S1x[1]*S1y[0]+(S0x[1]*S1y[0]-S1x[1]*S0y[0])/2.0f);
		J[ vp[k].ix   + nrow*(vp[k].iy+1)].z += vp[k].qvz * (S0x[0]*S0y[1]+S1x[0]*S1y[1]+(S0x[0]*S1y[1]-S1x[0]*S0y[1])/2.0f);
		J[ vp[k].ix+1 + nrow*(vp[k].iy+1)].z += vp[k].qvz * (S0x[1]*S0y[1]+S1x[1]*S1y[1]+(S0x[1]*S1y[1]-S1x[1]*S0y[1])/2.0f);
	}
}


/*********************************************************************************************
 
 Particle advance
 
 *********************************************************************************************/

void interpolate_fld( const t_vfld* restrict const E, const t_vfld* restrict const B, const int nrow, 
			  const t_part* restrict const part, t_vfld* restrict const Ep, t_vfld* restrict const Bp )
{
	register int i, j, ih, jh;
	register t_fld w1, w2, w1h, w2h;
	
	i = part->ix;
	j = part->iy;
	
	w1 = part->x;
	w2 = part->y;
	
	ih = (w1 < 0.5f)? -1 : 0;
	jh = (w2 < 0.5f)? -1 : 0;
	
	// w1h = w1 - 0.5f - ih;
	// w2h = w2 - 0.5f - jh;
	w1h = w1 + ((w1 < 0.5f) ? 0.5f : -0.5f);
	w2h = w2 + ((w2 < 0.5f) ? 0.5f : -0.5f);

	
	ih += i;
	jh += j;
	
	Ep->x = ( E[ih +     j *nrow].x * (1.0f - w1h) + E[ih+1 +     j*nrow].x * w1h ) * (1.0f -  w2 ) +
			( E[ih + (j +1)*nrow].x * (1.0f - w1h) + E[ih+1 + (j+1)*nrow].x * w1h ) * w2;

	Ep->y = ( E[i  +     jh*nrow].y * (1.0f -  w1) + E[i+1  +     jh*nrow].y * w1 ) * (1.0f - w2h ) +
			( E[i  + (jh+1)*nrow].y * (1.0f -  w1) + E[i+1  + (jh+1)*nrow].y * w1 ) * w2h;

	Ep->z = ( E[i  +     j *nrow].z * (1.0f - w1) + E[i+1 +     j*nrow].z * w1 ) * (1.0f - w2 ) +
			( E[i  + (j +1)*nrow].z * (1.0f - w1) + E[i+1 + (j+1)*nrow].z * w1 ) * w2;

	Bp->x = ( B[i  +     jh*nrow].x * (1.0f - w1) + B[i+1 +     jh*nrow].x * w1 ) * (1.0f - w2h ) +
			( B[i  + (jh+1)*nrow].x * (1.0f - w1) + B[i+1 + (jh+1)*nrow].x * w1 ) * w2h;

	Bp->y = ( B[ih +     j *nrow].y * (1.0f - w1h) + B[ih+1 +     j*nrow].y * w1h ) * (1.0f - w2 ) +
			( B[ih + (j +1)*nrow].y * (1.0f - w1h) + B[ih+1 + (j+1)*nrow].y * w1h ) * w2;

	Bp->z = ( B[ih +     jh*nrow].z * (1.0f - w1h) + B[ih+1 +     jh*nrow].z * w1h ) * (1.0f - w2h ) +
			( B[ih + (jh+1)*nrow].z * (1.0f - w1h) + B[ih+1 + (jh+1)*nrow].z * w1h ) * w2h;
}

int ltrim( t_part_data x )
{
	return ( x >= 1.0f ) - ( x < 0.0f );
}

// Wait for each write and save particles
void wait_save_particles(t_species* species_array, const int num_spec)
{
	int num_new_part[num_spec][NUM_ADJ];
	int num_new_part_spec[num_spec]; memset(num_new_part_spec, 0, num_spec * sizeof(int));

	const char moving_window = species_array->moving_window;

	// -1 because it was incremented before in spec_advance
	const int iter_num = species_array->iter - 1;

	// Segment offset multiplier, 0 if iteration is even, 2 if odd.
	const int seg_offset_mult = (iter_num % 2) == 0 ? 0 : 2;

	for (int spec_i = 0; spec_i < num_spec; spec_i++)
	{
		// printf("Receiving particles from species %d\n\n", spec_i); fflush(stdout);

		// Index of first received particle of this species on each segment.
		int spec_starting_index[NUM_ADJ];
		for (int dir = 0; dir < NUM_ADJ; dir++)
		{
			spec_starting_index[dir] = part_send_seg_size[dir] + seg_offset_mult * part_send_seg_size[dir];
		}
		
		for (int spec_i_2 = 0; spec_i_2 < spec_i; spec_i_2++)
		{
			for (int dir = 0; dir < NUM_ADJ; dir++)
			{
				spec_starting_index[dir] += num_new_part[spec_i_2][dir] + 1;
			}
		}

		// notif id depends in iteration num, if odd notif id = spec_id + num_spec, if even notif_id = spec_id
		const gaspi_notification_id_t notif_id = (iter_num % 2) == 0 ? spec_i : num_spec + spec_i;

		for (int dir = 0; dir < NUM_ADJ; dir++)
		{
			// Check if we will receive particles from this dir
			if ( !can_send_to_dir(moving_window, dir) )
				continue;

			gaspi_notification_id_t id;
			SUCCESS_OR_DIE( gaspi_notify_waitsome(
			dir,			// The segment id, = to direction
			notif_id,		// The notification id to wait for
			1,				// The number of notification ids this wait will accept, waiting for a specific write, so 1
			&id,			// Output parameter with the id of a received notification
			GASPI_BLOCK		// Timeout in milliseconds, wait until write is completed
			));

			gaspi_notification_t value;
			SUCCESS_OR_DIE( gaspi_notify_reset(dir, id, &value) );

			// printf("Starting index for spec %d dir %d is %d\n", spec_i, dir, spec_starting_index[dir]);
			// Number of particles received is saved in the ix field of the first particle of that species
			// -1 because of the fake particle
			const int num_part = particle_segments[dir][spec_starting_index[dir]].ix - 1;

			// printf("from dir %d fake part %d %d\n", dir, particle_segments[dir][spec_starting_index[dir]]. ix,particle_segments[dir][spec_starting_index[dir]].iy);

			num_new_part_spec[spec_i] += num_part;
			num_new_part[spec_i][dir] = num_part;
		}

		// if each write to this species was just a fake particle, continue
		if (num_new_part_spec[spec_i] == 0) continue;

		// printf("NEW PARTICLES!\n"); fflush(stdout);

		// realloc particle array if necessary
		check_part_buffer_size(&species_array[spec_i], num_new_part_spec[spec_i]);

		// Copy particles from each segment to the main particle array for this species
		int copy_index = species_array[spec_i].np;
		for (int dir = 0; dir < NUM_ADJ; dir++)
		{
			// Check if we received particles from this dir
			if ( !can_send_to_dir(moving_window, dir) )
				continue;

			const int num_part = num_new_part[spec_i][dir];
			const int starting_index = spec_starting_index[dir] + 1;

			// printf("dir: %d, num_part:%d\n", dir, num_part);
			// for (int i = 0; i < num_part; i++)
			// {
			// 	t_part part = particle_segments[dir][starting_index + i];
			// 	printf("New Part: ix:%d, iy:%d, x:%f, y:%f, ux:%f, uy:%f, uz:%f\n", part.ix, part.iy, part.x, part.y, part.ux, part.uy, part.uz); fflush(stdout);
			// }

			memcpy(&species_array[spec_i].part[copy_index], &particle_segments[dir][starting_index], num_part * sizeof(t_part));
			copy_index += num_part;
		}

		// update the number of particles of this species
		species_array[spec_i].np += num_new_part_spec[spec_i];
	}

	// // Print particles
	// for (int i = 0; i < num_spec; i++)
	// {
	// 	printf("part from species %d: %d\n", i, species_array[i].np);
	// 	for (int j = 0; j < species_array[i].np; j++)
	// 	{
	// 		t_part part = species_array[i].part[j];
	// 		printf("Part: ix:%d, iy:%d, x:%f, y:%f, ux:%f, uy:%f, uz:%f\n", part.ix, part.iy, part.x, part.y, part.ux, part.uy, part.uz); fflush(stdout);
	// 	}
	// }
}


// get the direction from which the particle left the proc zone, based on its cell coordinates
// If particle does no leave this proc, return -1
int get_part_seg_direction(const t_part* const part_pointer, const int nx_local[NUM_DIMS])
{
	const int ix = part_pointer->ix;
	const int iy = part_pointer->iy;

	// Left border
	if (ix < 0)
	{
		// Top left border
		if( iy < 0 )
			return UP_LEFT;
		
		// Bottom left border
		if ( iy >= nx_local[1] )
			return DOWN_LEFT;

		return LEFT;
	}
	
	// Right border
	if (ix >= nx_local[0])
	{
		// Top right border
		if( iy < 0 )
			return UP_RIGHT;
		
		// Bottom right border
		if ( iy >= nx_local[1] )
			return DOWN_RIGHT;

		return RIGHT;
	}

	// Top border
	if( iy < 0 )
		return UP;

	// Bottom border
	if ( iy >= nx_local[1] )
		return DOWN;

	return -1;
}

// Correct particle coords to coords relative to the new proc
void correct_coords(t_part* const part_pointer, const int dir)
{
	// Correct part x cell coord
	// If particle leaves through the left border
	if (dir == LEFT || dir == UP_LEFT || dir == DOWN_LEFT)
	{
		part_pointer->ix = neighbour_nx[dir][0] - 1;
	}
	// If particle leaves through the right border
	else if (dir == RIGHT || dir == UP_RIGHT || dir == DOWN_RIGHT)
	{
		part_pointer->ix = 0;
	}

	// Correct part y cell coord
	// If particle leaves through the bottom border
	if (dir == DOWN || dir == DOWN_LEFT || dir == DOWN_RIGHT)
	{
		part_pointer->iy = 0;
	}
	// If particle leaves through the top border
	else if (dir == UP || dir == UP_LEFT || dir == UP_RIGHT)
	{
		part_pointer->iy = neighbour_nx[dir][1] - 1;
	}
}

void send_spec(t_species* spec, int part_seg_write_index[NUM_ADJ], int num_part_to_send[][NUM_ADJ], const int num_spec)
{	
	const int spec_id = spec->id;

	// index of the fake particle on each direction for this species
	int fake_part_index[NUM_ADJ];
	
	// Add fake particle to the start of the particles of this species on each segment.
	// This particle will be used to transmit the number of particles on this write.
	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		// Check if we can send particles to this dir
		if ( !can_send_to_dir(spec->moving_window, dir) )
			continue;

		// these values are set to make this particle easy to identify on debugging situations
		const t_part fake_part =
		{
			//.ix = -42,
			.iy = -42,
			.x = 0.42,
			.y = 0.42,
			.ux = 0.42,
			.uy = 0.42,
			.uz = 0.42
		};

		// Save the index of the fake particle
		fake_part_index[dir] = part_seg_write_index[dir];

		// Copy particle to the correct send segment
		particle_segments[dir][ part_seg_write_index[dir]++ ] = fake_part;

		// Increment part send count for this segment, for this species
		num_part_to_send[spec_id][dir]++;
	}

	// Check if particle left the proc zone, if so, copy it to the correct send segment
	int i = 0;
	while ( i < spec->np )
	{
		
		// Get the direction this particle left the local simulation space from, also correct its coords relative to the new proc
		const int dir = get_part_seg_direction(&spec->part[i], spec->nx_local);

		// check if particle left the proc zone
		if ( dir != -1 )
		{
			// int old_x = spec->part[i].ix;
			// int old_y = spec->part[i].iy;

			// Correct part coords relative to new proc
			correct_coords(&spec->part[i], dir);

			// Copy particle to segment
			particle_segments[dir][ part_seg_write_index[dir]++ ] = spec->part[i];

			// if (old_x != particle_segments[dir][ part_seg_write_index[dir]-1 ].ix || old_y != particle_segments[dir][part_seg_write_index[dir]-1].iy)
			// {
			// 	printf("old part x:%d, y:%d\n", old_x, old_y);
			// 	printf("NEW part x:%d, y:%d\n\n", particle_segments[dir][ part_seg_write_index[dir]-1 ].ix, particle_segments[dir][ part_seg_write_index[dir]-1 ].iy);
			// }

			// Increment part send count for this segment, for this species
			num_part_to_send[spec_id][dir]++;

			// Remove particle
			spec->part[i] = spec->part[ --spec->np ];
			continue;
		}

		i++;
	}

	// Check if segment size is respected before sending
	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		// Compute the total number of particles sent to this dir
		unsigned int num_part_seg = 0;
		for (int spec_i_2 = 0; spec_i_2 <= spec_id; spec_i_2++)
		{
			num_part_seg += num_part_to_send[spec_i_2][dir];
		}

		// printf("ON TOTAL sending %5d particles to proc %d to dir %d\n", num_part_seg, neighbour_rank[dir], OPPOSITE_DIR(dir));
		// printf("We have room to send %d particles\n\n", part_send_seg_size[dir]);

		assert(num_part_seg <= (unsigned int) part_send_seg_size[dir]);
	}

	// notif id depends in iteration num, if odd notif id = spec_id + num_spec, if even notif_id = spec_id
	// iter - 1 because it was incremented in spec_advance
	const gaspi_notification_id_t notif_id = ( (spec->iter-1) % 2) == 0 ? spec_id : num_spec + spec_id;

	// Make sure it is safe to change the segment data
	SUCCESS_OR_DIE( gaspi_wait(Q_PARTICLES, GASPI_BLOCK) );

	// printf("Sending particles from species %d\n", spec_id); fflush(stdout);

	// Send particles on each segment
	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		// Check if we can send particles to this dir
		if ( !can_send_to_dir(spec->moving_window, dir) )
			continue;

		// if a particle leaves a proc zone by moving right, it will be written to the PAR_LEFT segment of the receiving proc
		int new_dir = OPPOSITE_DIR(dir);

		gaspi_offset_t local_offset = fake_part_index[dir] * sizeof(t_part); // in bytes
		gaspi_offset_t remote_offset = (fake_part_index[dir] + part_send_seg_size[dir]) * sizeof(t_part); // in bytes
		gaspi_size_t size = num_part_to_send[spec_id][dir] * sizeof(t_part); // in bytes

		// Save the number of particles sent (including fake part) to the ix field of the fake particle, for this species, for this dir
		particle_segments[dir][ fake_part_index[dir] ].ix = num_part_to_send[spec_id][dir];

		SUCCESS_OR_DIE( gaspi_write_notify(
		dir,					// The segment id where data is located.
		local_offset,			// The offset where the data is located.
		neighbour_rank[dir],	// The rank where to write and notify.
		new_dir, 				// The remote segment id to write the data to.
		remote_offset,			// The remote offset where to write to.
		size,					// The size of the data to write.
		notif_id,				// The notification id to use.
		1,						// The notification value used.
		Q_PARTICLES,			// The queue where to post the request.
		GASPI_BLOCK				// Timeout in milliseconds.
		));

		// if (num_part_to_send[spec_id][dir] == 1)
		// {
		// 	continue;
		// }
		
		// for (int i = 0; i < num_part_to_send[spec_id][dir]; i++)
		// {
		// 	t_part part = particle_segments[dir][fake_part_index[dir] + i];
		// 	printf("Sent particle ix:%d, iy:%d, x:%f, y:%f, ux:%f, uy:%f, uz:%f\n", part.ix, part.iy, part.x, part.y, part.ux, part.uy, part.uz);
		// }
		// printf("\n"); fflush(stdout);
	}
}

void spec_advance(t_species* spec, t_emf* emf, t_current* current)
{
	// uint64_t t0 = timer_ticks();
	
	const t_part_data tem   = 0.5 * spec->dt / spec->m_q;
	const t_part_data dt_dx = spec->dt / spec->dx[0];
	const t_part_data dt_dy = spec->dt / spec->dx[1];

	// Auxiliary values for current deposition
	const t_part_data qnx = spec->q * spec->dx[0] / spec->dt;
	const t_part_data qny = spec->q * spec->dx[1] / spec->dt;

	// Advance particles
	for (int i = 0; i < spec->np; i++)
	{
		t_vfld Ep, Bp;
		t_part_data utx, uty, utz;
		t_part_data ux, uy, uz, rg;
		t_part_data gtem, otsq;
		
		t_part_data x1, y1;
		
		int di, dj;
		float dx, dy;

		// Load particle momenta
		ux = spec->part[i].ux;
		uy = spec->part[i].uy;
		uz = spec->part[i].uz;

		// interpolate fields ONLY CHANGES LOCAL VARIABLES
		interpolate_fld( emf->E, emf->B, emf->nrow_local, &spec->part[i], &Ep, &Bp );
		
		// advance u using Boris scheme
		Ep.x *= tem;
		Ep.y *= tem;
		Ep.z *= tem;
		
		utx = ux + Ep.x;
		uty = uy + Ep.y;
		utz = uz + Ep.z;

		// Perform first half of the rotation
		gtem = tem / sqrtf( 1.0f + utx*utx + uty*uty + utz*utz );
		
		Bp.x *= gtem;
		Bp.y *= gtem;
		Bp.z *= gtem;

		otsq = 2.0f / ( 1.0f + Bp.x*Bp.x + Bp.y*Bp.y + Bp.z*Bp.z );

		ux = utx + uty*Bp.z - utz*Bp.y;
		uy = uty + utz*Bp.x - utx*Bp.z;
		uz = utz + utx*Bp.y - uty*Bp.x;
		
		// Perform second half of the rotation
		
		Bp.x *= otsq;
		Bp.y *= otsq;
		Bp.z *= otsq;
		
		utx += uy*Bp.z - uz*Bp.y;
		uty += uz*Bp.x - ux*Bp.z;
		utz += ux*Bp.y - uy*Bp.x;
		
		// Perform second half of electric field acceleration
		ux = utx + Ep.x;
		uy = uty + Ep.y;
		uz = utz + Ep.z;
		
		// Store new momenta
		spec->part[i].ux = ux;
		spec->part[i].uy = uy;
		spec->part[i].uz = uz;

		if (ux == 0.0f && uy == 0.0f && uz == 0.0f)
		{
			continue;
		}
		
		// push particle
		rg = 1.0f / sqrtf(1.0f + ux*ux + uy*uy + uz*uz);
				
		dx = dt_dx * rg * ux;
		dy = dt_dy * rg * uy;
		
		x1 = spec->part[i].x + dx; 
		y1 = spec->part[i].y + dy;
		
		di = ltrim(x1);
		dj = ltrim(y1);

		x1 -= di;
		y1 -= dj;
		
		t_part_data qvz = spec->q * uz * rg;

		dep_current_zamb(spec->part[i].ix, spec->part[i].iy, di, dj, 
						 spec->part[i].x , spec->part[i].y , dx, dy, 
						 qnx, qny, qvz, current );

		// Store results
		spec->part[i].x = x1;
		spec->part[i].y = y1;
		spec->part[i].ix += di;
		spec->part[i].iy += dj;
	}

	// Advance internal iteration number
	spec->iter += 1;
	// _spec_npush += spec->np;

	// Move simulation window if needed
	if ( spec->moving_window )
	{
		spec_move_window( spec );
	}
	
	//_spec_time += timer_interval_seconds( t0, timer_ticks() );
}

/*********************************************************************************************
 
 Charge Deposition
 
 *********************************************************************************************/


void spec_deposit_charge( const t_species* spec, t_part_data* charge )
{
	int i,j;
	
	// Charge array is expected to have 1 guard cell at the upper boundary
	int nrow = spec->nx[0] + 1;
	t_part_data q = spec->q;
	
	for (i=0; i<spec->np; i++) {
		int idx = spec->part[i].ix + nrow*spec->part[i].iy;
		t_fld w1, w2;
		
		w1 = spec->part[i].x;
		w2 = spec->part[i].y;
		
		charge[ idx            ] += ( 1.0f - w1 ) * ( 1.0f - w2 ) * q;
		charge[ idx + 1        ] += (        w1 ) * ( 1.0f - w2 ) * q;
		charge[ idx     + nrow ] += ( 1.0f - w1 ) * (        w2 ) * q;
		charge[ idx + 1 + nrow ] += (        w1 ) * (        w2 ) * q;
	}

	// Correct boundary values

	// x
	if ( ! spec->moving_window )
	{
		for (j = 0; j < spec->nx[1] + 1; j++)
		{
			charge[ 0 + j*nrow ] += charge[ spec->nx[0] + j*nrow ];
		}
	}
	
	// y - Periodic boundaries
	for (i = 0; i < spec->nx[0]+1; i++)
	{
		charge[ i + 0 ] += charge[ i + spec->nx[1] * nrow ];
	}

}

/*********************************************************************************************
 
 Diagnostics
 
 *********************************************************************************************/

void spec_rep_particles( const t_species *spec )
{
	
	t_zdf_file part_file;

	int i;
	
	const char * quants[] = {
		"x1","x2",
		"u1","u2","u3"
	};

	const char * units[] = {
		"c/\\omega_p", "c/\\omega_p",
		"c","c","c"
	};

	t_zdf_iteration iter = {
		.n = spec->iter,
		.t = spec->iter * spec->dt,
		.time_units = "1/\\omega_p"
	};

	// Allocate buffer for positions
	
	t_zdf_part_info info = {
		.name = (char *) spec->name,
		.nquants = 5,
		.quants = (char **) quants,
		.units = (char **) units,
		.np = spec ->np
	};

	// Create file and add description
	zdf_part_file_open( &part_file, &info, &iter, "PARTICLES" );

	// Add positions and generalized velocities
	size_t size = ( spec->np ) * sizeof( float );
	float* data = malloc( size );

	// x1
	for( i = 0; i < spec ->np; i++ )
		data[i] = ( spec->n_move + spec->part[i].ix + spec->part[i].x ) * spec->dx[0];
	zdf_part_file_add_quant( &part_file, quants[0], data, spec ->np );

	// x2
	for( i = 0; i < spec ->np; i++ )
		data[i] = (spec->part[i].iy + spec->part[i].y ) * spec->dx[1];
	zdf_part_file_add_quant( &part_file, quants[1], data, spec ->np );

	// ux
	for( i = 0; i < spec ->np; i++ ) data[i] = spec->part[i].ux;
	zdf_part_file_add_quant( &part_file, quants[2], data, spec ->np );

	// uy
	for( i = 0; i < spec ->np; i++ ) data[i] = spec->part[i].uy;
	zdf_part_file_add_quant( &part_file, quants[3], data, spec ->np );

	// uz
	for( i = 0; i < spec ->np; i++ ) data[i] = spec->part[i].uz;
	zdf_part_file_add_quant( &part_file, quants[4], data, spec ->np );

	free( data );

	zdf_close_file( &part_file );
}	


void spec_rep_charge( const t_species *spec )
{
	t_part_data *buf, *charge, *b, *c;
	size_t size;
	int i, j;
	
	// Add 1 guard cell to the upper boundary
	size = ( spec->nx[0] + 1 ) * ( spec->nx[1] + 1 ) * sizeof( t_part_data );
	charge = malloc( size );
	memset( charge, 0, size );
	
	// Deposit the charge
	spec_deposit_charge( spec, charge );
	
	// Compact the data to save the file (throw away guard cells)
	size = ( spec->nx[0] ) * ( spec->nx[1] );
	buf = malloc( size * sizeof( float ) );
	
	b = buf;
	c = charge;
	for( j = 0; j < spec->nx[1]; j++) {
		for ( i = 0; i < spec->nx[0]; i++ ) {
			b[i] = c[i];
		}
		b += spec->nx[0];
		c += spec->nx[0] + 1;
	}
	
	free( charge );

	t_zdf_grid_axis axis[2];
	axis[0] = (t_zdf_grid_axis) {
		.min = 0.0,
		.max = spec->box[0],
		.label = "x_1",
		.units = "c/\\omega_p"
	};

	axis[1] = (t_zdf_grid_axis) {
		.min = 0.0,
		.max = spec->box[1],
		.label = "x_2",
		.units = "c/\\omega_p"
	};

	t_zdf_grid_info info = {
		.ndims = 2,
		.label = "charge",
		.units = "n_e",
		.axis  = axis
	};

	info.nx[0] = spec->nx[0];
	info.nx[1] = spec->nx[1];

	t_zdf_iteration iter = {
		.n = spec->iter,
		.t = spec->iter * spec->dt,
		.time_units = "1/\\omega_p"
	};

	zdf_save_grid( buf, &info, &iter, spec->name );	


	free( buf );
}	


void spec_pha_axis( const t_species *spec, int i0, int np, int quant, float *axis )
{
	int i;
	
	switch (quant) {
		case X1:
			for (i = 0; i < np; i++) 
				axis[i] = ( spec->part[i0+i].x + spec->part[i0+i].ix ) * spec->dx[0];
			break;
		case X2:
			for (i = 0; i < np; i++) 
				axis[i] = ( spec->part[i0+i].y + spec->part[i0+i].iy ) * spec->dx[1];
			break;
		case U1:
			for (i = 0; i < np; i++) 
				axis[i] = spec->part[i0+i].ux;
			break;
		case U2:
			for (i = 0; i < np; i++) 
				axis[i] = spec->part[i0+i].uy;
			break;
		case U3:
			for (i = 0; i < np; i++) 
				axis[i] = spec->part[i0+i].uz;
			break;
	}
}

const char * spec_pha_axis_units( int quant ) {
	switch (quant) {
		case X1:
		case X2:
			return("c/\\omega_p");
			break;
		case U1:
		case U2:
		case U3:
			return("m_e c");
	}
	return("");
}


void spec_deposit_pha( const t_species *spec, const int rep_type,
			  const int pha_nx[], const float pha_range[][2], float* restrict buf )
{
	const int BUF_SIZE = 1024;
	float pha_x1[BUF_SIZE], pha_x2[BUF_SIZE];


	const int nrow = pha_nx[0];

	const int quant1 = rep_type & 0x000F;
	const int quant2 = (rep_type & 0x00F0)>>4;

	const float x1min = pha_range[0][0];
	const float x2min = pha_range[1][0];

	const float rdx1 = pha_nx[0] / ( pha_range[0][1] - pha_range[0][0] );
	const float rdx2 = pha_nx[1] / ( pha_range[1][1] - pha_range[1][0] );

	for ( int i = 0; i<spec->np; i+=BUF_SIZE ) {
		int np = ( i + BUF_SIZE > spec->np )? spec->np - i : BUF_SIZE;

		spec_pha_axis( spec, i, np, quant1, pha_x1 );
		spec_pha_axis( spec, i, np, quant2, pha_x2 );

		for ( int k = 0; k < np; k++ ) {

			float nx1 = ( pha_x1[k] - x1min ) * rdx1;
			float nx2 = ( pha_x2[k] - x2min ) * rdx2;

			int i1 = (int)(nx1 + 0.5f);
			int i2 = (int)(nx2 + 0.5f);

			float w1 = nx1 - i1 + 0.5f;
			float w2 = nx2 - i2 + 0.5f;

			int idx = i1 + nrow*i2;

			if ( i2 >= 0 && i2 < pha_nx[1] ) {

				if (i1 >= 0 && i1 < pha_nx[0]) {
					buf[ idx ] += (1.0f-w1)*(1.0f-w2)*spec->q;
				}

				if (i1+1 >= 0 && i1+1 < pha_nx[0] ) {
					buf[ idx + 1 ] += w1*(1.0f-w2)*spec->q;
				}
			}

			idx += nrow;
			if ( i2+1 >= 0 && i2+1 < pha_nx[1] ) {

				if (i1 >= 0 && i1 < pha_nx[0]) {
					buf[ idx ] += (1.0f-w1)*w2*spec->q;
				}

				if (i1+1 >= 0 && i1+1 < pha_nx[0] ) {
					buf[ idx + 1 ] += w1*w2*spec->q;
				}
			}

		}

	}
}

void spec_rep_pha( const t_species *spec, const int rep_type,
			  const int pha_nx[], const float pha_range[][2] )
{

	char const * const pha_ax_name[] = {"x1","x2","x3","u1","u2","u3"};
	char pha_name[64];

	// Allocate phasespace buffer
	float* restrict buf = malloc( pha_nx[0] * pha_nx[1] * sizeof( float ));
	memset( buf, 0, pha_nx[0] * pha_nx[1] * sizeof( float ));

	// Deposit the phasespace
	spec_deposit_pha( spec, rep_type, pha_nx, pha_range, buf );

	// save the data in hdf5 format
	int quant1 = rep_type & 0x000F;
	int quant2 = (rep_type & 0x00F0)>>4;

	const char * pha_ax1_units = spec_pha_axis_units(quant1);
	const char * pha_ax2_units = spec_pha_axis_units(quant2);

	sprintf( pha_name, "%s%s", pha_ax_name[quant1-1], pha_ax_name[quant2-1] );

	t_zdf_grid_axis axis[2];
	axis[0] = (t_zdf_grid_axis) {
		.min = pha_range[0][0],
		.max = pha_range[0][1],
		.label = (char *) pha_ax_name[ quant1 - 1 ],
		.units = (char *) pha_ax1_units
	};

	axis[1] = (t_zdf_grid_axis) {
		.min = pha_range[1][0],
		.max = pha_range[1][1],
		.label = (char *) pha_ax_name[ quant2 - 1 ],
		.units = (char *) pha_ax2_units
	};

	t_zdf_grid_info info = {
		.ndims = 2,
		.label = pha_name,
		.units = "a.u.",
		.axis  = axis
	};

	info.nx[0] = pha_nx[0];
	info.nx[1] = pha_nx[1];

	t_zdf_iteration iter = {
		.n = spec->iter,
		.t = spec->iter * spec->dt,
		.time_units = "1/\\omega_p"
	};

	zdf_save_grid( buf, &info, &iter, spec->name );

	// Free temp. buffer
	free( buf );

}

void spec_report( const t_species *spec, const int rep_type, 
				  const int pha_nx[], const float pha_range[][2] )
{
	
	switch (rep_type & 0xF000) {
		case CHARGE:
			spec_rep_charge( spec );
			break;

		case PHA:
			spec_rep_pha( spec, rep_type, pha_nx, pha_range );
			break;

		case PARTICLES:
			spec_rep_particles( spec );
			break;
	}
	
	
}
