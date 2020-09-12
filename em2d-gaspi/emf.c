/*
 *  emf.c
 *  zpic
 *
 *  Created by Ricardo Fonseca on 10/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "emf.h"
#include "zdf.h"
#include "timer.h"

#include <GASPI.h>

static double _emf_time = 0.0;

extern const int gc[NUM_DIMS][NUM_DIMS];
extern int proc_coords[NUM_DIMS];
extern int dims[NUM_DIMS];

extern int proc_block_low[NUM_DIMS];
extern int proc_block_high[NUM_DIMS];
extern int proc_block_size[NUM_DIMS];

extern int emf_seg_size[NUM_ADJ][NUM_DIMS];
extern int emf_cell_to_send_starting_coord[NUM_ADJ][NUM_DIMS];
extern int emf_cell_to_write_starting_coord[NUM_ADJ][NUM_DIMS];

extern gaspi_rank_t neighbour_rank[NUM_ADJ];

extern gaspi_rank_t proc_rank;

extern t_vfld* emf_segments[NUM_ADJ];
extern t_vfld* emf_b_segments[NUM_ADJ];
extern t_vfld* emf_e_segments[NUM_ADJ];

void print_emf_e(t_emf* emf)
{
	printf("EMF E X\n");
	for (int y = -gc[1][0]; y < emf->nx_local[1] + gc[1][1]; y++)
	{
		if (y == 0 || y == emf->nx_local[1])
		{
			printf("\n");
		}
		
		for (int x = -gc[0][0]; x < emf->nx_local[0] + gc[0][1]; x++)
		{
			if (x == 0 || x == emf->nx_local[0])
			{
				printf("    ");
			}
			
			printf("%f ", emf->E[y * emf->nrow_local + x].x);
		}
		printf("\n");
	}
	printf("EMF E Y\n");
	for (int y = -gc[1][0]; y < emf->nx_local[1] + gc[1][1]; y++)
	{
		if (y == 0 || y == emf->nx_local[1])
		{
			printf("\n");
		}
		
		for (int x = -gc[0][0]; x < emf->nx_local[0] + gc[0][1]; x++)
		{
			if (x == 0 || x == emf->nx_local[0])
			{
				printf("    ");
			}
			
			printf("%f ", emf->E[y * emf->nrow_local + x].y);
		}
		printf("\n");
	}
	printf("EMF E Z\n");
	for (int y = -gc[1][0]; y < emf->nx_local[1] + gc[1][1]; y++)
	{
		if (y == 0 || y == emf->nx_local[1])
		{
			printf("\n");
		}
		
		for (int x = -gc[0][0]; x < emf->nx_local[0] + gc[0][1]; x++)
		{
			if (x == 0 || x == emf->nx_local[0])
			{
				printf("    ");
			}
			
			printf("%f ", emf->E[y * emf->nrow_local + x].z);
		}
		printf("\n");
	}
	fflush(stdout);
}

void print_emf_b(t_emf* emf)
{
	printf("EMF B X\n");
	for (int y = -gc[1][0]; y < emf->nx_local[1] + gc[1][1]; y++)
	{
		if (y == 0 || y == emf->nx_local[1])
		{
			printf("\n");
		}
		
		for (int x = -gc[0][0]; x < emf->nx_local[0] + gc[0][1]; x++)
		{
			if (x == 0 || x == emf->nx_local[0])
			{
				printf("    ");
			}
			
			printf("%f ", emf->B[y * emf->nrow_local + x].x);
		}
		printf("\n");
	}
	printf("EMF B Y\n");
	for (int y = -gc[1][0]; y < emf->nx_local[1] + gc[1][1]; y++)
	{
		if (y == 0 || y == emf->nx_local[1])
		{
			printf("\n");
		}
		
		for (int x = -gc[0][0]; x < emf->nx_local[0] + gc[0][1]; x++)
		{
			if (x == 0 || x == emf->nx_local[0])
			{
				printf("    ");
			}
			
			printf("%f ", emf->B[y * emf->nrow_local + x].y);
		}
		printf("\n");
	}
	printf("EMF B Z\n");
	for (int y = -gc[1][0]; y < emf->nx_local[1] + gc[1][1]; y++)
	{
		if (y == 0 || y == emf->nx_local[1])
		{
			printf("\n");
		}
		
		for (int x = -gc[0][0]; x < emf->nx_local[0] + gc[0][1]; x++)
		{
			if (x == 0 || x == emf->nx_local[0])
			{
				printf("    ");
			}
			
			printf("%f ", emf->B[y * emf->nrow_local + x].z);
		}
		printf("\n");
	}
}

double emf_time( void )
{
	return _emf_time;
}

// size_x moving window iterations difference to normal iterations, perspective of the sender
inline int get_moving_iter_size_x_diff(const int dir)
{
	/* 
	int moving_size_diffs[NUM_ADJ][NUM_DIMS] = // {size_diff_x, size_diff_y}
	{
		//LEFT		CENTER		RIGHT
		{1, 0},		{0, 0},		{NULL, NULL},	// DOWN !!!
		{1, 0},					{NULL, NULL},	// CENTER
		{1, 0},		{0, 0},		{NULL, NULL}	// UP !!!
	};*/

	if ( dir == DOWN_LEFT || dir == LEFT || dir == UP_LEFT )
		return 1;
	
	return 0;
}

// starting_x moving window iterations difference to normal iterations, perspective of the receiver
inline int get_moving_iter_starting_write_x_diff(const int dir)
{
	/*
	int moving_starting_write_coord_diffs[NUM_ADJ][NUM_DIMS] = // {coord_x, coord_y}
	{
		//LEFT			CENTER		RIGHT
		{NULL, NULL},	{-1, 0},	{-1, 0},	// DOWN !!!
		{NULL, NULL},				{-1, 0},	// CENTER
		{NULL, NULL},	{-1, 0},	{-1, 0}		// UP !!!
	}; */

	if ( dir != DOWN_LEFT && dir != LEFT && dir != UP_LEFT )
		return -1;

	return 0;
}

/*	=== moving window edge proc diffs ===
// In moving window simulations
if (moving_window)
{
	// If proc is on the left edge of the simulation space
	if (proc_coords[0] == 0)
	{
		sizes[DOWN][0] += gc[0][0];
		starting_send_coord[DOWN][0] -= gc[0][0];
		starting_write_coord[DOWN][0] -= gc[0][0];

		sizes[UP][0] += gc[0][0];
		starting_send_coord[UP][0] -= gc[0][0];
		starting_write_coord[UP][0] -= gc[0][0];
	}

	// If proc is on the right edge of the simulation space
	if (proc_coords[0] == dims[0] - 1)
	{
		sizes[DOWN][0] += gc[0][1];

		sizes[UP][0] += gc[0][1];
	}
}*/

// x_size normal iterations difference to normal iterations on moving window simulations
inline int get_moving_window_edge_proc_size_x_diff(const int dir)
{	
	if ( dir == DOWN || dir == UP )
	{
		// If proc is on the left edge of the simulation space
		if (proc_coords[0] == 0)
		{
			return gc[0][0];
		}
		// If proc is on the right edge of the simulation space
		if (proc_coords[0] == dims[0] - 1)
		{
			return gc[0][1];
		}
	}

	return 0;
}

// starting send/write coord normal iterations difference on moving window simulations
inline int get_moving_window_edge_proc_coord_x_diff(const int dir)
{
	if ( dir == DOWN || dir == UP )
	{
		// If proc is on the left edge of the simulation space
		if (proc_coords[0] == 0)
		{
			return -gc[0][0];
		}
	}

	return 0;
}


/*********************************************************************************************
 
 Constructor / Destructor
 
 *********************************************************************************************/

void create_emf_segments(const int nx_local[NUM_DIMS])
{
	const int nxl0 = nx_local[0];
	const int nxl1 = nx_local[1];

	// size of the emf transmissions, in number of cells, perspective of the sender
	int sizes[NUM_ADJ][NUM_DIMS] = // {size_x, size_y}
	{
		//LEFT					CENTER				RIGHT
		{gc[0][1], gc[1][0]},	{nxl0, gc[1][0]},	{gc[0][0], gc[1][0]},	// DOWN !!!
		{gc[0][1], nxl1},							{gc[0][0], nxl1},		// CENTER
		{gc[0][1], gc[1][1]},	{nxl0, gc[1][1]},	{gc[0][0], gc[1][1]}	// UP !!!
	};

	// top left coord of the cells this proc will send to each direction, perspective of the sender
	int starting_send_coord[NUM_ADJ][NUM_DIMS] = // {coord_x, coord_y}
	{
		//LEFT			CENTER			RIGHT
		{0, nxl1 - 1},	{0, nxl1 - 1},	{nxl0 - 1, nxl1 - 1},	// DOWN !!!
		{0, 0},							{nxl0 - 1, 0},			// CENTER
		{0, 0},			{0, 0},			{nxl0 - 1, 0}			// UP !!!
	};

	// top left coord of the cells this proc will override with received cells from each direction, perspective of the receiver
	int starting_write_coord[NUM_ADJ][NUM_DIMS] = // {coord_x, coord_y}
	{
		//LEFT						CENTER				RIGHT
		{-gc[0][0], nxl1},			{0, nxl1},			{nxl0, nxl1},		// DOWN !!!
		{-gc[0][0], 0},									{nxl0, 0},			// CENTER
		{-gc[0][0], -gc[1][0]},		{0, -gc[1][0]},		{nxl0, -gc[1][0]}	// UP !!!
	};

	gaspi_pointer_t array;
	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		for (int dim = 0; dim < NUM_DIMS; dim++)
		{
			emf_seg_size[dir][dim] = sizes[dir][dim];
			emf_cell_to_send_starting_coord[dir][dim] = starting_send_coord[dir][dim];
			emf_cell_to_write_starting_coord[dir][dim] = starting_write_coord[dir][dim];
		}

		// Compute maximum segment size that will be needed for this direction
		// Segment will be used as follows:
		// Segment => [EMF_B Send zone, EMF_E Send zone, EMF_B Receive zone, EMF_E Receive zone]
		gaspi_size_t segment_size = (
			((sizes[dir][0] + gc[0][1]) * sizes[dir][1] * 2) + 								// EMF_B Send zone, EMF_E Send zone
			((sizes[OPPOSITE_DIR(dir)][0] + gc[0][1]) * sizes[OPPOSITE_DIR(dir)][1] * 2)	// EMF_B Receive zone, EMF_E Receive zone
		) * sizeof(t_vfld);

		// EMF size_x of each iteration deppends on several factors, for example:
		// If simulation is moving window or if the proc is on the left or right edges of the simulation space etc.
		// This difference is calculated on each iteration on the send_emf_gc and wait_save_emf_gc functions.
		// The difference is at most +gc[0][1], alocate enought space for this case

		SUCCESS_OR_DIE( gaspi_segment_alloc(DIR_TO_EMF_SEG_ID(dir), segment_size, GASPI_MEM_UNINITIALIZED) );
		SUCCESS_OR_DIE( gaspi_segment_register(DIR_TO_EMF_SEG_ID(dir), neighbour_rank[dir], GASPI_BLOCK) );

		SUCCESS_OR_DIE( gaspi_segment_ptr(DIR_TO_EMF_SEG_ID(dir), &array) );
		emf_segments[dir] = (t_vfld*) array;
	}
}

void emf_new(t_emf *emf, const int nx[NUM_DIMS], const int nx_local[NUM_DIMS], const t_fld box[NUM_DIMS], const float dt, const int moving_window)
{
	// Allocate arrays, initialized to 0

	size_t num_cells = (gc[0][0] + nx_local[0] + gc[0][1]) * (gc[1][0] + nx_local[1] + gc[1][1]);

	emf->E_buf = calloc(num_cells, sizeof(t_vfld));
	emf->B_buf = calloc(num_cells, sizeof(t_vfld));
	assert( emf->E_buf && emf->B_buf );

	emf->nrow = gc[0][0] + nx[0] + gc[0][1];
	emf->nrow_local = gc[0][0] + nx_local[0] + gc[0][1];

	// Make E and B point to local cell [0][0]
	emf->E = emf->E_buf + gc[0][0] + (gc[1][0] * emf->nrow_local);
	emf->B = emf->B_buf + gc[0][0] + (gc[1][0] * emf->nrow_local);


	create_emf_segments(nx_local);

	
	// store nx and gc values
	for(int i = 0; i < NUM_DIMS; i++)
	{
		emf->nx[i] = nx[i];
		emf->nx_local[i] = nx_local[i];
	}

	// store time step values
	emf->dt = dt;

	// Set cell sizes and box limits
	for(int i = 0; i < NUM_DIMS; i++)
	{
		emf->box[i] = box[i];
		emf->dx[i] = box[i] / nx[i];
		emf->box_local[i] = emf->dx[i] * emf->nx_local[i];
	}

	// Set time step
	emf->dt = dt;

	// Reset iteration number
	emf->iter = 0;

	// Set moving window information
	emf->moving_window = moving_window;
	emf->n_move = 0;
}
/*********************************************************************************************
 
 Laser Pulses
 
*********************************************************************************************/
 

t_fld gauss_phase( const t_emf_laser* const laser, const t_fld z, const t_fld r )
{
	t_fld z0   = laser -> omega0 * ( laser->W0 * laser->W0 ) / 2;
	t_fld rho2 = r*r;
	t_fld curv = rho2 * z / (z0*z0 + z*z);
	t_fld rWl2 = (z0*z0) / (z0*z0 + z*z);
	t_fld gouy_shift = atan2( z, z0 );
	
	return sqrt( sqrt(rWl2) ) * exp( -rho2 * rWl2/(laser->W0 * laser->W0) ) * cos( laser->omega0 * ( z + curv ) - gouy_shift );
}

t_fld lon_env( const t_emf_laser* const laser, const t_fld z )
{

	if ( z > laser->start ) {
		// Ahead of laser
		return 0.0;
	} else if ( z > laser->start - laser->rise ) {
		// Laser rise
		t_fld csi = z - laser->start;
		t_fld e = sin( M_PI_2 * csi / laser->rise );
		return e*e;
	} else if ( z > laser->start - (laser->rise + laser->flat) ) {
		// Flat-top
		return 1.0;
	} else if ( z > laser->start - (laser->rise + laser->flat + laser->fall) ) {
		// Laser fall
		t_fld csi = z - (laser->start - laser->rise - laser->flat - laser->fall);
		t_fld e = sin( M_PI_2 * csi / laser->fall );
		return e*e;
	}

	// Before laser
	return 0.0;
}

void div_corr_x( t_emf *emf, t_vfld* restrict E, t_vfld* restrict B)
{
	int i, j;
	
	double ex, bx;
	
	const int nrow = emf->nrow;
	const double dx_dy = emf->dx[0] / emf->dx[1];
	
	for (j = 0; j < emf->nx[1]; j++)
	{
		ex = 0.0;
		bx = 0.0;
		for (i = emf->nx[0] - 1; i >= 0; i--)
		{
			ex += dx_dy * (E[i+1 + j*nrow].y - E[i+1 + (j-1)*nrow ].y);
			E[i + j*nrow].x = ex;
			
			bx += dx_dy * (B[i + (j+1)*nrow].y - B[i + j*nrow ].y);
			B[i + j*nrow].x = bx;
		}
	}
}


void emf_add_laser( t_emf* const emf,  t_emf_laser*  laser )
{
	// Validate laser parameters
	if ( laser -> fwhm != 0 )
	{
		if ( laser -> fwhm <= 0 )
		{
			fprintf(stdout, "Invalid laser FWHM, must be > 0, aborting.\n" );
			exit(-1);
		}
		// The fwhm parameter overrides the rise/flat/fall parameters
		laser -> rise = laser -> fwhm;
		laser -> fall = laser -> fwhm;
		laser -> flat = 0.;
	}

	if ( laser -> rise <= 0 )
	{
		fprintf(stdout, "Invalid laser RISE, must be > 0, aborting.\n" );
		exit(-1);
	}

	if ( laser -> flat < 0 )
	{
		fprintf(stdout, "Invalid laser FLAT, must be >= 0, aborting.\n" );
		exit(-1);
	}

	if ( laser -> fall <= 0 )
	{
		fprintf(stdout, "Invalid laser FALL, must be > 0, aborting.\n" );
		exit(-1);
	}
	
	// Launch laser
	int i, j, nrow;
	
	t_fld r_center, z, z_2, r, r_2;
	t_fld amp, lenv, lenv_2, k;
	t_fld dx, dy;
	t_fld cos_pol, sin_pol;
	
	t_vfld* restrict E = emf -> E;
	t_vfld* restrict B = emf -> B;

	nrow = emf -> nrow;
	dx = emf -> dx[0];
	dy = emf -> dx[1];
	
	r_center = laser->axis;
	amp = laser->omega0 * laser->a0;
	
	cos_pol = cos( laser -> polarization );
	sin_pol = sin( laser -> polarization );
		
	switch (laser->type)
	{
		case PLANE:
			k = laser -> omega0;

			// i represents global x coord
			// While i is inside this procs x coord range (including guard cells)
			for (i = proc_block_low[0] - gc[0][0]; i <= proc_block_high[0] + gc[0][1] && i < emf->nx[0]; i++)
			{
				z = i * dx;
				z_2 = z + dx/2;
				
				lenv   = amp*lon_env( laser, z );
				lenv_2 = amp*lon_env( laser, z_2 );
				
				// j represents global y coord
				// While j is inside this procs y coord range (including guard cells)
				for (j = proc_block_low[1] - gc[1][0]; j <= proc_block_high[1] + gc[1][1] && j < emf->nx[1]; j++)
				{
					// Get local coords corresponding to the global coords i and j
					int x = i - proc_block_low[0];
					int y = j - proc_block_low[1];

					// E[x + y*nrow].x = 0.0
					E[x + y*nrow].y = +lenv * cos( k * z ) * cos_pol;
					E[x + y*nrow].z = +lenv * cos( k * z ) * sin_pol;

					// E[x + y*nrow].x = 0.0
					B[x + y*nrow].y = -lenv_2 * cos( k * z_2 ) * sin_pol;
					B[x + y*nrow].z = +lenv_2 * cos( k * z_2 ) * cos_pol;
				}
			}
			break;

		case GAUSSIAN:
		{
			size_t size = (gc[0][0] + emf->nx[0] + gc[0][1]) * (gc[1][0] + emf->nx[1] + gc[1][1]);

			// div_corr_x requires this proc to know other procs emf values.
			// Compute global laser emf and save it to these variables,
			// Then save relevant data on the local emf segments
			t_vfld* B_global_buff = calloc(size, sizeof(t_vfld));
			t_vfld* E_global_buff = calloc(size, sizeof(t_vfld));

			// Point to cell [0][0]
			t_vfld* restrict B_global = B_global_buff + gc[0][0] + gc[1][0] * emf->nrow;
			t_vfld* restrict E_global = E_global_buff + gc[0][0] + gc[1][0] * emf->nrow;

			// i represents global x coord
			for (i = 0; i < emf->nx[0]; i++)
			{
				z = i * dx;
				z_2 = z + dx/2;
				
				lenv   = amp*lon_env( laser, z );
				lenv_2 = amp*lon_env( laser, z_2 );
				
				// j represents to global y coord
				for (j = 0; j < emf->nx[1]; j++)
				{
					r = j * dy - r_center;
					r_2 = r + dy/2;

					// E_global[i + j*nrow].x += 0.0
					E_global[i + j*nrow].y = +lenv * gauss_phase( laser, z  , r_2 ) * cos_pol;
					E_global[i + j*nrow].z = +lenv * gauss_phase( laser, z  , r   ) * sin_pol;
					
					// B_global[i + j*nrow].x += 0.0
					B_global[i + j*nrow].y = -lenv_2 * gauss_phase( laser, z_2, r   ) * sin_pol;
					B_global[i + j*nrow].z = +lenv_2 * gauss_phase( laser, z_2, r_2 ) * cos_pol;
				}
			}

			div_corr_x(emf, E_global, B_global);

			const int nrow_local = emf->nrow_local; // Local nrow

			// j represents global y coord
			for (int j = proc_block_low[1]; j <= proc_block_high[1]; j++)
			{
				// Get local y coord corresponding to global coord j
				int y = j - proc_block_low[1];

				memcpy(&B[y * nrow_local], &B_global[proc_block_low[0] + j * nrow], proc_block_size[0] * sizeof(t_vfld));
				memcpy(&E[y * nrow_local], &E_global[proc_block_low[0] + j * nrow], proc_block_size[0] * sizeof(t_vfld));
			}
			
			free(B_global_buff);
			free(E_global_buff);
			
			break;
		}
		
		default:
			break;
	}
	
	// Set guard cell values
	emf_update_gc_gaspi(emf, 0);
	
	// printf("AFTER EMF LAZER GC UPDATE\n");
	// print_emf_e(emf);
	// print_emf_b(emf);
}

/*********************************************************************************************
 
 Diagnostics
 
 *********************************************************************************************/


void emf_report( const t_emf *emf, const char field, const char fc )
{
	int i, j;
	char vfname[3];

	// Choose field to save
	t_vfld * restrict f;
	switch (field) {
		case EFLD:
			f = emf->E;
			vfname[0] = 'E';
			break;
		case BFLD:
			f = emf->B;
			vfname[0] = 'B';
			break;
		default:
			fprintf(stdout, "Invalid field type selected, returning\n");
			return;
	}

	// Pack the information
	float * restrict const buf = malloc( emf->nx[0]*emf->nx[1]*sizeof(float) );
	float * restrict p = buf;
	switch (fc)
	{
		case 0:
			for( j = 0; j < emf->nx[1]; j++)
			{
				for ( i = 0; i < emf->nx[0]; i++ )
				{
					p[i] = f[i].x;
				}
				p += emf->nx[0];
				f += emf->nrow;
			}
			vfname[1] = '1';
			break;
		case 1:
			for( j = 0; j < emf->nx[1]; j++)
			{
				for ( i = 0; i < emf->nx[0]; i++ )
				{
					p[i] = f[i].y;
				}
				p += emf->nx[0];
				f += emf->nrow;
			}
			vfname[1] = '2';
			break;
		case 2:
			for( j = 0; j < emf->nx[1]; j++)
			{
				for ( i = 0; i < emf->nx[0]; i++ )
				{
					p[i] = f[i].z;
				}
				p += emf->nx[0];
				f += emf->nrow;
			}
			vfname[1] = '3';
			break;
		default:
			fprintf(stdout, "Invalid field component selected, returning\n");
			return;
	}
	vfname[2] = 0;

	t_zdf_grid_axis axis[2];
	axis[0] = (t_zdf_grid_axis) {
		.min = 0.0,
		.max = emf->box[0],
		.label = "x_1",
		.units = "c/\\omega_p"
	};

	axis[1] = (t_zdf_grid_axis) {
		.min = 0.0,
		.max = emf->box[1],
		.label = "x_2",
		.units = "c/\\omega_p"
	};

	t_zdf_grid_info info = {
		.ndims = 2,
		.label = vfname,
		.units = "m_e c \\omega_p e^{-1}",
		.axis = axis
	};

	info.nx[0] = emf->nx[0];
	info.nx[1] = emf->nx[1];

	t_zdf_iteration iter = {
		.n = emf->iter,
		.t = emf -> iter * emf -> dt,
		.time_units = "1/\\omega_p"
	};

	zdf_save_grid( buf, &info, &iter, "/home/bruno/zpic-out/gaspi/EMF" );
	// zdf_save_grid( buf, &info, &iter, "/home/pr1eja00/pr1eja17/zpic-out/gaspi/EMF" ); 

	// free local data
	free( buf );

}


/*********************************************************************************************
 
 Field solver
 
 *********************************************************************************************/

void yee_b( t_emf *emf, const float dt )
{
	t_vfld* const restrict B = emf -> B;
	const t_vfld* const restrict E = emf -> E;

	const t_fld dt_dx = dt / emf->dx[0];
	const t_fld dt_dy = dt / emf->dx[1];
	
	// Canonical implementation
	const int nrow = emf->nrow_local; // Local nrow

	// j = -gc[1][0]; j < nx_local + gc[1][1] - 1
	for (int j = -1; j <= emf->nx_local[1]; j++)
	{
		for (int i = -1; i <= emf->nx_local[0]; i++)
		{
			B[ i + j*nrow ].x += ( - dt_dy * ( E[i+(j+1)*nrow].z - E[i+j*nrow].z) );
			B[ i + j*nrow ].y += (   dt_dx * ( E[(i+1)+j*nrow].z - E[i+j*nrow].z) );
			B[ i + j*nrow ].z += ( - dt_dx * ( E[(i+1)+j*nrow].y - E[i+j*nrow].y) +
									 dt_dy * ( E[i+(j+1)*nrow].x - E[i+j*nrow].x) );
		}
	}
}


void yee_e( t_emf *emf, const t_current *current, const float dt )
{
	t_fld dt_dx = dt / emf->dx[0];
	t_fld dt_dy = dt / emf->dx[1];

	t_vfld* const restrict E = emf -> E;
	const t_vfld* const restrict B = emf -> B;
	const t_vfld* const restrict J = current -> J;
	
	// Canonical implementation
	const int nrow_e = emf->nrow_local;
	const int nrow_j = current->nrow_local;
	
	for (int j = 0; j <= emf->nx_local[1] + 1; j++)
	{
		for (int i = 0; i <= emf->nx_local[0] + 1; i++)
		{
			E[i+j*nrow_e].x += ( + dt_dy * ( B[i+j*nrow_e].z - B[i+(j-1)*nrow_e].z) ) - dt * J[i+j*nrow_j].x;  
			
			E[i+j*nrow_e].y += ( - dt_dx * ( B[i+j*nrow_e].z - B[(i-1)+j*nrow_e].z) ) - dt * J[i+j*nrow_j].y;  

			E[i+j*nrow_e].z += ( + dt_dx * ( B[i+j*nrow_e].y - B[(i-1)+j*nrow_e].y) - 
								   dt_dy * ( B[i+j*nrow_e].x - B[i+(j-1)*nrow_e].x) ) - dt * J[i+j*nrow_j].z;  
		}
	}
}

// OLD IMPLEMENTATION
void emf_update_gc( t_emf *emf )
{
	int i,j;
	const int nrow = emf->nrow_local; // Local nrow

	t_vfld* const restrict E = emf -> E;
	t_vfld* const restrict B = emf -> B;

	// For moving window don't update x boundaries
	if ( ! emf -> moving_window )
	{
		// x
		for (j = -gc[1][0]; j < emf->nx[1] + gc[1][1]; j++)
		{	
			// lower
			for (i = -gc[0][0]; i < 0; i++)
			{
				E[ i + j*nrow ].x = E[ emf->nx[0] + i + j*nrow ].x;
				E[ i + j*nrow ].y = E[ emf->nx[0] + i + j*nrow ].y;
				E[ i + j*nrow ].z = E[ emf->nx[0] + i + j*nrow ].z;

				B[ i + j*nrow ].x = B[ emf->nx[0] + i + j*nrow ].x;
				B[ i + j*nrow ].y = B[ emf->nx[0] + i + j*nrow ].y;
				B[ i + j*nrow ].z = B[ emf->nx[0] + i + j*nrow ].z;
			}

			// upper
			for (i = 0; i < gc[0][1]; i++)
			{
				E[ emf->nx[0] + i + j*nrow ].x = E[ i + j*nrow ].x;
				E[ emf->nx[0] + i + j*nrow ].y = E[ i + j*nrow ].y;
				E[ emf->nx[0] + i + j*nrow ].z = E[ i + j*nrow ].z;
				
				B[ emf->nx[0] + i + j*nrow ].x = B[ i + j*nrow ].x;
				B[ emf->nx[0] + i + j*nrow ].y = B[ i + j*nrow ].y;
				B[ emf->nx[0] + i + j*nrow ].z = B[ i + j*nrow ].z;
			}
		}
	}

	// y
	for (i = -gc[0][0]; i < emf->nx[0] + gc[0][1]; i++)
	{	
		// lower
		for (j = -gc[1][0]; j < 0; j++)
		{
			E[ i + j*nrow ].x = E[ i + (emf->nx[1] + j) * nrow ].x;
			E[ i + j*nrow ].y = E[ i + (emf->nx[1] + j) * nrow ].y;
			E[ i + j*nrow ].z = E[ i + (emf->nx[1] + j) * nrow ].z;
			
			B[ i + j*nrow ].x = B[ i + (emf->nx[1] + j) * nrow ].x;
			B[ i + j*nrow ].y = B[ i + (emf->nx[1] + j) * nrow ].y;
			B[ i + j*nrow ].z = B[ i + (emf->nx[1] + j) * nrow ].z;
		}
		
		// upper
		for (j = 0; j < gc[1][1]; j++)
		{
			E[ i + (emf->nx[1] + j) * nrow ].x = E[ i + j*nrow ].x;
			E[ i + (emf->nx[1] + j) * nrow ].y = E[ i + j*nrow ].y;
			E[ i + (emf->nx[1] + j) * nrow ].z = E[ i + j*nrow ].z;
			
			B[ i + (emf->nx[1] + j) * nrow ].x = B[ i + j*nrow ].x;
			B[ i + (emf->nx[1] + j) * nrow ].y = B[ i + j*nrow ].y;
			B[ i + (emf->nx[1] + j) * nrow ].z = B[ i + j*nrow ].z;
		}
	}
}

void send_emf_gc(t_emf* emf, const char moving_window_iter)
{
	const int nrow = emf->nrow_local; // Local nrow
	const int moving_window = emf->moving_window;

	const t_vfld* const restrict B = emf->B;
	const t_vfld* const restrict E = emf->E;

	// Make sure there are no uncompleted outgoing write requests
	SUCCESS_OR_DIE( gaspi_wait(Q_EMF, GASPI_BLOCK) );

	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		// For moving window simulations dont use pediodic boundaries for the left and right edge procs
		if ( !use_pediodic_boundaries(moving_window, dir) )
			continue;

		int starting_x = emf_cell_to_send_starting_coord[dir][0];
		int size_x = emf_seg_size[dir][0];

		const int starting_y = emf_cell_to_send_starting_coord[dir][1];
		const int size_y = emf_seg_size[dir][1];
		const int max_y = emf_cell_to_send_starting_coord[dir][1] + size_y;

		// On moving window simulations
		if (moving_window)
		{
			// Iterations where window will be moved
			if (moving_window_iter)
			{
				// Dont send anything to procs on the right
				if (dir == DOWN_RIGHT || dir == RIGHT || dir == UP_RIGHT)
					continue;
				
				// Apply propper size difference if window will be moved
				size_x += get_moving_iter_size_x_diff(dir);
			}
			// Iterations where the window is not moved
			else
			{
				// These only apply if proc is on the left or right simulation edges !!!				
				size_x += get_moving_window_edge_proc_size_x_diff(dir);
				starting_x += get_moving_window_edge_proc_coord_x_diff(dir);
			}
		}
		
		const size_t size_x_bytes = size_x * sizeof(t_vfld);

		// Copy EMF B data to segment
		int copy_index = 0;
		for (int y = starting_y; y < max_y; y++)
		{
			// Copy each line
			memcpy(&emf_segments[dir][copy_index], &B[starting_x + y * nrow], size_x_bytes);
			copy_index += size_x;
		}

		// Copy EMF E data to segment
		for (int y = starting_y; y < max_y; y++)
		{
			// Copy each line
			memcpy(&emf_segments[dir][copy_index], &E[starting_x + y * nrow], size_x_bytes);
			copy_index += size_x;
		}

		const int opposite_dir = OPPOSITE_DIR(dir);

		gaspi_segment_id_t local_segment = DIR_TO_EMF_SEG_ID(dir);
		gaspi_segment_id_t remote_segment = DIR_TO_EMF_SEG_ID(opposite_dir);

		gaspi_offset_t local_offset = 0; // in bytes
		gaspi_offset_t remote_offset = (emf_seg_size[opposite_dir][0] + gc[0][1]) * emf_seg_size[opposite_dir][1] * 2 * sizeof(t_vfld); // in bytes
		gaspi_size_t size = copy_index * sizeof(t_vfld); // in bytes

		// Send data
		SUCCESS_OR_DIE( gaspi_write_notify(
		local_segment,			// The segment id where data is located.
		local_offset, 			// The offset where the data is located.
		neighbour_rank[dir],	// The rank where to write and notify.
		remote_segment, 		// The remote segment id to write the data to.
		remote_offset,			// The remote offset where to write to.
		size,					// The size of the data to write.
		NOTIF_ID_EMF,			// The notification id to use.
		1,						// The notification value used.
		Q_EMF,					// The queue where to post the request.
		GASPI_BLOCK				// Timeout in milliseconds.
		));
	}
}

// Wait and save received EMF from each dir
void wait_save_emf_gc(t_emf* emf, const char moving_window_iter)
{
	const int nrow = emf->nrow_local; // Local nrow
	const int moving_window = emf->moving_window;

	t_vfld* const restrict B = emf->B;
	t_vfld* const restrict E = emf->E;

	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		// For moving window simulations dont use pediodic boundaries for the left and right edge procs
		if ( !use_pediodic_boundaries(moving_window, dir) )
			continue;

		const int opposite_dir = OPPOSITE_DIR(dir);

		int starting_x = emf_cell_to_write_starting_coord[dir][0];
		int size_x = emf_seg_size[opposite_dir][0];

		const int starting_y = emf_cell_to_write_starting_coord[dir][1];
		const int max_y = emf_cell_to_write_starting_coord[dir][1] + emf_seg_size[opposite_dir][1];

		// On moving window simulations
		if (moving_window)
		{
			// Iterations where window will be moved
			if (moving_window_iter)
			{
				// Will not receive anything from procs on the left
				if (dir == DOWN_LEFT || dir == LEFT || dir == UP_LEFT)
					continue;

				// Apply propper size and write coord differences if window will be moved
				size_x += get_moving_iter_size_x_diff(opposite_dir);
				starting_x += get_moving_iter_starting_write_x_diff(dir);
			}
			// Iterations where the window is not moved
			else
			{
				// These only apply if proc is on the left or right simulation edges !!!
				size_x += get_moving_window_edge_proc_size_x_diff(opposite_dir);
				starting_x += get_moving_window_edge_proc_coord_x_diff(dir);
			}
		}

		const size_t size_x_bytes = size_x * sizeof(t_vfld);

		// Wait for the write
		gaspi_notification_id_t id;
		SUCCESS_OR_DIE(gaspi_notify_waitsome(
			DIR_TO_EMF_SEG_ID(dir),
			NOTIF_ID_EMF,
			1,
			&id,
			GASPI_BLOCK
		));

		gaspi_notification_t value;
		SUCCESS_OR_DIE( gaspi_notify_reset(DIR_TO_EMF_SEG_ID(dir), id, &value) );

		// Copy EMF B data from segment
		int copy_index = (emf_seg_size[dir][0] + gc[0][1]) * emf_seg_size[dir][1] * 2;
		for (int y = starting_y; y < max_y; y++)
		{
			memcpy(&B[y * nrow + starting_x], &emf_segments[dir][copy_index], size_x_bytes);
			copy_index += size_x;
		}

		// Copy EMF E data from segment
		for (int y = starting_y; y < max_y; y++)
		{
			memcpy(&E[y * nrow + starting_x], &emf_segments[dir][copy_index], size_x_bytes);
			copy_index += size_x;
		}
	}
}

void emf_move_window( t_emf *emf )
{
	// printf("MOVING WINDOW NOW!\n"); fflush(stdout);

	int i, j;
	const int nrow = emf->nrow_local; // Local nrow

	const int nxl0 = emf->nx_local[0];
	const int nxl1 = emf->nx_local[1];

	t_vfld* const restrict B = emf->B;
	t_vfld* const restrict E = emf->E;

	const t_vfld zero_fld = {0.,0.,0.};

	// Shift data left 1 cell
	for (j = 0; j < nxl1; j++)
	{
		// Dont update the last row of cells, it will be overwritten later
		for (i = -gc[0][0]; i <= nxl0 - 2; i++)
		{
			B[ i + j*nrow ] = B[ i + 1 + j*nrow ];
			E[ i + j*nrow ] = E[ i + 1 + j*nrow ];
		}
	}

	// If proc is on the right edge of the simulation space
	if (proc_coords[0] == dims[0] - 1)
	{
		// Zero 3 leftmost cells on each line
		for (j = -gc[1][0]; j < nxl1 + gc[1][1]; j++)
		{
			for (i = nxl0 - 1; i < nxl0 + gc[0][1]; i++)
			{
				B[ i + j*nrow ] = zero_fld;
				E[ i + j*nrow ] = zero_fld;
			}
		}
	}
	
	// Increase moving window counter
	emf -> n_move++;
}

void emf_update_gc_gaspi(t_emf *emf, const char moving_window_iter)
{
	send_emf_gc(emf, moving_window_iter);

	// Move window if needed
	if(moving_window_iter)
		emf_move_window(emf);

	wait_save_emf_gc(emf, moving_window_iter);
}

void emf_advance( t_emf *emf, const t_current *current )
{
	uint64_t t0 = timer_ticks();
	const float dt = emf->dt;
	
	// Advance EM field using Yee algorithm modified for having E and B time centered
	// yee_b( emf, dt/2.0f ); // this is now done in sim_iter
	
	yee_e( emf, current, dt );

	yee_b( emf, dt/2.0f );
	
	// printf("BEFORE EMF GC UPDATE\n");
	// print_emf_e(emf);
	// print_emf_b(emf);

	// 1 if window will be moved on this iteration, 0 otherwise
	const char moving_window_iter = emf->moving_window && ( ((emf->iter + 1) * dt) > (emf->dx[0] * (emf->n_move + 1)) );

	// Update guard cells with new values and, if needed, move window
	emf_update_gc_gaspi(emf, moving_window_iter);

	// printf("AFTER EMF GC UPDATE\n");
	// print_emf_e(emf);
	// print_emf_b(emf);

	// Advance internal iteration number
	emf->iter++;
	
	// Update timing information
	// _emf_time += timer_interval_seconds(t0, timer_ticks());
}

void emf_get_energy( const t_emf *emf, double energy[] )
{
	int i,j;
	t_vfld* const restrict E = emf -> E;
	t_vfld* const restrict B = emf -> B;
	const int nrow = emf -> nrow;

	for(i = 0; i < 6; i++) energy[i] = 0;

	for(j = 0; i < emf -> nx[1]; j ++ ) {
		for( i = 0; i < emf -> nx[0]; i ++ ) {
			energy[0] += E[i + j*nrow].x * E[i + j*nrow].x;
			energy[1] += E[i + j*nrow].y * E[i + j*nrow].y;
			energy[2] += E[i + j*nrow].z * E[i + j*nrow].z;
			energy[3] += B[i + j*nrow].x * B[i + j*nrow].x;
			energy[4] += B[i + j*nrow].y * B[i + j*nrow].y;
			energy[5] += B[i + j*nrow].z * B[i + j*nrow].z;
		}
	}

	for( i = 0; i<6; i++) energy[i] *= 0.5 * emf -> dx[0] * emf -> dx[1];

}
