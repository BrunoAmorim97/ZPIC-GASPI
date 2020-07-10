/*
 *  current.c
 *  zpic
 *
 *  Created by Ricardo Fonseca on 12/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#include "current.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "zdf.h"

#include <GASPI.h>

extern int dims[NUM_DIMS];
extern const int gc[NUM_DIMS][NUM_DIMS];
extern int proc_coords[NUM_DIMS];
extern gaspi_rank_t neighbour_rank[NUM_ADJ];
extern gaspi_rank_t proc_rank;
extern gaspi_rank_t num_procs;

extern int neighbour_nx[NUM_ADJ][NUM_DIMS];

extern t_vfld* current_segments[NUM_ADJ];
extern int curr_send_size[NUM_ADJ][NUM_DIMS];
extern int curr_cell_to_send_starting_coord[NUM_ADJ][NUM_DIMS];

extern int curr_kernel_sizes[NUM_ADJ][NUM_DIMS];
extern int curr_kernel_send_coord[NUM_ADJ][NUM_DIMS];
extern int curr_kernel_write_coord[NUM_ADJ][NUM_DIMS];

#define NUM_KERNEL_X_DIRS 2
const int kernel_x_directions[NUM_KERNEL_X_DIRS] = {LEFT, RIGHT};

#define NUM_KERNEL_Y_DIRS 6
const int kernel_y_directions[NUM_KERNEL_Y_DIRS] = {DOWN_LEFT, DOWN, DOWN_RIGHT, UP_LEFT, UP , UP_RIGHT};

void print_local_current(t_current* current)
{
	t_vfld* J = current->J;
	int nrow = current->nrow_local;
	printf("CURRENT X\n");
	for (int y = -gc[1][0]; y < current->nx_local[1] + gc[1][1]; y++)
	{
		if (y == 0 || y == current->nx_local[1])
		{
			printf("\n");
		}
		
		for (int x = -gc[0][0]; x < current->nx_local[0] + gc[0][1]; x++)
		{
			if (x == 0 || x == current->nx_local[0])
			{
				printf("    ");
			}
			
			printf("%f ", J[y * nrow + x].x);
		}
		printf("\n");
	}
	printf("\n"); fflush(stdout);
	
	printf("CURRENT Y\n");
	for (int y = -gc[1][0]; y < current->nx_local[1] + gc[1][1]; y++)
	{
		if (y == 0 || y == current->nx_local[1])
		{
			printf("\n");
		}
		
		for (int x = -gc[0][0]; x < current->nx_local[0] + gc[0][1]; x++)
		{
			if (x == 0 || x == current->nx_local[0])
			{
				printf("    ");
			}
			
			printf("%f ", J[y * nrow + x].y);
		}
		printf("\n");
	}
	printf("\n"); fflush(stdout);

	printf("CURRENT Z\n");
	for (int y = -gc[1][0]; y < current->nx_local[1] + gc[1][1]; y++)
	{
		if (y == 0 || y == current->nx_local[1])
		{
			printf("\n");
		}
		
		for (int x = -gc[0][0]; x < current->nx_local[0] + gc[0][1]; x++)
		{
			if (x == 0 || x == current->nx_local[0])
			{
				printf("    ");
			}
			
			printf("%f ", J[y * nrow + x].z);
		}
		printf("\n");
	}
	printf("\n"); fflush(stdout);
}

void curr_set_moving_window(t_current* curr)
{
	curr->moving_window = 1;
}

// create the segments that will be used to send and receive current values from other procs
void create_current_segments(const int nx_local[NUM_DIMS], const int moving_window)
{
	const int nxl0 = nx_local[0];
	const int nxl1 = nx_local[1];

	// ===== Current transmission sizes and coords =====
	const int sizes[NUM_ADJ][NUM_DIMS] = // {size_x, size_y}
	{
		//LEFT												CENTER												RIGHT
		{gc[0][0] + gc[0][1], gc[1][0] + gc[1][1]},			{gc[0][0] + nxl0 + gc[0][1], gc[1][0] + gc[1][1]},	{gc[0][0] + gc[0][1], gc[1][0] + gc[1][1]},			// DOWN !!!
		{gc[0][0] + gc[0][1], gc[1][0] + nxl1 + gc[1][1]},														{gc[0][0] + gc[0][1], gc[1][0] + nxl1 + gc[1][1]},	// CENTER
		{gc[0][0] + gc[0][1], gc[1][0] + gc[1][1]},			{gc[0][0] + nxl0 + gc[0][1], gc[1][0] + gc[1][1]},	{gc[0][0] + gc[0][1], gc[1][0] + gc[1][1]}			// UP !!!
	};

	// top left coord of the cells this proc will send to each direction
	const int starting_local_coord[NUM_ADJ][NUM_DIMS] = // {cord_x, coord_y}
	{
		//LEFT					CENTER					RIGHT
		{-gc[0][0], nxl1 - 1},	{-gc[0][0], nxl1 - 1},	{nxl0 - 1, nxl1 - 1},	// DOWN !!!
		{-gc[0][0], -gc[1][0]}, 						{nxl0 - 1, -gc[1][0]},	// CENTER
		{-gc[0][0], -gc[1][0]},	{-gc[0][0], -gc[1][0]},	{nxl0 - 1, -gc[1][0]}	// UP !!!
	};


	// ===== Kernel GC sizes and coordinates =====
	int kernel_sizes[NUM_ADJ][NUM_DIMS] = // {size_x, size_y}
	{
		//LEFT					CENTER				RIGHT
		{gc[0][0], gc[1][1]},	{nxl0, gc[1][1]},	{gc[0][1], gc[1][1]},	// DOWN !!!
		{gc[0][0], nxl1},							{gc[0][1], nxl1},		// CENTER
		{gc[0][0], gc[1][0]},	{nxl0, gc[1][0]},	{gc[0][1], gc[1][0]}	// UP !!!
	};

	// top left coord of the cells this proc will send to each direction
	int kernel_starting_send_coord[NUM_ADJ][NUM_DIMS] = // {coord_x, coord_y}
	{
		//LEFT			CENTER			RIGHT
		{0, nxl1 - 1},	{0, nxl1 - 1},	{nxl0 - 1, nxl1 - 1},	// DOWN !!!
		{0, 0},							{nxl0 - 1, 0},			// CENTER
		{0, 0},			{0, 0},			{nxl0 - 1, 0}			// UP !!!
	};

	// top left coord of the cells this proc will override with received cells from each direction
	int kernel_starting_write_coord[NUM_ADJ][NUM_DIMS] = // {coord_x, coord_y}
	{
		//LEFT							CENTER					RIGHT
		{-gc[0][0], nxl1},				{0, nxl1},				{nxl0, nxl1},		// DOWN !!!
		{-gc[0][0], 0},											{nxl0, 0},			// CENTER
		{-gc[0][0], -gc[1][0]},			{0, -gc[1][0]},			{nxl0, -gc[1][0]}	// UP !!!
	};

	// In moving window simulations
	if (moving_window)
	{
		// If proc is on the left edge of the simulation space
		if (proc_coords[0] == 0)
		{
			kernel_sizes[DOWN][0] += gc[0][0];
			kernel_starting_send_coord[DOWN][0] -= gc[0][0];
			kernel_starting_write_coord[DOWN][0] -= gc[0][0];

			kernel_sizes[UP][0] += gc[0][0];
			kernel_starting_send_coord[UP][0] -= gc[0][0];
			kernel_starting_write_coord[UP][0] -= gc[0][0];
		}

		// If proc is on the right edge of the simulation space
		if (proc_coords[0] == dims[0] - 1)
		{
			kernel_sizes[DOWN][0] += gc[0][1];

			kernel_sizes[UP][0] += gc[0][1];
		}
	}
	
	gaspi_pointer_t array;
	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		gaspi_size_t size = sizeof(t_vfld);

		// Kernel gc communication will use the same segments as normal current cell update.
		// Extend segment size so that both communications have their own dedicated segment zone.
		// The gc zones have kernel_size*2 beacause they have alternating zones depending if smothing pass iteration
		// number is even or odd, this is to avoid possible data overwrite by a faster neighbour 
		gaspi_size_t size_kernel = sizeof(t_vfld);

		for (int dim = 0; dim < NUM_DIMS; dim++)
		{
			curr_send_size[dir][dim] = sizes[dir][dim];
			curr_cell_to_send_starting_coord[dir][dim] = starting_local_coord[dir][dim];

			curr_kernel_sizes[dir][dim] = kernel_sizes[dir][dim];
			curr_kernel_send_coord[dir][dim] = kernel_starting_send_coord[dir][dim];
			curr_kernel_write_coord[dir][dim] = kernel_starting_write_coord[dir][dim];

			size *= sizes[dir][dim];
			size_kernel *= kernel_sizes[dir][dim];
		}

		// Segment will be used as follows:
		// Segment => 	[Curr Send zone, Curr receive zone,
		// 				Curr even Kernel Send zone, Curr even Kernel Receive zone, Curr odd Kernel Send zone, Curr odd Kernel Receive zone]
		SUCCESS_OR_DIE( gaspi_segment_alloc(DIR_TO_CURR_SEG_ID(dir), (size  * 2) + (size_kernel * 4), GASPI_MEM_UNINITIALIZED) );
		SUCCESS_OR_DIE( gaspi_segment_register(DIR_TO_CURR_SEG_ID(dir), neighbour_rank[dir], GASPI_BLOCK) );

		SUCCESS_OR_DIE( gaspi_segment_ptr(DIR_TO_CURR_SEG_ID(dir), &array) );
		current_segments[dir] = (t_vfld*) array;
	}
}

//send current to neighbour procs
void send_current(t_current* current)
{
	const int nrow = current->nrow_local; // Local nrow
	const t_vfld* restrict const J = current -> J;

	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		// For moving window simulations dont use pediodic boundaries for the left and right edge procs
		if ( !can_send_gc(current->moving_window, dir) )
			continue;

		// printf("Sending current to dir %d, proc %d\n", dir, neighbour_rank[dir]); fflush(stdout);

		const int size_x = curr_send_size[dir][0];
		const size_t size_x_bytes = size_x * sizeof(t_vfld); // in bytes
		const int starting_x = curr_cell_to_send_starting_coord[dir][0];

		const int starting_y = curr_cell_to_send_starting_coord[dir][1];
		const int max_y = curr_cell_to_send_starting_coord[dir][1] + curr_send_size[dir][1];

		int copy_index = 0;
		// Copy data to segment
		for (int y = starting_y; y < max_y; y++)
		{
			memcpy(&current_segments[dir][copy_index], &J[starting_x + y * nrow], size_x_bytes);
			copy_index += size_x;
		}
		
		// Sending to opposite direction segment, if I send current values to the proc to my RIGHT,
		// I want those values to be written to the remote LEFT current segment
		gaspi_segment_id_t local_segment_id = DIR_TO_CURR_SEG_ID(dir);
		gaspi_segment_id_t remote_segment_id = DIR_TO_CURR_SEG_ID(OPPOSITE_DIR(dir));

		gaspi_size_t size = curr_send_size[dir][0] * curr_send_size[dir][1] * sizeof(t_vfld); // in bytes
		gaspi_offset_t local_offset = 0;
		gaspi_offset_t remote_offset = size; // in bytes

		// Send data
		SUCCESS_OR_DIE( gaspi_write_notify(
		local_segment_id,		// The segment id where data is located.
		local_offset, 			// The offset where the data is located.
		neighbour_rank[dir],	// The rank where to write and notify.
		remote_segment_id,		// The remote segment id to write the data to.
		remote_offset,			// The remote offset where to write to.
		size,					// The size of the data to write.
		NOTIF_ID_CURRENT,		// The notification id to use.
		1,						// The notification value used.
		Q_CURRENT,				// The queue where to post the request.
		GASPI_BLOCK				// Timeout in milliseconds.
		));
	}
}

// Also applies current smoothing if necessary
void wait_save_update_current(t_current* current)
{
	t_vfld* restrict const J = current->J;
	const int nrow = current->nrow_local; // Local nrow

	// printf("BEFORE CURRENT GC ADD\n"); fflush(stdout);
	// print_local_current(current);

	// Make sure there are no uncompleted outgoing writes
	SUCCESS_OR_DIE( gaspi_wait(Q_CURRENT, GASPI_BLOCK) );

	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		// For moving window simulations dont use pediodic boundaries for the left and right edge procs
		if ( !can_send_gc(current->moving_window, dir) )
			continue;

		// The cells this proc will receive from dir, are the same this proc has to send to that direction
		const int starting_x = curr_cell_to_send_starting_coord[dir][0];
		const int starting_y = curr_cell_to_send_starting_coord[dir][1];

		const int max_x = curr_cell_to_send_starting_coord[dir][0] + curr_send_size[dir][0];
		const int max_y = curr_cell_to_send_starting_coord[dir][1] + curr_send_size[dir][1];

		// printf("\nFrom dir %d\n", dir); fflush(stdout);

		int index = curr_send_size[dir][0] * curr_send_size[dir][1];

		// Wait for the write
		gaspi_notification_id_t id;
		SUCCESS_OR_DIE( gaspi_notify_waitsome(
		DIR_TO_CURR_SEG_ID(dir),	// The segment id
		NOTIF_ID_CURRENT,			// The notification id to wait for
		1,							// The number of notification ids this wait will accept, waiting for a specific write, so 1
		&id,						// Output parameter with the id of a received notification
		GASPI_BLOCK					// Timeout in milliseconds, wait until write is completed
		));

		gaspi_notification_t value;
		SUCCESS_OR_DIE( gaspi_notify_reset(DIR_TO_CURR_SEG_ID(dir), id, &value) );

		for (int y = starting_y; y < max_y; y++)
		{
			for (int x = starting_x; x < max_x; x++)
			{
				// printf("adding to cell x:%d y:%d, value.z:%f\n", x, y, current_segments[dir][index].z); fflush(stdout);

				J[x + y*nrow].x += current_segments[dir][index].x;
				J[x + y*nrow].y += current_segments[dir][index].y;
				J[x + y*nrow].z += current_segments[dir][index].z;

				index++;
			}
		}
	}

	// printf("AFTER CURRENT GC ADD\n");
	// print_local_current(current);

	// Smoothing
	current_smooth( current );

	// printf("AFTER SMOOTHING\n"); fflush(stdout);
	// print_local_current(current);

	current -> iter++;
}

void current_new(t_current *current, const int nx[NUM_DIMS], const int nx_local[NUM_DIMS], const t_fld box[NUM_DIMS], const float dt, const int moving_window)
{
	// Allocate local current array, innitialized to 0

	size_t num_cells = (gc[0][0] + nx_local[0] + gc[0][1]) * (gc[1][0] + nx_local[1] + gc[1][1]);
	current->J_size = num_cells * sizeof(t_vfld);

	current->J_buf = calloc(num_cells, sizeof(t_vfld));
	assert( current->J_buf );

	current -> nrow = gc[0][0] + nx[0] + gc[0][1];
	current -> nrow_local = gc[0][0] + nx_local[0] + gc[0][1];


	create_current_segments(nx_local, moving_window);


	// store nx and gc values
	for(int i = 0; i < NUM_DIMS; i++)
	{
		current->nx[i] = nx[i];
		current->nx_local[i] = nx_local[i];
	}
	
	// Make J point to local cell [0][0]
	current->buff_offset = gc[0][0] + (gc[1][0] * current->nrow_local); //offset not in bytes
	current->J = current->J_buf + current->buff_offset;
	
	// Set cell sizes and box limits
	for(int i = 0; i < NUM_DIMS; i++)
	{
		current -> box[i] = box[i];
		current -> dx[i] = box[i] / nx[i];
		current -> box_local[i] = current->dx[i] * current->nx_local[i];
	}

	// Clear smoothing options
	current -> smooth = (t_smooth) {
		.xtype = NONE,
		.ytype = NONE,
		.xlevel = 0,
		.ylevel = 0
	};

	// Initialize time information
	current -> iter = 0;
	current -> dt = dt;

	current -> moving_window = moving_window;
}

void current_zero( t_current *current )
{
	// zero field
	memset( current->J_buf, 0, current->J_size );
}

// OLD IMPLEMENTATION
void current_update( t_current *current )
{
	int i,j;
	const int nrow = current->nrow_local; // Local nrow
	t_vfld* restrict const J = current -> J;
	
	// x
	if ( ! current -> moving_window )
	{
		// 	 j = -1			j < nx + 2
		for (j = -gc[1][0]; j < current->nx_local[1] + gc[1][1]; j++)
		{
			// lower - add the values from upper boundary ( both gc and inside box )
			//	 i = -1 		i < 2
			for (i = -gc[0][0]; i < gc[0][1]; i++)
			{
				J[ i + j*nrow ].x += J[ current->nx_local[0] + i + j*nrow ].x;
				J[ i + j*nrow ].y += J[ current->nx_local[0] + i + j*nrow ].y;
				J[ i + j*nrow ].z += J[ current->nx_local[0] + i + j*nrow ].z;
			}
			
			// upper - just copy the values from the lower boundary
			//	 i = -1		    i < 2
			for (i = -gc[0][0]; i < gc[0][1]; i++)
			{
				J[ current->nx_local[0] + i + j*nrow ] = J[ i + j*nrow ];
			}
		}
	}
	
	// y
	//	 i = -1		   i < nx_local + 2
	for (i = -gc[0][0]; i < current->nx_local[0] + gc[0][1]; i++)
	{
		// lower - add the values from upper boundary ( both gc and inside box )
		//	 j = -1			j < 2
		for (j = -gc[1][0]; j < gc[1][1]; j++)
		{
			J[ i + j*nrow ].x += J[ i + (current->nx_local[1] + j) * nrow ].x;
			J[ i + j*nrow ].y += J[ i + (current->nx_local[1] + j) * nrow ].y;
			J[ i + j*nrow ].z += J[ i + (current->nx_local[1] + j) * nrow ].z;
		}
		
		// upper - just copy the values from the lower boundary
		//	 j = -1			j < 2
		for (j = -gc[1][0]; j < gc[1][1]; j++)
		{
			J[ i + (current->nx_local[1] + j) * nrow ] = J[ i + j * nrow ];
		}
	}

	// Smoothing (NOT USED ON WEIBEL)
	current_smooth( current );

	// print_local_current(current);

	current -> iter++;
}

void current_report( const t_current *current, const char jc )
{
	t_vfld *f;
	float *buf, *p;
	int i, j;
	char vfname[3];
		
	// Pack the information
	buf = malloc( current->nx[0]*current->nx[1]*sizeof(float) );
	p = buf;
	f = current->J;
	vfname[0] = 'J';
	
	switch (jc) {
		case 0:
			for( j = 0; j < current->nx[1]; j++) {
				for ( i = 0; i < current->nx[0]; i++ ) {
					p[i] = f[i].x;
				}
				p += current->nx[0];
				f += current->nrow;
			}
			vfname[1] = '1';
			break;
		case 1:
			for( j = 0; j < current->nx[1]; j++) {
				for ( i = 0; i < current->nx[0]; i++ ) {
					p[i] = f[i].y;
				}
				p += current->nx[0];
				f += current->nrow;
			}
			vfname[1] = '2';
			break;
		case 2:
			for( j = 0; j < current->nx[1]; j++) {
				for ( i = 0; i < current->nx[0]; i++ ) {
					p[i] = f[i].z;
				}
				p += current->nx[0];
				f += current->nrow;
			}
			vfname[1] = '3';
			break;
	}
	vfname[2] = 0;
	
	t_zdf_grid_axis axis[2];
	axis[0] = (t_zdf_grid_axis) {
		.min = 0.0,
		.max = current->box[0],
		.label = "x_1",
		.units = "c/\\omega_p"
	};

	axis[1] = (t_zdf_grid_axis) {
		.min = 0.0,
		.max = current->box[1],
		.label = "x_2",
		.units = "c/\\omega_p"
	};

	t_zdf_grid_info info = {
		.ndims = 2,
		.label = vfname,
		.units = "e \\omega_p^2 / c",
		.axis = axis
	};

	info.nx[0] = current->nx[0];
	info.nx[1] = current->nx[1];

	t_zdf_iteration iter = {
		.n = current->iter,
		.t = current -> iter * current -> dt,
		.time_units = "1/\\omega_p"
	};

	zdf_save_grid( buf, &info, &iter, "/home/bruno/zpic-out/gaspi/CURRENT" );
	// zdf_save_grid( buf, &info, &iter, "/home/pr1eja00/pr1eja17/zpic-out/gaspi/CURRENT-gaspi" );  
	
	// free local data
	free( buf );
}

/*
 * get_smooth_comp
 *  Gets the value of the compensator kernel for an n pass binomial kernel
 */

void get_smooth_comp( int n, t_fld* sa, t_fld* sb) {
	t_fld a,b,total;

	a = -1;
	b = (4.0 + 2.0*n)/n;
	total = 2*a + b;

	*sa = a / total;
	*sb = b / total;
}

void send_current_kernel_gc(t_current* current, const int num_dirs, const int dirs[], const int smoothing_pass_num)
{
	const int nrow = current->nrow_local; // Local nrow
	const int buff_offset = current->buff_offset; //offset not in bytes
	const int moving_window = current->moving_window;

	for (int dir_i = 0; dir_i < num_dirs; dir_i++)
	{
		const int dir = dirs[dir_i];

		// For moving window simulations dont use pediodic boundaries for the left and right edges
		if ( !can_send_gc(moving_window, dir) )
			continue;

		const int opposite_dir = OPPOSITE_DIR(dir);

		const int size_x = curr_kernel_sizes[opposite_dir][0];
		const int size_y = curr_kernel_sizes[opposite_dir][1];

		const int starting_x = curr_kernel_send_coord[dir][0];
		const int starting_y = curr_kernel_send_coord[dir][1];

		gaspi_segment_id_t remote_segment_id = DIR_TO_CURR_SEG_ID(opposite_dir);

		gaspi_segment_id_t local_segment_ids[size_y]; memset(local_segment_ids, CURRENT, size_y * sizeof(gaspi_segment_id_t));
		gaspi_segment_id_t remote_segment_ids[size_y]; memset(remote_segment_ids, remote_segment_id, size_y * sizeof(gaspi_segment_id_t));

		gaspi_offset_t local_offsets[size_y]; // in bytes
		local_offsets[0] = (buff_offset + starting_x + (starting_y * nrow)) * sizeof(t_vfld);
		for ( int i = 1; i < size_y; i++)
		{
			local_offsets[i] = local_offsets[i - 1] + (nrow * sizeof(t_vfld));
		}
		

		gaspi_offset_t remote_offsets[size_y]; // in bytes
		remote_offsets[0] = curr_send_size[dir][0] * curr_send_size[dir][1] * sizeof(t_vfld);

		// To avoid overwriting gc data of a slower neighbour, if smoothing pass iteration num is odd use alternate segment zone. 
		if (smoothing_pass_num % 2 == 1)
		{
			remote_offsets[0] += curr_kernel_send_coord[opposite_dir][0] * curr_kernel_send_coord[opposite_dir][1] * sizeof(t_vfld);
		}

		for ( int i = 1; i < size_y; i++)
		{
			remote_offsets[i] = remote_offsets[i - 1] + (size_x * sizeof(t_vfld));
		}


		// in bytes
		gaspi_size_t write_sizes[size_y]; for (int i = 0; i < size_y; i++) write_sizes[i] = size_x * sizeof(t_vfld);

		SUCCESS_OR_DIE(gaspi_write_list_notify(
			size_y,
			local_segment_ids,
			local_offsets,
			neighbour_rank[dir],
			remote_segment_ids,
			remote_offsets,
			write_sizes,
			remote_segment_id,
			NOTIF_ID_CURRENT_KERNEL,
			1,
			Q_CURRENT,
			GASPI_BLOCK
		));
	}
}

void wait_save_kernel_gc(t_current* current, const int num_dirs, const int dirs[], const int smoothing_pass_num)
{
	const int nrow = current->nrow_local; // Local nrow
	const int moving_window = current->moving_window;

	t_vfld* restrict const J = current -> J;

	// Make sure there are no uncompleted outgoing writes
	SUCCESS_OR_DIE(gaspi_wait(Q_CURRENT, GASPI_BLOCK));

	for (int dir_i = 0; dir_i < num_dirs; dir_i++)
	{
		const int dir = dirs[dir_i];

		// For moving window simulations dont use pediodic boundaries for the left and right edge procs
		if ( !can_send_gc(moving_window, dir) )
			continue;

		gaspi_notification_id_t id;
		SUCCESS_OR_DIE( gaspi_notify_waitsome(
		DIR_TO_CURR_SEG_ID(dir),	// The segment id
		NOTIF_ID_CURRENT_KERNEL,	// The notification id to wait for
		1,							// The number of notification ids this wait will accept, waiting for a specific write, so 1
		&id,						// Output parameter with the id of a received notification
		GASPI_BLOCK					// Timeout in milliseconds, wait until write is completed
		));

		gaspi_notification_t value;
		SUCCESS_OR_DIE( gaspi_notify_reset(DIR_TO_CURR_SEG_ID(dir), id, &value) );

		const int starting_x = curr_kernel_write_coord[dir][0];
		const int starting_y = curr_kernel_write_coord[dir][1];

		const int size_x = curr_kernel_sizes[dir][0];

		const int max_y = curr_kernel_write_coord[dir][1] + curr_kernel_sizes[dir][1];

		int copy_index = curr_send_size[dir][0] * curr_send_size[dir][1];

		// iteration number is odd, use alternative segment zone
		if (smoothing_pass_num % 2 == 1)
		{
			copy_index += curr_kernel_sizes[dir][0] * curr_kernel_sizes[dir][1];
		}

		for (int y = starting_y; y < max_y; y++)
		{
			memcpy(&J[starting_x + y * nrow], &current_segments[dir][copy_index], size_x * sizeof(t_vfld));
			copy_index += size_x;
		}
	}
}

void kernel_gc_update(t_current* current, const int num_kernel_directions, const int kernel_directions[], const int smoothing_pass_iter)
{
	send_current_kernel_gc(current, num_kernel_directions, kernel_directions, smoothing_pass_iter);

	wait_save_kernel_gc(current, num_kernel_directions, kernel_directions, smoothing_pass_iter);
}

void kernel_x( t_current* const current, const t_fld sa, const t_fld sb, const int smoothing_pass_iter)
{
	int i, j;
	t_vfld* restrict const J = current->J;
	const int nrow = current->nrow_local; // local nrow

	for( j = 0; j < current->nx_local[1]; j++)
	{	
		int idx = j * nrow;

		t_vfld fl = J[idx-1];
		t_vfld f0 = J[idx  ];

		for( i = 0; i < current->nx_local[0]; i++)
		{
			t_vfld fu = J[idx + i + 1];

			t_vfld fs;

			fs.x = sa * fl.x + sb * f0.x + sa * fu.x;
			fs.y = sa * fl.y + sb * f0.y + sa * fu.y;
			fs.z = sa * fl.z + sb * f0.z + sa * fu.z;

			J[idx + i] = fs;

			fl = f0;
			f0 = fu;
		}

	}

	// Update x boundaries
	kernel_gc_update(current, NUM_KERNEL_X_DIRS, kernel_x_directions, smoothing_pass_iter);

	// printf("AFTER X KERNEL\n");
	// print_local_current(current);
}

void kernel_y(t_current* const current, const t_fld sa, const t_fld sb, const int smoothing_pass_iter)
{
	t_vfld flbuf[ current->nx[0] ];
	t_vfld* restrict const J = current->J;
	const int nrow = current->nrow;

	int i, j;

	// buffer lower row
	for( i = 0; i < current -> nx[0]; i++)
	{
		flbuf[i] = J[i - nrow];
	}


	for( j = 0; j < current -> nx[1]; j++)
	{
		int idx = j * nrow;

		for( i = 0; i < current -> nx[0]; i++)
		{
			// Get lower, central and upper values
			t_vfld fl = flbuf[i];
			t_vfld f0 = J[ idx + i ];
			t_vfld fu = J[ idx + i + nrow ];

			// Store the value that will be overritten for use in the next row
			flbuf[i] = f0;

			// Convolution with kernel
			t_vfld fs;
			fs.x = sa * fl.x + sb * f0.x + sa * fu.x;
			fs.y = sa * fl.y + sb * f0.y + sa * fu.y;
			fs.z = sa * fl.z + sb * f0.z + sa * fu.z;

			// Store result
			J[idx + i] = fs;
		}
	}

	// Update y boundaries
	kernel_gc_update( current, NUM_KERNEL_Y_DIRS, kernel_y_directions, smoothing_pass_iter);

	// printf("AFTER Y KERNEL\n");
	// print_local_current(current);
}


void current_smooth( t_current* const current )
{
	// filter kernel [sa, sb, sa]
	t_fld sa, sb;

	int i;

	// x-direction filtering
	if ( current -> smooth.xtype != NONE )
	{
		// binomial filter
		sa = 0.25; sb = 0.5;
		for( i = 0; i < current -> smooth.xlevel; i++)
		{
			kernel_x( current, 0.25, 0.5, i);
		}

		// Compensator
		if ( current -> smooth.xtype == COMPENSATED )
		{
			get_smooth_comp( current -> smooth.xlevel, &sa, &sb );
			kernel_x( current, sa, sb, i);
		}
	}

	// y-direction filtering
	if ( current -> smooth.ytype != NONE )
	{
		// binomial filter
		sa = 0.25; sb = 0.5;
		for( i = 0; i < current -> smooth.xlevel; i++)
		{
			kernel_y( current, 0.25, 0.5, i);
		}

		// Compensator
		if ( current -> smooth.ytype == COMPENSATED )
		{
			get_smooth_comp( current -> smooth.ylevel, &sa, &sb );
			kernel_y( current, sa, sb, i);
		}
	}
}

