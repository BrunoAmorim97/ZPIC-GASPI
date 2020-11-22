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
extern bool is_on_edge[2];
extern gaspi_rank_t neighbour_rank[NUM_ADJ];
extern gaspi_rank_t proc_rank;
extern gaspi_rank_t num_procs;

extern int neighbour_nx[NUM_ADJ][NUM_DIMS];

extern t_vfld* current_segments[NUM_ADJ];
extern t_vfld* current_kernel_smoothing_segments[NUM_ADJ];

extern int curr_send_size[NUM_ADJ][NUM_DIMS];
extern int curr_cell_to_send_starting_coord[NUM_ADJ][NUM_DIMS];

extern int curr_kernel_size[NUM_ADJ][NUM_DIMS];
extern int curr_kernel_send_coord[NUM_ADJ][NUM_DIMS];
extern int curr_kernel_write_coord[NUM_ADJ][NUM_DIMS];

#define NUM_KERNEL_X_DIRS 2
const int kernel_x_directions[NUM_KERNEL_X_DIRS] = { LEFT, RIGHT };

#define NUM_KERNEL_Y_DIRS 6
const int kernel_y_directions[NUM_KERNEL_Y_DIRS] = { DOWN_LEFT, DOWN, DOWN_RIGHT, UP_LEFT, UP, UP_RIGHT };

extern t_current_reduce current_reduce_priv;
#pragma omp threadprivate(current_reduce_priv)

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
	printf("\n");
	fflush(stdout);

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
	printf("\n");
	fflush(stdout);

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
	printf("\n");
	fflush(stdout);
}



void alloc_private_current_reduce(t_current* current)
{
	// printf("alloc_private_current_reduce\n"); fflush(stdout);
	current_reduce_priv.J_buff = malloc(current->J_buff_size);
	current_reduce_priv.J_buff_size = current->J_buff_size;
	current_reduce_priv.J_buff_num_cells = current->J_buff_size / sizeof(t_vfld);

	// Make J point to cell [0][0]
	current_reduce_priv.J = current_reduce_priv.J_buff + gc[0][0] + (gc[1][0] * current->nrow_local);
}

t_current_reduce reset_thread_current()
{
	// printf("reset_thread_current\n"); fflush(stdout);
	// zero field
	memset(current_reduce_priv.J_buff, 0, current_reduce_priv.J_buff_size);

	return current_reduce_priv;
}

void add_thread_current(t_current_reduce current_out, t_current_reduce current_in)
{
	// printf("add_thread_current\n"); fflush(stdout);
	const unsigned int num_cells = current_in.J_buff_num_cells;

	const t_vfld* const restrict J_in = current_in.J_buff;
	t_vfld* const restrict J_out = current_out.J_buff;

	for (unsigned int i = 0; i < num_cells; i++)
	{
		// printf("J_out %f J_in %f, J_out + J_in %f\n", J_out[i].x, J_in[i].x, J_out[i].x + J_in[i].x);
		J_out[i].x += J_in[i].x;
		J_out[i].y += J_in[i].y;
		J_out[i].z += J_in[i].z;
	}
}



void curr_set_moving_window(t_current* curr)
{
	curr->moving_window = 1;
}

// create the segments that will be used to send and receive current values from other procs
void create_current_segments(const int nx_local[NUM_DIMS])
{
	const int nxl0 = nx_local[0];
	const int nxl1 = nx_local[1];

	// ===== Current transmission sizes and coords =====
	// size of the curr data transmissions, in number of cells
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

	gaspi_pointer_t pointer;
	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		for (int dim = 0; dim < NUM_DIMS; dim++)
		{
			curr_send_size[dir][dim] = sizes[dir][dim];
			curr_cell_to_send_starting_coord[dir][dim] = starting_local_coord[dir][dim];
		}

		const unsigned int curr_size_send = sizes[dir][0] * sizes[dir][1];
		const unsigned int curr_size_recv = sizes[OPPOSITE_DIR(dir)][0] * sizes[OPPOSITE_DIR(dir)][1];

		const gaspi_size_t seg_size = (curr_size_send + curr_size_recv) * sizeof(t_vfld);

		SUCCESS_OR_DIE(gaspi_segment_alloc(DIR_TO_CURR_SEG_ID(dir), seg_size, GASPI_MEM_UNINITIALIZED));
		SUCCESS_OR_DIE(gaspi_segment_register(DIR_TO_CURR_SEG_ID(dir), neighbour_rank[dir], GASPI_BLOCK));

		SUCCESS_OR_DIE(gaspi_segment_ptr(DIR_TO_CURR_SEG_ID(dir), &pointer));
		current_segments[dir] = (t_vfld*)pointer;
	}
}

// create the segments that will be used to send and receive current kernel smoothing data
void create_current_kernel_smoothing_segments(const int nx_local[NUM_DIMS], const int num_dirs, const int dirs[], const bool moving_window)
{
	const int nxl0 = nx_local[0];
	const int nxl1 = nx_local[1];

	// ===== Kernel smoothing GC sizes and coordinates =====
	// size of the curr kernel smoothing data transmissions, in number of cells, perspective of the sender
	int kernel_sizes[NUM_ADJ][NUM_DIMS] = // {size_x, size_y}
	{
		//LEFT					CENTER				RIGHT
		{gc[0][1], gc[1][0]},	{nxl0, gc[1][0]},	{gc[0][0], gc[1][0]},	// DOWN !!!
		{gc[0][1], nxl1},							{gc[0][0], nxl1}, 		// CENTER
		{gc[0][1], gc[1][1]},	{nxl0, gc[1][1]},	{gc[0][0], gc[1][1]}	// UP !!!
	};

	// top left coord of the cells this proc will send to each direction, perspective of the sender
	int kernel_starting_send_coord[NUM_ADJ][NUM_DIMS] = // {coord_x, coord_y}
	{
		//LEFT			CENTER			RIGHT
		{0, nxl1 - 1},	{0, nxl1 - 1},	{nxl0 - 1, nxl1 - 1},	// DOWN !!!
		{0, 0},							{nxl0 - 1, 0},			// CENTER
		{0, 0},			{0, 0},			{nxl0 - 1, 0}			// UP !!!
	};

	// top left coord of the cells this proc will override with received cells from each direction, perspective of the receiver
	int kernel_starting_write_coord[NUM_ADJ][NUM_DIMS] = // {coord_x, coord_y}
	{
		//LEFT							CENTER					RIGHT
		{-gc[0][0], nxl1},				{0, nxl1},				{nxl0, nxl1},		// DOWN !!!
		{-gc[0][0], 0},											{nxl0, 0},			// CENTER
		{-gc[0][0], -gc[1][0]},			{0, -gc[1][0]},			{nxl0, -gc[1][0]}	// UP !!!
	};

	// On moving window simulations
	if (moving_window)
	{
		// If proc is on the left edge of the simulation space
		if (is_on_edge[0])
		{
			kernel_sizes[DOWN][0] += gc[0][0];
			kernel_starting_send_coord[DOWN][0] -= gc[0][0];
			kernel_starting_write_coord[DOWN][0] -= gc[0][0];

			kernel_sizes[UP][0] += gc[0][0];
			kernel_starting_send_coord[UP][0] -= gc[0][0];
			kernel_starting_write_coord[UP][0] -= gc[0][0];
		}

		// If proc is on the right edge of the simulation space
		if (is_on_edge[1])
		{
			kernel_sizes[DOWN][0] += gc[0][1];

			kernel_sizes[UP][0] += gc[0][1];
		}
	}

	gaspi_pointer_t pointer;
	for (int dir_i = 0; dir_i < num_dirs; dir_i++)
	{
		const int dir = dirs[dir_i];

		for (int dim = 0; dim < NUM_DIMS; dim++)
		{
			curr_kernel_size[dir][dim] = kernel_sizes[dir][dim];
			curr_kernel_send_coord[dir][dim] = kernel_starting_send_coord[dir][dim];
			curr_kernel_write_coord[dir][dim] = kernel_starting_write_coord[dir][dim];
		}

		// Dont create segments that will not be used
		if (!can_talk_to_dir(moving_window, dir))
			continue;
		
		// The gc zones have size*2 because they have alternating zones depending if smothing pass iteration
		// number is even or odd, this is to avoid possible data overwrite by a faster neighbour
		const unsigned int kernel_size_send = kernel_sizes[dir][0] * kernel_sizes[dir][1];
		const unsigned int kernel_size_recv = kernel_sizes[OPPOSITE_DIR(dir)][0] * kernel_sizes[OPPOSITE_DIR(dir)][1];

		// Segment will be used as follows:
		// Segment => 	[Curr even Kernel Send zone, Curr even Kernel Receive zone, Curr odd Kernel Send zone, Curr odd Kernel Receive zone]
		const gaspi_size_t seg_size = (2 * (kernel_size_send + kernel_size_recv)) * sizeof(t_vfld);

		SUCCESS_OR_DIE(gaspi_segment_alloc(DIR_TO_CURR_KER_SEG_ID(dir), seg_size, GASPI_MEM_UNINITIALIZED));
		SUCCESS_OR_DIE(gaspi_segment_register(DIR_TO_CURR_KER_SEG_ID(dir), neighbour_rank[dir], GASPI_BLOCK));

		SUCCESS_OR_DIE(gaspi_segment_ptr(DIR_TO_CURR_KER_SEG_ID(dir), &pointer));
		current_kernel_smoothing_segments[dir] = (t_vfld*)pointer;
	}
}

void current_new(t_current* current, const int nx[NUM_DIMS], const int nx_local[NUM_DIMS], const t_fld box[NUM_DIMS], const float dt, const bool moving_window)
{
	// Allocate local current array, innitialized to 0

	size_t num_cells = (gc[0][0] + nx_local[0] + gc[0][1]) * (gc[1][0] + nx_local[1] + gc[1][1]);
	current->J_buff_size = num_cells * sizeof(t_vfld);

	current->J_buff = calloc(num_cells, sizeof(t_vfld));
	assert(current->J_buff);

	current->nrow = gc[0][0] + nx[0] + gc[0][1];
	current->nrow_local = gc[0][0] + nx_local[0] + gc[0][1];

	// store nx and gc values
	for (int i = 0; i < NUM_DIMS; i++)
	{
		current->nx[i] = nx[i];
		current->nx_local[i] = nx_local[i];
	}

	current->moving_window = moving_window;

	create_current_segments(nx_local);

	// Make J point to local cell [0][0]
	current->J = current->J_buff + gc[0][0] + (gc[1][0] * current->nrow_local);

	// Set cell sizes and box limits
	for (int i = 0; i < NUM_DIMS; i++)
	{
		current->box[i] = box[i];
		current->dx[i] = box[i] / nx[i];
	}

	// Clear smoothing options
	current->smooth = (t_smooth)
	{
		.xtype = NONE,
		.ytype = NONE,
		.xlevel = 0,
		.ylevel = 0
	};

	// Initialize time information
	current->iter = 0;
	current->dt = dt;
}

void curr_set_smooth(t_current* current, t_smooth* smooth)
{
	if (smooth->xtype != NONE)
	{
		if (smooth->xlevel <= 0)
		{
			printf("Invalid smooth level along x direction\n");
			exit(-1);
		}
		
		create_current_kernel_smoothing_segments(current->nx_local, NUM_KERNEL_X_DIRS, kernel_x_directions, current->moving_window);
	}

	if (smooth->ytype != NONE)
	{
		if (smooth->ylevel <= 0)
		{
			fprintf(stdout, "Invalid smooth level along y direction\n");
			exit(-1);
		}
		
		create_current_kernel_smoothing_segments(current->nx_local, NUM_KERNEL_Y_DIRS, kernel_y_directions, current->moving_window);
	}

	current->smooth = *smooth;
}

void current_zero(t_current* current)
{
	// zero field
	memset(current->J_buff, 0, current->J_buff_size);
}

// OLD IMPLEMENTATION
void current_update(t_current* current)
{
	int i, j;
	const int nrow = current->nrow_local; // Local nrow
	t_vfld* restrict const J = current->J;

	// x
	if (!current->moving_window)
	{
		// 	 j = -1			j < nx + 2
		for (j = -gc[1][0]; j < current->nx_local[1] + gc[1][1]; j++)
		{
			// lower - add the values from upper boundary ( both gc and inside box )
			//	 i = -1 		i < 2
			for (i = -gc[0][0]; i < gc[0][1]; i++)
			{
				J[i + j * nrow].x += J[current->nx_local[0] + i + j * nrow].x;
				J[i + j * nrow].y += J[current->nx_local[0] + i + j * nrow].y;
				J[i + j * nrow].z += J[current->nx_local[0] + i + j * nrow].z;
			}

			// upper - just copy the values from the lower boundary
			//	 i = -1		    i < 2
			for (i = -gc[0][0]; i < gc[0][1]; i++)
			{
				J[current->nx_local[0] + i + j * nrow] = J[i + j * nrow];
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
			J[i + j * nrow].x += J[i + (current->nx_local[1] + j) * nrow].x;
			J[i + j * nrow].y += J[i + (current->nx_local[1] + j) * nrow].y;
			J[i + j * nrow].z += J[i + (current->nx_local[1] + j) * nrow].z;
		}

		// upper - just copy the values from the lower boundary
		//	 j = -1			j < 2
		for (j = -gc[1][0]; j < gc[1][1]; j++)
		{
			J[i + (current->nx_local[1] + j) * nrow] = J[i + j * nrow];
		}
	}

	// Smoothing (NOT USED ON WEIBEL)
	current_smooth(current);

	// print_local_current(current);

	current->iter++;
}

// Send current to neighbour procs
void send_current(t_current* current)
{
	const int nrow = current->nrow_local; // Local nrow

	const t_vfld* restrict const J = current->J;

	// Make sure there are no uncompleted outgoing writes
	SUCCESS_OR_DIE(gaspi_wait(Q_CURRENT, GASPI_BLOCK));

	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		// For moving window simulations dont use pediodic boundaries for the left and right edge procs
		if (!can_talk_to_dir(current->moving_window, dir))
			continue;

		const int num_columns = curr_send_size[dir][0];
		const size_t row_size_bytes = num_columns * sizeof(t_vfld); // in bytes

		const int starting_column = curr_cell_to_send_starting_coord[dir][0];
		const int starting_row = curr_cell_to_send_starting_coord[dir][1];

		const int max_row = curr_cell_to_send_starting_coord[dir][1] + curr_send_size[dir][1];

		int copy_index = 0;
		// Copy data to segment
		for (int row = starting_row; row < max_row; row++)
		{
			memcpy(&current_segments[dir][copy_index], &J[starting_column + row * nrow], row_size_bytes);
			copy_index += num_columns;
		}

		// Sending to opposite direction segment, if I send current values to the proc to my RIGHT,
		// I want those values to be written to the remote LEFT current segment
		const gaspi_segment_id_t local_segment_id = DIR_TO_CURR_SEG_ID(dir);
		const gaspi_segment_id_t remote_segment_id = DIR_TO_CURR_SEG_ID(OPPOSITE_DIR(dir));

		const gaspi_size_t size = curr_send_size[dir][0] * curr_send_size[dir][1] * sizeof(t_vfld);			   // in bytes
		const gaspi_offset_t remote_offset = curr_send_size[dir][0] * curr_send_size[dir][1] * sizeof(t_vfld); // in bytes
		const gaspi_offset_t local_offset = 0;

		// Send data
		SUCCESS_OR_DIE(gaspi_write_notify(
			local_segment_id,	 // The segment id where data is located.
			local_offset,		 // The offset where the data is located.
			neighbour_rank[dir], // The rank where to write and notify.
			remote_segment_id,	 // The remote segment id to write the data to.
			remote_offset,		 // The remote offset where to write to.
			size,				 // The size of the data to write.
			NOTIF_ID_CURRENT,	 // The notification id to use.
			1,					 // The notification value used.
			Q_CURRENT,			 // The queue where to post the request.
			GASPI_BLOCK			 // Timeout in milliseconds.
		));
	}
}

void wait_save_update_current(t_current* current)
{
	// printf("BEFORE CURRENT GC ADD\n");
	// print_local_current(current);

	t_vfld* restrict const J = current->J;
	const int nrow = current->nrow_local; // Local nrow

	bool received_notif[NUM_ADJ];
	int num_expected_notifs = 0;

	// Check each dir for expected notifs
	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		// 1 if we are expecting a notif from dir, 0 otherwise
		const bool expecting_notif = can_talk_to_dir(current->moving_window, dir);

		received_notif[dir] = !expecting_notif;
		num_expected_notifs += expecting_notif;
	}

	int num_received_notif = 0;
	// While we have not received all the notifs we need
	while (num_received_notif != num_expected_notifs)
	{
		int dir = 0;

		// Get a dir that has an unprocessed write
		while (1)
		{
			// If we are expecting a notif from this dir
			if (!received_notif[dir])
			{
				gaspi_notification_id_t id;
				gaspi_return_t return_value;

				// Test if the notification has arrived
				SUCCESS_TIMEOUT_OR_DIE(return_value = gaspi_notify_waitsome(
					DIR_TO_CURR_SEG_ID(dir),	// The segment id
					NOTIF_ID_CURRENT,			// The notification id to wait for
					1,							// The number of notification ids this wait will accept, waiting for a specific write, so 1
					&id,						// Output parameter with the id of a received notification
					GASPI_BLOCK					// Timeout
				));

				// If this notification has arrived
				if (return_value == GASPI_SUCCESS)
				{
					received_notif[dir] = 1;
					num_received_notif++;

					// Reset notification
					gaspi_notification_t value;
					SUCCESS_OR_DIE(gaspi_notify_reset(DIR_TO_CURR_SEG_ID(dir), id, &value));

					break;
				}
			}
			// Try next dir
			dir = (dir + 1) % NUM_ADJ;
		}

		// Process received write
		// printf("Processing write from dir %d\n", dir); fflush(stdout);

		// The cells this proc will receive from dir, are the same this proc has to send to that direction
		const int starting_column = curr_cell_to_send_starting_coord[dir][0];
		const int starting_row = curr_cell_to_send_starting_coord[dir][1];

		const int max_column = curr_cell_to_send_starting_coord[dir][0] + curr_send_size[dir][0];
		const int max_row = curr_cell_to_send_starting_coord[dir][1] + curr_send_size[dir][1];

		int seg_index = curr_send_size[dir][0] * curr_send_size[dir][1];

		for (int y = starting_row; y < max_row; y++)
		{
			for (int x = starting_column; x < max_column; x++)
			{
				// printf("adding to cell x:%d y:%d, value.z:%f\n", x, y, current_segments[dir][seg_index].z); fflush(stdout);

				J[x + y * nrow].x += current_segments[dir][seg_index].x;
				J[x + y * nrow].y += current_segments[dir][seg_index].y;
				J[x + y * nrow].z += current_segments[dir][seg_index].z;

				seg_index++;
			}
		}

		// next dir
		dir = (dir + 1) % NUM_ADJ;
	}

	// printf("AFTER CURRENT GC ADD\n");
	// print_local_current(current);
}

void send_current_kernel_gc(t_current* current, const int num_dirs, const int dirs[], const int smoothing_pass_iter)
{
	const int nrow = current->nrow_local; // Local nrow
	const bool moving_window = current->moving_window;
	const t_vfld* restrict const J = current->J;

	// Make sure it is safe to modify the segment data
	SUCCESS_OR_DIE(gaspi_wait(Q_CURRENT_KERNEL, GASPI_BLOCK));

	for (int dir_i = 0; dir_i < num_dirs; dir_i++)
	{
		const int dir = dirs[dir_i];

		// For moving window simulations dont use pediodic boundaries for the left and right edges
		if (!can_talk_to_dir(moving_window, dir))
			continue;

		const int num_columns = curr_kernel_size[dir][0];
		const size_t column_size_bytes = num_columns * sizeof(t_vfld); // in bytes
		const int starting_column = curr_kernel_send_coord[dir][0];

		const int starting_row = curr_kernel_send_coord[dir][1];
		const int max_row = starting_row + curr_kernel_size[dir][1];

		const int opposite_dir = OPPOSITE_DIR(dir);

		int copy_index = 0;

		gaspi_offset_t remote_offset = curr_kernel_size[opposite_dir][0] * curr_kernel_size[opposite_dir][1];

		gaspi_notification_id_t notification_id = NOTIF_ID_CURRENT_KERNEL_EVEN;

		// If smoothing pass iteration num is odd, use alternative segment zone and corresponding notification id
		if (smoothing_pass_iter % 2 == 1)
		{
			copy_index += curr_kernel_size[dir][0] * curr_kernel_size[dir][1] +
				curr_kernel_size[opposite_dir][0] * curr_kernel_size[opposite_dir][1];

			remote_offset += curr_kernel_size[dir][0] * curr_kernel_size[dir][1] +
				curr_kernel_size[opposite_dir][0] * curr_kernel_size[opposite_dir][1];

			notification_id = NOTIF_ID_CURRENT_KERNEL_ODD;
		}

		// Local data offset is equivalent to initial copy index
		gaspi_offset_t local_offset = copy_index * sizeof(t_vfld); // in bytes
		remote_offset *= sizeof(t_vfld);

		// printf("Sending kernel gc to dir %d iter %d:\n", dir, smoothing_pass_iter);
		// printf("local_offset = %ld, size = %ld, remote_offset = %ld\n", local_offset, curr_kernel_size[dir][0] * curr_kernel_size[dir][1] * sizeof(t_vfld), remote_offset); fflush(stdout);
		// Copy data to segment
		for (int row = starting_row; row < max_row; row++)
		{
			memcpy(&current_kernel_smoothing_segments[dir][copy_index], &J[starting_column + row * nrow], column_size_bytes);

			// for (int i = copy_index; i < copy_index + num_columns ; i++)
			// {
			// 	printf("%f\n", current_kernel_smoothing_segments[dir][i].x); fflush(stdout);
			// }

			copy_index += num_columns;
		}

		gaspi_segment_id_t local_segment_id = DIR_TO_CURR_KER_SEG_ID(dir);
		gaspi_segment_id_t remote_segment_id = DIR_TO_CURR_KER_SEG_ID(opposite_dir);

		gaspi_size_t size = curr_kernel_size[dir][0] * curr_kernel_size[dir][1] * sizeof(t_vfld); // in bytes

		// Send data
		SUCCESS_OR_DIE(gaspi_write_notify(
			local_segment_id,		// The segment id where data is located.
			local_offset,			// The offset where the data is located.
			neighbour_rank[dir],	// The rank where to write and notify.
			remote_segment_id,		// The remote segment id to write the data to.
			remote_offset,			// The remote offset where to write to.
			size,					// The size of the data to write.
			notification_id,		// The notification id to use.
			1,						// The notification value used.
			Q_CURRENT_KERNEL,		// The queue where to post the request.
			GASPI_BLOCK				// Timeout in milliseconds.
		));
	}
}

void wait_save_kernel_gc(t_current* current, const int num_dirs, const int dirs[], const int smoothing_pass_iter)
{
	const int nrow = current->nrow_local; // Local nrow
	const bool moving_window = current->moving_window;

	t_vfld* restrict const J = current->J;

	for (int dir_i = 0; dir_i < num_dirs; dir_i++)
	{
		const int dir = dirs[dir_i];

		// For moving window simulations dont use pediodic boundaries for the left and right edge procs
		if (!can_talk_to_dir(moving_window, dir))
			continue;

		const int opposite_dir = OPPOSITE_DIR(dir);

		const int starting_column = curr_kernel_write_coord[dir][0];
		const int num_columns = curr_kernel_size[opposite_dir][0];

		const int starting_row = curr_kernel_write_coord[dir][1];
		const int max_row = curr_kernel_write_coord[dir][1] + curr_kernel_size[opposite_dir][1];

		int copy_index = curr_kernel_size[dir][0] * curr_kernel_size[dir][1];

		gaspi_notification_id_t notification_id = NOTIF_ID_CURRENT_KERNEL_EVEN;

		// If smoothing pass iteration number is odd, use alternative segment zone and corresponding notification id
		if (smoothing_pass_iter % 2 == 1)
		{
			copy_index += curr_kernel_size[opposite_dir][0] * curr_kernel_size[opposite_dir][1] +
				curr_kernel_size[dir][0] * curr_kernel_size[dir][1];

			notification_id = NOTIF_ID_CURRENT_KERNEL_ODD;
		}

		gaspi_notification_id_t id;
		SUCCESS_OR_DIE(gaspi_notify_waitsome(
			DIR_TO_CURR_KER_SEG_ID(dir), // The segment id
			notification_id,			 // The notification id to wait for
			1,							 // The number of notification ids this wait will accept, waiting for a specific write, so 1
			&id,						 // Output parameter with the id of a received notification
			GASPI_BLOCK					 // Timeout in milliseconds, wait until write is completed
		));

		gaspi_notification_t value;
		SUCCESS_OR_DIE(gaspi_notify_reset(DIR_TO_CURR_KER_SEG_ID(dir), id, &value));

		// printf("Received kernel gc from dir %d iter %d with value %d\n", dir, smoothing_pass_iter, value-1);
		// printf("local_offset = %ld, size = %ld\n", copy_index*sizeof(t_vfld), curr_kernel_size[opposite_dir][0] * curr_kernel_size[opposite_dir][1] * sizeof(t_vfld)); fflush(stdout);
		for (int row = starting_row; row < max_row; row++)
		{
			// for (int i = copy_index; i < copy_index + num_columns ; i++)
			// {
			// 	printf("%f\n", current_kernel_smoothing_segments[dir][i].x); fflush(stdout);
			// }

			memcpy(&J[starting_column + row * nrow], &current_kernel_smoothing_segments[dir][copy_index], num_columns * sizeof(t_vfld));
			copy_index += num_columns;
		}
	}
}

void current_report(const t_current* current, const char jc)
{
	t_vfld* f;
	float* buf, * p;
	int i, j;
	char vfname[3];

	// Pack the information
	buf = malloc(current->nx[0] * current->nx[1] * sizeof(float));
	p = buf;
	f = current->J;
	vfname[0] = 'J';

	switch (jc)
	{
	case 0:
		for (j = 0; j < current->nx[1]; j++)
		{
			for (i = 0; i < current->nx[0]; i++)
			{
				p[i] = f[i].x;
			}
			p += current->nx[0];
			f += current->nrow;
		}
		vfname[1] = '1';
		break;
	case 1:
		for (j = 0; j < current->nx[1]; j++)
		{
			for (i = 0; i < current->nx[0]; i++)
			{
				p[i] = f[i].y;
			}
			p += current->nx[0];
			f += current->nrow;
		}
		vfname[1] = '2';
		break;
	case 2:
		for (j = 0; j < current->nx[1]; j++)
		{
			for (i = 0; i < current->nx[0]; i++)
			{
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
	axis[0] = (t_zdf_grid_axis){
		.min = 0.0,
		.max = current->box[0],
		.label = "x_1",
		.units = "c/\\omega_p" };

	axis[1] = (t_zdf_grid_axis){
		.min = 0.0,
		.max = current->box[1],
		.label = "x_2",
		.units = "c/\\omega_p" };

	t_zdf_grid_info info = {
		.ndims = 2,
		.label = vfname,
		.units = "e \\omega_p^2 / c",
		.axis = axis };

	info.nx[0] = current->nx[0];
	info.nx[1] = current->nx[1];

	t_zdf_iteration iter = {
		.n = current->iter,
		.t = current->iter * current->dt,
		.time_units = "1/\\omega_p" };

	zdf_save_grid(buf, &info, &iter, "/home/bruno/zpic-out/gaspi/CURRENT");
	// zdf_save_grid( buf, &info, &iter, "/home/pr1eja00/pr1eja17/zpic-out/gaspi/CURRENT" );

	// free local data
	free(buf);
}

/*
 * get_smooth_comp
 *  Gets the value of the compensator kernel for an n pass binomial kernel
 */

void get_smooth_comp(int n, t_fld* sa, t_fld* sb)
{
	t_fld a, b, total;

	a = -1;
	b = (4.0 + 2.0 * n) / n;
	total = 2 * a + b;

	*sa = a / total;
	*sb = b / total;
}

void kernel_gc_update(t_current* current, const int num_kernel_directions, const int kernel_directions[], const int smoothing_pass_iter)
{
	send_current_kernel_gc(current, num_kernel_directions, kernel_directions, smoothing_pass_iter);

	wait_save_kernel_gc(current, num_kernel_directions, kernel_directions, smoothing_pass_iter);
}

void kernel_x(t_current* const current, const t_fld sa, const t_fld sb, const int smoothing_pass_iter)
{
	int i, j;
	t_vfld* restrict const J = current->J;
	const int nrow = current->nrow_local; // local nrow

	for (j = 0; j < current->nx_local[1]; j++)
	{
		int idx = j * nrow;

		t_vfld fl = J[idx - 1];
		t_vfld f0 = J[idx];

		for (i = 0; i < current->nx_local[0]; i++)
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
	t_vfld flbuf[current->nx[0]];
	t_vfld* restrict const J = current->J;
	const int nrow = current->nrow_local; // local nrow

	int i, j;

	// buffer lower row
	for (i = 0; i < current->nx[0]; i++)
	{
		flbuf[i] = J[i - nrow];
	}

	for (j = 0; j < current->nx[1]; j++)
	{
		int idx = j * nrow;

		for (i = 0; i < current->nx[0]; i++)
		{
			// Get lower, central and upper values
			t_vfld fl = flbuf[i];
			t_vfld f0 = J[idx + i];
			t_vfld fu = J[idx + i + nrow];

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
	kernel_gc_update(current, NUM_KERNEL_Y_DIRS, kernel_y_directions, smoothing_pass_iter);

	// printf("AFTER Y KERNEL\n");
	// print_local_current(current);
}

void current_smooth(t_current* const current)
{
	// filter kernel [sa, sb, sa]
	t_fld sa, sb;
	int i;

	// x-direction filtering
	if (current->smooth.xtype != NONE)
	{
		// binomial filter
		sa = 0.25;
		sb = 0.5;
		for (i = 0; i < current->smooth.xlevel; i++)
		{
			kernel_x(current, 0.25, 0.5, i);
		}

		// Compensator
		if (current->smooth.xtype == COMPENSATED)
		{
			get_smooth_comp(current->smooth.xlevel, &sa, &sb);
			kernel_x(current, sa, sb, i);
		}
	}

	// y-direction filtering
	if (current->smooth.ytype != NONE)
	{
		// binomial filter
		sa = 0.25;
		sb = 0.5;
		for (i = 0; i < current->smooth.xlevel; i++)
		{
			kernel_y(current, 0.25, 0.5, i);
		}

		// Compensator
		if (current->smooth.ytype == COMPENSATED)
		{
			get_smooth_comp(current->smooth.ylevel, &sa, &sb);
			kernel_y(current, sa, sb, i);
		}
	}
}
