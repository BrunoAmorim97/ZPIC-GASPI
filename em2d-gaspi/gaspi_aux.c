#include "gaspi_aux.h"

extern int proc_block_low[NUM_DIMS];
extern int proc_block_high[NUM_DIMS];
extern int proc_block_size[NUM_DIMS];

extern int proc_coords[NUM_DIMS];
extern int dims[NUM_DIMS];

extern gaspi_rank_t neighbour_rank[NUM_ADJ];
extern unsigned int neighbour_nx[NUM_ADJ][NUM_DIMS];

int cart_rank(int coords[NUM_DIMS])
{
	return coords[0] + (coords[1] * dims[0]);
}

void cart_coords(int rank, int coords[NUM_DIMS])
{
	coords[0] = rank % dims[0];
	coords[1] = rank / dims[0];
}

void create_dims(int num_procs)
{
	int num_factors;
	int* factors = get_factors(num_procs, &num_factors);

	assign_factors(factors, num_factors, dims);
	free(factors);
}


int* get_factors(int num, int* num_factors)
{
	if(num  < 2)
	{
		(*num_factors) = 0;
		return NULL;
	}

	int num_sqrt = ceil(sqrt(num));
	int size = ceil(log2(num));
	int* factors = (int *) malloc((unsigned) size * sizeof(int));

	int i = 0;

	// occurances of factor 2
	while((num % 2) == 0)
	{
		num /= 2;
		factors[i++] = 2;
	}

	int d;
	// occurances of uneven primes up to sqrt(num)
	for(d = 3; (num > 1) && (d <= num_sqrt); d += 2)
	{
		while((num % d) == 0)
		{
			num /= d;
			factors[i++] = d;
		}
	}

	//add last factor
	if(num != 1)
	{
		factors[i++] = num;
	}

	(*num_factors) = i;

	return factors;
}

void assign_factors(int* factors, int num_factors, int dims[NUM_DIMS])
{
	dims[0] = 1;
	dims[1] = 1;

	// assign factors from highest to lowest
	for (int i = num_factors - 1; i >= 0; i--)
	{
		// assign to dimention with the lowest number of divisions, prioritize divisions on the y axis
		if (dims[0] <= dims[1])
			dims[0] *= factors[i];
		
		else
			dims[1] *= factors[i];	
	}
}

void assign_proc_blocks(const int nx[NUM_DIMS])
{
	static char blocks_set = 0;

	if (blocks_set == 0)
	{
		//only need to assign the blocks once
		blocks_set = 1;


		proc_block_low[0] = BLOCK_LOW(proc_coords[0], dims[0], nx[0]);
		proc_block_low[1] = BLOCK_LOW(proc_coords[1], dims[1], nx[1]);

		proc_block_high[0] = BLOCK_HIGH(proc_coords[0], dims[0], nx[0]);
		proc_block_high[1] = BLOCK_HIGH(proc_coords[1], dims[1], nx[1]);

		proc_block_size[0] = BLOCK_SIZE(proc_coords[0], dims[0], nx[0]);
		proc_block_size[1] =  BLOCK_SIZE(proc_coords[1], dims[1], nx[1]);

		// printf("proc low x:%d y:%d\n", proc_block_low[0], proc_block_low[1]);
		// printf("proc high x:%d y:%d\n", proc_block_high[0], proc_block_high[1]);
		// printf("proc size x:%d y:%d\n", proc_block_size[0], proc_block_size[1]); fflush(stdout);
	}
}

// periodic boundary enforcer. If a coord leaves its allowed range, enforce a periodic boundary
inline int periodic_coord(int coord, int coord_max)
{
	if(coord < 0)
		return coord + coord_max;

	if(coord >= coord_max)
		return coord - coord_max;

	return coord;
}

// populate the neighbour_rank array with correct rank in that direction and save their simulation space size
void discover_neighbours(int proc_coords[NUM_DIMS], int dims[NUM_DIMS], int nx[NUM_DIMS])
{
	int adj_i = 0;

	for (int y = -1; y <= 1; y++)
	{
		for (int x = -1; x <= 1; x++)
		{
			if(y == 0 && x == 0)
				continue;

			int new_x = periodic_coord(x + proc_coords[0], dims[0]);
			int new_y = periodic_coord(y + proc_coords[1], dims[1]);

			int coords[NUM_DIMS] = {new_x, new_y};

			for (int dim = 0; dim < NUM_DIMS; dim++)
			{
				neighbour_nx[adj_i][dim] = BLOCK_SIZE(coords[dim], dims[dim], nx[dim]);
			}

			neighbour_rank[adj_i++] = cart_rank(coords);
			//printf("rank x:%d y:%d = %d\n", nx, ny, cart_rank(coords));
		}
	}

/* 	for (int dir = 0; dir < NUM_ADJ; dir++)
	{
		printf("In dir %d I have proc %d, he has x:%d y:%d\n", dir, neighbour_rank[dir], neighbour_nx[dir][0], neighbour_nx[dir][1]);
		fflush(stdout);
	} */
}


/* See macro
int get_opposite_dir(int dir)
{
	switch (dir)
	{
		case DOWN_LEFT:		return UP_RIGHT;
		case DOWN: 			return UP;
		case DOWN_RIGHT:	return UP_LEFT;
		case LEFT:			return RIGHT;
		case RIGHT:			return LEFT;
		case UP_LEFT:		return DOWN_RIGHT;
		case UP:			return DOWN;
		case UP_RIGHT:		return DOWN_LEFT;

		default:
			printf("ERROR: get_opposite_dir reached default."); fflush(stdout);
			return -1;
	}

	return abs( dir - (NUM_ADJ-1) );
}
*/

// returns 1 if proc can send/receive gc data to/from neighbour at direction dir, 0 otherwise
inline int can_send_gc(const int moving_window, const int dir)
{
	// restrictions only apply to moving window simulations
	if ( !moving_window )
		return 1;
	
	// if proc is on the left edge of the simulation space, dont send gc data to the left
	if ( (dir == UP_LEFT || dir == LEFT || dir == DOWN_LEFT) && proc_coords[0] == 0 )
		return 0;
	
	// if proc is on the right edge of the simulation space, dont send gc data to the right
	if ( (dir == UP_RIGHT || dir == RIGHT || dir == DOWN_RIGHT) && proc_coords[0] == dims[0] - 1 )
		return 0;
	
	return 1;
}