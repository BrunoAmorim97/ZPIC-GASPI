#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "success_or_die_fast.h"
#include <stdbool.h>

#define NUM_DIMS 2
#define ROOT 0

#define MAX_PPC_MULTIPLIER 3

#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_HIGH(id, p, n) - BLOCK_LOW(id, p, n) + 1)
// #define BLOCK_OWNER(index, p, n) (((p) * ((index) + 1) - 1) / (n))

#define NUM_ADJ 8
enum dir
{
	DOWN_LEFT = 0,
	DOWN = 1,
	DOWN_RIGHT = 2,
	LEFT = 3,
	RIGHT = 4,
	UP_LEFT = 5,
	UP = 6,
	UP_RIGHT = 7
};

#define DIR_TO_CURR_SEG_ID(dir) (dir + NUM_ADJ)
#define DIR_TO_CURR_KER_SEG_ID(dir) (dir + (2 * NUM_ADJ))
#define DIR_TO_EMF_SEG_ID(dir) (dir + (3 * NUM_ADJ))

#define OPPOSITE_DIR(dir) abs(dir - (NUM_ADJ - 1))

// Particle segment IDs equal the constant for that direction.
enum seg_id
{
	PAR_DOWN_LEFT = DOWN_LEFT,
	PAR_DOWN = DOWN,
	PAR_DOWN_RIGHT = DOWN_RIGHT,
	PAR_LEFT = LEFT,
	PAR_RIGHT = RIGHT,
	PAR_UP_LEFT = UP_LEFT,
	PAR_UP = UP,
	PAR_UP_RIGHT = UP_RIGHT,

	CURR_DOWN_LEFT = DIR_TO_CURR_SEG_ID(DOWN_LEFT),
	CURR_DOWN = DIR_TO_CURR_SEG_ID(DOWN),
	CURR_DOWN_RIGHT = DIR_TO_CURR_SEG_ID(DOWN_RIGHT),
	CURR_LEFT = DIR_TO_CURR_SEG_ID(LEFT),
	CURR_RIGHT = DIR_TO_CURR_SEG_ID(RIGHT),
	CURR_UP_LEFT = DIR_TO_CURR_SEG_ID(UP_LEFT),
	CURR_UP = DIR_TO_CURR_SEG_ID(UP),
	CURR_UP_RIGHT = DIR_TO_CURR_SEG_ID(UP_RIGHT),

	CURR_KER_DOWN_LEFT = DIR_TO_CURR_KER_SEG_ID(DOWN_LEFT),
	CURR_KER_DOWN = DIR_TO_CURR_KER_SEG_ID(DOWN),
	CURR_KER_DOWN_RIGHT = DIR_TO_CURR_KER_SEG_ID(DOWN_RIGHT),
	CURR_KER_LEFT = DIR_TO_CURR_KER_SEG_ID(LEFT),
	CURR_KER_RIGHT = DIR_TO_CURR_KER_SEG_ID(RIGHT),
	CURR_KER_UP_LEFT = DIR_TO_CURR_KER_SEG_ID(UP_LEFT),
	CURR_KER_UP = DIR_TO_CURR_KER_SEG_ID(UP),
	CURR_KER_UP_RIGHT = DIR_TO_CURR_KER_SEG_ID(UP_RIGHT),

	EMF_DOWN_LEFT = DIR_TO_EMF_SEG_ID(DOWN_LEFT),
	EMF_DOWN = DIR_TO_EMF_SEG_ID(DOWN),
	EMF_DOWN_RIGHT = DIR_TO_EMF_SEG_ID(DOWN_RIGHT),
	EMF_LEFT = DIR_TO_EMF_SEG_ID(LEFT),
	EMF_RIGHT = DIR_TO_EMF_SEG_ID(RIGHT),
	EMF_UP_LEFT = DIR_TO_EMF_SEG_ID(UP_LEFT),
	EMF_UP = DIR_TO_EMF_SEG_ID(UP),
	EMF_UP_RIGHT = DIR_TO_EMF_SEG_ID(UP_RIGHT),

	REPORTING,
	CURRENT_REPORT,
	EMF_B_REPORT,
	EMF_E_REPORT,
};

enum queue_id
{
	Q_PARTICLES,
	Q_CURRENT,
	Q_CURRENT_KERNEL,
	Q_EMF,
	Q_REPORTING
};

enum notification_id
{
	NOTIF_ID_CURRENT,
	NOTIF_ID_CURRENT_KERNEL_EVEN,
	NOTIF_ID_CURRENT_KERNEL_ODD,
	NOTIF_ID_EMF,
	NOTIF_ID_EMF_B,
	NOTIF_ID_EMF_E,
	NOTIF_ID_REPORTING
};

int *get_factors(int num, int *num_factors);
void assign_factors(int *factors, int num_factors, int dims[NUM_DIMS]);
void create_dims(int num_procs);
int cart_rank(int coords[NUM_DIMS]);
void cart_coords(int rank, int coords[NUM_DIMS]);
void assign_proc_blocks(int const nx[NUM_DIMS]);
void discover_neighbours(int proc_coords[NUM_DIMS], int dims[NUM_DIMS], int nx[NUM_DIMS]);
int periodic_coord(int coord, int coord_max);

bool can_talk_to_dir(const bool moving_window, const int dir);
int get_num_incoming_notifs(const bool moving_window);