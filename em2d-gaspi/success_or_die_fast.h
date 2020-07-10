#ifndef SUCCESS_OR_DIE_FAST_H
#define SUCCESS_OR_DIE_FAST_H

#include <GASPI.h>
#include <stdlib.h>


#define SUCCESS_OR_DIE(...) if ( __VA_ARGS__ != GASPI_SUCCESS) exit (EXIT_FAILURE);


#endif
