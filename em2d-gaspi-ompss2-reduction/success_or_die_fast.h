#ifndef SUCCESS_OR_DIE_FAST_H
#define SUCCESS_OR_DIE_FAST_H

#include <GASPI.h>
#include <stdlib.h>


#define SUCCESS_OR_DIE(...)								\
	do													\
	{													\
		const gaspi_return_t r = __VA_ARGS__;			\
														\
		if (r != GASPI_SUCCESS)							\
		{												\
			exit (EXIT_FAILURE);						\
		}												\
	} while (0)


#define SUCCESS_TIMEOUT_OR_DIE(...)						\
	do													\
	{													\
		const gaspi_return_t r = __VA_ARGS__;			\
														\
		if (r != GASPI_SUCCESS && r != GASPI_TIMEOUT)	\
		{												\
			exit (EXIT_FAILURE);						\
		}												\
	} while (0)


#endif
