#ifndef SUCCESS_OR_DIE_H
#define SUCCESS_OR_DIE_H

#include <GASPI.h>
#include <stdio.h>
#include <stdlib.h>

#define SUCCESS_OR_DIE(...)																			\
	do																								\
	{																								\
		const gaspi_return_t r = __VA_ARGS__;														\
																									\
		if (r != GASPI_SUCCESS)																		\
		{																							\
			char* error = malloc(128);																\
			gaspi_print_error(r, &error);															\
			printf ("Error: '%s' [%s:%i]: %i => %s\n", #__VA_ARGS__, __FILE__, __LINE__, r, error);	\
			exit (EXIT_FAILURE);																	\
		}																							\
	} while (0)


#define SUCCESS_TIMEOUT_OR_DIE(...)																	\
	do																								\
	{																								\
		const gaspi_return_t r = __VA_ARGS__;														\
																									\
		if (r != GASPI_SUCCESS && r != GASPI_TIMEOUT)												\
		{																							\
			char* error = malloc(128);																\
			gaspi_print_error(r, &error);															\
			printf ("Error: '%s' [%s:%i]: %i => %s\n", #__VA_ARGS__, __FILE__, __LINE__, r, error);	\
			exit (EXIT_FAILURE);																	\
		}																							\
	} while (0)

#endif
