# ZPIC-GASPI

Distributed implementation of the Particle-in-cell code [ZPIC](https://github.com/ricardo-fonseca/zpic), using the [GPI-2](http://www.gpi-site.com/) implementation of the [GASPI](http://www.gaspi.de/) standard as the communication API. A GASPI/OpenMP hybrid version is also available.

## Requirements

- GASPI installation (regular or MPI interoperability configuration (recommended)).

## Compilation and Execution

- Change input file in main.c to desired simulation input.

### GASPI MPI interoperability (recommended)
```
make
mpirun -n <NUM PROCS> zpic
```

### GASPI (non-MPI interoperability)

- Comment MPI references in main.c.
- Change compiler in Makefile to non-MPI compiler.
- Create a machine file with the names of the participating nodes (as seen [here](https://github.com/cc-hpc-itwm/GPI-2)).

```
make
gaspi_run -m <MACHINE FILE> zpic
```

### GASPI\OpenMP hybrid

Same as described before, force the number of threads by command-line argument:
```
(...) zpic <NUM THREADS PER PROC>
```