#include "Simulation.h"
#if defined(_OPENMP)
  #include <omp.h>
#endif

/* File description:
This file exits because my simulations use a lot of timers. Before I wrote
these functions, my simulation code repeated the same 8 or so lines of code
about 20 times. These functions therefore exist to reduce code repetition. */

TIME_TYPE Simulation::Get_Time(void) {
  /* This function returns the current time. What gets reported depends on
  whether or not the simulation is using OpenMP */

  #if defined(_OPENMP)
    return omp_get_wtime();
  #else
    return ((double)clock())/((double)CLOCKS_PER_SEC);
  #endif
} // TIME_TYPE Simulation::Get_Time(void) {



TIME_TYPE Simulation::Time_Since(TIME_TYPE time) {
  /* This function is used to compute the difference between the current time
  and the time argument. In other words, it calculates the time since the
  time stored in the time argument. The way that we do this, however, depends on
  whether OpenMP is defined or not. */

  #if defined(_OPENMP)
    return (omp_get_wtime() - time);
  #else
    return ( (((double)clock())/((double)CLOCKS_PER_SEC)) - time);
  #endif
} // TIME_TYPE Simulation::Time_Since(TIME_TYPE time) {
