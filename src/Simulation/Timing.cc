#include "Simulation.h"
#if defined(_OPENMP)
  #include <omp.h>
#endif

/* File description:
This file exits because my simulations use a lot of timers. Before I wrote
these functions, my simulation code repeated the same 8 or so lines of code
about 20 times. These functions therefore exist to reduce code repetition. */

TIME_TYPE Simulation::Get_Time(void) {
  #if defined(_OPENMP)
    return omp_get_wtime();
  #else
    return clock();
  #endif
} // TIME_TYPE Simulation::Get_Time(void) {



TIME_TYPE Simulation::Time_Since(TIME_TYPE time) {
  #if defined(_OPENMP)
    return (omp_get_wtime() - time);
  #else
    return (clock() - time);
  #endif
} // TIME_TYPE Simulation::Time_Since(TIME_TYPE time) {
