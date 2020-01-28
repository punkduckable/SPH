#if !defined(SAVE_SIMULATION_HEADER)
#define SAVE_SIMULATION_HEADER

#include "Classes.h"
#include <stdio.h>

namespace IO {
  void Save_Simulation(const Body * Bodies,
                       const unsigned Num_Bodies);
} // namespace IO {

#endif
