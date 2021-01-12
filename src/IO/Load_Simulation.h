#if !defined(LOAD_SIMULATION_HEADER)
#define LOAD_SIMULATION_HEADER

#define LOAD_MONITOR

#include "Classes.h"

namespace IO {
  void Load_Simulation(Body ** Bodies_Ptr,
                       unsigned & Num_Bodies);

  /* Why is this function not static? Because it needs to be a friend of the
  Particle class (it modifies members of the particle object that do not have
  setters)*/
  void Load_Particle(Particle & P_In,
                     std::ifstream & File);
} // namespace IO {

#endif
