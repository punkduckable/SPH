#if !defined(LOAD_SIMULATION_HEADER)
#define LOAD_SIMULATION_HEADER

#define LOAD_MONITOR

#include "Classes.h"

namespace IO {
  void Load_Simulation(Body ** Array_Ptr,
                       unsigned & Num_Bodies);

  void Load_Body(Body & Body_In);

  void Load_Particle(Particle & P_In,
                     std::ifstream & File);
} // namespace IO {

#endif
