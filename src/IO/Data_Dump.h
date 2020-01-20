#if !defined(DATA_DUMP_HEADER)
#define DATA_DUMP_HEADER

#define LOAD_MONITOR

#include "Classes.h"
#include <stdio.h>

namespace Data_Dump {
  void Save_Simulation(const Body * Bodies,
                       const unsigned Num_Bodies);

  void Save_Body(const Body & Body_In);

  void Save_Particle(const Particle & P_In,
                     FILE * File);

  void Load_Simulation(Body ** Array_Ptr,
                      unsigned & Num_Bodies);

  void Load_Body(Body & Body_In);

  void Load_Particle(Particle & P_In,
                     FILE * File);
} // namespace Data_Dump {

#endif
