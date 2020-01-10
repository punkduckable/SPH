#if !defined(DATA_DUMP_HEADER)
#define DATA_DUMP_HEADER

#include "Classes.h"
#include <stdio.h>

namespace Data_Dump {
  void Save_Simulation(const Body * Bodies,
                       const unsigned Num_Bodies);

  void Save_Body(const Body & Body_In);

  void Save_Particle(const Particle & P_In,
                     FILE * File);

  int Load_Simulation(Body ** Array_Ptr,
                      unsigned & Num_Bodies);

  int Load_Body(Body & Body_In);

  void Load_Particle(Particle & P_In,
                     FILE * File);
} // namespace Data_Dump {

#endif
