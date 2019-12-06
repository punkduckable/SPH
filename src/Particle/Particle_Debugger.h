#if !defined(PARTICLE_FILE_HEADER)
#define PARTICLE_FILE_HEADER

#include "Particle.h"

// If defined, PARTICLE_DEBUG adds the Visc and Force_Visc members to the
// Tensor class. These members are for debugging the particle class.
#define PARTICLE_DEBUG

namespace Particle_Debugger {
  List<std::string> Name_List;
  List<unsigned> File_Number_List;
  const unsigned File_Number_Max_Digits = 5;

  void Export_Particle_Forces(const Body & Particles);

  std::string Get_File_Name(const std::string & Str);
} // namespace VTK_File {


#endif
