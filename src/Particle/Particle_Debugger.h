#if !defined(PARTICLE_FILE_HEADER)
#define PARTICLE_FILE_HEADER

#include "Particle.h"

namespace Particle_Debugger {
  List<std::string> Name_List;
  List<unsigned int> File_Number_List;
  const unsigned int File_Number_Max_Digits = 5;

  void Export_Particle_Forces(const Particle_Array & Particles);

  std::string Get_File_Name(const std::string & Str);
} // namespace VTK_File {


#endif
