#if !defined(PARTICLE_FILE_HEADER)
#define PARTICLE_FILE_HEADER

namespace Particle_Debugger_File {
  using std::string;

  unsigned int File_Number = 0;
  const unsigned int File_Number_Max_Digits = 5;

  void Export_Pariticle_Properties(const unsigned int Num_Particles, const Particle * Particles);
  void Get_File_Name(string & Str);
} // namespace VTK_File {

#endif
