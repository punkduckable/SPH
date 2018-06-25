#if !defined(VTK_FILE_HEADER)
#define VTK_FILE_HEADER

namespace VTK_File {
  using std::string;

  unsigned int File_Number = 0;
  const unsigned int File_Number_Max_Digits = 5;

  bool Append_Digits(unsigned int N, string & Str);
  void Export_Pariticle_Positions(const unsigned int Num_Particles, const Particle * Particles);
  void Get_File_Name(string & Str);
} // namespace VTK_File {

#endif
