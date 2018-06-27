#if !defined(VTK_FILE_HEADER)
#define VTK_FILE_HEADER

namespace VTK_File {
  using std::string;

  unsigned int File_Number = 0;
  const unsigned int File_Number_Max_Digits = 5;

  void Get_File_Name(string & Str);
  void Add_Point_Data(FILE * File, char * Weight_Name, unsigned int Num_Particles, double * Data);
  void Export_Pariticle_Positions(const unsigned int Num_Particles, const Particle * Particles);
} // namespace VTK_File {

#endif
