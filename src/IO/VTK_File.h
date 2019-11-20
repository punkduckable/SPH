#if !defined(VTK_FILE_HEADER)
#define VTK_FILE_HEADER

namespace VTK_File {
  // List to keep track of which arrays we've seen before and how many times
  // we've printed data from each
  List<std::string> Name_List;
  List<unsigned int> File_Number_List;
  const unsigned int File_Number_Max_Digits = 5;

  std::string Get_File_Name(const std::string & Str);

  void Add_Point_Data(FILE * File,
                      char * Weight_Name,
                      unsigned int Num_Particles,
                      double * Data);

  void Export_Particle_Positions(const Particle_Array & Particles);
} // namespace VTK_File {

#endif