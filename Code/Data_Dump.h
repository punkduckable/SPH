#if !defined(DATA_DUMP_HEADER)
#define DATA_DUMP_HEADER

namespace Data_Dump {
  void Print_Particle_Array_To_File(const Particle_Array & Particles);

  void Print_Particle_To_File(const Particle & P_In,
                              FILE * File);

  int Load_Saved_Data(Particle_Array ** Array_Ptr,
                      unsigned int & Num_Bodies);

  int Load_Particle_Array_From_File(Particle_Array & Particles);

  void Load_Particle_From_File(Particle & P_In,
                               FILE * File);
} // namespace Data_Dump {

#endif
