#if !defined(DATA_DUMP_HEADER)
#define DATA_DUMP_HEADER

namespace Data_Dump {
  void Save_Simulation(const Particle_Array * Arrays,
                       const unsigned int Num_Arrays);

  void Save_Particle_Array(const Particle_Array & Particles);

  void Save_Particle(const Particle & P_In,
                     FILE * File,
                     const bool Is_Cuboid);

  int Load_Simulation(Particle_Array ** Array_Ptr,
                      unsigned int & Num_Bodies);

  int Load_Particle_Array_From_File(Particle_Array & Particles);

  void Load_Particle_From_File(Particle & P_In,
                               FILE * File);
} // namespace Data_Dump {

#endif
