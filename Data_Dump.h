#if !defined(DATA_DUMP_HEADER)
#define DATA_DUMP_HEADER

namespace Data_Dump {
  void Load_From_File(Particle * Particles,
                      unsigned int & Num_Particles);

  void Print_Data_To_File(const Particle * Particles,
                          const unsigned int Num_Particles);

  void Print_Particle_To_File(const Particle & P_In,
                              FILE * File);

  void Load_Particle_From_File(Particle * Particles,
                               const FILE * File);
} // namespace Data_Dump {

#endif
