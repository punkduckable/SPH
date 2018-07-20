#if !defined(PARTICLE_FRIENDS_HEADER)
#define PARTICLE_FRIENDS_HEADER

namespace Particle_Helpers {
  //////////////////////////////////////////////////////////////////////////////
  // Neighbor methods. Function definitions are in Particle_Neighbors.c
  void Find_Neighbors(const unsigned int Num_Particles,    // Generate neighbor list for every particle in 'Partilces' array
                      Particle * Particles);

  void Find_Neighbors_Box(Particle & P_In,                 // Generate neighbor list for a 'box' or 'cuboid' geometry
                          Particle * Particles,
                          const unsigned X_SIDE_LENGTH,
                          const unsigned Y_SIDE_LENGTH,
                          const unsigned Z_SIDE_LENGTH);

  bool Are_Neighbors(const Particle & P1,                  // Returns true P1 and P2 are neighbors, false otherwise
                     const Particle & P2);

  void Remove_Neighbor(Particle & P_In,                    // Removes a paricular one of P_In's neighbors.
                      const unsigned int Remove_Neighbor_ID,
                      const Particle * Particles);
  //////////////////////////////////////////////////////////////////////////////
  // Update methods. Function definitions are in Particle_Update.c
  void Update_P(Particle & P_In,                           // Updates P_In's Second Piola-Kirchhoff stress tensor
                Particle * Particles,
                const double dt);

  void Update_x(Particle & P_In,                           // Updates P_In's spacial position
                const Particle * Particles,
                const double dt);

  //////////////////////////////////////////////////////////////////////////////
  // Damage methods. Function definitions are in Particle_Damage.c
  void Remove_Damaged_Particle(Particle & P_In,
                               Particle * Particles);

  //////////////////////////////////////////////////////////////////////////////
  // Contact methods. Function definitions are in Particle_Contact.c

  double K = 100;                                                              //        : N/(mm^2)

  void Contact(Particle * Body_A, const unsigned int Num_Particles_A,
               Particle * Body_B, const unsigned int Num_Particles_B,
               const double h_In);

} // namespace Particle_Helpers {

#endif
