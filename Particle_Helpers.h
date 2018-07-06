#if !defined(PARTICLE_FRIENDS_HEADER)
#define PARTICLE_FRIENDS_HEADER

namespace Particle_Helpers {
  //////////////////////////////////////////////////////////////////////////////
  // Neighbor methods. Function definitions are in Particle_Neighbors.c
  void Find_Neighbors(const unsigned int Num_Particles,
                      Particle * Particles);               // Generate neighbor list for every particle in 'Partilces' array

  void Find_Neighbors_Box(Particle & P_In, Particle * Particles);// Generate neighbor list for a 'box' or 'cuboid' geometry

  bool Are_Neighbors(const Particle & P1,
                            const Particle & P2);          // Returns true P1 and P2 are neighbors, false otherwise

  //////////////////////////////////////////////////////////////////////////////
  // Update methods. Function definitions are in Particle_Update.c
  void Update_P(Particle & P_In,
                Particle * Particles,
                const double dt);                          // Updates P_In's Second Piola-Kirchhoff stress tensor

  void Update_x(Particle & P_In,
                const Particle * Particles,
                const double dt);                          // Updates P_In's spacial position

  //////////////////////////////////////////////////////////////////////////////
  // Damage methods. Function definitions are in Particle_Damage.c
  void Remove_Damaged_Particle(Particle & P_In, Particle * Particles);

  //////////////////////////////////////////////////////////////////////////////
  // Contact methods. Function definitions are in Particle_Contact.c

  double K = 10;

  void Contact(Particle * Body_A, const unsigned int Num_Particles_A,
               Particle * Body_B, const unsigned int Num_Particles_B,
               const double h_particle, const double IPS);

} // namespace Particle_Helpers {

#endif
