#if !defined(PARTICLE_FRIENDS_HEADER)
#define PARTICLE_FRIENDS_HEADER

namespace Particle_Helpers {
  //////////////////////////////////////////////////////////////////////////////
  // Neighbor methods. Function definitions are in Particle_Neighbors.c
  void Find_Neighbors(Particle_Array & Particles);         // Generate neighbor list for every particle in 'Partilces' array

  void Find_Neighbors_Box(Particle & P_In,                 // Generate neighbor list for a 'box' or 'cuboid' geometry
                          Particle_Array & Particles);

  bool Are_Neighbors(const Particle_Array & Particles,     // Returns true P1 and P2 are neighbors, false otherwise
                     const unsigned int i,
                     const unsigned int j);

  void Remove_Neighbor(Particle & P_In,                    // Removes a paricular one of P_In's neighbors.
                      const unsigned int Remove_Neighbor_ID,
                      const Particle_Array & Particles);
  //////////////////////////////////////////////////////////////////////////////
  // Update methods. Function definitions are in Particle_Update.c
  void Update_P(Particle_Array & Particles,               // Updates P_In's Second Piola-Kirchhoff stress tensor
               const double dt);

  void Update_x(Particle_Array & Particles,               // Updates P_In's spacial position
                const double dt);

  //////////////////////////////////////////////////////////////////////////////
  // Damage methods. Function definitions are in Particle_Damage.c
  void Remove_Damaged_Particle(Particle & P_In,
                               Particle_Array & Particles);

  //////////////////////////////////////////////////////////////////////////////
  // Contact methods. Function definitions are in Particle_Contact.c

  const double K = 400;                                                               //        : N/(mm^2)

  void Contact(Particle_Array & Body_A,
               Particle_Array & Body_B);

} // namespace Particle_Helpers {

#endif
