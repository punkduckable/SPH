#if !defined(PARTICLE_HELPERS_HEADER)
#define PARTICLE_HELPERS_HEADER

#include "Namespaces.h"
#include "Particle.h"
#include "Body/Body.h"
#include "List.h"

namespace Particle_Helpers {
  //////////////////////////////////////////////////////////////////////////////
  // Update methods. Function definitions are in Particle_Update.c
  void Update_P(Body & Particles,               // Updates P_In's Second Piola-Kirchhoff stress tensor
                const double dt);

  void Update_x(Body & Particles,               // Updates P_In's spacial position
                const double dt);


  //////////////////////////////////////////////////////////////////////////////
  // Damage methods. Function definitions are in Particle_Damage.c
  void Remove_Damaged_Particle(Particle & P_In,
                               Body & Particles);


  //////////////////////////////////////////////////////////////////////////////
  // Contact methods. Function definitions are in Particle_Contact.c

  const double K = 400;                                                        //        : N/(mm^2)

  void Contact(Body & Body_A,
               Body & Body_B);

  unsigned Times_Printed_Net_External_Force = 0;

  void Print_Net_External_Force(const Body & Particles,
                                const unsigned l);

} // namespace Particle_Helpers {

#endif
