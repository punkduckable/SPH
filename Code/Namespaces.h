#if !defined(NAMESPACES_HEADER)
#define NAMESPACES_HEADER

/* Thi file exists to serve as a sort of 'forward declaration' for all the
namespaces in my program. Remember that namespaces are open.

Why does this file need to exist?
Good question, the reason is scoping. This file was created because of issues
with the 'Simulation' namespace. The simulation namespace modifies static
Particle class members. Because of this, these parts of the namespace needed
to be defined AFTER the particle class was defined (so Simulation.h comes
after Particle.h). The issue is, many of the functions in the Simulation
namespace are friends of the Particle class. Thus, some of the Simulation
class functions needed to be declared (prototyped) before the particle class
definition. My resolution to this issue was to make a new file that can be used
to house the parts of each namespace that needs to come before everything else
(such as the Run_Simulation prototype for the Simulation namespace)*/

namespace Simulation {
  void Run_Simulation(void);
  void Set_Static_Particle_Members(void);
}
namespace Data_Dump {}
namespace FEB_File {}
namespace Particle_Helpers {}
namespace Quick_Math {}
namespace Particle_Debugger {}
namespace OP_Count {}
namespace VTK_File {}

#endif
