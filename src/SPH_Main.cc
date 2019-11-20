#include <iostream>
#include <math.h>
#include <string>
#include <cstring>
#include <random>
#include <unistd.h>
#include <omp.h>

// Type definitions, definitions
#define CLOCKS_PER_MS (CLOCKS_PER_SEC/1000.)               // Conversion from cpu cycles to msec
typedef signed char Byte;
typedef unsigned char uByte;
#if defined(_OPENMP)
  #define clock_t double
  #define clock() omp_get_wtime()
#else
  #include <time.h>

#endif

// Header files
#include "Classes.h"
#include "Namespaces.h"
#include "Materials.h"
#include "List.h"
#include "SPH_Diagnostics.h"
#include "VTK_File.h"
#include "Data_Dump.h"
#include "FEB_File.h"
#include "Quick_Math.h"
#include "Vector.h"
#include "Tensor.h"
#include "Particle_Helpers.h"
#include "Particle.h"
#include "Particle_Array.h"
#include "Tests.h"
#include "Simulation.h"

// Source files
#include "SPH_Diagnostics.cc"
#include "VTK_File.cc"
#include "Vector.cc"
#include "Tensor.cc"
#include "Particle.cc"
#include "Particle_Neighbors.cc"
#include "Particle_Update.cc"
#include "Particle_Damage.cc"
#include "Particle_Contact.cc"
#include "Particle_Array.cc"
#include "Data_Dump.cc"
#include "FEB_File.cc"
#include "Tests.cc"
#include "Simulation.cc"

int main() {
  //////////////////////////////////////////////////////////////////////////////
  // Tests
  //Vector_Tests();
  //Tensor_Tests();
  //List_Tests();
  //Particle_Tests();
  //Timing_Tests();

  Simulation::Run_Simulation();

  return 0;
} // int main() {