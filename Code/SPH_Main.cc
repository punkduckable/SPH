#include <iostream>
#include <math.h>
#include <time.h>
#include <string>
#include <cstring>
#include <random>
#include <unistd.h>

// Type definitions
typedef signed char Byte;
typedef unsigned char uByte;

// Header files
#include "Classes.h"
#include "Namespaces.h"
#include "SPH_Diagnostics.h"
#include "VTK_File.h"
#include "Data_Dump.h"
#include "FEB_File.h"
#include "Quick_Math.h"
#include "Vector.h"
#include "Tensor.h"
#include "Particle_Helpers.h"
#include "Particle.h"
#include "List.h"
#include "Tests.h"
#include "Simulation.h"

// Prototypes
Tensor Dyadic_Product(const Vector & V1,const Vector & V2);

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
#include "Data_Dump.cc"
#include "FEB_File.cc"
#include "Tests.cc"
#include "Simulation.cc"

////////////////////////////////////////////////////////////////////////////////
// Function definitions

Tensor Dyadic_Product(const Vector & V1,const Vector & V2) {
  Tensor S;

  /* Assign the elements of our dyadic product using nested for loop. Note that
     We only cycle through the columns. Let S denote the dyadic product of V1
     and v2. The jth column of S is equal to V1*V2[j] (scale the vector V1 by
     the jth component of V2)
  */

  /* Unrolled loop (runs slower than 1 rolled loop for some reason)
  T[3*0 + 0] = V1[0]*V2[0];            // i = 0, j = 0
  T[3*0 + 1] = V1[0]*V2[1];            // i = 0, j = 1
  T[3*0 + 2] = V1[0]*V2[2];            // i = 0, j = 2

  T[3*1 + 0] = V1[1]*V2[0];            // i = 1, j = 0
  T[3*1 + 1] = V1[1]*V2[1];            // i = 1, j = 1
  T[3*1 + 2] = V1[1]*V2[2];            // i = 1, j = 2

  T[3*2 + 0] = V1[2]*V2[0];            // i = 2, j = 0
  T[3*2 + 1] = V1[2]*V2[1];            // i = 2, j = 1
  T[3*2 + 2] = V1[2]*V2[2];            // i = 2, j = 2
  */

  // Old loop. (works better than 9 statements with O2 optimization)
  for(int i = 0; i < 3; i++) {
    S[3*i + 0] = V1[i]*V2[0];
    S[3*i + 1] = V1[i]*V2[1];
    S[3*i + 2] = V1[i]*V2[2];
  } //   for(int i = 0; i < 3; i++)

  //OP_Count::Dyadic_Product++;                    // Increment operator count (See SPH Diagnostics)

  return S;
} // Tensor Dyatic_Product(const Vector & V2,const Vector & V2) {

int main() {
  //////////////////////////////////////////////////////////////////////////////
  // Tests
  //Vector_Tests();
  //Tensor_Tests();
  //List_Tests();
  //Particle_Tests();
  //Timing_Tests();

  Simulation::Run_Simulation();

  Vector * X = NULL;
  unsigned int Num_Nodes;
  FEB_File::Read_FEB_File("Needle.feb", &X, Num_Nodes);

  printf("%p\n",X);
  printf("%u\n",Num_Nodes);
  for(unsigned int i = 0; i < Num_Nodes; i++) {
    X[i].Print();
  }

  return 0;
} // int main() {
