#if !defined(PARTICLE_SOURCE)
#define PARTICLE_SOURCE

#include "Particle.h"
#include "Tensor.h"
#include "Vector.h"
#include "List.h"

////////////////////////////////////////////////////////////////////////////////
// Constructors and destructor

Particle::Particle(void) {
  //printf("Particle default constructor \n");
  Has_Neighbors = false;
  Num_Neighbors = 0;
  Vol = 0;                                                                     //        : mm^3
  Mass = 0;                                                                    //        : g

  // Now randomly set critical stress
  unsigned seed = std::rand();
  std::default_random_engine generator (seed);
  std::normal_distribution<double> distribution(1.3,.02);
  Stretch_Critical = distribution(generator);
} // Particle::Particle(void) {

Particle::Particle(const Particle & P_In) {
  // Do nothing copy constructor
  printf("Bad! You tried using the particle copy constructor. Don't do that!\n");
} // Particle::Particle(const Particle & P_In) {

Particle::~Particle(void) {
  //printf("Removing particle\n");

  // Note, we should only free the memory if it has been allocated.
  if(Has_Neighbors == true) {
    delete [] R;                                                               //        : mm
    delete [] Mag_R;                                                           //        : mm
    delete [] W;                                                               //        : unitless
    delete [] Grad_W;                                                          //        : mm^-1
    delete [] Neighbor_IDs;
  }
} // Particle::~Particle(void) {



////////////////////////////////////////////////////////////////////////////////
//  Particle equality

Particle & Particle::operator=(const Particle & P_In) {
  // Do nothing = operator overload.
  printf("You can't use equality with particles!\n");
  return *this;
}



////////////////////////////////////////////////////////////////////////////////
// Particle setup methods:
// Set Neighbors

void Particle::Set_Neighbors(const unsigned int N, const unsigned int * Neighbor_ID_Array, const Particle * Particles) {
  /* First check if this particle already has neighbors. This function should
  only be called if the neighbors have not been set. The reason for this is
  that this method allocates pointers. If the pointers have already been set,
  then allocating them again will cause a memory leak.  */
  if(Has_Neighbors == true) {
    printf("Neigbors already set, returning.\n");
    return;
  }

  /* Next, check that N > 0. if N = 0, then there are no neighbors to add. */
  if(N == 0) {
    printf("You didn't supply any neighbors!\n");
    return;
  }

  // Set Num_Neighbors using input
  Num_Neighbors = N;

  // Allocate memory for the Dynamic arrays
  Neighbor_IDs = new unsigned int[Num_Neighbors];
  R = new Vector[Num_Neighbors];                                               //        : mm
  Mag_R = new double[Num_Neighbors];                                           //        : mm
  W = new double[Num_Neighbors];                                               //        : unitless
  Grad_W = new Vector[Num_Neighbors];                                          //        : mm^-1

  /* Now that we know our neighbors IDs, we can figure out everything that we
  want to know about them. We an set the Neighbor_IDs, r, R, W, and Grad_W
  members. These can be used to calculate the shape matrix (and its inverse)! */

  int Neighbor_ID;                               // Keep track of current particle
  double V_j;                                    // Volume of jth neighbor               : mm^3
  Tensor A{0,0,0,
           0,0,0,
           0,0,0};                               // Shape Tensor (zero initialized)      : unitless

  for(unsigned int j = 0; j < Num_Neighbors; j++) {
    Neighbor_ID = Neighbor_ID_Array[j];          // Get Neighbor ID (index in Particles array)
    Neighbor_IDs[j] = Neighbor_ID;               // Set jth element of Neighbor_IDs member

    // Calculate displacement vectors
    R[j] = Particles[Neighbor_ID].X - X;         // Reference displacement vector        : mm
    Mag_R[j] = R[j].Magnitude();                 // |R[j]|                               : mm

    // Calculate shape function, shape function gradient for jth neighbor
    W[j] = Shape_Function_Amp*(h - Mag_R[j])*(h - Mag_R[j])*(h - Mag_R[j]);                             //        :unitless
    Grad_W[j] = -3*Shape_Function_Amp*((h - Mag_R[j])*(h - Mag_R[j]))*(R[j] / Mag_R[j]); //    : mm^-1

    // Add in the Current Neighbor's contribution to the Shape tensor
    V_j = Particles[Neighbor_ID].Vol;            // Neighbor Volume                      : mm^3
    A += Dyadic_Product((V_j*Grad_W[j]), R[j]);                                //        : unitless
  } // for(unsigned int j = 0; j < N; j++) {

  // Now we can calculate A^(-1) from A.
  A_Inv = A^(-1);                                                              //        : unitless

  // Now that neighbors have been set, we set 'Has_Neighbors' to true
  Has_Neighbors = true;
} // void Particle::Set_Neighbors(const unsigned int N, const unsigned int *Neighbor_Id_List, const Particle *Particles) {



////////////////////////////////////////////////////////////////////////////////
// Printing functions

void Particle::Print(void) const {
  // Print basic particle parameters.
  printf("X:   ");
  X.Print();
  //printf("x:   ");
  //x.Print();
  //printf("vel: ");
  //vel.Print();
  //printf("F:   \n");
  //F.Print();
  //printf("P:   \n");
  //P.Print();
  //printf("A^(-1)\n");
  //A_Inv.Print();
  printf("Num Neighbors: %d\n",Num_Neighbors);
  printf("F_Int = ");
  Force_Int.Print();
  printf("F_Visc = ");
  //Force_Visc.Print();                          // For debugging
  //printf("F_Hg = ");
  Force_Hg.Print();
  printf("\n");

  // If we have neighbors, print neighbor information
  if(Has_Neighbors == true) {
    //unsigned int i;                    // Loop index variable

    /* Print neighbor ID's
    printf("Neighbor ID's  : {");
    for(i = 0; i < Num_Neighbors-1; i++) {
      printf("%5d, ",Neighbor_IDs[i]);
    } // for(i = 0; i < Num_Neighbors-1; i++) {
    printf("%5d } \n", Neighbor_IDs[Num_Neighbors-1]);  // */

    /* Print Grad_W magnitudes
    printf("%p\n",Grad_W);
    printf("|Grad_W|       : {");
    for(unsigned int i = 0; i < Num_Neighbors-1; i++) {
      printf("%5.3f, ",Magnitude(Grad_W[i]));
    } // for(unsigned int i = 0; i < Num_Neighbors-1; i++) {
    printf("%5.3f } \n", Magnitude(Grad_W[Num_Neighbors-1])); // */

  } // if(Has_Neighbors == true) {
} // void Particle::Print(void) const {

void Print(const Particle & P_In) {
  P_In.Print();
} // void Print(const Particle & P_In) {

#endif
