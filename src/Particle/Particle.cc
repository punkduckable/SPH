#include "Particle.h"
#include <random>

////////////////////////////////////////////////////////////////////////////////
// Constructors and destructor

Particle::Particle(void) {
  //printf("Particle default constructor \n");
  Neighbors_Are_Set = false;
  Num_Neighbors = 0;
  Vol = 0;                                                                     //        : mm^3
  Mass = 0;                                                                    //        : g

  // Now randomly set critical stress
  unsigned seed = std::rand();
  std::default_random_engine generator (seed);
  std::normal_distribution<double> distribution(1.3,.00);
  Stretch_Critical = distribution(generator);
} // Particle::Particle(void) {



Particle::~Particle(void) {
  //printf("Removing particle\n");

  // Note, we should only free the memory if it has been allocated.
  if(Neighbors_Are_Set == true) {
    delete [] R;                                                               //        : mm Vectro
    delete [] Mag_R;                                                           //        : mm
    delete [] W;                                                               //        : unitless
    delete [] Grad_W;                                                          //        : 1/mm Vector
    delete [] Neighbor_IDs;
  } // if(Neighbors_Are_Set == true) {
} // Particle::~Particle(void) {





////////////////////////////////////////////////////////////////////////////////
// Printing functions

void Particle::Print(void) const {
  // Print basic particle parameters.
  printf("X:   ");
  (*this).X.Print();
  printf("x:   ");
  (*this).x.Print();
  printf("vel: ");
  (*this).V.Print();

  printf("F[0]:   \n");
  (*this).F[0].Print();
  printf("F[1]:   \n");
  (*this).F[1].Print();
  printf("P:   \n");
  (*this).P.Print();
  printf("A^(-1)\n");
  (*this).A_Inv.Print();

  printf("F_Int = ");
  (*this).Force_Int.Print();
  #if defined(PARTICLE_DEBUG)
    printf("F_Visc = ");
    (*this).Force_Visc.Print();
  #endif
  printf("F_Hg = ");
  (*this).Force_HG.Print();
  printf("\n");

  // If we have neighbors, print neighbor information
  if(Neighbors_Are_Set == true) {
    printf("Num Neighbors: %d\n",Num_Neighbors);
    //unsigned i;                              // Loop index variable

    /* Print neighbor ID's */
    printf("Neighbor ID's  : {");
    for(unsigned i = 0; i < (*this).Num_Neighbors-1; i++) {
      printf("%5d, ",(*this).Neighbor_IDs[i]);
    } // for(i = 0; i < Num_Neighbors-1; i++) {
    printf("%5d } \n", (*this).Neighbor_IDs[Num_Neighbors-1]); // */

    /* Print Grad_W magnitudes */
    printf("%p\n",Grad_W);
    printf("|Grad_W|       : {");
    for(unsigned i = 0; i < (*this).Num_Neighbors-1; i++) {
      printf("%5.3f, ",Magnitude((*this).Grad_W[i]));
    } // for(unsigned i = 0; i < Num_Neighbors-1; i++) {
    printf("%5.3f } \n", Magnitude((*this).Grad_W[Num_Neighbors-1])); // */
  } // if(Neighbors_Are_Set == true) {
} // void Particle::Print(void) const {
