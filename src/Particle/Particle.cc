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
// Getters

void Particle::Set_ID(const unsigned ID_In) { ID = ID_In; }
void Particle::Set_Mass(const double Mass_In) { Mass = Mass_In; }
void Particle::Set_Vol(const double Vol_In) { Vol = Vol_In; }
void Particle::Set_Radius(const double Radius_In) { Radius = Radius_In; }

void Particle::Set_X(const Vector & X_In) { X = X_In; }
void Particle::Set_x(const Vector & x_In) { x = x_In; }
void Particle::Set_V(const Vector & V_In) { V = V_In; }
void Particle::Set_a(const Vector & a_In) { a = a_In; }
void Particle::Set_D(const double D_In) { D = D_In; }





////////////////////////////////////////////////////////////////////////////////
// Getters

unsigned Particle::Get_ID(void) const { return ID; }
double Particle::Get_Mass(void) const { return Mass; }
double Particle::Get_Vol(void) const { return Vol; }
double Particle::Get_Radius(void) const { return Radius; }

const Vector & Particle::Get_X(void) const { return X; }
const Vector & Particle::Get_x(void) const { return x; }
const Vector & Particle::Get_V(void) const { return V; }
const Vector & Particle::Get_a(void) const { return a; }
const Tensor & Particle::Get_P(void) const { return P; }
const Tensor & Particle::Get_F(const unsigned i) const {
  assert(i < 2);
  return F[i];
} // const Tensor & Particle::Get_F(const unsigned i) const {

const Vector & Particle::Get_Force_Friction(void) const { return Force_Friction; }
const Vector & Particle::Get_Force_Contact(void) const { return Force_Contact; }

double Particle::Get_Stretch_M(void) const { return Stretch_M; }
double Particle::Get_Stretch_H(void) const { return Stretch_H; }
double Particle::Get_Stretch_Critical(void) const { return Stretch_Critical; }
double Particle::Get_D(void) const { return D; }

unsigned Particle::Get_Num_Neighbors(void) const { return Num_Neighbors; }
unsigned Particle::Get_Neighbor_IDs(unsigned i) const {
  if(i < Num_Neighbors) { return Neighbor_IDs[i]; }
  else {
    char Error_Message_Buffer[500];
    sprintf(Error_Message_Buffer,
            "Bad Neighbor Index exception: Thrown by Particle::Get_Neighbor_IDs\n"
            "Particle %u has %u neighbors. You requested neighbor %u\n",
            ID, Num_Neighbors, i);
    throw Bad_Neighbor_Index(Error_Message_Buffer);
  } // else
} // unsigned Particle::Get_Neighbor_IDs(unsigned i) const {





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



void Print(const Particle & P_In) { P_In.Print(); };         // Calls Print method
