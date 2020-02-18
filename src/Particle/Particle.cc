#include "Particle.h"
#include "Body/Body.h"
#include "Errors.h"
#include <random>
#include <assert.h>

////////////////////////////////////////////////////////////////////////////////
// Constructors and destructor

Particle::Particle(void) {
  (*this).Neighbors_Are_Set = false;
  (*this).Num_Neighbors = 0;
  (*this).Mass_Set = false;
  (*this).Volume_Set = false;
  (*this).Radius_Set = false;
  for(unsigned i = 0; i < 3; i++) { (*this).Has_BC[i] = false; }

  // Now randomly set critical stress
  unsigned seed = std::rand();
  std::default_random_engine generator (seed);
  std::normal_distribution<double> distribution(1.3,.00);
  Stretch_Critical = distribution(generator);
} // Particle::Particle(void) {



Particle::~Particle(void) {

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
// Constructors and destructor

void Particle::Apply_BCs(void) {
  // This function applies the BCs as specified in the BC Vector.
  for(unsigned i = 0; i < 3; i++) {
    if((*this).Has_BC[i] == true) { (*this).V[i] = BC[i]; }
  } // for(unsigned i = 0; i < 3; i++) {
} // void Particle::Apply_BCs(void) {






////////////////////////////////////////////////////////////////////////////////
// Setters

void Particle::Set_ID(const unsigned ID_In) { ID = ID_In; }

void Particle::Set_Mass(const double Mass_In) {
  assert(Mass_In != 0 );
  (*this).Mass = Mass_In;
  (*this).Mass_Set = true;
} // void Particle::Set_Mass(const double Mass_In) {
void Particle::Set_Volume(const double Volume_In) {
  assert(Volume_In != 0);
  (*this).Volume = Volume_In;
  (*this).Volume_Set = true;
} // void Particle::Set_Volume(const double Volume_In) {
void Particle::Set_Radius(const double Radius_In) {
  assert(Radius_In != 0);
  (*this).Radius = Radius_In;
  (*this).Radius_Set = true;
} // void Particle::Set_Radius(const double Radius_In) {

void Particle::Set_X(const Vector & X_In) { X = X_In; }
void Particle::Set_x(const Vector & x_In) { x = x_In; }
void Particle::Set_V(const Vector & V_In) { V = V_In; }
void Particle::Set_a(const Vector & a_In) { a = a_In; }

void Particle::Set_D(const double D_In) { D = D_In; }

void Particle::Set_BC(const unsigned Component, const double Value) {
  assert(Component < 3);

  // Now set the BC
  BC[Component] = Value;
  Has_BC[Component] = true;
} // void Particle::Set_BC(const unsigned Component, const double Value) {



////////////////////////////////////////////////////////////////////////////////
// Getters

unsigned Particle::Get_ID(void) const { return ID; }

double Particle::Get_Mass(void) const {
  assert((*this).Mass_Set == true);
  return (*this).Mass;
} // double Particle::Get_Mass(void) const {
double Particle::Get_Volume(void) const {
  assert((*this).Volume_Set == true);
  return (*this).Volume;
} // double Particle::Get_Volume(void) const {
double Particle::Get_Radius(void) const {
  assert((*this).Radius_Set == true);
  return (*this).Radius;
} // double Particle::Get_Radius(void) const {

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

bool Particle::Get_Has_BC(unsigned Component) const {
  assert(Component < 3);
  return Has_BC[Component];
} // bool Particle::Get_Has_BC(unsigned Component) const {

double Particle::Get_BC(unsigned Component) const {
  assert(Component < 3);

  if(Has_BC[Component] == false) {
    char Buf[500];
    sprintf(Buf,
            "No BC Exception: Thrown by Particle::Get_BC\n"
            "You requested the %u component of the particle with ID %u. However, the\n"
            "Particle with ID %u has no BC in the %u component direction\n",
            Component,
            (*this).ID,
            (*this).ID,
            Component);
    throw No_BC(Buf);
  } // if(Has_BC[Component] == false) {

  return BC[Component];
} // double Particle::Get_BC(unsigned Component) const {





////////////////////////////////////////////////////////////////////////////////
// Printing functions

void Particle::Print(void) const {
  // Print basic particle parameters.
  printf("ID:       %u\n",(*this).ID);

  printf("Mass   (Set = %u):  %lf\n", (*this).Mass_Set,   (*this).Mass);
  printf("Volume (Set = %u):  %lf\n", (*this).Volume_Set, (*this).Volume);
  printf("Radius (Set = %u):  %lf\n", (*this).Radius_Set, (*this).Radius);

  printf("X:   "); (*this).X.Print();
  printf("x:   "); (*this).x.Print();
  printf("V:   "); (*this).V.Print();
  printf("a:   "); (*this).a.Print();

  printf("P:   \n");
  (*this).P.Print();
  printf("F[0]:   \n");
(*this).F[0].Print();
  printf("F[1]:   \n");
  (*this).F[1].Print();
  printf("A^(-1)\n");
  (*this).A_Inv.Print();

  printf("F_Int = "); (*this).Force_Int.Print();
  printf("F_Hg = "); (*this).Force_HG.Print();
  printf("F_Contact = "); (*this).Force_Contact.Print();
  printf("F_Friction = "); (*this).Force_Friction.Print();
  printf("\n");

  printf("Stretch_M:        %lf\n",(*this).Stretch_M);
  printf("Stretch_H:        %lf\n",(*this).Stretch_H);
  printf("Stretch_Critical: %lf\n",(*this).Stretch_Critical);
  printf("D:                %lf\n",(*this).D);

  #if defined(PARTICLE_DEBUG)
    printf("F_Visc = "); (*this).Force_Visc.Print();
    printf("Visc: \n"); (*this).Visc.Print();
  #endif

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
    printf("Grad_W address: %p\n",Grad_W);
    printf("|Grad_W|       : {");
    for(unsigned i = 0; i < (*this).Num_Neighbors-1; i++) {
      printf("%5.3f, ",Magnitude((*this).Grad_W[i]));
    } // for(unsigned i = 0; i < Num_Neighbors-1; i++) {
    printf("%5.3f } \n", Magnitude((*this).Grad_W[Num_Neighbors-1])); // */
  } // if(Neighbors_Are_Set == true) {

  printf("Has BC:    <%u %u %u>\n", (*this).Has_BC[0], (*this).Has_BC[1], (*this).Has_BC[2]);
  printf("BC:        <%6.3lf %6.3lf %6.3lf>\n", (*this).BC[0], (*this).BC[1], (*this).BC[2]);
} // void Particle::Print(void) const {



void Print(const Particle & P_In) { P_In.Print(); };         // Calls Print method
