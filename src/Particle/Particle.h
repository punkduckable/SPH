#if !defined(PARTICLE_HEADER)
#define PARTICLE_HEADER

#include "Classes.h"
#include "Namespaces.h"
#include "Body/Body.h"
#include "IO/Data_Dump.h"
#include "Simulation/Simulation.h"
#include "Particle/Particle_Helpers.h"
#include "Particle_Debugger.h"
#include "Vector/Vector.h"
#include "Tensor/Tensor.h"
#include "Errors.h"

class Particle {
  private:
    // Cuboid parameters (should be deleted if not working with a rectangle)
    unsigned i,j,k;                                        // ijk position of particle in the cuboid

    // Particle ID
    unsigned ID;                                           // Particle ID in the Particle's array.

    // Particle dimensions (Mass, Vol, etc...)
    double Mass;                                           // Particle Mass              : g
    double Vol;                                            // Particle volume            : mm^3
    double Radius;                                         // Particle radius            : mm

    // Particle dynamics variables
    Vector X{0,0,0};                                       // reference position         : mm Vector
    Vector x{0,0,0};                                       // Particle's current position. x_i at start of iteration (x_i+1 at end)        : mm Vector
    Vector V{0,0,0};                                       // Particle's velocity. v_i+1/2 at start of iteration v_i+3/2 at end (Leap Frog): mm/s Vector
    Vector a{0,0,0};                                       // Particle's acceleration. a_i at start of iteration (a_i+1 at end)            : mm/s^2 Vector

    Tensor P{0,0,0,                                        // First Piola-Kirchhoff stress tensor  : Mpa Tensor
             0,0,0,
             0,0,0};
    Tensor F[2] = {{1,0,0,                                 // deformation gradients      : unitless Tensor
                    0,1,0,
                    0,0,1},

                   {1,0,0,
                    0,1,0,
                    0,0,1}};

    // Forces acting on the particle
    Vector Force_Int{0,0,0};                               // Internal Force vector      : N Vector
    Vector Force_Contact{0,0,0};                           // Contact Force Vector       : N Vector
    Vector Force_Friction{0,0,0};                          // Frictional force Vector (from contact)    : N Vector
    Vector Force_HG{0,0,0};                                // Hour-glass force           : N Vector

    #if defined(PARTICLE_DEBUG)
      Vector Force_Visc{0,0,0};
      Tensor Visc{0,0,0,
                  0,0,0,
                  0,0,0};
    #endif

    // Damage parameters
    double Stretch_H = 0;                                  // Historical max stretch     : unitless
    double Stretch_M = 0;                                  // Current max stretch        : unitless
    double Stretch_Critical;                               // If Principle stretch excedes Critical then particle incures damage : unitless
    double D = 0;                                          // Damage parameter           : unitless

    // Neighbor variables
    unsigned Neighbors_Are_Set = false;                    // True if the particle has neighbors, false otherwise
    unsigned Num_Neighbors;                                // Keeps track of number of Neighbors
    unsigned *Neighbor_IDs;                                // Dynamic array. Stores neighbor ID's
    Vector *R;                                             // Dynamic array. Stores neighbor reference displacement                        : mm
    double *Mag_R;                                         // Dynamic array. Stores magnitude of reference displacement to each neighbor   : mm
    double *W;                                             // Dynamic array. Stores shape function value for each neighbor                 : 1/(mm^3)
    Vector *Grad_W;                                        // Dynamic array. Stores Gradient of the Shape function at each neighbor        : 1/(mm^4)
    Tensor A_Inv;                                          // Inverse of shape tensor                                                      : unitless

  public:
    // Constructors, destructor
    Particle(void);                                        // Default constructor
    Particle(const Particle & P_In) = delete;              // Copy constructor

    ~Particle(void);                                       // Destructor

    // Particle equality
    Particle & operator=(const Particle & P_In) = delete;  // Defines P1 = P2

    // Methods to set particle properties
    void Set_ijk(const unsigned i_in, const unsigned j_in, const unsigned k_in) {
      i = i_in;
      j = j_in;
      k = k_in;
    } // void Set_ijk(const unsigned i_in, const unsigned j_in, const unsigned k_in) {
    void Set_ID(const unsigned ID_In) { ID = ID_In; }
    void Set_Mass(const double Mass_In) { Mass = Mass_In; }                    //        : g
    void Set_Vol(const double Vol_In) { Vol = Vol_In; }                        //        : mm^3
    void Set_Radius(const double Radius_In) { Radius = Radius_In; }            //        : mm
    void Set_X(const Vector & X_In) { X = X_In; }          // Set ref position           : mm Vector
    void Set_x(const Vector & x_In) { x = x_In; }          // Set spacial position       : mm Vector
    void Set_V(const Vector & V_In) { V = V_In; }          // Set particle's velocity    : mm/s Vector
    void Set_a(const Vector & a_In) { a = a_In; }          // Set particle's acel        : mm/s^2 Vector
    void Set_D(const double D_In) { D = D_In; }            // Set Damage parameter       : unitless


    // Methods to get particle properties
    unsigned Get_i(void) const { return i; }
    unsigned Get_j(void) const { return j; }
    unsigned Get_k(void) const { return k; }
    unsigned Get_ID(void) const { return ID; }
    double Get_Mass(void) const { return Mass; }
    double Get_Vol(void) const { return Vol; }
    double Get_Radius(void) const { return Radius; }
    const Vector & Get_X(void) const { return X; }
    const Vector & Get_x(void) const { return x; }
    const Vector & Get_V(void) const { return V; }
    const Vector & Get_a(void) const { return a; }
    const Tensor & Get_P(void) const { return P; }
    const Tensor & Get_F(const unsigned i) const { return F[i]; }
    const Vector & Get_Force_Friction(void) const { return Force_Friction; }
    const Vector & Get_Force_Contact(void) const { return Force_Contact; }
    double Get_Stretch_M(void) const { return Stretch_M; }
    double Get_Stretch_H(void) const { return Stretch_H; }
    double Get_Stretch_Critical(void) const { return Stretch_Critical; }
    double Get_D(void) const { return D; }
    unsigned Get_Num_Neighbors(void) const { return Num_Neighbors; }
    unsigned Get_Neighbor_IDs(unsigned i) const {
      if(i < Num_Neighbors)
        return Neighbor_IDs[i];

      printf("Requested neighbor ID is out of bounds! Num_Neighbors = %u, requested index = %u\n",Num_Neighbors, i);
      return 0;
    } // unsigned Get_Neighbor_IDs(unsigned i) const {

      // body method
    // Update P
    friend void Particle_Helpers::Update_P(Body & Particles,                   // Updates Second Piola-Kirchhoff stress tensor for each particle in a Body
                                           const double dt);

       // body method
    // Update x
    friend void Particle_Helpers::Update_x(Body & Particles,                   // Updates spatial position for each particle in a Body
                                           const double dt);

    // Neighbor friends (other neighbor methods are in Particle_Neighbors.cc)
    void Set_Neighbors(const unsigned N,                                       // Set Neighbors
                       const unsigned * Neighbor_ID_Array,
                       const Body & Particles);

      // body method
    friend bool Particle_Helpers::Are_Neighbors(const Body & Particles,
                                                const unsigned i,
                                                const unsigned j);

      // body method
    friend void Particle_Helpers::Find_Neighbors_Cuboid(Particle & P_In,
                                                        Body & Particles);

      // body method
    friend void Particle_Helpers::Remove_Neighbor(Particle & P_In,
                                                  const unsigned Remove_Neighbor_ID,
                                                  const Body & Particles);

       // body method
    // Damage method friends
    friend void Particle_Helpers::Remove_Damaged_Particle(Particle & P_In,
                                                          Body & Particles);

    // Contact method friends
    friend void Particle_Helpers::Contact(Body & Body_A,
                                          Body & Body_B);

    // Other friends
    friend void Particle_Tests(void);
    friend void Simulation::Run_Simulation(void);
      // body method
    friend void Particle_Debugger::Export_Particle_Forces(const Body & Particles);
      // body method
    friend int Data_Dump::Load_Body(Body & Particles);
    friend void Data_Dump::Load_Particle(Particle & P_In,
                                         FILE * File,
                                         const bool Is_Cuboid);

    // Printing function
    void Print(void) const;                                  // Print's info about particle (mostly for testing)
}; // class Particle {

void Print(const Particle & P_In) { P_In.Print(); };         // Calls Print method

#endif
