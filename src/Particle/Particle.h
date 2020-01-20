#if !defined(PARTICLE_HEADER)
#define PARTICLE_HEADER

#include "Classes.h"
#include "Vector/Vector.h"
#include "Tensor/Tensor.h"
#include "Simulation/Simulation.h"
#include "IO/Save_Simulation.h"
#include "IO/Load_Simulation.h"


class Particle {
  friend Body;

  private:
    unsigned ID;

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

    // BC variables.
    bool Has_BC[3];                                        // ith component is true if this particle has a BC in the ith direction.
    double BC[3];                                          // velocity BC       : mm/s Vector

  public:
    ////////////////////////////////////////////////////////////////////////////
    // Constructors, destructor
    Particle(void);                                        // Default constructor
    Particle(const Particle & P_In) = delete;              // Copy constructor

    ~Particle(void);                                       // Destructor



    ////////////////////////////////////////////////////////////////////////////
    // Operator Overloading

    Particle & operator=(const Particle & P_In) = delete;  // deleted = operator



    ////////////////////////////////////////////////////////////////////////////
    // Other functions

    void Apply_BCs(void);                                  // Applies boundary conditions




    ////////////////////////////////////////////////////////////////////////////
    // Setters
    void Set_ID(const unsigned ID_In);
    void Set_Mass(const double Mass_In);                                       //        : g
    void Set_Vol(const double Vol_In);                                         //        : mm^3
    void Set_Radius(const double Radius_In);                                   //        : mm

    void Set_X(const Vector & X_In);                       // Set ref position           : mm Vector
    void Set_x(const Vector & x_In);                       // Set spacial position       : mm Vector
    void Set_V(const Vector & V_In);                       // Set particle's velocity    : mm/s Vector
    void Set_a(const Vector & a_In);                       // Set particle's acel        : mm/s^2 Vector

    void Set_D(const double D_In);                         // Set Damage parameter       : unitless

    void Set_BC(const unsigned Component,                  // Set a component of the particles BC
                const double Value);


    ////////////////////////////////////////////////////////////////////////////
    // Getters
    unsigned Get_ID(void) const;                                               //
    double Get_Mass(void) const;                                               //        : g
    double Get_Vol(void) const;                                                //        : mm^3
    double Get_Radius(void) const;                                             //        : mm

    const Vector & Get_X(void) const;                                          //        : mm Vector
    const Vector & Get_x(void) const;                                          //        : mm Vector
    const Vector & Get_V(void) const;                                          //        : mm/s Vector
    const Vector & Get_a(void) const;                                          //        : mm/s^2 Vector
    const Tensor & Get_P(void) const;                                          //        : Mpa Tensor
    const Tensor & Get_F(const unsigned i) const;                              //        : unitless Tensor

    const Vector & Get_Force_Friction(void) const;                             //        : N Vector
    const Vector & Get_Force_Contact(void) const;                              //        : N Vector

    double Get_Stretch_M(void) const;                                          //        : unitless
    double Get_Stretch_H(void) const;                                          //        : unitless
    double Get_Stretch_Critical(void) const;                                   //        : unitless
    double Get_D(void) const;                                                  //        : unitless

    unsigned Get_Num_Neighbors(void) const;
    unsigned Get_Neighbor_IDs(unsigned i) const;

    bool Get_Has_BC(unsigned Component) const;             // true if particle has a BC in the direction of specified Component
    double Get_BC(unsigned Component) const;               // returns BC in the Component direction (if the particle has one)



  ////////////////////////////////////////////////////////////////////////////
  // Friends, Printing
  friend void Particle_Tests(void);
  friend void Simulation::Run_Simulation(void);
  friend void Simulation::Apply_Box_Particle_BCs(Particle & P_In, Vector BC);
  friend void Simulation::Apply_General_BCs(Body & Body_In, Array<Boundary_Condition> & BCs_In);
  friend void IO::Load_Body(Body & Body_In);
  friend void IO::Load_Particle(Particle & P_In,
                                  FILE * File);

  // Printing function
  void Print(void) const;                                  // Print's info about particle (mostly for testing)
}; // class Particle {

void Print(const Particle & P_In);

#endif
