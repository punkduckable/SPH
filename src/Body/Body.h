#if !defined(BODY_HEADER)
#define BODY_HEADER

#define DAMAGE_MONITOR

#include "Materials.h"
#include "Classes.h"
#include "IO/Data_Dump.h"
#include <string>

#if !defined(PI)
  #define PI 3.1415926535897932384626
#endif

class Body {
  private:
    ////////////////////////////////////////////////////////////////////////////
    // Body type parameters

    // Cuboid parameters.
    bool Is_Cuboid = false;
    unsigned X_SIDE_LENGTH = 0;
    unsigned Y_SIDE_LENGTH = 0;
    unsigned Z_SIDE_LENGTH = 0;

    // Boundary parameters
    bool Is_Boundary = false;


    ////////////////////////////////////////////////////////////////////////////
    // Body parameters

    // Particle array name
    std::string Name;

    // The array of particles
    bool Particles_Set_Up = false;
    unsigned Num_Particles = 0;
    Particle * Particles = nullptr;

    // Kernel parameters
    double Inter_Particle_Spacing = 0;                     //                            : mm
    unsigned Support_Radius = 0;                           // Support radius in Inter Particle Spacings's    : unitless
    double h = 0;                                          // Support radius in mm's     : mm
    double Shape_Function_Amplitude = 0;                   // Shape function Amplitude   : 1/(mm^6)

    // Material Parameters
    Materials::Material Body_Material;

    // Viscosity parameters
    double mu;                                             // Viscosity                  : Mpa*s
    unsigned char F_Index = 0;                             // Keeps track of which F (in each partilce) is the 'new' one

    // Hourglass (Hg) correction parameters
    double alpha;                                          // Hg control parameter       : unitless

    // Damage paramaters
    double Tau;                                            // Damage rate parameter (see eq 26)
    bool Damageable = 1;                                   // If true, allows this Body to take damage.

    // First time step?
    bool First_Time_Step = true;



    ////////////////////////////////////////////////////////////////////////////
    // Contact parameters
    static double K;                             // Set in Body.h                        : N/(mm^2)



    ////////////////////////////////////////////////////////////////////////////
    // Printing Parameters
    unsigned Times_Printed_Net_External_Force = 0;
    unsigned Times_Printed_Particle_Forces = 0;



    ////////////////////////////////////////////////////////////////////////////
    // Private methods
    void Set_h(const double h_In);

  public:
    ////////////////////////////////////////////////////////////////////////////
    // Constructors, destructor
    Body(void);                                  // default constructor
    Body(const unsigned Num_Particles_In);       // generate array constructor
    Body(const Body & Ar_In) = delete;           // deleted copy constructor

    ~Body(void);                                 // destructor



    ////////////////////////////////////////////////////////////////////////////
    // Operator overloading
    Body & operator=(Body & Ar_In) = delete;     // deleted = operator
    Particle & operator[](const unsigned i);
    const Particle & operator[](const unsigned i) const;



    ////////////////////////////////////////////////////////////////////////////
    // Neighbor methods.
    // These functions are defined in Neighbors.cc
    void Set_Neighbors(const unsigned i,         // Set Neighbors
                       const unsigned Num_Neighbors_In,
                       const unsigned * Neighbor_ID_Array);
    void Find_Neighbors(void);                   // Generate neighbor list for every particle in 'Particles' array
    void Find_Neighbors_Cuboid(void);            // Generate neighbor list for 'cuboid' geometry
    bool Are_Neighbors(const unsigned i,         // Returns true if particles i,j are neighbors
                       const unsigned j) const;
    void Remove_Neighbor(const unsigned i,       // Removes one of particle i's neighbors
                         const unsigned Remove_Neighbor_ID);



    //////////////////////////////////////////////////////////////////////////////
    // Update methods.
    // These function are defined in Update.cc
    void Update_P(const double dt);             // Update Second Piola-Kirchhoff stress tensor of each particle in a Body
    void Update_x(const double dt);             // Updates spacial position of each particle in a Body



    //////////////////////////////////////////////////////////////////////////////
    // Damage methods.
    // This function is defined in Damage.cc
    void Remove_Damaged_Particle(const unsigned i);



    //////////////////////////////////////////////////////////////////////////////
    // Contact methods.
    // This Function is defined in Contact.cc
    static void Contact(Body & Body_A,
                        Body & Body_B);



    ////////////////////////////////////////////////////////////////////////////
    // Setters
    void Set_Num_Particles(const unsigned Num_Particles_In);
    void Set_Name(const std::string & S_In);

    void Set_Inter_Particle_Spacing(const double IPS);
    void Set_Support_Radius(const unsigned SR_In);

    void Set_Material(const Materials::Material & Mat_In);
    void Set_mu(const double mu_In);
    void Set_alpha(const double alpha_In);

    void Set_Tau(const double Tau_In);
    void Set_Damageable(const bool D_In);

    void Set_Cuboid_Dimensions(const Vector & Dimensions);

    void Set_Boundary(const bool Boundary_In);

    void Set_First_Time_Step(const bool First_In);

    void Set_F_Index(const unsigned char i);
    void Increment_F_Index(void);



    ////////////////////////////////////////////////////////////////////////////
    // Getters
    unsigned Get_Num_Particles(void) const;
    std::string Get_Name(void) const;

    double Get_Inter_Particle_Spacing(void) const;
    unsigned Get_Support_Radius(void) const;

    double Get_h(void) const;
    double Get_Shape_Function_Amplitude(void) const;
    Materials::Material Get_Material(void) const;
    double Get_Lame(void) const;
    double Get_mu0(void) const;
    double Get_mu(void) const;
    double Get_E(void) const;
    double Get_density(void) const;
    double Get_alpha(void) const;

    unsigned char Get_F_Index(void) const;

    double Get_Tau(void) const;
    bool Get_Damagable(void) const;

    bool Get_Cuboid(void) const;
    unsigned Get_X_SIDE_LENGTH(void) const;
    unsigned Get_Y_SIDE_LENGTH(void) const;
    unsigned Get_Z_SIDE_LENGTH(void) const;

    bool Get_Boundary(void) const;

    bool Get_First_Time_Step(void) const;



    ////////////////////////////////////////////////////////////////////////////
    // Printing methods
    void Print_Net_External_Force(const unsigned time_step);
    void Print_Parameters(void) const;
    void Print_Particle_Forces(void);



    ////////////////////////////////////////////////////////////////////////////
    // Friends
    friend void Data_Dump::Save_Simulation(const Body * Arrays,
                                           const unsigned Num_Arrays);
    friend int Data_Dump::Load_Simulation(Body ** Array_Ptr,
                                          unsigned & Num_Bodies);
}; // class Body {

#endif
