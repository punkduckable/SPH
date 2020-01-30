#if !defined(BODY_HEADER)
#define BODY_HEADER

#define DAMAGE_MONITOR
#define IO_MONITOR

#include "Materials.h"
#include "Classes.h"
#include "IO/Save_Simulation.h"
#include "IO/Load_Simulation.h"
#include <string>

#if !defined(PI)
  #define PI 3.1415926535897932384626
#endif

class Body {
  private:
    ////////////////////////////////////////////////////////////////////////////
    // Body type parameters

    // Box parameters.
    bool Is_Box = false;
    unsigned X_SIDE_LENGTH = 0;
    unsigned Y_SIDE_LENGTH = 0;
    unsigned Z_SIDE_LENGTH = 0;

    // Is Fixed parameters
    bool Is_Fixed = false;


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
    bool Is_Damageable = 1;                                // If true, allows this Body to take damage.

    // Time step information
    bool First_Time_Step = true;
    unsigned Time_Steps_Per_Update;                   // The number of time steps between when P and F are updated



    ////////////////////////////////////////////////////////////////////////////
    // Contact parameters
    static double K;                             // Set in Body.h                        : N/(mm^2)



    ////////////////////////////////////////////////////////////////////////////
    // Printing Parameters
    unsigned Times_Printed_Net_External_Force = 0;
    unsigned Times_Printed_Particle_Forces = 0;
    unsigned Times_Printed_Particle_Positions = 0;



    ////////////////////////////////////////////////////////////////////////////
    // Private methods
    void Set_h(const double h_In);               // Defined in Body.cc

    void Add_Point_Data(FILE * File,             // Defined in Body_IO.cc
                        char * Weight_Name,      // Helps Export_Particle_Positions.
                        unsigned Num_Particles,
                        double * Data) const;



  public:
    ////////////////////////////////////////////////////////////////////////////
    // Constructors, destructor
    // Defined in Body.cc
    Body(void);                                  // default constructor
    Body(const unsigned Num_Particles_In);       // generate array constructor
    Body(const Body & Ar_In) = delete;           // deleted copy constructor

    ~Body(void);                                 // destructor



    ////////////////////////////////////////////////////////////////////////////
    // Operator overloading
    // Defined in Body.cc
    Body & operator=(Body & Ar_In) = delete;     // deleted = operator
    Particle & operator[](const unsigned i);
    const Particle & operator[](const unsigned i) const;



    ////////////////////////////////////////////////////////////////////////////
    // Boundary Conditions
    // Defined in Body.cc

    void Apply_BCs(void);                        // Applies BCs for each particle in the Body



    ////////////////////////////////////////////////////////////////////////////
    // Neighbor methods.
    // Defined in Neighbors.cc
    void Set_Neighbors(const unsigned i,         // Set Neighbors + Neighbor dependent members (using Set_Neighbor_Dependent_Members)
                       const Array<unsigned> & Neighbor_IDs_In);
    void Set_Neighbor_Dependent_Members(const unsigned i); // Sets Particle[i]'s members that can't be set up without a Neighbor_List
    void Find_Neighbors(void);                   // Generate neighbor list for every particle in 'Particles' array
    void Find_Neighbors_Box(void);               // Generate neighbor list for 'Box' geometry
    bool Are_Neighbors(const unsigned i,         // Returns true if particles i,j are neighbors
                       const unsigned j) const;
    void Remove_Neighbor(const unsigned i,       // Removes one of particle i's neighbors
                         const unsigned Remove_Neighbor_ID);



    //////////////////////////////////////////////////////////////////////////////
    // Update methods.
    // Defined in Update.cc
    void Update_P(const double dt);             // Update Second Piola-Kirchhoff stress tensor of each particle in a Body
    void Update_x(const double dt);             // Updates spacial position of each particle in a Body



    //////////////////////////////////////////////////////////////////////////////
    // Damage methods.
    // Defined in Damage.cc
    void Remove_Damaged_Particle(const unsigned i);



    //////////////////////////////////////////////////////////////////////////////
    // Contact methods.
    // Defined in Contact.cc
    static void Contact(Body & Body_A,
                        Body & Body_B);



    ////////////////////////////////////////////////////////////////////////////
    // Setters
    // Defined in Body.cc
    void Set_Num_Particles(const unsigned Num_Particles_In);
    void Set_Name(const std::string & S_In);

    void Set_Inter_Particle_Spacing(const double IPS);
    void Set_Support_Radius(const unsigned SR_In);

    void Set_Material(const Materials::Material & Mat_In);
    void Set_mu(const double mu_In);
    void Set_alpha(const double alpha_In);

    void Set_Tau(const double Tau_In);
    void Set_Is_Damageable(const bool D_In);

    void Set_Box_Dimensions(const unsigned Dim_x,
                            const unsigned Dim_y,
                            const unsigned Dim_z);

    void Set_Is_Fixed(const bool Is_Fixed_In);

    void Set_First_Time_Step(const bool First_In);
    void Set_Time_Steps_Per_Update(const unsigned Steps_In);

    void Set_F_Index(const unsigned char i);
    void Increment_F_Index(void);



    ////////////////////////////////////////////////////////////////////////////
    // Getters
    // Defined in Body.cc
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
    bool Get_Is_Damageable(void) const;

    bool Get_Is_Box(void) const;
    unsigned Get_X_SIDE_LENGTH(void) const;
    unsigned Get_Y_SIDE_LENGTH(void) const;
    unsigned Get_Z_SIDE_LENGTH(void) const;

    bool Get_Is_Fixed(void) const;

    bool Get_First_Time_Step(void) const;
    unsigned Get_Time_Steps_Per_Update(void) const;



    ////////////////////////////////////////////////////////////////////////////
    // Printing methods
    // Defined in Body_IO.cc
    void Print_Parameters(void) const;                               // Prints to command line
    void Export_Net_External_Force(const unsigned time_step);        // Prints to file
    void Export_Particle_Forces(void);                               // Prints to file
    void Export_Particle_Positions(void);                            // Prints to file


    ////////////////////////////////////////////////////////////////////////////
    // Friends
    friend void IO::Save_Simulation(const Body * Arrays,
                                    const unsigned Num_Arrays);
    friend void IO::Load_Simulation(Body ** Array_Ptr,
                                    unsigned & Num_Bodies);
}; // class Body {

#endif
