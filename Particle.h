#if !defined(PARTICLE_HEADER)
#define PARTICLE_HEADER

class Particle {
  private:
    // Cuboid parameters (should be deleted if not working with a rectangle)
    static double Inter_Particle_Spacing;                  //                            : mm
    unsigned int i,j,k;                                    // ijk position of particle in the cuboid

    // Kernel parameters
    static double h;                                       // Support radius             : mm
    static double Shape_Function_Amp;                      // Shape function Amplitude   : 1/mm^3

    // Strain energy function parameters
    static double Lame;                                    // Lame parameter             : Mpa
    static double mu0;                                     // Shear modulus              : Mpa

    // Viscosity parameters
    static double mu;                                      // Viscosity                  : Mpa*s

    // Hourglass (Hg) correction parameters
    static double E;                                       // Hourglass stiffness        : Mpa
    static double alpha;                                   // Hg control parameter       : unitless

    // Particle ID
    unsigned int ID;                                       // Particle ID in the Particle's array.

    // Particle dimensions (Mass, Vol, etc...)
    double Mass;                                           // Particle Mass              : g
    double Vol;                                            // Particle volume            : mm^3
    double Radius;                                         // Particle radius            : mm

    // Particle dynamics variables
    Vector X{0,0,0};                                       // reference position         : mm Vector
    Vector x{0,0,0};                                       // Particle's current position. x_i at start of iteration (x_i+1 at end)        :  mm Vector
    Vector V{0,0,0};                                       // Particle's velocity. v_i+1/2 at start of iteration v_i+3/2 at end (Leap Frog):  mm/s Vector

    unsigned int First_Time_Step = true;                   // True if we're on first time step. Tells us to use Forward Euler to get initial velocity (leap frog)
    Tensor P{0,0,0,                                        // First Piola-Kirchhoff stress tensor  : Mpa Tensor
             0,0,0,
             0,0,0};
    Tensor F{1,0,0,                                        // deformation gradient       : unitless Tensor
             0,1,0,
             0,0,1};

    // Forces acting on the particle
    Vector Force_Int{0,0,0};                               // Internal Force vector       : N Vector
    Vector Force_Contact{0,0,0};                           // Contact Force Vector        : N Vector
    Vector Force_HG{0,0,0};                                // Hour-glass force            : N Vector

    //Vector Force_Visc{0,0,0};                              // For debugging
    //Tensor Visc{0,0,0,                                     // For debugging
    //            0,0,0,
    //            0,0,0};

    // Damage parameters
    double Stretch_H = 0;                                  // Historical max stretch     : unitless
    double Stretch_M = 0;                                  // Current max stretch        : unitless
    double Stretch_Critical;                               // If Principle stretch excedes Critical then particle incures damage : unitless
    double D = 0;                                          // Damage parameter           : unitless
    static double Tau;                                     // Damage rate parameter (see eq 26)

    // Neighbor variables
    unsigned int Neighbors_Are_Set = false;                // True if the particle has neighbors, false otherwise
    unsigned int Num_Neighbors;                            // Keeps track of number of Neighbors
    unsigned int *Neighbor_IDs;                            // Dynamic array that stores neighbor ID's (arry index's for Particle array in main file)
    Vector *R;                                             // Dynamic array that stores neighbor reference displacement
    double *Mag_R;                                         // Dynamic array that stores the magnitude of each neighbor's reference displacement
    double *W;                                             // Dynamic array that stores shape function value for each neighbor
    Vector *Grad_W;                                        // Dynamic array that stores Gradient of the Shape function at each neighbors
    Tensor A_Inv;                                          // Inverse of shape tensor

  public:
    // Constructors, destructor
    Particle(void);                                        // Default constructor
    Particle(const Particle & P_In);                       // Copy constructor (performs a deep copy)

    ~Particle(void);                                       // Destructor

    // Particle equality
    Particle & operator=(const Particle & P_In);           // Defines P1 = P2 (performs a deep copy)

    // Particle setup methods
    void Set_Mass(const double Mass_In) { Mass = Mass_In; }                    //        : g
    void Set_Vol(const double Vol_In) { Vol = Vol_In; }                        //        : mm^3
    void Set_Radius(const double Radius_In) { Radius = Radius_In; }            //        : mm
    void Set_X(const Vector & X_In) { X = X_In; }          // Set ref position           : mm Vector
    void Set_x(const Vector & x_In) { x = x_In; }          // Set spacial position       : mm Vector
    void Set_V(const Vector & V_In) { V = V_In; }             // Set particle's velocity    : mm/s Vector
    void Set_Neighbors(const unsigned int N,
                       const unsigned int * Neighbor_ID_Array,
                       const Particle * Particles);        // Set Neighbors

    // Update P
    friend void Particle_Helpers::Update_P(Particle & P_In,
                                           Particle * Particles,
                                           const double dt);                   // Updates P_In's Second Piola-Kirchhoff stress tensor

    // Update x
    friend void Particle_Helpers::Update_x(Particle & P_In,
                                           const Particle * Particles,
                                           const double dt);                   // Updates P_In's spacial position

    // Neighbor methods
    friend bool Particle_Helpers::Are_Neighbors(const Particle & P1,
                                                const Particle & P2);          // Returns true P1 and P2 are neighbors, false otherwise

    friend void Particle_Helpers::Find_Neighbors(const unsigned int Num_Particles,
                                                 Particle * Particles);        // Generate neighbor list for every particle in 'Partilces' array

    friend void Particle_Helpers::Find_Neighbors_Box(Particle & P_In,
                                                     Particle * Particles);

    friend void Particle_Helpers::Remove_Neighbor(Particle & P_In,
                                                  const unsigned int Remove_Neighbor_ID,
                                                  const Particle * Particles);

    // Damage methods
    friend void Particle_Helpers::Remove_Damaged_Particle(Particle & P_In,
                                                          Particle * Particles);

    // Contact methods
    friend void Particle_Helpers::Contact(Particle * Body_A,
                                          const unsigned int Num_Particles_A,
                                          Particle * Body_B,
                                          const unsigned int Num_Particles_B,
                                          const double h_In);

    // Other friends
    friend void Particle_Tests(void);

    friend void Particle_Debugger::Export_Pariticle_Forces(const unsigned int Num_Particles,
                                                           const Particle * Particles);

    friend void Data_Dump::Print_Data_To_File(const Particle * Particles,
                                              const unsigned int Num_Particles);

    friend void Data_Dump::Print_Particle_To_File(const Particle & P_In,
                                                  FILE * File);

    friend Particle * Data_Dump::Load_Data_From_File(unsigned int & Num_Particles);

    friend void Data_Dump::Load_Particle_From_File(Particle & P_In,
                                                   FILE * File);

    friend void VTK_File::Export_Pariticle_Positions(const unsigned int Num_Particles,
                                                     const Particle * Particles);

  // Printing function
  void Print(void) const;                                  // Print's info about particle (mostly for testing)

};

void Print(const Particle & P_In);                         // Calls Print method

#endif
