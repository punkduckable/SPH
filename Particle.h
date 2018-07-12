#if !defined(PARTICLE_HEADER)
#define PARTICLE_HEADER

class Particle {
  private:
    // Cuboid paramaters (should be deleted if not working with a rectangle)
    const static double Inter_Particle_Spacing;            //                            : mm
    unsigned int ijk[3];                                   // ijk position of particle in the cuboid

    // Kernel paramaters
    const static double h;                                 // Support radius             : mm
    const static double Shape_Function_Amp;                // Shape function Amplitude   : mm^-3

    // Strain energy function parameters
    const static double Lame;                              // Lame paramater             : Mpa
    const static double mu0;                               // Shear modulus              : Mpa

    // Viscosity paramaters
    const static double mu;                                // Viscosity                  : Mpa*s

    // Hourglass (Hg) correction parameters
    const static double E;                                 // Hourglass stiffness        : Mpa
    const static double alpha;                             // Hg control parameter       : unitless

    // Particle dimensions (Mass, Vol, etc...)
    double Mass;                                           // Particle Mass              : g
    double Vol;                                            // Particle volume            : mm^3
    double Radius;                                         // Particle radius            : mm

    // Neighbor variables
    bool Has_Neighbors = false;                            // True if the particle has neighbors, false otherwise
    unsigned int Num_Neighbors;                            // Keeps track of number of Neighbors
    unsigned int *Neighbor_IDs;                            // Dynamic array that stores neighbor ID's (arry index's for Particle array in main file)
    unsigned int ID;                                       // Particle ID in the Particle's array.

    Vector *R;                                             // Dynamic array that stores neighbor reference displacement
    double *Mag_R;                                         // Dynamic array that stores the magnitude of each neighbor's reference displacement
    double *W;                                             // Dynamic array that stores shape function value for each neighbor
    Vector *Grad_W;                                        // Dynamic array that stores Gradient of the Shape function at each neighbors
    Tensor A_Inv;                                          // Inverse of shape tensor

    // Particle dynamics variables
    Vector X{0,0,0};                                       // reference position         : mm
    Vector x{0,0,0};                                       // Particle's current position. x_i at start of iteration (x_i+1 at end)        :  mm
    Vector vel{0,0,0};                                     // Particle's velocity. v_i+1/2 at start of iteration v_i+3/2 at end (Leap Frog):  mm/s

    bool First_Time_Step = true;                           // True if we're on first time step. Tells us to use Forward Euler to get initial velocity (leap frog)
    Tensor P{0,0,0,
             0,0,0,
             0,0,0};                                       // First Piola-Kirchhoff stress tensor  : Mpa
    Tensor F{1,0,0,
             0,1,0,
             0,0,1};                                       // deformation gradient       : unitless

    // Damage parameters
    double Stretch_H = 0;                                  // Historical max stretch     : unitless
    double Stretch_M = 0;                                  // Current max stretch        : unitless
    double Stretch_Critical;                               // If Principle stretch excedes Critical then particle incures damage : unitless
    double D = 0;                                          // Damage parameter           : unitless
    const static double Tau;                               // Damage rate paramater (see eq 26)

    // Forces acting on the particle
    Vector Force_Int{0,0,0};                               // Internal Force vector       : N
    Vector Force_Contact{0,0,0};                           // Contact Force Vector        : N
    Vector Force_Hg{0,0,0};                                // Hour-glass force            : N

    Vector Force_Visc{0,0,0};                              // For debugging
    Tensor Visc{0,0,0,
                0,0,0,
                0,0,0};                                    // For debugging

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
    void Set_X(const Vector & X_In) { X = X_In; }          // Set ref position           : mm
    void Set_x(const Vector & x_In) { x = x_In; }          // Set spacial position       : mm
    void Set_vel(const Vector & vel_In) { vel = vel_In; }  // Set particle's velocity    : mm/s
    void Set_Neighbors(const unsigned int N,
                       const unsigned int * Neighbor_ID_Array,
                       const Particle * Particles);        // Set Neighbors

    // Update P
    friend void Particle_Helpers::Update_P(Particle & P_In,
                         Particle * Particles,
                         const double dt);                 // Updates P_In's Second Piola-Kirchhoff stress tensor

    // Update x
    friend void Particle_Helpers::Update_x(Particle & P_In,
                         const Particle * Particles,
                         const double dt); // Updates P_In's spacial position

    // Neighbor methods
    friend bool Particle_Helpers::Are_Neighbors(const Particle & P1,
                              const Particle & P2);        // Returns true P1 and P2 are neighbors, false otherwise
    friend void Particle_Helpers::Find_Neighbors(const unsigned int Num_Particles,
                               Particle * Particles);      // Generate neighbor list for every particle in 'Partilces' array

    friend void Particle_Helpers::Find_Neighbors_Box(Particle & P_In, Particle * Particles);

    // Damage methods
    friend void Particle_Helpers::Remove_Damaged_Particle(Particle & P_In, Particle * Particles);

    // Contact methods
    friend void Particle_Helpers::Contact(Particle * Body_A, const unsigned int Num_Particles_A,
                                          Particle * Body_B, const unsigned int Num_Particles_B,
                                          const double h_In);

    // Other friends
    friend void Particle_Tests(void);
    friend void Particle_Debugger::Export_Pariticle_Properties(const unsigned int Num_Particles, const Particle * Particles);
    friend void VTK_File::Export_Pariticle_Positions(const unsigned int Num_Particles, const Particle * Particles);

  // Printing function
  void Print(void) const;                                  // Print's info about particle (mostly for testing)

};

void Print(const Particle & P_In);                         // Calls Print method

#endif
