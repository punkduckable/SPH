#if !defined(PARTICLE_HEADER)
#define PARTICLE_HEADER

class Particle {
  private:
    // Kernel paramaters
    static double h;                                       // Support radius             : mm
    static double Shape_Function_Amp;                      // Shape function Amplitude   : mm^-3

    // Strain energy function parameters
    static double Lame;                                    // Lame paramater             : Mpa
    static double mu0;                                     // Shear modulus              : Mpa

    // Viscosity paramaters
    static double mu;                                      // Viscosity                  : Mpa*s

    // Hourglass (Hg) correction parameters
    static double E;                                       // Hourglass stiffness        : Mpa
    static double alpha;                                   // Hg control parameter       : unitless

    // Particle dimensions (Mass, Vol, etc...)
    double Mass;                                           // Particle Mass              : g
    double Vol;                                            // Particle volume            : mm^3

    // Neighbor variables
    bool Has_Neighbors = false;                            // True if the particle has neighbors, false otherwise
    unsigned int Num_Neighbors;                            // Keeps track of number of Neighbors
    unsigned int *Neighbor_IDs;                            // Dynamic array that stores neighbor ID's (arry index's for Particle array in main file)

    Vector *R;                                             // Dynamic array that stores neighbor reference displacement
    double *Mag_R;                                         // Dynamic array that stores the magnitude of each neighbor's reference displacement

    double *W;                                             // Dynamic array that stores shape function value for each neighbor
    Vector *Grad_W;                                        // Dynamic array that stores Gradient of the Shape function at each neighbors
    //Vector *Grad_W_Tilde;                                // Dynamic array that stores Grad_W_Tilde
    Tensor A_Inv;                                          // Inverse of shape tensor

    // Particle dynamics variables
    Vector X{0,0,0};                                       // reference position         : mm
    Vector x{0,0,0};                                       // Particle's current position. x_i at start of iteration (x_i+1 at end)        :  mm
    Vector vel{0,0,0};                                     // Particle's velocity. v_i+1/2 at start of iteration v_i+3/2 at end (Leap Frog)    :  mm/s

    bool First_Iteration = true;                           // True if we're on first time step. Tells us to use Forward Euler to get initial velocity (leap frog)
    Tensor P{0,0,0,
             0,0,0,
             0,0,0};                                       // First Piola-Kirchhoff stress tensor  : Mpa
    Tensor F{1,0,0,
             0,1,0,
             0,0,1};                                       // deformation gradient       : unitless
    double Max_Stretch = 0;
    const double Critical_Stretch = 1.3;

    // Forces acting on the particle
    Vector Force_Int{0,0,0};                              // Internal Force vector       : N
    Vector Force_Ext{0,0,0};                              // External/body force         : N
    Vector Force_Hg{0,0,0};                               // Hour-glass force            : N

  public:
    // Constructors, destructor
    Particle(void);                                        // Default constructor
    Particle(const Particle & P_In);                       // Copy constructor (performs a deep copy)

    ~Particle(void);                                       // Destructor

    // Particle equality
    Particle & operator=(const Particle & P_In);           // Defines P1 = P2 (performs a deep copy)

    // Particle setup methods
    void Set_Mass(const double Mass_In) { Mass = Mass_In; }// Set particle's mass        : g
    void Set_Vol(const double Vol_In) { Vol = Vol_In; }    // Set particle's volume      : mm^3
    void Set_X(const Vector & X_In) { X = X_In; }          // Set ref position           : mm
    void Set_x(const Vector & x_In) { x = x_In; }          // Set spacial position       : mm
    void Set_vel(const Vector & vel_In) { vel = vel_In; }  // Set particle's velocity    : mm/s
    void Set_Neighbors(const unsigned int N,
                       const unsigned int * Neighbor_ID_Array,
                       const Particle * Particles);        // Set Neighbors

    // Friend functions
    friend void Update_P(Particle & P_In,
                         const Particle * Particles,
                         const double dt);                 // Updates P_In's Second Piola-Kirchhoff stress tensor
    friend void Update_Particle_Position(Particle & P_In,
                                         const Particle * Particles,
                                         const double dt); // Updates P_In's x (spacial position)
    friend bool Are_Neighbors(const Particle & P1,
                              const Particle & P2);        // Returns true P1 and P2 are neighbors, false otherwise


    // Temp friends... Should be removed!
    friend void Particle_Tests(void);
    friend void Generate_Neighbor_Lists_Box(const unsigned int Num_Particles, Particle * Particles,
                                    const unsigned int num_x, const unsigned int num_y, const unsigned int num_z,
                                    const unsigned int Support_Radius);
    friend void Particle_Debugger_File::Export_Pariticle_Properties(const unsigned int Num_Particles, const Particle * Particles);
    friend void VTK_File::Export_Pariticle_Positions(const unsigned int Num_Particles, const Particle * Particles);

  // Printing function
  void Print(void) const;                                  // Print's info about particle (mostly for testing)

};

void Print(const Particle & P_In);                         // Calls Print method

void Generate_Neighbor_Lists(const unsigned int Num_Particles,
                             Particle * Particles);        // Generate neighbor list for every particle in 'Partilces' array

void Generate_Neighbor_Lists_Box(const unsigned int Num_Particles, Particle * Particles,
                                const unsigned int num_x, const unsigned int num_y, const unsigned int num_z,
                                const unsigned int Support_Radius);  // Generate neighbor list for a 'box' or 'cuboid' geometry

#endif
