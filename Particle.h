#if !defined(_PARTICLE_HEADER)
#define _PARTICLE_HEADER

class Particle {
  private:
    // Kernel paramaters
    static double h;                                    // Specifies connection radius (same for every particle)
    static double A;                                    // Scaling constant for spikey kernel

    // Strain energy function parameters
    static double k1;                                   // Constant in strain energy function
    static double k2;                                   // Constant in strain energy function
    static double mu0;                                  // Constant in strain energy function

    // Viscosity paramaters
    static double mu;                                   // Viscosity (in viscohyperelastic term)

    // Hourglass correction parameters
    static double E;                                    // Hourglass stiffness
    static double alpha;                                // Coefficient to contols amplitue of hourclass correction

    // Mass parameters
    static double density;                              // Particle density
    static Vector M;                                    // Fiber orientation vector

    double V;                                           // Particle volume
    double M;                                           // Particle Mass

    bool Has_Neighbors = false;                         // True if the particle has neighbors, false otherwise
    unsigned int Num_Neighbors;                         // Keeps track of number of Neighbors
    unsigned int *Neighbor_List;                        // Dynamic array that stores neighbor ID's (arry index's for Particle array in main file)

    Vector *r;                                          // Stores neighbor spacial positions
    Vector *R;                                          // Stores neighbor reference positions

    Vector *Grad_W;                                     // Dynamic array that stores Grad_W
    Vector *Grad_W_Tilde;                               // Dynamic array that stores Grad_W_Tilde
    double Calc_W(const Vector & Rj);                   // Returns value of W (kernel function) for a given displacement vector

    Vector x;                                           // Particle's current position: x_i at start of iteration (x_i+1 at end)
    Vector v;                                           // Particle's velocity at half time step (Leap-Frog method): v_i+1/2 at start of iteration (v_i+3/2 at end)
    Vector X;                                           // reference position

    bool First_Iteration = true;                        // True if we're on first time step. Tells us to use Forward Euler to get initial velocity (leap frog)
    Tensor P;                                           // First Piola-Kirchhoff stress tensor (needed to update position)
    Tensor F;                                           // deformation gradient (needed to calculate (d/dt)F)

  public:
    Particle(void);                                     // Default constructor
    Particle(const Particle & P_In);                    // Copy constructor (performs a deep copy)

    ~Particle(void);                                    // Destructor

    void Set_X(const Vector & X_In);                    // Set ref position
    void Set_V(const double V_In);
    void Set_Neighbors(const unsigned int N,
                      const unsigned int *Neighbor_Id_List,
                      const Particle *Particles);       // Set Neighbors

    void Set_P(const Tensor & P_In);
    void Set_F(const Tensor & F_In);

    Particle & operator=(const Particle & P_In);        // Defines P1 = P2 (performs a deep copy)

    // Friend functions
    friend void Update_P(Particle & P_In,
                         const Particle *Particles,
                         const double dt);

   friend void Update_Particle_Position(Particle & P_In,
                                        const Particle *Particles,
                                        const double dt);
};

#endif
