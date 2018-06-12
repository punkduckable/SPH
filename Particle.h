#if !defined(_PARTICLE_HEADER)
#define _PARTICLE_HEADER

class Particle {
  private:
    static double h;                                    // Specifies connection radius (same for every particle)
    static double A;                                    // Scaling constant for spikey kernel
    static double k1;                                   // Constant in strain energy function
    static double k2;                                   // Constant in strain energy function
    static double mu0;                                  // Constant in strain energy function
    static double mu;                                   // Viscosity (in viscohyperelastic term)
    static double E;                                    // Hourglass stiffness
    static double alpha;                                // Coefficient to contols amplitue of hourclass correction
    static double rho;                                  // Particle density
    static Vector M;                                    // Fiber orientation vector

    double V;

    bool Neighbors_Set;
    unsigned int Num_Neighbors;                         // Keeps track of number of Neighbors
    unsigned int *Neighbor_List;                        // Dynamic array that stores neighbor ID's (arry index's for Particle array in main file)

    Vector *r;                                          // Stores neighbor spacial positions
    Vector *R;                                          // Stores neighbor reference positions

    Vector *Grad_W;                                     // Dynamic array that stores Grad_W
    Vector *Grad_W_Tilde;                               // Dynamic array that stores Grad_W_Tilde
    double Calc_W(const Vector & Rj);                   // Returns value of W (kernel function) for a given displacement vector

    Vector x;                                           // Particle's current position
    Vector X;                                           // reference position
    bool First_Step = true;
    Vector v;                                           // Particle's velocity at half time step (Leap-Frog method)
    Tensor P;                                           // First Piola-Kirchhoff stress tensor

    Tensor F;                                           // deformation gradient (this needs to be a member so that we can calculate F')

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
