#if !defined(_PARTICLE_HEADER)
#define _PARTICLE_HEADER

class Particle {
  private:
    static double h;                                    // Specifies connection radius (same for every particle)
    static double k1;
    static double k2;
    static double mu0;
    static double mu;
    static double E;
    static double alpha;
    static double rho;
    static Vector M;

    double V;

    bool Neighbors_Set;
    unsigned int Num_Neighbors;                         // Keeps track of number of Neighbors
    unsigned int *Neighbor_List;                        // Dynamic array that stores neighbor ID's (arry index's for Particle array in main file)
    Vector *r;                                           // Stores neighbor spacial positions
    Vector *R;                                           // Stores neighbor reference positions
    Vector *Grad_W;                                     // Dynamic array that stores Grad_W
    Vector *Grad_W_Tilde;                               // Dynamic array that stores Grad_W_Tilde

    Vector x;                                           // Particle's current position
    Vector x_new;                                       // Particle's position after current time step.
    Vector X;                                           // reference position
    Tensor P;                                           // First Piola-Kirchhoff stress tensor

    Tensor F;                                           // deformation gradient (needed to calculate F')

    double Calc_W(const Vector & Rj);

  public:
    Particle(void);                                     // Default constructor
    Particle(const Particle & P_In);                    // Copy constructor (performs a deep copy)

    ~Particle(void);                                    // Destructor

    void Set_X(const Vector & X_In);                      // Set ref position
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
