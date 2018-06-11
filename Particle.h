#if !defined(_PARTICLE_HEADER)
#define _PARTICLE_HEADER

class Particle {
  private:
    double V = 1;
    static double h;                                    // Specifies connection radius (same for every particle)

    unsigned int Num_Neighbors;                         // Keeps track of number of Neighbors
    unsigned int *Neighbor_List;                        // Dynamic array that stores neighbor ID's (arry index's for Particle array in main file)
    Vector *Grad_W_tilde;                               // Dynamic array that stores Grad_W_Tilde for each particle

    Vector x;                                           // Particle's current position
    Vector x_new;                                       // Particle's position after current time step
    Vector X;                                           // reference position

    Vector Calc_Grad_W_Tilde(const Vector Rj) const;    // Calculate Grad_W_Tilde: Grad_W_Tilde = A^(-1)*Grad_W
    Vector Calc_Grad_W(const Vector Rj) const;          // Calculates Grad_W for a given Rj

  public:
    Particle(void);                                     // Default constructor
    Particle(const Vector & X_In);                      // Vector constructor (accept ref position)
    Particle(const Particle & P_In);                    // Copy constructor (performs a deep copy)

    ~Particle(void);                                    // Destructor

    void Set_X(const Vector X_In);                      // Set ref position
    Vector Get_X(void) const;                           // Returns ref position
    Vector Get_x(void) const;                           // Returns current position

    void Set_Neighbors(const unsigned int N,
                       const unsigned int *Neighbor_List_In,
                       const Vector *X_Neighbors);      // Set Neighbors

    Particle & operator=(const Particle & P_In);        // Defines P1 = P2 (performs a deep copy)

    // Friend functions
    friend void Update_Particle_Position(const Particle & P_In,
                         const Particle *Particles,
                         const double dt);
};

#endif
