#if !defined(_PARTICLE_HEADER)
#define _PARTICLE_HEADER

class Particle {
  private:
    double V = 1;
    static double h;                                    // Specifies connection radius (same for every particle)

    unsigned int Num_Neighbors;                         // Keeps track of number of Neighbors
    unsigned int *Neighbor_List;                        // Dynamic array that stores neighbor ID's (arry index's for Particle array in main file)
    Vector *Grad_W_tilde;                               // Dynamic array that stores Grad_W_Tilde for each particle

    Vector x;                                           // current position
    Vector X;                                           // reference position
    Vector Force;                                       // Force vector
    Vector acceleration;                                // acceleration vector

    Tensor F;                                           // Deformation gradient
    Tensor A_Inv;                                       // Inverse of shape tensor
    Tensor P;                                           // First Poila-Kirchoff stress tensor
    Tensor S;                                           // Second Poila-Kirchoff stress tensor

    Vector Calc_Grad_W(const Vector Rj) const;          // Calculates Grad_W for a given Rj
    void Calc_A_Inv(void);                              // Calculates A_Inverse (can only be called once Neighbors are set)
    Vector Calc_Grad_W_Tilde(const Vector Rj) const;    // Calculate Grad_W_Tilde: Grad_W_Tilde = A^(-1)*Grad_W

    void Update_F(void);
    void Update_P(void);
    void Update_Force(void);
    void Update_acceleration(void);

  public:
    Particle(void);                                     // Default constructor
    Particle(const Vector & X_In);                      // Vector constructor (accept ref position)
    Particle(const Particle & P_In);                    // Copy constructor

    ~Particle(void);                                    // Destructor

    void Set_X(const Vector X_In);                      // Set ref position
    Vector Get_X(void) const;                           // Returns ref position
    Vector Get_x(void) const;                           // Returns current position

    void Set_Neighbors(const unsigned int N,
                       const unsigned int *Neighbor_List_In,
                       const Vector *X_Neighbors);      // Set Neighbors

    void Update_x(const Vector *x_Neighbors, const double dt);
};

#endif
