#if !defined(_PARTICLE_HEADER)
#define _PARTICLE_HEADER

class Particle {
  private:
    double V = 1;

    unsigned int *Neighbor_List;
    unsigned int Num_Neighbors;

    Vector x;
    Vector X;
    Vector Force;
    Vector acceleration;

    Tensor F;
    Tensor A_Inv;
    Tensor P;
    Tensor S;

    Vector Grad_W(const Vector Rj, const double h) const;
    void Calculate_A_Inv(const double h);
    Vector Grad_W_Tilde(const Vector Rj, const double h) const;
    void Update_F(void);
    void Update_P(void);
    void Update_Force(void);
    void Update_acceleration(void);

  public:
    Particle(void);
    Particle(Vector X_In);

    void Set_X(const Vector X_In);
    Vector Get_X(void) const;
    Vector Get_x(void) const;

    void Set_Neighbor_List(const unsigned int *List, const unsigned int N);

    void Update_x(const double dt);
};

#endif
