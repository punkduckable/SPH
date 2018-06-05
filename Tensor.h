#if !defined(_TENSOR_HEADER)
#define _TENSOR_HEADER

class Tensor {
  public:
    double T[9];

    Tensor(void);                                 // Default constructor
    Tensor(double t11, double t12, double t13,
           double t21, double t22, double t23,
           double t31, double t32, double t33);   // Component constructor

    Tensor operator+(const Tensor S_In) const;
    Tensor operator*(const Tensor S_In) const;
    Vector operator*(const Vector V_In) const;
    Tensor operator=(const double S_In[9]);
    Tensor operator=(const Tensor S_In);
    double& operator()(const uByte row, const uByte col);
    double operator()(const uByte row, const uByte col) const;
    void Print(void);
}; // class Tensor {

#endif
