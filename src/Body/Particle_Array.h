#if !defined PARTICLE_ARRAY_HEADER
#define PARTICLE_ARRAY_HEADER

#if !defined(PI)
  #define PI 3.1415926535897932384626
#endif

class Particle_Array {
  private:
    // Particle array name
    std::string Name;

    // Cuboid parameters.
    bool Is_Cuboid = false;
    unsigned int X_SIDE_LENGTH = 0;
    unsigned int Y_SIDE_LENGTH = 0;
    unsigned int Z_SIDE_LENGTH = 0;

    // Boundary parameters
    bool Is_Boundary = false;

    // The array of particles
    unsigned int Num_Particles = 0;
    Particle * Array = NULL;

    // Kernel parameters
    double Inter_Particle_Spacing = 0;                     //                            : mm
    unsigned int Support_Radius = 0;                       // Support radius in Inter Particle Spacints's    : unitless
    double h = 0;                                          // Support radius in mm's     : mm
    double Shape_Function_Amplitude = 0;                   // Shape function Amplitude   : 1/(mm^6)

    // Material Parameters
    Materials::Material Array_Material;

    // Viscosity parameters
    double mu;                                             // Viscosity                  : Mpa*s
    unsigned char F_Index = 0;                             // Keeps track of which F (in each partilce) is the 'new' one

    // Hourglass (Hg) correction parameters
    double alpha;                                          // Hg control parameter       : unitless

    // Damage paramaters
    double Tau;                                            // Damage rate parameter (see eq 26)
    bool Damageable = 1;                                   // If true, allows this Particle_Array to take damage.

    // First time step?
    bool First_Time_Step = true;

    // Private methods
    void Set_h(const double h_in) {
      h = h_in;
      Shape_Function_Amplitude =  15./(PI*h*h*h*h*h*h);
    } // void Set_h(const double h_in) {

  public:
    // Constructors, destructor
    Particle_Array(void);                                  // default constructor
    Particle_Array(const unsigned int Num_Particles_In);   // generate array constructor
    Particle_Array(const Particle_Array & Ar_In) {         // do nothing copy constructor
      printf("Copy constructor is NOT enabled for Particle Arrays\n");
    }
    ~Particle_Array(void);                                 // destructor

    // Operator overloading
    Particle_Array & operator=(Particle_Array & Ar_In) {   // do nothing = operator
      printf("Equality is NOT defined for Particle arrays\n");
      return Ar_In;
    }
    Particle & operator[](const unsigned int i) {
      return Array[i];
    }
    const Particle & operator[](const unsigned int i) const {
      return Array[i];
    }


    // Set methods
    void Set_Num_Particles(const unsigned int Num_Particles_In);
    void Set_Name(const std::string & S_In) { Name = S_In; }

    void Set_Inter_Particle_Spacing(const double IPS);
    void Set_Support_Radius(const unsigned int SR_In);

    void Set_Material(const Materials::Material & Mat_In) { Array_Material = Mat_In; }
    void Set_mu(const double mu_In) { mu = mu_In; }
    void Set_alpha(const double alpha_In) { alpha = alpha_In; }

    void Set_Tau(const double Tau_In) { Tau = Tau_In; }
    void Set_Damageable(const bool D_In) { Damageable = D_In; }

    void Set_Cuboid_Dimensions(const Vector & Dimensions);

    void Set_Boundary(const bool Boundary_In) { Is_Boundary = Boundary_In; }

    void Set_First_Time_Step(const bool First_In) { First_Time_Step = First_In; }

    void Set_F_Index(const unsigned char i) {
      if(i > 1) {
        printf("Error! F_Counter must be 0 or 1!!! You supplied %u\n", i);
        return;
      } // if(i > 1) {
      F_Index = i;
    } // void Set_F_Counter(const unsigned char i) {

    void Increment_F_Index(void) {
      if(F_Index == 0)
        F_Index++;
      else
        F_Index = 0;
    } // void Increment_F_Counter(void) {

    // Getters
    unsigned int Get_Num_Particles(void) const { return Num_Particles; }
    std::string Get_Name(void) const { return Name; }

    double Get_Inter_Particle_Spacing(void) const { return Inter_Particle_Spacing; }
    unsigned int Get_Support_Radius(void) const { return Support_Radius; }

    double Get_h(void) const { return h; }
    double Get_Shape_Function_Amplitude(void) const { return Shape_Function_Amplitude; }
    Materials::Material Get_Material(void) const { return Array_Material; }
    double Get_Lame(void) const { return Array_Material.Lame; }
    double Get_mu0(void) const { return Array_Material.mu0; }
    double Get_mu(void) const { return mu; }
    double Get_E(void) const { return Array_Material.E; }
    double Get_density(void) const { return Array_Material.density; }
    double Get_alpha(void) const { return alpha; }
    unsigned char Get_F_Index(void) const { return F_Index; }

    double Get_Tau(void) const { return Tau; }
    bool Get_Damagable(void) const { return Damageable; }

    bool Get_Cuboid(void) const { return Is_Cuboid; }
    unsigned int Get_X_SIDE_LENGTH(void) const { return X_SIDE_LENGTH; }
    unsigned int Get_Y_SIDE_LENGTH(void) const { return Y_SIDE_LENGTH; }
    unsigned int Get_Z_SIDE_LENGTH(void) const { return Z_SIDE_LENGTH; }

    bool Get_Boundary(void) const { return Is_Boundary; }

    bool Get_First_Time_Step(void) const { return First_Time_Step; }

    // Other Methods
    void Print_Parameters(void) const;
}; // class Particle_Array {

#endif
