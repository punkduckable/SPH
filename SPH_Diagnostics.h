#if !defined(PARTICLE_FILE_HEADER)
#define PARTICLE_FILE_HEADER

namespace Particle_Debugger {
  using std::string;

  unsigned int File_Number = 0;
  const unsigned int File_Number_Max_Digits = 5;

  void Export_Pariticle_Forces(const unsigned int Num_Particles,
                               const Particle * Particles);
                               
  void Get_File_Name(string & Str);
} // namespace VTK_File {

namespace OP_Count {
  //////////////////////////////////////////////////////////////////////////////
  // Operation count variables

  // Tensors
  unsigned int T_Default_Constructor = 0;
  unsigned int T_Component_Constructor = 0;
  unsigned int T_Copy_Constructor = 0;
  unsigned int T_Equality = 0;
  unsigned int T_T_Addition = 0;
  unsigned int T_T_Subtraction = 0;
  unsigned int T_T_Multiplication = 0;
  unsigned int T_V_Multiplication = 0;
  unsigned int T_S_Multiplication = 0;
  unsigned int T_S_Division = 0;
  unsigned int Compound_T_T_Addition = 0;
  unsigned int Compound_T_T_Subtraction = 0;
  unsigned int Compound_T_T_Multiplication = 0;
  unsigned int T_Inverse = 0;
  unsigned int T_Determinant = 0;
  unsigned int T_Transpose = 0;
  unsigned int T_Dot_Product = 0;

  // Vectors
  unsigned int V_Default_Constructor = 0;
  unsigned int V_Component_Constructor = 0;
  unsigned int V_Copy_Constructor = 0;
  unsigned int V_Equality = 0;
  unsigned int V_V_Addition = 0;
  unsigned int V_V_Subtraction  = 0;
  unsigned int V_S_Multiplication = 0;
  unsigned int V_S_Division = 0;
  unsigned int Compound_V_V_Addition = 0;
  unsigned int Compound_V_V_Subtraction = 0;
  unsigned int Compound_V_S_Multiplication = 0;
  unsigned int V_Magnitude = 0;
  unsigned int V_Dot_Product = 0;

  // Other
  unsigned int Dyadic_Product = 0;

  //////////////////////////////////////////////////////////////////////////////
  // Functions

  void Reset_Counts(void);
  void Print_Counts(void);
} // namespace OP_Count {

#endif
