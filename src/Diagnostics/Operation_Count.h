#if !defined(OPERATION_COUNT_HEADER)
#define OPERATION_COUNT_HEADER

//#define OPERATION_COUNT

namespace OP_Count {
  /* Floating point Operation counters */
  extern unsigned long Addition;
  extern unsigned long Subtraction;
  extern unsigned long Multiplication;
  extern unsigned long Division;
  extern unsigned long Modulus;
  extern unsigned long Sqrt;
  extern unsigned long Exp;
  extern unsigned long Cos;
  extern unsigned long Acos;

  /* Functions */
  void Print(void);
} // namespace OP_Count {

#endif
