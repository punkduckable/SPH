#if !defined(QUICK_MATH_HEADER)
#define QUICK_MATH_HEADER

#if !defined(PI_SQUARED)
#define PI_SQUARED 9.86960440109
#endif

namespace Quick_Math {
  double Acos(const double x) {
    return (-0.69813170079773212 * x * x - 0.87266462599716477) * x + 1.5707963267948966;
  } // double Acos(const double x) {

  double cos(double x) {
    if(x > PI)
      x -= 2*PI;
    else if(x < -PI)
      x += 2*PI;

    return 1. + ((x*x)/(3.*PI_SQUARED))*(-14. + (8./PI_SQUARED)*(x*x));
  } // double cos(double x) {
} // namespace Quick_Math {

#endif
