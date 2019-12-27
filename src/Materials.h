#if !defined MATERIALS_HEADER
#define MATERAIALS_HEADER

#include <string>

namespace Materials {
  // Material structure
  struct Material {
    std::string Name;
    double Lame;                                           // Lame parameter             : Mpa
    double mu0;                                            // Shear modulus              : Mpa
    double E;                                              // Hourglass stiffness        : Mpa
    double density;                                        // Material density           : g/(mm^3)
  }; // struct Material {

  // Materials
  const Material Default               = { "Default",
                                           1.125,          // Lame parameter             : Mpa
                                           .275,           // Shear modulus              : Mpa
                                           0.770982,       // Young's modulus (E)        : Mpa
                                           .01};           // Material density           : g/(mm^3)

  const Material Stainless_Steel       = { "Stainless Steel",
                                           104374.0,        // Lame parameter             : Mpa
                                           86000.0,        // Shear modulus              : Mpa
                                           195000.0,            // Young's modulus (E)        : Mpa
                                           .08};           // Material density           : g/(mm^3)

} // namespace Materials {

#endif
