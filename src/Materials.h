#if !defined MATERIALS_HEADER
#define MATERIALS_HEADER

#include <string>

namespace Materials {
  // Material structure
  struct Material {
    double Lame;                                           // Lame parameter             : Mpa
    double mu0;                                            // Shear modulus              : Mpa
    double E;                                              // Young's Modulus            : Mpa
    double density;                                        // Material density           : g/(mm^3)
  }; // struct Material {

  // Materials
  const Material Default               = { 1.125,          // Lame parameter             : Mpa
                                           .275,           // Shear modulus              : Mpa
                                           0.770982,       // Young's modulus (E)        : Mpa
                                           .01};           // Material density           : g/(mm^3)

  const Material Old_Needle            = { 104.374,        // Lame parameter             : Mpa
                                           75.5814,        // Shear modulus              : Mpa
                                           195.0,          // Young's modulus (E)        : Mpa
                                           .08};           // Material density           : g/(mm^3)

  const Material Stainless_Steel       = { 104374.0,       // Lame parameter             : Mpa
                                           86000.0,        // Shear modulus              : Mpa
                                           195000.0,       // Young's modulus (E)        : Mpa
                                           .08};           // Material density           : g/(mm^3)

} // namespace Materials {

#endif
