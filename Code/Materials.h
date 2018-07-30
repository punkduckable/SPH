#if !defined MATERIALS_HEADER
#define MATERAIALS_HEADER

namespace Materials {
  // Material structure
  struct Material {
    std::string Name;
    double Lame;                                           // Lame parameter             : Mpa
    double mu0;                                            // Shear modulus              : Mpa
    double E;                                              // Hourglass stiffness        : Mpa
  }; // struct Material {

  // Materials
  const Material Default               = { "Default",
                                           1.125,          // Lame parameter             : Mpa
                                           .275,           // Shear modulus              : Mpa
                                           0.770982};      // Young's modulus (E)        : Mpa

  const Material Stainless_Steel       = { "Stainless Steel",
                                           104.374,        // Lame parameter             : Mpa
                                           75.5814,        // Shear modulus              : Mpa
                                           195};           // Young's modulus (E)        : Mpa
} // namespace Materials {

#endif
