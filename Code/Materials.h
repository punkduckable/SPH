#if !defined MATERIALS_HEADER
#define MATERAIALS_HEADER

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
                                           1};             // Material density           : g/(mm^3)

  const Material Stainless_Steel       = { "Stainless Steel",
                                           104.374,        // Lame parameter             : Mpa
                                           75.5814,        // Shear modulus              : Mpa
                                           195,            // Young's modulus (E)        : Mpa
                                           8.0};           // Material density           : g/(mm^3)

} // namespace Materials {

#endif
