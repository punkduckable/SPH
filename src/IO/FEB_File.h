#if !defined(FEB_HEADER)
#define FEB_HEADER

#include "Classes.h"
#include <string>

namespace FEB_File {
  int Read_FEB_File(const std::string & File_Name, Vector ** X_Ptr, unsigned & Num_Nodes);
} // namespace FEB_File {

#endif
