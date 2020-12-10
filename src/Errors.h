#if !defined(ERRORS_HEADER)
#define ERRORS_HEADER

/* File Description:
This file houses all of the exception classes for my SPH code. */

/* Base esception class
All exceptions are based off of this class. */
class Exception {
  static const unsigned Max_Msg_Length = 512;

  private:
    char Message[Max_Msg_Length];
  public:
    Exception(const char* Message_In) {
      /* Note: We do not allow the last character to be populated. By default,
      C initializes char arrays to 0. Thus, not allowing the last cell of
      Message to be populated ensures that Message will always be a
      null-terminated char array. */
      for(unsigned i = 0; i < Max_Msg_Length-1; i++) {
        if(Message_In[i] == '\0') { break; }
        else { Message[i] = Message_In[i]; }
      } // for(unsigned i = 0; i < Max_Msg_Length-1; i++) {
    } // Exception(const char* Message_In) {

    const char* what(void) const { return Message; }
}; // class Exception {





////////////////////////////////////////////////////////////////////////////////
/* General exceptions
These exceptions are not specific to any particular class. */

class Divide_By_Zero : public Exception {
  public:
    Divide_By_Zero(const char* Message_In) : Exception(Message_In) {}
}; // class Divide_By_Zero : public Exception {





////////////////////////////////////////////////////////////////////////////////
/* List exceptions
List_Exception: Base List exception class

Empty_List: This is thrown whenever the user attempts to pop or access (via
Front/Back) elements of an empty list. */

class List_Exception : public Exception {
  public:
    List_Exception(const char* Message_In) : Exception(Message_In) {}
}; // class List_Exception : public Exception {


class Empty_List : public List_Exception {
  public:
    Empty_List(const char* Message_In) : List_Exception(Message_In) {}
}; // class Empty_List : public Exception {





////////////////////////////////////////////////////////////////////////////////
/* Array exceptions
Array_Exception: Base Array exception class

Array_Length_Not_Set: This is thrown whenever the user attempts to access array
information before the array length has been set set up.

Array_Length_Already_Set: This is thrown whenever the user attemtps to set the
Length of an Array that already has a Length. */

class Array_Exception : public Exception {
  public:
    Array_Exception(const char* Message_In) : Exception(Message_In) {}
}; // class Array_Exception : public Exception {


class Array_Length_Not_Set : public Array_Exception {
  public:
    Array_Length_Not_Set(const char* Message_In) : Array_Exception(Message_In) {}
}; // class Array_Length_Not_Set : public Exception {


class Array_Length_Already_Set : public Array_Exception {
  public:
    Array_Length_Already_Set(const char* Message_In) : Array_Exception(Message_In) {}
}; // class Array_Length_Already_Set : public Exception {





////////////////////////////////////////////////////////////////////////////////
/* Tensor exceptions
Tensor_Exception: Base Tensor exception class.

Singular_Matrix: This is thrown whenever the user tries to find the inverse of
a singular matrix (Determinant = 0) */

class Tensor_Exception : public Exception {
  public:
    Tensor_Exception(const char* Message_In) : Exception(Message_In) {}
}; // class Tensor_Exception : public Exception {


class Singular_Matrix : public Tensor_Exception {
  public:
    Singular_Matrix(const char* Message_In) : Tensor_Exception(Message_In) {}
}; // class Singular_Matrix : public Tensor_Exception {


class Undefined_Exponent : public Tensor_Exception {
  public:
    Undefined_Exponent(const char* Message_In) : Tensor_Exception(Message_In) {}
}; // class Undefined_Exponent : public Tensor_Exception {





////////////////////////////////////////////////////////////////////////////////
/* Particle exceptions
Particle_Exception: Base Particle exception class

No_BC: This is thrown whenever the user tries to access the ith component of a
particle's boundary condition, but the particle has no boundary condition for
that component. */

class Particle_Exception : public Exception {
  public:
    Particle_Exception(const char* Message_In) : Exception(Message_In) {}
}; // class Particle_Exception : Public Exception {


class No_BC : public Particle_Exception {
  public:
    No_BC(const char* Message_In) : Particle_Exception(Message_In) {}
}; // class No_BC : public Particle_Exception {





////////////////////////////////////////////////////////////////////////////////
/* Body exceptions
Body_Exception: Base Particle exception class

Not_A_Box: This is thrown whenever a non-box body tries to do something that
only boxes can do (call functions that only work for boxes, for example). */

class Body_Exception : public Exception {
  public:
    Body_Exception(const char* Message_In) : Exception(Message_In) {}
}; // class Body_Exception : Public Exception {


class Not_A_Box : public Body_Exception {
  public:
    Not_A_Box(const char* Message_In) : Body_Exception(Message_In) {}
}; // class Not_A_Box : public Body_Exception {





////////////////////////////////////////////////////////////////////////////////
/* IO Exceptions
IO_Exception: Base input/output exception class

Cant_Open_File: This exception is thrown whenever the program can not open a
file that it's supposed to be opening.

Bad_Read: This exception is thrown whenever the program reads in something that
does not match what it was excpecting. */

class IO_Exception : public Exception {
  public:
    IO_Exception(const char* Message_In) : Exception(Message_In) {}
}; // class IO_Exception : public Exception {


class Cant_Open_File : public IO_Exception {
  public:
    Cant_Open_File(const char* Message_In) : IO_Exception(Message_In) {}
}; // class Cant_Open_File : public IO_Exception {


class Bad_Read : public IO_Exception {
  public:
    Bad_Read(const char* Message_In) : IO_Exception(Message_In) {}
}; // class Bad_Read : public IO_Exception {





////////////////////////////////////////////////////////////////////////////////
/* Simulation exceptions
Simulation_Exception: Base Simulation_Exception class.

Bad_Body_Setup: This exception is thrown whenever the bodies are set up
improperly (such as when a body is designated as both a box and from file) */

class Simulation_Exception : public Exception {
  public:
    Simulation_Exception(const char* Message_In) : Exception(Message_In) {}
}; // class Simulation_Exception : public Exception {


class Bad_Body_Setup : public Simulation_Exception {
  public:
    Bad_Body_Setup(const char* Message_In) : Simulation_Exception(Message_In) {}
}; // class Bad_Body_Setup : Public Simulation_Exception {

#endif
