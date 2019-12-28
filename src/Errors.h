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
  public: Divide_By_Zero(const char* Message_In) : Exception(Message_In) {}
}; // class Divide_By_Zero : public Exception {



////////////////////////////////////////////////////////////////////////////////
/* Tensor exceptions
Tensor_Exception: Base Tensor exception class.

Zero_Determinant: This is thrown whenever it is impossible to calculate some
value because a tensor's determinant is zero. */

class Tensor_Exception : public Exception {
  public: Tensor_Exception(const char* Message_In) : Exception(Message_In) {}
}; // class Tensor_Exception : public Exception {


class Zero_Determinant : public Tensor_Exception {
  public: Zero_Determinant(const char* Message_In) : Tensor_Exception(Message_In) {}
}; // class Zero_Determinant : public Tensor_Exception {


class Undefined_Exponent : public Tensor_Exception {
  public: Undefined_Exponent(const char* Message_In) : Tensor_Exception(Message_In) {}
}; // class Undefined_Exponent : public Tensor_Exception {



////////////////////////////////////////////////////////////////////////////////
/* Particle exceptions
Particle_Exception: Base Particle exception class

Bad_Neighbor_Index: This is thrown whenever the user tries to access a particle's
jth neighbor, when j is greater than the particle's number of neighbors */

class Particle_Exception : public Exception {
  public: Particle_Exception(const char* Message_In) : Exception(Message_In) {}
}; // class Particle_Exception : Public Exception {


class Bad_Neighbor_Index : public Particle_Exception {
  public: Bad_Neighbor_Index(const char* Message_In) : Particle_Exception(Message_In) {}
}; // class Bad_Neighbor_ID : public Particle_Exception {

#endif
