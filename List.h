#if !defined(_LIST_HEADER)
#define _LIST_HEADER

class List {
  private:
    struct Node {
      unsigned int Particle_ID;
      Node *Next_Node;
      Node *Previous_Node;
    };

    Node* Front;                         // Start of the list
    Node* End;                           // End of the list
    int Num_Nodes;                       // Number of nodes in the list

  public:
    List();                              // Default constructor
    List(const List & L_In);             // Copy constructor (to prevent inadvertent shallow copying)
    ~List();                             // Destructor

    List & operator=(const List & L_In); // Explicit = operator overload (to prevent inadvertent shallow copying)

    void Add_End(const int ID_In);       // Add an item to the front of the list
    void Add_Front(const int ID_In);     // Add an item to the end of the list
    int Remove_End(void);                // Remove an item from the end of the list
    int Remove_Front(void);              // Add an item to the front of the list

    int Node_Count(void) const;          // Returns number of nodes
    void Cycle_Nodes(void) const;      // For testing
};

#endif
