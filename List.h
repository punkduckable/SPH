#if !defined(LIST_HEADER)
#define LIST_HEADER

class List {
  private:
    struct Node {
      unsigned int Particle_ID;
      Node *Next_Node;
      Node *Previous_Node;
    };

    Node* Beginning;                     // Start of the list
    Node* End;                           // End of the list
    int Num_Nodes;                       // Number of nodes in the list

  public:
    // Constructors, destructor
    List();                              // Default constructor
    List(const List & L_In);             // Copy constructor (to prevent inadvertent shallow copying)
    ~List();                             // Destructor

    // List equality
    List & operator=(const List & L_In); // Explicit = operator overload (to prevent inadvertent shallow copying)

    // Method to add or remove items from the list
    void Add_Back(const int ID_In);      // Add an item to the start of the list
    void Add_Front(const int ID_In);     // Add an item to the end of the list
    int Remove_Back(void);               // Remove an item from the end of the list
    int Remove_Front(void);              // Add an item to the start of the list

    // Methods to get list info
    int Node_Count(void) const;          // Returns number of nodes
    void Print_Node_Info(void) const;    // For testing
};

#endif
