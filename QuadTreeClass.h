#include "body.h"
#include "Vector3D.h"

extern auto Fg( Body, Body ) -> Vector3D;

class QuadTreeNode
{

 private:  // these attributes are only visible by members of the class
  
    // location+size parameters of each node (center)
    static const double LENGTH; 	   // = 1.0 for DEPTH = 0
    static const double HORIZ_OFFSET;	   // = 0.5*LENGTH for DEPTH = 0
    static const double VERT_OFFSET;	   // "                        "

    // quadtree subnode int ID's
    static constexpr int SUBNODE_A_ID { 1 };
    static constexpr int SUBNODE_B_ID { 2 };
    static constexpr int SUBNODE_C_ID { 3 };
    static constexpr int SUBNODE_D_ID { 4 };
    
    // depth of the quadtree
    static const int DEPTH;
    static constexpr int MAX_DEPTH { 10 };

    // threshold of far-ness for force calculations
    static constexpr double THRESHOLD { 1.0 };

 public:
  
    // the total mass and CoM for all bodies inserted
    // below this node in the tree
    Body BODY;

    // the four subnodes / quadrants
    QuadTreeNode* SUBNODE_A = NULL;
    QuadTreeNode* SUBNODE_B = NULL;
    QuadTreeNode* SUBNODE_C = NULL;
    QuadTreeNode* SUBNODE_D = NULL;

    // node boundaries
    static auto NODE_XMAX { HORIZ_OFFSET + 0.5*LENGTH };
    static auto NODE_XMIN { HORIZ_OFFSET - 0.5*LENGTH };
    static auto NODE_YMAX { VERT_OFFSET + 0.5*LENGTH };
    static auto NODE_YMIN { VERT_OFFSET - 0.5*LENGTH };

    // used to initialize root node parameters
    bool IS_ROOT_NODE = true;

/* ----------------- QuadTree Methods ------------------ */

    // return the subnode a body is in
    auto subnode ( Body temp ) -> int {
    
      // cerr << "Searching which node the body is in..." << endl;
    	
      // individual body x and y positions
      auto x { temp.POSITION.Getx() };
      auto y { temp.POSITION.Gety() };

      // get subnode boundaries
      auto minX { this->NODE_XMIN };    auto minY { this->NODE_YMIN };
      auto midX { this->HORIZ_OFFSET }; auto midY { this->VERT_OFFSET };
      auto maxX { this->NODE_XMAX };    auto maxY { this->NODE_YMAX };
      
      // 1 = first quadrant (A), 2 = second quadrant (B), etc.
      if       ( ( x >= midX ) && ( x <= maxX ) &&
	         ( y >= midY ) && ( y <= maxY ) ) {
	return( this->SUBNODE_A_ID );
      }else if ( ( x >= minX ) && ( x < midX ) &&
		 ( y >= midY ) && ( y <= maxY ) ) {
	return( this->SUBNODE_B_ID );
      }else if ( ( x >= minX ) && ( x <= midX ) &&
		 ( y >= minY ) && ( y < midY ) ) {
	return( this->SUBNODE_C_ID );
      }else if ( ( x > midX ) && ( x <= maxX ) &&
		 ( y >= minY ) && ( y < midY ) ) {
	return( this->SUBNODE_D_ID );
      }else{
	cerr << "WARNING: Attempted to find containing subnode for a body but the body is not in a subnode !" << endl;
	return( 0 );
      }
    }
	
    // Methods to create new subnodes when necessary:
    // define new parameters of subnodes
    auto createSubnode ( int subnode_ID ) -> void {
      
      auto subnode { new QuadTreeNode() };      // create node
      subnode->IS_ROOT_NODE = false;	 // subnodes are not a root node
      subnode->DEPTH = this->DEPTH + 1;	 // depth increases by 1
      subnode->LENGTH = this->LENGTH/2.; // length of subnode is half of parent node

      // find offset based off of subnode ID
      if       ( subnode_ID == 1 ) {
	subnode->HORIZ_OFFSET = this->HORIZ_OFFSET + 0.25*this->LENGTH;
	subnode->VERT_OFFSET = this->VERT_OFFSET + 0.25*this->LENGTH;	
      }else if ( subnode_ID == 2 ) {
	subnode->HORIZ_OFFSET = this->HORIZ_OFFSET - 0.25*this->LENGTH;
	subnode->VERT_OFFSET = this->VERT_OFFSET + 0.25*this->LENGTH;	
      }else if ( subnode_ID == 3 ) {
	subnode->HORIZ_OFFSET = this->HORIZ_OFFSET - 0.25*this->LENGTH;
	subnode->VERT_OFFSET = this->VERT_OFFSET - 0.25*this->LENGTH;	
      }else if ( subnode_ID == 4 ) {
	subnode->HORIZ_OFFSET = this->HORIZ_OFFSET + 0.25*this->LENGTH;
	subnode->VERT_OFFSET = this->VERT_OFFSET - 0.25*this->LENGTH;	
      }else{
	cerr << "WARNING: subnode offsets not set !" << endl;
      }

    }
	
    // checks if a body is in the current node
    auto in_Node ( Body temp ) -> bool {

      auto x{ temp.POSITION.Getx() };
      auto y{ temp.POSITION.Gety() };
		
      if ( ( x >= NODE_XMIN ) && ( x <= NODE_XMAX ) &&
	   ( y >= NODE_YMIN ) && ( y <= NODE_YMAX ) ) return true;
      return false;  // not in node if above conditions not met
      
    }
	
    // insert into appropriate subnode
    auto insert_into_subnode ( Body temp ) -> void {
			
      if       ( subnode(temp) == 1 ) {	   // body is in 1st quadrant
	// create 1st quadrant node if it doesn't already exist
	if ( !this->SUBNODE_A ) createSubnode( this->SUBNODE_A_ID );
	// insert body into 1st quadrant node
	this->SUBNODE_A->insert( temp );
      }else if ( subnode(temp) == 2 ) {    // body is in 2nd quadrant
	if ( !this->SUBNODE_B ) createSubnode( this->SUBNODE_B_ID );
	this->SUBNODE_B->insert( temp );
      }else if ( subnode(temp) == 3 ) {    // body is in 3rd quadrant
	if ( !this->SUBNODE_C ) createSubnode( this->SUBNODE_C_ID );
	this->SUBNODE_C->insert( temp );
      }else if ( subnode(temp) == 4 ) {    // body is in 4th quadrant
	if ( !this->SUBNODE_D ) createSubnode( this->SUBNODE_D_ID );
	this->SUBNODE_D->insert( temp );
      }else{
	cerr << "WARNING: attempted to create appropriate subnode for a body, but the body is not in any subnodes !" << endl;
	return;
      }
      
    }

    // check if current node is internal or external
    auto is_internal() -> bool {
      // if current node has existing subnodes, it is internal
      if ( this->A ) return true;
      if ( this->B ) return true;
      if ( this->C ) return true;
      if ( this->D ) return true;
      return false;	// otherwise, this is an external node
    }

    ///
    /// this is the main function for constructing the QuadTree
    ///
    // We begin at the root node and insert bodies one by one,
    // inserting them into subnodes whenever appropriate.
    // This process is carried out for each node, traveling along
    // each branch until the body is in an external node.
    
    // insert body into this node
    auto insert ( Body next ) -> void {
    // consider 'the next' body in the system to insert into the tree

      // prevent subdividing past max depth
      if ( this->DEPTH > QuadTreeNode::MAX_DEPTH ) return;
		
      // initialize values if this is a root node
      if ( this->IS_ROOT_NODE ) {
	this->DEPTH = 0;	      // 1st level
	this->LENGTH = 1.0;	      // max length
	this->HORIZ_OFFSET = 0.0;     // node at origin
	this->VERT_OFFSET = 0.0;
      }

      // do not insert if the body is outside this node
      if( !in_Node( next ) ) return;
		
      // if above checks out, then follow inserting algorithm below
		
      // if there's no body in this node, put the next body in here
      if ( this->BODY.MASS == 0 ) {
	this->BODY = next;
	return;
      }

      if ( is_internal() ) {	// if this is an internal node
			
	// update this node's values before inserting the next body
			
	// update the node's body location (the center of mass)
	this->BODY.POSITION = center_of_mass( next );
	// update node's body mass (all masses in and below this node)
	this->BODY.MASS += next.MASS;

	// recursively insert the next body into appropriate subnode
	insert_into_subnode( next );	 
			   	
      } else { 	// if this is an external node

	// update this node's values
	this->BODY.POSITION = center_of_mass( next );
	this->BODY.MASS += next.MASS;

	// since this is an external node,
	// then this node contains two bodies, so we must
	// recursively insert this node's (individual) body
	// and the next body into their appropriate subnodes
	insert_into_subnode( this->BODY );
	insert_into_subnode( next );
 
      }

    }

    // finds center of mass between this node's body and another body
    auto center_of_mass ( Body temp ) -> Vector3D {
      return { ( this->BODY.MASS*this->BODY.POSITION
		 + temp.MASS*temp.POSITION )
	  * ( 1. / ( this->BODY.MASS + temp.MASS ) ) };
    }
    
    // now we calculate the forces on each body using nodes
    // starting from the root
	
    // acummulate the net force on the next body
    Vector3D netforce_on(Body next) { 
   		
      auto Fnet{ Vector3D() };
      auto dvdt{ Vector3D() };
    	
      // ignore the next body if it's outside the tree
      if( ( this->IS_ROOT_NODE ) && ( !in_Node( next ) ) ) {
	return(Fnet);
      }

      // if the current node is external and contains a body other than the next body
      if ( ( !is_internal() ) && ( this->BODY.MASS != next.MASS ) ) {
			
	// calculate the forces and add it to the netforce
	Fnet += Fg( this->BODY, next );
						
      } else if ( is_internal() ) {    // if current node is internal
			
	// check how far the next body is from the node body's position
	auto dr{ next.POSITION - this->BODY.POSITION };
	auto dr_mag{ dr.GetMagnitude() };
						
	// this is the parameter to decide what is 'far'
	auto dist_ratio{ (this->LENGTH) / dr_mag };
								
	if ( dist_ratio < QuadTreeNode::THRESHOLD ){ 	// far enough
								
	  // add gravitational force due to node's body mass
	  Fnet += Fg( this->BODY, next );
								
	} else{	     // too close, check subnode forces
						
	  // check if subnode exists and if it does, recursively
	  // move down subnodes until dist_ratio < threshold or
	  // until an external node is reached
				
	  if ( this->SUBNODE_A ) {		
	    Fnet += this->SUBNODE_A->netforce_on( next );
	  }
	  if ( this->SUBNODE_B ) {		
	    Fnet += this->SUBNODE_B->netforce_on( next );
	  }
	  if ( this->SUBNODE_C ) {		
	    Fnet += this->SUBNODE_C->netforce_on( next );
	  }
	  if ( this->SUBNODE_D ) {		
	    Fnet += this->SUBNODE_D->netforce_on( next );
	  }
	        		
	}
	
      }
					
      // return the accumulated net force	
      return( Fnet );
	
    }
    
    // delete the quadtree after use
    void free() {
    
      if ( this->SUBNODE_A ) {   // if this subnode exists, delete its
	this->SUBNODE_A->free(); // subnodes recursively until a NULL
	delete this->SUBNODE_A;	 // subnode, then delete everything
      }			         // starting from the bottom
      if (this->SUBNODE_B) {
	this->SUBNODE_B->free();
	delete this->SUBNODE_B;
      }
      if (this->SUBNODE_C) {
	this->SUBNODE_C->free();
	delete this->SUBNODE_C;
      }
      if (this->SUBNODE_D) {
	this->SUBNODE_D->free();
	delete this->SUBNODE_D;
      }
      return;
    }

    // method to print the tree values
    void print() {
      cerr << "Tree: " << endl;
      cerr << "      mass: " << BODY.MASS << endl;
      cerr << "      position: " << BODY.POSITION.printt() << endl;
      cerr << "      -----------" << endl;

      if ( this->SUBNODE_A ) {
	cerr << "          A";
	this->SUBNODE_A->print();
      }
      if ( this->SUBNODE_B ) {
	cerr << "          B";
	this->SUBNODE_B->print();
      }
      if ( this->SUBNODE_C ) {
	cerr << "          C";
	this->SUBNODE_C->print();
      }
      if ( this->SUBNODE_D ) {
	cerr << "          D";
	this->SUBNODE_D->print();
      }
		
    }

};
