#define threshold 0.5 	// threshold for determining "far-ness"
// 0.5 for best results, higher means faster 
// #define G 6.67e-11	// m^3 kg^-1 s^-2, gravitational constant
//#define G 1.56700e-9	// ly^3 Msun^-1 (100yr)^-2, gravitational constant
#define max_depth 15	// max depth of quadtree
//#define scale 2.5*384.4e6 	// m
//#define scale 90		// ly

//#define notes false		// turn notes on/off

struct Body {
    double mass;
    Vector3D position;
    Vector3D velocity;

    // makes a new body
Body(double m, Vector3D pos, Vector3D vel) : mass(m), position(pos), velocity(vel) {};

Body() : mass(0), position(0, 0, 0), velocity(0, 0, 0) {};
};

struct QuadTreeNode;

struct QuadTreeNode {

/* ---------- Stored Properties of this node ----------- */

    // the total mass and CoM for all bodies inserted
    // below this node in the tree
    Body body;
    
    // rescale positions to fit in nodes
    // note: will need to have all bodies in
    // at least positive coordinates to be included
    // body.position.update(body.position.Getx()/scale, 
    // 						body.position.Gety()/scale);
    // scale defined in program

    // the four subnodes / quadrants
    QuadTreeNode* A = NULL;
    QuadTreeNode* B = NULL;
    QuadTreeNode* C = NULL;
    QuadTreeNode* D = NULL;
    
    // location+size parameters of each node (center)
    double length; 			// = 1.0 for depth = 0
    double horiz_offset;	// = 0.5*length for depth = 0
    double vert_offset;		// "                        "
    
    // used to initialize root node parameters
    bool root_node = true;
    
    // depth of the quadtree
    int depth;
  	
/* ----------------------------------------------------- */

/* ----------------- QuadTree Methods ------------------ */

    // search which subnode a body is in
    int Subnode(Body temp){
    
		// cerr << "Searching which node the body is in..." << endl;
    	
    	// individual body x and y positions
		double x = temp.position.Getx();
		double y = temp.position.Gety();
		
		// define subnode boundaries
		double node_xmax = this->horiz_offset + 0.5*this->length;
		double node_xmin = this->horiz_offset - 0.5*this->length;
		double node_ymax = this->vert_offset + 0.5*this->length;
		double node_ymin = this->vert_offset - 0.5*this->length;
		
		// position of edges between quadrants
		double AB = node_xmax - 0.5*this->length;	// = CD
		double BC = node_ymax - 0.5*this->length;	// = AD
			
		// 1 = first quadrant (A), 2 = second quadrant (B), etc.
		if 	 	 ((x >= AB) && (x <= node_xmax) && 
			  	  (y >= BC) && (y <= node_ymax)){return(1);
		}else if ((x >= node_xmin) && (x < AB) && 
			  	  (y >= BC) && (y <= node_ymax)){return(2);
		}else if ((x >= node_xmin) && (x <= AB) && 
			  	  (y >= node_ymin) && (y < BC)){return(3);
		}else if ((x > AB) && (x <= node_xmax) && 
			  	  (y >= node_ymin) && (y < BC)){return(4);
		}else{ 
			return(0);
		}
    }
	
	// Methods to create new subnodes when necessary:
    // define new parameters of subnodes
    void createA(){			// 1st quadrant / subnode
		this->A = new QuadTreeNode();				// create node
		this->A->root_node = false;					// sub nodes are not the root node
		this->A->depth = this->depth + 1;			// depth increases by 1
		this->A->length = this->length/2.;			// length of subnode is half of parent node
		this->A->horiz_offset = this->horiz_offset + 0.25*this->length;
		this->A->vert_offset = this->vert_offset + 0.25*this->length;
    }
	    
    void createB() {   		// 2nd quadrant / subnode
		this->B = new QuadTreeNode();				// create node
		this->B->root_node = false;					// sub nodes are not the root node
		this->B->depth = this->depth + 1;			// depth increases by 1
		this->B->length = this->length/2.;			// length of subnode is half of parent node
		this->B->horiz_offset = this->horiz_offset - 0.25*this->length;
		this->B->vert_offset = this->vert_offset + 0.25*this->length;
    }
	    
    void createC() {		// 3rd quadrant / subnode
		this->C = new QuadTreeNode();				// create node
		this->C->root_node = false;					// sub nodes are not the root node
		this->C->depth = this->depth + 1;			// depth increases by 1
		this->C->length = this->length/2.;			// length of subnode is half of parent node
		this->C->horiz_offset = this->horiz_offset - 0.25*this->length;
		this->C->vert_offset = this->vert_offset - 0.25*this->length;
    }
	    
    void createD() {	    // 4th quadrant / subnode
	this->D = new QuadTreeNode();				// create node
	this->D->root_node = false;					// sub nodes are not the root node
	this->D->depth = this->depth + 1;			// depth increases by 1
	this->D->length = this->length/2.;			// length of subnode is half of parent node
	this->D->horiz_offset = this->horiz_offset + 0.25*this->length;
	this->D->vert_offset = this->vert_offset - 0.25*this->length;
    }
	
	// checks if a body is in the current node
	bool in_Node(Body temp){

		double x = temp.position.Getx();
		double y = temp.position.Gety();
			
		// define tree boundaries
		double node_xmax = this->horiz_offset + 0.5*this->length;
		double node_xmin = this->horiz_offset - 0.5*this->length;
		double node_ymax = this->vert_offset + 0.5*this->length;
		double node_ymin = this->vert_offset - 0.5*this->length;
		
		if ((x >= node_xmin) && (x <= node_xmax) &&
			(y >= node_ymin) && (y <= node_ymax)) return true;
		return false;
	}
	
    // insert body into appropriate subnode
    void insert_Subnode(Body temp){
			
		// insert body recursively
		if (Subnode(temp) == 1) {			// body in 1st quadrant
			if (!this->A) createA();		// create 1st quadrant node if it doesn't already exist
			this->A->insert(temp);			// insert
		}
		else if (Subnode(temp) == 2) {		// body in 2nd quadrant
			if (!this->B) createB();		// create 2nd quadrant node if it doesn't already exist
			this->B->insert(temp);			// insert
		}
		else if (Subnode(temp) == 3) {		// body in 3rd quadrant
			if (!this->C) createC();		// create 3rd quadrant node if it doesn't already exist
			this->C->insert(temp);			// insert
		}
		else if (Subnode(temp) == 4){		// body in 4th quadrant
			if (!this->D) createD();		// create 4th quadrant node if it doesn't already exist
			this->D->insert(temp);			// insert
		}else return;
    }

    // check if current node is internal or external
    bool is_internal(){
    	// if current node has existing subnodes, it is internal
		if (this->A) return true;
		if (this->B) return true;
		if (this->C) return true;
		if (this->D) return true;
		return false;	// otherwise, this is an external node
    }

    // this is the main function for constructing the QuadTree
    // insert body into this node
    void insert(Body next) {		// consider the next body in the system
									// to insert into the tree

		if (this->depth > max_depth) return;	// prevent subdividing past max depth
		
		// initialize values if this is a root node
		if (this->root_node){
			this->depth = 0;			// 1st level
			this->length = scale;		// max length
			this->horiz_offset = 0.0; 	// 0.5*this->length; for positive quadrant
			this->vert_offset = 0.0; 	// '									 '
		}
		
		if(!in_Node(next)) return;		// do not insert if the body is outside this node
		
		// if above checks out, then follow inserting algorithm below
		
		// if there's no body in this node, put the next body in here
		if (this->body.mass == 0) {
			this->body = next;
			return;
		}

		if (is_internal()) {	// if this is an internal node
			
			// update this node's values before inserting the next body
			
			// update the node's body location (the center of mass)
			this->body.position = (this->body.mass*this->body.position 
						+ next.mass*next.position)*(1./(this->body.mass + next.mass));
			// update the node's body mass (total mass of masses in and below this node)
			this->body.mass += next.mass;

			// recursively insert the next body into appropriate subnode
			insert_Subnode(next);	 
			   	
		} else { 	// if this is an external node
			// since the current node's body mass is not 0, then
			// this node contains a body and so we must
			// recursively insert this node's (individual) body
			// and the next body into their appropriate subnodes
			insert_Subnode(this->body);
			insert_Subnode(next);
				
			// then update this node's body position (center of mass)
			this->body.position = (this->body.mass*this->body.position 
						+ next.mass*next.position)*(1./(this->body.mass + next.mass));
			// and update the node's body mass (total mass)
			this->body.mass += next.mass;
		}
		
		// return nothing if above is not satisfied for any reason
		return;

    }

	// now we calculate the forces on each body using nodes
	// starting from the root

    // calculates gravitational force
    Vector3D Fg(Body temp){

		if(notes) cerr << "Calculating grav force..." << endl;

		double M1, M2, dist;
		Vector3D  Dr, F;
			
		M1 = this->body.mass;
		M2 = temp.mass;

		// this takes the position vector difference
		Dr =  temp.position - this->body.position;
		// and this is the magnitude
		dist = pow(Dr.GetMagnitude(),2.) + 0.0001;			// small epsilon for softening
		
		if(notes) cerr << "M1, M2, Dist = " << M1 << ", " << M2 << ", " << pow(dist,0.5) << endl;

		// here is the gravitational force
		F = (-G*M1*M2/(pow(dist,1.5))) * Dr;
			
		if(notes) cerr << "Fg = " << F.printt() <<  endl;
			
		return(F);
		
    }
	
    // acummulate the net force on the next body
    Vector3D netforce_on(Body next) { 
   		
   		if(notes) cerr << "Finding net force on mass " << next.mass << endl;
   		if(notes) cerr << "Current node: " << this->depth << ", " << this->body.mass << endl;
   		
    	Vector3D Fnet(0.0, 0.0, 0.0);
    	Vector3D dvdt;
    	
    	// ignore the next body if it's outside the tree
    	if((this->root_node) && (!in_Node(next))) {
    		if(notes) cerr << "Not in tree" << endl;
    		return(Fnet);
    	}

		// if the current node is external and not containing only the next body
		if ((!is_internal()) && (this->body.mass != next.mass)){
			
			if(notes) cerr << "External calc..." << endl;
   			if(notes) cerr << "Current node: " << this->depth << ", " << this->body.mass << endl;
			
			// simply calculate the grav force by two bodies and add it to the netforce
			Fnet = Fnet + Fg(next);
			
			if(notes) cerr << "External: Fnet = " << Fnet.printt() << endl;
			
		} else if (is_internal()){		// if current node is internal
		
			if(notes) cerr << "Internal calc..." << endl;
			if(notes) cerr << "Current node: " << this->depth << ", " << this->body.mass << endl;
			
			// check how far the next body is from the node body's center of mass
			Vector3D Dr;
			double dist;
			Dr =  next.position - this->body.position;
			dist = Dr.GetMagnitude();
			
			if(notes) cerr << "Dist = " << dist << endl;
			
			// this is the parameter to decide what is 'far'
			double dist_ratio = (this->length) / dist;
			
			if(notes) cerr << "Dist ratio = " << dist_ratio << endl;
					
			if (dist_ratio < threshold){ 	// far enough
				
				if(notes) cerr << "Far enough" << endl;
				
				// add gravitational force due to node's body mass
				Fnet = Fnet + Fg(next);
				
				if(notes) cerr << "Internal, far enough: Fnet = " << Fnet.printt() << endl;
				
			} else{		// too close, check subnode forces
			
				if(notes) cerr << "Too close, checking subnodes" << endl;
			
				// check if subnode exists and if it does, recursively
				// move down subnodes until dist_ratio < threshold or
				// until an external node is reached
				
				if (this->A){
				
					if(notes) cerr << "Checking subnode A" << endl;
				
					Fnet = Fnet + this->A->netforce_on(next);
				}
				if (this->B){
									
					if(notes) cerr << "Checking subnode B" << endl;
				
					Fnet = Fnet + this->B->netforce_on(next);
				}
				if (this->C){
									
					if(notes) cerr << "Checking subnode C" << endl;
				
					Fnet = Fnet + this->C->netforce_on(next);
				}
				if (this->D){
									
					if(notes) cerr << "Checking subnode D" << endl;
				
					Fnet = Fnet + this->D->netforce_on(next);
				}
				
				if(notes) cerr << "Finished checking subnodes" << endl;
				
			}
		}
		
		if(notes) cerr << "Accumulated force = " << Fnet.printt() << endl;
			
		// return the accumulated net force	
		return(Fnet);
	
    }
    
    // delete the quadtree after use
    void free() {
    
    	if (this->A) {					// if this subnode exists, delete
			this->A->free();			// its subnodes (keeps checking) until a
			delete this->A;				// NULL subnode is reached, then delete everything
    	}								// starting from the bottom
    	if (this->B) {
			this->B->free();
			delete this->B;
    	}
    	if (this->C) {
			this->C->free();
			delete this->C;
    	}
    	if (this->D) {
			this->D->free();
			delete this->D;
    	}
    	return;
    }

	// method to print the tree values
    void print() {
		cerr << "Tree: " << endl;
		cerr << "      mass: " << body.mass << endl;
		cerr << "      position: " << body.position.Getx() << ", " << body.position.Gety() 
						   << ", " << body.position.Getz() << endl;
		cerr << "      -----------" << endl;

		if (this->A) {
			cerr << "          A";
			this->A->print();
		}
		if (this->B) {
			cerr << "          B";
			this->B->print();
		}
		if (this->C) {
			cerr << "          C";
			this->C->print();
		}
		if (this->D) {
			cerr << "          D";
			this->D->print();
		}
		
    }

};
