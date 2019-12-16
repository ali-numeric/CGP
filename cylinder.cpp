// === Code Name: CGP/FEM  Project                                                            
// === Author: Ali A. Kashefi																  	
// === Language: C++                                                                          

// === Flow Past A Cylinder
// === Class Name: Cylinder														   
// === cylinder.cpp

#define _USE_MATH_DEFINES
#include <cmath>
#include "cylinder.h"

bool compareAngle(const cylinder & a , const cylinder & b);

// -------------------- Constructor --------------------

Cylinder :: Cylinder(double Viscosity , double Density , double Time , double Time_Step , string S , double Alfa , double penalty_term , double Beta_q , int level , int n_center)
	   : Computation(Viscosity , Density , Time , Time_Step , S , Alfa , penalty_term , Beta_q , level , n_center){
		
	//https://msdn.microsoft.com/en-us/library/4hwaceh6.aspx

	// Re = (U_stream)(D)/(dynamic_viscosiy) = (1)(1)/(dynamic_viscosiy) 
	
	pi = M_PI ;

	dynamic_viscosiy = Viscosity ;
	density = Density ;

	U_stream = 1.0 ;
	H = 32.0 ;
	L_Box = 38.0 ;
	D = 1.00 ;

	X_Center = 8.0 ;
	Y_Center = 16.0 ;
}

// -------------------- Initialize Velocity --------------------

void Cylinder :: initializeVelocity(){
	
	// I.C.

	ifstream fin ("MeshP1.txt");
	
	int Node_Number , vain ;
	
	fin >> Node_Number ;
	fin >> vain ;

	double x , y , z ;
	
	int counter = 0 ; 

	for (int i = 0 ; i < Node_Number ; i ++){
		
		fin >> x ;  
		fin >> y ;
		fin >> z ;
		
		if (y == 0.0){
		velocity_u[0][counter] = U_stream ;
		velocity_v[0][counter] = 0.0 ; 
		velocity[0][counter] = velocity_u[0][counter] ;
		velocity[0][counter + Node_Number] = velocity_v[0][counter] ; counter ++ ; continue ;}
		
		if (x == 0.0 || y == H ){
		velocity_u[0][counter] = U_stream ;
		velocity_v[0][counter] = 0.0 ; 
		velocity[0][counter] = velocity_u[0][counter] ;
		velocity[0][counter + Node_Number] = velocity_v[0][counter] ; counter ++ ; continue ;}

		velocity_u[0][counter] = 0.0 ;
		velocity_v[0][counter] = 0.0 ; 
		velocity[0][counter] = velocity_u[0][counter] ;
		velocity[0][counter + Node_Number] = velocity_v[0][counter] ;
		
		counter ++ ;
	}

	fin.close() ;

	// B.C.
	
	Mesh Initializr("MeshP1.txt" , 1 , N_center) ;
	Initializr.generateMesh() ;

	for(int k = 0 ; k < time_step ; k ++){
		for (int i = 0 ; i < Dirichlet_bound_number_velocity ; i ++){
			for (int j = 0 ; j < Initializr.boundary_nodes_counter ; j ++){
				
				if (Initializr.boundary_nodes[j].node_number == Dirichlet_bound_velocity[i]){
			
				x = Initializr.boundary_nodes[j].x_coordinate ;
				y = Initializr.boundary_nodes[j].y_coordinate ;
			
				double di = sqrt((x-X_Center)*(x-X_Center) + (y-Y_Center)*(y-Y_Center)) ;
			
				if ( di < 1.5 ){
				intermediate_u_BC[k][i] = 0.0 ; 
				intermediate_v_BC[k][i] = 0.0 ; continue ;}

				intermediate_u_BC[k][i] = U_stream ; 
				intermediate_v_BC[k][i] = 0.0 ;
				
				}
			}
		}
	}
	
	Initializr.Destroy() ;
	//Initializr.~Mesh() ;
}

// -------------------- Initialize Pressure--------------------

void Cylinder ::  initializePressure(){

	// I.C.

	ifstream fin ("MeshP1.txt") ;
	
	int Node_Number , vain ;
	
	fin >> Node_Number ;
	fin >> vain ;
	
	double x , y , z ;
	
	int counter = 0 ; 

	for (int i = 0 ; i < Node_Number ; i ++){
		
		fin >> x ;  
		fin >> y ;
		fin >> z ;
		
		pressure[0][counter] = 0.0 ;
		pressure_exp[counter] = 0.0 ;

		counter ++ ;
	}

	fin.close() ;
}

// -------------------- Initialize Force --------------------

void Cylinder :: initializeForce(){

	ifstream fin ("MeshP1.txt") ;
	
	int Node_Number , vain ;
	
	fin >> Node_Number ;	fin >> vain ;

	double x , y , z ;

	for(int k = 0 ; k < time_step ; k ++){
		for (int i = 0 ; i < Node_Number ; i ++){
		
			fin >> x ;  
			fin >> y ;
			fin >> z ;

			RHS_Force[k][i] = 0.0 ; //F(x) ;
			RHS_Force[k][i + velocity_node_number] = 0.0 ; //F(y) 
		}
	}

	fin.close() ;
}

// -------------------- Construct Nuemann & Dirichlet Bounds --------------------

void Cylinder :: constructNuemannDirichletBounds(){

	Mesh BC("MeshP1.txt" , 1 , N_center) ;
	BC.generateMesh() ;
		
	Dirichlet_bound_number_pressure = 0 ; //BC.boundary_nodes_counter  ; 
	Nuemann_bound_number_pressure = 0 ; //BC.boundary_nodes_counter ;
	
	for(int i = 0 ; i < BC.boundary_nodes_counter ; i ++){
	
			if (BC.boundary_nodes[i].x_coordinate == L_Box){
				Dirichlet_bound_number_pressure ++ ;
			}
		}

	Nuemann_bound_number_pressure = BC.boundary_nodes_counter - Dirichlet_bound_number_pressure ; 

	Dirichlet_bound_number_velocity = Nuemann_bound_number_pressure ;
	Nuemann_bound_number_velocity = Dirichlet_bound_number_pressure ;
	
	Dirichlet_bound_pressure = new int[Dirichlet_bound_number_pressure] ; 
	Dirichlet_bound_velocity = new int[Dirichlet_bound_number_velocity] ;
	Nuemann_bound_pressure = new int[Nuemann_bound_number_pressure] ;
	Nuemann_bound_velocity = new int[Nuemann_bound_number_velocity] ;
	
	int counter = 0 ;
	int tally = 0 ; 

	for (int k = 0 ; k < BC.boundary_nodes_counter ; k ++){
			
			if (BC.boundary_nodes[k].x_coordinate == L_Box ){
				Dirichlet_bound_pressure[counter] = BC.boundary_nodes[k].node_number ;
				counter ++ ; continue ;
			}
				Nuemann_bound_pressure[tally] = BC.boundary_nodes[k].node_number ;
				tally ++ ;
		}

	for (int i = 0 ; i < Nuemann_bound_number_velocity ; i ++){
			
			Nuemann_bound_velocity[i] = Dirichlet_bound_pressure[i] ;		
		}
	
	for (int i = 0 ; i < Nuemann_bound_number_pressure ; i ++){
			
			Dirichlet_bound_velocity[i] = Nuemann_bound_pressure[i] ;	
		}
	
	//Determine Nodes on the cylinder
	Node_Cyl = 0 ;
	for (int i = 0 ; i < BC.boundary_nodes_counter ; i ++){
			
		double xx = BC.boundary_nodes[i].x_coordinate ;
		double yy = BC.boundary_nodes[i].y_coordinate ;
		
		double di = sqrt((xx-X_Center)*(xx-X_Center) + (yy-Y_Center)*(yy-Y_Center)) ;

		if (di > 1.50){ continue ;}
		Node_Cyl ++ ;
	}
	
	Cyl_Inf = new cylinder [Node_Cyl] ;
	int numbering = 0 ;
	for (int i = 0 ; i < BC.boundary_nodes_counter ; i ++){
		
		double xx = BC.boundary_nodes[i].x_coordinate ;
		double yy = BC.boundary_nodes[i].y_coordinate ;
		double di = sqrt((xx-X_Center)*(xx-X_Center) + (yy-Y_Center)*(yy-Y_Center)) ;

		if (di > 1.50){ continue ;}
		
		Cyl_Inf[numbering].x = xx ;
		Cyl_Inf[numbering].y = yy ;
		Cyl_Inf[numbering].node = BC.boundary_nodes[i].node_number ;
		numbering ++ ;
	}

	for (int i = 0 ; i < Node_Cyl ; i ++){

		if ((Cyl_Inf[i].x < X_Center) && (Cyl_Inf[i].y >= Y_Center)){
			
			Cyl_Inf[i].Si = (2.0/D)*(Cyl_Inf[i].y - Y_Center) ; 
			Cyl_Inf[i].Co = (2.0/D)*(X_Center - Cyl_Inf[i].x) ;
			Cyl_Inf[i].teta = asin(Cyl_Inf[i].Si) ;}

		if ((Cyl_Inf[i].x >= X_Center) && (Cyl_Inf[i].y > Y_Center)){

			Cyl_Inf[i].Si = (2.0/D)*(-Y_Center + Cyl_Inf[i].y) ; 
			Cyl_Inf[i].Co = (2.0/D)*(-(-X_Center + Cyl_Inf[i].x)) ;
			Cyl_Inf[i].teta = pi - asin(Cyl_Inf[i].Si) ;}

		if ((Cyl_Inf[i].x > X_Center) && (Cyl_Inf[i].y <= Y_Center)){
			
			Cyl_Inf[i].Si = (2.0/D)*(-(Y_Center - Cyl_Inf[i].y)) ; 
			Cyl_Inf[i].Co = (2.0/D)*(-(Cyl_Inf[i].x - X_Center)) ;
			Cyl_Inf[i].teta = -asin(Cyl_Inf[i].Si) + pi ;}

		if ((Cyl_Inf[i].x <= X_Center) && (Cyl_Inf[i].y < Y_Center)){
			
			Cyl_Inf[i].Si = (2.0/D)*(-(Y_Center - Cyl_Inf[i].y)) ; 
			Cyl_Inf[i].Co = (2.0/D)*(X_Center - Cyl_Inf[i].x) ;
			Cyl_Inf[i].teta = asin(Cyl_Inf[i].Si) + 2.0*pi ;}
	}
	
	sort(Cyl_Inf , Cyl_Inf + Node_Cyl , compareAngle) ;
	
	Cyl_Inf[Node_Cyl - 1].d_teta = 0.0 ;
	for (int i = 0 ; i < Node_Cyl - 1 ; i ++){
		Cyl_Inf[i].d_teta = Cyl_Inf[i+1].teta - Cyl_Inf[i].teta ; 
	}

	Cyl_Inf[Node_Cyl - 1].d_teta = Cyl_Inf[0].teta + (2.0*pi - Cyl_Inf[Node_Cyl - 1].teta) ;
	
	// Information
	ofstream fout ("OUT.txt") ;
	fout << Node_Cyl  << endl ;
	fout << endl ;
	for (int i = 0 ; i < Node_Cyl ; i ++){
		fout << "Node " << Cyl_Inf[i].node << endl ;
		fout << "x " << Cyl_Inf[i].x << endl ;
		fout << "y " << Cyl_Inf[i].y << endl ;
		fout << "d_teta "<< Cyl_Inf[i].d_teta << endl ;
		fout << "teta "<< Cyl_Inf[i].teta << endl ;
		fout << "Sin "<< Cyl_Inf[i].Si << endl ;
		fout << "Cos "<< Cyl_Inf[i].Co << endl ;
		fout << "-----------------------" << endl ;
	}
	// Information
	
	BC.Destroy() ;
	//BC.~Mesh() ;
	
	if (cycle == 1){ constructNuemannDirichletBoundsVkMinus1() ; }
	if (cycle == 2){ constructNuemannDirichletBoundsVkMinus2() ; }
	if (cycle == 3){ constructNuemannDirichletBoundsVkMinus3() ; }

	return ;
}

// -------------------- Construct Nuemann & Dirichlet Bounds V(k - 1) --------------------

void Cylinder :: constructNuemannDirichletBoundsVkMinus1(){
	
	Mesh BC("MeshP1C1.txt" , 1 , N_center) ;
	BC.generateMesh() ;
	
	Dirichlet_bound_number_pressure_Minus1 = 0 ; //BC.boundary_nodes_counter  ; 
	Nuemann_bound_number_pressure_Minus1 = 0 ; //BC.boundary_nodes_counter ;
	
	for(int i = 0 ; i < BC.boundary_nodes_counter ; i ++){
	
			if (BC.boundary_nodes[i].x_coordinate == L_Box ){
				// Dirichlet Pressure
				Dirichlet_bound_number_pressure_Minus1 ++ ;
			}
		}

	Nuemann_bound_number_pressure_Minus1 = BC.boundary_nodes_counter - Dirichlet_bound_number_pressure_Minus1 ; 

	Dirichlet_bound_pressure_Minus1 = new int[Dirichlet_bound_number_pressure_Minus1] ; 
	Nuemann_bound_pressure_Minus1 = new int[Nuemann_bound_number_pressure_Minus1] ;
	
	int counter = 0 ;
	int tally = 0 ; 

	for (int k = 0 ; k < BC.boundary_nodes_counter ; k ++){
			
			if (BC.boundary_nodes[k].x_coordinate == L_Box ){
				Dirichlet_bound_pressure_Minus1[counter] = BC.boundary_nodes[k].node_number ;
				counter ++ ; continue ;
			}
				Nuemann_bound_pressure_Minus1[tally] = BC.boundary_nodes[k].node_number ;
				tally ++ ;
		}
	
	BC.Destroy() ;
	//BC.~Mesh() ;
}

// -------------------- Construct Nuemann & Dirichlet Bounds V(k - 2) --------------------

void Cylinder :: constructNuemannDirichletBoundsVkMinus2(){
	
	Mesh BC("MeshP1C2.txt" , 1 , N_center) ;
	BC.generateMesh() ;
	
	Dirichlet_bound_number_pressure_Minus2 = 0 ; //BC.boundary_nodes_counter  ; 
	Nuemann_bound_number_pressure_Minus2 = 0 ; //BC.boundary_nodes_counter ;
	
	for(int i = 0 ; i < BC.boundary_nodes_counter ; i ++){
	
			if (BC.boundary_nodes[i].x_coordinate == L_Box ){
				// Dirichlet Pressure
				Dirichlet_bound_number_pressure_Minus2 ++ ;
			}
		}

	Nuemann_bound_number_pressure_Minus2 = BC.boundary_nodes_counter - Dirichlet_bound_number_pressure_Minus2 ; 

	Dirichlet_bound_pressure_Minus2 = new int[Dirichlet_bound_number_pressure_Minus2] ; 
	Nuemann_bound_pressure_Minus2 = new int[Nuemann_bound_number_pressure_Minus2] ;
	
	int counter = 0 ;
	int tally = 0 ; 

	for (int k = 0 ; k < BC.boundary_nodes_counter ; k ++){
			
			if (BC.boundary_nodes[k].x_coordinate == L_Box ){
				Dirichlet_bound_pressure_Minus2[counter] = BC.boundary_nodes[k].node_number ;
				counter ++ ; continue ;
			}
				Nuemann_bound_pressure_Minus2[tally] = BC.boundary_nodes[k].node_number ;
				tally ++ ;
		}
	
	BC.Destroy() ;
	//BC.~Mesh() ;
}

// -------------------- Construct Nuemann & Dirichlet Bounds V(k - 3) --------------------

void Cylinder :: constructNuemannDirichletBoundsVkMinus3(){
		
	Mesh BC("MeshP1C3.txt" , 1 , N_center) ;
	BC.generateMesh() ;
	
	Dirichlet_bound_number_pressure_Minus3 = 0 ; //BC.boundary_nodes_counter  ; 
	Nuemann_bound_number_pressure_Minus3 = 0 ; //BC.boundary_nodes_counter ;
	
	for(int i = 0 ; i < BC.boundary_nodes_counter ; i ++){
	
			if (BC.boundary_nodes[i].x_coordinate == L_Box ){
				// Dirichlet Pressure
				Dirichlet_bound_number_pressure_Minus3 ++ ;
			}
		}

	Nuemann_bound_number_pressure_Minus3 = BC.boundary_nodes_counter - Dirichlet_bound_number_pressure_Minus3 ; 

	Dirichlet_bound_pressure_Minus3 = new int[Dirichlet_bound_number_pressure_Minus3] ; 
	Nuemann_bound_pressure_Minus3 = new int[Nuemann_bound_number_pressure_Minus3] ;
	
	int counter = 0 ;
	int tally = 0 ; 

	for (int k = 0 ; k < BC.boundary_nodes_counter ; k ++){
			
			if (BC.boundary_nodes[k].x_coordinate == L_Box ){
				Dirichlet_bound_pressure_Minus3[counter] = BC.boundary_nodes[k].node_number ;
				counter ++ ; continue ;
			}
				Nuemann_bound_pressure_Minus3[tally] = BC.boundary_nodes[k].node_number ;
				tally ++ ;
		}
	
	BC.Destroy() ;
	//BC.~Mesh() ;
}

// -------------------- Distructor --------------------

Cylinder :: ~Cylinder(){

	cout << "The class is destroyed." << endl ;
}

//==================================
//==================================
//===== Out of Class functions =====
//==================================
//==================================

bool compareAngle(const cylinder & a , const cylinder & b){

   if (a.teta < b.teta) return true;
   return false ;
}
