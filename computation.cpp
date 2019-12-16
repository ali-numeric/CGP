// === Code Name: CGP/FEM  Project             
// === Author: Ali A. Kashefi																 	
// === Language: C++                                                                         

// === Class Name: Computation														  
// === computation.cpp

#include "computation.h"

bool compareMiddle(const middle & a , const middle & b);
bool compareDistance(const interspace & a , const interspace & b);
bool compareIindex(const AAIJ & a , const AAIJ & b);
int checkSparsity(AAIJ L[] , int J , int K , int Spicy , int NonZeroApp);

// -------------------- Constructor -------------------

Computation :: Computation(double Viscosity , double Density , double Time , double Time_Step , string s , double Alfa , double penalty_term , double Beta_q , int level , int n_center){
		
		mesh_type = true ;	
		
		cycle = level ;	
		N_center = n_center ;

		time = Time ;
		delta_t =  Time_Step ;

		time_step = (int)(time/delta_t) + 1 ;
		
		ifstream fin ("MeshP1.txt") ;
		fin >> pressure_node_number ;

		velocity_node_number = pressure_node_number ;
		
		if (mesh_type == true){
			
			ifstream fin ("MeshP1.txt") ;
			fin >> pressure_node_number ;
			fin.close() ;

			velocity_node_number = pressure_node_number ;
	
			ifstream finC1 ("MeshP1C1.txt") ;
			finC1 >> pressure_node_number_minus_1 ;
			finC1.close() ;

			ifstream finC2 ("MeshP1C2.txt") ;
			finC2 >> pressure_node_number_minus_2 ;
			finC2.close() ;

			ifstream finC3 ("MeshP1C3.txt") ;
			finC3 >> pressure_node_number_minus_3 ;
			finC3.close() ;
		}
		
		if (cycle == 1){velocity_coarse_number = pressure_node_number_minus_1 ;} // Just for P1-P1
		if (cycle == 2){velocity_coarse_number = pressure_node_number_minus_2 ;} // Just for P1-P1
		if (cycle == 3){velocity_coarse_number = pressure_node_number_minus_3 ;} // Just for P1-P1

		viscosity = Viscosity ; 
		density = Density ;
		r = penalty_term ;
		r2 = penalty_term ;
		beta_q = Beta_q  ;
		alfa = Alfa ;

		pressure_exp = new double[pressure_node_number] ;
		
		// We only store the last two steps.

		int memory = 2 ;

		velocity  = new double *[memory] ;
			for(int j = 0 ; j < memory ; j++ ){
				velocity [j] = new double [2 * velocity_node_number] ;
			}

		velocity_u  = new double *[memory] ;
			for(int j = 0 ; j < memory ; j++ ){
				velocity_u [j] = new double[velocity_node_number] ;
			}

		velocity_v  = new double *[memory] ;
			for(int j = 0 ; j < memory ; j++ ){
				velocity_v [j] = new double[velocity_node_number] ;
			}

		pressure  = new double *[memory] ;
			for(int j = 0 ; j < memory  ; j++ ){
				pressure [j] = new double[pressure_node_number] ;
			}
		
		RHS_Force = new double *[memory] ;
		for(int j = 0 ; j < memory ; j++ ){
				RHS_Force[j] = new double[2*velocity_node_number] ;
			}

		intermediate_velocity = new double [2 * velocity_node_number] ;
		
		if (cycle == 1){ intermediate_velocity_coarse_grid_L_1 = new double[2*velocity_coarse_number] ;}
		if (cycle == 2){ intermediate_velocity_coarse_grid_L_2 = new double[2*velocity_coarse_number] ;}
		if (cycle == 3){ intermediate_velocity_coarse_grid_L_3 = new double[2*velocity_coarse_number] ;}
}

// -------------------- Update RHS Pressure Part--------------------

void Computation :: updateRHSPressurePart(){ // current n

	// first update pressure
	
	updatePressureExp() ;
	productAb( Gradient_aa , Gradient_ai , Gradient_aj , pressure_exp , RHS_pressure_part , 1.0 , 2*velocity_node_number ) ; 
}

// -------------------- Update RHS Mass Part --------------------

void Computation :: updateRHSMassPart(){ //current n
		
	double sum = 0 ;

	for (int i = 0 ; i < 2*velocity_node_number ; i ++){

		for (int k = Mass_Velocity_aj[i] ; k < Mass_Velocity_aj[i+1] ; k ++ ){
       
				sum +=  Mass_Velocity_aa[k]*(1.0*velocity[concise/*current_step*/ - 1][Mass_Velocity_ai[k]]) ;
			}

		RHS_mass_part[i] = (1.0/delta_t)*sum  ;
		sum = 0 ;
	}
}

// -------------------- Update RHS --------------------

void Computation :: updateRHS(){ // current n

	// Pressure Part , Mass Part

	updateRHSPressurePart() ; // n is the current time step. In fact (n+1)
	updateRHSMassPart() ;	// n is the current time step. In fact (n+1)

	for (int i = 0 ; i < 2 * velocity_node_number ; i ++)
		RHS[i] = (RHS_mass_part[i]) ; //+RHS_Force[concise/*current_step*/][i] ;
}

// -------------------- Break RHS --------------------

void Computation :: breakRHS(){ // current n

	for (int i = 0 ; i < velocity_node_number ; i ++){

		RHS_U[i] = RHS[i] ; 
		RHS_V[i] = RHS[i + velocity_node_number] ; 
	}
}

// -------------------- Update Pressure --------------------

void Computation :: updatePressure(){

	for (int i = 0 ; i < pressure_node_number ; i ++)	
		pressure[concise/*current_step*/][i] = fi_fine_grid[i] ;
}

// -------------------- Updat Pressure Exp --------------------

void Computation :: updatePressureExp(){

	for (int i = 0 ; i < pressure_node_number ; i ++)
		pressure_exp[i] = 0.0 ;
}

// -------------------- Calculate RHS Poisson In Fine Grid --------------------

void Computation :: calculateRHSPoissonInFineGrid(){
	
	productAb( Divergence_aa , Divergence_ai , Divergence_aj , intermediate_velocity , RHS_Poisson_fine , (beta_q*density)/delta_t , pressure_node_number) ;
}

// -------------------- Calculate RHS Poisson In Coarse Grid V(k - 1) --------------------

void Computation :: calculateRHSPoissonInCoarseGridVkMinus1(){

	productAb( Divergence_V_k_minus_1_aa , Divergence_V_k_minus_1_ai , Divergence_V_k_minus_1_aj , intermediate_velocity_coarse_grid_L_1 , RHS_Poisson_coarse , (beta_q*density)/delta_t , pressure_node_number_minus_1) ;
}

// -------------------- Calculate RHS Poisson In Coarse Grid V(k - 2) --------------------

void Computation :: calculateRHSPoissonInCoarseGridVkMinus2(){
	
	productAb( Divergence_V_k_minus_2_aa , Divergence_V_k_minus_2_ai , Divergence_V_k_minus_2_aj , intermediate_velocity_coarse_grid_L_2 , RHS_Poisson_coarse , (beta_q*density)/delta_t , pressure_node_number_minus_2) ;
}

// -------------------- Calculate RHS Poisson In Coarse Grid V(k - 3) --------------------

void Computation :: calculateRHSPoissonInCoarseGridVkMinus3(){
	
	productAb( Divergence_V_k_minus_3_aa , Divergence_V_k_minus_3_ai , Divergence_V_k_minus_3_aj , intermediate_velocity_coarse_grid_L_3 , RHS_Poisson_coarse , (beta_q*density)/delta_t , pressure_node_number_minus_3) ;
}

// -------------------- Calculate Global Laplace Operator Matrix In Coarse Grid --------------------

void Computation :: calculateGlobalPressureLaplaceOperatorMatrixInCoarseGridPaper(int level){

//Level 1

	if (level == 1){

	int non_zero_approximation = 100*pressure_node_number_minus_1 ;
	AAIJ *Lp ;
	Lp = new AAIJ[non_zero_approximation] ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		Lp[i].aa = 0.0 ;
		Lp[i].i = 0 ;
		Lp[i].j = 0 ;}

	int spicy = 0 ;
	
	Mesh Vk_minus_1("MeshP1C1.txt" , 1 , N_center) ;
	Vk_minus_1.generateMesh() ;
		
	for (int i = 0 ; i < Vk_minus_1.number_of_elements ; i ++){

		int *node_array ;
		node_array = new int [3] ;

		Shape element(Vk_minus_1.domain[i]) ;
		element.setOperators() ;

		node_array[0] = Vk_minus_1.domain[i].node_1 ;
	    node_array[1] = Vk_minus_1.domain[i].node_2 ;
	    node_array[2] = Vk_minus_1.domain[i].node_3 ;

		for (int j = 0 ; j < 3 ; j ++){
			for (int k = 0 ; k < 3 ; k++){ 
				
				double value = 0.0 ;
				value =	element.getPressureLaplace ( element.mapping (node_array[j]), element.mapping (node_array[k])) ; 
				spicy = updateAAIJ(Lp , value , node_array[j] - 1 , node_array[k] - 1 , spicy , non_zero_approximation) ;
				
				}
			}
			
			element.Destroy() ;
			element.~Shape() ; 
			delete [] node_array ;
	}

	Vk_minus_1.Destroy() ;
	//Vk_minus_1.~Mesh() ;
	
	// Applying B.C.
	for (int j = 0 ; j < Dirichlet_bound_number_pressure_Minus1 ; j ++ ){
		for (int i = 0 ; i < spicy + 1 ; i ++ ){		
			if (Lp[i].i == Dirichlet_bound_pressure_Minus1[j] - 1){
				Lp[i].aa = 0.0 ;}
		}
	}
	
	for (int j = 0 ; j < Dirichlet_bound_number_pressure_Minus1 ; j ++){
		
			Lp[spicy].aa = 1.0 ;
			Lp[spicy].i = Dirichlet_bound_pressure_Minus1[j] - 1 ;
			Lp[spicy].j = Dirichlet_bound_pressure_Minus1[j] - 1 ;
			spicy ++ ;
		}

	int counter = 0 ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		if (Lp[i].aa != 0.0){counter ++ ;}
	}

	AAIJ *real_Lp ;
	real_Lp = new AAIJ[counter] ;
	int tally = 0 ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		if (Lp[i].aa != 0.0){

			real_Lp[tally].aa = Lp[i].aa ;
			real_Lp[tally].i = Lp[i].i ;
			real_Lp[tally].j = Lp[i].j ;
			tally ++ ;
		}
	}
	
	delete [] Lp ;

	Laplace_Pressure_V_k_minus_1_non_zeros_elements = counter ;
	Laplace_Pressure_V_k_minus_1_aa = new double [Laplace_Pressure_V_k_minus_1_non_zeros_elements] ;
	Laplace_Pressure_V_k_minus_1_ai = new int [Laplace_Pressure_V_k_minus_1_non_zeros_elements] ;
	Laplace_Pressure_V_k_minus_1_aj = new int [pressure_node_number_minus_1 + 1] ;

	convertAAIJtoCSR(pressure_node_number_minus_1 , Laplace_Pressure_V_k_minus_1_non_zeros_elements , real_Lp , Laplace_Pressure_V_k_minus_1_aa , Laplace_Pressure_V_k_minus_1_ai , Laplace_Pressure_V_k_minus_1_aj) ;
	
	delete [] real_Lp ;
	return ;

	}

//Level 2

	if (level == 2){
	
	int non_zero_approximation = 100*pressure_node_number_minus_2 ;
	AAIJ *Lp ;
	Lp = new AAIJ[non_zero_approximation] ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		Lp[i].aa = 0.0 ;
		Lp[i].i = 0 ;
		Lp[i].j = 0 ;}

	int spicy = 0 ;
	
	Mesh Vk_minus_2("MeshP1C2.txt" , 1 , N_center) ;
	Vk_minus_2.generateMesh() ;
		
	for (int i = 0 ; i < Vk_minus_2.number_of_elements ; i ++){

		int *node_array ;
		node_array = new int [3] ;

		Shape element(Vk_minus_2.domain[i]) ;
		element.setOperators() ;

		node_array[0] = Vk_minus_2.domain[i].node_1 ;
	    node_array[1] = Vk_minus_2.domain[i].node_2 ;
	    node_array[2] = Vk_minus_2.domain[i].node_3 ;

		for (int j = 0 ; j < 3 ; j ++){
			for (int k = 0 ; k < 3 ; k++){ 
				
				double value = 0.0 ;
				value =	element.getPressureLaplace ( element.mapping (node_array[j]), element.mapping (node_array[k])) ; 
				spicy = updateAAIJ(Lp , value , node_array[j] - 1 , node_array[k] - 1 , spicy , non_zero_approximation) ;
				
				}
			}
			
			element.Destroy() ;
			element.~Shape() ; 
			delete [] node_array ;
	}

	Vk_minus_2.Destroy() ;
	//Vk_minus_1.~Mesh() ;
	
	// Applying B.C.
	for (int j = 0 ; j < Dirichlet_bound_number_pressure_Minus2 ; j ++ ){
		for (int i = 0 ; i < spicy + 1 ; i ++ ){		
			if (Lp[i].i == Dirichlet_bound_pressure_Minus2[j] - 1){
				Lp[i].aa = 0.0 ;}
		}
	}
	
	for (int j = 0 ; j < Dirichlet_bound_number_pressure_Minus2 ; j ++){
		
			Lp[spicy].aa = 1.0 ;
			Lp[spicy].i = Dirichlet_bound_pressure_Minus2[j] - 1 ;
			Lp[spicy].j = Dirichlet_bound_pressure_Minus2[j] - 1 ;
			spicy ++ ;
		}

	int counter = 0 ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		if (Lp[i].aa != 0.0){counter ++ ;}
	}

	AAIJ *real_Lp ;
	real_Lp = new AAIJ[counter] ;
	int tally = 0 ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		if (Lp[i].aa != 0.0){

			real_Lp[tally].aa = Lp[i].aa ;
			real_Lp[tally].i = Lp[i].i ;
			real_Lp[tally].j = Lp[i].j ;
			tally ++ ;
		}
	}
	
	delete [] Lp ;

	Laplace_Pressure_V_k_minus_2_non_zeros_elements = counter ;
	Laplace_Pressure_V_k_minus_2_aa = new double [Laplace_Pressure_V_k_minus_2_non_zeros_elements] ;
	Laplace_Pressure_V_k_minus_2_ai = new int [Laplace_Pressure_V_k_minus_2_non_zeros_elements] ;
	Laplace_Pressure_V_k_minus_2_aj = new int [pressure_node_number_minus_2 + 1] ;

	convertAAIJtoCSR(pressure_node_number_minus_2 , Laplace_Pressure_V_k_minus_2_non_zeros_elements , real_Lp , Laplace_Pressure_V_k_minus_2_aa , Laplace_Pressure_V_k_minus_2_ai , Laplace_Pressure_V_k_minus_2_aj) ;
	
	delete [] real_Lp ;
	return ;
		
	}
	
//Level 3

	if (level == 3){

	int non_zero_approximation = 100*pressure_node_number_minus_3 ;
	AAIJ *Lp ;
	Lp = new AAIJ[non_zero_approximation] ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		Lp[i].aa = 0.0 ;
		Lp[i].i = 0 ;
		Lp[i].j = 0 ;}

	int spicy = 0 ;
	
	Mesh Vk_minus_3("MeshP1C3.txt" , 1 , N_center) ;
	Vk_minus_3.generateMesh() ;
		
	for (int i = 0 ; i < Vk_minus_3.number_of_elements ; i ++){

		int *node_array ;
		node_array = new int [3] ;

		Shape element(Vk_minus_3.domain[i]) ;
		element.setOperators() ;

		node_array[0] = Vk_minus_3.domain[i].node_1 ;
	    node_array[1] = Vk_minus_3.domain[i].node_2 ;
	    node_array[2] = Vk_minus_3.domain[i].node_3 ;

		for (int j = 0 ; j < 3 ; j ++){
			for (int k = 0 ; k < 3 ; k++){ 
				
				double value = 0.0 ;
				value =	element.getPressureLaplace ( element.mapping (node_array[j]), element.mapping (node_array[k])) ; 
				spicy = updateAAIJ(Lp , value , node_array[j] - 1 , node_array[k] - 1 , spicy , non_zero_approximation) ;
				
				}
			}
			
			element.Destroy() ;
			element.~Shape() ; 
			delete [] node_array ;
	}

	Vk_minus_3.Destroy() ;
	//Vk_minus_1.~Mesh() ;
	
	// Applying B.C.
	for (int j = 0 ; j < Dirichlet_bound_number_pressure_Minus3 ; j ++ ){
		for (int i = 0 ; i < spicy + 1 ; i ++ ){		
			if (Lp[i].i == Dirichlet_bound_pressure_Minus3[j] - 1){
				Lp[i].aa = 0.0 ;}
		}
	}
	
	for (int j = 0 ; j < Dirichlet_bound_number_pressure_Minus3 ; j ++){
		
			Lp[spicy].aa = 1.0 ;
			Lp[spicy].i = Dirichlet_bound_pressure_Minus3[j] - 1 ;
			Lp[spicy].j = Dirichlet_bound_pressure_Minus3[j] - 1 ;
			spicy ++ ;
		}

	int counter = 0 ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		if (Lp[i].aa != 0.0){counter ++ ;}
	}

	AAIJ *real_Lp ;
	real_Lp = new AAIJ[counter] ;
	int tally = 0 ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		if (Lp[i].aa != 0.0){

			real_Lp[tally].aa = Lp[i].aa ;
			real_Lp[tally].i = Lp[i].i ;
			real_Lp[tally].j = Lp[i].j ;
			tally ++ ;
		}
	}
	
	delete [] Lp ;

	Laplace_Pressure_V_k_minus_3_non_zeros_elements = counter ;
	Laplace_Pressure_V_k_minus_3_aa = new double [Laplace_Pressure_V_k_minus_3_non_zeros_elements] ;
	Laplace_Pressure_V_k_minus_3_ai = new int [Laplace_Pressure_V_k_minus_3_non_zeros_elements] ;
	Laplace_Pressure_V_k_minus_3_aj = new int [pressure_node_number_minus_3 + 1] ;

	convertAAIJtoCSR(pressure_node_number_minus_3 , Laplace_Pressure_V_k_minus_3_non_zeros_elements , real_Lp , Laplace_Pressure_V_k_minus_3_aa , Laplace_Pressure_V_k_minus_3_ai , Laplace_Pressure_V_k_minus_3_aj) ;
	
	delete [] real_Lp ;
	return ;

	}
}

// -------------------- Calculate Global Divergence Operator Matrix In Coarse Grid --------------------

void Computation :: calculateGlobalDivergenceOperatorMatrixInCoarseGridPaper(int level){

	if (level == 1){

	int velocity_coarse_node_number = pressure_node_number_minus_1 ; // It is true in case: P1-P1
	int non_zero_approximation_Di = 100*2*velocity_coarse_node_number ;
	AAIJ *Di ;

	Di = new AAIJ[non_zero_approximation_Di] ;
 
	for (int i = 0 ; i < non_zero_approximation_Di ; i ++){
		Di[i].aa = 0.0 ;
		Di[i].i = 0 ;
		Di[i].j = 0 ;}

	int spicy_Di = 0 ;
	
	Mesh Vk_minus_1("MeshP1C1.txt" , 1 , N_center) ;
	Vk_minus_1.generateMesh() ;
	N = Vk_minus_1.number_of_elements ;
		
	for (int i = 0 ; i < Vk_minus_1.number_of_elements ; i ++){

		int *node_array ;
		node_array = new int [3] ;

		Shape element(Vk_minus_1.domain[i]) ;
		element.setOperators() ;

		node_array[0] = Vk_minus_1.domain[i].node_1 ;
	    node_array[1] = Vk_minus_1.domain[i].node_2 ;
	    node_array[2] = Vk_minus_1.domain[i].node_3 ;

		for (int j = 0 ; j < 3 ; j ++){
			for (int k = 0 ; k < 3 ; k++){ 
				
				double value = 0.0 ;
				value = element.getH( element.mapping (node_array[j]), element.mapping(node_array[k])) ;
				spicy_Di = updateAAIJ(Di , value , node_array[j] - 1 , node_array[k] - 1 , spicy_Di , non_zero_approximation_Di) ;
				
				value = 0.0 ;
				value = element.getJ( element.mapping (node_array[j]), element.mapping(node_array[k])) ;
				spicy_Di = updateAAIJ(Di , value , node_array[j] - 1 , node_array[k] - 1 + velocity_coarse_node_number , spicy_Di , non_zero_approximation_Di) ;
			}
		}
			
			element.Destroy() ;
			element.~Shape() ; 
			delete [] node_array ;
	}

	Vk_minus_1.Destroy() ;
	//Vk.~Mesh() ;

	int counter = 0 ;
	for (int i = 0 ; i < non_zero_approximation_Di ; i ++){
		if (Di[i].aa != 0.0){counter ++ ;}
	}

	AAIJ *real_Di ;
	real_Di = new AAIJ[counter] ;
	int tally = 0 ;
	for (int i = 0 ; i < non_zero_approximation_Di ; i ++){
		if (Di[i].aa != 0.0){

			real_Di[tally].aa = Di[i].aa ;
			real_Di[tally].i = Di[i].i ;
			real_Di[tally].j = Di[i].j ;
			tally ++ ;
		}
	}
	
	delete [] Di ;
	
	Divergence_V_k_minus_1_non_zeros_elements = counter ;
	Divergence_V_k_minus_1_aa = new double [Divergence_V_k_minus_1_non_zeros_elements] ;
	Divergence_V_k_minus_1_ai = new int [Divergence_V_k_minus_1_non_zeros_elements] ;
	Divergence_V_k_minus_1_aj = new int [pressure_node_number_minus_1 + 1] ;

	convertAAIJtoCSR(pressure_node_number_minus_1 , Divergence_V_k_minus_1_non_zeros_elements , real_Di , Divergence_V_k_minus_1_aa , Divergence_V_k_minus_1_ai , Divergence_V_k_minus_1_aj) ;
	
	delete [] real_Di ;
	return ;

	}

	if (level == 2){
		
	int velocity_coarse_node_number = pressure_node_number_minus_2 ; // It is true in case: P1-P1
	int non_zero_approximation_Di = 100*2*velocity_coarse_node_number ;
	AAIJ *Di ;

	Di = new AAIJ[non_zero_approximation_Di] ;
 
	for (int i = 0 ; i < non_zero_approximation_Di ; i ++){
		Di[i].aa = 0.0 ;
		Di[i].i = 0 ;
		Di[i].j = 0 ;}

	int spicy_Di = 0 ;
	
	Mesh Vk_minus_2("MeshP1C2.txt" , 1 , N_center) ;
	Vk_minus_2.generateMesh() ;
	N = Vk_minus_2.number_of_elements ;
		
	for (int i = 0 ; i < Vk_minus_2.number_of_elements ; i ++){

		int *node_array ;
		node_array = new int [3] ;

		Shape element(Vk_minus_2.domain[i]) ;
		element.setOperators() ;

		node_array[0] = Vk_minus_2.domain[i].node_1 ;
	    node_array[1] = Vk_minus_2.domain[i].node_2 ;
	    node_array[2] = Vk_minus_2.domain[i].node_3 ;

		for (int j = 0 ; j < 3 ; j ++){
			for (int k = 0 ; k < 3 ; k++){ 
				
				double value = 0.0 ;
				value = element.getH( element.mapping (node_array[j]), element.mapping(node_array[k])) ;
				spicy_Di = updateAAIJ(Di , value , node_array[j] - 1 , node_array[k] - 1 , spicy_Di , non_zero_approximation_Di) ;
				
				value = 0.0 ;
				value = element.getJ( element.mapping (node_array[j]), element.mapping(node_array[k])) ;
				spicy_Di = updateAAIJ(Di , value , node_array[j] - 1 , node_array[k] - 1 + velocity_coarse_node_number , spicy_Di , non_zero_approximation_Di) ;
			}
		}
			
			element.Destroy() ;
			element.~Shape() ; 
			delete [] node_array ;
	}

	Vk_minus_2.Destroy() ;
	//Vk.~Mesh() ;

	int counter = 0 ;
	for (int i = 0 ; i < non_zero_approximation_Di ; i ++){
		if (Di[i].aa != 0.0){counter ++ ;}
	}

	AAIJ *real_Di ;
	real_Di = new AAIJ[counter] ;
	int tally = 0 ;
	for (int i = 0 ; i < non_zero_approximation_Di ; i ++){
		if (Di[i].aa != 0.0){

			real_Di[tally].aa = Di[i].aa ;
			real_Di[tally].i = Di[i].i ;
			real_Di[tally].j = Di[i].j ;
			tally ++ ;
		}
	}
	
	delete [] Di ;
	
	Divergence_V_k_minus_2_non_zeros_elements = counter ;
	Divergence_V_k_minus_2_aa = new double [Divergence_V_k_minus_2_non_zeros_elements] ;
	Divergence_V_k_minus_2_ai = new int [Divergence_V_k_minus_2_non_zeros_elements] ;
	Divergence_V_k_minus_2_aj = new int [pressure_node_number_minus_2 + 1] ;

	convertAAIJtoCSR(pressure_node_number_minus_2 , Divergence_V_k_minus_2_non_zeros_elements , real_Di , Divergence_V_k_minus_2_aa , Divergence_V_k_minus_2_ai , Divergence_V_k_minus_2_aj) ;
	
	delete [] real_Di ;
	return ;

	}

	if (level == 3){

	int velocity_coarse_node_number = pressure_node_number_minus_3 ; // It is true in case: P1-P1
	int non_zero_approximation_Di = 100*2*velocity_coarse_node_number ;
	AAIJ *Di ;

	Di = new AAIJ[non_zero_approximation_Di] ;
 
	for (int i = 0 ; i < non_zero_approximation_Di ; i ++){
		Di[i].aa = 0.0 ;
		Di[i].i = 0 ;
		Di[i].j = 0 ;}

	int spicy_Di = 0 ;
	
	Mesh Vk_minus_3("MeshP1C3.txt" , 1 , N_center) ;
	Vk_minus_3.generateMesh() ;
	N = Vk_minus_3.number_of_elements ;
		
	for (int i = 0 ; i < Vk_minus_3.number_of_elements ; i ++){

		int *node_array ;
		node_array = new int [3] ;

		Shape element(Vk_minus_3.domain[i]) ;
		element.setOperators() ;

		node_array[0] = Vk_minus_3.domain[i].node_1 ;
	    node_array[1] = Vk_minus_3.domain[i].node_2 ;
	    node_array[2] = Vk_minus_3.domain[i].node_3 ;

		for (int j = 0 ; j < 3 ; j ++){
			for (int k = 0 ; k < 3 ; k++){ 
				
				double value = 0.0 ;
				value = element.getH( element.mapping (node_array[j]), element.mapping(node_array[k])) ;
				spicy_Di = updateAAIJ(Di , value , node_array[j] - 1 , node_array[k] - 1 , spicy_Di , non_zero_approximation_Di) ;
				
				value = 0.0 ;
				value = element.getJ( element.mapping (node_array[j]), element.mapping(node_array[k])) ;
				spicy_Di = updateAAIJ(Di , value , node_array[j] - 1 , node_array[k] - 1 + velocity_coarse_node_number , spicy_Di , non_zero_approximation_Di) ;
			}
		}
			
			element.Destroy() ;
			element.~Shape() ; 
			delete [] node_array ;
	}

	Vk_minus_3.Destroy() ;
	//Vk.~Mesh() ;

	int counter = 0 ;
	for (int i = 0 ; i < non_zero_approximation_Di ; i ++){
		if (Di[i].aa != 0.0){counter ++ ;}
	}

	AAIJ *real_Di ;
	real_Di = new AAIJ[counter] ;
	int tally = 0 ;
	for (int i = 0 ; i < non_zero_approximation_Di ; i ++){
		if (Di[i].aa != 0.0){

			real_Di[tally].aa = Di[i].aa ;
			real_Di[tally].i = Di[i].i ;
			real_Di[tally].j = Di[i].j ;
			tally ++ ;
		}
	}
	
	delete [] Di ;
	
	Divergence_V_k_minus_3_non_zeros_elements = counter ;
	Divergence_V_k_minus_3_aa = new double [Divergence_V_k_minus_3_non_zeros_elements] ;
	Divergence_V_k_minus_3_ai = new int [Divergence_V_k_minus_3_non_zeros_elements] ;
	Divergence_V_k_minus_3_aj = new int [pressure_node_number_minus_3 + 1] ;

	convertAAIJtoCSR(pressure_node_number_minus_3 , Divergence_V_k_minus_3_non_zeros_elements , real_Di , Divergence_V_k_minus_3_aa , Divergence_V_k_minus_3_ai , Divergence_V_k_minus_3_aj) ;
	
	delete [] real_Di ;
	return ;

	}
}
	
// -------------------- Calculate RHS Real Velocity --------------------

void Computation :: calculateRHSRealVelocity(){ // Using in Eq3
	
	double sum1 = 0 ;
	double sum2 = 0 ; 

	for (int i = 0 ; i < 2 * velocity_node_number  ; i ++){

		for (int k = Mass_Velocity_aj[i] ; k < Mass_Velocity_aj[i+1] ; k ++ )
			sum1 +=  Mass_Velocity_aa[k]*intermediate_velocity[Mass_Velocity_ai[k]] ;
		
		for (int k = Gradient_aj[i] ; k < Gradient_aj[i+1] ; k ++ )
			sum2 +=  Gradient_aa[k]*fi_fine_grid[Gradient_ai[k]] ;

		RHS_real_velocity[i] = sum1 - (delta_t/beta_q)*sum2 ; 
		sum1 = 0 ;
		sum2 = 0 ;
	}
}

// -------------------- Apply Mixed BC to RHS of Poisson --------------------

void Computation :: applyMixedBCtoRHSofPoissonVk(){

	for (int j = 0 ; j < Dirichlet_bound_number_pressure ; j ++)
		RHS_Poisson_fine[Dirichlet_bound_pressure[j] - 1] = 0.0 ;
}

// -------------------- Apply Mixed BC to RHS of Poisson V(k - 1) --------------------
		
void Computation :: applyMixedBCtoRHSofPoissonVkMinus1(){
	
	for (int j = 0 ; j < Dirichlet_bound_number_pressure_Minus1 ; j ++)
		RHS_Poisson_coarse[Dirichlet_bound_pressure_Minus1[j] - 1] = 0.0 ;
}

// -------------------- Apply Mixed BC to RHS of Poisson V(k - 2) --------------------
		
void Computation :: applyMixedBCtoRHSofPoissonVkMinus2(){
	
	for (int j = 0 ; j < Dirichlet_bound_number_pressure_Minus2 ; j ++)
		RHS_Poisson_coarse[Dirichlet_bound_pressure_Minus2[j] - 1] = 0.0 ;
}

// -------------------- Apply Mixed BC to RHS of Poisson V(k - 3) --------------------

void Computation :: applyMixedBCtoRHSofPoissonVkMinus3(){
	
	for (int j = 0 ; j < Dirichlet_bound_number_pressure_Minus3 ; j ++)
		RHS_Poisson_coarse[Dirichlet_bound_pressure_Minus3[j] - 1] = 0.0 ;
}

// -------------------- Apply Mixed B.C. to the RHS of the Momentum Equation --------------------

void Computation :: applyMixedBCtoRHSofMomentum(){

	for (int j = 0 ; j < Dirichlet_bound_number_velocity ; j ++)
		RHS[Dirichlet_bound_velocity[j] - 1] = intermediate_u_BC[current_step][j] ;
		
	for (int j = 0 ; j < Dirichlet_bound_number_velocity ; j ++)
		RHS[(Dirichlet_bound_velocity[j] - 1) + velocity_node_number] = intermediate_v_BC[current_step][j] ; 
}

// -------------------- Solve Eq1 --------------------

void Computation :: solveEq1(int itr_max , int rpt_max){

	int maximum_iteration = itr_max ;
	int maximum_repeat = rpt_max ;

	double tol_abs;
	double tol_rel;

//  Set the initial solution estimation

    tol_abs = 1.0E-08;
    tol_rel = 1.0E-08;

	for (int i = 0 ; i < 2*velocity_node_number ; i ++)
		intermediate_velocity[i] = 0.0 ; 
		
    pmgmres_ilu_cr ( 2*velocity_node_number , LHS_non_zeros_elements , LHS_Velocity_aj , LHS_Velocity_ai , LHS_Velocity_aa , intermediate_velocity , RHS , maximum_iteration , maximum_repeat, tol_abs , tol_rel );
	
	return ;
}

// -------------------- Solve Eq1 Primary In X Direction --------------------

void Computation :: solveEq1primary(int itr_max , int rpt_max){

	int maximum_iteration = itr_max ;
	int maximum_repeat = rpt_max ;

	double tol_abs;
	double tol_rel;

//  Set the initial solution estimation

    tol_abs = 1.0E-08;
    tol_rel = 1.0E-08;

	for (int i = 0 ; i < velocity_node_number ; i ++)
		intermediateX[i] = 0.0 ;
	
    pmgmres_ilu_cr (velocity_node_number , LHS_non_zeros_elements , LHS_Velocity_aj , LHS_Velocity_ai , LHS_Velocity_aa , intermediateX , RHS_U , maximum_iteration , maximum_repeat, tol_abs , tol_rel );
	
	for (int i = 0 ; i < velocity_node_number ; i ++)
		intermediate_velocity[i] = intermediateX[i] ; 
	
	return ;
}

// -------------------- SolveEq1 Secondary In Y Direction --------------------

void Computation :: solveEq1secondary(int itr_max, int rpt_max){
	
	int maximum_iteration = itr_max ;
	int maximum_repeat = rpt_max ;

	double tol_abs;
	double tol_rel;

//  Set the initial solution estimation

    tol_abs = 1.0E-08;
    tol_rel = 1.0E-08;

	for (int i = 0 ; i < velocity_node_number ; i ++)
		intermediateX[i] = 1.0 ;
	
    pmgmres_ilu_cr (velocity_node_number , LHS_non_zeros_elements , LHS_Velocity_aj , LHS_Velocity_ai , LHS_Velocity_aa , intermediateX , RHS_V , maximum_iteration , maximum_repeat, tol_abs , tol_rel );
	
	for (int i = 0 ; i < velocity_node_number ; i ++)
		intermediate_velocity[i + velocity_node_number] = intermediateX[i] ; 
	
	return ;
}

// -------------------- Solve Eq2 --------------------

void Computation :: solveEq2 (int itr_max , int rpt_max){

	int maximum_iteration = itr_max ;
	int maximum_repeat = rpt_max ;

	double tol_abs;
	double tol_rel;

//  Set the initial solution estimate.
  
    tol_abs = 1.0E-08;
    tol_rel = 1.0E-08;

	for (int i = 0 ; i < pressure_node_number ; i ++)
		fi_fine_grid[i] = 0.0 ; 

    pmgmres_ilu_cr ( pressure_node_number , Laplace_Pressure_V_k_non_zeros_elements , Laplace_Pressure_V_k_aj , Laplace_Pressure_V_k_ai , Laplace_Pressure_V_k_aa , fi_fine_grid , RHS_Poisson_fine , maximum_iteration , maximum_repeat, tol_abs , tol_rel );
	
	return ;
}

// -------------------- solve Eq2 Vk Minus 1 --------------------

void Computation :: solveEq2VkMinus1(int itr_max, int rpt_max){
	
	int maximum_iteration = itr_max ;
	int maximum_repeat = rpt_max ;

	double tol_abs;
	double tol_rel;

//  Set the initial solution estimate.
  
    tol_abs = 1.0E-08;
    tol_rel = 1.0E-08;
	
	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++)
		fi_coarse_grid_L_1[i] = 0.0 ; 
	
    pmgmres_ilu_cr ( pressure_node_number_minus_1 , Laplace_Pressure_V_k_minus_1_non_zeros_elements , Laplace_Pressure_V_k_minus_1_aj , Laplace_Pressure_V_k_minus_1_ai , Laplace_Pressure_V_k_minus_1_aa , fi_coarse_grid_L_1 , RHS_Poisson_coarse , maximum_iteration , maximum_repeat, tol_abs , tol_rel );

	return ;
}

// -------------------- solve Eq2 Vk Minus 2 --------------------

void Computation :: solveEq2VkMinus2(int itr_max, int rpt_max){
	
	int maximum_iteration = itr_max ;
	int maximum_repeat = rpt_max ;

	double tol_abs;
	double tol_rel;
  
    tol_abs = 1.0E-08;
    tol_rel = 1.0E-08;

	for (int i = 0 ; i < pressure_node_number_minus_2 ; i ++)
		fi_coarse_grid_L_2[i] = 0.0 ; 

	pmgmres_ilu_cr ( pressure_node_number_minus_2 , Laplace_Pressure_V_k_minus_2_non_zeros_elements , Laplace_Pressure_V_k_minus_2_aj , Laplace_Pressure_V_k_minus_2_ai , Laplace_Pressure_V_k_minus_2_aa , fi_coarse_grid_L_2 , RHS_Poisson_coarse , maximum_iteration , maximum_repeat, tol_abs , tol_rel );
	
	return ;
}

// -------------------- solve Eq2 Vk Minus 3 --------------------

void Computation :: solveEq2VkMinus3(int itr_max, int rpt_max){
	
	int maximum_iteration = itr_max ;
	int maximum_repeat = rpt_max ;

	double tol_abs;
	double tol_rel;
  
    tol_abs = 1.0E-08;
    tol_rel = 1.0E-08;
	
	for (int i = 0 ; i < pressure_node_number_minus_3 ; i ++)
		fi_coarse_grid_L_3[i] = 0.0 ; 

    pmgmres_ilu_cr ( pressure_node_number_minus_3 , Laplace_Pressure_V_k_minus_3_non_zeros_elements , Laplace_Pressure_V_k_minus_3_aj , Laplace_Pressure_V_k_minus_3_ai , Laplace_Pressure_V_k_minus_3_aa , fi_coarse_grid_L_3 , RHS_Poisson_coarse , maximum_iteration , maximum_repeat, tol_abs , tol_rel );
		
	return ;
}

// -------------------- solve Eq3 --------------------

void Computation :: solveEq3(int itr_max , int rpt_max){
	
	int maximum_iteration = itr_max ;
	int maximum_repeat = rpt_max ;
   
	double tol_abs;
	double tol_rel;

    tol_abs = 1.0E-08;
    tol_rel = 1.0E-08;

	for (int i = 0 ; i < 2*velocity_node_number ; i ++ )
		U_GMRES[i] =  intermediate_velocity[i] ;

    pmgmres_ilu_cr ( 2*velocity_node_number , Mass_Prime_Velocity_non_zeros_elements , Mass_Prime_Velocity_aj , Mass_Prime_Velocity_ai , Mass_Prime_Velocity_aa , U_GMRES , RHS_real_velocity , maximum_iteration , maximum_repeat, tol_abs , tol_rel );

	return ;
}

// -------------------- Transfer GMRES to Velocity Vectors --------------------

void Computation :: transferGMREStoVectors(){
	
	for(int i = 0 ; i < 2*velocity_node_number ; i ++){
			velocity[concise/*current_step*/][i] = U_GMRES[i] ;
			
		if ( i < velocity_node_number ){ velocity_u[concise/*current_step*/][i] = velocity[concise/*current_step*/][i]  ;}
		if ( velocity_node_number <= i ){ velocity_v[concise/*current_step*/][i - velocity_node_number] = velocity[concise/*current_step*/][i]  ;}			 
		
	 }
}

// -------------------- Construct Global Gradient Operator Matrix Super --------------------

void Computation :: constructGlobalGradientOperatorMatrixSuper(){

	int non_zero_approximation_Gr = 12*2*velocity_node_number ;
	AAIJ *Gr ;

	Gr = new AAIJ[non_zero_approximation_Gr] ;
 
	for (int i = 0 ; i < non_zero_approximation_Gr ; i ++){
		Gr[i].aa = 0.0 ;
		Gr[i].i = 0 ;
		Gr[i].j = 0 ;}

	int spicy_Gr = 0 ;
	
	Mesh Vk("MeshP1.txt" , 1 , N_center) ;
	Vk.generateMesh() ;
	N = Vk.number_of_elements ;
		
	for (int i = 0 ; i < Vk.number_of_elements ; i ++){

		int *node_array ;
		node_array = new int [3] ;

		Shape element(Vk.domain[i]) ;
		element.setOperators() ;

		node_array[0] = Vk.domain[i].node_1 ;
	    node_array[1] = Vk.domain[i].node_2 ;
	    node_array[2] = Vk.domain[i].node_3 ;

		for (int j = 0 ; j < 3 ; j ++){
			for (int k = 0 ; k < 3 ; k++){ 
				
				double value = 0.0 ;
				value = element.getH( element.mapping (node_array[j]), element.mapping(node_array[k])) ;
				spicy_Gr = updateAAIJ(Gr , value , node_array[j] - 1 , node_array[k] - 1 , spicy_Gr , non_zero_approximation_Gr) ;
				
				value = 0.0 ;
				value = element.getJ( element.mapping (node_array[j]), element.mapping(node_array[k])) ;
				spicy_Gr = updateAAIJ(Gr , value , node_array[j] - 1 + velocity_node_number , node_array[k] - 1 , spicy_Gr , non_zero_approximation_Gr) ;
			}
		}
			
			element.Destroy() ;
			element.~Shape() ; 
			delete [] node_array ;
	}

	Vk.Destroy() ;
	//Vk.~Mesh() ;

	int counter = 0 ;
	for (int i = 0 ; i < non_zero_approximation_Gr ; i ++){
		if (Gr[i].aa != 0.0){counter ++ ;}
	}

	AAIJ *real_Gr ;
	real_Gr = new AAIJ[counter] ;
	int tally = 0 ;
	for (int i = 0 ; i < non_zero_approximation_Gr ; i ++){
		if (Gr[i].aa != 0.0){

			real_Gr[tally].aa = Gr[i].aa ;
			real_Gr[tally].i = Gr[i].i ;
			real_Gr[tally].j = Gr[i].j ;
			tally ++ ;
		}
	}
	
	delete [] Gr ;
	
	Gradient_non_zeros_elements = counter ;
	Gradient_aa = new double [Gradient_non_zeros_elements] ;
	Gradient_ai = new int [Gradient_non_zeros_elements] ;
	Gradient_aj = new int [2*velocity_node_number + 1] ;

	convertAAIJtoCSR(2*velocity_node_number , Gradient_non_zeros_elements , real_Gr , Gradient_aa , Gradient_ai , Gradient_aj) ;
	
	delete [] real_Gr ;
}

// -------------------- Construct Global Divergence Operator Matrix Super --------------------

void Computation :: constructGlobalDivergenceOperatorMatrixSuper(){
	
	int non_zero_approximation_Di = 100*2*velocity_node_number ;
	AAIJ *Di ;

	Di = new AAIJ[non_zero_approximation_Di] ;
 
	for (int i = 0 ; i < non_zero_approximation_Di ; i ++){
		Di[i].aa = 0.0 ;
		Di[i].i = 0 ;
		Di[i].j = 0 ;}

	int spicy_Di = 0 ;
	
	Mesh Vk("MeshP1.txt" , 1 , N_center) ;
	Vk.generateMesh() ;
	N = Vk.number_of_elements ;
		
	for (int i = 0 ; i < Vk.number_of_elements ; i ++){

		int *node_array ;
		node_array = new int [3] ;

		Shape element(Vk.domain[i]) ;
		element.setOperators() ;

		node_array[0] = Vk.domain[i].node_1 ;
	    node_array[1] = Vk.domain[i].node_2 ;
	    node_array[2] = Vk.domain[i].node_3 ;

		for (int j = 0 ; j < 3 ; j ++){
			for (int k = 0 ; k < 3 ; k++){ 
				
				double value = 0.0 ;
				value = element.getH( element.mapping (node_array[j]), element.mapping(node_array[k])) ;
				spicy_Di = updateAAIJ(Di , value , node_array[j] - 1 , node_array[k] - 1 , spicy_Di , non_zero_approximation_Di) ;
				
				value = 0.0 ;
				value = element.getJ( element.mapping (node_array[j]), element.mapping(node_array[k])) ;
				spicy_Di = updateAAIJ(Di , value , node_array[j] - 1 , node_array[k] - 1 + velocity_node_number , spicy_Di , non_zero_approximation_Di) ;
			}
		}
			
			element.Destroy() ;
			element.~Shape() ; 
			delete [] node_array ;
	}

	Vk.Destroy() ;
	//Vk.~Mesh() ;

	int counter = 0 ;
	for (int i = 0 ; i < non_zero_approximation_Di ; i ++){
		if (Di[i].aa != 0.0){counter ++ ;}
	}

	AAIJ *real_Di ;
	real_Di = new AAIJ[counter] ;
	int tally = 0 ;
	for (int i = 0 ; i < non_zero_approximation_Di ; i ++){
		if (Di[i].aa != 0.0){

			real_Di[tally].aa = Di[i].aa ;
			real_Di[tally].i = Di[i].i ;
			real_Di[tally].j = Di[i].j ;
			tally ++ ;
		}
	}
	
	delete [] Di ;
	
	Divergence_non_zeros_elements = counter ;
	Divergence_aa = new double [Divergence_non_zeros_elements] ;
	Divergence_ai = new int [Divergence_non_zeros_elements] ;
	Divergence_aj = new int [pressure_node_number + 1] ;
	
	convertAAIJtoCSR(pressure_node_number , Divergence_non_zeros_elements , real_Di , Divergence_aa , Divergence_ai , Divergence_aj) ;
	
	delete [] real_Di ;
}

// -------------------- Construct Laplace Pressure Vk Super --------------------

void Computation :: constructGlobalLaplacePressureVkSuper(){

	int non_zero_approximation = 100*pressure_node_number ;
	AAIJ *Lp ;
	Lp = new AAIJ[non_zero_approximation] ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		Lp[i].aa = 0.0 ;
		Lp[i].i = 0 ;
		Lp[i].j = 0 ;}

	int spicy = 0 ;
	
	Mesh Vk("MeshP1.txt" , 1 , N_center) ;
	Vk.generateMesh() ;
	N = Vk.number_of_elements ;
		
	for (int i = 0 ; i < Vk.number_of_elements ; i ++){

		int *node_array ;
		node_array = new int [3] ;

		Shape element(Vk.domain[i]) ;
		element.setOperators() ;

		node_array[0] = Vk.domain[i].node_1 ;
	    node_array[1] = Vk.domain[i].node_2 ;
	    node_array[2] = Vk.domain[i].node_3 ;

		for (int j = 0 ; j < 3 ; j ++){
			for (int k = 0 ; k < 3 ; k++){ 
				
				double value = 0.0 ;
				value =	element.getPressureLaplace ( element.mapping (node_array[j]), element.mapping (node_array[k])) ; 
				spicy = updateAAIJ(Lp , value , node_array[j] - 1 , node_array[k] - 1 , spicy , non_zero_approximation) ;
				
				}
			}
			
			element.Destroy() ;
			element.~Shape() ; 
			delete [] node_array ;
	}

	Vk.Destroy() ;
	
	// Applying B.C.
	for (int j = 0 ; j < Dirichlet_bound_number_pressure ; j ++ ){
		for (int i = 0 ; i < spicy + 1 ; i ++ ){		
			if (Lp[i].i == Dirichlet_bound_pressure[j] - 1){
				Lp[i].aa = 0.0 ;}
		}
	}
	
	for (int j = 0 ; j < Dirichlet_bound_number_pressure ; j ++){
		
			Lp[spicy].aa = 1.0 ;
			Lp[spicy].i = Dirichlet_bound_pressure[j] - 1 ;
			Lp[spicy].j = Dirichlet_bound_pressure[j] - 1 ;
			spicy ++ ;
		}

	int counter = 0 ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		if (Lp[i].aa != 0.0){counter ++ ;}
	}

	AAIJ *real_Lp ;
	real_Lp = new AAIJ[counter] ;
	int tally = 0 ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		if (Lp[i].aa != 0.0){

			real_Lp[tally].aa = Lp[i].aa ;
			real_Lp[tally].i = Lp[i].i ;
			real_Lp[tally].j = Lp[i].j ;
			tally ++ ;
		}
	}
	
	delete [] Lp ;

	Laplace_Pressure_V_k_non_zeros_elements = counter ;
	Laplace_Pressure_V_k_aa = new double [Laplace_Pressure_V_k_non_zeros_elements] ;
	Laplace_Pressure_V_k_ai = new int [Laplace_Pressure_V_k_non_zeros_elements] ;
	Laplace_Pressure_V_k_aj = new int [pressure_node_number + 1] ;

	convertAAIJtoCSR(pressure_node_number , Laplace_Pressure_V_k_non_zeros_elements , real_Lp , Laplace_Pressure_V_k_aa , Laplace_Pressure_V_k_ai , Laplace_Pressure_V_k_aj) ;
	
	delete [] real_Lp ;
}

// -------------------- Construct Velocity Mass Matrix Super --------------------

void Computation :: constructGlobalVelocityMassMatrixSuper(){
	
	int non_zero_approximation = 25*2*velocity_node_number ;
	AAIJ *Mv ;
	Mv = new AAIJ[non_zero_approximation] ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		Mv[i].aa = 0.0 ;
		Mv[i].i = 0 ;
		Mv[i].j = 0 ;}

	int spicy = 0 ;
	
	Mesh Vk("MeshP1.txt" , 1 , N_center);
	Vk.generateMesh() ;
	N = Vk.number_of_elements ;
		
	for (int i = 0 ; i < Vk.number_of_elements ; i ++){

		int *node_array ;
		node_array = new int [3] ;

		Shape element(Vk.domain[i]) ;
		element.setOperators() ;

		node_array[0] = Vk.domain[i].node_1 ;
	    node_array[1] = Vk.domain[i].node_2 ;
	    node_array[2] = Vk.domain[i].node_3 ;

		for (int j = 0 ; j < 3 ; j ++){
			for (int k = 0 ; k < 3 ; k++){ 
				 		
				double value = 0.0 ;
				value =	density*element.getVelocityMass( element.mapping (node_array[j]), element.mapping (node_array[k]))  ; 
				spicy = updateAAIJ(Mv , value , node_array[j] - 1 , node_array[k] - 1 , spicy , non_zero_approximation) ;
				spicy = updateAAIJ(Mv , value , node_array[j] - 1 + velocity_node_number , node_array[k] - 1 + velocity_node_number , spicy , non_zero_approximation) ;
				}
			}
			
			element.Destroy() ;
			element.~Shape() ; 
			delete [] node_array ;
	}

	Vk.Destroy() ;
	//Vk.~Mesh() ;

	int counter = 0 ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		if (Mv[i].aa != 0.0){counter ++ ;}
	}

	AAIJ *real_Mv ;
	real_Mv = new AAIJ[counter] ;
	int tally = 0 ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		if (Mv[i].aa != 0.0){

			real_Mv[tally].aa = Mv[i].aa ;
			real_Mv[tally].i = Mv[i].i ;
			real_Mv[tally].j = Mv[i].j ;
			tally ++ ;
		}
	}

	Mass_Velocity_non_zeros_elements = counter ; 
	Mass_Velocity_aa = new double [Mass_Velocity_non_zeros_elements] ;
	Mass_Velocity_ai = new int [Mass_Velocity_non_zeros_elements] ;
	Mass_Velocity_aj = new int [2*velocity_node_number + 1] ;
	
	convertAAIJtoCSR(2*velocity_node_number , Mass_Velocity_non_zeros_elements , real_Mv , Mass_Velocity_aa , Mass_Velocity_ai , Mass_Velocity_aj) ;

	delete [] real_Mv ;
	
	// Applying B.C.
	for (int j = 0 ; j < Dirichlet_bound_number_velocity ; j ++ ){
		for (int i = 0 ; i < spicy + 1 ; i ++ ){		
			if (Mv[i].i == Dirichlet_bound_velocity[j] - 1 || Mv[i].i == Dirichlet_bound_velocity[j] - 1 + velocity_node_number){
				Mv[i].aa = 0.0 ;}
		}
	}
	
	for (int j = 0 ; j < Nuemann_bound_number_velocity ; j ++ ){
		for (int i = 0 ; i < spicy + 1 ; i ++ ){		
			if (Mv[i].i == Nuemann_bound_velocity[j] - 1 || Mv[i].i == Nuemann_bound_velocity[j] - 1 + velocity_node_number){
				Mv[i].aa = 0.0 ;}
		}
	}

	for (int j = 0 ; j < Dirichlet_bound_number_velocity ; j ++){
		
			Mv[spicy].aa = 1.0 ;
			Mv[spicy].i = Dirichlet_bound_velocity[j] - 1 ;
			Mv[spicy].j = Dirichlet_bound_velocity[j] - 1 ;
			spicy ++ ;

			Mv[spicy].aa = 1.0 ;
			Mv[spicy].i = Dirichlet_bound_velocity[j] - 1 + velocity_node_number ;
			Mv[spicy].j = Dirichlet_bound_velocity[j] - 1 + velocity_node_number ;
			spicy ++ ;
		}

	for (int j = 0 ; j < Nuemann_bound_number_velocity ; j ++){

			Mv[spicy].aa = 1.0 ;
			Mv[spicy].i = Nuemann_bound_velocity[j] - 1  ;
			Mv[spicy].j = Nuemann_bound_velocity[j] - 1  ;
			spicy ++ ;

			Mv[spicy].aa = 1.0 ;
			Mv[spicy].i = Nuemann_bound_velocity[j] - 1 + velocity_node_number ;
			Mv[spicy].j = Nuemann_bound_velocity[j] - 1 + velocity_node_number ;
			spicy ++ ;
	
		}
	
	counter = 0 ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		if (Mv[i].aa != 0.0){counter ++ ;}
	}
	
	AAIJ *prime_Mv ;
	prime_Mv = new AAIJ[counter] ;
	tally = 0 ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		if (Mv[i].aa != 0.0){

			prime_Mv[tally].aa = Mv[i].aa ;
			prime_Mv[tally].i = Mv[i].i ;
			prime_Mv[tally].j = Mv[i].j ;
			tally ++ ;
		}
	}
	
	delete [] Mv ;

	Mass_Prime_Velocity_non_zeros_elements = counter ; 
	Mass_Prime_Velocity_aa = new double [Mass_Prime_Velocity_non_zeros_elements] ;
	Mass_Prime_Velocity_ai = new int [Mass_Prime_Velocity_non_zeros_elements] ;
	Mass_Prime_Velocity_aj = new int [2*velocity_node_number + 1] ;
	
	convertAAIJtoCSR(2*velocity_node_number , Mass_Prime_Velocity_non_zeros_elements , prime_Mv , Mass_Prime_Velocity_aa , Mass_Prime_Velocity_ai , Mass_Prime_Velocity_aj) ;
	
	delete [] prime_Mv ;
}

// -------------------- Apply Mixed BC To RHS of Correction Equ --------------------

void Computation :: applyMixedBCToRHSofCorrectionEqu(){
	
	for (int j = 0 ; j < Dirichlet_bound_number_velocity ; j ++)
		RHS_real_velocity[Dirichlet_bound_velocity[j] - 1] = intermediate_u_BC[current_step][j] ;	
		
	for (int j = 0 ; j < Dirichlet_bound_number_velocity ; j ++)
		RHS_real_velocity[(Dirichlet_bound_velocity[j] - 1) + velocity_node_number] = intermediate_v_BC[current_step][j] ; 
}

// -------------------- Construct Global Laplace Velocity Matrix Paper --------------------

void Computation :: constructGlobalLaplaceVelocityMatrixPaper(){

	int non_zero_approximation = 20*2*velocity_node_number ;
	AAIJ *Lhs ;
	Lhs = new AAIJ[non_zero_approximation] ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		Lhs[i].aa = 0.0 ;
		Lhs[i].i = 0 ;
		Lhs[i].j = 0 ;}

	int spicy = 0 ;

	Mesh Vk("MeshP1.txt" , 1 , N_center) ;
	Vk.generateMesh() ;
	N = Vk.number_of_elements ;
		
	for (int i = 0 ; i < Vk.number_of_elements ; i ++){

		int *node_array ;
		node_array = new int [3] ;

		Shape element(Vk.domain[i]) ;
		element.setOperators() ;

		node_array[0] = Vk.domain[i].node_1 ;
	    node_array[1] = Vk.domain[i].node_2 ;
	    node_array[2] = Vk.domain[i].node_3 ;

		for (int j = 0 ; j < 3 ; j ++){
			for (int k = 0 ; k < 3 ; k++){ 
				 	
					double value = 0.0 ;
					
					value = -viscosity*element.getVelocityLaplace( element.mapping (node_array[j]), element.mapping (node_array[k]))  ;
					value += (beta_q/delta_t)*density*element.getVelocityMass( element.mapping (node_array[j]), element.mapping (node_array[k]))  ; 

					spicy = updateAAIJ(Lhs , value , node_array[j] - 1 , node_array[k] - 1 , spicy , non_zero_approximation) ;
				}
			}
			
			element.Destroy() ;
			element.~Shape() ; 
			delete [] node_array ;
	}

	Vk.Destroy() ;
	//Vk.~Mesh() ;

	// Applying B.C.	
	for (int j = 0 ; j < Dirichlet_bound_number_velocity ; j ++ ){
		for (int i = 0 ; i < spicy + 1 ; i ++ ){		
			if (Lhs[i].i == Dirichlet_bound_velocity[j] - 1){
				Lhs[i].aa = 0.0 ;}
		}
	}
	
	for (int j = 0 ; j < Dirichlet_bound_number_velocity ; j ++){
		
			Lhs[spicy].aa = 1.0 ;
			Lhs[spicy].i = Dirichlet_bound_velocity[j] - 1 ;
			Lhs[spicy].j = Dirichlet_bound_velocity[j] - 1 ;
			spicy ++ ;
		}

	int counter = 0 ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		if (Lhs[i].aa != 0.0){counter ++ ;}
	}

	AAIJ *real_Lhs ;
	real_Lhs = new AAIJ[counter] ;
	int tally = 0 ;
	for (int i = 0 ; i < non_zero_approximation ; i ++){
		if (Lhs[i].aa != 0.0){

			real_Lhs[tally].aa = Lhs[i].aa ;
			real_Lhs[tally].i = Lhs[i].i ;
			real_Lhs[tally].j = Lhs[i].j ;
			tally ++ ;
		}
	}
	
	delete [] Lhs ;
	
	LHS_non_zeros_elements = counter ;
	LHS_Velocity_aa = new double [LHS_non_zeros_elements] ;
	LHS_Velocity_ai = new int [LHS_non_zeros_elements] ;
	LHS_Velocity_aj = new int [velocity_node_number + 1] ;

	convertAAIJtoCSR(velocity_node_number , LHS_non_zeros_elements , real_Lhs , LHS_Velocity_aa , LHS_Velocity_ai , LHS_Velocity_aj) ;
	
	delete [] real_Lhs ;

	Reserved_LHS_aa = new double[LHS_non_zeros_elements] ; 
	for(int i =0 ; i < LHS_non_zeros_elements ; i ++)
		Reserved_LHS_aa[i] = LHS_Velocity_aa[i] ; 
}

// -------------------- Construct Global Force Vector --------------------

void Computation :: constructGlobalForceVector(){

	for (int k = 0 ; k < concise + 1 /*time_step*/ ; k ++)
		productAb( Mass_Velocity_aa , Mass_Velocity_ai , Mass_Velocity_aj , RHS_Force[k] , RHS_Force[k] , 1.0/density , 2*velocity_node_number) ; 			
}

// -------------------- Construct Global Convective Matrix --------------------

void Computation :: constructGlobalConvectiveMatrix(){
	
	Mesh Vk("MeshP1.txt" , 1 , N_center) ;
	Vk.generateMesh() ;
	N = Vk.number_of_elements ;
	N_element = Vk.number_of_elements ;

	NonLinear = new convectionIntegration[N] ;
	Coefficient = new LocalGlobal[N] ;

	for (int i = 0 ; i < N ; i ++){
		NonLinear[i].createK() ;
		NonLinear[i].createL() ;
	}

	for (int i = 0 ; i < Vk.number_of_elements ; i ++){

		int *node_array ;
		node_array = new int [3] ;

		Shape element(Vk.domain[i]) ;
		element.setOperators() ;
		
		node_array[0] = Vk.domain[i].node_1 ;
	    node_array[1] = Vk.domain[i].node_2 ;
	    node_array[2] = Vk.domain[i].node_3 ;
		
		double *u_components ;
		double *v_components ;
		u_components = new double [3] ;
		v_components = new double [3] ;
		
		Coefficient[i].Local0 = element.mapping(node_array[0]) ;
		Coefficient[i].Global0 = node_array[0] ;
		Coefficient[i].Local1 = element.mapping(node_array[1]) ;
		Coefficient[i].Global1 = node_array[1] ;
		Coefficient[i].Local2 = element.mapping(node_array[2]) ;
		Coefficient[i].Global2 = node_array[2] ;

		for (int q = 0 ; q < 3 ; q ++){
			for (int l = 0 ; l < 3 ; l ++){
				for (int p = 0 ; p < 3 ; p ++){

					NonLinear[i].K[p][q][l] = element.getNonLinearU(element.mapping(node_array[p]), element.mapping(node_array[q]) , element.mapping(node_array[l])) ; 
					NonLinear[i].L[p][q][l] = element.getNonLinearV(element.mapping(node_array[p]), element.mapping(node_array[q]) , element.mapping(node_array[l])) ;
				}
			}
		}
		
		double nonlinear = 0.0 ; 
		
			element.Destroy() ;
			element.~Shape() ; 
			delete [] node_array ;
			delete [] u_components ;
			delete [] v_components ;	
	}

	Vk.Destroy() ;
	//Vk.~Mesh() ;
}

// -------------------- Construct Nonlinear Paper --------------------

void Computation :: constructNonlinearPaper(){

	non_zero_convection = LHS_non_zeros_elements ; 

	aN = new NonlinearSuper[non_zero_convection] ;
	for (int i = 0 ; i < non_zero_convection ; i ++){ aN[i].ak = 0 ; aN[i].counter = 0 ; aN[i].boundary = true ;}

	// PART 2
	for (int i = 0 ; i < N_element ; i ++){

		int *node_array ;
		node_array = new int [3] ;

		node_array[0] = Coefficient[i].Global0 ;
		node_array[1] = Coefficient[i].Global1 ;
		node_array[2] = Coefficient[i].Global2 ;

		for (int j = 0 ; j < 3 ; j ++){
			for (int k = 0 ; k < 3 ; k++){ 

					if (checkBoundaryPoint(node_array[j] - 1)){	
						continue ;
					}
					
					int s = searchCSR(node_array[j] - 1 , node_array[k] - 1 , LHS_Velocity_ai , LHS_Velocity_aj , LHS_non_zeros_elements) ;
					aN[s].ak ++ ;
					aN[s].boundary = false ;

				}
			}
			
			delete [] node_array ;
	}	

	// PART 3

	for (int i = 0 ; i < non_zero_convection ; i ++){ aN[i].create() ;}

	// PART 4
	
	for (int i = 0 ; i < N_element ; i ++){

		int *node_array ;
		node_array = new int [3] ;

		node_array[0] = Coefficient[i].Global0 ;
		node_array[1] = Coefficient[i].Global1 ;
		node_array[2] = Coefficient[i].Global2 ;

		for (int j = 0 ; j < 3 ; j ++){
			for (int k = 0 ; k < 3 ; k++){ 

					if (checkBoundaryPoint(node_array[j] - 1)){continue ;}
					
					int s = searchCSR(node_array[j] - 1 , node_array[k] - 1 , LHS_Velocity_ai , LHS_Velocity_aj , LHS_non_zeros_elements) ;
					aN[s].f[aN[s].counter] = i ; // element number 
					aN[s].X[aN[s].counter] = j ; // J index
					aN[s].Z[aN[s].counter] = k ; // K index
					aN[s].counter ++ ;

				}
			}
			
			delete [] node_array ;
	}
}

// -------------------- Update Convection Advanced Super --------------------

void Computation :: updateConvectionAdvancedSuper(double Virtual_Velocity[]){

	for (int i = 0 ; i < non_zero_convection ; i ++){
		if (aN[i].boundary == true){continue ;}

		for (int k = 0 ; k < aN[i].ak ; k ++){
		
		int *node_array ;
		node_array = new int [3] ;

		node_array[0] = Coefficient[aN[i].f[k]].Global0 ;
		node_array[1] = Coefficient[aN[i].f[k]].Global1 ;
		node_array[2] = Coefficient[aN[i].f[k]].Global2 ;

		double *u_components ;
		double *v_components ;
		u_components = new double [3] ;
		v_components = new double [3] ;
		
		u_components[Coefficient[aN[i].f[k]].Local0] = Virtual_Velocity[Coefficient[aN[i].f[k]].Global0 - 1] ;
		u_components[Coefficient[aN[i].f[k]].Local1] = Virtual_Velocity[Coefficient[aN[i].f[k]].Global1 - 1] ;
		u_components[Coefficient[aN[i].f[k]].Local2] = Virtual_Velocity[Coefficient[aN[i].f[k]].Global2 - 1] ;

		v_components[Coefficient[aN[i].f[k]].Local0] = Virtual_Velocity[Coefficient[aN[i].f[k]].Global0 - 1 + velocity_node_number] ; 
		v_components[Coefficient[aN[i].f[k]].Local1] = Virtual_Velocity[Coefficient[aN[i].f[k]].Global1 - 1 + velocity_node_number] ;
		v_components[Coefficient[aN[i].f[k]].Local2] = Virtual_Velocity[Coefficient[aN[i].f[k]].Global2 - 1 + velocity_node_number] ;
		
		double nonlinear = 0.0 ; 

		for (int p = 0 ; p < 3 ; p ++){
					 
			nonlinear += (u_components[p]*NonLinear[aN[i].f[k]].K[p][aN[i].X[k]][aN[i].Z[k]] + v_components[p]*NonLinear[aN[i].f[k]].L[p][aN[i].X[k]][aN[i].Z[k]])  ;}

			LHS_Velocity_aa[i] += nonlinear ;
			nonlinear = 0 ;
			
			delete [] node_array ;
			delete [] u_components ;
			delete [] v_components ;
 		}
	}
}

// -------------------- Check Boundary Points --------------------

bool Computation :: checkBoundaryPoint(int i){

	for (int j = 0 ; j < Dirichlet_bound_number_velocity ; j ++){
		
		if(i == Dirichlet_bound_velocity[j] - 1){return true ;}
	}

	return false ;
}

// -------------------- Go To Convection Advanced --------------------

void Computation :: goToConvectionAdvanced(){ 
	
	for (int i = 0 ; i < LHS_non_zeros_elements ; i ++)
		LHS_Velocity_aa[i] = Reserved_LHS_aa[i] ;

	updateConvectionAdvancedSuper(velocity[concise/*current_step*/ - 1]) ;
}

// -------------------- Construct The Rest Of The Matrices & Vectors --------------------

void Computation :: constructTheRestOfTheMatricesVectors(){

		fi_fine_grid = new double [pressure_node_number] ;

		for (int i = 0 ; i < pressure_node_number ; i ++)	
			fi_fine_grid[i] = 0.0 ; 

		if (cycle == 1){

			fi_coarse_grid_L_1 = new double [pressure_node_number_minus_1 /*(int)pow(coarsen_number/2.0 + 1.0 , 2.0)*/] ;
			for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
				fi_coarse_grid_L_1[i] = 0.0 ;
			}
		}
		
		if (cycle == 2){

			fi_coarse_grid_L_1 = new double [pressure_node_number_minus_1 /*(int)pow(coarsen_number/2.0 + 1.0 , 2.0)*/] ;
			fi_coarse_grid_L_2 = new double [pressure_node_number_minus_2 /*(int)pow(coarsen_number/4.0 + 1.0 , 2.0)*/] ;
			for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
				fi_coarse_grid_L_1[i] = 0.0 ;
			}

			for (int i = 0 ; i < pressure_node_number_minus_2 ; i ++)
				fi_coarse_grid_L_2[i] = 0.0 ;
		}

		if (cycle == 3){

			fi_coarse_grid_L_1 = new double [pressure_node_number_minus_1/*(int)pow(coarsen_number/2.0 + 1.0 , 2.0)*/] ;
			fi_coarse_grid_L_2 = new double [pressure_node_number_minus_2/*(int)pow(coarsen_number/4.0 + 1.0 , 2.0)*/] ;
			fi_coarse_grid_L_3 = new double [pressure_node_number_minus_3/*(int)pow(coarsen_number/8.0 + 1.0 , 2.0)*/] ;
			for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
				fi_coarse_grid_L_1[i] = 0.0 ;
			}

			for (int i = 0 ; i < pressure_node_number_minus_2 ; i ++)
				fi_coarse_grid_L_2[i] = 0.0 ;

			for (int i = 0 ; i < pressure_node_number_minus_3 ; i ++)
				fi_coarse_grid_L_3[i] = 0.0 ;
		}
		
		RHS_Poisson_fine = new double [pressure_node_number] ; 
		RHS_real_velocity = new double [2 * velocity_node_number] ;
		
		double l = pow(2.0 , (double)cycle) ;
		if (mesh_type == false) {RHS_Poisson_coarse = new double [(int)(pow(coarsen_number/l + 1.0 , 2.0))] ;}
		if (mesh_type == true) {
			if (cycle == 1){RHS_Poisson_coarse = new double [pressure_node_number_minus_1] ;}
			if (cycle == 2){RHS_Poisson_coarse = new double [pressure_node_number_minus_2] ;}
			if (cycle == 3){RHS_Poisson_coarse = new double [pressure_node_number_minus_3] ;}
		}

		U_GMRES = new double [2*velocity_node_number] ;

		RHS = new double [2 * velocity_node_number] ;
		RHS_U = new double [velocity_node_number] ;
		RHS_V = new double [velocity_node_number] ;

		RHS_pressure_part = new double [2 * velocity_node_number] ;
		RHS_mass_part = new double [2 * velocity_node_number] ; 

		U = new double[velocity_node_number] ;
		V = new double[velocity_node_number] ;

		intermediateX = new double[velocity_node_number] ;
		stress = new double [2*velocity_node_number] ;
		Lift = new double [time_step] ;
		Drag = new double [time_step] ;
		PressureDrag = new double [time_step] ;
		ViscousDrag = new double [time_step] ;
		PressureLift = new double [time_step] ;
		ViscousLift = new double [time_step] ;
}

// -------------------- Solve The Problem In Advanced Version Paper --------------------

void Computation :: solveTheProblemAdvancedPaper(int level){
	
	if (cycle == 0){solveTheProblemAdvancedVkPaper() ; return ;}
	if (cycle == 1){solveTheProblemAdvancedVkMinus1Paper() ; return ;}
	if (cycle == 2){solveTheProblemAdvancedVkMinus2Paper() ; return ;}
	if (cycle == 3){solveTheProblemAdvancedVkMinus3Paper() ; return ;}
}

// -------------------- solve The Problem Advanced Vk Paper --------------------

void Computation :: solveTheProblemAdvancedVkPaper(){

	StopWatch preprocessing2(1) ;

	int itr_max1 = 2*velocity_node_number - 5 ; 
	int rpt_max1 = 20 ;

	int itr_max2 = pressure_node_number - 2 ;  
	int rpt_max2 = 20 ;

	int itr_max3 = 2*velocity_node_number - 5 ; 
	int rpt_max3 = 20 ;

	//** Preprocessing **//

	preprocessing2.startTime() ;
	constructGlobalConvectiveMatrix() ;
	preprocessing2.stopTime() ;
	preprocessing2.resetTime() ;
	
	preprocessing2.startTime() ;
	constructGlobalGradientOperatorMatrixSuper() ;
	preprocessing2.stopTime() ;
	preprocessing2.resetTime() ;
	
	preprocessing2.startTime() ;
	constructGlobalDivergenceOperatorMatrixSuper() ;
	preprocessing2.stopTime() ;
	preprocessing2.resetTime() ;
	
	preprocessing2.startTime() ;
	constructGlobalLaplacePressureVkSuper() ;
	preprocessing2.stopTime() ;
	preprocessing2.resetTime() ;
	
	preprocessing2.startTime() ;
	constructGlobalVelocityMassMatrixSuper() ;  
	preprocessing2.stopTime() ;
	preprocessing2.resetTime() ;
	
	preprocessing2.startTime() ;
	constructGlobalLaplaceVelocityMatrixPaper() ; 
	preprocessing2.stopTime() ;
	preprocessing2.resetTime() ;
	
	preprocessing2.startTime() ;
	constructTheRestOfTheMatricesVectors() ;
	constructNonlinearPaper() ; 
	preprocessing2.stopTime() ;
	preprocessing2.resetTime() ;
	preprocess2 = preprocessing2.finalTime() ;
	 
	// Processing //
	// Start the iterations

	StopWatch Adv(1) ;
	StopWatch Poisson(1) ;

	current_step = 1 ;
	concise = 1 ;

	while( current_step < time_step ){

		cout << "Current Step: " << current_step << endl ; cout << endl ;
		
		goToConvectionAdvanced() ;
		updateRHS() ; 
		applyMixedBCtoRHSofMomentum() ; 
		breakRHS() ;

		//Eq1
		
		Adv.startTime() ;
		solveEq1primary(itr_max1/2 , rpt_max1/2) ;
		solveEq1secondary(itr_max1/2 , rpt_max1/2) ;
		Adv.stopTime() ;
		Adv.resetTime() ;
		
		//Eq2

		calculateRHSPoissonInFineGrid() ;
		applyMixedBCtoRHSofPoissonVk() ;

		Poisson.startTime() ;

		solveEq2(itr_max2 , rpt_max2) ;
		
		Poisson.stopTime() ;
		Poisson.resetTime() ;
		
		//Eq3

		calculateRHSRealVelocity() ;
		applyMixedBCToRHSofCorrectionEqu() ;
		
		solveEq3(itr_max3 , rpt_max3) ;
		
		transferGMREStoVectors() ;

		//Eq4
		updatePressure() ;

//		computeDragAndLift() ;

		exchangeVariables() ;
		current_step ++ ;
		updatePressureExp() ;
		concise = 1 ; // current_step ; 
	}

	advection = Adv.finalTime() ; 
	poissonic = Poisson.finalTime() ;

	cout << "The job is finished." << endl ;

	return ;

	// current_step ; // physcially n in formulas, computationally n - 1
}

// -------------------- solve The Problem Advanced V(k - 1) Paper --------------------

void Computation :: solveTheProblemAdvancedVkMinus1Paper(){

	StopWatch preprocessing2(1) ;

	int itr_max1 = 2*velocity_node_number - 5 ; 
	int rpt_max1 = 20 ;

	int itr_max2_Vk_1 = pressure_node_number_minus_1 - 2 ; 
	int rpt_max2_Vk_1 = 20 ;

	int itr_max3 = 2*velocity_node_number - 5 ; 
	int rpt_max3 = 20 ;

	//** Preprocessing **//

	preprocessing2.startTime() ;
	constructGlobalConvectiveMatrix() ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;

	preprocessing2.startTime() ;
	constructGlobalGradientOperatorMatrixSuper() ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;
	
	preprocessing2.startTime() ;
	constructGlobalVelocityMassMatrixSuper() ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;
	
	preprocessing2.startTime() ;
	constructGlobalLaplaceVelocityMatrixPaper() ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;
	
	preprocessing2.startTime() ;
	constructTheRestOfTheMatricesVectors() ;
	constructNonlinearPaper() ; 
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;

	preprocessing2.startTime() ;
	// Set Multi-grid Tools //
	StopWatch MG(1) ;
	MG.startTime() ;

	if (mesh_type == true){constructMapFineCoarseNonUniformGrid() ;}
	
	MG.stopTime() ;
    MG.resetTime() ;
	multigrid = MG.finalTime() ;
	
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;

	preprocessing2.startTime() ;
	calculateGlobalPressureLaplaceOperatorMatrixInCoarseGridPaper(cycle) ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;

	preprocessing2.startTime() ;
	calculateGlobalDivergenceOperatorMatrixInCoarseGridPaper(cycle) ;			
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;
	preprocess2 = preprocessing2.finalTime() ;
	
	// Processing //
	// Start the iterations

	StopWatch Adv(1) ;
	StopWatch Poisson(1) ;
	StopWatch Inj(1) ;
	StopWatch Pro(1) ;

	current_step = 1 ;
	concise = 1 ;

	while( current_step < time_step ){

		cout << "Current Step: " << current_step << endl ; cout << endl ;

		goToConvectionAdvanced() ;	
		updateRHS() ;
		applyMixedBCtoRHSofMomentum() ;
		breakRHS() ;

		//Eq1

		Adv.startTime() ;
		solveEq1primary(itr_max1/2 , rpt_max1/2) ;
		solveEq1secondary(itr_max1/2 , rpt_max1/2) ;
		Adv.stopTime() ;
		Adv.resetTime() ;

		// Multigrid //

		// Injection //

		Inj.startTime() ;

		injectIntermediateVelocityVktoVkMinus1() ;
		
		Inj.stopTime() ;
		Inj.resetTime() ;

		//Eq2

		calculateRHSPoissonInCoarseGridVkMinus1() ;
		applyMixedBCtoRHSofPoissonVkMinus1() ;

		Poisson.startTime() ;

		solveEq2VkMinus1(itr_max2_Vk_1 , itr_max2_Vk_1) ;
		
		Poisson.stopTime() ;
		Poisson.resetTime() ;

		// Prolongation //

		Pro.startTime() ;

		prolongationPressureVkMinus1toVk() ; 
		
		Pro.stopTime() ;
		Pro.resetTime() ;

		//Eq3
		
		calculateRHSRealVelocity() ;
		applyMixedBCToRHSofCorrectionEqu() ;
		
		solveEq3(itr_max3 , rpt_max3) ;

		transferGMREStoVectors() ;

		//Eq4
		updatePressure() ;
		
		// computeDragAndLift() ;

		exchangeVariables() ;
		current_step ++ ;
		updatePressureExp() ;
		concise = 1 ; // current_step ;
	}

	advection = Adv.finalTime() ; 
	poissonic = Poisson.finalTime() ;
	inject = Inj.finalTime() ; 
	prolong = Pro.finalTime() ;
	multigrid = MG.finalTime() ;

	cout << " The job is finished. " << endl ;

	return ;

	// current_step ; // physcially n in formulas, computationally n - 1
}

// -------------------- solve The Problem Advanced V(k - 2) Paper --------------------

void Computation :: solveTheProblemAdvancedVkMinus2Paper(){

	int itr_max1 = 2*velocity_node_number - 5 ; 
	int rpt_max1 = 20 ;

	int itr_max2_Vk_2 = pressure_node_number_minus_2 - 2 ; 
	int rpt_max2_Vk_2 = 20 ;

	int itr_max3 = 2*velocity_node_number - 5 ; 
	int rpt_max3 = 20 ;

	//** Preprocessing **//

	StopWatch preprocessing2(1) ;
	preprocessing2.startTime() ;
	constructGlobalConvectiveMatrix() ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;
	
	preprocessing2.startTime() ;
	constructGlobalGradientOperatorMatrixSuper() ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;

	preprocessing2.startTime() ;
	constructGlobalVelocityMassMatrixSuper() ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;
	
	preprocessing2.startTime() ;
	constructGlobalLaplaceVelocityMatrixPaper() ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;
	
	preprocessing2.startTime() ;
	constructTheRestOfTheMatricesVectors() ;
	constructNonlinearPaper() ; 
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;

	preprocessing2.startTime() ;
	// Set Multi-grid Tools //
	StopWatch MG(1) ;
	MG.startTime() ;

	if (mesh_type == true)  {constructMapFineCoarseNonUniformGrid() ;}
	
	MG.stopTime() ;
    MG.resetTime() ;
	multigrid = MG.finalTime() ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;

	preprocessing2.startTime() ;
	calculateGlobalPressureLaplaceOperatorMatrixInCoarseGridPaper(cycle) ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;

	preprocessing2.startTime() ;
	calculateGlobalDivergenceOperatorMatrixInCoarseGridPaper(cycle) ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;
	preprocess2 = preprocessing2.finalTime() ;
	
	// Processing //
	// Start the iterations

	StopWatch Adv(1) ;
	StopWatch Poisson(1) ;
	StopWatch Inj(1) ;
	StopWatch Pro(1) ;

	current_step = 1 ;
	concise = 1 ;

	while( current_step < time_step ){

		cout << "Current Step: " << current_step << endl ; cout << endl ;

		goToConvectionAdvanced() ;	
		updateRHS() ;
		applyMixedBCtoRHSofMomentum() ;
		breakRHS() ;

		//Eq1

		Adv.startTime() ;
		solveEq1primary(itr_max1/2 , rpt_max1/2) ;
		solveEq1secondary(itr_max1/2 , rpt_max1/2) ;
		Adv.stopTime() ;
		Adv.resetTime() ;

		// Multigrid //

		// Injection //

		Inj.startTime() ;

		injectIntermediateVelocityVktoVkMinus2() ;
		
		Inj.stopTime() ;
		Inj.resetTime() ;

		//Eq2

		calculateRHSPoissonInCoarseGridVkMinus2() ;
		applyMixedBCtoRHSofPoissonVkMinus2() ;

		Poisson.startTime() ;

		solveEq2VkMinus2(itr_max2_Vk_2 , itr_max2_Vk_2) ;
		
		Poisson.stopTime() ;
		Poisson.resetTime() ;

		// Prolongation //

		Pro.startTime() ;

		prolongationPressureVkMinus2toVk() ; 
		
		Pro.stopTime() ;
		Pro.resetTime() ;

		//Eq3
		
		calculateRHSRealVelocity() ;
		applyMixedBCToRHSofCorrectionEqu() ;

		solveEq3(itr_max3 , rpt_max3) ;

		transferGMREStoVectors() ;

		//Eq4
		updatePressure() ;
			
		// computeDragAndLift() ;

		exchangeVariables() ;
		current_step ++ ;
		updatePressureExp() ;
		concise = 1 ; // current_step ;
	}

	advection = Adv.finalTime() ; 
	poissonic = Poisson.finalTime() ;
	inject = Inj.finalTime() ; 
	prolong = Pro.finalTime() ;
	multigrid = MG.finalTime() ;

	cout << " The job is finished. " << endl ;

	return ;

	// current_step ; // physcially n in formulas, computationally n - 1
}

// -------------------- solve The Problem Advanced V(k - 3) Paper --------------------

void Computation :: solveTheProblemAdvancedVkMinus3Paper(){

	int itr_max1 = 2*velocity_node_number - 5 ; 
	int rpt_max1 = 20 ;

	int itr_max2_Vk_3 = pressure_node_number_minus_3 - 2 ; 
	int rpt_max2_Vk_3 = 20 ;

	int itr_max3 = 2*velocity_node_number - 5 ; 
	int rpt_max3 = 20 ;

	//** Preprocessing **//

	StopWatch preprocessing2(1) ;
	preprocessing2.startTime() ;
	constructGlobalConvectiveMatrix() ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;

	preprocessing2.startTime() ;
	constructGlobalGradientOperatorMatrixSuper() ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;

	preprocessing2.startTime() ;
	constructGlobalVelocityMassMatrixSuper() ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;

	preprocessing2.startTime() ;
	constructGlobalLaplaceVelocityMatrixPaper() ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;
	
	preprocessing2.startTime() ;
	constructTheRestOfTheMatricesVectors() ;
	constructNonlinearPaper() ; 
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;

	preprocessing2.startTime() ;
	// Set Multi-grid Tools //
	StopWatch MG(1) ;
	MG.startTime() ;

	if (mesh_type == true)  {constructMapFineCoarseNonUniformGrid() ;}
	
	MG.stopTime() ;
    MG.resetTime() ;
	multigrid = MG.finalTime() ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;

	preprocessing2.startTime() ;
	calculateGlobalPressureLaplaceOperatorMatrixInCoarseGridPaper(cycle) ;
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;

	preprocessing2.startTime() ;
	calculateGlobalDivergenceOperatorMatrixInCoarseGridPaper(cycle) ;		
	preprocessing2.stopTime() ;
    preprocessing2.resetTime() ;
	preprocess2 = preprocessing2.finalTime() ;
	
	// Processing //
	// Start the iterations

	StopWatch Adv(1) ;
	StopWatch Poisson(1) ;
	StopWatch Inj(1) ;
	StopWatch Pro(1) ;

	current_step = 1 ;
	concise = 1 ;

	while( current_step < time_step ){

		cout << "Current Step: " << current_step << endl ; cout << endl ;

		goToConvectionAdvanced() ;	
		updateRHS() ;
		applyMixedBCtoRHSofMomentum() ;
		breakRHS() ;

		//Eq1

		Adv.startTime() ;
		solveEq1primary(itr_max1/2 , rpt_max1/2) ;
		solveEq1secondary(itr_max1/2 , rpt_max1/2) ;
		Adv.stopTime() ;
		Adv.resetTime() ;

		// Multigrid //

		// Injection //

		Inj.startTime() ;

		injectIntermediateVelocityVktoVkMinus3() ;
		
		Inj.stopTime() ;
		Inj.resetTime() ;

		//Eq2

		calculateRHSPoissonInCoarseGridVkMinus3() ;
		applyMixedBCtoRHSofPoissonVkMinus3() ;

		Poisson.startTime() ;

		solveEq2VkMinus3(itr_max2_Vk_3 , itr_max2_Vk_3) ;
		
		Poisson.stopTime() ;
		Poisson.resetTime() ;

		// Prolongation //

		Pro.startTime() ;

		prolongationPressureVkMinus3toVk() ; 
		
		Pro.stopTime() ;
		Pro.resetTime() ;

		//Eq3
		
		calculateRHSRealVelocity() ;
		applyMixedBCToRHSofCorrectionEqu() ;
		
		solveEq3(itr_max3 , rpt_max3) ;

		transferGMREStoVectors() ;

		//Eq4
		updatePressure() ;
		
		// computeDragAndLift() ;

		exchangeVariables() ;
		current_step ++ ;
		updatePressureExp() ;
		concise = 1 ; // current_step ;
	}

	advection = Adv.finalTime() ; 
	poissonic = Poisson.finalTime() ;
	inject = Inj.finalTime() ; 
	prolong = Pro.finalTime() ;
	multigrid = MG.finalTime() ;

	cout << " The job is finished. " << endl ;

	return ;

	// current_step ; // physcially n in formulas, computationally n - 1
}

// -------------------- Exchange Variables --------------------

void Computation :: exchangeVariables(){
	
	for (int i = 0 ; i < velocity_node_number ; i ++){
		
		velocity[0][i] = velocity_u[0][i] = velocity_u[1][i] ;
		velocity[0][i + velocity_node_number] = velocity_v[0][i] = velocity_v[1][i] ;
	}

	for (int i = 0 ; i < pressure_node_number ; i ++)
		pressure[0][i] = pressure[1][i] ;
}

// -------------------- Product A*b --------------------

void Computation :: productAb(double A_aa[] , int A_ai[] , int A_aj[] , double b[] , double result[] , double coefficient , int N){

	double sum = 0 ;

	for (int i = 0 ; i < N ; i ++){
		for (int k = A_aj[i] ; k < A_aj[i+1] ; k ++ ){
       
		 sum +=  A_aa[k]*b[A_ai[k]] ;
      }
	result[i] = coefficient*sum ;
	sum = 0 ;
   }
}

// -------------------- Construct Time Dependent BC--------------------

void Computation :: constructTimeDependentBC(){
	
	Mesh Initializr("MeshP1.txt" , 1 , N_center) ;
	Initializr.generateMesh() ;
	
	intermediate_u_BC_N = Dirichlet_bound_number_velocity ;
	intermediate_v_BC_N = Dirichlet_bound_number_velocity ;
	
	Initializr.Destroy() ;
	//Initializr.~Mesh() ;

	intermediate_u_BC = new double *[time_step] ;
	for(int j = 0 ; j < time_step ; j ++)
		intermediate_u_BC[j] = new double[intermediate_u_BC_N] ;

	intermediate_v_BC = new double *[time_step] ;
	for(int j = 0 ; j < time_step ; j ++)
		intermediate_v_BC[j] = new double[intermediate_v_BC_N] ;
}

// -------------------- Update CSR --------------------

void Computation :: updateCSR(int r , int s , double variable , double A_aa[] , int A_ai[] , int A_aj[] , int A_non_zeros_elements){

	// A(r , s) = variable
	// Note: Indices start from zero!

	for (int k = 0 ; k < A_non_zeros_elements ; k ++){
		
		if((A_ai[k] == s) && (A_aj[r] < k + 1 )){
			
			//A_aa[k] = variable ; // Normal Update
			A_aa[k] += variable ; //Finite Element Format, Using the format, A_aa[k] should be zero!
			break ;
		}
	}
}

// -------------------- Search CSR--------------------

int Computation :: searchCSR(int r , int s , int A_ai[] , int A_aj[] , int A_non_zeros_elements){

	// A(r , s) = variable
	// Note: Indices start from zero!

	for (int k = 0 ; k < A_non_zeros_elements ; k ++){
		
		if((A_ai[k] == s) && (A_aj[r] < k + 1 )){
			
			return k ;
			//A_aa[k] = variable ; // Normal Update
			//A_aa[k] += variable ; //Finite Element Format, Using the format, A_aa[k] should be zero!
		}
	}
}

// -------------------- zero CSR --------------------

void Computation :: zeroCSR(double A_aa[] , int A_non_zeros_elements){
	
	for (int k = 0 ; k < A_non_zeros_elements ; k ++)
		A_aa[k] = 0.0 ; 
}

// -------------------- Convert AAIJ to CSR --------------------

void Computation :: convertAAIJtoCSR(int size_i , int non_zeros_elements , AAIJ A[] , double AA_aa[] , int AA_ai[] , int AA_aj[]){

	// Note: Indices start from zero!

	sort (A , A + non_zeros_elements , compareIindex) ;

	for (int i = 0 ; i < non_zeros_elements ; i ++){
		
		AA_aa[i] = A[i].aa ;
		AA_ai[i] = A[i].j ;
	}

	AA_aj[0] = 0 ;
	int counter = 0 ;
	int tally = 1 ;

	while(counter < non_zeros_elements ){
	
		if( A[counter].i < A[counter + 1].i){

			AA_aj[tally] = counter + 1 ;
			tally ++ ;
		}
		
		counter ++ ;
	}

	AA_aj[size_i] = non_zeros_elements ;
}

// -------------------- Update AAIJ --------------------

int Computation :: updateAAIJ( AAIJ A[] , double value , int I , int J , int spicy , int non_zero_approximation){

	int checker = checkSparsity(A , I , J , spicy , non_zero_approximation) ;
	
	if (checker != -1){

	A[checker].aa += value ;
	A[checker].i = I ;
	A[checker].j = J ;
	return spicy ; }
				
	if (checker == -1){

	A[spicy].aa += value ;
	A[spicy].i = I ;
	A[spicy].j = J ;
	spicy ++ ;
	return spicy ;}
}

// -------------------- Impose AAIJ --------------------

int Computation :: imposeAAIJ( AAIJ A[] , double value , int I , int J , int spicy , int non_zero_approximation){

	int checker = checkSparsity(A , I , J , spicy , non_zero_approximation) ;
	
	if (checker != -1 && value == 1.0){

	A[checker].aa = value ;
	A[checker].i = I ;
	A[checker].j = J ;
	return spicy ; }
				
	if (checker == -1 && value == 1.0){

	A[spicy].aa += value ;
	A[spicy].i = I ;
	A[spicy].j = J ;
	spicy ++ ;
	return spicy ;}

	if (checker != -1 && value == 0.0){

	A[checker].aa = value ;
	A[checker].i = I ;
	A[checker].j = J ;
	return spicy ; }

	if (checker == -1 && value == 0.0){
	return spicy ;}
}

// -------------------- Compute Residual --------------------

double Computation :: residual(double v_previous[] , double v_current[] , int N_size){
	
	double second_norm = 0 ;
	double u_previous_norm = 0 ;
	double u_current_norm = 0 ;

	for (int i = 0 ; i < N_size ; i ++){
	
		u_previous_norm = sqrt(pow(v_previous[i] , 2.0) + pow(v_previous[i+N_size] , 2.0)) ;
		u_current_norm = sqrt(pow(v_current[i] , 2.0) + pow(v_current[i+N_size] , 2.0)) ;
		second_norm += pow(u_current_norm - u_previous_norm , 2.0 ) ;	
	}

	return sqrt (second_norm) ; 
}

// -------------------- Compute Stress --------------------

void Computation :: computeStress(){
	 
	double *R1 , *R2 , *result ;
	R1 = new double [2*velocity_node_number] ;
	R2 = new double [2*velocity_node_number] ;
	result = new double [2*velocity_node_number] ;
	
	for (int i = 0 ; i < 2*velocity_node_number ; i ++){
		
		stress[i] = 1.0 ;
		R1[i] = 0.0 ;
		R2[i] = 0.0 ;
		result[i] = 0.0 ;
	}

	productAb(Gradient_aa , Gradient_ai , Gradient_aj , velocity_u[1] , R1 , 1.0 , 2*velocity_node_number) ;
	productAb(Gradient_aa , Gradient_ai , Gradient_aj , velocity_v[1] , R2 , 1.0 , 2*velocity_node_number) ;
	
	for (int i = 0 ; i < velocity_node_number ; i ++){
		result[i] = R2[i] + R1[i + velocity_node_number] ;
		result[i + velocity_node_number] = result[i] ;
	}

	delete [] R1 ;
	delete [] R2 ;

    pmgmres_ilu_cr ( 2*velocity_node_number , Mass_Velocity_non_zeros_elements , Mass_Velocity_aj , Mass_Velocity_ai , Mass_Velocity_aa , stress , result , 2*velocity_node_number - 10 , 20 , pow(10.0 , -8.0) , pow(10.0 , -8.0));
		
	delete [] result ;
}

// -------------------- Compute Drag and Lift --------------------

void Computation :: computeDragAndLift(){

	computeStress() ;
	
	double sum = 0.0 ;
	for (int i = 0 ; i < Node_Cyl ; i ++){
		
		sum += (pressure[1][Cyl_Inf[i].node - 1]*Cyl_Inf[i].Co*(Cyl_Inf[i].d_teta) + viscosity*stress[Cyl_Inf[i].node - 1]*Cyl_Inf[i].Si*(Cyl_Inf[i].d_teta)) ;	
		}

	Drag[current_step] =  sum ;

	sum = 0.0 ;
	for (int i = 0 ; i < Node_Cyl ; i ++){
		
		sum += (-pressure[1][Cyl_Inf[i].node - 1]*Cyl_Inf[i].Si*(Cyl_Inf[i].d_teta) + viscosity*stress[Cyl_Inf[i].node - 1]*Cyl_Inf[i].Co*(Cyl_Inf[i].d_teta)) ;	
		}

	Lift[current_step] = sum ;
	
	sum = 0.0 ;
	for (int i = 0 ; i < Node_Cyl ; i ++){
		
		sum += (pressure[1][Cyl_Inf[i].node - 1]*Cyl_Inf[i].Co*(Cyl_Inf[i].d_teta)) ;	
		}

	PressureDrag[current_step] =  sum ;
	
	sum = 0.0 ;
	for (int i = 0 ; i < Node_Cyl ; i ++){
		
		sum += (viscosity*stress[Cyl_Inf[i].node - 1]*Cyl_Inf[i].Si*(Cyl_Inf[i].d_teta)) ;	
		}

	ViscousDrag[current_step] =  sum ;
	
	sum = 0.0 ;
	for (int i = 0 ; i < Node_Cyl ; i ++){
		
		sum += (-pressure[1][Cyl_Inf[i].node - 1]*Cyl_Inf[i].Si*(Cyl_Inf[i].d_teta)) ;	
		}

	PressureLift[current_step] = sum ;
	
	sum = 0.0 ;
	for (int i = 0 ; i < Node_Cyl ; i ++){
		
		sum += (viscosity*stress[Cyl_Inf[i].node - 1]*Cyl_Inf[i].Co*(Cyl_Inf[i].d_teta)) ;	
		}

	ViscousLift[current_step] = sum ;
}

// -------------------- Print Velocity Contour TecPlot --------------------

void Computation :: printVelocityContourTecPlot(double desired_time){
	
	ifstream fin ("MeshP1.txt") ;
	
	int Node_N , Element_N ;
	
	fin >> Node_N ;
	fin >> Element_N ;

	int step = concise ; // int(desired_time/delta_t) ;
	string x_c ;
	ostringstream convert ;
	convert << desired_time ;
	x_c = convert.str() ;
	string s1 , s2 , s3 , s4 , s5 , s6 ;

	s1 = "velocityU" + x_c + ".txt" ;
	s2 = "velocityV" + x_c + ".txt" ;
	s3 = "velocity" + x_c + ".txt" ;
	s4 = "IntermeidateVelocityMagnitude" + x_c + ".txt" ;
	s5 = "IntermeidateVelocityU" + x_c + ".txt" ;
	s6 = "IntermeidateVelocityV" + x_c + ".txt" ;
	
	ofstream fout1 (s1.c_str()) ;
	ofstream fout2 (s2.c_str()) ;
	ofstream fout3 (s3.c_str()) ;
	ofstream fout4 (s4.c_str()) ;
	ofstream fout5 (s5.c_str()) ;
	ofstream fout6 (s6.c_str()) ;
	
	fout1 << "title=Velocity U" << endl ;
	fout1 << "VARIABLES=X, Y, P " << endl ;
	fout1 << "ZONE N = " << Node_N << ",E = " << Element_N << " ,DATAPACKING=POINT, ZONETYPE=FETRIANGLE " << endl ; 
	
	fout2 << "title=Velocity V" << endl ;
	fout2 << "VARIABLES=X, Y, P " << endl ;
	fout2 << "ZONE N = " << Node_N << ",E = " << Element_N << " ,DATAPACKING=POINT, ZONETYPE=FETRIANGLE " << endl ; 
	
	fout3 << "title=Velocity" << endl ;
	fout3 << "VARIABLES=X, Y, P " << endl ;
	fout3 << "ZONE N = " << Node_N << ",E = " << Element_N << " ,DATAPACKING=POINT, ZONETYPE=FETRIANGLE " << endl ; 
	
	fout4 << "title=Intermediate Velocity" << endl ;
	fout4 << "VARIABLES=X, Y, P " << endl ;
	fout4 << "ZONE N = " << Node_N << ",E = " << Element_N << " ,DATAPACKING=POINT, ZONETYPE=FETRIANGLE " << endl ; 
	
	fout5 << "title=Intermediate Velocity U" << endl ;
	fout5 << "VARIABLES=X, Y, P " << endl ;
	fout5 << "ZONE N = " << Node_N << ",E = " << Element_N << " ,DATAPACKING=POINT, ZONETYPE=FETRIANGLE " << endl ; 
	
	fout6 << "title=Intermediate Velocity V" << endl ;
	fout6 << "VARIABLES=X, Y, P " << endl ;
	fout6 << "ZONE N = " << Node_N << ",E = " << Element_N << " ,DATAPACKING=POINT, ZONETYPE=FETRIANGLE " << endl ; 
	
	double x_coordinate , y_coordinate , z_coordinate , extra ;
	
	int index1 , index2 , index3 ; 

	for (int i = 0 ; i < Node_N  ; i ++){
	
			fin >> x_coordinate ;
			fin >> y_coordinate ;
			fin >> z_coordinate ;
			
			fout1 << x_coordinate << "  " << y_coordinate << "  " << velocity_u[step][i] << endl ;		
			fout2 << x_coordinate << "  " << y_coordinate << "  " << velocity_v[step][i] << endl ;
			fout3 << x_coordinate << "  " << y_coordinate << "  " << sqrt(pow(velocity_u[step][i] , 2.0) + pow(velocity_v[step][i] , 2.0)) << endl ;
			fout4 << x_coordinate << "  " << y_coordinate << "  " << sqrt (pow(intermediate_velocity[i] , 2.0) + pow(intermediate_velocity[i + velocity_node_number] , 2.0)) << endl ;		
			fout5 << x_coordinate << "  " << y_coordinate << "  " << intermediate_velocity[i] << endl ;		
			fout6 << x_coordinate << "  " << y_coordinate << "  " << intermediate_velocity[i + velocity_node_number] << endl ;		
		}
	
	for (int i = 0 ; i < Element_N  ; i ++){

			fin >> extra ;
			fin >> index1 ;
			fin >> index2 ;
			fin >> index3 ;
			
			fout1 << index1 + 1 - N_center  << "  " << index2 + 1 - N_center  << "  " << index3 + 1 - N_center << endl ;
			fout2 << index1 + 1 - N_center  << "  " << index2 + 1 - N_center  << "  " << index3 + 1 - N_center << endl ;
			fout3 << index1 + 1 - N_center  << "  " << index2 + 1 - N_center  << "  " << index3 + 1 - N_center << endl ;
			fout4 << index1 + 1 - N_center  << "  " << index2 + 1 - N_center  << "  " << index3 + 1 - N_center << endl ;
			fout5 << index1 + 1 - N_center  << "  " << index2 + 1 - N_center  << "  " << index3 + 1 - N_center << endl ;
			fout6 << index1 + 1 - N_center  << "  " << index2 + 1 - N_center  << "  " << index3 + 1 - N_center << endl ;
		}

	fin.close() ;
	return ;
}

// -------------------- Print Pressure Contour TecPlot --------------------

void Computation :: printPressureContourTecPlot(double desired_time){
	
	ifstream fin ("MeshP1.txt") ;
	
	int Node_N , Element_N ;
	
	fin >> Node_N ;
	fin >> Element_N ;

	int step = concise ; // int(desired_time/delta_t) ;
	
	string x_c ;
	ostringstream convert ;
	convert << desired_time ;
	x_c = convert.str() ;
	string s1 , s2 , s3 ;

	s1 = "PressureT" + x_c + ".txt" ;
	s2 = "MediatePressureT" + x_c + ".txt" ;
	s3 = "GradientPressureT" + x_c + ".txt" ;

	ofstream fout1 (s1.c_str()) ;
	ofstream fout2 (s2.c_str()) ;
	ofstream fout3 (s3.c_str()) ;

	fout1 << "title=Pressure Field" << endl ;
	fout1 << "VARIABLES=X, Y, P " << endl ;
	fout1 << "ZONE N = " << Node_N << ",E = " << Element_N << " ,DATAPACKING=POINT, ZONETYPE=FETRIANGLE " << endl ; 
	
	fout2 << "title=Intermediate Pressure Field" << endl ;
	fout2 << "VARIABLES=X, Y, P " << endl ;
	fout2 << "ZONE N = " << Node_N << ",E = " << Element_N << " ,DATAPACKING=POINT, ZONETYPE=FETRIANGLE " << endl ; 
	
	double x_coordinate , y_coordinate , z_coordinate , extra ;
	
	int index1 , index2 , index3 ; 

	for (int i = 0 ; i < Node_N  ; i ++){
	
			fin >> x_coordinate ;
			fin >> y_coordinate ;
			fin >> z_coordinate ;
			
			fout1 << x_coordinate << "  " << y_coordinate << "  " << pressure[step][i] << endl ;		
			fout2 << x_coordinate << "  " << y_coordinate << "  " << fi_fine_grid[i] << endl ;		
		}
	
	for (int i = 0 ; i < Element_N  ; i ++){

			fin >> extra ;
			fin >> index1 ;
			fin >> index2 ;
			fin >> index3 ;
			
			fout1 << index1 + 1 - N_center  << "  " << index2 + 1 - N_center  << "  " << index3 + 1 - N_center << endl ;
			fout2 << index1 + 1 - N_center  << "  " << index2 + 1 - N_center  << "  " << index3 + 1 - N_center << endl ;
		}
	
	double *grad_p ;
	grad_p = new double [2*velocity_node_number] ;
	productAb(Gradient_aa , Gradient_ai , Gradient_aj , pressure[1] , grad_p , 1.0 , 2*velocity_node_number ) ;
	
	for (int i = 0 ; i < 2*Node_N ; i ++){
		fout3 << i << "  " <<  grad_p[i] << endl ; }

	fin.close() ;
	return ;
}

// -------------------- Print Stream Lines TecPlot --------------------

void Computation :: printStreamLinesTecPlot(double desired_time){
	
	ifstream fin ("MeshP1.txt") ;
	
	int Node_N , Element_N ;
	
	fin >> Node_N ;
	fin >> Element_N ;
	
	int step = concise ; // int(desired_time/delta_t) ;
	
	string x_c ;
	ostringstream convert ;
	convert << desired_time ;
	x_c = convert.str() ;
	string s1 ;

	s1 = "streamlinesT" + x_c + ".txt" ;
	ofstream fout1 (s1.c_str()) ;

	fout1 << "VARIABLES = X , Y, U, V" << endl ;
	fout1 << "ZONE N = " << Node_N << ",E = " << Element_N << " ,DATAPACKING=POINT, ZONETYPE=FETRIANGLE " << endl ; 
		
	double x_coordinate , y_coordinate , z_coordinate , extra ;
	
	int index1 , index2 , index3 ; 

	for (int i = 0 ; i < Node_N  ; i ++){
	
			fin >> x_coordinate ;
			fin >> y_coordinate ;
			fin >> z_coordinate ;
			
			fout1 << x_coordinate << "  " << y_coordinate << "  " << velocity_u[step][i] << "  " << velocity_v[step][i] << endl ;		
			
		}
	
	for (int i = 0 ; i < Element_N  ; i ++){

			fin >> extra ;
			fin >> index1 ;
			fin >> index2 ;
			fin >> index3 ;
			
			fout1 << index1 + 1 - N_center  << "  " << index2 + 1 - N_center  << "  " << index3 + 1 - N_center << endl ;
		}

	fin.close() ;
}

// -------------------- Print Vorticity Tecplot --------------------

void Computation :: printVorticityTecPlot(double desired_time){
	
	double *vorticity , *R1 , *R2 , *result ;
	vorticity = new double [2*velocity_node_number] ;
	R1 = new double [2*velocity_node_number] ;
	R2 = new double [2*velocity_node_number] ;
	result = new double [2*velocity_node_number] ;
	
	for (int i = 0 ; i < 2*velocity_node_number ; i ++){
		vorticity[i] = 1.0 ;
		R1[i] = 0.0 ;
		R2[i] = 0.0 ;
		result[i] = 0.0 ;
	}

	productAb(Gradient_aa , Gradient_ai , Gradient_aj , velocity_u[1] , R1 , 1.0 , 2*velocity_node_number) ;
	productAb(Gradient_aa , Gradient_ai , Gradient_aj , velocity_v[1] , R2 , 1.0 , 2*velocity_node_number) ;
	
	for (int i = 0 ; i < velocity_node_number ; i ++){
		result[i] = R2[i] - R1[i + velocity_node_number] ;
		result[i + velocity_node_number] = result[i] ;
	}

	delete [] R1 ;
	delete [] R2 ;

    pmgmres_ilu_cr ( 2*velocity_node_number , Mass_Velocity_non_zeros_elements , Mass_Velocity_aj , Mass_Velocity_ai , Mass_Velocity_aa , vorticity , result , 2*velocity_node_number - 10 , 20 , pow(10.0 , -8.0) , pow(10.0 , -8.0));
		
	delete [] result ;

	ifstream fin ("MeshP1.txt") ;
	
	int Node_N , Element_N ;
	
	fin >> Node_N ;
	fin >> Element_N ;
	
	int step = concise ; // int(desired_time/delta_t) ;
	
	string x_c ;
	ostringstream convert ;
	convert << desired_time ;
	x_c = convert.str() ;
	string s1 ;

	s1 = "vorticityT" + x_c + ".txt" ;
	ofstream fout1 (s1.c_str()) ;

	fout1 << "VARIABLES = X , Y, W" << endl ;
	fout1 << "ZONE N = " << Node_N << ",E = " << Element_N << " ,DATAPACKING=POINT, ZONETYPE=FETRIANGLE " << endl ; 
		
	double x_coordinate , y_coordinate , z_coordinate , extra ;
	
	int index1 , index2 , index3 ; 

	for (int i = 0 ; i < Node_N  ; i ++){
	
			fin >> x_coordinate ;
			fin >> y_coordinate ;
			fin >> z_coordinate ;
			
			fout1 << x_coordinate << "  " << y_coordinate << "  " << vorticity[i] << endl ;		
			
		}
	
	for (int i = 0 ; i < Element_N  ; i ++){

			fin >> extra ;
			fin >> index1 ;
			fin >> index2 ;
			fin >> index3 ;
			
			fout1 << index1 + 1 - N_center  << "  " << index2 + 1 - N_center  << "  " << index3 + 1 - N_center << endl ;
		}

	fin.close() ;
	
	delete [] vorticity ;
}

// -------------------- Print Drag And Lift --------------------

void Computation :: printDragAndLift(){
	
	ofstream fout1("Drag.txt") ;
	ofstream fout2("Lift.txt") ;
	ofstream fout3("PressureDrag.txt") ;
	ofstream fout4("ViscousDrag.txt") ;
	ofstream fout5("PressureLift.txt") ;
	ofstream fout6("ViscousLift.txt") ;

	for (int i = 0 ; i < time_step ; i ++){
		
		fout1 << ((double)i)*delta_t << "  " << Drag[i] << endl ;
		fout2 << ((double)i)*delta_t << "  " << Lift[i] << endl ;
		fout3 << ((double)i)*delta_t << "  " << PressureDrag[i] << endl ;
		fout4 << ((double)i)*delta_t << "  " << ViscousDrag[i] << endl ;
		fout5 << ((double)i)*delta_t << "  " << PressureLift[i] << endl ;
		fout6 << ((double)i)*delta_t << "  " << ViscousLift[i] << endl ;
	}
}

// -------------------- Print Variable Function --------------------

void Computation :: printVariableFunction(double x , double y , double desired_time){
	
	ifstream fin ("MeshP1.txt") ;
	
	int Node_N , Element_N ;
	
	fin >> Node_N ;
	fin >> Element_N ;
	
	interspace *Function ;
	Function = new interspace[Node_N] ;

	double x_coordinate , y_coordinate , z_coordinate , extra ;
	int index1 , index2 , index3 ; 

	for (int i = 0 ; i < Node_N  ; i ++){
			
			fin >> x_coordinate ;
			fin >> y_coordinate ;
			fin >> z_coordinate ;

			double D = sqrt(pow( x - x_coordinate , 2.0) + pow( y - y_coordinate , 2.0)) ;
			Function[i].D = D ;
			Function[i].index = i ;
		}

	sort(Function , Function + Node_N , compareDistance) ;

	fin.close() ;

	int step = int(desired_time/delta_t) ;
	
	ofstream fout1 ("UTime.txt") ;
	ofstream fout2 ("VTime.txt") ;

	fout1 << "VARIABLES = t , U" << endl ;
	fout2 << "VARIABLES = t , V" << endl ;		

	for (int k = 0 ; k < step + 1 ; k ++){
		
		double time = k*delta_t ;
		
		fout1 << time << "  " << velocity_u[k][Function[0].index] << endl ;
		fout2 << time << "  " << velocity_v[k][Function[0].index] << endl ;	
	}
}

// -------------------- Print Time Processing --------------------

void Computation :: printTimeProcessing(float T){

	preprocess1 = T ;

	ofstream fout ("Time.txt") ;
	
	if (cycle == 0){
		
	float Full = advection + poissonic + preprocess2 + preprocess1 ;

	fout << "Level: " << cycle << endl ;
	fout << endl ;
	fout << "Full Time: " << Full << endl ;
	fout << endl ;
	fout << "Preprocess: " <<  preprocess2 + preprocess1 << endl ;
	fout << endl ;
	fout << "Process: " << advection + poissonic << endl ;
	fout << endl ;
	fout << "Poisson: " << poissonic << endl ;
	fout << endl ;
	fout << "Advection & Velocity Correction: " << advection << endl ;
	fout << endl ;
	fout << "Mapp: " << 0 << endl ;
	fout << endl ;
	fout << "Injection: " << 0  << endl ;
	fout << endl ;
	fout << "Prolongation: " << 0  << endl ;
	fout << endl ;
	fout << "Multigrid Tool: " << 0  << endl ;
	
	fout << endl ;
	fout << " ************** " << endl ;
	fout << endl ;

	fout << "% Preprocess: " <<  ((preprocess2 + preprocess1)/Full)*100.00 << endl ;
	fout << endl ;
	fout << "% Process: " << ((advection + poissonic)/Full)*100.00   << endl ;
	fout << endl ;
	fout << "% Poisson: " << (poissonic/Full)*100.00 << endl ;
	fout << endl ;
	fout << "% Advection & Velocity Correction: " << (advection/Full)*100.00 << endl ;
	fout << endl ;
	fout << "% Mapp: " << 0 << endl ;
	fout << endl ;
	fout << "% Injection: " << 0  << endl ;
	fout << endl ;
	fout << "% Prolongation: " << 0 << endl ;
	
	}

	if (cycle != 0){

	float Full = (advection + poissonic + inject + prolong + preprocess1 + preprocess2) ;
	
	fout << "Level: " << cycle << endl ;
	fout << endl ;
	fout << "Full Time: " << Full << endl ;
	fout << endl ;
	fout << "Preprocess: " <<  preprocess2 + preprocess1 << endl ;
	fout << endl ;
	fout << "Process: " << advection + poissonic + inject + prolong   << endl ;
	fout << endl ;
	fout << "Poisson: " << poissonic << endl ;
	fout << endl ;
	fout << "Advection & Velocity Correction: " << advection << endl ;
	fout << endl ;
	fout << "Mapp: " << inject + prolong << endl ;
	fout << endl ;
	fout << "Injection: " << inject  << endl ;
	fout << endl ;
	fout << "Prolongation: " << prolong  << endl ;
	fout << endl ;
	fout << "Multigrid: " << multigrid  << endl ;
	
	fout << endl ;
	fout << " ************** " << endl ;
	fout << endl ;

	fout << "% Preprocess: " <<  ((preprocess2 + preprocess1)/Full)*100.00 << endl ;
	fout << endl ;
	fout << "% Process: " << ((advection + poissonic + inject + prolong)/Full)*100.00   << endl ;
	fout << endl ;
	fout << "% Poisson: " << (poissonic/Full)*100.00 << endl ;
	fout << endl ;
	fout << "% Advection & Velocity Correction: " << (advection/Full)*100.00 << endl ;
	fout << endl ;
	fout << "% Mapp: " << ((inject + prolong)/(Full))*100.00 << endl ;
	fout << endl ;
	fout << "% Injection: " << (inject/Full)*100.00  << endl ;
	fout << endl ;
	fout << "% Prolongation: " << (prolong/Full)*100.00  << endl ;
	
	}
}

// ===============================
// === Coarse Grid Projection  ===
// ===============================

// -------------------- Construct Map Fine Coarse Non Uniform Grid --------------------

void Computation :: constructMapFineCoarseNonUniformGrid(){
		
	if (cycle == 1){constructMapFineCoarseNonUniformGridVkMinus1() ; return ;}
	if (cycle == 2){constructMapFineCoarseNonUniformGridVkMinus2() ; return ;}
	if (cycle == 3){constructMapFineCoarseNonUniformGridVkMinus3() ; return ;}
	
	return ;
}

// -------------------- Construct Map Fine Coarse Non Uniform Grid V(k - 1) --------------------

void Computation :: constructMapFineCoarseNonUniformGridVkMinus1(){

	//  Cycle == 1

	// Restriction
	 
	Coordinate *artificial_c1_mesh ;
	artificial_c1_mesh = new Coordinate[pressure_node_number_minus_1] ; 

	double z_stuff = 0.0 ;

	ifstream finC1("MeshP1C1.txt") ;
	finC1 >> z_stuff ;
	finC1 >> z_stuff ;
	
	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
		
		artificial_c1_mesh[i].node_number = i + 1 ;
		finC1 >> artificial_c1_mesh[i].x_coordinate ;  
		finC1 >> artificial_c1_mesh[i].y_coordinate ;
		finC1 >> z_stuff ;			
	}
	
	finC1.close() ;
	
	Coordinate *artificial_f_mesh ;
	artificial_f_mesh = new Coordinate[pressure_node_number] ; 

	ifstream finF("MeshP1.txt") ;
	finF >> z_stuff ;
	finF >> z_stuff ;
	
	for (int i = 0 ; i < pressure_node_number ; i ++){
		
		artificial_f_mesh[i].node_number = i + 1 ;
		finF >> artificial_f_mesh[i].x_coordinate ;  
		finF >> artificial_f_mesh[i].y_coordinate ;
		finF >> z_stuff ;			
	}
	
	finF.close() ;
	
	mapp_V_k_minus_1 = new restriction[pressure_node_number_minus_1] ;
	int counter = 0 ; 

	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
		for (int j = 0 ; j < pressure_node_number ; j ++){
			
			if ((artificial_f_mesh[j].x_coordinate == artificial_c1_mesh[i].x_coordinate) && (artificial_f_mesh[j].y_coordinate == artificial_c1_mesh[i].y_coordinate)){
				
				mapp_V_k_minus_1[counter].N_coarse = artificial_c1_mesh[i].node_number ;
				mapp_V_k_minus_1[counter].N_fine = artificial_f_mesh[j].node_number ;
				counter ++ ;
				continue ;
			}
		}	
	}
	
	// Prolongation
	
	index_V_k_minus_1 = new prolongation[pressure_node_number] ;
	
	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){

		index_V_k_minus_1[i].i = mapp_V_k_minus_1[i].N_coarse ; 
		index_V_k_minus_1[i].j = mapp_V_k_minus_1[i].N_coarse ; 
		index_V_k_minus_1[i].k = mapp_V_k_minus_1[i].N_fine ; 
	}

	Mesh Coarse("MeshP1C1.txt" , 1 , N_center) ;
	Coarse.generateMesh() ;
	
	middle *imaginary ;
	int imaginary_length = 3 * Coarse.number_of_elements ; 
	imaginary = new middle[imaginary_length] ;
	
	counter = 0 ; 

	for (int i = 0 ; i < Coarse.number_of_elements ; i ++){
	
		imaginary[counter].x = (Coarse.domain[i].x1_global + Coarse.domain[i].x2_global)/2.0 ; 
		imaginary[counter].y = (Coarse.domain[i].y1_global + Coarse.domain[i].y2_global)/2.0 ;
		imaginary[counter].N1 = Coarse.domain[i].node_1 ;
		imaginary[counter].N2 = Coarse.domain[i].node_2 ;

		counter ++ ;

		imaginary[counter].x = (Coarse.domain[i].x1_global + Coarse.domain[i].x3_global)/2.0 ; 
		imaginary[counter].y = (Coarse.domain[i].y1_global + Coarse.domain[i].y3_global)/2.0 ;
		imaginary[counter].N1 = Coarse.domain[i].node_1 ;
		imaginary[counter].N2 = Coarse.domain[i].node_3 ;

		counter ++ ;

		imaginary[counter].x = (Coarse.domain[i].x3_global + Coarse.domain[i].x2_global)/2.0 ; 
		imaginary[counter].y = (Coarse.domain[i].y3_global + Coarse.domain[i].y2_global)/2.0 ;
		imaginary[counter].N1 = Coarse.domain[i].node_3 ;
		imaginary[counter].N2 = Coarse.domain[i].node_2 ;

		counter ++ ;
	}

	Coarse.Destroy() ;
	//Coarse.~Mesh() ;
	
	sort (imaginary , imaginary + imaginary_length , compareMiddle) ;
	
	double x_test , y_test ;
	x_test = imaginary[0].x ;
	y_test = imaginary[0].y ;
	
	middle *total_imaginary ; 
	total_imaginary = new middle[pressure_node_number - pressure_node_number_minus_1] ;
	
	total_imaginary[0].x = imaginary[0].x ;
	total_imaginary[0].y = imaginary[0].y ;
	total_imaginary[0].N1 = imaginary[0].N1 ;
	total_imaginary[0].N2 = imaginary[0].N2 ;

	int tally = 1 ; 

	for (int i = 1 ; i < imaginary_length ; i ++){
		
		if ((x_test == imaginary[i].x) && (y_test == imaginary[i].y)){continue ;}

		total_imaginary[tally].x = imaginary[i].x ;
		total_imaginary[tally].y = imaginary[i].y ;
		total_imaginary[tally].N1 = imaginary[i].N1 ;
		total_imaginary[tally].N2 = imaginary[i].N2 ;
		tally ++ ;
		x_test = imaginary[i].x ; 
		y_test = imaginary[i].y ;
	}

	interspace *Distance ;
	Distance = new interspace[pressure_node_number - pressure_node_number_minus_1] ;
	
	counter = pressure_node_number_minus_1 ;
	int repetitious = 0 ;
	
	for (int j = 0 ; j < pressure_node_number ; j ++){
		
		double x = artificial_f_mesh[j].x_coordinate ; 
		double y = artificial_f_mesh[j].y_coordinate ;
		
		for (int k = 0 ; k < pressure_node_number_minus_1 ; k ++){
			if (x == artificial_c1_mesh[k].x_coordinate && y == artificial_c1_mesh[k].y_coordinate){repetitious = 1 ; break ;}
		}

		if (repetitious == 1){ repetitious = 0 ; continue ;}

		for(int i = 0 ; i < pressure_node_number - pressure_node_number_minus_1 ; i ++){
			
			Distance[i].D = sqrt (pow( x - total_imaginary[i].x , 2.0) + pow( y - total_imaginary[i].y , 2.0)) ; 
			Distance[i].index = i ;
		}

		sort(Distance , Distance + (pressure_node_number - pressure_node_number_minus_1) , compareDistance) ;
		
		index_V_k_minus_1[counter].i = total_imaginary[Distance[0].index].N1 ;
		index_V_k_minus_1[counter].j = total_imaginary[Distance[0].index].N2 ;
		index_V_k_minus_1[counter].k = artificial_f_mesh[j].node_number ;
		counter ++ ;
	}

	delete [] total_imaginary ;
	delete [] Distance ;
	delete [] imaginary ;
	
	delete [] artificial_f_mesh ;
	delete [] artificial_c1_mesh ;

	// For test now

	Pro = new double *[pressure_node_number] ;
	for (int i = 0 ; i < pressure_node_number ; i ++)
		Pro[i] = new double[pressure_node_number_minus_1] ;
	
	Res = new double *[pressure_node_number_minus_1] ;
	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++)
		Res[i] = new double[pressure_node_number] ;
	
	for (int i = 0 ; i < pressure_node_number ; i ++){
		for (int j = 0 ; j < pressure_node_number_minus_1 ; j ++){
			Pro[i][j] = 0.0 ;		
		}
	}

	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
		for (int j = 0 ; j < pressure_node_number ; j ++){
			Res[i][j] = 0.0 ;
		}
	}
	
	for (int i = 0 ; i < pressure_node_number ; i ++){

			//fi_fine_grid[index_V_k_minus_1[i].k - 1] = (fi_coarse_grid_L_1[index_V_k_minus_1[i].i - 1] + fi_coarse_grid_L_1[index_V_k_minus_1[i].j - 1])/2.0 ;
			if ( index_V_k_minus_1[i].i == index_V_k_minus_1[i].j ){ Pro[index_V_k_minus_1[i].k - 1][index_V_k_minus_1[i].i - 1] = 1.0 ; continue ;}
			Pro[index_V_k_minus_1[i].k - 1][index_V_k_minus_1[i].i - 1] = 0.5 ;
			Pro[index_V_k_minus_1[i].k - 1][index_V_k_minus_1[i].j - 1] = 0.5 ;
		}

	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
		for (int j = 0 ; j < pressure_node_number ; j ++){
			Res[i][j] = Pro[j][i] ;
		}
	}

	double sume = 0 ; 

	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
		for (int j = 0 ; j < pressure_node_number ; j ++){
			sume += Res[i][j] ;
		}
		for (int j = 0 ; j < pressure_node_number ; j ++){
		Res[i][j] = Res[i][j]/sume ;
		}
		sume = 0.0 ;
	}
}

// -------------------- Construct Map Fine Coarse Non Uniform Grid V(k - 2) --------------------

void Computation :: constructMapFineCoarseNonUniformGridVkMinus2(){

	//  Cycle == 1

	// Restriction
	 
	Coordinate *artificial_c1_mesh ;
	artificial_c1_mesh = new Coordinate[pressure_node_number_minus_1] ; 

	double z_stuff = 0.0 ;

	ifstream finC1("MeshP1C1.txt") ;
	finC1 >> z_stuff ;
	finC1 >> z_stuff ;
	
	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
		
		artificial_c1_mesh[i].node_number = i + 1 ;
		finC1 >> artificial_c1_mesh[i].x_coordinate ;  
		finC1 >> artificial_c1_mesh[i].y_coordinate ;
		finC1 >> z_stuff ;			
	}
	
	finC1.close() ;
	
	Coordinate *artificial_f_mesh ;
	artificial_f_mesh = new Coordinate[pressure_node_number] ; 

	ifstream finF("MeshP1.txt") ;
	finF >> z_stuff ;
	finF >> z_stuff ;
	
	for (int i = 0 ; i < pressure_node_number ; i ++){
		
		artificial_f_mesh[i].node_number = i + 1 ;
		finF >> artificial_f_mesh[i].x_coordinate ;  
		finF >> artificial_f_mesh[i].y_coordinate ;
		finF >> z_stuff ;			
	}
	
	finF.close() ;
	
	mapp_V_k_minus_1 = new restriction[pressure_node_number_minus_1] ;
	int counter = 0 ; 

	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
		for (int j = 0 ; j < pressure_node_number ; j ++){
			
			if ((artificial_f_mesh[j].x_coordinate == artificial_c1_mesh[i].x_coordinate) && (artificial_f_mesh[j].y_coordinate == artificial_c1_mesh[i].y_coordinate)){
				
				mapp_V_k_minus_1[counter].N_coarse = artificial_c1_mesh[i].node_number ;
				mapp_V_k_minus_1[counter].N_fine = artificial_f_mesh[j].node_number ;
				counter ++ ;
				continue ;
			}
		}	
	}
	
	// Prolongation
	
	index_V_k_minus_1 = new prolongation[pressure_node_number] ;
	
	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){

		index_V_k_minus_1[i].i = mapp_V_k_minus_1[i].N_coarse ; 
		index_V_k_minus_1[i].j = mapp_V_k_minus_1[i].N_coarse ; 
		index_V_k_minus_1[i].k = mapp_V_k_minus_1[i].N_fine ; 
	}

	Mesh Coarse("MeshP1C1.txt" , 1 , N_center) ;
	Coarse.generateMesh() ;
	
	middle *imaginary ;
	int imaginary_length = 3 * Coarse.number_of_elements ; 
	imaginary = new middle[imaginary_length] ;
	
	counter = 0 ; 

	for (int i = 0 ; i < Coarse.number_of_elements ; i ++){
	
		imaginary[counter].x = (Coarse.domain[i].x1_global + Coarse.domain[i].x2_global)/2.0 ; 
		imaginary[counter].y = (Coarse.domain[i].y1_global + Coarse.domain[i].y2_global)/2.0 ;
		imaginary[counter].N1 = Coarse.domain[i].node_1 ;
		imaginary[counter].N2 = Coarse.domain[i].node_2 ;

		counter ++ ;

		imaginary[counter].x = (Coarse.domain[i].x1_global + Coarse.domain[i].x3_global)/2.0 ; 
		imaginary[counter].y = (Coarse.domain[i].y1_global + Coarse.domain[i].y3_global)/2.0 ;
		imaginary[counter].N1 = Coarse.domain[i].node_1 ;
		imaginary[counter].N2 = Coarse.domain[i].node_3 ;

		counter ++ ;

		imaginary[counter].x = (Coarse.domain[i].x3_global + Coarse.domain[i].x2_global)/2.0 ; 
		imaginary[counter].y = (Coarse.domain[i].y3_global + Coarse.domain[i].y2_global)/2.0 ;
		imaginary[counter].N1 = Coarse.domain[i].node_3 ;
		imaginary[counter].N2 = Coarse.domain[i].node_2 ;

		counter ++ ;
	}

	Coarse.Destroy() ;
	//Coarse.~Mesh() ;
	
	sort (imaginary , imaginary + imaginary_length , compareMiddle) ;
	
	double x_test , y_test ;
	x_test = imaginary[0].x ;
	y_test = imaginary[0].y ;
	
	middle *total_imaginary ; 
	total_imaginary = new middle[pressure_node_number - pressure_node_number_minus_1] ;
	
	total_imaginary[0].x = imaginary[0].x ;
	total_imaginary[0].y = imaginary[0].y ;
	total_imaginary[0].N1 = imaginary[0].N1 ;
	total_imaginary[0].N2 = imaginary[0].N2 ;

	int tally = 1 ; 

	for (int i = 1 ; i < imaginary_length ; i ++){
		
		if ((x_test == imaginary[i].x) && (y_test == imaginary[i].y)){continue ;}

		total_imaginary[tally].x = imaginary[i].x ;
		total_imaginary[tally].y = imaginary[i].y ;
		total_imaginary[tally].N1 = imaginary[i].N1 ;
		total_imaginary[tally].N2 = imaginary[i].N2 ;
		tally ++ ;
		x_test = imaginary[i].x ; 
		y_test = imaginary[i].y ;
	}

	interspace *Distance ;
	Distance = new interspace[pressure_node_number - pressure_node_number_minus_1] ;
	
	counter = pressure_node_number_minus_1 ;
	int repetitious = 0 ;
	
	for (int j = 0 ; j < pressure_node_number ; j ++){
		
		double x = artificial_f_mesh[j].x_coordinate ; 
		double y = artificial_f_mesh[j].y_coordinate ;
		
		for (int k = 0 ; k < pressure_node_number_minus_1 ; k ++){
			if (x == artificial_c1_mesh[k].x_coordinate && y == artificial_c1_mesh[k].y_coordinate){repetitious = 1 ; break ;}
		}

		if (repetitious == 1){ repetitious = 0 ; continue ;}

		for(int i = 0 ; i < pressure_node_number - pressure_node_number_minus_1 ; i ++){
			
			Distance[i].D = sqrt (pow( x - total_imaginary[i].x , 2.0) + pow( y - total_imaginary[i].y , 2.0)) ; 
			Distance[i].index = i ;
		}

		sort(Distance , Distance + (pressure_node_number - pressure_node_number_minus_1) , compareDistance) ;
		
		index_V_k_minus_1[counter].i = total_imaginary[Distance[0].index].N1 ;
		index_V_k_minus_1[counter].j = total_imaginary[Distance[0].index].N2 ;
		index_V_k_minus_1[counter].k = artificial_f_mesh[j].node_number ;
		counter ++ ;
	}

	delete [] total_imaginary ;
	delete [] Distance ;
	delete [] imaginary ;

	//  Cycle == 2
	
	// Restriction
	 
	Coordinate *artificial_c2_mesh ;
	artificial_c2_mesh = new Coordinate[pressure_node_number_minus_2] ; 

	ifstream finC2("MeshP1C2.txt") ;
	finC2 >> z_stuff ;
	finC2 >> z_stuff ;
	
	for (int i = 0 ; i < pressure_node_number_minus_2 ; i ++){
		
		artificial_c2_mesh[i].node_number = i + 1 ;
		finC2 >> artificial_c2_mesh[i].x_coordinate ;  
		finC2 >> artificial_c2_mesh[i].y_coordinate ;
		finC2 >> z_stuff ;			
	}
	
	finC2.close() ;
	
	mapp_V_k_minus_2 = new restriction[pressure_node_number_minus_2] ;
	counter = 0 ; 

	for (int i = 0 ; i < pressure_node_number_minus_2 ; i ++){
		for (int j = 0 ; j < pressure_node_number ; j ++){
			
			if ((artificial_f_mesh[j].x_coordinate == artificial_c2_mesh[i].x_coordinate) && (artificial_f_mesh[j].y_coordinate == artificial_c2_mesh[i].y_coordinate)){
				
				mapp_V_k_minus_2[counter].N_coarse = artificial_c2_mesh[i].node_number ;
				mapp_V_k_minus_2[counter].N_fine = artificial_f_mesh[j].node_number ;
				counter ++ ;
				continue ;
			}
		}	
	}
	
	// Prolongation

	index_V_k_minus_2 = new prolongation[pressure_node_number_minus_1] ;
	
	restriction *mapp_V_k_minus_2_local ; 
	mapp_V_k_minus_2_local = new restriction[pressure_node_number_minus_2] ;
	
	counter = 0 ; 

	for (int i = 0 ; i < pressure_node_number_minus_2 ; i ++){
		for (int j = 0 ; j < pressure_node_number_minus_1 ; j ++){
			
			if ((artificial_c1_mesh[j].x_coordinate == artificial_c2_mesh[i].x_coordinate) && (artificial_c1_mesh[j].y_coordinate == artificial_c2_mesh[i].y_coordinate)){
				
				mapp_V_k_minus_2_local[counter].N_coarse = artificial_c2_mesh[i].node_number ;
				mapp_V_k_minus_2_local[counter].N_fine = artificial_c1_mesh[j].node_number ;
				counter ++ ;
				continue ;
			}
		}	
	}
	
	for (int i = 0 ; i < pressure_node_number_minus_2 ; i ++){

		index_V_k_minus_2[i].i = mapp_V_k_minus_2_local[i].N_coarse ;  
		index_V_k_minus_2[i].j = mapp_V_k_minus_2_local[i].N_coarse ; 
		index_V_k_minus_2[i].k = mapp_V_k_minus_2_local[i].N_fine ;
	}

	Mesh CoarseC2("MeshP1C2.txt" , 1 , N_center) ;
	CoarseC2.generateMesh() ;
	
	middle *imaginaryC2 ;
	int imaginary_length_C2 = 3 * CoarseC2.number_of_elements ; 
	imaginaryC2 = new middle[imaginary_length_C2] ;
	
	counter = 0 ; 

	for (int i = 0 ; i < CoarseC2.number_of_elements ; i ++){
	
		imaginaryC2[counter].x = (CoarseC2.domain[i].x1_global + CoarseC2.domain[i].x2_global)/2.0 ; 
		imaginaryC2[counter].y = (CoarseC2.domain[i].y1_global + CoarseC2.domain[i].y2_global)/2.0 ;
		imaginaryC2[counter].N1 = CoarseC2.domain[i].node_1 ;
		imaginaryC2[counter].N2 = CoarseC2.domain[i].node_2 ;

		counter ++ ;

		imaginaryC2[counter].x = (CoarseC2.domain[i].x1_global + CoarseC2.domain[i].x3_global)/2.0 ; 
		imaginaryC2[counter].y = (CoarseC2.domain[i].y1_global + CoarseC2.domain[i].y3_global)/2.0 ;
		imaginaryC2[counter].N1 = CoarseC2.domain[i].node_1 ;
		imaginaryC2[counter].N2 = CoarseC2.domain[i].node_3 ;

		counter ++ ;

		imaginaryC2[counter].x = (CoarseC2.domain[i].x3_global + CoarseC2.domain[i].x2_global)/2.0 ; 
		imaginaryC2[counter].y = (CoarseC2.domain[i].y3_global + CoarseC2.domain[i].y2_global)/2.0 ;
		imaginaryC2[counter].N1 = CoarseC2.domain[i].node_3 ;
		imaginaryC2[counter].N2 = CoarseC2.domain[i].node_2 ;

		counter ++ ;
	}

	CoarseC2.Destroy() ;
	//CoarseC2.~Mesh() ;
	
	sort (imaginaryC2 , imaginaryC2 + imaginary_length_C2 , compareMiddle) ;
	
	x_test = imaginaryC2[0].x ;
	y_test = imaginaryC2[0].y ;
	
	middle *total_imaginaryC2 ; 
	total_imaginaryC2 = new middle[pressure_node_number_minus_1 - pressure_node_number_minus_2] ;
	
	total_imaginaryC2[0].x = imaginaryC2[0].x ;
	total_imaginaryC2[0].y = imaginaryC2[0].y ;
	total_imaginaryC2[0].N1 = imaginaryC2[0].N1 ;
	total_imaginaryC2[0].N2 = imaginaryC2[0].N2 ;

	tally = 1 ; 

	for (int i = 1 ; i < imaginary_length_C2 ; i ++){
		
		if ((x_test == imaginaryC2[i].x) && (y_test == imaginaryC2[i].y)){continue ;}

		total_imaginaryC2[tally].x = imaginaryC2[i].x ;
		total_imaginaryC2[tally].y = imaginaryC2[i].y ;
		total_imaginaryC2[tally].N1 = imaginaryC2[i].N1 ;
		total_imaginaryC2[tally].N2 = imaginaryC2[i].N2 ;
		tally ++ ;
		x_test = imaginaryC2[i].x ; 
		y_test = imaginaryC2[i].y ;
	}

	interspace *DistanceC2 ;
	DistanceC2 = new interspace[pressure_node_number_minus_1 - pressure_node_number_minus_2] ;
	
	counter = pressure_node_number_minus_2 ;
	repetitious = 0 ;
	
	for (int j = 0 ; j < pressure_node_number_minus_1 ; j ++){
		
		double x = artificial_c1_mesh[j].x_coordinate ; 
		double y = artificial_c1_mesh[j].y_coordinate ;
		
		for (int k = 0 ; k < pressure_node_number_minus_2 ; k ++){
			if (x == artificial_c2_mesh[k].x_coordinate && y == artificial_c2_mesh[k].y_coordinate){repetitious = 1 ; break ;}
		}

		if (repetitious == 1){ repetitious = 0 ; continue ;}

		for(int i = 0 ; i < pressure_node_number_minus_1 - pressure_node_number_minus_2 ; i ++){
			
			DistanceC2[i].D = sqrt (pow( x - total_imaginaryC2[i].x , 2.0) + pow( y - total_imaginaryC2[i].y , 2.0)) ; 
			DistanceC2[i].index = i ;
		}

		sort(DistanceC2 , DistanceC2 + (pressure_node_number_minus_1 - pressure_node_number_minus_2) , compareDistance) ;
		
		index_V_k_minus_2[counter].i = total_imaginaryC2[DistanceC2[0].index].N1 ;
		index_V_k_minus_2[counter].j = total_imaginaryC2[DistanceC2[0].index].N2 ;
		index_V_k_minus_2[counter].k = artificial_c1_mesh[j].node_number ;
		counter ++ ;
	}
	
	delete [] mapp_V_k_minus_2_local ;
	delete [] total_imaginaryC2 ;
	delete [] imaginaryC2 ;
	delete [] DistanceC2 ;
	
	delete [] artificial_f_mesh ;
	delete [] artificial_c1_mesh ;
	delete [] artificial_c2_mesh ;

	// For test now

	Pro = new double *[pressure_node_number] ;
	for (int i = 0 ; i < pressure_node_number ; i ++)
		Pro[i] = new double[pressure_node_number_minus_1] ;
	
	Res = new double *[pressure_node_number_minus_1] ;
	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++)
		Res[i] = new double[pressure_node_number] ;
	
	for (int i = 0 ; i < pressure_node_number ; i ++){
		for (int j = 0 ; j < pressure_node_number_minus_1 ; j ++){
			Pro[i][j] = 0.0 ;		
		}
	}

	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
		for (int j = 0 ; j < pressure_node_number ; j ++){
			Res[i][j] = 0.0 ;
		}
	}
	
	for (int i = 0 ; i < pressure_node_number ; i ++){

			//fi_fine_grid[index_V_k_minus_1[i].k - 1] = (fi_coarse_grid_L_1[index_V_k_minus_1[i].i - 1] + fi_coarse_grid_L_1[index_V_k_minus_1[i].j - 1])/2.0 ;
			if ( index_V_k_minus_1[i].i == index_V_k_minus_1[i].j ){ Pro[index_V_k_minus_1[i].k - 1][index_V_k_minus_1[i].i - 1] = 1.0 ; continue ;}
			Pro[index_V_k_minus_1[i].k - 1][index_V_k_minus_1[i].i - 1] = 0.5 ;
			Pro[index_V_k_minus_1[i].k - 1][index_V_k_minus_1[i].j - 1] = 0.5 ;
		}

	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
		for (int j = 0 ; j < pressure_node_number ; j ++){
			Res[i][j] = Pro[j][i] ;
		}
	}

	double sume = 0 ; 

	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
		for (int j = 0 ; j < pressure_node_number ; j ++){
			sume += Res[i][j] ;
		}
		for (int j = 0 ; j < pressure_node_number ; j ++){
		Res[i][j] = Res[i][j]/sume ;
		}
		sume = 0.0 ;
	}

	virtual_velocity = new double[2*velocity_node_number] ;\
}

// -------------------- Construct Map Fine Coarse Non Uniform Grid V(k - 3) --------------------

void Computation :: constructMapFineCoarseNonUniformGridVkMinus3(){

	//  Cycle == 1

	// Restriction
	 
	Coordinate *artificial_c1_mesh ;
	artificial_c1_mesh = new Coordinate[pressure_node_number_minus_1] ; 

	double z_stuff = 0.0 ;

	ifstream finC1("MeshP1C1.txt") ;
	finC1 >> z_stuff ;
	finC1 >> z_stuff ;
	
	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
		
		artificial_c1_mesh[i].node_number = i + 1 ;
		finC1 >> artificial_c1_mesh[i].x_coordinate ;  
		finC1 >> artificial_c1_mesh[i].y_coordinate ;
		finC1 >> z_stuff ;			
	}
	
	finC1.close() ;
	
	Coordinate *artificial_f_mesh ;
	artificial_f_mesh = new Coordinate[pressure_node_number] ; 

	ifstream finF("MeshP1.txt") ;
	finF >> z_stuff ;
	finF >> z_stuff ;
	
	for (int i = 0 ; i < pressure_node_number ; i ++){
		
		artificial_f_mesh[i].node_number = i + 1 ;
		finF >> artificial_f_mesh[i].x_coordinate ;  
		finF >> artificial_f_mesh[i].y_coordinate ;
		finF >> z_stuff ;			
	}
	
	finF.close() ;
	
	mapp_V_k_minus_1 = new restriction[pressure_node_number_minus_1] ;
	int counter = 0 ; 

	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){
		for (int j = 0 ; j < pressure_node_number ; j ++){
			
			if ((artificial_f_mesh[j].x_coordinate == artificial_c1_mesh[i].x_coordinate) && (artificial_f_mesh[j].y_coordinate == artificial_c1_mesh[i].y_coordinate)){
				
				mapp_V_k_minus_1[counter].N_coarse = artificial_c1_mesh[i].node_number ;
				mapp_V_k_minus_1[counter].N_fine = artificial_f_mesh[j].node_number ;
				counter ++ ;
				continue ;
			}
		}	
	}
	
	// Prolongation
	
	index_V_k_minus_1 = new prolongation[pressure_node_number] ;
	
	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++){

		index_V_k_minus_1[i].i = mapp_V_k_minus_1[i].N_coarse ; 
		index_V_k_minus_1[i].j = mapp_V_k_minus_1[i].N_coarse ; 
		index_V_k_minus_1[i].k = mapp_V_k_minus_1[i].N_fine ; 
	}

	Mesh Coarse("MeshP1C1.txt" , 1 , N_center) ;
	Coarse.generateMesh() ;
	
	middle *imaginary ;
	int imaginary_length = 3 * Coarse.number_of_elements ; 
	imaginary = new middle[imaginary_length] ;
	
	counter = 0 ; 

	for (int i = 0 ; i < Coarse.number_of_elements ; i ++){
	
		imaginary[counter].x = (Coarse.domain[i].x1_global + Coarse.domain[i].x2_global)/2.0 ; 
		imaginary[counter].y = (Coarse.domain[i].y1_global + Coarse.domain[i].y2_global)/2.0 ;
		imaginary[counter].N1 = Coarse.domain[i].node_1 ;
		imaginary[counter].N2 = Coarse.domain[i].node_2 ;

		counter ++ ;

		imaginary[counter].x = (Coarse.domain[i].x1_global + Coarse.domain[i].x3_global)/2.0 ; 
		imaginary[counter].y = (Coarse.domain[i].y1_global + Coarse.domain[i].y3_global)/2.0 ;
		imaginary[counter].N1 = Coarse.domain[i].node_1 ;
		imaginary[counter].N2 = Coarse.domain[i].node_3 ;

		counter ++ ;

		imaginary[counter].x = (Coarse.domain[i].x3_global + Coarse.domain[i].x2_global)/2.0 ; 
		imaginary[counter].y = (Coarse.domain[i].y3_global + Coarse.domain[i].y2_global)/2.0 ;
		imaginary[counter].N1 = Coarse.domain[i].node_3 ;
		imaginary[counter].N2 = Coarse.domain[i].node_2 ;

		counter ++ ;
	}

	Coarse.Destroy() ;
	//Coarse.~Mesh() ;
	
	sort (imaginary , imaginary + imaginary_length , compareMiddle) ;
	
	double x_test , y_test ;
	x_test = imaginary[0].x ;
	y_test = imaginary[0].y ;
	
	middle *total_imaginary ; 
	total_imaginary = new middle[pressure_node_number - pressure_node_number_minus_1] ;
	
	total_imaginary[0].x = imaginary[0].x ;
	total_imaginary[0].y = imaginary[0].y ;
	total_imaginary[0].N1 = imaginary[0].N1 ;
	total_imaginary[0].N2 = imaginary[0].N2 ;

	int tally = 1 ; 

	for (int i = 1 ; i < imaginary_length ; i ++){
		
		if ((x_test == imaginary[i].x) && (y_test == imaginary[i].y)){continue ;}

		total_imaginary[tally].x = imaginary[i].x ;
		total_imaginary[tally].y = imaginary[i].y ;
		total_imaginary[tally].N1 = imaginary[i].N1 ;
		total_imaginary[tally].N2 = imaginary[i].N2 ;
		tally ++ ;
		x_test = imaginary[i].x ; 
		y_test = imaginary[i].y ;
	}

	interspace *Distance ;
	Distance = new interspace[pressure_node_number - pressure_node_number_minus_1] ;
	
	counter = pressure_node_number_minus_1 ;
	int repetitious = 0 ;
	
	for (int j = 0 ; j < pressure_node_number ; j ++){
		
		double x = artificial_f_mesh[j].x_coordinate ; 
		double y = artificial_f_mesh[j].y_coordinate ;
		
		for (int k = 0 ; k < pressure_node_number_minus_1 ; k ++){
			if (x == artificial_c1_mesh[k].x_coordinate && y == artificial_c1_mesh[k].y_coordinate){repetitious = 1 ; break ;}
		}

		if (repetitious == 1){ repetitious = 0 ; continue ;}

		for(int i = 0 ; i < pressure_node_number - pressure_node_number_minus_1 ; i ++){
			
			Distance[i].D = sqrt (pow( x - total_imaginary[i].x , 2.0) + pow( y - total_imaginary[i].y , 2.0)) ; 
			Distance[i].index = i ;
		}

		sort(Distance , Distance + (pressure_node_number - pressure_node_number_minus_1) , compareDistance) ;
		
		index_V_k_minus_1[counter].i = total_imaginary[Distance[0].index].N1 ;
		index_V_k_minus_1[counter].j = total_imaginary[Distance[0].index].N2 ;
		index_V_k_minus_1[counter].k = artificial_f_mesh[j].node_number ;
		counter ++ ;
	}

	delete [] total_imaginary ;
	delete [] Distance ;
	delete [] imaginary ;

	//  Cycle == 2
	
	// Restriction
	 
	Coordinate *artificial_c2_mesh ;
	artificial_c2_mesh = new Coordinate[pressure_node_number_minus_2] ; 

	ifstream finC2("MeshP1C2.txt") ;
	finC2 >> z_stuff ;
	finC2 >> z_stuff ;
	
	for (int i = 0 ; i < pressure_node_number_minus_2 ; i ++){
		
		artificial_c2_mesh[i].node_number = i + 1 ;
		finC2 >> artificial_c2_mesh[i].x_coordinate ;  
		finC2 >> artificial_c2_mesh[i].y_coordinate ;
		finC2 >> z_stuff ;			
	}
	
	finC2.close() ;
	
	mapp_V_k_minus_2 = new restriction[pressure_node_number_minus_2] ;
	counter = 0 ; 

	for (int i = 0 ; i < pressure_node_number_minus_2 ; i ++){
		for (int j = 0 ; j < pressure_node_number ; j ++){
			
			if ((artificial_f_mesh[j].x_coordinate == artificial_c2_mesh[i].x_coordinate) && (artificial_f_mesh[j].y_coordinate == artificial_c2_mesh[i].y_coordinate)){
				
				mapp_V_k_minus_2[counter].N_coarse = artificial_c2_mesh[i].node_number ;
				mapp_V_k_minus_2[counter].N_fine = artificial_f_mesh[j].node_number ;
				counter ++ ;
				continue ;
			}
		}	
	}
	
	// Prolongation

	index_V_k_minus_2 = new prolongation[pressure_node_number_minus_1] ;
	
	restriction *mapp_V_k_minus_2_local ; 
	mapp_V_k_minus_2_local = new restriction[pressure_node_number_minus_2] ;
	
	counter = 0 ; 

	for (int i = 0 ; i < pressure_node_number_minus_2 ; i ++){
		for (int j = 0 ; j < pressure_node_number_minus_1 ; j ++){
			
			if ((artificial_c1_mesh[j].x_coordinate == artificial_c2_mesh[i].x_coordinate) && (artificial_c1_mesh[j].y_coordinate == artificial_c2_mesh[i].y_coordinate)){
				
				mapp_V_k_minus_2_local[counter].N_coarse = artificial_c2_mesh[i].node_number ;
				mapp_V_k_minus_2_local[counter].N_fine = artificial_c1_mesh[j].node_number ;
				counter ++ ;
				continue ;
			}
		}	
	}
	
	for (int i = 0 ; i < pressure_node_number_minus_2 ; i ++){

		index_V_k_minus_2[i].i = mapp_V_k_minus_2_local[i].N_coarse ;  
		index_V_k_minus_2[i].j = mapp_V_k_minus_2_local[i].N_coarse ; 
		index_V_k_minus_2[i].k = mapp_V_k_minus_2_local[i].N_fine ;
	}

	Mesh CoarseC2("MeshP1C2.txt" , 1 , N_center) ;
	CoarseC2.generateMesh() ;
	
	middle *imaginaryC2 ;
	int imaginary_length_C2 = 3 * CoarseC2.number_of_elements ; 
	imaginaryC2 = new middle[imaginary_length_C2] ;
	
	counter = 0 ; 

	for (int i = 0 ; i < CoarseC2.number_of_elements ; i ++){
	
		imaginaryC2[counter].x = (CoarseC2.domain[i].x1_global + CoarseC2.domain[i].x2_global)/2.0 ; 
		imaginaryC2[counter].y = (CoarseC2.domain[i].y1_global + CoarseC2.domain[i].y2_global)/2.0 ;
		imaginaryC2[counter].N1 = CoarseC2.domain[i].node_1 ;
		imaginaryC2[counter].N2 = CoarseC2.domain[i].node_2 ;

		counter ++ ;

		imaginaryC2[counter].x = (CoarseC2.domain[i].x1_global + CoarseC2.domain[i].x3_global)/2.0 ; 
		imaginaryC2[counter].y = (CoarseC2.domain[i].y1_global + CoarseC2.domain[i].y3_global)/2.0 ;
		imaginaryC2[counter].N1 = CoarseC2.domain[i].node_1 ;
		imaginaryC2[counter].N2 = CoarseC2.domain[i].node_3 ;

		counter ++ ;

		imaginaryC2[counter].x = (CoarseC2.domain[i].x3_global + CoarseC2.domain[i].x2_global)/2.0 ; 
		imaginaryC2[counter].y = (CoarseC2.domain[i].y3_global + CoarseC2.domain[i].y2_global)/2.0 ;
		imaginaryC2[counter].N1 = CoarseC2.domain[i].node_3 ;
		imaginaryC2[counter].N2 = CoarseC2.domain[i].node_2 ;

		counter ++ ;
	}

	CoarseC2.Destroy() ;
	//CoarseC2.~Mesh() ;
	
	sort (imaginaryC2 , imaginaryC2 + imaginary_length_C2 , compareMiddle) ;
	
	x_test = imaginaryC2[0].x ;
	y_test = imaginaryC2[0].y ;
	
	middle *total_imaginaryC2 ; 
	total_imaginaryC2 = new middle[pressure_node_number_minus_1 - pressure_node_number_minus_2] ;
	
	total_imaginaryC2[0].x = imaginaryC2[0].x ;
	total_imaginaryC2[0].y = imaginaryC2[0].y ;
	total_imaginaryC2[0].N1 = imaginaryC2[0].N1 ;
	total_imaginaryC2[0].N2 = imaginaryC2[0].N2 ;

	tally = 1 ; 

	for (int i = 1 ; i < imaginary_length_C2 ; i ++){
		
		if ((x_test == imaginaryC2[i].x) && (y_test == imaginaryC2[i].y)){continue ;}

		total_imaginaryC2[tally].x = imaginaryC2[i].x ;
		total_imaginaryC2[tally].y = imaginaryC2[i].y ;
		total_imaginaryC2[tally].N1 = imaginaryC2[i].N1 ;
		total_imaginaryC2[tally].N2 = imaginaryC2[i].N2 ;
		tally ++ ;
		x_test = imaginaryC2[i].x ; 
		y_test = imaginaryC2[i].y ;
	}

	interspace *DistanceC2 ;
	DistanceC2 = new interspace[pressure_node_number_minus_1 - pressure_node_number_minus_2] ;
	
	counter = pressure_node_number_minus_2 ;
	repetitious = 0 ;
	
	for (int j = 0 ; j < pressure_node_number_minus_1 ; j ++){
		
		double x = artificial_c1_mesh[j].x_coordinate ; 
		double y = artificial_c1_mesh[j].y_coordinate ;
		
		for (int k = 0 ; k < pressure_node_number_minus_2 ; k ++){
			if (x == artificial_c2_mesh[k].x_coordinate && y == artificial_c2_mesh[k].y_coordinate){repetitious = 1 ; break ;}
		}

		if (repetitious == 1){ repetitious = 0 ; continue ;}

		for(int i = 0 ; i < pressure_node_number_minus_1 - pressure_node_number_minus_2 ; i ++){
			
			DistanceC2[i].D = sqrt (pow( x - total_imaginaryC2[i].x , 2.0) + pow( y - total_imaginaryC2[i].y , 2.0)) ; 
			DistanceC2[i].index = i ;
		}

		sort(DistanceC2 , DistanceC2 + (pressure_node_number_minus_1 - pressure_node_number_minus_2) , compareDistance) ;
		
		index_V_k_minus_2[counter].i = total_imaginaryC2[DistanceC2[0].index].N1 ;
		index_V_k_minus_2[counter].j = total_imaginaryC2[DistanceC2[0].index].N2 ;
		index_V_k_minus_2[counter].k = artificial_c1_mesh[j].node_number ;
		counter ++ ;
	}
	
	delete [] mapp_V_k_minus_2_local ;
	delete [] total_imaginaryC2 ;
	delete [] imaginaryC2 ;
	delete [] DistanceC2 ;
	
	// Cycle == 3
	
	// Restriction
	 
	Coordinate *artificial_c3_mesh ;
	artificial_c3_mesh = new Coordinate[pressure_node_number_minus_3] ; 

	ifstream finC3("MeshP1C3.txt") ;
	finC3 >> z_stuff ;
	finC3 >> z_stuff ;
	
	for (int i = 0 ; i < pressure_node_number_minus_3 ; i ++){
		
		artificial_c3_mesh[i].node_number = i + 1 ;
		finC3 >> artificial_c3_mesh[i].x_coordinate ;  
		finC3 >> artificial_c3_mesh[i].y_coordinate ;
		finC3 >> z_stuff ;			
	}
	
	finC3.close() ;
	
	mapp_V_k_minus_3 = new restriction[pressure_node_number_minus_3] ;
	counter = 0 ; 

	for (int i = 0 ; i < pressure_node_number_minus_3 ; i ++){
		for (int j = 0 ; j < pressure_node_number ; j ++){
			
			if ((artificial_f_mesh[j].x_coordinate == artificial_c3_mesh[i].x_coordinate) && (artificial_f_mesh[j].y_coordinate == artificial_c3_mesh[i].y_coordinate)){
				
				mapp_V_k_minus_3[counter].N_coarse = artificial_c3_mesh[i].node_number ;
				mapp_V_k_minus_3[counter].N_fine = artificial_f_mesh[j].node_number ;
				counter ++ ;
				continue ;
			}
		}	
	}

	// Prolongation
	
	index_V_k_minus_3 = new prolongation[pressure_node_number_minus_2] ;
	
	restriction *mapp_V_k_minus_3_local ; 
	mapp_V_k_minus_3_local = new restriction[pressure_node_number_minus_3] ;
	
	counter = 0 ; 

	for (int i = 0 ; i < pressure_node_number_minus_3 ; i ++){
		for (int j = 0 ; j < pressure_node_number_minus_2 ; j ++){
			
			if ((artificial_c2_mesh[j].x_coordinate == artificial_c3_mesh[i].x_coordinate) && (artificial_c2_mesh[j].y_coordinate == artificial_c3_mesh[i].y_coordinate)){
				
				mapp_V_k_minus_3_local[counter].N_coarse = artificial_c3_mesh[i].node_number ;
				mapp_V_k_minus_3_local[counter].N_fine = artificial_c2_mesh[j].node_number ;
				counter ++ ;
				continue ;
			}
		}	
	}

	for (int i = 0 ; i < pressure_node_number_minus_3 ; i ++){

		index_V_k_minus_3[i].i = mapp_V_k_minus_3_local[i].N_coarse ; 
		index_V_k_minus_3[i].j = mapp_V_k_minus_3_local[i].N_coarse ; 
		index_V_k_minus_3[i].k = mapp_V_k_minus_3_local[i].N_fine ;
	}

	Mesh CoarseC3("MeshP1C3.txt" , 1 , N_center) ;
	CoarseC3.generateMesh() ;
	
	middle *imaginaryC3 ;
	int imaginary_length_C3 = 3 * CoarseC3.number_of_elements ; 
	imaginaryC3 = new middle[imaginary_length_C3] ;
	
	counter = 0 ; 

	for (int i = 0 ; i < CoarseC3.number_of_elements ; i ++){
	
		imaginaryC3[counter].x = (CoarseC3.domain[i].x1_global + CoarseC3.domain[i].x2_global)/2.0 ; 
		imaginaryC3[counter].y = (CoarseC3.domain[i].y1_global + CoarseC3.domain[i].y2_global)/2.0 ;
		imaginaryC3[counter].N1 = CoarseC3.domain[i].node_1 ;
		imaginaryC3[counter].N2 = CoarseC3.domain[i].node_2 ;

		counter ++ ;

		imaginaryC3[counter].x = (CoarseC3.domain[i].x1_global + CoarseC3.domain[i].x3_global)/2.0 ; 
		imaginaryC3[counter].y = (CoarseC3.domain[i].y1_global + CoarseC3.domain[i].y3_global)/2.0 ;
		imaginaryC3[counter].N1 = CoarseC3.domain[i].node_1 ;
		imaginaryC3[counter].N2 = CoarseC3.domain[i].node_3 ;

		counter ++ ;

		imaginaryC3[counter].x = (CoarseC3.domain[i].x3_global + CoarseC3.domain[i].x2_global)/2.0 ; 
		imaginaryC3[counter].y = (CoarseC3.domain[i].y3_global + CoarseC3.domain[i].y2_global)/2.0 ;
		imaginaryC3[counter].N1 = CoarseC3.domain[i].node_3 ;
		imaginaryC3[counter].N2 = CoarseC3.domain[i].node_2 ;

		counter ++ ;
	}

	CoarseC3.Destroy() ;
	//CoarseC3.~Mesh() ;

	sort (imaginaryC3 , imaginaryC3 + imaginary_length_C3 , compareMiddle) ;
	
	x_test = imaginaryC3[0].x ;
	y_test = imaginaryC3[0].y ;
	
	middle *total_imaginaryC3 ; 
	total_imaginaryC3 = new middle[pressure_node_number_minus_2 - pressure_node_number_minus_3] ;
	
	total_imaginaryC3[0].x = imaginaryC3[0].x ;
	total_imaginaryC3[0].y = imaginaryC3[0].y ;
	total_imaginaryC3[0].N1 = imaginaryC3[0].N1 ;
	total_imaginaryC3[0].N2 = imaginaryC3[0].N2 ;

	tally = 1 ; 

	for (int i = 1 ; i < imaginary_length_C3 ; i ++){
		
		if ((x_test == imaginaryC3[i].x) && (y_test == imaginaryC3[i].y)){continue ;}

		total_imaginaryC3[tally].x = imaginaryC3[i].x ;
		total_imaginaryC3[tally].y = imaginaryC3[i].y ;
		total_imaginaryC3[tally].N1 = imaginaryC3[i].N1 ;
		total_imaginaryC3[tally].N2 = imaginaryC3[i].N2 ;
		tally ++ ;
		x_test = imaginaryC3[i].x ; 
		y_test = imaginaryC3[i].y ;
	}

	interspace *DistanceC3 ;
	DistanceC3 = new interspace[pressure_node_number_minus_2 - pressure_node_number_minus_3] ;
	
	counter = pressure_node_number_minus_3 ;
	repetitious = 0 ;
	
	for (int j = 0 ; j < pressure_node_number_minus_2 ; j ++){
		
		double x = artificial_c2_mesh[j].x_coordinate ; 
		double y = artificial_c2_mesh[j].y_coordinate ;
		
		for (int k = 0 ; k < pressure_node_number_minus_3 ; k ++){
			
			if (x == artificial_c3_mesh[k].x_coordinate && y == artificial_c3_mesh[k].y_coordinate){repetitious = 1 ; break ;}
		}

		if (repetitious == 1){ repetitious = 0 ; continue ;}

		for(int i = 0 ; i < pressure_node_number_minus_2 - pressure_node_number_minus_3 ; i ++){
			
			DistanceC3[i].D = sqrt (pow( x - total_imaginaryC3[i].x , 2.0) + pow( y - total_imaginaryC3[i].y , 2.0)) ; 
			DistanceC3[i].index = i ;
		}

		sort(DistanceC3 , DistanceC3 + (pressure_node_number_minus_2 - pressure_node_number_minus_3) , compareDistance) ;
		
		index_V_k_minus_3[counter].i = total_imaginaryC3[DistanceC3[0].index].N1 ;
		index_V_k_minus_3[counter].j = total_imaginaryC3[DistanceC3[0].index].N2 ;
		index_V_k_minus_3[counter].k = artificial_c2_mesh[j].node_number ;
		counter ++ ;
	}

	delete [] total_imaginaryC3 ;
	delete [] DistanceC3 ;
	delete [] imaginaryC3 ;
	delete [] mapp_V_k_minus_3_local ;

	delete [] artificial_f_mesh ;
	delete [] artificial_c1_mesh ;
	delete [] artificial_c2_mesh ;
	delete [] artificial_c3_mesh ;
}

// -------------------- Inject Intermediate Velocity Vk to Vk - 1 --------------------

void Computation :: injectIntermediateVelocityVktoVkMinus1(){
	
	for (int i = 0 ; i < pressure_node_number_minus_1 /*array_length*/  ; i ++){
		
			intermediate_velocity_coarse_grid_L_1[mapp_V_k_minus_1[i].N_coarse - 1] = intermediate_velocity[mapp_V_k_minus_1[i].N_fine - 1] ;
			intermediate_velocity_coarse_grid_L_1[mapp_V_k_minus_1[i].N_coarse - 1 + velocity_coarse_number] = intermediate_velocity[mapp_V_k_minus_1[i].N_fine - 1 + velocity_node_number] ;			
		}
}

// -------------------- Inject Intermediate Velocity V(k) to V(k - 2) --------------------

void Computation :: injectIntermediateVelocityVktoVkMinus2(){
	
	for (int i = 0 ; i < pressure_node_number_minus_2 /*array_length*/  ; i ++){

		 intermediate_velocity_coarse_grid_L_2[mapp_V_k_minus_2[i].N_coarse - 1] = intermediate_velocity[mapp_V_k_minus_2[i].N_fine - 1] ;			
		 intermediate_velocity_coarse_grid_L_2[mapp_V_k_minus_2[i].N_coarse - 1 + velocity_coarse_number] = intermediate_velocity[mapp_V_k_minus_2[i].N_fine - 1 + velocity_node_number] ;			
		}
}

// -------------------- Inject Intermediate Velocity V(k) to V(k - 3) --------------------

void Computation :: injectIntermediateVelocityVktoVkMinus3(){
	
	for (int i = 0 ; i < pressure_node_number_minus_3 /*array_length*/  ; i ++){

		 intermediate_velocity_coarse_grid_L_3[mapp_V_k_minus_3[i].N_coarse - 1] = intermediate_velocity[mapp_V_k_minus_3[i].N_fine - 1] ;
		 intermediate_velocity_coarse_grid_L_3[mapp_V_k_minus_3[i].N_coarse - 1 + velocity_coarse_number] = intermediate_velocity[mapp_V_k_minus_3[i].N_fine - 1 + velocity_node_number] ;
	  }
}

// -------------------- Prolongate Pressure V(k - 1) to V(k) --------------------

void Computation :: prolongationPressureVkMinus1toVk(){
	
	for (int i = 0 ; i < pressure_node_number ; i ++)
		fi_fine_grid[index_V_k_minus_1[i].k - 1] = (1.000)*(fi_coarse_grid_L_1[index_V_k_minus_1[i].i - 1] + fi_coarse_grid_L_1[index_V_k_minus_1[i].j - 1])/2.0 ;		
}

// -------------------- Prolongate Pressure V(k - 2) to V(k) --------------------

void Computation :: prolongationPressureVkMinus2toVk(){
	
	for (int i = 0 ; i < pressure_node_number_minus_1 ; i ++)
		fi_coarse_grid_L_1[index_V_k_minus_2[i].k - 1] = (fi_coarse_grid_L_2[index_V_k_minus_2[i].i - 1] + fi_coarse_grid_L_2[index_V_k_minus_2[i].j - 1])/2.0 ;
		
	prolongationPressureVkMinus1toVk() ;
}

// -------------------- Prolongate Pressure V(k - 3) to V(k) --------------------

void Computation :: prolongationPressureVkMinus3toVk(){
	
	for (int i = 0 ; i < pressure_node_number_minus_2 ; i ++)
			fi_coarse_grid_L_2[index_V_k_minus_3[i].k - 1] = (fi_coarse_grid_L_3[index_V_k_minus_3[i].i - 1] + fi_coarse_grid_L_3[index_V_k_minus_3[i].j - 1])/2.0 ;
	
	prolongationPressureVkMinus2toVk() ;
}

// -------------------- Distructor --------------------

Computation :: ~Computation(){

	cout << "The Computation Class is destroyed." << endl ;
}

// ===============================
// === Out of Class Functions  ===
// ===============================

// -------------------- Make a Comparison Middle Points --------------------

bool compareMiddle(const middle & a , const middle & b){
	
   if (a.x < b.x) return true;
   if (b.x < a.x) return false;

   if (a.y < b.y) return true;
   if (b.y < a.y) return false;
	
   return false ;
}

// -------------------- Find the Minimum Distance --------------------

bool compareDistance(const interspace & a , const interspace & b){

	return a.D < b.D ;
}

// -------------------- Compare Iindex --------------------

bool compareIindex(const AAIJ & a , const AAIJ & b){

   if (a.i < b.i) return true;
   if (b.i < a.i) return false;

   if (a.j < b.j) return true;
   if (b.j < a.j) return false;
	
   return false ;
}

// -------------------- Check Sparsity of AAIJ --------------------

int checkSparsity(AAIJ L[] , int J , int K , int Spicy , int NonZeroApp){

	for (int i = 0 ; i < Spicy ; i ++){
		if (L[i].i == J && L[i].j == K){ return i ;	}		
	}
	return -1 ;
}
