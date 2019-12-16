// === Code Name: CGP/FEM  Project                                                            ===
// === Author: Ali A. Kashefi																  ===	
// === Language: C++                                                                          ===

// === Class Name: Computation
// === computation.h
// === Description:
// The class assembles the elements of a domain and manages LHSs &
// RHSs of the governing equations.

#ifndef COMPUTATION_H
#define COMPUTATION_H

#include "stopwatch.h"
#include "shape.h"
#include "mesh.h"
#include "mgmres.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

using namespace std ;

struct NonlinearSuper{

	int *f , *X , *Z ;
	int ak ;
	int counter ;
	bool boundary ;

	void create(){
		f = new int [ak] ;
		X = new int [ak] ;
		Z = new int [ak] ;
	}
};

struct AAIJ{
	
	double aa ;
	int i , j ;
};

struct restriction{

	// fine --> coarse
	int N_fine , N_coarse ;
};

struct prolongation{
	
	// coarse --> fine
	int k ; // index of fine
	int i , j ; // index of coarse
};

struct middle{
	
	// N1*-----middle-----*N2
	int N1 , N2 ;
	// x = (x(N1) + x(N2))/2.0
	// y = (y(N1) + y(N2))/2.0
	double x , y ; 
};

struct interspace{
	
	int index ;
	double D ; // D = distance between two points	
};

struct LocalGlobal{
	
	int Local0 , Global0 , Local1 , Global1 , Local2 , Global2 ; 
};

struct convectionIntegration{ // For P1
	
	// N(i , j) = u1*K1ij + u2*K2ij + u3*K3ij + v1*L1ij + v2*L2ij + v3*L3ij

	double ***K ;
	double ***L ;

	void createK(){
		K = new double **[3] ;
		for (int i = 0 ; i < 3 ; i ++){
			K [i]= new double *[3] ;
			for (int j = 0 ; j < 3 ; j ++)
				K[i][j] = new double[3] ;

		}
	}

	void createL(){
		L = new double **[3] ;
		for (int i = 0 ; i < 3 ; i ++){
			L [i]= new double *[3] ;
			for (int j = 0 ; j < 3 ; j ++)
				L[i][j] = new double[3] ;

		}
	}
};

struct cylinder{ // Determine nodes on the cylinder
	
	double x , y ;
	int node ;
	double Si , Co ; // Sin and Cos
	double teta ;
	double d_teta ;
};

class Computation{

	protected :

		int velocity_node_number , pressure_node_number ;
		int pressure_node_number_minus_1 , pressure_node_number_minus_2 , pressure_node_number_minus_3 ;  
		int coarsen_number ;
		int velocity_coarse_number ;

		double **velocity_u , **velocity_v , **velocity ; // **velocity = {v} + {u} in each time
		double *intermediate_velocity, *intermediate_velocity_coarse_grid_L_1 , *intermediate_velocity_coarse_grid_L_2 , *intermediate_velocity_coarse_grid_L_3 , *intermediateX ;
		double **pressure , *pressure_exp ;
		
		double delta_t ;

		double *fi_fine_grid ;
		double *fi_coarse_grid_L_1 , *fi_coarse_grid_L_2 , *fi_coarse_grid_L_3 ; // The answer of L(Fi) = Gi , will be assigned to fi_coarse_grid in each level
 
		double *U_GMRES ; 

		double *U , *V ;
		
		double *RHS_pressure_part , *RHS_mass_part , *RHS , *RHS_U , *RHS_V ;

		double *RHS_Poisson_fine , *RHS_Poisson_coarse ;
		double *RHS_real_velocity ;
		
		double **RHS_Force ; 
		double *stress ; 

		double *global_force ;

		double viscosity , density , beta_q ;
		double r ; // penalty parameter, it is really an integer
		double r2 ;
		double alfa ; // It is 1.0 or zero.

		double time ;
		int time_step ; //t_size = (int)(final_time/delta_t) + 1 ;
		int current_step ; // physcially n in formulas, computationally n - 1

		int N ; //element numbers in the domain
		int N_element ; //element numbers in the domain 
		int cycle ;

		double *LHS_Velocity_aa ;
		int *LHS_Velocity_ai , *LHS_Velocity_aj ;
		int LHS_non_zeros_elements ;

		double *Mass_Velocity_aa ; // It belongs to a block mass velocity matrix not a single mass velocity matrix
		int *Mass_Velocity_ai , *Mass_Velocity_aj ;
		int Mass_Velocity_non_zeros_elements ;

		// Mass Prime is for applying Dirichlet B.C. to the Equ. 3 (velocity correction)
		double *Mass_Prime_Velocity_aa ; // It belongs to a block mass velocity matrix not a single mass velocity matrix
		int *Mass_Prime_Velocity_ai , *Mass_Prime_Velocity_aj ;
		int Mass_Prime_Velocity_non_zeros_elements ;
		
		double *Gradient_aa ; 
		int *Gradient_ai , *Gradient_aj ;
		int Gradient_non_zeros_elements ;

		double *Divergence_aa ; 
		int *Divergence_ai , *Divergence_aj ;
		int Divergence_non_zeros_elements ;
		
		double *Reserved_LHS_aa ; // All LHS parts except the non-linear term

		double *Laplace_Pressure_V_k_aa ;
		int *Laplace_Pressure_V_k_ai , *Laplace_Pressure_V_k_aj ;
		int Laplace_Pressure_V_k_non_zeros_elements ;

		double *Laplace_Pressure_V_k_minus_1_aa ;
		int *Laplace_Pressure_V_k_minus_1_ai , *Laplace_Pressure_V_k_minus_1_aj ;
		int Laplace_Pressure_V_k_minus_1_non_zeros_elements ;

		double *Laplace_Pressure_V_k_minus_2_aa ;
		int *Laplace_Pressure_V_k_minus_2_ai , *Laplace_Pressure_V_k_minus_2_aj ;
		int Laplace_Pressure_V_k_minus_2_non_zeros_elements ;

		double *Laplace_Pressure_V_k_minus_3_aa ;
		int *Laplace_Pressure_V_k_minus_3_ai , *Laplace_Pressure_V_k_minus_3_aj ;
		int Laplace_Pressure_V_k_minus_3_non_zeros_elements ;

		double *Divergence_V_k_minus_1_aa ; 
		int *Divergence_V_k_minus_1_ai , *Divergence_V_k_minus_1_aj ;
		int Divergence_V_k_minus_1_non_zeros_elements ;

		double *Divergence_V_k_minus_2_aa ; 
		int *Divergence_V_k_minus_2_ai , *Divergence_V_k_minus_2_aj ;
		int Divergence_V_k_minus_2_non_zeros_elements ;

		double *Divergence_V_k_minus_3_aa ; 
		int *Divergence_V_k_minus_3_ai , *Divergence_V_k_minus_3_aj ;
		int Divergence_V_k_minus_3_non_zeros_elements ;
		
		convectionIntegration *NonLinear ;
		LocalGlobal *Coefficient ;

		restriction *mapp_V_k_minus_1 ;
		restriction *mapp_V_k_minus_2 ;
		restriction *mapp_V_k_minus_3 ;
		
		prolongation *index_V_k_minus_1 ;
		prolongation *index_V_k_minus_2 ;
		prolongation *index_V_k_minus_3 ;
		
		//double ** pressure_BC ;
		double ** intermediate_u_BC ;
		double ** intermediate_v_BC ;

		int intermediate_u_BC_N ; // Number of velocity nodes in B.C.
		int intermediate_v_BC_N ; // Number of velocity nodes in B.C.

		int *Dirichlet_bound_pressure , *Dirichlet_bound_velocity ;
		int Dirichlet_bound_number_pressure , Dirichlet_bound_number_velocity  ;
		int Dirichlet_bound_number_pressure_Minus1 , Dirichlet_bound_number_pressure_Minus2 , Dirichlet_bound_number_pressure_Minus3 ;  
		int *Dirichlet_bound_pressure_Minus1 , *Dirichlet_bound_pressure_Minus2 , *Dirichlet_bound_pressure_Minus3 ;  

		int *Nuemann_bound_pressure , *Nuemann_bound_velocity ;
		int Nuemann_bound_number_pressure , Nuemann_bound_number_velocity ;
		int Nuemann_bound_number_pressure_Minus1 , Nuemann_bound_number_pressure_Minus2 , Nuemann_bound_number_pressure_Minus3 ;  
		int *Nuemann_bound_pressure_Minus1 , *Nuemann_bound_pressure_Minus2 , *Nuemann_bound_pressure_Minus3 ;   

		float advection , poissonic , inject , prolong , preprocess1 , preprocess2 , multigrid ;
		
		bool mesh_type ; // structured --> true, unstructured --> false 
		int N_center ; // the number of circles inside the domain

		double *function ; //for test
		double *vector ; // for test

		double *virtual_velocity ; // for test now
		double *virtual_coarse ; // for test now

		NonlinearSuper *aN ;
		int non_zero_convection ;
		int concise ;
		cylinder *Cyl_Inf ; int Node_Cyl ;
		double *Drag , *Lift , *PressureDrag , *ViscousDrag , *PressureLift , *ViscousLift ;

		double **Pro ; // for test now
		double **Res ; // for test now

	public :
		
		Computation(double Viscosity , double Density , double Time , double Time_Step , string s , double Alfa = 1.0 , double penalty_term = 0.0 , double Beta_q = 3.0/2.0 , int level = 0 , int n_center = 0) ;

		// ==== Pre-Processing ====
		
		void constructGlobalGradientOperatorMatrixSuper() ;
		void constructGlobalDivergenceOperatorMatrixSuper() ;
		void constructGlobalLaplacePressureVkSuper() ;
		void constructGlobalVelocityMassMatrixSuper() ;
		void constructNonlinearPaper() ;
		void constructGlobalLaplaceVelocityMatrixPaper() ;
		void calculateGlobalPressureLaplaceOperatorMatrixInCoarseGridPaper(int level) ;
		void calculateGlobalDivergenceOperatorMatrixInCoarseGridPaper(int level) ;
		void constructGlobalForceVector() ;
		void constructGlobalConvectiveMatrix() ;
		void constructTheRestOfTheMatricesVectors() ;
		void constructTimeDependentBC() ;

		// ==== Processing ====
		
		void solveTheProblemAdvancedPaper(int level = 0) ;
		void solveTheProblemAdvancedVkPaper() ; 
		void solveTheProblemAdvancedVkMinus1Paper() ;
		void solveTheProblemAdvancedVkMinus2Paper() ;
		void solveTheProblemAdvancedVkMinus3Paper() ;
		void breakRHS() ;
		void exchangeVariables() ;

		void updateRHSPressurePart() ; // n is the current time step. In fact (n+1)
		void updateRHSMassPart() ; // n is the current time step. In fact (n+1)
		void updateRHS() ;
		void updatePressure() ;
		void updatePressureExp() ;
		void calculateRHSPoissonInFineGrid() ;
		void calculateRHSPoissonInCoarseGridVkMinus1() ;
		void calculateRHSPoissonInCoarseGridVkMinus2() ;
		void calculateRHSPoissonInCoarseGridVkMinus3() ;
		void calculateRHSRealVelocity() ;
		void goToConvectionAdvanced() ;
		void updateConvectionAdvancedSuper(double Virtual_Velocity[]) ;
		
		void solveEq1(int itr_max , int rpt_max) ;
		void solveEq1primary(int itr_max , int rpt_max) ;
		void solveEq1secondary(int itr_max, int rpt_max) ;
		void solveEq2(int itr_max , int rpt_max) ;
		void solveEq2VkMinus1(int itr_max , int rpt_max) ;
		void solveEq2VkMinus2(int itr_max , int rpt_max) ;
		void solveEq2VkMinus3(int itr_max , int rpt_max) ;
		void solveEq3(int itr_max , int rpt_max) ;
		
		void applyMixedBCtoRHSofPoissonVk() ;
		void applyMixedBCtoRHSofPoissonVkMinus1() ;
		void applyMixedBCtoRHSofPoissonVkMinus2() ;
		void applyMixedBCtoRHSofPoissonVkMinus3() ;
		void applyMixedBCtoRHSofMomentum() ;
		void applyMixedBCToRHSofCorrectionEqu() ;
		
		void transferGMREStoVectors() ;
		void productAb(double A_aa[] , int A_ai[] , int A_aj[] , double b[] , double result[] , double coefficient , int N) ;
		void updateCSR(int r , int s , double variable , double A_aa[] , int A_ai[] , int A_aj[] , int A_non_zero_elements) ;
		int searchCSR(int r , int s , int A_ai[] , int A_aj[] , int A_non_zeros_elements) ;
		void zeroCSR(double A_aa[] , int A_non_zeros_elements) ;
		bool checkBoundaryPoint(int i) ;
		void convertAAIJtoCSR(int size_i , int non_zeros_elements , AAIJ A[] , double AA_aa[] , int AA_ai[] , int AA_aj[]);
		int updateAAIJ(AAIJ A[] , double value , int I , int J , int spicy , int non_zero_approximation) ;
		int imposeAAIJ( AAIJ A[] , double value , int I , int J , int spicy , int non_zero_approximation) ;
		
		double residual(double v_previous[] , double v_current[] , int N_size) ;
		
		// ==== Post-Processing ====
		
		void printVelocityContourTecPlot(double desired_time) ;
		void printPressureContourTecPlot(double desired_time) ;
		void printVorticityTecPlot(double desired_time) ;
		void printStreamLinesTecPlot(double desired_time) ;
		void printVariableFunction(double x , double y , double desired_time) ;
		void printDragAndLift() ;
		void printTimeProcessing(float T) ;
	
		void computeStress() ;
		void computeDragAndLift() ;

		//==== Coarse Grid Projection Method ====

		void constructMapFineCoarseNonUniformGrid() ;
		void constructMapFineCoarseNonUniformGridVkMinus1() ;
		void constructMapFineCoarseNonUniformGridVkMinus2() ;
		void constructMapFineCoarseNonUniformGridVkMinus3() ;
		
		void injectIntermediateVelocityVktoVkMinus1() ;
		void injectIntermediateVelocityVktoVkMinus2() ;
		void injectIntermediateVelocityVktoVkMinus3() ;

		void prolongationPressureVkMinus1toVk() ;
		void prolongationPressureVkMinus2toVk() ;
		void prolongationPressureVkMinus3toVk() ;

		~Computation() ;
};

#endif




