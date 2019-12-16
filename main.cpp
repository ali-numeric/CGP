// Ali Kashefi

// The C++ code for the simulation of the 2D
// incompressible Navier-Stokes equations using
// incremental/non-incremental pressure correction
// schemes with an unstructured finite-element discretization
// (with Gmsh as a mesh generator) and
// a semi-implicit time integration method.

#define _USE_MATH_DEFINES
#include <cmath>
#include "shape.h"
#include "computation.h"
#include "cylinder.h"

using namespace std ;

int main() {
		
	int N_center = 1 ;

	double Time = 100.00 ;
	double delta_T = 0.05 ;

	double viscosity = 0.01 ;
	double density = 1.0 ;
	
	string S = "New Version" ;
	
	double alfa = 1.0 ;
	double beta_q = 1.0 ;
	
	int level = 0 ;

	StopWatch preprocessing1(1) ;
	preprocessing1.startTime() ;
		
	Cylinder C(viscosity , density , Time , delta_T , S , alfa , 0.0 , beta_q , level , N_center) ;
	C.constructNuemannDirichletBounds() ;
	C.constructTimeDependentBC() ;	
	C.initializePressure() ;
	C.initializeVelocity() ;

	preprocessing1.stopTime() ;
	preprocessing1.resetTime() ;
    float preprocess1 = preprocessing1.finalTime() ;

	C.solveTheProblemAdvancedPaper(level) ;
	C.printPressureContourTecPlot(Time) ;
	C.printVelocityContourTecPlot(Time) ;
	C.printStreamLinesTecPlot(Time) ;
	C.printTimeProcessing(preprocess1) ;
	C.printVorticityTecPlot(Time) ;
	
	return 0 ;
}
