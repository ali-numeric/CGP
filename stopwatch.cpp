// === Code Name: CGP/FEM  Project                                                           
// === Author: Ali A. Kashefi																  	
// === Language: C++                                                                         

// === Class Name: StopWatch														  
// === stopwatch.cpp

#include "stopwatch.h"
#include <iostream>
using namespace std ;

// -------------------- Constructor -------------------

StopWatch :: StopWatch(int random){
	
	full = 0 ;
}

// -------------------- Start Time -------------------

void StopWatch :: startTime(){

	start = clock() ;
}

// -------------------- Stop Time -------------------

void StopWatch :: stopTime(){

	stop = clock() ;
}

// -------------------- Output Time -------------------

int StopWatch :: computeCPUTime(){
	
	return stop - start ;
}

// -------------------- Reset Time -------------------

void StopWatch :: resetTime(){

	full += computeCPUTime() ;
}

// -------------------- Final Time -------------------

float StopWatch :: finalTime(){
	
	return final = ((float)full)/CLOCKS_PER_SEC ;

	cout << ((float)full)/CLOCKS_PER_SEC << endl ;
}

// -------------------- Destruction -------------------

StopWatch :: ~StopWatch(){

}
