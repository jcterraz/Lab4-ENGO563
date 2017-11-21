/*
Author: Juan Carlos Terrazas Borbon
ID:10130921
Last Update: 11/20/2017
Course: ENGO 553
Lab 4:*/

// Include Libraries needed
#include<Eigen/Dense>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

// Include Files created
#include "LeastSquares.h"
#include "files.h"
#include "Lab4.h"

int main(){
	// Perform Lab 1 adjustments
	MatrixXd Res, P, A, Qv_P1, Corr;
	lab_1(Res, P, A);

	// Part 2 Question 1: statistical tests using data snooping method
	Qv_P1 = snooping_method(Res, P, A, 1, 23.68,2.99);
	
	// Part 2 Question 2: correlation
	Corr = correlation_coefficient(Qv_P1);

	// Part 4

	return 0;
}