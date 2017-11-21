/*
Author: Juan Carlos Terrazas Borbon
ID:10130921
Last Update: 11/20/2017
Course: ENGO 553
Lab 4:*/
#ifndef LAB_1
#define LAB_1

// include needed libraries, files and constants
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <vector>

#include "LeastSquares.h"
#include "files.h"

using namespace std;
using namespace Eigen;
const double PI = 3.141592653589793;

void lab_1(MatrixXd &Res, MatrixXd &P, MatrixXd &A);

MatrixXd snooping_method(MatrixXd v, MatrixXd P, MatrixXd A, double apriori, double Chi, double K);

MatrixXd correlation_coefficient(MatrixXd Q_v);


#endif