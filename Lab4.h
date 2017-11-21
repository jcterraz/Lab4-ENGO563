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

void least_squares(MatrixXd &Res, MatrixXd &P, MatrixXd &A, vector<angles> ang_data, vector<distances> dist_data, vector<coordinates> coords_data, double std_ang, double std_dist);

MatrixXd snooping_method(MatrixXd v, MatrixXd P, MatrixXd A, double apriori, double Chi, double K, bool &check, int &obs_del);

MatrixXd correlation_coefficient(MatrixXd Q_v);


#endif