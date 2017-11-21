/*
Author: Juan Carlos Terrazas Borbon
Last Update: 09/19/2017

Purpose:
This Header file contains all the functions needed to perform a Least
Squeares Adjustment
*/

#ifndef LEAST_H
#define LEAST_H

// Needed Libraries and files needed
#include<fstream>
#include<iostream>
#include<Eigen/Dense>
#include<string>
#include<vector>
#include<cmath>

#include "files.h";
using namespace std;
using namespace Eigen;

/* Function: compute_A_matrix
Purpose: This Function modifies the size of a given A matrix and computes its values.
Inputs:																			Type
- A: Matrix A that would be resize and filled with values						[MatrixXd]
- angles: vector of struct of the angles observations							[vector<struct>]
- dist: vector of struct of the distances observations							[vector<struct>]
- coords: vector of struc of the given coordinates								[vector<struct>]
- unk: number of unkwons														[int]
*/
void compute_A_matrix(MatrixXd& A, vector<angles> angles, vector<distances> dist, vector<coordinates> coords, int unk);

/* Function: compute_P_matrix
Purpose: This Function returns a matrix which is the weight matrix (P).
Inputs:																			Type
- ang_data: vector of struct of the angles observations							[vector<struct>]
- dist_data: vector of struct of the distances observations						[vector<struct>]
- std_ang: given standard deviation of angles observations						[double]
- std_dist: given standard deviation of distances observations					[double]
*/
MatrixXd compute_P_matrix(vector<angles> ang_data, vector<distances> dist_data, double std_ang, double std_dist);

/* Function: delta_angle
Purpose: This Function returns a VectorXd, which contains the delta values of the angles observations
which are passed to the A matrix.

Inputs:																			Type
- temp_angl: angle observation													[struct]
- coords: vector of struct of the given coordinates								[vector<struct>]
- id: unknown point id (P1 or P2)												[string]
*/
VectorXd  delta_angle(angles temp_angl, vector<coordinates> coords, string id);

/* Function: delta_dist
Purpose: This Function returns a VectorXd, which contains the delta values of the distances observations
which are passed to the A matrix.

Inputs:																			Type
- temp_dist: distance observation												[struct]
- coords: vector of struct of the given coordinates								[vector<struct>]
- id: unknown point id (P1 or P2)												[string]
*/
VectorXd delta_dist(distances temp_dist, vector<coordinates> coords, string id);

/* Function: compute_w_matrix
Purpose: This Function returns a matrix which is the w matrix needed for least squares adjustment.
Inputs:																			Type
- ang_data: vector of struct of the angles observations							[vector<struct>]
- dist: vector of struct of the distances observations							[vector<struct>]
- coord_data: vector of struct of the given coordinates							[vector<struct>]
*/
MatrixXd compute_w_matrix(vector<angles> ang_data, vector<distances> dist, vector<coordinates> coord_data);

/* Function: calculate_angle
Purpose: This Function returns the calculated angle using the observation id and the given coordinates

Inputs:																			Type
- temp_angl: angle observation													[struct]
- coords: vector of struct of the given coordinates								[vector<struct>]
*/
double calculate_angle(angles temp_angl, vector<coordinates> coords);

/* Function: calculate_dist
Purpose: This Function returns the calculated distance using the observation id and the given coordinates

Inputs:																			Type
- temp_dist: distance observation												[struct]
- coords: vector of struct of the given coordinates								[vector<struct>]
*/
double calculate_dist(distances temp_dist, vector<coordinates> coords);

#endif