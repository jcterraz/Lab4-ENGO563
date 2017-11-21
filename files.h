/*
Author: Juan Carlos Terrazas Borbon
ID:10130921
Last Update: 11/20/2017
Course: ENGO 553
Lab 4:*/

#ifndef FILES_H
#define FILES_H
// Include needed libraries
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

// Struct for angles observations
struct angles{
	string id1;
	string id2;
	string id3;
	double degrees;
	double minutes;
	double seconds;
};

// Struct for distances observations
struct distances{
	string id1;
	string id2;
	double distance;
};

// Struct for given coordinates
struct coordinates{
	string id;
	double x;
	double y;
};

/* Function: read_ang_files
Purpose: This Function has the purpose of reading angles files and storing their data.
Inputs:																			Type
-file_name: file name or directory of file										[string]
*/
vector<angles> read_ang_files(string file_name);

/* Function: read_dist_files
Purpose: This Function has the purpose of reading distance files and storing their data.
Inputs:																			Type
-file_name: file name or directory of file										[string]
*/
vector<distances> read_dist_files(string file_name);

/* Function: read_coords_files
Purpose: This Function has the purpose of reading coordinates files and storing their data.
Inputs:																			Type
-file_name: file name or directory of file										[string]
*/
vector<coordinates> read_coords_files(string file_name);

/* Function: output_matrix
Purpose: This Function has the purpose of creating and output text file and output an specific matrix
Inputs:																			Type
-file_name: file name or directory of file										[string]
-Mat: Matrix to be outputed														[MatrixXd]
*/
void output_matrix(string file_name, MatrixXd Mat);

#endif