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
	// Read needed files with observations and coordinates
	vector<angles> ang_data;
	ang_data = read_ang_files("Angles.txt");

	vector<distances> dist_data;
	dist_data = read_dist_files("Distances.txt");

	vector<coordinates> coords_data;
	coords_data = read_coords_files("Coordinates.txt");

	// Given STDs from observations
	double std_ang = 1.5;
	double std_dist = 2;

	// Perform Lab 1 adjustments
	MatrixXd Res, P, A, Qv_P1, Corr;
	least_squares(Res, P, A, ang_data, dist_data, coords_data, std_ang, std_dist);

	// Part 2 Question 1: statistical tests using data snooping method
	bool check = false;
	int obs_del;
	while (check == false)
	{
		Qv_P1 = snooping_method(Res, P, A, 1, 23.68, 2.99, check, obs_del);

		if (obs_del + 1 < ang_data.size())
		{
			ang_data.erase(ang_data.begin() + obs_del);
		}
		else
		{
			dist_data.erase(dist_data.begin() + (obs_del - ang_data.size()));
		}
		Res.resize(0,0);
		P.resize(0,0);
		A.resize(0,0);
		Qv_P1.resize(0,0); 
		Corr.resize(0,0);
		least_squares(Res, P, A, ang_data, dist_data, coords_data, std_ang, std_dist);
	}
	
	// Part 2 Question 2: correlation
	Corr = correlation_coefficient(Qv_P1);
	output_matrix("Correlation.txt", Corr);

	// Part 4

	return 0;
}