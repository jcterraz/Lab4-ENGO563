/*
Author: Juan Carlos Terrazas Borbon
ID:10130921
Last Update: 11/21/2017
Course: ENGO 563
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

int main() {
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
	least_squares(Res, P, A, ang_data, dist_data, coords_data, std_ang, std_dist, 1);

	// Part 2 Question 1: statistical tests using data snooping method
	bool check = false;
	int obs_del;
	int iter = 0;
	while (check == false)
	{
		iter++;
		Qv_P1 = snooping_method(Res, P, A, 1, 23.68, 2.99, check, obs_del);

		Corr = correlation_coefficient(Qv_P1);
		output_matrix("Correlation_" + std::to_string(iter)+".txt", Corr);
		if (check == true)
		{
			break;
		}

		if (obs_del + 1 < ang_data.size())
		{
			ang_data.erase(ang_data.begin() + obs_del);
		}
		else
		{
			dist_data.erase(dist_data.begin() + (obs_del - ang_data.size()));
		}

		Res.resize(0, 0);
		P.resize(0, 0);
		A.resize(0, 0);
		Qv_P1.resize(0, 0);
		Corr.resize(0, 0);
		
		least_squares(Res, P, A, ang_data, dist_data, coords_data, std_ang, std_dist, 1);
	}

	// Part 2 Question 4: variance values of the measurements
	bool check_p4 = false;
	double tolerance = 0.001;

	// Calculate A-Posteriori
	double Apos = (Res.transpose() * P * Res)(0, 0) / (Res.size() - 4);
	
	// Separate P, V and A		
	MatrixXd P_11 = P.block(0, 0, ang_data.size(), ang_data.size());
	MatrixXd P_22 = P.block(ang_data.size(), ang_data.size(), dist_data.size(), dist_data.size());

	MatrixXd A_1 = A.block(0, 0, ang_data.size(), 4);
	MatrixXd A_2 = A.block(ang_data.size(), 0, dist_data.size(), 4);

	MatrixXd V_1 = Res.block(0, 0, ang_data.size(), 1);
	MatrixXd V_2 = Res.block(ang_data.size(), 0, dist_data.size(), 1);

	int iterations = 0;
	// Iteration section
	while (check_p4 == false)
	{
		iterations++;

		MatrixXd N = A.transpose() * P * A;

		MatrixXd N_1 = A_1.transpose() * P_11 * A_1;
		MatrixXd N_2 = A_2.transpose() * P_22 * A_2;

		MatrixXd W(2, 1);
		W(0, 0) = (V_1.transpose() * P_11 * V_1)(0,0);
		W(1, 0) = (V_2.transpose() * P_22 * V_2)(0, 0);

		MatrixXd S(2, 2);
		S(0, 0) = ang_data.size() - (2 * (N_1 * N.inverse()).trace()) + (N_1 * N.inverse()*N_1 * N.inverse()).trace();
		S(0, 1) = S(1, 0) = (N_2 * N.inverse() * N_1 * N.inverse()).trace();
		S(1,1) = dist_data.size() - (2 * (N_2 * N.inverse()).trace()) + (N_2 * N.inverse()*N_2 * N.inverse()).trace();

		MatrixXd Theta = S.inverse() * W;

		// Update
		P_11 = (Apos / Theta(0, 0)) * P_11;
		P_22 = (Apos / Theta(1, 0)) * P_22;
		
		for (int i = 0; i < P_11.rows(); i++)
		{
			P(i, i) = P_11(i, i);
		}

		for (int i = 0; i < P_22.rows(); i++)
		{
			P(i + P_11.rows(), i + P_11.rows()) = P_22(i, i);
		}
		
		// Check
		if (abs(Theta(0, 0) - Theta(1, 0)) < tolerance && abs(Theta(0, 0) - Apos) < tolerance)
		{
			check_p4 = true;
		}
	}

	double var_ang = Apos / P_11(0, 0);
	double var_dist = Apos / P_22(0, 0);

	cout << endl << "Number of iterations: " << iterations << endl;
	cout << "Aposteriori:" << Apos << endl;
	cout << "Variance of angle observations (rad^2): " << var_ang << endl;
	cout << "Variance of distance observations (m^2): " << var_dist << endl;
	cout << "Standard Deviation of angle observations (sec): " << sqrt(var_ang)*180/PI*3600 << endl;
	cout << "Standard Deviation of distance observations (cm): " << sqrt(var_dist)*100 << endl;

	return 0;
}