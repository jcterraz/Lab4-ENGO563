/*
Author: Juan Carlos Terrazas Borbon
ID:10130921
Last Update: 11/20/2017
Course: ENGO 553
Lab 4:*/

#include "Lab4.h";
void lab_1(MatrixXd &Res, MatrixXd &P, MatrixXd &A)
{
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

// Obtain weight matrix
P = compute_P_matrix(ang_data, dist_data, std_ang, std_dist);

MatrixXd S, W, Unk;
bool tresh = false;
while (tresh == false)
{
	// Compute the A and P matrix for Least Squares
	int num_unk = 4; // Given Unknowns

	compute_A_matrix(A, ang_data, dist_data, coords_data, num_unk);

	// Estimate the coordinates of points P1 and P2 and their variance-covariance matrix
	W = compute_w_matrix(ang_data, dist_data, coords_data);

	// Compute Delta
	S = -1 * (A.transpose()*P*A).inverse()*A.transpose()*P*W;
	int check = 0;

	// Check for iteration
	for (int i = 0; i < S.rows(); i++)
	{
		if (abs(S(i, 0)) > 0.0005)
		{
			check += 1;
		}
	}

	if (check == 0)
		tresh = true;


	// Obtain Unkwons and update
	MatrixXd X_unk;
	X_unk.resize(4, 1);
	X_unk(0, 0) = coords_data.at(3).x;
	X_unk(1, 0) = coords_data.at(3).y;
	X_unk(2, 0) = coords_data.at(4).x;
	X_unk(3, 0) = coords_data.at(4).y;

	Unk = X_unk + S;

	coords_data.at(3).x = Unk(0, 0);
	coords_data.at(3).y = Unk(1, 0);
	coords_data.at(4).x = Unk(2, 0);
	coords_data.at(4).y = Unk(3, 0);
}
output_matrix("Lab1_P.txt", P);
output_matrix("Lab1_A.txt", A);
output_matrix("Lab1_Unk.txt", Unk);

// Obtain Residuals
Res = A*S + W;

// Create an Observation Matrix
MatrixXd l;
l.resize(ang_data.size() + dist_data.size(), 1);

for (int i = 0; i < (ang_data.size() + dist_data.size()); i++)
{
	if (i < ang_data.size())
		l(i, 0) = ((ang_data.at(i).degrees) + (ang_data.at(i).minutes / 60) + ang_data.at(i).seconds / 3600)*PI / 180;
	else
		l(i, 0) = dist_data.at(i - 12).distance;
}

// Obtain Residuals
MatrixXd l_hat = l + Res;
};

MatrixXd snooping_method(MatrixXd v, MatrixXd P, MatrixXd A, double apriori, double Chi, double K)
{
	// Global Test
	double T = (v.transpose() * P * v)(0, 0) / apriori;
	cout << endl << "Global Test: " << T << endl;

	// Local Test 
	MatrixXd w(v.rows(), 1), Q_v;
	if (T > Chi)
	{
		cout << "Global Test Failed: Blunder Detected" << endl;
		Q_v = P.inverse() - A*(A.transpose() * P * A).inverse() * A.transpose();

		for (int i = 0; i < v.rows(); i++)
		{
			w(i, 0) = v(i, 0) / (sqrt(apriori)*sqrt(Q_v(i, i)));
			if (abs(w(i, 0)) > K)
			{
				cout << "Observation " << i + 1 << " failed" << endl;
			}
		}
		output_matrix("w.txt", w);
	}
	return Q_v;
};

MatrixXd correlation_coefficient(MatrixXd Q_v)
{
	MatrixXd corr(Q_v.rows(), Q_v.cols());

	for (int i = 0; i < Q_v.rows(); i++)
	{
		for (int j = 0; j < Q_v.cols(); j++)
		{
			corr(i,j) = Q_v(i, j) / (sqrt(Q_v(i, i)) * sqrt(Q_v(j, j)));
		}
	}

	return corr;
};