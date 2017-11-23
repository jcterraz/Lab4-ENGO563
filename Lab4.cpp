/*
Author: Juan Carlos Terrazas Borbon
ID:10130921
Last Update: 11/21/2017
Course: ENGO 563
Lab 4:*/

#include "Lab4.h";
void least_squares(MatrixXd &Res, MatrixXd &P, MatrixXd &A, vector<angles> ang_data, vector<distances> dist_data, vector<coordinates> coords_data, double std_ang, double std_dist, double apriori)
{

// Obtain weight matrix
P = compute_P_matrix(ang_data, dist_data, std_ang, std_dist, apriori);

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
		l(i, 0) = dist_data.at(i - ang_data.size()).distance;
}

// Obtain Residuals
MatrixXd l_hat = l + Res;
};

MatrixXd snooping_method(MatrixXd v, MatrixXd P, MatrixXd A, double apriori, double Chi, double K, bool &check, int &obs_del)
{
	// Global Test
	double T = (v.transpose() * P * v)(0, 0) / apriori;
	cout << endl << "Global Test: " << T << endl;

	// Local Test 
	MatrixXd w(v.rows(), 1), Q_v;
	bool global_check = false;
	int max = 0;
	bool local_check = false;
	if (T > Chi)
	{
		cout << "Global Test Failed: Blunder Detected" << endl;
		Q_v = P.inverse() - A*(A.transpose() * P * A).inverse() * A.transpose();
		
		for (int i = 0; i < v.rows(); i++)
		{
			w(i, 0) = v(i, 0) / (sqrt(apriori)*sqrt(Q_v(i, i)));
			if (i > 0)
			{
				if (abs(w(i, 0)) > abs(w(max, 0)))
				{
					max = i;
				}
			}
		}
		obs_del = max;
		cout << "Highest blunder in observation: " << max + 1 << " with a blunder of " << w(max,0)<< endl;
	}
	else
	{
		cout << "Global Test Successfull" << endl;
		Q_v = P.inverse() - A*(A.transpose() * P * A).inverse() * A.transpose();
		check = true;
		for (int i = 0; i < v.rows(); i++)
		{
			w(i, 0) = v(i, 0) / (sqrt(apriori)*sqrt(Q_v(i, i)));

			if (abs(w(i, 0)) > K)
			{
				cout << "Local Test in observation " << i + 1 << "did not pass test with a blunder of: " << w(i, 0) << endl;
				local_check = true;
			}
			if (i > 0)
			{
				if (abs(w(i, 0)) > abs(w(max, 0)))
				{
					max = i;
				}
			}
		}
		if (local_check == false)
		{
			cout << "All observations passed local test" << endl;
		}
		cout << "Highest blunder in observation: " << max + 1 << " with a blunder of " << w(max, 0) << endl;
		output_matrix("Blunders_Final.txt", w);
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