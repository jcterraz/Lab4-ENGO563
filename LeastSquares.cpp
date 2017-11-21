/*
Author: Juan Carlos Terrazas Borbon
ID:10130921
Last Update: 11/21/2017
Course: ENGO 563
Lab 1:*/
#include "LeastSquares.h"
const double PI = 3.141592653589793;

void compute_A_matrix(MatrixXd& A, vector<angles> ang_data, vector<distances> dist, vector<coordinates> coords, int unk)
{
	A.resize(ang_data.size() + dist.size(), unk);

	for (int i = 0; i < ang_data.size(); i++)
	{
		if (ang_data.at(i).id1 == "P1" || ang_data.at(i).id2 == "P1" || ang_data.at(i).id3 == "P1")
		{	
			VectorXd temp = delta_angle(ang_data.at(i), coords, "P1");
			A(i, 0) = temp(0);
			A(i, 1) = temp(1);
		}
		else
		{
			A(i, 0) = 0;
			A(i, 1) = 0;
		}

		if (ang_data.at(i).id1 == "P2" || ang_data.at(i).id2 == "P2" || ang_data.at(i).id3 == "P2")
		{
			VectorXd temp = delta_angle(ang_data.at(i), coords, "P2");
			A(i, 2) = temp(0);
			A(i, 3) = temp(1);
		}
		else
		{
			A(i, 2) = 0;
			A(i, 3) = 0;
		}
	}

	for (int i = 0; i < dist.size(); i++)
	{
		if (dist.at(i).id1 == "P1" || dist.at(i).id2 == "P1")
		{
			VectorXd temp = delta_dist(dist.at(i), coords, "P1");
			A(i + ang_data.size(), 0) = temp(0);
			A(i + ang_data.size(), 1) = temp(1);
		}
		else
		{
			A(i + ang_data.size(), 0) = 0;
			A(i + ang_data.size(), 1) = 0;
			
		}
		if (dist.at(i).id1 == "P2" || dist.at(i).id2 == "P2")
		{
			VectorXd temp = delta_dist(dist.at(i), coords, "P2");
			A(i + ang_data.size(), 2) = temp(0);
			A(i + ang_data.size(), 3) = temp(1);
		}
		else
		{
			A(i + ang_data.size(), 2) = 0;
			A(i + ang_data.size(), 3) = 0;
		}
	}
};

VectorXd  delta_angle(angles temp_angl, vector<coordinates> coords, string id)
{
	VectorXd temp(2);

	coordinates temp_i, temp_j, temp_k;

	for (int z = 0; z < coords.size(); z++)
	{
		if (coords.at(z).id == temp_angl.id2)
		{
			temp_i = coords.at(z);
		}
		
		else if (coords.at(z).id == temp_angl.id1)
		{
			temp_j = coords.at(z);
		}
		
		else if (coords.at(z).id == temp_angl.id3)
		{
			temp_k = coords.at(z);
		}
	}

	if (temp_angl.id2 == id)
	{
		temp(0) = ((temp_k.y - temp_i.y) / (pow(temp_k.x - temp_i.x, 2) + pow(temp_k.y - temp_i.y, 2)))
			- ((temp_j.y - temp_i.y) / (pow(temp_j.x - temp_i.x, 2) + pow(temp_j.y - temp_i.y, 2)));
		temp(1) = -1*(temp_k.x - temp_i.x) / (pow(temp_k.x - temp_i.x, 2) + pow(temp_k.y - temp_i.y, 2))
			+ ((temp_j.x - temp_i.x) / (pow(temp_j.x - temp_i.x, 2) + pow(temp_j.y - temp_i.y, 2)));
	}
	else if (temp_angl.id1 == id)
	{
		temp(0) = (temp_j.y - temp_i.y) / (pow(temp_j.x - temp_i.x, 2) + pow(temp_j.y - temp_i.y, 2));
		temp(1) = -1*(temp_j.x - temp_i.x) / (pow(temp_j.x - temp_i.x, 2) + pow(temp_j.y - temp_i.y, 2));
	}
	else if (temp_angl.id3 == id)
	{
		temp(0) = -1*(temp_k.y - temp_i.y) / (pow(temp_k.x - temp_i.x, 2) + pow(temp_k.y - temp_i.y, 2));
		temp(1) = (temp_k.x - temp_i.x) / (pow(temp_k.x - temp_i.x, 2) + pow(temp_k.y - temp_i.y, 2));
	}

	return temp;
};

VectorXd delta_dist(distances temp_dist, vector<coordinates> coords, string id)
{
	VectorXd temp(2);

	coordinates temp_i, temp_j;

	for (int z = 0; z < coords.size(); z++)
	{
		if (coords.at(z).id == temp_dist.id1)
		{
			temp_i = coords.at(z);
		}

		else if (coords.at(z).id == temp_dist.id2)
		{
			temp_j = coords.at(z);
		}

	}

	if (temp_dist.id1 == id)
	{
		temp(0) = (temp_i.x - temp_j.x) / sqrt(pow(temp_i.x - temp_j.x, 2) + pow(temp_i.y - temp_j.y, 2));
		temp(1) = (temp_i.y - temp_j.y) / sqrt(pow(temp_i.x - temp_j.x, 2) + pow(temp_i.y - temp_j.y, 2));
	}
	else if (temp_dist.id2 == id)
	{
		temp(0) = -1 * (temp_i.x - temp_j.x) / sqrt(pow(temp_i.x - temp_j.x, 2) + pow(temp_i.y - temp_j.y, 2));
		temp(1) = -1 * (temp_i.y - temp_j.y) / sqrt(pow(temp_i.x - temp_j.x, 2) + pow(temp_i.y - temp_j.y, 2));
	}

	return temp;
};

MatrixXd compute_P_matrix(vector<angles> ang_data, vector<distances> dist_data, double std_ang, double std_dist, double apriori)
{
	MatrixXd C_l = MatrixXd::Zero((ang_data.size() + dist_data.size()), (ang_data.size() + dist_data.size()));
	for (int i = 0; i < (ang_data.size() + dist_data.size()); i++)
	{
		if (i < ang_data.size())
			C_l(i, i) = pow(std_ang/3600/180*PI,2);
		else
			C_l(i, i) = pow(std_dist/100,2);
	}
	MatrixXd P = pow(apriori,2) * C_l.inverse();
	return P;
};


MatrixXd compute_w_matrix(vector<angles> ang_data, vector<distances> dist, vector<coordinates> coord_data)
{
	MatrixXd W;
	W.resize(ang_data.size() + dist.size(), 1);
	for (int i = 0; i < ang_data.size(); i++)
	{
		W(i, 0) = calculate_angle(ang_data.at(i), coord_data) - ((ang_data.at(i).degrees * 3600) + (ang_data.at(i).minutes * 60) + ang_data.at(i).seconds)/3600/180*PI;
	}
	for (int i = 0; i < dist.size(); i++)
	{
		W(i + ang_data.size(), 0) = (calculate_dist(dist.at(i), coord_data) - dist.at(i).distance);
	}
	return W;
}

double calculate_angle(angles temp_angl, vector<coordinates> coords)
{
	coordinates temp_i, temp_j, temp_k;

	for (int z = 0; z < coords.size(); z++)
	{
		if (coords.at(z).id == temp_angl.id2)
		{
			temp_i = coords.at(z);
		}

		else if (coords.at(z).id == temp_angl.id1)
		{
			temp_j = coords.at(z);
		}

		else if (coords.at(z).id == temp_angl.id3)
		{
			temp_k = coords.at(z);
		}
	}

	double angle_1 = atan2((temp_k.y - temp_i.y), (temp_k.x - temp_i.x));
	if (angle_1 < 0)
	{
		angle_1 = angle_1 + (PI * 2);
	}

	double angle_2 = atan2((temp_j.y - temp_i.y), (temp_j.x - temp_i.x));
	if (angle_2 < 0)
	{
		angle_2 = angle_2 + (PI * 2);
	}

	double angle;
	if (angle_1 > angle_2)
		angle = abs(angle_1 - angle_2);
	else
		angle = abs((PI * 2) - (angle_2 - angle_1));

	return angle;
}

double calculate_dist(distances temp_dist, vector<coordinates> coords)
{
	coordinates temp_i, temp_j;

	for (int z = 0; z < coords.size(); z++)
	{
		if (coords.at(z).id == temp_dist.id1)
		{
			temp_i = coords.at(z);
		}

		else if (coords.at(z).id == temp_dist.id2)
		{
			temp_j = coords.at(z);
		}

	}

	double dist = sqrt(pow((temp_j.y - temp_i.y), 2) + pow(temp_j.x - temp_i.x, 2));

	return dist;
}