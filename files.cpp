/*
Author: Juan Carlos Terrazas Borbon
ID:10130921
Last Update: 11/20/2017
Course: ENGO 553
Lab 4:*/

#include "files.h"


vector<angles> read_ang_files(string file_name)
{
	ifstream file;
	file.open(file_name);

	// Check if file opens
	if (file.fail())
	{
		// Output error message in case of not opening the file and exit program
		cout << "ERROR: there was an error opening your angles file." << endl;
		system("pause");
		exit(1);
	}

	else
	{
		//Output message that reading the file has been successful
		cout << "Succesfully opened " << file_name << endl;

		// Create specific vectors and input the data where it corresponds.
		// Angles file
		vector<angles> angl;
		
		while (!file.eof())
		{
			angles temp;
			file >> temp.id1 >> temp.id2 >> temp.id3 >> temp.degrees >> temp.minutes >> temp.seconds;
			angl.push_back(temp);
		}
		file.close();
		return angl;
	}
};

vector<distances> read_dist_files(string file_name)
{
	ifstream file;
	file.open(file_name);

	// Check if file opens
	if (file.fail())
	{
		// Output error message in case of not opening the file and exit program
		cout << "ERROR: there was an error opening your distance file." << endl;
		system("pause");
		exit(1);
	}

	else
	{
		//Output message that reading the file has been successful
		cout << "Succesfully opened " << file_name << endl;

		// Distance file
		vector<distances> dist;
		distances temp;
		while (!file.eof())
		{
			file >> temp.id1 >> temp.id2 >> temp.distance;
			dist.push_back(temp);
		}
		file.close();
		return dist;
	}
};

vector<coordinates> read_coords_files(string file_name)
{
	ifstream file;
	file.open(file_name);

	// Check if file opens
	if (file.fail())
	{
		// Output error message in case of not opening the file and exit program
		cout << "ERROR: there was an error opening your coordinates file." << endl;
		system("pause");
		exit(1);
	}

	else
	{
		//Output message that reading the file has been successful
		cout << "Succesfully opened " << file_name << endl;

	// Coordinates file
		vector<coordinates> coords;
		coordinates temp;
		while (!file.eof())
		{
			file >> temp.id >> temp.x >> temp.y;
			coords.push_back(temp);
		}
		file.close();
		return coords;
	}
};


void output_matrix(string file_name, MatrixXd Mat)
{
	ofstream outfile;
	outfile.open(file_name);
	
	if (outfile.fail())
	{
		// Output error message in case of not opening the file and exit program
		cout << "ERROR: there was an error opening" << file_name << endl;
		system("pause");
		exit(1);
	}
	else
	{
		//Output message that reading the file has been successful
		cout << "Succesfully opened " << file_name << endl;

		outfile.precision(10);

		int rows = Mat.rows();
		int cols = Mat.cols();

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				outfile << fixed << Mat(i, j) << '\t';

				if (j == (cols - 1))
					outfile << endl;
			}
		}
	}

	outfile.close();
	return;
}