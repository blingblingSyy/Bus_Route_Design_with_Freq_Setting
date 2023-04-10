#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <typeinfo>
#include <algorithm>
 
using namespace std;


/*Import data from the csv file and store the data to a 2D array*/
double** DataImport()	 //relative path
{
	string fname = "tin shui wei time matrix.CSV";
	vector<vector<string>> content;
	vector<string> row;
	string line, word;
 
	fstream file (fname, ios::in);
	if(file.is_open())
	{
		while(getline(file, line))
		{
			row.clear();
 
			stringstream str(line);     //input flow
 
			while (getline(str, word, ','))
			{
				stringstream ss(word);
				string snum;
				while (getline(ss, snum, '"'))
				{
					stringstream ss1(snum);
					string num;				
					while (getline(ss1, num, '}'))
					{
						row.push_back(num);	
					}
				}
            }
			content.push_back(row);
		}
	}
	else
		cout<<"Could not open the file" << endl;

	int size = content.size() - 1;
	double** arr = new double* [size];
	for (int i = 0; i < size; i++)
	{
		arr[i] = new double [size];
		for (int j = 0; j < size; j++)
		{
			arr[i][j] = stod(content[i+1][j+1]);
		}
	}

	// for (int i = 0; i < size; i++)
	// {
	// 	for (int j = 0; j < size; j++)
	// 	{
	// 		cout << arr[i][j] << " ";
	// 	}
	// 	cout << endl;
	// }

	return arr;
}


// int main()
// {
// 	DataImport();
// 	return 0;
// }