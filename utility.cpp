#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<vector>
#include<algorithm>
#include <numeric> 
using namespace std;



/*find the max number in an array*/
template <typename T>
T FindMax(T arr[], int length)
{
	T num_max = arr[0];
    int i;

	for (int i = 0; i < length; i++)
	{
		num_max = (num_max <= arr[i]) ? arr[i] : num_max;
	}
	return num_max;
}



/*decide whether two elements exist in a same array*/
bool Pass2Element(int arr[], int length, int p, int des)
{
	bool exist1, exist2;
	
	if ((!(p == 0 || des == 0)) && (p != des))
	{
		exist1 = (find(arr, arr + length, p) != arr + length);
		exist2 = (find(arr, arr + length, des) != arr + length);

		if (p == 29 || des == 29)
		{
			return exist1 || exist2;
		}
		return exist1 && exist2;
	}
	return false;
}



/*calculate the frequency of bus service*/
double CalFreq(int v, double t)
{
	double freq;
	freq = (double) v / (2*t);
	return freq;
}



/*find the last stop id in a route*/
int FindLastId(int* route, int len)
{
	int last_id = ((find(route, route+len, 0) - route) - 1);
	return last_id;
}



/*a checking mechanism: drop duplicate node in an array*/
void DropDuplicate(vector<int> &route)
{
    auto end = route.end();
    for (auto it = route.begin(); it != end; it++) 
    {
        end = remove(it + 1, end, *it);
    }
 
    route.erase(end, route.end());
}



/*calculate the different consecutive node pairs between two routes*/
int DiffPairs(int *route1, int *route2, int size1, int size2)
{
	vector<vector<int>> pair1_vec(size1-1);
	vector<vector<int>> pair2_vec(size2-1);
	int i, j;
	int diff_num = size1 + size2 - 2;

	for (i = 0; i < size1-1; i++)
	{
		pair1_vec[i].push_back(route1[i]);
		pair1_vec[i].push_back(route1[i+1]);
	}

	for (j = 0; j < size2-1; j++)
	{
		pair2_vec[j].push_back(route2[j]);
		pair2_vec[j].push_back(route2[j+1]);
	}

	for (i = 0; i < size1-1; i++)
	{
		for (j = 0; j < size2-1; j++)
		{
			if (pair1_vec[i] == pair2_vec[j])
			{
				diff_num -= 2;
			}
		}
	}

	return diff_num;	//if two routes are the same, their distance will be 0
}


                                                             
/*calculate the average of a vector*/
double avg(vector<double> vec)
{
	float average = accumulate(vec.begin(), vec.end(), 0.0)/vec.size();
	return average;
}                                                      



/*calculate the standard deviation of a vector*/
double stdev(vector<double> vec, double mean)
{
	double stdev = 0;

	if (vec.size() >= 1)
	{
		for (int i = 0; i < vec.size(); i++)
		{
			stdev += pow((vec[i] - mean), 2);
		}
		stdev = sqrt(stdev / (vec.size()-1));
	}
	else
	{
		stdev = 0;
	}

	return stdev;
}



/*calculate the t-test of two vectors of data (assume two samples have different variances)*/
double ttest(double mean1, double sd1, double mean2, double sd2, double n1, double n2)
{
   double t_test = (mean1 - mean2) / sqrt(((sd1 * sd1) / n1) + ((sd2 * sd2) / n2));
   return t_test;
}



/*calculate the degree of freedom for two t-test samples (assumed with different variances)*/
double dof(double mean1, double sd1, double mean2, double sd2, double n1, double n2)
{
   double nom;
   double denom;

   nom = pow(((sd1 * sd1) / n1) + ((sd2 * sd2) / n2), 2);
   denom = pow(((sd1 * sd1) / n1), 2) / (n1 - 1) + pow(((sd2 * sd2) / n2), 2) / (n2 - 1);

   return nom / denom;
}
