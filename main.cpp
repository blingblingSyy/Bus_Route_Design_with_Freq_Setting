#include<iostream>
#include <cstdio>
#include <ctime>
#include <fstream>
#include "genetic_algorithm.cpp"

using namespace std;


int main()
{
    srand((unsigned)time(NULL));
	int run_time = 10;	//set new parameters: run times
	int gen_num = 100;  //set new parameters: iteration times
	
	vector<double> transfer1;
	vector<double> transfer2;
	vector<double> travel_time1;
	vector<double> travel_time2;
	vector<double> obj_v1;
	vector<double> obj_v2;


	/*1. start running the simultaneous genetic algorithm*/
	ofstream outfile1;
	outfile1.open("result.dat", ios::out | ios::trunc);

	clock_t start1;
    double duration1;
    start1 = clock(); // get current time

	Genetic_Algorithm GA1;
	GA1.Set_Run_Times(run_time);
	GA1.Set_Gens_Num(gen_num);
	GA1.Set_Simultaneous(true);
	GA1.GA();

	outfile1 << "* Results of Simultaneous Genetic Algorithm: "<< endl;
	for (int i = 0; i < GA1.best_inds.size(); i++)
	{
		transfer1.push_back(GA1.best_inds[i].transfer);
		travel_time1.push_back(GA1.best_inds[i].travel_time);
		obj_v1.push_back(GA1.best_inds[i].obj_value);
		outfile1 << "Best Individual of the " << i+1 << "-th run time: ";
		outfile1 << "transfer times = " << GA1.best_inds[i].transfer << ", ";
		outfile1 << "total travel time = " << GA1.best_inds[i].travel_time << ", ";
		outfile1 << "objective value = " << GA1.best_inds[i].obj_value << endl;
	}
	outfile1 << endl;

	double mean_transfer1 = avg(transfer1);
	double sd_transfer1 = stdev(transfer1, mean_transfer1);
	double mean_time1 = avg(travel_time1);
	double sd_time1 = stdev(travel_time1, mean_time1);
	double mean_obj1 = avg(obj_v1);
	double sd_obj1 = stdev(obj_v1, mean_obj1);
	outfile1 << "The average transfer is " << mean_transfer1 << endl;
	outfile1 << "The standard deviation of transfer is " << sd_transfer1 << endl;
	outfile1 << "The average travel time is " << mean_time1 << endl;
	outfile1 << "The standard deviation of travel time is " << sd_time1 << endl;
	outfile1 << "The average objective value is " << mean_obj1 << endl;
	outfile1 << "The standard deviation of objective value is " << sd_obj1 << endl;
	outfile1 << endl;

	duration1 = (clock() - start1 ) / (double) CLOCKS_PER_SEC;
    outfile1 << "Algorithm took an average of "<< double(duration1)/run_time << " seconds per run." << endl;
	outfile1 << "----------------------------------------------------------------------------------" << endl;
	outfile1 << endl;


	/*2. start running the sequential genetic algorithm*/
	clock_t start2;
    double duration2;
    start2 = clock(); // get current time

	Genetic_Algorithm GA2;
	GA2.Set_Run_Times(run_time);
	GA2.Set_Gens_Num(gen_num);
	GA2.Set_Simultaneous(false);
	GA2.GA();

	outfile1 << "* Results of Sequential Genetic Algorithm: "<< endl;
	for (int i = 0; i < GA2.best_inds.size(); i++)
	{
		transfer2.push_back(GA2.best_inds[i].transfer);
		travel_time2.push_back(GA2.best_inds[i].travel_time);
		obj_v2.push_back(GA2.best_inds[i].obj_value);
		outfile1 << "Best Individual of the " << i+1 << "-th run time: ";
		outfile1 << "transfer times = " << GA2.best_inds[i].transfer << ", ";
		outfile1 << "total travel time = " << GA2.best_inds[i].travel_time << ", ";
		outfile1 << "objective value = " << GA2.best_inds[i].obj_value << endl;
	}
	outfile1 << endl;

	double mean_transfer2 = avg(transfer2);
	double sd_transfer2 = stdev(transfer2, mean_transfer2);
	double mean_time2 = avg(travel_time2);
	double sd_time2 = stdev(travel_time2, mean_time2);
	double mean_obj2 = avg(obj_v2);
	double sd_obj2 = stdev(obj_v2, mean_obj2);
	outfile1 << "The average transfer is " << mean_transfer2 << endl;
	outfile1 << "The standard deviation of transfer is " << sd_transfer2 << endl;
	outfile1 << "The average travel time is " << mean_time2 << endl;
	outfile1 << "The standard deviation of travel time is " << sd_time2 << endl;
	outfile1 << "The average objective value is " << mean_obj2 << endl;
	outfile1 << "The standard deviation of objective value is " << sd_obj2 << endl;
	outfile1 << endl;

	duration2 = (clock() - start2 ) / (double) CLOCKS_PER_SEC;
    outfile1 << "Algorithm took an average of "<< double(duration2)/run_time << " seconds per run." << endl;
	outfile1 << "----------------------------------------------------------------------------------" << endl;
	outfile1 << endl;

	double ttest_objv = ttest(mean_obj1, sd_obj1, mean_obj2, sd_obj2, obj_v1.size(), obj_v2.size());
	double dof_objv = dof(mean_obj1, sd_obj1, mean_obj2, sd_obj2, obj_v1.size(), obj_v2.size());
	outfile1 << "* t-stat value of the objective of the two algorithms: " << ttest_objv << endl;
	outfile1 << "degree of freedom of the two samples = " << dof_objv << endl;

	outfile1.close();


	system("pause");
    return 0;
}

