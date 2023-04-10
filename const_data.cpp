#include "data_import.cpp"

/*PARAMETERS OF ALGORITHM SETTING*/
const int P = 20;                //population size (number of individuals for each population);
const int O = 16; 				 //number of offspring generated from genetic operators

/*PARAMETERS OF PROBLEM SETTING*/
const int Z = 29;                //number of nodes without considering the dummy node 0
const int U = 23;                //number of stops (including terminals) within TSW (demand locations)
const int V = 5;                 //number of destination nodes in urban area (24-28)
const int Y = 7;                 //number of bus terminals
const int C = 1;                 //number of interchange (node L)
const double S = 1.5;            //average time for stopping at a node, i.e., 1.5 minutes
const int W = 176;               //maximum bus fleet size, i.e., 176
const int R = 10;                //maximum number of routes on the bus network, i.e., 10
const double FMIN = 0.08;        //minimum frequency of a route, 4.8 buses per hour = 0.08 buses per minute
const int SMAX = 8;              //maximum number of intermediate stops on the TSW portion of a route, i.e., 8
const int TMAX = 35;             //maximum travel time from the bus terminal in TSW to TLT interchange, including stopping time, i.e., 35 minutes


/*set of nodes without considering the dummy node 0 (Z=U+V+C)*/
const int Z_SET[Z] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29}; 

/*set of stops (including terminals) within TSW (demand locations)*/
const int U_SET[U] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};

/*set of destination nodes in urban area (24-28)*/ 
const int V_SET[V] = {24,25,26,27,28}; 

/*set of bus terminals*/
const int Y_SET[Y] = {1,7,9,14,16,20,23}; 

/*set of interchange (node L)*/
const int C_SET[C] = {29};      //use 29 instead of "L" for simplicity


/*travel demand from node i to destination e (hourly demands during peak hours)*/
const int demands[U][V] =
{
	{595,372,256,192,441},
    {167,99,97,46,162},
    {147,100,97,57,164},
    {106,57,53,29,92},
    {313,187,196,114,301},
	{270,194,179,95,336},
    {349,192,178,94,307},
    {313,191,179,97,345},
    {298,158,124,71,254},
    {104,62,48,33,101},
	{60,36,35,21,70},
    {485,337,287,142,488},
    {547,264,226,158,423},
    {196,120,92,62,177},
    {316,203,158,82,276},
	{784,425,375,258,629},
    {87,52,50,30,81},
    {237,158,147,79,211},
    {107,63,57,30,90},
    {186,98,77,55,147},
    {113,58,57,32,83},
    {104,63,51,35,87},
    {638,461,369,197,617}
};


/*TSW-TLT-DESTINATION: in-vehicle travel time (in minutes) on the shortest path between node i and j*/
double** travel_time1 = DataImport();


/*TLT-DESTINATION: in-vehicle travel time (in minutes) on the shortest path between node i and j*/
double travel_time2[V] = {28.5, 50.5, 57.5, 50.0, 29.0};
const double max_t2 = 57.5;