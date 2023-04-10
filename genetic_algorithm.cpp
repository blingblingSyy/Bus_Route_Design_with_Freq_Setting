#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<vector>
#include<algorithm>
#include <Windows.h>
#include "utility.cpp"
#include "const_data.cpp"

using namespace std;


/*DEFINITION OF CLASS Genetic_Algorithm*/
class Genetic_Algorithm
{
protected:
	double CROSS_RATE = 0.8;  								//crossover rate
	double MUTA_RATE = 0.2;	 								//mutation rate
	double ROUTE_CPROB = 0.5;        						//rate of selecting route crossover operator 
	double STOP_CPROB = 0.5;         						//rate of selecting stop crossover operator 
	double INSERT_MPROB = 0.4;       						//rate of selecting insert mutation operator 
	double DELETE_MPROB = 0.4;       						//rate of selecting delete mutation operator 
	double SWAP_MPROB = 0.1;         						//rate of selecting swap mutation operator 
	double TRANSFER_MPROB = 0.1;     						//rate of selecting transfer mutation operator 
	int IT = 2000;             								//iteration times(generations) of a GA algorithm
	int RT = 5;				 								//run times for each GA
	int B1 = 60;                     						//weight for the number of transfers 
	const int B2 = 1;                   					//weight for the total travel time 
	double Freq_Penalty = 500;		 						//a non-negative parameter divided by the frequency to increase the objective value
	int Time_Penalty = 200;	 		 						//a penalty added to the average travel time for solution that leave some nodes unvisited	
	bool diversity = true;									//to use diversity control mechanism
	bool sim = true;										//simultaneous genetic algorithm (false = sequential genetic algorithm)

public:
	void Set_Cross_Rate(double rate);						//set the crossover trigger rate for each individual
	void Set_Muta_Rate(double rate);						//set the mutation trigger rate for each individual
	void Set_Cross_Combine_Rate(double route_rate, double stop_rate);
															//set the combinative rate of crossover methods
	void Set_Muta_Combine_Rate(double insert_rate, double delete_rate, double swap_rate, double transfer_rate);
															//set the combinative rate of mutation methods 
	void Set_Gens_Num(int generations);						//set the number of iterations/generations for each run time
	void Set_Run_Times(int run_times);						//set algorithm run times
	void Set_B1_Weight(int weight);							//set B1 weight
	void Set_Freq_Penalty(double penl);						//set the frequency penalty value
	void Set_Time_Penalty(int penl);						//set the time penalty value
	void Set_Diversity(bool div);							//decide whether to use diversity control mechanism
	void Set_Simultaneous(bool simul);						//decide whether to use simultaneous or sequential algorithm

protected:
	int population[P][R][SMAX+2] = {0};						//initial population of #P individuals with #R routes, each route consisting of at most TMAX+1 stops
	int iter;												//an iteration variable
	int run;												//a run time variable
	int bus_array1[P][R]= {0};								//the bus allocation of the initial population
	double bus_fn1[P][R] = {0.0};							//the service frequency of the initial population
	double Fit1[P];											//the fitness value for each solution in the initial population
	int bus_array2[O][R] = {0};								//the bus allocation of the repaired mutated population
	double bus_fn2[O][R] = {0.0};							//the service frequency of the repaired mutated population
	double Fit2[O];											//the fitness value fo each solution in the mutated population
	double Rec[P];											//the reciprocal of fitness value for each solution in the initial population to select parents
	double pr[P];											//the probability of selecting each solution by roulette wheel selection
	int roul_population[P][R][SMAX+2] = {0};				//select the parents for crossover through routelette wheel selection
	int cross_population[O][R][2*SMAX+2] = {0};				//the crossover children population
	int muta_population[O][R][2*SMAX+2] = {0};				//the mutated children population
	int repair_population[O][R][SMAX+2] = {0};				//the repaired mutated children population
	int survival_population[P][R][SMAX+2] = {0};			//the mixed population through survival selection
	struct Best_Individual;									//store the results of the best individual of each iteration

public:
	vector<Best_Individual> best_inds;						//store the results of the best individuals in all iterations

public:
	void GA();												//the main body of the simultaneous genetic algorithm
	void Init_Population(int (*pop)[R][SMAX+2], int pop_size = P);
															//Initialization
	double Trip_Time(int* input_route, int len);			//calculate the trip time of a route
	bool Direct_Bus(int n1, int n2, int (*ind)[SMAX+2]);	//0-1 whether there is direct bus service from node n1 to node n2 in an individual (1: no direct bus service) 
	double Travel_Tn(int p, int des, int *route, int route_len);			
															//calculate the travel time of route n from node n1 to node n2
	double Avg_Travel_T(int p, int des, double* fn, int (*ind)[SMAX+2]);	
															//calculate the average travel time of passenges from node p to destination among all 10 routes
	double Total_Transfer(int (*ind)[SMAX+2]);				//calculate the total transfer times of a solution
	double Total_Travel_T(double* fn, int (*ind)[SMAX+2]);	//calculate the total travel time of a solution
	void Freq_Setting(int *p, int sum_w, int low, int up);	//Initialize the bus frequency setting for each route in an individual
	double Freq_Heuristic(int (*input_ind)[SMAX+2], int bus_array[10], double bus_fn[10]);
															//Frequency Setting Heuristic to improve total travel time
	double Obj_Value_Sim(int (*input_ind)[SMAX+2], int bus_array[10], double bus_fn[10]);
	double Obj_Value_Seq(int (*input_ind)[SMAX+2]);
															//calculate the objective value of a given individual
	double Fitness(double obj_v); 							//Fitness Evaluation
	void Parents_Selection();   							//Roulette Wheel Selection approach to select parents for crossover
	void Route_Crossover(int (*father)[SMAX+2], int (*mother)[SMAX+2], int (*route_child_ind1)[2*SMAX+2], int (*route_child_ind2)[2*SMAX+2]);
															//route crossover operator, need to initialize the route child individual
	void Stop_Crossover(int (*father)[SMAX+2], int (*mother)[SMAX+2], int (*stop_child_ind1)[2*SMAX+2], int (*stop_child_ind2)[2*SMAX+2]);
															//stop crossover operator, need to initialize the stop child individual
	void Insert_Mutation(int (*ind)[2*SMAX+2]);  			//insert mutation operator
	void Delete_Mutation(int (*ind)[2*SMAX+2]);  			//delete mutation operator
	void Swap_Mutation(int (*ind)[2*SMAX+2]);  				//swap mutation operator
	void Transfer_Mutation(int (*ind)[2*SMAX+2]);  			//transfer mutation operator
	void Crossover();										//integrated crossover operator
	void Mutation();										//integrated mutation operator
	void Stop_Seq_Improve(int *route, int route_len);		//stop sequence improvement heuristic
	void Repair(int *route, int route_len);					//solution Repairing to satisfy the constraints of intermediate stops and travel time
	void Offspring_Generate();  							//offspring generation by genetic operator (crossover and mutation)
	double Diversity_Control(int (*ind1)[SMAX+2], int (*ind2)[SMAX+2], double param_c = 0.08, double param_a = 0.002);
															//a new diversity control mechanism - calculate new hamming distance 
	bool Diff_Individuals(int (*ind1)[SMAX+2], int (*ind2)[SMAX+2]);
															//find whether two individuals are the same
	void Mix_Population(int mix_population[P+O][R][SMAX+2], double mix_fit[P+O], int mix_bus_array[P+O][R], 
						double mix_bus_fn[P+O][R], double sort_mix_fit[P+O], int sort_fit_preindex[P+O]);
															//mix the initial population and the repaired population
	void Best_Solution(int best_ind[R][SMAX+2], double best_fit, int bus_arr[R], double bus_fn[R]);
															//store and present the best individual of each iteration
	void Freq_Adjust(int bus_arr[R], double bus_fn[R]);		//adjust the bus allocation if frequency constraint is not satisfied
	void Survival_Selection();								//prepare new generation for the next iteration
};




/*DEFINITION OF struct Best_Individual: store the results of the best individual of each iteration*/
struct Genetic_Algorithm::Best_Individual									
{
	int (*solution)[SMAX+2];
	double fitness;
	int *bus_array;
	double *bus_fn;
	double transfer;
	double travel_time;
	double obj_value;
};	



/*DEFINITION OF GA(): the main body of the simulteneous genetic algorithm*/
void Genetic_Algorithm::GA()
{
	for (run = 0; run < RT; run ++)
	{
		cout << "Run time " << run << ": " << endl;
		/*1. call function Init_Population() to generate initial population*/
		cout << "1. create initial population" << endl;
		Init_Population(population, P);

		// for (int i = 0; i < P; i++)
		// {
		// 	cout << "Individual " << i+1 << ": " << endl;
		// 	for (int j = 0; j < R; j++)
		// 	{
		// 		cout << "Route " << j << ": " << " "; 
		// 		for (int k = 0; k < SMAX+2; k++)
		// 		{
		// 			cout << population[i][j][k] << " ";
		// 		}
		// 		cout << endl;
		// 	}
		// }

		/*2. call function Fitness() to calculate the fitness value for each solution of the initial population*/
		cout << "2. calculate fitness value for initial population" << endl;
		if (sim)
		{
			for (int i = 0; i < P; i++)
			{
				double obj_value1 = Obj_Value_Sim(population[i], bus_array1[i], bus_fn1[i]);
				Fit1[i] = Fitness(obj_value1);
			}
		}
		else
		{
			for (int i = 0; i < P; i++)
			{
				double obj_value1 = Obj_Value_Seq(population[i]);
				Fit1[i] = Fitness(obj_value1);
				Freq_Heuristic(population[i], bus_array1[i], bus_fn1[i]);
			}			
		}


		/*Iterate for #IT generations*/
		for (iter = 0; iter < IT; iter++)
		{
			/*3. call function Parents_Selection() to select #O parents for variation with Roulette Wheel Selection*/
			cout << "iteration " << iter << ": 3. roulette wheel selection" << endl;
			Parents_Selection();

			/*4. call function Offspring_Generate() to generate children with genetic operators*/
			cout << "iteration " << iter << ": 4. generate offspring with genetic operator" << endl;
			Offspring_Generate();

			/*5. call function Fitness() to calculate fitness value for each solution of the offspring population*/ 
			cout << "iteration " << iter << ": 5. calculate fitness value for offspring" << endl;
			if (sim)
			{
				for (int i = 0; i < O; i++)
				{
					double obj_value2 = Obj_Value_Sim(repair_population[i], bus_array2[i], bus_fn2[i]);
					Fit2[i] = Fitness(obj_value2);
				}
			}
			else
			{
				for (int i = 0; i < O; i++)
				{
					double obj_value2 = Obj_Value_Seq(repair_population[i]);
					Fit2[i] = Fitness(obj_value2);
					Freq_Heuristic(repair_population[i], bus_array2[i], bus_fn2[i]);
				}			
			}			

			/*6. call function Survival_Selection() to prepare new population for the next iteration*/
			cout << "iteration " << iter << ": 6. find the best solution and prepare new population for next iteration" << endl;
			Survival_Selection();
		}
	}
}



/*DEFINITION OF Init_Population(): generate the initial generation*/
void Genetic_Algorithm::Init_Population(int (*pop)[R][SMAX+2], int pop_size /*= P*/)
{
	int i, j, k;
	int start_node_id;
	int end_node_id;
	int inter_node_id;
	int start_node;
	int end_node;
	int inter_node;
	int pos;
	double inc_t;
	double travel_t;
	double tlt_des_t;
	double tsw_tlt_t;
	bool time_exceed = false;  //set time_exceed = false

	//srand((unsigned)time(NULL));

	/*generate a sequence of stops for each route in each solution of a population*/
	for (i = 0; i < pop_size; i++)
	{
		for (j = 0; j < R; j++)
		{
			/*1. randomly select a terminal node and a destination node for route i*/
			start_node_id = rand() % Y;
			end_node_id = rand() % V;
			start_node = Y_SET[start_node_id];
			end_node = V_SET[end_node_id];

			vector<int> temp_stop;
			temp_stop.push_back(start_node);
			temp_stop.push_back(end_node);

			travel_t = travel_time1[temp_stop[0]-1][temp_stop[1]-1];   //travel time of the current route
			tlt_des_t = travel_time2[end_node - 24];   //TLT-DES travel time
			tsw_tlt_t = travel_t - tlt_des_t;	 //TSW-TLT travel time
			time_exceed = (tsw_tlt_t > TMAX)? true : false;   //evaluate the value of time_exceed
			
			/*2. randomly select a node that is not appeared on the route*/
			while (temp_stop.size()-2 < SMAX && (!time_exceed))
			{
				inter_node_id = rand() % U;
				inter_node = U_SET[inter_node_id];
				while (find(temp_stop.begin(), temp_stop.end(), inter_node) != temp_stop.end())
				{
					inter_node_id = rand() % U;
					inter_node = U_SET[inter_node_id];
				}
				
				/*3. find a position to insert the node to minimize the total travel time*/
				double min_t = TMAX + max_t2;
				for (k = 1; k < temp_stop.size(); k++)
				{
					inc_t = (travel_time1[temp_stop[k-1]-1][inter_node-1] + travel_time1[inter_node-1][temp_stop[k]-1] - 
							 travel_time1[temp_stop[k-1]-1][temp_stop[k]-1]) + S;
					if (inc_t < min_t)
					{
						min_t = inc_t;
						pos = k;
					}
				}
				
				temp_stop.insert(temp_stop.begin()+pos, inter_node);
				tsw_tlt_t += inc_t;	 //TSW-TLT travel time
				time_exceed = (tsw_tlt_t > TMAX)? true : false;   //evaluate the value of time_exceed	
			}
			for (k = 0; k < temp_stop.size(); k++)
			{
				pop[i][j][k] = temp_stop[k]; //population[i][j][k]: the k-th gene (stop) of the j-th route of the i-th individual
			}
		}
	}
}



/*DEFINITION OF Trip_Time(): calculate the single-round trip time of a route*/
double Genetic_Algorithm::Trip_Time(int* input_route, int len)		//input is a route
{
	int i;
	double trip_t = 0;

	for (i = 0; i < len-1; i++)
	{
		if (input_route[i+1] != 0)
		{
			trip_t += travel_time1[input_route[i]-1][input_route[i+1]-1] + S;
		}
	}
	return trip_t;  //consider the stopping time at the interchange which is implicitly expressed in a route
}



/*DEFINITION OF Direct_Bus(): whether there is direct bus service from node n1 to node n2 in given 10 routes (1: no direct bus service)*/
bool Genetic_Algorithm::Direct_Bus(int n1, int n2, int (*ind)[SMAX+2])
{
	int i;
	bool direct = 1;

	for (i = 0; i < R; i++)
	{
		direct *= (1 - Pass2Element(ind[i], (SMAX+2), n1, n2));
	}
	return direct;
}



/*DEFINITION OF Travel_Tn(): calculate the travel time of route n from node n1 to node n2*/
double Genetic_Algorithm::Travel_Tn(int p, int des, int *route, int route_len)
{
	int i;
	int n1_index, n2_index;
	double travel_tn = 0;

	if (p == 29)
	{
		n2_index = des - 24;					
		travel_tn = travel_time2[n2_index];
	}
	else if (des == 29)
	{
		n1_index = find(route, route + route_len, p) - route;
		n2_index = FindLastId(route, route_len);
		int des_id = route[n2_index] - 24;		
		for (i = n1_index; i < n2_index; i++)
		{
			travel_tn += travel_time1[route[i]-1][route[i+1]-1] + S;
		}
		travel_tn = travel_tn - travel_time2[des_id] - S;
	}
	else
	{
		n1_index = find(route, route + route_len, p) - route;
		n2_index = find(route, route + route_len, des) - route;

		if (n1_index != n2_index)
		{
			for (i = n1_index; i < n2_index; i++)
			{
				travel_tn += travel_time1[route[i]-1][route[i+1]-1] + S;
			}

			if (n2_index != FindLastId(route, route_len))
			{
				travel_tn -= S;
			}
		}
		else
		{
			cerr << "Please input different nodes!" << endl;
		}
	}
	return travel_tn;
}



/*DEFINITION OF Avg_Travel_T: calculate the average travel time of passenges from node n1 to destination among all 10 routes*/
double Genetic_Algorithm::Avg_Travel_T(int p, int des, double* fn, int (*ind)[SMAX+2])  //n2 inputs a destination node
{
	int i;
	double avg_travel_t = 0;
	double fn_sum1 = 0, fn_sum2 = 0;
	double tn_sum1 = 0, tn_sum2 = 0;

	if (!(Direct_Bus(p, des, ind)))	//direct = 1: no direct bus service; direct = 0: direct bus service
	{
		for (i = 0; i < R; i++)  //for all routes
		{
			bool p_d = Pass2Element(ind[i], SMAX+2, p, des);
			fn_sum1 += fn[i] * p_d;
			if (p_d)
			{
				tn_sum1 += fn[i] * p_d * Travel_Tn(p, des, ind[i], SMAX+2);
			}
		}
		avg_travel_t = (tn_sum1 + 1) / fn_sum1;
	}
	else
	{
		for (i = 0; i < R; i++)  //for all routes
		{
			bool p_t = Pass2Element(ind[i], SMAX+2, p, 29);
			bool t_d = Pass2Element(ind[i], SMAX+2, 29, des);
			fn_sum1 += fn[i] * p_t;  //all node will pass through the interchange 
			fn_sum2 += fn[i] * t_d;
			if (p_t)
			{
				tn_sum1 += fn[i] * Travel_Tn(p, 29, ind[i], SMAX+2) * p_t;
			}
			if (t_d)
			{
				tn_sum2 += fn[i] * Travel_Tn(29, des, ind[i], SMAX+2) * t_d;
			}
		}
		if (!(fn_sum1 && fn_sum2))
		{
			avg_travel_t = Time_Penalty;
		}
		else
		{
			avg_travel_t = (tn_sum1 + 1) / fn_sum1 + (tn_sum2 + 1) / fn_sum2;		
		}
	}
	return avg_travel_t;
}



/*DEFINITION OF Total_Transfer(): calculate the total transfer times of a solution*/
double Genetic_Algorithm::Total_Transfer(int (*ind)[SMAX+2])
{
	int i, j;
	double total_transfer = 0;
	bool direct[U][V];

	/*Decide whether there's direct service between any demand points and destination*/
	for (i = 0; i < U; i++)
	{
		for (j = 0; j < V; j++)
		{
			direct[i][j] = Direct_Bus(i+1, j+24, ind);
			total_transfer += demands[i][j] * direct[i][j];
		}
	}

	return total_transfer;
}



/*DEFINITION OF Total_Travel_T(): calculate the original total travel time of a solution*/
double Genetic_Algorithm::Total_Travel_T(double* fn, int (*ind)[SMAX+2])
{
	int i, j;
	double total_travel_t = 0;

	for (i = 0; i < U; i++)
	{
		for (j = 0; j < V; j++)
		{
			total_travel_t += demands[i][j] * Avg_Travel_T(i+1, j+24, fn, ind);
		}
	}
	return total_travel_t;
}



/*DEFINITION OF FreqSetting(): a Frequencing Setting Initialization Heuristic*/
void Genetic_Algorithm::Freq_Setting(int p[10], int sum_w, int low=50, int up=100)  //prefer more even distribution
{
	int i;
	int sum_ori = 0;
	int sum_new = 0;
	int sum_diff = 0;

	//srand(time(0));		//set time seed

	for(i = 0; i < 10; i++)
	{
		p[i]= rand() % (up-low+1) + up;		//generate random numbers ranging from lower to upper bound 
		sum_ori += p[i];
	}

	for(i = 0; i < 10; i++)
	{
		p[i]= (p[i]/sum_ori) * sum_w;
		sum_new += p[i];
	}

	sum_diff = sum_w - sum_new;

	for(i = 0; i < sum_diff; i++)
	{
		p[rand() % 10] += 1;
	}
}



/*DEFINITION OF Freq_Heuristic(): Frequency Setting Heuristic to calculate the total improved travel time*/
double Genetic_Algorithm::Freq_Heuristic(int (*input_ind)[SMAX+2], int bus_array[10], double bus_fn[10])	//input_ind is an individual with 10 routes
{
	int i, j;
	double bus_fn_copy[10];
	double repair_bus_fn[10];
	double new_travel_time;
	double route_trip_time[10];

	/*1. Initialize a solution by randomly allocating buses to 10 routes*/
	Freq_Setting(bus_array, W);

	for (i = 0; i < R; i++)
	{
		route_trip_time[i] = Trip_Time(input_ind[i], SMAX+2);
		bus_fn[i] = CalFreq(bus_array[i], route_trip_time[i]);
	}

	double total_travel_time = Total_Travel_T(bus_fn, input_ind);

	copy(bus_fn, bus_fn+10, bus_fn_copy);

	for (i = 0; i < R-1; i++)
	{
		for (j = i+1; j < R; j++)
		{
			/*2. move 1 bus from route i to route j and evaluate the objective*/
			bus_fn_copy[i] = CalFreq(bus_array[i]-1, route_trip_time[i]);
			bus_fn_copy[j] = CalFreq(bus_array[j]+1, route_trip_time[j]);
			new_travel_time = Total_Travel_T(bus_fn_copy, input_ind);

			/*3. if objective value is improved, then continue*/
			if (new_travel_time < total_travel_time)
			{
				total_travel_time = new_travel_time;
				bus_array[i] -= 1;
				bus_array[j] += 1;
				bus_fn[i] = bus_fn_copy[i];
				bus_fn[j] = bus_fn_copy[j];
				continue;
			}
			/*4. if not, move 1 bus from route j to route i and evaluate the objective*/
			else
			{
				bus_fn_copy[i] = CalFreq(bus_array[i]+1, route_trip_time[i]);
				bus_fn_copy[j] = CalFreq(bus_array[j]-1, route_trip_time[j]);
				new_travel_time = Total_Travel_T(bus_fn_copy, input_ind);

				/*5. if objective value is improved, then continue*/
				if (new_travel_time < total_travel_time)
				{
					total_travel_time = new_travel_time;
					bus_array[i] += 1;
					bus_array[j] -= 1;
					bus_fn[i] = bus_fn_copy[i];
					bus_fn[j] = bus_fn_copy[j];
					continue;
				}
			}
		}
	}
	
	/*6. add penalty to the route whose frequency is smaller than the minimum frequency*/
	for (i = 0; i < R; i++)
	{
		repair_bus_fn[i] = bus_fn[i];
		if (bus_fn[i] < FMIN)
		{
			repair_bus_fn[i] = (double) (bus_fn[i] / Freq_Penalty);
		}
	}
	double repair_travel_time = Total_Travel_T(repair_bus_fn, input_ind);

	return repair_travel_time;
}



/*DEFINITION OF Obj_Value(): calculate the objective value of a given individual in simultaneous genetic algorithm*/
double Genetic_Algorithm::Obj_Value_Sim(int (*input_ind)[SMAX+2], int bus_array[10], double bus_fn[10])
{
	int i, j;
	double obj_v;
	double total_transfer = Total_Transfer(input_ind);
	double improved_travel_t = Freq_Heuristic(input_ind, bus_array, bus_fn);
	obj_v = B1 * total_transfer + B2 * improved_travel_t;

	return obj_v;
}



/*DEFINITION OF Obj_Value(): calculate the objective value of a given individual in sequential genetic algorithm*/
double Genetic_Algorithm::Obj_Value_Seq(int (*input_ind)[SMAX+2])
{
	int i, j;
	double obj_v;
	double total_transfer = Total_Transfer(input_ind);
	obj_v = total_transfer;

	return obj_v;
}


/*DEFINITION OF Fitness(): Frequency Setting Heuristic and Fitness Evaluation*/
double Genetic_Algorithm::Fitness(double obj_v)
{
	double fitness;

	if (obj_v > 0)
	{
		fitness = 1/obj_v;
	}
	else
	{
		fitness = 0;
	}
	return fitness;
}



/*DEFINITION OF Parents_Selection(): Roulette Wheel Selection approach to select parents for crossover*/
void Genetic_Algorithm::Parents_Selection()
{
	int i, j, k, h;
	double sum = 0;
	double SUM = 0;
	double r;
	
	/*1. calculate the sum of all fitness */
	for (i = 0; i < P; i++)
	{
		SUM += Fit1[i];
	}

	/*2. calculate the accumulated probability of the fitness value*/
	for (i = 0; i < P; i++)
	{
		sum += Fit1[i];
		pr[i] = sum / SUM;
	}
	
	/*3. select the crossover parents population by roulette wheel selection*/
	for (i = 0; i < P; i++)		
	{
		r = (double)rand() / RAND_MAX;

		for (int j = 0; j < P; j++)
		{
			if (r <= pr[j])		//if the random number is less than the accumulated probability of individual j
			{
				for (k = 0; k < R; k++)
				{
					for (h = 0; h < SMAX+2; h++)
					{
						roul_population[i][k][h] = population[j][k][h];
					}
				}
				break;
			}
		}
	}
}



/*DEFINITION OF Route_Crossover(): route crossover operator*/
void Genetic_Algorithm::Route_Crossover(int (*father)[SMAX+2], int (*mother)[SMAX+2], int (*route_child_ind1)[2*SMAX+2], int (*route_child_ind2)[2*SMAX+2])
{
	int i, j;
	int cut_point1 = 0, cut_point2 = 0;

	/*1. randomly generate two crossover points*/
	while (cut_point1 >= cut_point2)
	{
		cut_point1 = rand() % (R-1);
		cut_point2 = rand() % R;		//cut_point2 > cut_point1
	}

	/*2. get the middle crossover part of the two children*/
	for (i = cut_point1; i <= cut_point2; i++)
	{
		for (j = 0; j < SMAX+2; j++)
		{
			route_child_ind1[i][j] = father[i][j];	//ind1 takes father's middle part
			route_child_ind2[i][j] = mother[i][j];	//ind2 takes mother's middle part
		}
	}

	/*3. fill the front and end part of the two children*/
	if (cut_point1 == 0)
	{
		if (cut_point2 < R-1)
		{
			for (i = cut_point2+1; i < R; i++)
			{
				for (j = 0; j < SMAX+2; j++)
				{
					route_child_ind1[i][j] = mother[i][j];
					route_child_ind2[i][j] = father[i][j];
				}
			}
		}
	}
	else if (cut_point2 == R-1)
	{
		for (i = 0; i < cut_point1; i++)
		{
			for (j = 0; j < SMAX+2; j++)
			{
				route_child_ind1[i][j] = mother[i][j];
				route_child_ind2[i][j] = father[i][j];
			}
		}
	}
	else
	{
		for (i = 0; i < cut_point1; i++)
		{
			for (j = 0; j < SMAX+2; j++)
			{
				route_child_ind1[i][j] = mother[i][j];
				route_child_ind2[i][j] = father[i][j];
			}
		}
		for (i = cut_point2+1; i < R; i++)
		{
			for (j = 0; j < SMAX+2; j++)
			{
				route_child_ind1[i][j] = mother[i][j];
				route_child_ind2[i][j] = father[i][j];
			}
		}		
	}
}



/*DEFINITION OF Stop_Crossover(): stop crossover operator*/
void Genetic_Algorithm::Stop_Crossover(int (*father)[SMAX+2], int (*mother)[SMAX+2], int (*stop_child_ind1)[2*SMAX+2], int (*stop_child_ind2)[2*SMAX+2])
{
	int i, j;
	int cut_point_f1 = 1, cut_point_f2 = 0;
	int cut_point_m1 = 1, cut_point_m2 = 0;
	int rand_route_id1, rand_route_id2;
	int *stop_father, *stop_mother;
	int last_father_node_id, last_mother_node_id;
	int last_father_node, last_mother_node;
	vector<int> mother_candidates;
	vector<int> stop_mother_copy;
	vector<int> stop_father_copy;

	/*1. one route is randomly selected from the father*/
	while(mother_candidates.size() == 0)	//select a father which has the same destination in a mother
	{
		rand_route_id1 = rand() % R;
		stop_father = father[rand_route_id1];	//select a route as the father of the stop crossover
		last_father_node_id = FindLastId(stop_father, SMAX+2);
		last_father_node = stop_father[last_father_node_id];
		if (last_father_node_id <= 1)
			continue;

		/*2. another route with the same destination as the first route is randomly selected from mother*/
		for (i = 0; i < R; i++)
		{
			last_mother_node_id = FindLastId(mother[i], SMAX+1);
			last_mother_node = mother[i][last_mother_node_id];
			if ((last_mother_node == last_father_node) && (last_mother_node_id >= 2))
			{
				mother_candidates.push_back(i);
			}
		}
	}
	rand_route_id2 = mother_candidates[rand() % (mother_candidates.size())];	//randomly pick a route id from the mother candidates
	stop_mother = mother[rand_route_id2];	//select the route as the mother of the stop crossover
	last_mother_node_id = FindLastId(stop_mother, SMAX+1);
	last_mother_node = stop_mother[last_mother_node_id];

	/*copy the selected father and mother to vectors to conduct the stop crossover operator*/
	for (i = 0; i < last_father_node_id+1; i++)
	{
		stop_father_copy.push_back(stop_father[i]);
	}

	for (i = 0; i < last_mother_node_id+1; i++)
	{
		stop_mother_copy.push_back(stop_mother[i]);
	}

	/*3. randomly generate two crossover points in father and mother*/
	while (cut_point_f1 > cut_point_f2)
	{
		cut_point_f1 = rand() % (last_father_node_id-1) + 1;	//the starting terminal can't be involved in the stop crossover operation
		cut_point_f2 = rand() % (last_father_node_id-1) + 1;	
	}
	
	while (cut_point_m1 > cut_point_m2)
	{
		cut_point_m1 = rand() % (last_mother_node_id-1) + 1;
		cut_point_m2 = rand() % (last_mother_node_id-1) + 1;
	}

	/*4. the selected intermediate stops of the two routes are exchanged*/
	vector<int> sub_father {&stop_father_copy[cut_point_f1], &stop_father_copy[cut_point_f2+1]};
	stop_father_copy.erase(stop_father_copy.begin()+cut_point_f1, stop_father_copy.begin()+cut_point_f2+1);
	stop_father_copy.insert(stop_father_copy.begin()+cut_point_f1, stop_mother_copy.begin()+cut_point_m1, stop_mother_copy.begin()+cut_point_m2+1);
	
	stop_mother_copy.erase(stop_mother_copy.begin()+cut_point_m1, stop_mother_copy.begin()+cut_point_m2+1);
	stop_mother_copy.insert(stop_mother_copy.begin()+cut_point_m1, sub_father.begin(), sub_father.end());

	/*5. duplicated stops in the same route are eliminated*/
	DropDuplicate(stop_father_copy);
	DropDuplicate(stop_mother_copy);

	/*6. assign value to stop crossover children*/
	for (i = 0; i < R; i++)
	{
		if (i == rand_route_id1)
		{
			for (j = 0; j < stop_father_copy.size(); j++)
			{
				stop_child_ind1[i][j] = stop_father_copy[j];
			}			
		}
		else
		{
			for (j = 0; j < SMAX+2; j++)
			{
				stop_child_ind1[i][j] = father[i][j];
			}
		}
	}

	for (i = 0; i < R; i++)
	{
		if (i == rand_route_id2)
		{
			for (j = 0; j < stop_mother_copy.size(); j++)
			{
				stop_child_ind2[i][j] = stop_mother_copy[j];
			}			
		}
		else
		{
			for (j = 0; j < SMAX+2; j++)
			{
				stop_child_ind2[i][j] = mother[i][j];
			}
		}
	}
}



/*DEFINITION OF Insert_Mutation(): insert mutation operator*/
void Genetic_Algorithm::Insert_Mutation(int (*ind)[2*SMAX+2])
{
	int i;
	int rand_route;
	int rand_node;
	int route_len;
	int rand_pos;
	vector<int> route_vec;

	rand_route = rand() % R;
	route_len = FindLastId(ind[rand_route], 2*SMAX+2) + 1;
	
	/*insert random stop before random position in the random route*/
	for (i = 0; i < route_len; i++)
	{
		route_vec.push_back(ind[rand_route][i]);
		ind[rand_route][i] = 0;
	}

	rand_pos = rand() % route_len;		//rand_pos = rand() % (route_len-1) + 1;

	if (rand_pos == 0)
	{
		rand_node = Y_SET[rand() % Y];
	}
	else
	{
		rand_node = rand() % U + 1;
	}

	route_vec.insert(route_vec.begin() + rand_pos, rand_node);
	DropDuplicate(route_vec);	//a checking mechanism

	for (i = 0; i < route_vec.size(); i++)
	{
		ind[rand_route][i] = route_vec[i];
	}
}



/*DEFINITION OF Delete_Mutation(): delete mutation operator*/
void Genetic_Algorithm::Delete_Mutation(int (*ind)[2*SMAX+2])
{
	int i;
	int rand_route;
	int route_len;
	int rand_pos;
	vector<int> route_vec;

	rand_route = rand() % R;
	route_len = FindLastId(ind[rand_route], 2*SMAX+2) + 1;

	if (route_len >= 3)
	{
		for (i = 0; i < route_len; i++)
		{
			route_vec.push_back(ind[rand_route][i]);
			ind[rand_route][i] = 0;
		}

		rand_pos = rand() % (route_len-2) + 1;	//will not delete the starting node and the destination node
		route_vec.erase(route_vec.begin() + rand_pos);
		DropDuplicate(route_vec);	//a checking mechanism

		for (i = 0; i < route_vec.size(); i++)
		{
			ind[rand_route][i] = route_vec[i];
		}
	}
}


/*DEFINITION OF Swap_Mutation(): swap mutation operator*/
void Genetic_Algorithm::Swap_Mutation(int (*ind)[2*SMAX+2])
{
	int i;
	int rand_route1 = 0, rand_route2 = 0;
	int route_len1, route_len2;
	int rand_pos1, rand_pos2;
	vector<int> route_vec1; 
	vector<int> route_vec2;

	while (rand_route1 == rand_route2)
	{
		rand_route1 = rand() % R;
		rand_route2 = rand() % R;
	}
	
	route_len1 = FindLastId(ind[rand_route1], 2*SMAX+2) + 1;
	route_len2 = FindLastId(ind[rand_route2], 2*SMAX+2) + 1;

	/*two nodes in each route are of the same type*/
	if (route_len2 <= 2)
	{
		double prob = (double)rand() / RAND_MAX;
		if (prob <= 0.5)
		{
			rand_pos1 = 0;
			rand_pos2 = 0;
		}
		else
		{
			rand_pos1 = route_len1 - 1;
			rand_pos2 = route_len2 - 1;
		}
	}
	else if (route_len2 >= 3)
	{
		rand_pos1 = rand() % route_len1;
		if (rand_pos1 == 0)
		{
			rand_pos2 = 0;
		}
		else if (rand_pos1 == route_len1 - 1)
		{
			rand_pos2 = route_len2 - 1;
		}
		else
		{
			rand_pos2 = rand() % (route_len2 - 2) + 1;
		}
	}

	/*assign value to vectors to operate swap operation*/
	for (i = 0; i < route_len1; i++)
	{
		route_vec1.push_back(ind[rand_route1][i]);
		ind[rand_route1][i] = 0;
	}

	for (i = 0; i < route_len2; i++)
	{
		route_vec2.push_back(ind[rand_route2][i]);
		ind[rand_route2][i] = 0;
	}

	swap_ranges(route_vec1.begin()+rand_pos1, route_vec1.begin()+rand_pos1+1, route_vec2.begin()+rand_pos2);
	DropDuplicate(route_vec1);	//a checking mechanism
	DropDuplicate(route_vec2);	//a checking mechanism

	for (i = 0; i < route_vec1.size(); i++)
	{
		ind[rand_route1][i] = route_vec1[i];
	}

	for (i = 0; i < route_vec2.size(); i++)
	{
		ind[rand_route2][i] = route_vec2[i];
	}
}



/*DEFINITION OF Transfer_Mutation(): transfer mutation operator*/
void Genetic_Algorithm::Transfer_Mutation(int (*ind)[2*SMAX+2])
{
	int i;
	int rand_route1 = 0, rand_route2 = 0;
	int route_len1, route_len2;
	int rand_pos1, rand_pos2;
	vector<int> route_vec1; 
	vector<int> route_vec2;

	while (rand_route1 == rand_route2)
	{
		rand_route1 = rand() % R;
		rand_route2 = rand() % R;
	}

	route_len1 = FindLastId(ind[rand_route1], 2*SMAX+2) + 1;
	route_len2 = FindLastId(ind[rand_route2], 2*SMAX+2) + 1;

	if (route_len2 >= 3)
	{
		/*assign value to vectors to operate swap operation*/
		for (i = 0; i < route_len1; i++)
		{
			route_vec1.push_back(ind[rand_route1][i]);
			ind[rand_route1][i] = 0;
		}

		for (i = 0; i < route_len2; i++)
		{
			route_vec2.push_back(ind[rand_route2][i]);
			ind[rand_route2][i] = 0;
		}

		/*moves an intermediate stop on one route to another*/
		rand_pos1 = rand() % (route_len1-1) + 1;
		rand_pos2 = rand() % (route_len2-2) + 1;

		route_vec1.insert(route_vec1.begin() + rand_pos1, route_vec2[rand_pos2]);
		route_vec2.erase(route_vec2.begin() + rand_pos2);
		DropDuplicate(route_vec1);	//a checking mechanism
		DropDuplicate(route_vec2);	//a checking mechanism

		for (i = 0; i < route_vec1.size(); i++)
		{
			ind[rand_route1][i] = route_vec1[i];
		}

		for (i = 0; i < route_vec2.size(); i++)
		{
			ind[rand_route2][i] = route_vec2[i];
		}
	}
}



/*DEFINITION OF Crossover(): integrated crossover operator*/
void Genetic_Algorithm::Crossover()
{
	int i, j, k;
	int chromoN1 = 0, chromoN2 = 0;
	int Z1 = 0;
	int Z2 = 1;
	int father[R][SMAX+2] = {0};
	int mother[R][SMAX+2] = {0};

	/*1. initialize cross_population*/
	for (i = 0; i < O; i++)
	{
		for (j = 0; j < R; j++)
		{
			for (k = 0; k < 2*SMAX+2; k++)
			{
				cross_population[i][j][k] = 0;
			}
		}
	}

	for (i = 0; i < O/2; i++)		//crossover times: O/2
	{
		/*2. select the crossover parents*/
		chromoN1 = rand() % P;		//select #O crossover parents to generate #O offspring
		chromoN2 = rand() % P;		
		while (chromoN1 == chromoN2)
		{
			chromoN2 = rand() % P;
		}
		//cout << "chromoN1 " << chromoN1 << endl;
		//cout << "chromoN2 " << chromoN2 << endl;

		for (j = 0; j < R; j++)
		{
			for (k = 0; k < SMAX+2; k++)
			{
				father[j][k] = roul_population[chromoN1][j][k];
				mother[j][k] = roul_population[chromoN2][j][k];
			}
		}
		//cout << "crossover time: " << i << endl;
		
		/*3. conduct route crossover and stop crossover with probability*/
		double r = (double)rand() / RAND_MAX;
		//cout << "r " << r << endl; 

		if (r <= CROSS_RATE)
		{
			double sr = (double) rand() / RAND_MAX;
			//cout << "sr " << sr << endl;
			
			if (sr <= ROUTE_CPROB)
			{
				/*conduct route crossover*/
				//cout << "route crossover" << endl;
				Route_Crossover(father, mother, cross_population[Z1], cross_population[Z2]);
			}
			else if (sr <= STOP_CPROB + ROUTE_CPROB)
			{
				/*conduct stop crossover*/
				//cout << "stop_crossover" << endl;
				Stop_Crossover(father, mother, cross_population[Z1], cross_population[Z2]);
			}
		}
		else
		{
			/*do not crossover parents*/
			//cout << "not crossover" << endl;
			for (j = 0; j < R; j++)
			{
				for (k = 0; k < SMAX+2; k++)
				{
					cross_population[Z1][j][k] = father[j][k];
					cross_population[Z2][j][k] = mother[j][k];
				}
			}
		}
		Z1 += 2;
		Z2 += 2;
	}
	//cout << "crossover finish" << endl;
}



/*DEFINITION OF Mutation(): integrated mutation operator*/
void Genetic_Algorithm::Mutation()
{
	int i, j, k; 

	/*1. initialize muta_population as cross_population*/
	for (i = 0; i < O; i++)
	{
		for (j = 0; j < R; j++)
		{
			for (k = 0; k < 2*SMAX+2; k++)
			{
				muta_population[i][j][k] = cross_population[i][j][k];
			}
		}
	}

	for (i = 0; i < O; i++)
	{
		/*2. conduct mutation with probability*/
		double r = double(rand()) / RAND_MAX;
		//cout << "r " << r << endl;

		if (r <= MUTA_RATE)
		{
			double mr = double(rand()) / RAND_MAX;
			//cout << "mr " << mr << endl;
			if (mr <= INSERT_MPROB)
			{
				//cout << "insert mutation" << endl;
				Insert_Mutation(muta_population[i]);
			}
			else if (mr <= INSERT_MPROB + DELETE_MPROB)
			{
				//cout << "delete mutation" << endl;
				Delete_Mutation(muta_population[i]);
			}
			else if (mr <= INSERT_MPROB + DELETE_MPROB + SWAP_MPROB)
			{
				//cout << "swap mutation" << endl;
				Swap_Mutation(muta_population[i]);
			}
			else if (mr <= INSERT_MPROB + DELETE_MPROB + SWAP_MPROB + TRANSFER_MPROB)
			{
				//cout << "transfer mutation" << endl;
				Transfer_Mutation(muta_population[i]);
			}
		}
	}
}



/*DEFINITION OF Stop_Seq_Improve(): stop sequence improvement heuristic*/
void Genetic_Algorithm::Stop_Seq_Improve(int *route, int route_len)
{
	double trip_t = Trip_Time(route, route_len);
	int stop_len = FindLastId(route, route_len) + 1;

	if (stop_len >= 4)
	{
		for (int i = 1; i < stop_len-2; i++)
		{
			for (int j = i+1; j < stop_len-1; j++)
			{
				/*exchange i and j*/
				iter_swap(route + i, route + j);		

				double tmp_trip_t = Trip_Time(route, route_len);
				if (tmp_trip_t < trip_t)
				{
					trip_t = tmp_trip_t;	//if the trip time is reduced, continue the exchange
				}
				else			//otherwise, undo the exchange
				{
					iter_swap(route + i, route + j);	
				}
			}
		}		
	}
}



/*DEFINITION OF Repair(): solution Repairing to satisfy the constraints of intermediate stops and travel time*/
void Genetic_Algorithm::Repair(int *route, int route_len)	//1. maximum number of intermediate stops and 2. maximum allowable travel time
{
	int i, j, k;
	int stop_len = FindLastId(route, route_len) + 1;
	double tsw_tlt_t = Travel_Tn(route[0], 29, route, route_len);
	vector<int> route_vec;
	double dec_t;
	int pos;

	for (k = 0; k < stop_len; k++)
	{
		route_vec.push_back(route[k]);
		route[k] = 0;
	}

	while ((route_vec.size() > SMAX + 2) || (tsw_tlt_t > TMAX))	
	{
		double max_t = 0;
		for (i = 1; i < route_vec.size()-1; i++)
		{
			dec_t = (travel_time1[route_vec[i-1]-1][route_vec[i]-1] + travel_time1[route_vec[i]-1][route_vec[i+1]-1] - 
					 travel_time1[route_vec[i-1]-1][route_vec[i+1]-1]) + S;
			if (dec_t > max_t)
			{
				max_t = dec_t;
				pos = i;
			}
		}
		route_vec.erase(route_vec.begin() + pos);
		tsw_tlt_t -= max_t;
	}

	for (i = 0; i < route_vec.size(); i++)
	{
		route[i] = route_vec[i];
	}
}



/*DEFINITION OF Offspring_Generate(): generate offspring through genetic operator*/
void Genetic_Algorithm::Offspring_Generate()
{
	int i, j, k;

	Crossover();
	//cout << "crossover over" << endl;

	Mutation();
	//cout << "mutation over" << endl;

	// for (i = 0; i < O; i++)
	// {
	// 	cout << "individual " << i << endl;
	// 	for (j = 0; j < R; j++)
	// 	{
	// 		cout << "route " << j << ": ";
	// 		for (k = 0; k < 2*SMAX+2; k++)
	// 		{
	// 			cout << muta_population[i][j][k] << " ";
	// 		}
	// 		cout << endl;
	// 	}
	// }
	
	for (i = 0; i < O; i++)
	{
		for (j = 0; j < R; j++)
		{
			//cout << "Stop Seq Improve: " << i << "-" << j << endl;
			Stop_Seq_Improve(muta_population[i][j], 2*SMAX+2);
			
			// for (k = 0; k < 2*SMAX+2; k++)
			// {
			// 	cout << muta_population[i][j][k] << " ";
			// }
			// cout << endl;

			//cout << "Repair: " << i << "-" << j << endl;
			Repair(muta_population[i][j], 2*SMAX+2);
			//int route[2*SMAX+2] = {9, 10, 15, 14, 19, 23, 18, 17, 8, 6, 3, 12, 21, 26, 0, 0, 0, 0};
			//Repair(route, 2*SMAX+2);
						
			for (k = 0; k < SMAX+2; k++)
			{
				repair_population[i][j][k] = muta_population[i][j][k];
			}
		}
	}
}



/*DEFINITION OF Diversity_Control(): calculate new hamming distance between an individual and the best individual*/
double Genetic_Algorithm::Diversity_Control(int (*ind1)[SMAX+2], int (*ind2)[SMAX+2], double param_c /*= 0.08*/, double param_a /*= 0.002*/)
{
	int i;
	int stop_len1;
	int stop_len2;
	int hamming_dist = 0;
	int total_pairs = 0;
	double select_prob;
	
	for (i = 0; i < R; i++)
	{
		stop_len1 = FindLastId(ind1[i], SMAX+2) + 1;
		stop_len2 = FindLastId(ind2[i], SMAX+2) + 1;
		hamming_dist += DiffPairs(ind1[i], ind2[i], stop_len1, stop_len2);
		total_pairs += stop_len1 + stop_len2 - 2;
	}

	select_prob = pow(((1-param_c) * (double(hamming_dist)/total_pairs) + param_c), param_a);

	return select_prob;
}



/*DEFINITION OF Diff_Individuals(): find whether two individuals are the same*/
bool Genetic_Algorithm::Diff_Individuals(int (*ind1)[SMAX+2], int (*ind2)[SMAX+2])
{
	int i, j;
	for (i = 0; i < R; i++)
	{
		for (j = 0; j < SMAX+2; j++)
		{
			if (ind1[i][j] != ind2[i][j])
			{
				return true;
			}
		}
	}
	return false;
}



void Genetic_Algorithm::Mix_Population(int mix_population[P+O][R][SMAX+2], double mix_fit[P+O], int mix_bus_array[P+O][R], 
									   double mix_bus_fn[P+O][R], double sort_mix_fit[P+O], int sort_fit_preindex[P+O])
{
	int i, j, k;
	
	/*1. mix population and repair population*/
	for (i = 0; i < (P + O); i++)
	{
		if (i < P)
		{
			for (j = 0; j < R; j++)
			{
				for (k = 0; k < SMAX+2; k++)
				{
					mix_population[i][j][k] = population[i][j][k];
				}
				mix_bus_array[i][j] = bus_array1[i][j];
				mix_bus_fn[i][j] = bus_fn1[i][j];
			}
			mix_fit[i] = Fit1[i];
		}
		else
		{
			for (j = 0; j < R; j++)
			{
				for (k = 0; k < SMAX+2; k++)
				{
					mix_population[i][j][k] = repair_population[i - P][j][k];
				}
				mix_bus_array[i][j] = bus_array2[i - P][j];
				mix_bus_fn[i][j] = bus_fn2[i - P][j];
			}
			mix_fit[i] = Fit2[i - P];
		}
	}

	/*2. sort mix fitness from large to small*/	
	for (i = 0; i < P+O; i++)
	{
		sort_mix_fit[i] = mix_fit[i];
	}
	
	sort(sort_mix_fit, sort_mix_fit + (P + O), greater<double>());

	for (i = 0; i < P+O; i++)
	{
		sort_fit_preindex[i] = (find(mix_fit, mix_fit + (P + O), sort_mix_fit[i]) - mix_fit);
	}
}



/*DEFINITION OF Freq_Adjust(): adjust bus allocation to satisfy the minimum frequency constraint*/
void Genetic_Algorithm::Freq_Adjust(int bus_arr[R], double bus_fn[R])
{
	int i, j;
	double round_trip_tn[R];

	for (i = 0; i < R; i++)
	{
		round_trip_tn[i] = double(bus_arr[i]) / bus_fn[i];
	}

	for (i = 0; i < R; i++)
	{
		while (bus_fn[i] < FMIN)
		{
			int max_id = max_element(bus_fn, bus_fn + R) - bus_fn;
			if (bus_fn[max_id] < FMIN)
			{
				break;
			}
            bus_arr[max_id] -= 1;
			bus_arr[i] += 1;
			bus_fn[max_id] = double(bus_arr[max_id]) / round_trip_tn[max_id];
			bus_fn[i] = double(bus_arr[i]) / round_trip_tn[i];		
		}
	}
}



/*DEFINITION OF Best_Solution(): report the best solution of the iteration (including the parents and the children)*/
void Genetic_Algorithm::Best_Solution(int best_ind[R][SMAX+2], double best_fit, int bus_arr[R], double bus_fn[R])
{
	int i, j, k;
	
	double transfer = Total_Transfer(best_ind);
	Freq_Adjust(bus_arr, bus_fn);
	double total_travel_time = Total_Travel_T(bus_fn, best_ind);
	double obj_v = B1 * transfer + B2 * total_travel_time;

	struct Best_Individual best = {best_ind, best_fit, bus_arr, bus_fn, transfer, total_travel_time, obj_v};
	if (iter == IT-1)
	{
		best_inds.push_back(best);
	}

	cout << "----------------------------" << endl;
	cout << "Best Solution of iteration " << iter << ": " << endl;	//report the best solution
	for (j = 0; j < R; j++)
	{
		cout << "Route " << j << ": ";
		int best_stop_len = FindLastId(best_ind[j], SMAX+2) + 1;
		for (k = 0; k < best_stop_len-1; k++)
		{
			cout << best_ind[j][k] << "->";
		}
		cout << best_ind[j][best_stop_len-1] << endl;
	}

	cout << "Bus Allocation of the Best Solution: " << endl;

	for (j = 0; j < R; j++)
	{
		cout << bus_arr[j] << " ";
	}
	cout << endl;
	
	cout << "Route Frequency (bus/min) of the Best Solution: " << endl;
	for (j = 0; j < R; j++)
	{
		cout << bus_fn[j] << " ";
	}
	cout << endl;

	cout << "Fitness value is " << best_fit << endl;
	cout << "Total Transfer is " << transfer << endl;
	cout << "Total Travel Time is (min)" << total_travel_time << endl;
	cout << "Objective Value is " << obj_v << endl;
	cout << "----------------------------" << endl;
}



/*DEFINITION OF Survival_Selection(): prepare new generation for the initial population of the  next iteration*/
void Genetic_Algorithm::Survival_Selection()
{
	int i, j, k, h;
	int mix_population[P+O][R][SMAX+2] = {0};
	double mix_fit[P+O] = {0};
	int mix_bus_array[P+O][R] = {0};
	double mix_bus_fn[P+O][R] = {0};
	double sort_mix_fit[P+O] = {0};
	int sort_fit_preindex[P+O] = {0};
	bool diff_ind[P+O];
	double diversity_pr[P+O];
	double survival_fit[P] = {0};
	bool diff;
	double pr;
	int sort_ind_id;
	int select_id;
	int survival_id;
	int gap_num;

	Mix_Population(mix_population, mix_fit, mix_bus_array, mix_bus_fn, sort_mix_fit, sort_fit_preindex);
	Best_Solution(mix_population[sort_fit_preindex[0]], sort_mix_fit[0], mix_bus_array[sort_fit_preindex[0]], mix_bus_fn[sort_fit_preindex[0]]);

	/*1. select the best individual in M for survival into the next generation*/
	for (j = 0; j < R; j++)
	{
		for (k = 0; k < SMAX+2; k++)
		{
			survival_population[0][j][k] = mix_population[sort_fit_preindex[0]][j][k];
		}
	}
	survival_fit[0] = sort_mix_fit[0];
	
	/*1. eliminate the duplicate individuals*/
	diff_ind[0] = true;
	for (i = 1 ; i < P+O; i++)
	{
		diff = Diff_Individuals(mix_population[sort_fit_preindex[i-1]], mix_population[sort_fit_preindex[i]]);
		if (diff)
		{
			diff_ind[i] = true;
		}
		else
		{
			diff_ind[i] = false;
		}
	}

	/*2. calculate the selection probability by diversity control */
	for (i = 0; i < P+O; i++)
	{
		diversity_pr[i] = Diversity_Control(mix_population[sort_fit_preindex[i]], survival_population[0]);
	}
	
	/*3. Select individuals in the mix population in the order of their fitness values based on the diversity control probability*/
	survival_id = 1;	//control survival population's index	
	select_id = 1;		//control selection from mix population
	while ((survival_id < P) && (select_id < P+O))	//generate P survival population
	{
		sort_ind_id = sort_fit_preindex[select_id];
		if (diff_ind[select_id])
		{
			if (diversity)
			{
				pr = double(rand()) / RAND_MAX;
				if (pr <= diversity_pr[select_id])
				{
					for (k = 0; k < R; k++)
					{
						for (h = 0; h < SMAX+2; h++)
						{
							survival_population[survival_id][k][h] = mix_population[sort_ind_id][k][h];
						}
					}
					survival_fit[survival_id] = sort_mix_fit[select_id];
					survival_id++;
				}
			}
			else if (!(diversity))
			{
				for (k = 0; k < R; k++)
				{
					for (h = 0; h < SMAX+2; h++)
					{
						survival_population[survival_id][k][h] = mix_population[sort_ind_id][k][h];
					}
				}
				survival_fit[survival_id] = sort_mix_fit[select_id];	
				survival_id++;			
			}
		}
		select_id++;
	}

	/*4. new individuals randomly generated to fill the gap between the population size and the survived individual numbers*/
	gap_num = P - survival_id;
	if (gap_num > 0)
	{
		int gap_population[gap_num][R][SMAX+2] = {0};
		int gap_bus_array[gap_num][R] = {0};
		double gap_bus_fn[gap_num][R] = {0};

		Init_Population(gap_population, gap_num);
		
		for (i = 0; i < gap_num; i++)
		{
			for (j = 0; j < R; j++)
			{
				for (k = 0; k < SMAX+2; k++)
				{
					survival_population[i + survival_id][j][k] = gap_population[i][j][k];
				}
			}
			
			double obj_value;
			if (sim)
			{
				obj_value = Obj_Value_Sim(gap_population[i], gap_bus_array[i], gap_bus_fn[i]);
			}
			else
			{
				obj_value = Obj_Value_Seq(gap_population[i]);
			}

			survival_fit[i + survival_id] = Fitness(obj_value);
		}
	}

	/*5. assign mix population to population and calculate the new fitness value*/
	for (i = 0; i < P; i++)
	{
		for (j = 0; j < R; j++)
		{
			for (k = 0; k < SMAX+2; k++)
			{
				population[i][j][k] = survival_population[i][j][k];
			}
		}
		Fit1[i] = survival_fit[i];
	}
	
	/*9. report the solutions found in this iteration*/
	// cout << "Survived Population:" << endl;
	// for (i = 0; i < P; i++)
	// {
	// 	cout << "Individual" << i << ": " << endl;
	// 	for (j = 0; j < R; j++)
	// 	{
	// 		cout << "Route " << j << ": ";
	// 		for (k = 0; k < SMAX+2; k++)
	// 		{
	// 			cout << population[i][j][k] << " ";
	// 		}
	// 		cout << endl;
	// 	}
	// 	cout << endl;
	// }
}



/*DEFINITION OF Set_Cross_Rate(): set the crossover trigger rate for each individual*/
void Genetic_Algorithm::Set_Cross_Rate(double rate)
{
	CROSS_RATE = rate;
}



/*DEFINITION OF Set_Muta_Rate(): set the mutation trigger rate for each individual*/
void Genetic_Algorithm::Set_Muta_Rate(double rate)
{
	MUTA_RATE = rate;
}



/*DEFINITION OF Set_Cross_Combine_Rate(): set the combinative rate of crossover methods*/
void Genetic_Algorithm::Set_Cross_Combine_Rate(double route_rate, double stop_rate)
{
	if (route_rate + stop_rate == 1)
	{
		ROUTE_CPROB = route_rate;
		STOP_CPROB = stop_rate;
	}
	else
		cerr << "Rates should be sumed to 1!";
}



/*DEFINITION OF Set_Muta_Combine_Rate(): set the combinative rate of mutation methods*/
void Genetic_Algorithm::Set_Muta_Combine_Rate(double insert_rate, double delete_rate, double swap_rate, double transfer_rate)
{
	if (insert_rate + delete_rate + swap_rate + transfer_rate == 1)
	{
		INSERT_MPROB = insert_rate;
		DELETE_MPROB = delete_rate;
		SWAP_MPROB = swap_rate;
		TRANSFER_MPROB = transfer_rate;
	}
	else
		cerr << "Rates should be sumed to 1!";
}



/*DEFINITION OF Set_Gens_Num(): set the number of iterations/generations for each run time*/
void Genetic_Algorithm::Set_Gens_Num(int generations)
{
	IT = generations;
}



/*DEFINITION OF Set_Run_Times(): set algorithm run times*/
void Genetic_Algorithm::Set_Run_Times(int run_times)
{
	RT = run_times;
}



/*DEFINITION OF Set_B1_Weight(): set B1 weight*/
void Genetic_Algorithm::Set_B1_Weight(int weight)
{
	if ((B1 >= 1) && (B1 <= 200))
		B1 = weight;
	else
		cerr << "B1 should vary from 1 to 200!";
}



/*DEFINITION OF Set_Freq_Penalty(): set the frequency penalty value*/
void Genetic_Algorithm::Set_Freq_Penalty(double penl)
{
	Freq_Penalty = penl;
}



/*DEFINITION OF Set_Time_Penalty(): set the time penalty value*/
void Genetic_Algorithm::Set_Time_Penalty(int penl)
{
	Time_Penalty = penl;
}



/*DEFINITION OF Set_Diversity(): decide whether to use diversity control mechanism*/
void Genetic_Algorithm::Set_Diversity(bool div)
{
	diversity = div;
}						




/*DEFINITION OF Set_Simultaneous(): decide whether to use simultaneous or sequential algorithm*/
void Genetic_Algorithm::Set_Simultaneous(bool simul)
{
	sim = simul;
}			