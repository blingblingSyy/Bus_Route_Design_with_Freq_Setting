# Project Title
- The code reproduces the paper: A Simultaneous Bus Route Design and Frequency Setting Problem for Tin Shui Wai, Hong Kong (author: W.Y.Szeto, Yongzhong Wu).
- Copyright by Yingying Sun

## Project Description
- Introduced by the paper, the codes design a specific GA to simultaneously implement the route design problem and the frequency setting problem. The code also compares the sequential GA method with the simultaneous method. The difference of the two methods is illustrated with a student's t test.

- The result with 10 runs, each of which contains 500 generations, respectively for the simultaneous GA and the sequential GA, is involved in the "result.dat" file.

### Solution Representation
- Different from traditional route design problem, the problem takes into account the combination of routes and the stops sequence at the same time. Therefore, different genetic operators are designed to pertube the routes combination and sequences of stops at the same time.

- All routes will pass through the same interchange stop (TLT) before arrive at the destination. The interchange stop is implicitly presented in routes.

- As the crossover operators may lead to unequal length of stops, an multi-dimensional array with maximum possible stop length is introduced to contain a specific route.

### Fitness Evaluation
- The fitness values are measured differently between the simultaneous GA and the sequential GA. Therefore, it is not allowed to compare the two methods' fitness values. Instead, we compare the objective value of the two methods.

- The simultaneous GA considers minimizing the transfer times and the total travel time, in order to improve the route design problem and the frequency setting problem at the same time.

- The sequential GA first considers minimizing the transfer times to design the route, only after which the frequency setting heuristic is implemented once to optimize and obtain the total travel time. 

### Genetic Operators
- Two crossover operators are designed, which are route crossover operator and stop crossover operator.

- Four mutation operators are designed, which are insert muation operator, delete mutation operator, swap mutation operator and transfer mutation operator. Different from the paper, the insert mutation operator allows inserting a bus terminal stop at the front of a route, in order to increase the diversity of the starting nodes of each route.

### Penalty Parameters
- Frequency Parameter: reduce the frequency by a non-negative parameter to gradually throw away solutions with infeasible frequency. Set with trial and error. 

- Time Parameter: some solutions may leave some nodes unvisited all the time. If this occurs, add a large time penalty to the average travel time to reduce the fitness values. Set with trial and error.

- Both penalties are applied in simultaneous GA. Since sequencial GA does not consider total travel time in the target, the penalties do not work in sequential GA.

### Frequency Adjustment Heuristic
- The code improves the paper by a Frequency Adjustment Heuristic algorithm.

- Although a frequency penalty is applied, there still exists possibility that the frequency in the best solution cannot reach the minimum requirment, or the variance of frequency between each route is large. Therefore, the Frequency Adjustment Heuristic is further applied to revise the frequency setting of the best solution.

## Getting Started
- The code must be compiled as a C++11 application.