#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <vector>
#include <climits>
#include <ctime>
#include "nr3.h"
#include "ludcmp.h"

#define NUM_STATES 65
#define NUM_ACTIONS 4
#define PRINT_UTILITY_PRECISION 2

#define N 0
#define E 1
#define S 2
#define W 3


/*

TODO

epsilon print out in policy iteration?

fix return values and interactions with global arrays -- have things print and construct new utilities/policies
	IN each iteration function, rather than working with globals at all (needs work on init, then -- send a pointer)

fix print out format for correct double precision (2, currently 3?)
fix prototype orders and parameters based on above

add matrix library (armadillo library) and correct construction of new utilities in policy iteration

*/

double T[NUM_STATES][NUM_ACTIONS][NUM_STATES];
double R[NUM_STATES];

enum Iter {Value, Policy};
static const std::string iterStrings[] = {"Value Iteration", "Policy Iteration"};

void initMDP(double, double, double, double);
void valueIteration(double, double, double, double, double);
void policyIteration(double, double, double, double, double);

std::string action(int);
bool extractPolicy(std::vector<double> &, std::vector<int> &);

void printResults(double, int, Iter, double, double, double, double, double, std::vector<double>, std::vector<int>);

int main(int argc, char* argv[])
{

	// help/cmd errors
	if (argc != 8) {
		std::cout << "Incorrect Parameters. Exiting." << std::endl
					<< "1: Discount Rate (Double)" << std::endl
					<< "2: Max error state (Double)" << std::endl
					<< "3: Key loss probability (Double)" << std::endl
					<< "4: Positive terminal reward (Double)" << std::endl
					<< "5: Negative terminal reward (Double)" << std::endl
					<< "6: Step cost (double)" << std::endl
					<< "7: Iteration type? (v or p)" << std::endl;
		return -1;
	}

	// reading in args for problem parameters
	double discount = atof(argv[1]);
	double epsilon = atof(argv[2]);
	double keyLoss = atof(argv[3]);
	double posTerminal = atof(argv[4]);
	double negTerminal = atof(argv[5]);
	double stepCost = atof(argv[6]);

	std::string iterArg(argv[7]);

	Iter iter;
	if (iterArg.compare("v") == 0) 
		iter = Iter::Value;
	else if (iterArg.compare("p") == 0)
		iter = Iter::Policy;
	else {
		std::cout << "Incorrect iteration type specified (not v or p)" << std::endl;
		return 1;
	}

	initMDP(negTerminal, posTerminal, stepCost, keyLoss);

	if (iter == Iter::Value)
		valueIteration(discount, epsilon, posTerminal, negTerminal, stepCost);
	else if (iter == Iter::Policy)
		policyIteration(discount, epsilon, posTerminal, negTerminal, stepCost);

	return 0;
}

// synchronous, in-place
void valueIteration(double discount, double epsilon, double posTerminal, 
					double negTerminal, double stepCost)
{

	clock_t start = clock();
	int numIter = 0;

	std::vector<double> utility;
	std::vector<int> policy;

	for (int s = 0; s < NUM_STATES; s++) {
		utility.push_back(0);
		policy.push_back(N);
	}
	double delta = 0;

	do {
		delta = 0;
		for (int s = 0; s < NUM_STATES; s++) {

			double uPS = 0;
			double aMaxVal = INT_MIN;

			for (int a = 0; a < NUM_ACTIONS; a++) {

				double currSum = 0;
				for (int sP = 0; sP < NUM_STATES; sP++) {
					currSum += T[s][a][sP] * utility[sP];
				}

				if (currSum > aMaxVal) {
					policy[s] = a;
					aMaxVal = currSum;
				}
			}

			uPS = R[s] + discount * aMaxVal;

			if (std::abs(uPS - utility[s]) > delta)
				delta = std::abs(uPS - utility[s]);

			utility[s] = uPS;
		}

		numIter++;

	} while (delta >= epsilon * (1-discount) / discount);

	clock_t end = clock();
	double solTime = (double)(end - start) / CLOCKS_PER_SEC;

	printResults(solTime, numIter, Iter::Value, stepCost,
				 discount, epsilon, posTerminal, negTerminal,
				 utility, policy);

	return;
}

void policyIteration(double discount, double epsilon, double posTerminal, 
					double negTerminal, double stepCost)
{

	clock_t start = clock();
	int numIter = 0;

	//Initialize policy and utility vectors
	//In order to solve the systems of equations, we use the Numerical Recipes package
	//which asks for a matrix and two vectors of "Doubs", a data structure unique to Numerical Recipes
	std::vector<int> policy;
	VecDoub utilityVD(NUM_STATES);

	// set inital random (all north)
	for (int s = 0; s < NUM_STATES; s++) {
		policy.push_back(N);
		utilityVD[s] = N;
	}

	std::vector<double> utility;
	for(int s=0; s < NUM_STATES; s++){
		utility.push_back(s);
	}


	for(int s=0; s < NUM_STATES; s++){
		utility[s] = utilityVD[s];
	}


	bool policyChange = true;

	while(policyChange) {
		policyChange = false;

		
		//Initialize matrix of Left side of the system of equations
		MatDoub allLCoeffs(NUM_STATES,NUM_STATES);

		//Initialize the right side of the equation (b part of Ax=b)
		VecDoub currRCoeffs(NUM_STATES);

		// iterate through every state
		for (int s = 0; s < NUM_STATES; s++) {

			currRCoeffs[s] = R[s];

			std::vector<double> currLCoeffs;
			for (int sP = 0; sP < NUM_STATES; sP++) {
				allLCoeffs[s][sP] = -1*T[s][policy[s]][sP] * discount;

				if(s == sP){
					allLCoeffs[s][sP] ++;
				}
			}
			
		}

		//Performs Lower Upper Decomposition to simplify matrix
		//(Factors the matrix as a product of lower and upper triangular matrices)
		LUdcmp alu(allLCoeffs);

		//Solves the system of equations to find the utility of each state
		alu.solve(currRCoeffs,utilityVD);

		//Give utility the values of utilityVD
		for(int s=0; s < NUM_STATES; s++){
			utility[s] = utilityVD[s];
		}


		// UTILITY + SOLUTION TO MATRIX CALCULATION W/ COEFFS
		policyChange = extractPolicy(utility, policy);
		numIter++;
	}

	clock_t end = clock();
	double solTime = (double)(end - start)/ CLOCKS_PER_SEC;

	printResults(solTime, numIter,  Iter::Policy, stepCost,
			 	 discount, epsilon, posTerminal, negTerminal,
				 utility, policy);

}

bool extractPolicy(std::vector<double> &utilities, std::vector<int> &policy)
{

	bool policyChange = false;

	std::vector<int> newPolicy;
	// go through every state and fight the max weighted average of utilities given a move
	for (int s = 0; s < NUM_STATES; s++) {
		int newAction = N;
		double aMaxVal = INT_MIN;

		for (int a = 0; a < NUM_ACTIONS; a++) {
			double aCurrVal = 0;

			for (int sP = 0; sP < NUM_STATES; sP++) {
				aCurrVal += T[s][a][sP] * utilities[sP];
			}

			if (aCurrVal > aMaxVal) {
				newAction = a;
				aMaxVal = aCurrVal;
			}
		}

		if (newAction != policy[s])
			policyChange = true;

		policy[s] = newAction;
	}

	return policyChange;
}

void printResults(double solTime, int numIter, Iter iter, double stepCost,
				  double discount, double epsilon, double posTerminal, double negTerminal,
				  std::vector<double> utility, std::vector<int> policy)
{
	
	// formateString is a C-style format string to use with Java's printf-wannabe
	// method; the format string specifies what the output should look like, including
	// format specifications for values, and the actual items to be printed are
	// the arguments to printf that come after the format string.  in the following,
	// if PRINT_UTILITY_PRECISION is 2, the format string would be:
	//
	//    "%s%2d%s %.2f %s    "
	//
	// This means that the output will be:
	//    a string,
	//    an integer that should be printed in 2 spaces, 
	//    a string,
	//    a space (spaces in the format string are printed literally),
	//    a floating-point number printed in 5 spaces with 
	//          PRINT_UTILITY_PRECISION digits after the decimal point, 
	//    a space,
	//    a string, and
	//    4 spaces.
	//
	// The arguments that come after specify *what* string, *what* integer, etc.

	std::cout << std::endl << std::endl << std::fixed << std::setprecision(PRINT_UTILITY_PRECISION);

	for (int s = 58 ; s <= 64 ; s += 2)
		std::cout << "(" << std::setw(2) << s << std::setw(1) << ") " << std::setw(5) << utility[s] << " (" << action(policy[s]) << ")    "; 
	std::cout << std::endl;
	
	for (int s = 59 ; s <= 63 ; s += 2)
		std::cout << "(" << std::setw(2) << s << std::setw(1) << ") " << std::setw(5) << utility[s] << " (" << action(policy[s]) << ")    "; 
	std::cout << std::endl << std::endl;
	
	
	for (int s = 50 ; s <= 56 ; s += 2)
		std::cout << "(" << std::setw(2) << s << std::setw(1) << ") " << std::setw(5) << utility[s] << " (" << action(policy[s]) << ")    "; 
	std::cout << std::endl;
	
	for (int s = 51 ; s <= 57 ; s += 2)
		std::cout << "(" << std::setw(2) << s << std::setw(1) << ") " << std::setw(5) << utility[s] << " (" << action(policy[s]) << ")    "; 
	std::cout << std::endl << std::endl;

	for (int s = 40 ; s <= 48 ; s += 2)
		std::cout << "(" << std::setw(2) << s << std::setw(1) << ") " << std::setw(5) << utility[s] << " (" << action(policy[s]) << ")    "; 
	std::cout << std::endl;

	for (int s = 41 ; s <= 49 ; s += 2)
		std::cout << "(" << std::setw(2) << s << std::setw(1) << ") " << std::setw(5) << utility[s] << " (" << action(policy[s]) << ")    "; 
	std::cout << std::endl << std::endl;	
	
	for (int s = 30 ; s <= 38 ; s += 2)
		std::cout << "(" << std::setw(2) << s << std::setw(1) << ") " << std::setw(5) << utility[s] << " (" << action(policy[s]) << ")    "; 
	std::cout << std::endl;
	
	for (int s = 31 ; s <= 39 ; s += 2)
		std::cout << "(" << std::setw(2) << s << std::setw(1) << ") " << std::setw(5) << utility[s] << " (" << action(policy[s]) << ")    "; 
	std::cout << std::endl << std::endl;	
	
	for (int s = 0 ; s <= 14 ; s += 2)
		std::cout << "(" << std::setw(2) << s << std::setw(1) << ") " << std::setw(5) << utility[s] << " (" << action(policy[s]) << ")    "; 
	std::cout << std::endl;
	
	for (int s = 1 ; s <= 15 ; s += 2)
		std::cout << "(" << std::setw(2) << s << std::setw(1) << ") " << std::setw(5) << utility[s] << " (" << action(policy[s]) << ")    "; 
	std::cout << std::endl << std::endl << std::endl;
	
	std::cout << "    ";
	for (int s = 16 ; s <= 28 ; s += 2)
		std::cout << "(" << std::setw(2) << s << std::setw(1) << ") " << std::setw(5) << utility[s] << " (" << action(policy[s]) << ")    "; 
	std::cout << std::endl;
	
	std::cout << "    ";
	for (int s = 17 ; s <= 29 ; s += 2)
		std::cout << "(" << std::setw(2) << s << std::setw(1) << ") " << std::setw(5) << utility[s] << " (" << action(policy[s]) << ")    "; 
	std::cout << std::endl << std::endl;

	std::cout << std::fixed << std::setprecision(1) << "Solution Technique: " << iterStrings[iter] << std::endl << std::endl
			  << "Discount Factor = " << setprecision(6) << discount << std::endl
			  << "Max Error in State Utilities = " << epsilon << std::endl
			  << "Positive Reward = " << posTerminal << std::endl
			  << "Negative Reward = " << negTerminal << std::endl
			  << "Step Cost = " << stepCost << std::endl << std::endl
			  << "# Iterations: " << numIter << std::endl
			  << "Solution Time: " << std::setprecision(8) << solTime << " seconds" << std::endl
			  << "EXITING" << std::endl << std::endl << std::endl;

	return;
}
    
std::string action(int a)
{
	switch (a) {
	    
	   case N: return "N";
	    
	   case E: return "E";
	    
	   case S: return "S";
	    
	   case W: return "W";

	   default: return "X";
	}
}

void initMDP(double negTerminal, double posTerminal, double stepCost, double keyLoss)
{

	// initializing rewards
	for (int s = 0; s < NUM_STATES; s++) {
		R[s] = stepCost;
	}

	// reset the rewards for terminal states
 	R[44] = negTerminal;
 	R[45] = negTerminal;

 	R[48] = negTerminal;
 	R[49] = negTerminal;

 	R[28] = posTerminal;


	// initialize all transition probabilities to 0.0
	for (int s1 = 0 ; s1 < NUM_STATES ; s1++)
	    for (int a = 0 ; a < NUM_ACTIONS ; a++)
			for (int s2 = 0 ; s2 < NUM_STATES ; s2++)
		    	T[s1][a][s2] = 0.0;

	// reset those transition probabilities that are NOT 0
	T[0][N][0] = 0.1;
	T[0][N][30] = 0.8;
	T[0][N][2] = 0.1;

	T[0][E][30] = 0.1;
	T[0][E][2] = 0.8;
	T[0][E][0] = 0.1;

	T[0][S][2] = 0.1;
	T[0][S][0] = 0.9;

	T[0][W][0] = 0.9;
	T[0][W][30] = 0.1;


	T[1][N][1] = 0.1;
	T[1][N][31] = 0.8;
	T[1][N][3] = 0.1;

	T[1][E][31] = 0.1;
	T[1][E][3] = 0.8;
	T[1][E][1] = 0.1;

 	T[1][S][3] = 0.1;
 	T[1][S][1] = 0.9;

	T[1][W][1] = 0.9;
	T[1][W][31] = 0.1;


	T[2][N][0] = 0.1;
	T[2][N][32] = 0.8;
	T[2][N][4] = 0.1;

	T[2][E][32] = 0.1;
	T[2][E][4] = 0.8;
	T[2][E][2] = 0.1;

	T[2][S][4] = 0.1;
	T[2][S][2] = 0.8;
	T[2][S][0] = 0.1;

	T[2][W][2] = 0.1;
	T[2][W][0] = 0.8;
	T[2][W][32] = 0.1;


	T[3][N][1] = 0.1;
	T[3][N][33] = 0.8;
	T[3][N][5] = 0.1;

	T[3][E][33] = 0.1;
	T[3][E][5] = 0.8;
	T[3][E][3] = 0.1;

	T[3][S][5] = 0.1;
	T[3][S][3] = 0.8;
	T[3][S][1] = 0.1;

	T[3][W][3] = 0.1;
	T[3][W][1] = 0.8;
	T[3][W][33] = 0.1;


	T[4][N][2] = 0.1;
	T[4][N][34] = 0.8;
	T[4][N][6] = 0.1;

	T[4][E][34] = 0.1;
	T[4][E][6] = 0.8;
	T[4][E][4] = 0.1;

	T[4][S][6] = 0.1;
	T[4][S][4] = 0.8;
	T[4][S][2] = 0.1;

	T[4][W][4] = 0.1;
	T[4][W][2] = 0.8;
	T[4][W][34] = 0.1;


	T[5][N][3] = 0.1;
	T[5][N][35] = 0.8;
	T[5][N][7] = 0.1;

	T[5][E][35] = 0.1;
	T[5][E][7] = 0.8;
	T[5][E][5] = 0.1;

	T[5][S][7] = 0.1;
	T[5][S][5] = 0.8;
	T[5][S][3] = 0.1;

	T[5][W][5] = 0.1;
	T[5][W][3] = 0.8;
	T[5][W][35] = 0.1;


	T[6][N][4] = 0.1;
	T[6][N][36] = 0.8;
	T[6][N][8] = 0.1;

	T[6][E][36] = 0.1;
	T[6][E][8] = 0.8;
	T[6][E][6] = 0.1;

	T[6][S][8] = 0.1;
	T[6][S][6] = 0.8;
	T[6][S][4] = 0.1;

	T[6][W][6] = 0.1;
	T[6][W][4] = 0.8;
	T[6][W][36] = 0.1;


	T[7][N][5] = 0.1;
	T[7][N][37] = 0.8;
	T[7][N][9] = 0.1;

	T[7][E][37] = 0.1;
	T[7][E][9] = 0.8;
	T[7][E][7] = 0.1;

	T[7][S][9] = 0.1;
	T[7][S][7] = 0.8;
	T[7][S][5] = 0.1;

	T[7][W][7] = 0.1;
	T[7][W][5] = 0.8;
	T[7][W][37] = 0.1;


	T[8][N][6] = 0.1;
	T[8][N][38] = 0.8;
	T[8][N][10] = 0.1;

	T[8][E][38] = 0.1;
	T[8][E][10] = 0.8;
	T[8][E][8] = 0.1;

	T[8][S][10] = 0.1;
	T[8][S][8] = 0.8;
	T[8][S][6] = 0.1;

	T[8][W][8] = 0.1;
	T[8][W][6] = 0.8;
	T[8][W][38] = 0.1;


	T[9][N][7] = 0.1;
	T[9][N][39] = 0.8;
	T[9][N][11] = 0.1;

	T[9][E][39] = 0.1;
	T[9][E][11] = 0.8;
	T[9][E][9] = 0.1;

	T[9][S][11] = 0.1;
	T[9][S][9] = 0.8;
	T[9][S][7] = 0.1;

	T[9][W][9] = 0.1;
	T[9][W][7] = 0.8;
	T[9][W][39] = 0.1;


	T[10][N][8] = 0.1;
	T[10][N][10] = 0.8;
	T[10][N][12] = 0.1;

	T[10][E][10] = 0.2;
	T[10][E][12] = 0.8;

	T[10][S][12] = 0.1;
	T[10][S][10] = 0.8;
	T[10][S][8] = 0.1;

	T[10][W][10] = 0.2;
	T[10][W][8] = 0.8;


	T[11][N][9] = 0.1;
	T[11][N][11] = 0.8;
	T[11][N][13] = 0.1;

	T[11][E][11] = 0.2;
	T[11][E][13] = 0.8;

	T[11][S][13] = 0.1;
	T[11][S][11] = 0.8;
	T[11][S][9] = 0.1;

	T[11][W][11] = 0.2;
	T[11][W][9] = 0.8;


	T[12][N][10] = 0.1;
	T[12][N][12] = 0.8;
	T[12][N][14] = 0.1;

	T[12][E][12] = 0.2;
	T[12][E][14] = 0.8;

	T[12][S][14] = 0.1;
	T[12][S][12] = 0.8;
	T[12][S][10] = 0.1;

	T[12][W][12] = 0.2;
	T[12][W][10] = 0.8;


	T[13][N][11] = 0.1;
	T[13][N][13] = 0.8;
	T[13][N][15] = 0.1;

	T[13][E][13] = 0.2;
	T[13][E][15] = 0.8;

	T[13][S][15] = 0.1;
	T[13][S][13] = 0.8;
	T[13][S][11] = 0.1;

	T[13][W][13] = 0.2;
	T[13][W][11] = 0.8;


	T[14][N][12] = 0.1;
	T[14][N][14] = 0.8;
	T[14][N][16] = 0.1;

	T[14][E][14] = 0.2;
	T[14][E][16] = 0.8;

	T[14][S][16] = 0.1;
	T[14][S][14] = 0.8;
	T[14][S][12] = 0.1;

	T[14][W][14] = 0.2;
	T[14][W][12] = 0.8;


	T[15][N][13] = 0.1;
	T[15][N][15] = 0.8;
	T[15][N][17] = 0.1;

	T[15][E][15] = 0.2;
	T[15][E][17] = 0.8;

	T[15][S][17] = 0.1;
	T[15][S][15] = 0.8;
	T[15][S][13] = 0.1;

	T[15][W][15] = 0.2;
	T[15][W][13] = 0.8;


	T[16][N][14] = 0.1;
	T[16][N][16] = 0.8;
	T[16][N][18] = 0.1;

	T[16][E][16] = 0.2;
	T[16][E][18] = 0.8;

	T[16][S][18] = 0.1;
	T[16][S][16] = 0.8;
	T[16][S][14] = 0.1;

	T[16][W][16] = 0.2;
	T[16][W][14] = 0.8;


	T[17][N][15] = 0.1;
	T[17][N][17] = 0.8;
	T[17][N][19] = 0.1;

	T[17][E][17] = 0.2;
	T[17][E][19] = 0.8;

	T[17][S][19] = 0.1;
	T[17][S][17] = 0.8;
	T[17][S][15] = 0.1;

	T[17][W][17] = 0.2;
	T[17][W][15] = 0.8;


	T[18][N][16] = 0.1;
	T[18][N][18] = 0.8;
	T[18][N][20] = 0.1;

	T[18][E][18] = 0.2;
	T[18][E][20] = 0.8;

	T[18][S][20] = 0.1;
	T[18][S][18] = 0.8;
	T[18][S][16] = 0.1;

	T[18][W][18] = 0.2;
	T[18][W][16] = 0.8;


	T[19][N][17] = 0.1;
	T[19][N][19] = 0.8;
	T[19][N][21] = 0.1;

	T[19][E][19] = 0.2;
	T[19][E][21] = 0.8;

	T[19][S][21] = 0.1;
	T[19][S][19] = 0.8;
	T[19][S][17] = 0.1;

	T[19][W][19] = 0.2;
	T[19][W][17] = 0.8;



	T[20][N][18] = 0.1;
	T[20][N][20] = 0.8;
	T[20][N][22] = 0.1;

	T[20][E][20] = 0.2;
	T[20][E][22] = 0.8;

	T[20][S][22] = 0.1;
	T[20][S][20] = 0.8;
	T[20][S][18] = 0.1;

	T[20][W][20] = 0.2;
	T[20][W][18] = 0.8;


	T[21][N][19] = 0.1;
	T[21][N][21] = 0.8;
	T[21][N][23] = 0.1;

	T[21][E][21] = 0.2;
	T[21][E][23] = 0.8;

	T[21][S][23] = 0.1;
	T[21][S][21] = 0.8;
	T[21][S][19] = 0.1;

	T[21][W][21] = 0.2;
	T[21][W][19] = 0.8;


	T[22][N][20] = 0.1;
	T[22][N][22] = 0.8;
	T[22][N][24] = 0.1;

	T[22][E][22] = 0.2;
	T[22][E][24] = 0.8;

	T[22][S][24] = 0.1;
	T[22][S][22] = 0.8;
	T[22][S][20] = 0.1;

	T[22][W][22] = 0.2;
	T[22][W][20] = 0.8;


	T[23][N][21] = 0.1;
	T[23][N][23] = 0.8;
	T[23][N][25] = 0.1;

	T[23][E][23] = 0.2;
	T[23][E][25] = 0.8;

	T[23][S][25] = 0.1;
	T[23][S][23] = 0.8;
	T[23][S][21] = 0.1;

	T[23][W][23] = 0.2;
	T[23][W][21] = 0.8;


	T[24][N][22] = 0.1;
	T[24][N][24] = 0.8;
	T[24][N][26] = 0.1;

	T[24][E][24] = 0.2;
	T[24][E][26] = 0.8;

	T[24][S][26] = 0.1;
	T[24][S][24] = 0.8;
	T[24][S][22] = 0.1;

	T[24][W][24] = 0.2;
	T[24][W][22] = 0.8;


	T[25][N][23] = 0.1;
	T[25][N][25] = 0.8;
	T[25][N][27] = 0.1;

	T[25][E][25] = 0.2;
	T[25][E][27] = 0.8;

	T[25][S][27] = 0.1;
	T[25][S][25] = 0.8;
	T[25][S][23] = 0.1;

	T[25][W][25] = 0.2;
	T[25][W][23] = 0.8;


	T[26][N][24] = 0.1;
	T[26][N][26] = 0.8;
	T[26][N][28] = 0.1;

	T[26][E][26] = 0.2;
	T[26][E][28] = 0.8;

	T[26][S][28] = 0.1;
	T[26][S][26] = 0.8;
	T[26][S][24] = 0.1;

	T[26][W][26] = 0.2;
	T[26][W][24] = 0.8;


	T[27][N][25] = 0.1;
	T[27][N][27] = 0.8;
	T[27][N][29] = 0.1;

	T[27][E][27] = 0.2;
	T[27][E][29] = 0.8;

	T[27][S][29] = 0.1;
	T[27][S][27] = 0.8;
	T[27][S][25] = 0.1;

	T[27][W][27] = 0.2;
	T[27][W][25] = 0.8;


	// no transitions from states 28 and 29


	T[30][N][30] = 0.1;
	T[30][N][40] = 0.8 * (1.0 - keyLoss);
	T[30][N][41] = 0.8 * keyLoss;
	T[30][N][32] = 0.1;

	T[30][E][40] = 0.1 * (1.0 - keyLoss);
	T[30][E][41] = 0.1 * keyLoss;
	T[30][E][32] = 0.8;
	T[30][E][0] = 0.1;

	T[30][S][32] = 0.1;
	T[30][S][0] = 0.8;
	T[30][S][30] = 0.1;

	T[30][W][0] = 0.1;
	T[30][W][30] = 0.8;
	T[30][W][40] = 0.1 * (1.0 - keyLoss);
	T[30][W][41] = 0.1 * keyLoss;


	T[31][N][31] = 0.1;
	T[31][N][41] = 0.8;
	T[31][N][33] = 0.1;

	T[31][E][41] = 0.1;
	T[31][E][33] = 0.8;
	T[31][E][1] = 0.1;

	T[31][S][33] = 0.1;
	T[31][S][1] = 0.8;
	T[31][S][31] = 0.1;

	T[31][W][1] = 0.1;
	T[31][W][31] = 0.8;
	T[31][W][41] = 0.1;


	T[32][N][30] = 0.1;
	T[32][N][42] = 0.8;
	T[32][N][34] = 0.1;

	T[32][E][42] = 0.1;
	T[32][E][34] = 0.8;
	T[32][E][2] = 0.1;

	T[32][S][34] = 0.1;
	T[32][S][2] = 0.8;
	T[32][S][30] = 0.1;

	T[32][W][2] = 0.1;
	T[32][W][30] = 0.8;
	T[32][W][42] = 0.1;


	T[33][N][31] = 0.1;
	T[33][N][43] = 0.8;
	T[33][N][35] = 0.1;

	T[33][E][43] = 0.1;
	T[33][E][35] = 0.8;
	T[33][E][3] = 0.1;

	T[33][S][35] = 0.1;
	T[33][S][3] = 0.8;
	T[33][S][31] = 0.1;

	T[33][W][3] = 0.1;
	T[33][W][31] = 0.8;
	T[33][W][43] = 0.1;


	T[34][N][32] = 0.1;
	T[34][N][44] = 0.8;
	T[34][N][36] = 0.1;

	T[34][E][44] = 0.1;
	T[34][E][36] = 0.8;
	T[34][E][4] = 0.1;

	T[34][S][36] = 0.1;
	T[34][S][4] = 0.8;
	T[34][S][32] = 0.1;

	T[34][W][4] = 0.1;
	T[34][W][32] = 0.8;
	T[34][W][44] = 0.1;


	T[35][N][33] = 0.1;
	T[35][N][45] = 0.8;
	T[35][N][37] = 0.1;

	T[35][E][45] = 0.1;
	T[35][E][37] = 0.8;
	T[35][E][5] = 0.1;

	T[35][S][37] = 0.1;
	T[35][S][5] = 0.8;
	T[35][S][33] = 0.1;

	T[35][W][5] = 0.1;
	T[35][W][33] = 0.8;
	T[35][W][45] = 0.1;


	T[36][N][34] = 0.1;
	T[36][N][46] = 0.8;
	T[36][N][38] = 0.1;

	T[36][E][46] = 0.1;
	T[36][E][38] = 0.8;
	T[36][E][6] = 0.1;

	T[36][S][38] = 0.1;
	T[36][S][6] = 0.8;
	T[36][S][34] = 0.1;

	T[36][W][6] = 0.1;
	T[36][W][34] = 0.8;
	T[36][W][46] = 0.1;


	T[37][N][35] = 0.1;
	T[37][N][47] = 0.8;
	T[37][N][39] = 0.1;

	T[37][E][47] = 0.1;
	T[37][E][39] = 0.8;
	T[37][E][7] = 0.1;

	T[37][S][39] = 0.1;
	T[37][S][7] = 0.8;
	T[37][S][35] = 0.1;

	T[37][W][7] = 0.1;
	T[37][W][35] = 0.8;
	T[37][W][47] = 0.1;


	T[38][N][36] = 0.1;
	T[38][N][48] = 0.8;
	T[38][N][38] = 0.1;

	T[38][E][48] = 0.1;
	T[38][E][38] = 0.8;
	T[38][E][8] = 0.1;

	T[38][S][38] = 0.1;
	T[38][S][8] = 0.8;
	T[38][S][36] = 0.1;

	T[38][W][8] = 0.1;
	T[38][W][36] = 0.8;
	T[38][W][48] = 0.1;


	T[39][N][37] = 0.1;
	T[39][N][49] = 0.8;
	T[39][N][39] = 0.1;

	T[39][E][49] = 0.1;
	T[39][E][39] = 0.8;
	T[39][E][9] = 0.1;

	T[39][S][39] = 0.1;
	T[39][S][9] = 0.8;
	T[39][S][37] = 0.1;

	T[39][W][9] = 0.1;
	T[39][W][37] = 0.8;
	T[39][W][49] = 0.1;


	T[40][N][40] = 0.1 * (1.0 - keyLoss);
	T[40][N][41] = 0.1 * keyLoss;
	T[40][N][50] = 0.8;
	T[40][N][42] = 0.1;

	T[40][E][50] = 0.1;
	T[40][E][42] = 0.8;
	T[40][E][30] = 0.1;

	T[40][S][42] = 0.1;
	T[40][S][30] = 0.8;
	T[40][S][40] = 0.1 * (1.0 - keyLoss);
	T[40][S][41] = 0.1 * keyLoss;

	T[40][W][30] = 0.1;
	T[40][W][40] = 0.8 * (1.0 - keyLoss);
	T[40][W][41] = 0.8 * keyLoss;
	T[40][W][50] = 0.1;


	T[41][N][41] = 0.1;
	T[41][N][51] = 0.8;
	T[41][N][43] = 0.1;

	T[41][E][51] = 0.1;
	T[41][E][43] = 0.8;
	T[41][E][31] = 0.1;

	T[41][S][43] = 0.1;
	T[41][S][31] = 0.8;
	T[41][S][41] = 0.1;

	T[41][W][31] = 0.1;
	T[41][W][41] = 0.8;
	T[41][W][51] = 0.1;


	T[42][N][40] = 0.1 * (1.0 - keyLoss);
	T[42][N][41] = 0.1 * keyLoss;
	T[42][N][52] = 0.8;
	T[42][N][44] = 0.1;

	T[42][E][52] = 0.1;
	T[42][E][44] = 0.8;
	T[42][E][32] = 0.1;

 	T[42][S][44] = 0.1;
	T[42][S][32] = 0.8;
  	T[42][S][40] = 0.1 * (1.0 - keyLoss);
  	T[42][S][41] = 0.1 * keyLoss;
	
	T[42][W][32] = 0.1;
	T[42][W][40] = 0.8 * (1.0 - keyLoss);
	T[42][W][41] = 0.8 * keyLoss;
	T[42][W][52] = 0.1;


	T[43][N][41] = 0.1;
	T[43][N][53] = 0.8;
	T[43][N][45] = 0.1;

	T[43][E][53] = 0.1;
	T[43][E][45] = 0.8;
	T[43][E][33] = 0.1;

	T[43][S][45] = 0.1;
	T[43][S][33] = 0.8;
	T[43][S][41] = 0.1;

	T[43][W][33] = 0.1;
	T[43][W][41] = 0.8;
	T[43][W][53] = 0.1;

	
	// no transitions from states 44 and 45


	T[46][N][44] = 0.1;
	T[46][N][56] = 0.8;
	T[46][N][48] = 0.1;

	T[46][E][56] = 0.1;
	T[46][E][48] = 0.8;
	T[46][E][36] = 0.1;

	T[46][S][48] = 0.1;
	T[46][S][36] = 0.8;
	T[46][S][44] = 0.1;

	T[46][W][36] = 0.1;
	T[46][W][44] = 0.8;
	T[46][W][56] = 0.1;


	T[47][N][45] = 0.1;
	T[47][N][57] = 0.8;
	T[47][N][49] = 0.1;

	T[47][E][57] = 0.1;
	T[47][E][49] = 0.8;
	T[47][E][37] = 0.1;

	T[47][S][49] = 0.1;
	T[47][S][37] = 0.8;
	T[47][S][45] = 0.1;

	T[47][W][37] = 0.1;
	T[47][W][45] = 0.8;
	T[47][W][57] = 0.1;


	// no transitions from states 48 and 49


	T[50][N][50] = 0.1;
	T[50][N][58] = 0.8;
	T[50][N][52] = 0.1;

	T[50][E][58] = 0.1;
	T[50][E][52] = 0.8;
	T[50][E][40] = 0.1 * (1.0 - keyLoss);
	T[50][E][41] = 0.1 * keyLoss;

	T[50][S][52] = 0.1;
	T[50][S][40] = 0.8 * (1.0 - keyLoss);
	T[50][S][41] = 0.8 * keyLoss;
	T[50][S][50] = 0.1;

	T[50][W][40] = 0.1 * (1.0 - keyLoss);
	T[50][W][41] = 0.1 * keyLoss;
	T[50][W][50] = 0.8;
	T[50][W][58] = 0.1;


	T[51][N][51] = 0.1;
	T[51][N][59] = 0.8;
	T[51][N][53] = 0.1;

	T[51][E][59] = 0.1;
	T[51][E][53] = 0.8;
	T[51][E][41] = 0.1;

	T[51][S][53] = 0.1;
	T[51][S][41] = 0.8;
	T[51][S][51] = 0.1;

	T[51][W][41] = 0.1;
	T[51][W][51] = 0.8;
	T[51][W][59] = 0.1;


	T[52][N][50] = 0.1;
	T[52][N][60] = 0.8;
	T[52][N][54] = 0.1;

	T[52][E][60] = 0.1;
	T[52][E][54] = 0.8;
	T[52][E][42] = 0.1;

	T[52][S][54] = 0.1;
	T[52][S][42] = 0.8;
	T[52][S][50] = 0.1;

	T[52][W][42] = 0.1;
	T[52][W][50] = 0.8;
	T[52][W][60] = 0.1;


	T[53][N][51] = 0.1;
	T[53][N][61] = 0.8;
	T[53][N][55] = 0.1;

	T[53][E][61] = 0.1;
	T[53][E][55] = 0.8;
	T[53][E][43] = 0.1;

	T[53][S][55] = 0.1;
	T[53][S][43] = 0.8;
	T[53][S][51] = 0.1;

	T[53][W][43] = 0.1;
	T[53][W][51] = 0.8;
	T[53][W][61] = 0.1;


	T[54][N][52] = 0.1;
	T[54][N][62] = 0.8;
	T[54][N][56] = 0.1;

	T[54][E][62] = 0.1;
	T[54][E][56] = 0.8;
	T[54][E][44] = 0.1;

	T[54][S][56] = 0.1;
	T[54][S][44] = 0.8;
	T[54][S][52] = 0.1;

	T[54][W][44] = 0.1;
	T[54][W][52] = 0.8;
	T[54][W][62] = 0.1;


	T[55][N][53] = 0.1;
	T[55][N][63] = 0.8;
	T[55][N][57] = 0.1;

	T[55][E][63] = 0.1;
	T[55][E][57] = 0.8;
	T[55][E][45] = 0.1;

	T[55][S][57] = 0.1;
	T[55][S][45] = 0.8;
	T[55][S][53] = 0.1;

	T[55][W][45] = 0.1;
	T[55][W][53] = 0.8;
	T[55][W][63] = 0.1;


	T[56][N][54] = 0.1;
	T[56][N][64] = 0.8;
	T[56][N][56] = 0.1;

	T[56][E][64] = 0.1;
	T[56][E][56] = 0.8;
	T[56][E][46] = 0.1;

	T[56][S][56] = 0.1;
	T[56][S][46] = 0.8;
	T[56][S][54] = 0.1;

	T[56][W][46] = 0.1;
	T[56][W][54] = 0.8;
	T[56][W][64] = 0.1;


 	T[57][N][55] = 0.1;
 	T[57][N][64] = 0.8;
 	T[57][N][57] = 0.1;

	T[57][E][64] = 0.1;
	T[57][E][57] = 0.8;
	T[57][E][47] = 0.1;

	T[57][S][57] = 0.1;
	T[57][S][47] = 0.8;
	T[57][S][55] = 0.1;

	T[57][W][47] = 0.1;
	T[57][W][55] = 0.8;
	T[57][W][64] = 0.1;


	T[58][N][58] = 0.9;
	T[58][N][60] = 0.1;

	T[58][E][58] = 0.1;
	T[58][E][60] = 0.8;
	T[58][E][50] = 0.1;

	T[58][S][60] = 0.1;
	T[58][S][50] = 0.8;
	T[58][S][58] = 0.1;

	T[58][W][50] = 0.1;
	T[58][W][58] = 0.9;


	T[59][N][59] = 0.9;
	T[59][N][61] = 0.1;

	T[59][E][59] = 0.1;
	T[59][E][61] = 0.8;
	T[59][E][51] = 0.1;

	T[59][S][61] = 0.1;
	T[59][S][51] = 0.8;
	T[59][S][59] = 0.1;

	T[59][W][51] = 0.1;
	T[59][W][59] = 0.9;



	T[60][N][58] = 0.1;
	T[60][N][60] = 0.8;
	T[60][N][62] = 0.1;

	T[60][E][60] = 0.1;
	T[60][E][62] = 0.8;
	T[60][E][52] = 0.1;

	T[60][S][62] = 0.1;
	T[60][S][52] = 0.8;
	T[60][S][58] = 0.1;

	T[60][W][52] = 0.1;
	T[60][W][58] = 0.8;
	T[60][W][60] = 0.1;


	T[61][N][59] = 0.1;
	T[61][N][61] = 0.8;
	T[61][N][63] = 0.1;

	T[61][E][61] = 0.1;
	T[61][E][63] = 0.8;
	T[61][E][53] = 0.1;

	T[61][S][63] = 0.1;
	T[61][S][53] = 0.8;
	T[61][S][59] = 0.1;

	T[61][W][53] = 0.1;
	T[61][W][59] = 0.8;
	T[61][W][61] = 0.1;


	T[62][N][60] = 0.1;
	T[62][N][62] = 0.8;
	T[62][N][64] = 0.1;

	T[62][E][62] = 0.1;
	T[62][E][64] = 0.8;
	T[62][E][54] = 0.1;

	T[62][S][64] = 0.1;
	T[62][S][54] = 0.8;
	T[62][S][60] = 0.1;

	T[62][W][54] = 0.1;
	T[62][W][60] = 0.8;
	T[62][W][62] = 0.1;


	T[63][N][61] = 0.1;
	T[63][N][63] = 0.8;
	T[63][N][64] = 0.1;

	T[63][E][63] = 0.1;
	T[63][E][64] = 0.8;
	T[63][E][55] = 0.1;

	T[63][S][64] = 0.1;
	T[63][S][55] = 0.8;
	T[63][S][61] = 0.1;

	T[63][W][55] = 0.1;
	T[63][W][61] = 0.8;
	T[63][W][63] = 0.1;


	T[64][N][62] = 0.1;
	T[64][N][64] = 0.9;

	T[64][E][64] = 0.9;
	T[64][E][56] = 0.1;

	T[64][S][64] = 0.1;
	T[64][S][56] = 0.8;
	T[64][S][62] = 0.1;

	T[64][W][56] = 0.1;
	T[64][W][62] = 0.8;
	T[64][W][64] = 0.1;

	return;
}