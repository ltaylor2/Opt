#include <vector>
#include <iostream>

struct Variable {
	int name;
	double probability;
};

enum SolutionType { naive };

// prototypes
void readSSFile(std::vector<Variable>*, std::vector<std::vector<int>>);
void solve(SolutionType, std::vector<Variable>*, std::vector<std::vector<int>>*, std::vector<bool>*, double*);

int main(int argc, char* argv[]) {

	if (argc != 3) {
		std::cout << "Invalid Arguments. Exiting." << std::endl;
		return 1;
	}

	SolutionType directions;
	if (argv[1] == 'n') {
		directions = SolutionType::naive;
	}
	else {
		std::cout << "---" << argv[1] << "---" << std::endl;
		std::cout << "Incorrect solving directions. Exiting." << std::endl;
		return 1;
	}

	std::string fileName = argv[2];
	
	// data structures
	std::vector<Variable> variables;
	std::vector<bool> assignments;

	std::vector<std::vector<int>> clauses;

	double solutionProb;
	// file IO

	std::cout << "Hello World!" << std::endl;

	solve(directions, &variables, &clauses, &assignments, &solutionProb);
	// call appropriate solver



	// report results	

	return 0;
}

void readSSFile(std::vector<Variable>* variables, 
				std::vector<std::vector<int>>* clauses)
{
	// read in the file

	// fill the argument vectors
	
	return;
}

void solve(SolutionType directions,
		   std::vector<Variable>* variables,
		   std::vector<std::vector<int>>* clauses,
		   std::vector<bool>* assignments,
		   double* solutionProb)
{
	if (directions == SolutionType::naive)
		std::cout << "testing complete" << std::endl;
	else
		std::cout << "whoops" << std::endl;
	return;
}