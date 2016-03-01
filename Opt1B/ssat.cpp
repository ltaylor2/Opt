#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

enum SolutionType { naive };

// prototypes
int readSSATFile(std::string fileName, std::vector<double>*, std::vector<std::vector<int>>*);
void solve(SolutionType, std::vector<double>*, std::vector<std::vector<int>>*, std::vector<bool>*, double*);

int main(int argc, char* argv[]) {

	// read in cmd line arguments
	if (argc != 3) {
		std::cout << "Invalid Arguments (" << argc << ") Exiting." << std::endl;
		return 1;
	}

	SolutionType directions;
	if (std::string(argv[1]).compare("n") == 0) {
		directions = SolutionType::naive;
	}
	else {
		std::cout << "---" << argv[1] << "---" << std::endl;
		std::cout << "Incorrect solving directions. Exiting." << std::endl;
		return 1;
	}

	std::string fileName = std::string(argv[2]);

	// data structures
	std::vector<double> variables;
	std::vector<bool> assignments;

	std::vector<std::vector<int>> clauses;

	double solutionProb;
	// file IO
	if (readSSATFile(fileName, &variables, &clauses) == 1) {
		return 1;
	}

	solve(directions, &variables, &clauses, &assignments, &solutionProb);
	// call appropriate solver





	// report results	

	return 0;
}

int readSSATFile(std::string fileName,
				std::vector<double>* variables, 
				std::vector<std::vector<int>>* clauses)
{
	std::cout << "Reading in " << fileName << std::endl;

	std::ifstream file(fileName);

	if (!file) {
		std::cout << "Failed to open file. Exiting." << std::endl;
		return 1;
	}

	std::string line;
	while(getline(file, line)) {
		if (line.compare("variables") == 0) {
			getline(file, line);
			while(line.compare("") != 0) {
				std::stringstream ss(line);
				int varName = 0;
				double varValue = 0.0;
				ss >> varName >> varValue;
				variables->push_back(varValue);
				getline(file, line);
			}
		}
		if (line.compare("clauses") == 0) {
			getline(file, line);
			while(line.compare("") != 0) {
				std::vector<int> clause;
				std::stringstream ss(line);
				int literal = variables->size() + 1;
				ss >> literal;
				while (literal != 0) {
					clause.push_back(literal);
					ss >> literal;
					std::cout << literal << "   ";
				}
				std::cout << "\n";
				clauses->push_back(clause);
				getline(file,line);
			}
		}
	}

	return 0;
}

void solve(SolutionType directions,
		   std::vector<double>* variables,
		   std::vector<std::vector<int>>* clauses,
		   std::vector<bool>* assignments,
		   double* solutionProb)
{
}