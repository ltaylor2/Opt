#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

enum SolutionType { naive };

// prototypes
int readSSATFile(std::string fileName, std::vector<double>*, std::vector<std::vector<int>>*, std::vector<std::vector<int>>*);
double solve(SolutionType, std::vector<double>*, std::vector<std::vector<int>>*, std::vector<int>*, std::vector<int>, std::vector<std::vector<int>>*);

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
	std::vector<int> assignments;
	std::vector<std::vector<int>> clauses;
	std::vector<int> clauseSats;
	std::vector<std::vector<int>> varsByClause;

	// file IO
	if (readSSATFile(fileName, &variables, &clauses, &varsByClause) == 1) {
		return 1;
	}

	// fill assignments to match variables
	for (unsigned int i = 0; i < variables.size(); i++) {
		assignments.push_back(0);
	}
	// fill clauses for satisfaction
	for (unsigned int i = 0; i < clauses.size(); i++) {
		clauseSats.push_back(0);
	}

	double solutionProb = solve(directions, &variables, &clauses, &clauseSats, assignments, &varsByClause);
	std::cout << "Solution is:  " << solutionProb << std::endl;

	return 0;
}

int readSSATFile(std::string fileName,
				std::vector<double>* variables, 
				std::vector<std::vector<int>>* clauses,
				std::vector<std::vector<int>>* varsByClause)
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
				varsByClause->push_back(std::vector<int>());
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
					int litSign = 1;
					if (literal < 0)
						litSign = -1;
					varsByClause->at(abs(literal) - 1).push_back(clauses->size() * litSign);
					ss >> literal;
				}
				clauses->push_back(clause);

				getline(file,line);
			}
		}
	}

	return 0;
}

double solve(SolutionType directions,
		   std::vector<double>* variables,
		   std::vector<std::vector<int>>* clauses,
		   std::vector<int>* clauseSats,
		   std::vector<int> assignments,
		   std::vector<std::vector<int>>* varsByClause)
{
	// loop through clause satisfaction, checking for fully satisfied or any unsatisfied
	bool allSat = true;
	for (unsigned int i = 0; i < clauseSats->size(); i++) {
		if (clauseSats->at(i) == -1)		// any unsatisfied means that the SSAT is unsatisfied
			return 0.0;
		else if (clauseSats->at(i) == 0) {	// a non-satisfied (but not unsatisfied) clause
			allSat = false;
		}
	}
	if (allSat)			// no -1's or 0's, all are satisfied
		return 1.0;

	// there is guaranteed to be a 0 in assignments, because if there was not we would have retunred from allSat == TRUE
	int nextVarIndex = std::distance(assignments.begin(), std::find(assignments.begin(), assignments.end(), 0));
	// trying false
	assignments[nextVarIndex] = -1;
	std::vector<std::vector<int>> falseCopy(*clauses);	// a copy so we don't have to revert changes
	std::vector<int> falseSats(*clauseSats);
	// satisfy clauses
	for (unsigned int c = 0; c < falseCopy.size(); c++) {
		if (falseSats[c] == 1)
			continue;
		for (unsigned int l = 0; l < falseCopy[c].size(); l++) {
			if (falseCopy[c][l] == (nextVarIndex + 1) * assignments[nextVarIndex]) {
				falseSats[c] = 1;
			} else if (falseCopy[c][l] == (nextVarIndex + 1) * assignments[nextVarIndex] * -1) {
				if (falseCopy[c].size() == 1)
					falseSats[c] = -1;
				else {
					falseCopy[c].erase(falseCopy[c].begin() + l);	
					l--;
				}
			}
		}
	}
	double probSatFalse = solve(directions, variables, &falseCopy, &falseSats, assignments, varsByClause);

	// trying true
	assignments[nextVarIndex] = 1;
	std::vector<std::vector<int>> trueCopy(*clauses);	// a copy so we don't have to revert changes
	std::vector<int> trueSats(*clauseSats);
	// satisfy clauses
	for (unsigned int c = 0; c < trueCopy.size(); c++) {
		if (trueSats[c] == 1)
			continue;
		for (unsigned int l = 0; l < trueCopy[c].size(); l++) {
			if (trueCopy[c][l] == (nextVarIndex + 1) * assignments[nextVarIndex]) {
				trueSats[c] = 1;
			} else if (trueCopy[c][l] == (nextVarIndex + 1) * assignments[nextVarIndex] * -1) {
				if (trueCopy[c].size() == 1)
					trueSats[c] = -1;
				else {
					trueCopy[c].erase(trueCopy[c].begin() + l);
					l--;
				}
			}
		}
	}
	double probSatTrue = solve(directions, variables, &trueCopy, &trueSats, assignments, varsByClause);

	if (variables->at(nextVarIndex) == -1) { 	// v is a choice variable
		return std::max(probSatFalse, probSatTrue);
	}
	// v is a chance variable
	return probSatTrue * variables->at(nextVarIndex) + probSatFalse * (1 - variables->at(nextVarIndex));
}