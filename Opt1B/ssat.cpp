#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>

enum SolutionType { naive, unit };

// prototypes
int readSSATFile(std::string fileName, std::vector<double>*, std::vector<std::vector<int>>*, std::vector<std::vector<int>>*);
double solve(SolutionType, std::vector<double>*, std::vector<std::vector<int>>*, std::vector<int>*, std::vector<int>, std::vector<std::vector<int>>*);
void satisfyClauses(int, std::vector<std::vector<int>>*, std::vector<int>*, std::vector<int>*);

int main(int argc, char* argv[])
{
	// read in cmd line arguments
	if (argc != 3) {
		std::cout << "Invalid Arguments (" << argc << ") Exiting." << std::endl;
		return 1;
	}

	SolutionType directions;
	if (std::string(argv[1]).compare("n") == 0)
		directions = SolutionType::naive;
	else if (std::string(argv[1]).compare("u") == 0)
		directions = SolutionType::unit;
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

	std::cout << "Beginning to solve!" << std::endl;
	clock_t start = clock();
	double solutionProb = solve(directions, &variables, &clauses, &clauseSats, assignments, &varsByClause);
	clock_t end = clock();

	double solveTime  = (double)(end-start) / CLOCKS_PER_SEC;
	std::cout << "Solution is: " << solutionProb << " (found in " << solveTime << " seconds)" << std::endl;

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


	// find the next unit clause (regardless of ordering)
	if (directions == SolutionType::unit) {
		std::vector<std::vector<int>> unitCopy(*clauses);
		std::vector<int> unitSats(*clauseSats);

		int unitVar	= 0;
		for (unsigned int c = 0; c < unitCopy.size(); c++) {
			if (unitCopy[c].size() == 1 && unitSats[c] != 1) {
				unitVar = unitCopy[c][0];
				break;
			}
		}
		// if you've found a unit clause, proceed
		if (unitVar != 0) {
			assignments[abs(unitVar)-1] = unitVar / abs(unitVar);
			satisfyClauses(abs(unitVar)-1, &unitCopy, &unitSats, &assignments);
			double probSatUnit = solve(directions, variables, &unitCopy, &unitSats, assignments, varsByClause);
			if (variables->at(abs(unitVar)-1) == -1)
				return probSatUnit;
			else {
				if (assignments[abs(unitVar)-1] == 1)
					return probSatUnit * variables->at(abs(unitVar)-1);
				else
					return probSatUnit * (1 - variables->at(abs(unitVar)-1));
			}
		}
	}

	// there is guaranteed to be a 0 in assignments, because if there was not we would have retunred from allSat == TRUE
	int nextVarIndex = std::distance(assignments.begin(), std::find(assignments.begin(), assignments.end(), 0));
	// trying false
	assignments[nextVarIndex] = -1;
	std::vector<std::vector<int>> falseCopy(*clauses);	// a copy so we don't have to revert changes
	std::vector<int> falseSats(*clauseSats);
	// satisfy clauses
	satisfyClauses(nextVarIndex, &falseCopy, &falseSats, &assignments);
	double probSatFalse = solve(directions, variables, &falseCopy, &falseSats, assignments, varsByClause);

	// trying true
	assignments[nextVarIndex] = 1;
	std::vector<std::vector<int>> trueCopy(*clauses);	// a copy so we don't have to revert changes
	std::vector<int> trueSats(*clauseSats);
	// satisfy clauses
	satisfyClauses(nextVarIndex, &trueCopy, &trueSats, &assignments);
	double probSatTrue = solve(directions, variables, &trueCopy, &trueSats, assignments, varsByClause);

	if (variables->at(nextVarIndex) == -1) { 	// v is a choice variable
		return std::max(probSatFalse, probSatTrue);
	}
	// v is a chance variable
	return probSatTrue * variables->at(nextVarIndex) + probSatFalse * (1 - variables->at(nextVarIndex));
}

void satisfyClauses(int varIndex, std::vector<std::vector<int>>* clauses, std::vector<int>* sats, std::vector<int>* assignments)
{
	for (unsigned int c = 0; c < clauses->size(); c++) {
		if (sats->at(c) == 1)
			continue;
		for (unsigned int l = 0; l < clauses->at(c).size(); l++) {
			if (clauses->at(c).at(l) == (varIndex + 1) * assignments->at(varIndex))
				sats->at(c) = 1;
			else if (clauses->at(c).at(l) == (varIndex + 1) * assignments->at(varIndex) * -1) {
				if (clauses->at(c).size() == 1)
					sats->at(c) = -1;
				else {
					clauses->at(c).erase(clauses->at(c).begin() + l);	
					l--;
				}
			}
		}
	}

}