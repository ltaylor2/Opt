#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>

enum SolutionType { naive, unit, pure };

// prototypes
int readSSATFile(std::string fileName, std::vector<double>*, std::vector<std::vector<int>>*, std::vector<std::vector<int>>*);
double solve(SolutionType, std::vector<double>*, std::vector<std::vector<int>>*, std::vector<int>*, std::vector<int>, std::vector<std::vector<int>>*);
void satisfyClauses(int, std::vector<std::vector<int>>*, std::vector<int>*, std::vector<int>*, std::vector<std::vector<int>>*);

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
	else if (std::string(argv[1]).compare("p") == 0)
		directions = SolutionType::pure;
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
		std::vector<std::vector<int>> unitVBC(*varsByClause);

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
			satisfyClauses(abs(unitVar)-1, &unitCopy, &unitSats, &assignments, &unitVBC);
			double probSatUnit = solve(directions, variables, &unitCopy, &unitSats, assignments, &unitVBC);
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

	// find pure choice variables
	if (directions == SolutionType::pure) {
		std::vector<std::vector<int>> pureCopy(*clauses);
		std::vector<int> pureSats(*clauseSats);
		std::vector<std::vector<int>> pureVBC(*varsByClause);

		int pureVar = -1;
		for (int l = 0; l < (signed int)pureVBC.size(); l++) {
			if (variables->at(l) != -1 || assignments[l] == 0)	// only looking at unassigned choice variables
				continue;
			for (unsigned int c = 1; c < varsByClause->at(l).size(); c++) {
				pureVar = l;
				if (!((pureVBC[l][c-1] < 0 && pureVBC[l][c] < 0) || (pureVBC[l][c-1] > 0 && pureVBC[l][c] > 0)) || pureVBC[l].size() <= 1) {
					pureVar = -1;	// know it's not pure
					break;		// so move on to the next variable
				}
			}

			if (pureVar == l)	// found one
				break;			// use it as the pure
		}

		if (pureVar != -1) {	// now assign pure var correctly, check for satisfaction, etc.
			assignments[pureVar] = pureVBC[pureVar][0] / abs(pureVBC[pureVar][0]);
			satisfyClauses(pureVar, &pureCopy, &pureSats, &assignments, &pureVBC);
			return solve(directions, variables, &pureCopy, &pureSats, assignments, &pureVBC);
		}
	}

	// there is guaranteed to be a 0 in assignments, because if there was not we would have retunred from allSat == TRUE
	int nextVarIndex = std::distance(assignments.begin(), std::find(assignments.begin(), assignments.end(), 0));
	// trying false
	assignments[nextVarIndex] = -1;
	std::vector<std::vector<int>> falseCopy(*clauses);	// a copy so we don't have to revert changes
	std::vector<int> falseSats(*clauseSats);
	std::vector<std::vector<int>> falseVBC(*varsByClause);	// a copy so we don't have to revert changes

	// satisfy clauses
	satisfyClauses(nextVarIndex, &falseCopy, &falseSats, &assignments, &falseVBC);
	double probSatFalse = solve(directions, variables, &falseCopy, &falseSats, assignments, &falseVBC);

	// trying true
	assignments[nextVarIndex] = 1;
	std::vector<std::vector<int>> trueCopy(*clauses);	// a copy so we don't have to revert changes
	std::vector<int> trueSats(*clauseSats);
	std::vector<std::vector<int>> trueVBC(*varsByClause);	// a copy so we don't have to revert changes

	// satisfy clauses
	satisfyClauses(nextVarIndex, &trueCopy, &trueSats, &assignments, &trueVBC);
	double probSatTrue = solve(directions, variables, &trueCopy, &trueSats, assignments, &trueVBC);

	if (variables->at(nextVarIndex) == -1) { 	// v is a choice variable
		return std::max(probSatFalse, probSatTrue);
	}
	// v is a chance variable
	return probSatTrue * variables->at(nextVarIndex) + probSatFalse * (1 - variables->at(nextVarIndex));
}

void satisfyClauses(int varIndex, std::vector<std::vector<int>>* clauses, std::vector<int>* sats, std::vector<int>* assignments, std::vector<std::vector<int>>* varsByClause)
{
	for (int c = 0; c < (signed int)clauses->size(); c++) {
		if (sats->at(c) == 1) {
			continue;
		}

		for (unsigned int l = 0; l < clauses->at(c).size(); l++) {
			if (clauses->at(c).at(l) == (varIndex + 1) * assignments->at(varIndex)) {
				sats->at(c) = 1;
				// now that a clause is satisfied, erase that clause from all the variable's list of active clauses
				for (unsigned int v = 0; v < varsByClause->size(); v++) {
					std::vector<int>::iterator it = std::find(varsByClause->at(v).begin(), varsByClause->at(v).end(), c * assignments->at(varIndex));

					if (it != varsByClause->at(v).end()) {
						varsByClause->at(v).erase(it);
					}
				}
			}
			else if (clauses->at(c).at(l) == (varIndex + 1) * assignments->at(varIndex) * -1) {
				varsByClause->at(varIndex).erase(std::find(varsByClause->at(varIndex).begin(), varsByClause->at(varIndex).end(), c * assignments->at(varIndex) * -1));

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