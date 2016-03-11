/* 
 * File:   ssat.cpp
 * Authosr: Liam Taylor and Henry Daniels-Koch
 * 
 * Created on February 27, 2016
 */
//ssat.cpp.
//      Assignment1b.c++-Source code for Planning as Satisfiability.
//FUNCTIONAL DESCRIPTION
//      This program in C++ fulfills the requirements of Assignment 1b in
//      CSCI3460 Spring 2016
//
//      1. 
//Notice
//      Copyright (C) Feburary 27, 2016 to March 9, 2016
//      Liam Taylor and Henry Daniels-Koch All Rights Reserved.
//


#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <climits>

//Specifies which solution the user would like
enum SolutionType { naive, unit, pure, both, hOne, hTwo, hThree };

//
//Reads in a file with an ssat problem and fills the vector of variables, vector of clauses, and vector of variables that
//each contain a vector of the clause #s they appear in
//
int readSSATFile(std::string fileName, std::vector<double>*, std::vector<std::vector<int>>*, std::vector<std::vector<int>>*);

//Solves the SSAT problem
double solve(SolutionType, std::vector<double>*, std::vector<std::vector<int>>, std::vector<int>, std::vector<int>*, std::vector<std::vector<int>>);

//Makes sure all clauses are satisfied
void satisfyClauses(int, std::vector<std::vector<int>>*, std::vector<int>*, std::vector<int>*, std::vector<std::vector<int>>*);


//Main
int main(int argc, char* argv[])
{
    //
    //Reads in the command line arguments
    //User did not put in exactly 3 arguments
    if (argc != 3) {
		std::cout << "Invalid Arguments (" << argc << "). Need [directions] [filetype] -- Exiting." << std::endl;
		return 1;
    }

    SolutionType directions;

    if (std::string(argv[1]).compare("n") == 0)
		directions = SolutionType::naive;
    else if (std::string(argv[1]).compare("u") == 0)
		directions = SolutionType::unit;
    else if (std::string(argv[1]).compare("p") == 0)
		directions = SolutionType::pure;
    else if (std::string(argv[1]).compare("b") == 0)
		directions = SolutionType::both;
	else if (std::string(argv[1]).compare("1") == 0)
		directions = SolutionType::hOne;
	else if (std::string(argv[1]).compare("2") == 0)
		directions = SolutionType::hTwo;
	else if (std::string(argv[1]).compare("3") == 0)
		directions = SolutionType::hThree;
    else {
		std::cout << "---" << argv[1] << "---" << std::endl;
		std::cout << "Incorrect solving directions. Exiting." << std::endl;
		return 1;
    }

    //Get filename of user input ssat file
    std::string fileName = std::string(argv[2]);

    //
    //Initialize data structures to hold variables and clauses
    //
	
    //Vector of all variables in order with their probabilities
    std::vector<double> variables;

    //Vector of variable assignments where:
    //-1 is false
    //0 is unassigned
    //1 is true
    std::vector<int> assignments;

    //Vector of all clauses (where each clause is a vector)
    std::vector<std::vector<int>> clauses;

    //Vector that shows whether each clause is satisfied where:
    //-1 is unsatisfied
    //0 un unassigned
    //1 is satisfied
    std::vector<int> clauseSats;

    //Vector of variables that each contain a vector of the clause #s they appear in
    std::vector<std::vector<int>> varsByClause;

    //Read file in and assign values to variables and clauses
    //If file could not be opened, return 1
    if (readSSATFile(fileName, &variables, &clauses, &varsByClause) == 1) {
		return 1;
    }

    //Fill assignments such that each variable has an unassigned value
    for (unsigned int i = 0; i < variables.size(); i++) {
		assignments.push_back(0);
    }
	
    //Fill clauseSats such that each clause is unassigned a satisfaction value
    for (unsigned int i = 0; i < clauses.size(); i++) {
		clauseSats.push_back(0);
    }

    //Start solving the SSAT Problem and time it
    std::cout << "Beginning to solve!" << std::endl;
    clock_t start = clock();
    double solutionProb = solve(directions, &variables, clauses, clauseSats, &assignments, varsByClause);
    clock_t end = clock();

    double solveTime  = (double)(end-start) / CLOCKS_PER_SEC;
    std::cout << "Solution is: " << solutionProb << " (found in " << solveTime << " seconds)" << std::endl;

    return 0;
}

//
//Reads in a file with an ssat problem and fills the vector of variables, array of clauses, and
// vector of variables that each contain a vector of the clause #s they appear in
//
int readSSATFile(std::string fileName,
		 std::vector<double>* variables, 
		 std::vector<std::vector<int>>* clauses,
		 std::vector<std::vector<int>>* varsByClause)
{
    std::cout << "Reading in " << fileName << std::endl;

    //Does this construct an ifstream object?
    std::ifstream file(fileName);

    //File could not be opened
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

		//Loop has reached the "variables" line in the file
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

				    varsByClause->at(abs(literal) - 1).push_back(((int)clauses->size() + 1) * litSign);
				    ss >> literal;
				}

				clauses->push_back(clause);

				getline(file,line);
		    }
		}
    }

    return 0;
}

//Solve SSAT problem
double solve(SolutionType directions,
	     std::vector<double>* variables,
	     std::vector<std::vector<int>> clauses,
	     std::vector<int> clauseSats,
	     std::vector<int>* assignments,
	     std::vector<std::vector<int>> varsByClause)
{

    bool allSat = true;

    for (unsigned int i = 0; i < clauseSats.size(); i++) {
		if (clauseSats[i] == -1)	        
		    return 0.0;
		else if (clauseSats[i] == 0)
		    allSat = false;
    }

    if (allSat)	       
		return 1.0;

    //User wants solution to execute unit clause propogation
    if (directions == SolutionType::unit || directions == SolutionType::both
    	 || directions == SolutionType::hOne || directions == SolutionType::hTwo || directions == SolutionType::hThree) {
		int unitVar	= 0;

		for (unsigned int c = 0; c < clauses.size(); c++) {
		    if (clauses[c].size() == 1 && clauseSats[c] == 0) {
				unitVar = clauses[c][0];
				break;
		    }
		}
	    
		//Unit clause found
		if (unitVar != 0) {

		    assignments->at(abs(unitVar)-1) = unitVar / abs(unitVar);

		    satisfyClauses(abs(unitVar)-1, &clauses, &clauseSats, assignments, &varsByClause);
		    double probSatUnit = solve(directions, variables, clauses, clauseSats, assignments, varsByClause);

		    if (variables->at(abs(unitVar)-1) == -1)
				return probSatUnit;
		    else {
				if (assignments->at(abs(unitVar)-1) == 1)
				    return probSatUnit * variables->at(abs(unitVar)-1);
				else
				    return probSatUnit * (1 - variables->at(abs(unitVar)-1));
		    }
		}
    }

    //User wants to eliminate pure variables before solving SSAT
    if (directions == SolutionType::pure || directions == SolutionType::both
    	 || directions == SolutionType::hOne || directions == SolutionType::hTwo || directions == SolutionType::hThree) {

		int pureVar = -1;

		for (int l = 0; l < (signed int)varsByClause.size(); l++) {

		    if (variables->at(l) != -1 || assignments->at(l) != 0)	// only looking at unassigned choice variables
				continue;

		    for (unsigned int c = 1; c < varsByClause[l].size(); c++) {
				pureVar = l;

				// Variable is not pure
				if (varsByClause[l].size() < 1 || !((varsByClause[l][c-1] < 0 && varsByClause[l][c] < 0) || (varsByClause[l][c-1] > 0 && varsByClause[l][c] > 0))) {
				    pureVar = -1;
				    break;
				}
		    }

		    //Variable is pure
		    if (pureVar == l) 
				break;
		}

		//Found pure variable
		if (pureVar != -1) {	// now assign pure var correctly, check for satisfaction, etc.
		    assignments->at(pureVar) = varsByClause[pureVar][0] / abs(varsByClause[pureVar][0]);

		    satisfyClauses(pureVar, &clauses, &clauseSats, assignments, &varsByClause);
		    return solve(directions, variables, clauses, clauseSats, assignments, varsByClause);
		}
    }

    //There is guaranteed to be a 0 in assignments, because if there was not we would have retunred from allSat == TRUE
    int nextVarIndex = std::distance(assignments->begin(), std::find(assignments->begin(), assignments->end(), 0));
    
    if (directions == SolutionType::hOne) {
    	int minClauseLength = INT_MAX;

   		// decide the current block
   		bool choiceVar = (variables->at(nextVarIndex) == -1);

   		// now run through every variable in the current block
   		bool currBlock = true;
   		for (int i = nextVarIndex; currBlock; i++) {

   			// make sure we're still in the block
   			if (i+1 == varsByClause.size())
   				currBlock = false;
   			else if (choiceVar) {
   				if (variables->at(i+1) != -1)
   					currBlock = false;
   			} else {
   				if (variables->at(i+1) == -1)
   					currBlock = false;
   			}

   			if (assignments->at(i) != 0)
   				continue;

   			for (unsigned int c = 0; c < varsByClause[i].size(); c++) {
   				int currLength = clauses[abs(varsByClause[i][c]) - 1].size();

   				if (currLength < minClauseLength) {
   					minClauseLength = currLength;
   					nextVarIndex = i;
   				}
   			}
   		}
    }
    if (directions == SolutionType::hTwo) {
    	int maxCount = 0;

    	bool choiceVar = (variables->at(nextVarIndex) == -1);

    	bool currBlock = true;

    	for (int i = nextVarIndex; currBlock; i++) {
    		int currCount = 0;

    		if (i+1 == varsByClause.size()
    			|| (choiceVar && (variables->at(i+1) != -1))
    			|| (!choiceVar && (variables->at(i+1) == 1)))
    			currBlock = false;

    		if (assignments->at(i) != 0)
    			continue;

    		currCount = varsByClause[nextVarIndex].size();

    		if (currCount > maxCount) {
    			maxCount = currCount;
    			nextVarIndex = i;
    		}
    	}
    }
    
    if (directions == SolutionType::hThree) {
    	double maxCount = 0.0;

   		// decide the current block
   		bool choiceVar = (variables->at(nextVarIndex) == -1);

   		// now run through every variable in the current block
   		bool currBlock = true;

    	for (int i = nextVarIndex; currBlock; i++) {
    		double currPosCount = 0.0;
    		double currNegCount = 0.0;

    		// make sure we're still in the block
    		if (i+1 == varsByClause.size())
    			currBlock = false;
    		else if (choiceVar) {
    			if (variables->at(i+1) != -1)
    				currBlock = false;
    		} else {
    			if (variables->at(i+1) == -1)
    				currBlock = false;
    		}

    		if (assignments->at(i) != 0)
    			continue;

    		for (unsigned int n = 0; n < varsByClause[i].size(); n++) {
    			if ((varsByClause[i][n]  / abs(varsByClause[i][n])) == 1)
    				currPosCount++;
    			else
    				currNegCount++;
    		}

    		if (variables->at(i) != -1) {
    			currPosCount *= variables->at(i);
    			currNegCount *= (1 - variables->at(i));
    		}

    		if (assignments->at(i) == 0 && std::max(currPosCount, currNegCount) > maxCount ) {
    			maxCount = std::max(currPosCount, currNegCount);
    			nextVarIndex = i;
    		}
    	}
    }

    // trying false
    assignments->at(nextVarIndex) = -1;

    std::vector<std::vector<int>> falseClauses(clauses);
    std::vector<int> falseSats(clauseSats);
    std::vector<int> falseAssignments(*assignments);
    std::vector<std::vector<int>> falseVBC(varsByClause);

    satisfyClauses(nextVarIndex, &falseClauses, &falseSats, &falseAssignments, &falseVBC);
    double probSatFalse = solve(directions, variables, falseClauses, falseSats, &falseAssignments, falseVBC);

    // trying true
    assignments->at(nextVarIndex) = 1;

    // satisfy clauses
    satisfyClauses(nextVarIndex, &clauses, &clauseSats, assignments, &varsByClause);
    double probSatTrue = solve(directions, variables, clauses, clauseSats, assignments, varsByClause);

    if (variables->at(nextVarIndex) == -1) { 	// v is a choice variable
		return std::max(probSatFalse, probSatTrue);
    }
    
    // v is a chance variable
    return probSatTrue * variables->at(nextVarIndex) + probSatFalse * (1 - variables->at(nextVarIndex));
}

//Checks if all clauses are satisfied
void satisfyClauses(int varIndex, std::vector<std::vector<int>>* clauses, std::vector<int>* sats, std::vector<int>* assignments, std::vector<std::vector<int>>* varsByClause)
{
    for (int c = 0; c < (signed int)clauses->size(); c++) {
		if (sats->at(c) == 1)
		    continue;

		for (unsigned int l = 0; l < clauses->at(c).size(); l++) {
		    if (clauses->at(c)[l] == (varIndex + 1) * assignments->at(varIndex)) {
				sats->at(c) = 1;

				for (unsigned int v = 0; v < varsByClause->size(); v++) {

				    std::vector<int>::iterator it = std::find(varsByClause->at(v).begin(), varsByClause->at(v).end(), (c + 1)  * assignments->at(varIndex));

				    if (it != varsByClause->at(v).end())
						varsByClause->at(v).erase(it);
				}
		    }
		    else if (clauses->at(c)[l] == (varIndex + 1) * assignments->at(varIndex) * -1) {
				varsByClause->at(varIndex).erase(std::find(varsByClause->at(varIndex).begin(), varsByClause->at(varIndex).end(), (c + 1) * assignments->at(varIndex) * -1));

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
