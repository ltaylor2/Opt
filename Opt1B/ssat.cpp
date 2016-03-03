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

//Specifies which solution the user would like
enum SolutionType { naive, unit, pure, both, hOne, hTwo, hThree };

// Prototype functions

//
//Reads in a file with an ssat problem and fills the vector of variables, vector of clauses, and vector of variables that
//each contain a vector of the clause #s they appear in
//
int readSSATFile(std::string fileName, std::vector<double>*, std::vector<std::vector<int>>*, std::vector<std::vector<int>>*);

//Solves the SSAT problem
double solve(SolutionType, std::vector<double>*, std::vector<std::vector<int>>*, std::vector<int>*, std::vector<int>, std::vector<std::vector<int>>*);

//Makes sure all clauses are satisfied
void satisfyClauses(int, std::vector<std::vector<int>>*, std::vector<int>*, std::vector<int>*, std::vector<std::vector<int>>*);

//Main
int main(int argc, char* argv[])
{
  //
  //Reads in the command line arguments
  //User did not put in exactly 3 arguments
  if (argc != 3) {
    std::cout << "Invalid Arguments (" << argc << ") Exiting." << std::endl;
    return 1;
  }

  //
  //Determines which solution type the user selected
  //
  SolutionType directions;

  //User selected naive solution
  if (std::string(argv[1]).compare("n") == 0)
    directions = SolutionType::naive;

  //User selected solution with unit clause propogation
  else if (std::string(argv[1]).compare("u") == 0)
    directions = SolutionType::unit;

  //User selected solution with pure variable elimination
  else if (std::string(argv[1]).compare("p") == 0)
    directions = SolutionType::pure;

  //User selected solution with both unit clause propogation and pure variable elimination
  else if (std::string(argv[1]).compare("b") == 0)
    directions = SolutionType::both;

  //User did not put in a valid solution type
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
  double solutionProb = solve(directions, &variables, &clauses, &clauseSats, assignments, &varsByClause);
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

  //Line a line in the file
  std::string line;

  //Iteratue until the loop reaches the end of the file
  while(getline(file, line)) {

    //Loop has reached the "variables" line in the file
    if (line.compare("variables") == 0) {

      //Skip a line
      getline(file, line);

      //Iterate until all the variables have been read in
      while(line.compare("") != 0) {

	//Convert the line into a string stream object
	std::stringstream ss(line);

	//Name of Variable
	int varName = 0;

	//Probabiltiy of variable being true
	double varValue = 0.0;

	//Assign values to the name of the variable and its probability from the file
	ss >> varName >> varValue;

	//Fill variable vector
	variables->push_back(varValue);

	//Fill vector of variables that each contain a vector of the clause #s they appear in
	varsByClause->push_back(std::vector<int>());

	//Skip a line
	getline(file, line);
      }
    }

    //Loop has reached the "variables" line in the file
    if (line.compare("clauses") == 0) {

      //Skip line
      getline(file, line);

      //Iterate until all the clauses have been read in
      while(line.compare("") != 0) {

	//Create a clause vector for each line
	std::vector<int> clause;

	//Create stringstream object
	std::stringstream ss(line);

	//Create int to hold each variable within the clause
	int literal = variables->size() + 1;
	ss >> literal;

	//Iterate until the end of the clause
	while (literal != 0) {

	  //Fill the clause vector
	  clause.push_back(literal);

	  //Create an int to hold the assignment of each variable
	  int litSign = 1;

	  //Variable is false
	  if (literal < 0)
	    litSign = -1;

	  //Fill vector of variables with the clause #s they appear in???
	  varsByClause->at(abs(literal) - 1).push_back(clauses->size() * litSign);
	  ss >> literal;
	}

	//Fill clauses vector with each clause
	clauses->push_back(clause);

	//Skip a line
	getline(file,line);
      }
    }
  }

  return 0;
}

//Solve SSAT problem
double solve(SolutionType directions,
	     std::vector<double>* variables,
	     std::vector<std::vector<int>>* clauses,
	     std::vector<int>* clauseSats,
	     std::vector<int> assignments,
	     std::vector<std::vector<int>>* varsByClause)
{

  
  //True if all clauses are satisfied
  bool allSat = true;

  //Loop through clause satisfactions, checking for fully satisfied or any unsatisfied
  
  for (unsigned int i = 0; i < clauseSats->size(); i++) {

    //A clause can never be satisfied and thus the entire SSAT problem cannot be satisfied
    if (clauseSats->at(i) == -1)	        
      return 0.0;

    //Non-satisfied (but not unsatisfied) clause
    else if (clauseSats->at(i) == 0) {  
      allSat = false;
    }
  }
  //All clauses are satisfied 
  if (allSat)	       
    return 1.0;


  //User wants solution to execute unit clause propogation
  if (directions == SolutionType::unit || directions == SolutionType::both) {

    //Create a vector that is a copy of all the clauses
    std::vector<std::vector<int>> unitCopy(*clauses);

    //Create a vector that is a copy of all the clause satisfaction values
    std::vector<int> unitSats(*clauseSats);

    //Create a vector that is a copy of the vector of variables that each contain a vector of the clause #s they appear in
    std::vector<std::vector<int>> unitVBC(*varsByClause);

    int unitVar	= 0;

    //Loop through all the clauses to look for unit clauses
    for (unsigned int c = 0; c < unitCopy.size(); c++) {

      //Clause is a unit clause (size 1) and is not satisfied
      if (unitCopy[c].size() == 1 && unitSats[c] != 1) {

	//Assign unitVar the value of the unit clause variable
	unitVar = unitCopy[c][0];
	break;
      }
    }
    
    //Unit clause found
    if (unitVar != 0) {

      //Fill assignments vector with variable assignment (-1 or 1)
      assignments[abs(unitVar)-1] = unitVar / abs(unitVar);

      //With new state, go through all clauses and determine if they are satisfied or can be shortened
      satisfyClauses(abs(unitVar)-1, &unitCopy, &unitSats, &assignments, &unitVBC);

      //Solve SSAT with unit clause propogation already conducted
      double probSatUnit = solve(directions, variables, &unitCopy, &unitSats, assignments, &unitVBC);

      //Unit variable is a choice variable
      if (variables->at(abs(unitVar)-1) == -1)
	return probSatUnit;
      else {

	//Unit variable is a probabilistic  variable that must be true
	if (assignments[abs(unitVar)-1] == 1)
	  return probSatUnit * variables->at(abs(unitVar)-1);

	//Unit variable is a probabilistic variable that must be false
	else
	  return probSatUnit * (1 - variables->at(abs(unitVar)-1));
      }
    }
  }

  //User wants to eliminate pure variables before solving SSAT
  if (directions == SolutionType::pure || directions == SolutionType::both) {

    //Create a vector that is a copy of all the clauses
    std::vector<std::vector<int>> pureCopy(*clauses);

    //Create a vector that is a copy of all the clause satisfaction values
    std::vector<int> pureSats(*clauseSats);

    //Create a vector that is a copy of the vector of variables that each contain a vector of the clause #s they appear in
    std::vector<std::vector<int>> pureVBC(*varsByClause);

    //Initialize an unassigned pure variable
    int pureVar = -1;

    //Loop through all the variables in pureVBC
    for (int l = 0; l < (signed int)pureVBC.size(); l++) {

      //Variable is a chance variable or already assigned
      if (variables->at(l) != -1 || assignments[l] != 0)	// only looking at unassigned choice variables
	continue;

      //Loop through each instance of each variable in its respective clauses
      for (unsigned int c = 1; c < varsByClause->at(l).size(); c++) {

	//Assign the pure variable
	pureVar = l;

	//Variable is not pure
	if (!((pureVBC[l][c-1] < 0 && pureVBC[l][c] < 0) || (pureVBC[l][c-1] > 0 && pureVBC[l][c] > 0)) || pureVBC[l].size() <= 1) {
	  pureVar = -1;

	  //Move on to next variable
	  break;
	}
      }

      //Variable is pure
      if (pureVar == l) 
	break;
    }

    //Found pure variable
    if (pureVar != -1) {	// now assign pure var correctly, check for satisfaction, etc.

      //Assign pure variable final value in assignments vector
      assignments[pureVar] = pureVBC[pureVar][0] / abs(pureVBC[pureVar][0]);

      //Check if all 3clauses are satisfied
      satisfyClauses(pureVar, &pureCopy, &pureSats, &assignments, &pureVBC);
      
      return solve(directions, variables, &pureCopy, &pureSats, assignments, &pureVBC);
    }
  }

  //There is guaranteed to be a 0 in assignments, because if there was not we would have retunred from allSat == TRUE
  
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

//Checks if all clauses are satisfied
void satisfyClauses(int varIndex, std::vector<std::vector<int>>* clauses, std::vector<int>* sats, std::vector<int>* assignments, std::vector<std::vector<int>>* varsByClause)
{
  //Loops through all of the clauses
  for (int c = 0; c < (signed int)clauses->size(); c++) {

    //Clause is satisfied
    if (sats->at(c) == 1)
      continue;

    //Clause might still be satisfied, Let's check

    //Loop through each clause
    for (unsigned int l = 0; l < clauses->at(c).size(); l++) {

      //Clause is satisfied, mark this in sats vector
      if (clauses->at(c)[l] == (varIndex + 1) * assignments->at(varIndex)) {
	sats->at(c) = 1;

	// Now that a clause is satisfied, erase that clause from the variable's list of active clauses
	for (unsigned int v = 0; v < varsByClause->size(); v++) {

	  //I'm confused
	  std::vector<int>::iterator it = std::find(varsByClause->at(v).begin(), varsByClause->at(v).end(), c * assignments->at(varIndex));

	  // Found variable from vector of 
	  if (it != varsByClause->at(v).end()) {
	    varsByClause->at(v).erase(it);
	  }
	}
      }
      else if (clauses->at(c)[l] == (varIndex + 1) * assignments->at(varIndex) * -1) {
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
