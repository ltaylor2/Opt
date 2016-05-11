README.txt

Liam Taylor and Henry Daniels-Koch
Optimization and Uncertainty
May 11, 2016
Assignment5 -- SSAT Algorithm

The program initializes a markov decision process test problem (three terminal states and
one state variable, a key item that must be held for the positive terminal state),
and, given input parameters, solves that problem with either Policy Iteration or
Value Iteration techiques. Upon finding a solution, the program will output the 
calculated utility for each state and the current optimal policy, in addition to
run-time and iteration information as well as a record of parameter settings.

IMPORTANT NOTES: The full contents of nr3.h and ludcmp.h are taken from the Numerical Recipes library by S. Teukolsky and W Press. LT and HDK did not alter the contents of these files in any way. The printResults and initMDP functions were adapted from java code by SM

This program can be compiled with the command:
g++ -std=c++11 main.cpp -o mdp

The program can then be run with:
./mdp [discount] [epsilon] [keyLoss] [posTerm] [negTerm] [stepCost] [solutionType]

PARAMETERS:
	discount (NUMERIC): The factor that decreases reward impact across steps in both policy
				and value iteration

	epsilon (NUMERIC): The error factor that determines the stop condition range during
				value iteration

	keyLoss (NUMERIC): The chance of losing the key item in the key loss state square 
				of the specific MDP problem initialized here
	
	posTerm (NUMERIC): The positive terminal reward given for a terminal state
				requiring a key item in this MDP problem

	negTerm (NUMERIC): The negative terminal reward given for the two negative 
				states near the key retrieval state in this MDP problem

	stepCost (NUMERIC): The cost of taking each action in this MDP (influencing the rewards)

	solutionType (CHAR): p OR v, where
				p = Policy iteration 
				v = Value iteration
