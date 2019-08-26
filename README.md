# opt-simplex-tableau

A Julia script for extracting the optimal simplex tableau for a given integer linear program in SEF. Currently only works for .mps files and uses IBM's CPLEX Solver. 

The problem is assumed to be in SEF. The tableau that is returned is for the original instance with slack variables (added by the solver). That is, given an instance in the form  
 	A x = b  
 	l <= x <= u  
 The optimal tableau for the following augmented instance is returned:  
 	[A | I] [x // s] <= b  
 	l <= x <= u  
 	s >= 0  
  
Command line arguments: FILE_NAME
