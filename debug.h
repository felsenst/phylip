/* Version 4.0. (c) Copyright 2012-2013 by the University of Washington.
   Written by Bob Giansiracusa for debugging purposes.
   This file contains declarations used only for debugging, not to be included in distributed systems. */


#ifndef _DEBUG_H_
#define _DEBUG_H_


// These two functions are NOT called in source-code.  They are available for calling manually in GDB.
//
void seenode(node * pq);
void seetree(tree * pt);


// This function IS called in source-code, when debugging.
//
// Prints debugging information whether tree is valid or erroneous.
// Checks for erroneous tree data structures and (if so found) aborts.
//
void treechecker(tree * pt, const char * filename, const unsigned int linenum, const char * funcname);


#endif /* _DEBUG_H_ */


// End.
