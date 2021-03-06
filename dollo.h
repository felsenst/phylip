/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


/*
   dollo.h: included in dollop, dolmove & dolpenny
*/


#ifndef OLDC
/* function prototypes */
void   correct(tree *, node *, long, boolean, bitptr);
void   fillin(node *);
void   postorder(node *);
void   count(long *, bitptr, steptr, steptr);
void   filltrav(node *);
void   hyprint(struct htrav_vars *, boolean *, bitptr, Char *);
void   hyptrav(tree *, node *, boolean *, bitptr, long, boolean,
                 Char *pointptr, gbit *, bitptr, bitptr);
void   hypstates(tree *, long, boolean, Char *,
                   node *, gbit *, bitptr, bitptr);
void   drawline(long, double, node *);
void   printree(double, boolean, node *);
void   writesteps(boolean, boolean, steptr);
/* function prototypes */
#endif


// End.
