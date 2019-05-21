/* Version 4.0. (c) Copyright 2012-2013 by the University of Washington.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#include "phylip.h"

boolean usejtt;
boolean usepmb;

double  * probcat;
double  * rate;
double  * rrate;

long *  enterorder;

long    categs;
long    njumble;
long    rcategs;
long    sites;

steptr  aliasweight;

tree *curtree, *bestree, *bestree2, *priortree;

/* Variables introduced to allow for protein probability calculations   */
long   max_num_sibs;            /* maximum number of siblings used in a */
                                /* nuview calculation.  determines size */
                                /* final size of pmatrices              */
double *eigmat;                 /* eig matrix variable                  */
double **probmat;               /* prob matrix variable                 */
double ****dpmatrix;            /* derivative of pmatrix                */
double ****ddpmatrix;           /* derivative of xpmatrix               */
double *****pmatrices;          /* matrix of probabilities of protein   */
                                /* conversion.  The 5 subscripts refer  */
                                /* to sibs, rcategs, categs, final and  */
                                /* initial states, respectively.        */
double freqaa[20];              /* amino acid frequencies               */


extern double jtteigmat[];
extern double jttprobmat[20][20];
extern double pmbeigmat[20];
extern double pmbprobmat[20][20];
extern double pameigmat[];
extern double pamprobmat[20][20];


void    prom_alloc_pmatrix(long sib);
void    prom_allocrest(void);
void    prom_clean_up(void);
void    prom_free_all_x(long nonodes, pointarray treenode);
void    prom_free_pmatrix(long sib);
void    prom_init_protmats(void);
void    prom_makeweights(void);
void    prom_reallocsites(void);


// End.
