/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Dan Fineman, Joseph Felsenstein, Mike Palczewski, Hisashi Horino,
   Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee
   is charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "cons.h"


long output_scheme, trees2;

extern long tree_pairing;
extern double *lengths;
extern hashtype hashp;

/* The following extern's refer to things declared in cons.c */

typedef enum { SYMMETRIC, BSD, RF } distance_type;

distance_type  dtype;
extern Char outfilename[FNMLNGTH], intreename[FNMLNGTH], intree2name[FNMLNGTH], outtreename[FNMLNGTH];

extern long numopts, outgrno, col;
extern long maxgrp;               /* max. no. of groups in all trees found  */

extern boolean trout, firsttree, noroot, outgropt, didreroot, prntsets, progress, treeprint, goteof;
extern group_type **grouping, **grping2, **group2;/* to store groups found  */
extern long **order, **order2, lasti;
extern group_type *fullset;
extern long tipy;

extern double **timesseen, **tmseen2, **times2;
extern double trweight, ntrees;

tree * treeForNodes;    // treedist doesn't use trees quite the way other
                        // programs do -- it stores info with tree nodes,
                        // though, so we use this tree to create them

#ifndef OLDC
/* function prototpes */
void   assign_tree(group_type **, pattern_elm ***, long, long *);
boolean group_is_null(group_type **, long);
void   compute_distances(pattern_elm ***, long, long);
void   free_patterns(pattern_elm ***, long);
void   produce_square_matrix(long, long *);
void   produce_full_matrix(long, long, long *);
void   output_submenu(void);
void   pairing_submenu(void);

void   read_second_file(pattern_elm ***, long, long);
void   getoptions(void);
void   assign_lengths(double **lengths, pattern_elm ***pattern_array, long tree_index);
void   print_header(long trees_in_1, long trees_in_2);
void   output_distances(long trees_in_1, long trees_in_2);
void   output_long_distance(long diffl, long tree1, long tree2, long trees_in_1, long trees_in_2);
void output_matrix_long(long diffl, long tree1, long tree2, long trees_in_1, long trees_in_2);
void output_matrix_double(double diffl, long tree1, long tree2, long trees_in_1, long trees_in_2);
void   output_double_distance(double diffd, long tree1, long tree2, long trees_in_1, long trees_in_2);
long   symetric_diff(group_type **tree1, group_type **tree2, long ntree1, long ntree2, long patternsz1, long patternsz2);
double bsd_tree_diff(group_type **tree1, group_type **tree2, long ntree1, long ntree2,
                     double* lengths1, double *lengths2, long patternsz1, long patternsz2);
void   tree_diff(group_type **tree1, group_type **tree2, double *lengths1, double* lengths2, long patternsz1,
                 long patternsz2, long ntree1, long ntree2, long trees_in_1, long trees_in_2);
void   print_line_heading(long tree);
int    get_num_columns(void);
void   print_matrix_heading(long tree, long maxtree);
void   doinit(void);
void   treedistrun(void);
void   treedist(char *intreename, char *intree2name, char *outfilename, char *outfileopt, char *DistType,
                char *DistKind, int Rooted, char *OutputKind, int PrintInd);
/* function prototpes */
#endif


void assign_lengths(double **lengths, pattern_elm ***pattern_array, long tree_index)
{
  *lengths = pattern_array[0][tree_index]->length;
}


void assign_tree(group_type **treeN, pattern_elm ***pattern_array, long tree_index, long *pattern_size)
{ /* set treeN to be the tree_index-th tree in pattern_elm */
  long i;

  for (i = 0 ; i < setsz ; i++) {
    treeN[i] = pattern_array[i][tree_index]->apattern;
  }
  *pattern_size = *pattern_array[0][tree_index]->patternsize;
}  /* assign_tree */


boolean group_is_null(group_type **treeN, long index)
{
  /* Check to see if a given index to a tree array points to an empty
     group */
  long i;

  for (i = 0 ; i < setsz ; i++)
    if (treeN[i][index] != (group_type) 0)
      return false;

  /* If we've gotten this far, then the index is to an empty group in
     the tree. */
  return true;
}  /* group_is_null */


double bsd_tree_diff(group_type **tree1, group_type **tree2,
                     long ntree1, long ntree2, double *lengths1,
                     double* lengths2, long patternsz1, long patternsz2)
{
  /* Compute the difference between 2 given trees. Return
     that value as a double. */

  long index1, index2;
  double return_value = 0;
  boolean match_found;
  long i;

  if (group_is_null (tree1, 0) || group_is_null (tree2, 0)) {
    printf ("Error computing tree difference between tree %ld and tree %ld.\n", ntree1, ntree2);
    exxit(-1);
  }

  for (index1 = 0; index1 < patternsz1; index1++) {
    if (!group_is_null (tree1, index1)) {
      if (lengths1[index1] == -1) {
        printf("\nError: tree %ld is missing a length from at least one branch\n", ntree1);
        printf("Did you instead intend to compute the symmetric difference?\n\n");
        exxit(-1);
      }
    }
  }

  for (index2 = 0; index2 < patternsz2; index2++) {
    if (!group_is_null (tree2, index2)) {
      if (lengths2[index2] == -1) {
        printf("Error: tree %ld is missing a length from at least one branch\n", ntree2);
        exxit(-1);
      }
    }
  }

  for (index1 = 0 ; index1 < patternsz1; index1++) {
    /* For every element in the first tree, see if there's
       a match to it in the second tree. */
    match_found = false;

    if (group_is_null (tree1, index1)) {
      /* When we've gone over all the elements in tree1, greater
         number of elements in tree2 will constitute that much more
         of a difference... */
      while (! group_is_null (tree2, index1)) {
        if (dtype == BSD)
          return_value+= pow(lengths1[index1], 2);
        else  /* RF case */
          return_value+= fabs(lengths1[index1]);
        index1++;
      }
      break;
    }

    for (index2 = 0 ; index2 < patternsz2 ; index2++) {
      /* For every element in the second tree, see if any match
         the current element in the first tree. */
      if (group_is_null (tree2, index2)) {
        /* When we've gone over all the elements in tree2 */
        match_found = false;
        break;
      } else {
        /* Tentatively set match_found; will be changed later if
           neccessary. . . */
        match_found = true;

        for (i = 0 ; i < setsz ; i++) {
          /* See if we've got a match, */
          if (tree1[i][index1] != tree2[i][index2])
            match_found = false;
        }

        if (match_found == true) {
          break;
        }
      }
    }

    if (match_found == false) {
      if (dtype == BSD)
        return_value+= pow(lengths1[index1], 2);
      else  /* RF case */
        return_value+= fabs(lengths1[index1]);
    }
  }

  for (index2 = 0 ; index2 < patternsz2 ; index2++) {
    /* For every element in the second tree, see if there's
       a match to it in the first tree. */
    match_found = false;
    if (group_is_null (tree2, index2)) {
      /* When we've gone over all the elements in tree2, greater
         number of elements in tree1 will constitute that much more
         of a difference... */

      while (! group_is_null (tree1, index2)) {
        if (dtype == BSD)
          return_value+= pow(lengths2[index2], 2);
        else  /* RF case */
          return_value+= fabs(lengths2[index2]);
        index2++;
      }
      break;
    }

    for (index1 = 0 ; index1 < patternsz1 ; index1++) {
      /* For every element in the first tree, see if any match
         the current element in the second tree. */
      if (group_is_null (tree1, index1)) {
        /* When we've gone over all the elements in tree2 */
        match_found = false;
        break;
      } else {
        /* Tentatively set match_found; will be changed later if
           neccessary. . . */
        match_found = true;

        for (i = 0 ; i < setsz ; i++) {
          /* See if we've got a match, */
          if (tree2[i][index2] != tree1[i][index1])
            match_found = false;
        }

        if (match_found == true) {
          if (dtype == BSD)
            return_value += pow(lengths1[index1] - lengths2[index2], 2);
          else  /* RF case */
            return_value += fabs(lengths1[index1] - lengths2[index2]);
          break;
        }
      }
    }

    if (match_found == false) {
      if (dtype == BSD)
        return_value+= pow(lengths2[index2], 2);
      else  /* RF case */
        return_value+= fabs(lengths2[index2]);
    }
  }
  if (return_value > 0.0) {
    if (dtype == BSD)
      return_value = sqrt(return_value);
    else  /* RF case */
      return_value = fabs(return_value);
  }
  else
    return_value = 0.0;
  return return_value;
}


long symetric_diff(group_type **tree1, group_type **tree2,
                   long ntree1, long ntree2, long patternsz1, long patternsz2)
{
  /* Compute the symmetric difference between 2 given trees. Return
     that value as a long. */

  long index1, index2, return_value = 0;
  boolean match_found;
  long i;

  (void)ntree1;                         // RSGnote: Parameter never used.
  (void)ntree2;                         // RSGnote: Parameter never used.

  if (group_is_null (tree1, 0) || group_is_null (tree2, 0)) {
    printf ("Error computing tree difference.\n");
    return 0;
  }

  for (index1 = 0 ; index1 < patternsz1 ; index1++) {
    /* For every element in the first tree, see if there's
       a match to it in the second tree. */
    match_found = false;
    if (group_is_null (tree1, index1)) {
      /* When we've gone over all the elements in tree1, greater
         number of elements in tree2 will constitute that much more
         of a difference... */

      while (! group_is_null (tree2, index1)) {
        return_value++;
        index1++;
      }
      break;
    }

    for (index2 = 0 ; index2 < patternsz2 ; index2++) {
      /* For every element in the second tree, see if any match
         the current element in the first tree. */
      if (group_is_null (tree2, index2)) {
        /* When we've gone over all the elements in tree2 */
        match_found = false;
        break;
      } else {
        /* Tentatively set match_found; will be changed later if
           neccessary. . . */
        match_found = true;

        for (i = 0 ; i < setsz ; i++) {
          /* See if we've got a match, */
          if (tree1[i][index1] != tree2[i][index2])
            match_found = false;
        }

        if (match_found == true) {
          /* If the previous loop ran from 0 to setsz without setting
             match_found to false, */
          break;
        }
      }
    }

    if (match_found == false) {
      return_value++;
    }
  }
  return return_value;
}  /* symetric_diff */


void output_double_distance(double diffd, long tree1, long tree2, long trees_in_1, long trees_in_2)
{
  switch (tree_pairing) {
    case ADJACENT_PAIRS:
      if (output_scheme == VERBOSE ) {
        fprintf (outfile, "Trees %ld and %ld:    %e\n", tree1, tree2, diffd);
      } else if (output_scheme == SPARSE) {
        fprintf (outfile, "%ld %ld %e\n", tree1, tree2, diffd);
      }
      break;

    case ALL_IN_FIRST:
      if (output_scheme == VERBOSE) {
        fprintf (outfile, "Trees %ld and %ld:    %e\n", tree1, tree2, diffd);
      } else if (output_scheme == SPARSE) {
        fprintf (outfile, "%ld %ld %e\n", tree1, tree2, diffd );
      } else if ((output_scheme == FULL_MATRIX) ||
                 (output_scheme == COMPUTER_READABLE_MATRIX)) {
        output_matrix_double(diffd, tree1, tree2, trees_in_1, trees_in_2);
      }
      break;

    case CORR_IN_1_AND_2:
      if (output_scheme == VERBOSE) {
        fprintf (outfile, "Tree pair %ld:    %e\n", tree1, diffd);
      } else if (output_scheme == SPARSE) {
        fprintf (outfile, "%ld %e\n", tree1, diffd);
      }
      break;

    case ALL_IN_1_AND_2:
      if (output_scheme == VERBOSE )
        fprintf (outfile, "Trees %ld and %ld:    %e\n", tree1, tree2, diffd);
      else if (output_scheme == SPARSE)
        fprintf (outfile, "%ld %ld %e\n", tree1, tree2, diffd);
      else if ((output_scheme == FULL_MATRIX) || 
               (output_scheme == COMPUTER_READABLE_MATRIX)) {
        output_matrix_double(diffd, tree1, tree2, trees_in_1, trees_in_2);
      }
      break;
  }
} /* output_double_distance */


void print_matrix_heading(long tree, long maxtree)
{
  long i;

  if ( tree_pairing == ALL_IN_1_AND_2 ) {
    fprintf(outfile, "\n\nFirst\\  Second tree file:\n");
    fprintf(outfile, "tree  \\\n");
    fprintf(outfile, "file:  \\");
  } else fprintf(outfile, "\n\n      ");

  for ( i = tree ;  i <= maxtree ; i++ ) {
    if ( dtype == SYMMETRIC )
      fprintf(outfile, "%5ld ", i);
    else
      fprintf(outfile, "    %7ld ", i);
  }
  fprintf(outfile, "\n");
  if ( tree_pairing == ALL_IN_1_AND_2 )
    fprintf(outfile, "        \\");
  else
    fprintf(outfile, "      \\");
  for ( i = tree ;  i <= maxtree ; i++ ) {
    if ( dtype == SYMMETRIC )
      fprintf(outfile, "------");
    else fprintf(outfile, "------------");
  }
}


void print_line_heading(long tree)
{
  if ( tree_pairing == ALL_IN_1_AND_2 )
    fprintf(outfile, "\n%4ld    |", tree);
  else fprintf(outfile, "\n%5ld |", tree);

}


void output_matrix_double(double diffl, long tree1, long tree2, long trees_in_1, long trees_in_2)
{
  if ( tree1 == 1 && ((tree2 - 1) % get_num_columns() == 0 || tree2 == 1 ))
  {
    if ( (tree_pairing == ALL_IN_FIRST && tree2 + get_num_columns() - 1 < trees_in_1) ||
         (tree_pairing == ALL_IN_1_AND_2 && tree2 + get_num_columns() - 1 < trees_in_2))
    {
      if ( output_scheme == FULL_MATRIX )
        print_matrix_heading(tree2, tree2 + get_num_columns() - 1);
    }
    else
    {
      if ( tree_pairing == ALL_IN_FIRST)
        print_matrix_heading(tree2, trees_in_1);
      else
        print_matrix_heading(tree2, trees_in_2);
    }
  }
  if ( (tree2 - 1) % get_num_columns() == 0 || tree2 == 1)
  {
    print_line_heading(tree1);
  }
  fprintf(outfile, " %9g  ", diffl);
  if ((tree_pairing == ALL_IN_FIRST && tree1 == trees_in_1 && tree2 == trees_in_1) ||
      (tree_pairing == ALL_IN_1_AND_2 && tree1 == trees_in_1 && tree2 == trees_in_2))
    fprintf(outfile, "\n\n\n");
} /* output_matrix_double */


void output_matrix_long(long diffl, long tree1, long tree2, long trees_in_1, long trees_in_2)
{
  if ( tree1 == 1 && ((tree2 - 1) % get_num_columns() == 0 || tree2 == 1 ))
  {
    if ( (tree_pairing == ALL_IN_FIRST && tree2 + get_num_columns() - 1 < trees_in_1) ||
         (tree_pairing == ALL_IN_1_AND_2 && tree2 + get_num_columns() - 1 < trees_in_2))
    {
      print_matrix_heading(tree2, tree2 + get_num_columns() - 1);
    }
    else
    {
      if ( tree_pairing == ALL_IN_FIRST)
        print_matrix_heading(tree2, trees_in_1);
      else
        print_matrix_heading(tree2, trees_in_2);
    }
  }
  if ( (tree2 - 1) % get_num_columns() == 0 || tree2 == 1)
  {
    print_line_heading(tree1);
  }
  fprintf(outfile, "%4ld  ", diffl);
  if ((tree_pairing == ALL_IN_FIRST && tree1 == trees_in_1 && tree2 == trees_in_1) ||
      (tree_pairing == ALL_IN_1_AND_2 && tree1 == trees_in_1 && tree2 == trees_in_2))
    fprintf(outfile, "\n\n\n");
} /* output_matrix_long */


void output_long_distance(long diffl, long tree1, long tree2, long trees_in_1, long trees_in_2)
{
  switch (tree_pairing)
  {
    case ADJACENT_PAIRS:
      if (output_scheme == VERBOSE )
      {
        fprintf (outfile, "Trees %ld and %ld:    %ld\n", tree1, tree2, diffl);
      }
      else if (output_scheme == SPARSE)
      {
        fprintf (outfile, "%ld %ld %ld\n", tree1, tree2, diffl);
      }
      break;

    case ALL_IN_FIRST:
      if (output_scheme == VERBOSE)
      {
        fprintf (outfile, "Trees %ld and %ld:    %ld\n", tree1, tree2, diffl);
      }
      else if (output_scheme == SPARSE)
      {
        fprintf (outfile, "%ld %ld %ld\n", tree1, tree2, diffl );
      }
      else if ((output_scheme == FULL_MATRIX) || 
               (output_scheme == COMPUTER_READABLE_MATRIX))
      {
        output_matrix_long(diffl, tree1, tree2, trees_in_1, trees_in_2);
      }
      break;

    case CORR_IN_1_AND_2:
      if (output_scheme == VERBOSE)
      {
        fprintf (outfile, "Tree pair %ld:    %ld\n", tree1, diffl);
      }
      else if (output_scheme == SPARSE)
      {
        fprintf (outfile, "%ld %ld\n", tree1, diffl);
      }
      break;

    case ALL_IN_1_AND_2:
      if (output_scheme == VERBOSE)
        fprintf (outfile, "Trees %ld and %ld:    %ld\n", tree1, tree2, diffl);
      else if (output_scheme == SPARSE)
        fprintf (outfile, "%ld %ld %ld\n", tree1, tree2, diffl);
      else if ((output_scheme == FULL_MATRIX) || 
               (output_scheme == COMPUTER_READABLE_MATRIX))
      {
        output_matrix_long(diffl, tree1, tree2, trees_in_1, trees_in_2);
      }
      break;
  }
}


void tree_diff(group_type **tree1, group_type **tree2, double *lengths1,
               double* lengths2, long patternsz1, long patternsz2,
               long ntree1, long ntree2, long trees_in_1, long trees_in_2)
{
  long diffl;
  double diffd;

  switch (dtype)
  {
    case SYMMETRIC:
      diffl = symetric_diff (tree1, tree2, ntree1, ntree2, patternsz1, patternsz2);
      diffl += symetric_diff (tree2, tree1, ntree1, ntree2, patternsz2, patternsz1);
      output_long_distance(diffl, ntree1, ntree2, trees_in_1, trees_in_2);
      break;
    case BSD:
    case RF:
      diffd = bsd_tree_diff(tree1, tree2, ntree1, ntree2, lengths1, lengths2, patternsz1, patternsz2);
      output_double_distance(diffd, ntree1, ntree2, trees_in_1, trees_in_2);
      break;
  }
} /* tree_diff */


int get_num_columns(void)
{
  if (output_scheme == COMPUTER_READABLE_MATRIX)
    return (10*trees2);
  else if ( dtype == SYMMETRIC )
        return 10;
    else return 7;
} /* get_num_columns */


void compute_distances(pattern_elm ***pattern_array, long trees_in_1, long trees_in_2)
{
  /* Compute symmetric distances between arrays of trees */
  // RSGnote: Variable "diff_index" formerly initialized and never used; removed.
  long  tree_index, end_tree, index1, index2, index3;
  group_type **treeA, **treeB;
  long patternsz1, patternsz2;
  double *length1 = NULL, *length2 = NULL;
  int num_columns = get_num_columns();

  index1 = 0;
  /* Put together space for treeA and treeB */
  treeA = (group_type **) Malloc (setsz * sizeof (group_type *));
  treeB = (group_type **) Malloc (setsz * sizeof (group_type *));
  print_header(trees_in_1, trees_in_2);
  switch (tree_pairing)
  {
    case ADJACENT_PAIRS:
      end_tree = trees_in_1 - 1;
      for (tree_index = 0 ; tree_index < end_tree ; tree_index += 2) {
        /* For every tree, compute the distance between it and the tree
           at the next location; do this in both directions */
        assign_tree (treeA, pattern_array, tree_index, &patternsz1);
        assign_tree (treeB, pattern_array, tree_index + 1, &patternsz2);
        assign_lengths(&length1, pattern_array, tree_index);
        assign_lengths(&length2, pattern_array, tree_index + 1);
        tree_diff (treeA, treeB, length1, length2, patternsz1, patternsz2, tree_index+1, tree_index+2, trees_in_1, trees_in_2);
        if (tree_index + 2 == end_tree)
          printf("\nWARNING: extra tree at the end of input tree file.\n");
      }
      break;

    case ALL_IN_FIRST:
      end_tree   = trees_in_1;

      if ((output_scheme == FULL_MATRIX) || 
               (output_scheme == COMPUTER_READABLE_MATRIX))
      {
        for (index1 = 0 ; index1 < end_tree ; index1++)
        {
          /* For every tree, compute the distance between it and every
             other tree in that file. */
          assign_tree (treeA, pattern_array, index1, &patternsz1);
          assign_lengths(&length1, pattern_array, index1);

          for (index2 = 0 ; index2 < end_tree ; index2++)
          {
            assign_tree (treeB, pattern_array, index2, &patternsz2);
            assign_lengths(&length2, pattern_array, index2);
            tree_diff (treeA, treeB, length1, length2, patternsz1, patternsz2, index1 + 1, index2 + 1, trees_in_1, trees_in_2);
          }
        }
      }
      else
      {
        for ( index3 = 0 ; index3 < trees_in_1 ; index3 += num_columns)
        {
          for ( index1 = 0 ; index1 < trees_in_1 ; index1++)
          {
            assign_tree (treeA, pattern_array, index1, &patternsz1);
            assign_lengths(&length1, pattern_array, index1);
            for ( index2 = index3 ;
                  index2 < index3 + num_columns && index2 < trees_in_1 ;
                  index2++)
            {
              assign_tree (treeB, pattern_array, index2, &patternsz2);
              assign_lengths(&length2, pattern_array, index2);
              tree_diff (treeA, treeB, length1, length2, patternsz1, patternsz2, index1 + 1, index2 + 1, trees_in_1, trees_in_2);
            }
          }
        }
      }
      break;

    case CORR_IN_1_AND_2:
      if (trees_in_1 != trees_in_2)
      {
        /* Set end tree to the smaller of the two totals. */
        end_tree = trees_in_1 > trees_in_2 ? trees_in_2 : trees_in_1;

        /* Print something out to the outfile and to the terminal. */
        fprintf (outfile, "\n\n");
        fprintf (outfile, "*** Warning: differing number of trees in first and second\n");
        fprintf (outfile, "*** tree files.  Only computing %ld pairs.\n\n", end_tree);
        printf ("\n *** Warning: differing number of trees in first and second\n");
        printf (" *** tree files.  Only computing %ld pairs.\n\n", end_tree);
      }
      else
        end_tree = trees_in_1;

      for (tree_index = 0 ; tree_index < end_tree ; tree_index++)
      {
        /* For every tree, compute the distance between it and the
           tree at the parallel location in the other file; do this in
           both directions */

        assign_tree (treeA, pattern_array, tree_index, &patternsz1);
        assign_lengths(&length1, pattern_array, tree_index);
        /* (tree_index + trees_in_1) will be the corresponding tree in
           the second file. */
        assign_tree (treeB, pattern_array, tree_index + trees_in_1, &patternsz2);
        assign_lengths(&length2, pattern_array, tree_index + trees_in_1);
        tree_diff (treeA, treeB, length1, length2, patternsz1, patternsz2, tree_index + 1, 0, trees_in_1, trees_in_2);
      }
      break;

    case ALL_IN_1_AND_2:
      end_tree = trees_in_1 + trees_in_2;

      if ((output_scheme == FULL_MATRIX) || 
               (output_scheme == COMPUTER_READABLE_MATRIX))
      {
        for (tree_index = 0 ; tree_index < trees_in_1 ; tree_index++)
        {
          /* For every tree in the first file, compute the distance
             between it and every tree in the second file. */
          assign_tree (treeA, pattern_array, tree_index, &patternsz1);
          assign_lengths(&length1, pattern_array, tree_index);
          for (index2 = trees_in_1 ; index2 < end_tree ; index2++)
          {
            assign_tree (treeB, pattern_array, index2, &patternsz2);
            assign_lengths(&length2, pattern_array, index2);
            tree_diff(treeA, treeB, length1, length2, patternsz1, patternsz2, tree_index + 1 , index2 - trees_in_1 + 1, trees_in_1, trees_in_2);
          }
        }
      }
      else
      {
        for ( index3 = trees_in_1 ; index3 < end_tree ; index3 += num_columns)
        {
          for ( index1 = 0 ; index1 < trees_in_1 ; index1++) {
            assign_tree (treeA, pattern_array, index1, &patternsz1);
            assign_lengths(&length1, pattern_array, index1);
            for ( index2 = index3 ;
                  index2 < index3 + num_columns && index2 < end_tree ;
                  index2++)
            {
              assign_tree (treeB, pattern_array, index2, &patternsz2);
              assign_lengths(&length2, pattern_array, index2);
              tree_diff (treeA, treeB, length1, length2, patternsz1, patternsz2, index1 + 1, index2 - trees_in_1 + 1, trees_in_1, trees_in_2);
            }
          }
        }
      }
      break;
  }
  /* Free up treeA and treeB */
  free (treeA);
  free (treeB);
}  /* compute_distances */


void free_patterns(pattern_elm ***pattern_array, long total_trees)
{
  long i, j;
  long end_pattern = total_trees - 1;

  /* Free each pattern array, */
  for (i=0 ; i < setsz ; i++)
  {
    for (j = 0 ; j < end_pattern ; j++)   /* debug -- one off? */
    {
      free (pattern_array[i][j]->apattern);
      free (pattern_array[i][j]->patternsize);
    }
    free (pattern_array[i]);
  }
  free (pattern_array);
}  /* free_patterns */


void print_header(long trees_in_1, long trees_in_2)
{
  long end_tree;

  (void)trees_in_2;                     // RSGnote: Parameter never used.

  switch (tree_pairing)
  {
    case ADJACENT_PAIRS:
      end_tree = trees_in_1 - 1;

      if (output_scheme == VERBOSE)
      {
        fprintf(outfile, "\nTree distance program, version %s\n\n", VERSION);
        if (dtype == BSD)
          fprintf (outfile, "Branch score distances between adjacent pairs of trees:\n\n");
        else
        {
          if (dtype == SYMMETRIC)
            fprintf (outfile, "Symmetric differences between adjacent pairs of trees:\n\n");
          else
            fprintf (outfile, "Robinson-Foulds distances between adjacent pairs of trees:\n\n");
        }
      }
      else if ( output_scheme != SPARSE)
        printf ("Error -- cannot output adjacent pairs into a full matrix.\n");
      break;

    case ALL_IN_FIRST:
      end_tree   = trees_in_1;

      if (output_scheme == VERBOSE)
      {
        fprintf(outfile, "\nTree distance program, version %s\n\n", VERSION);
        if (dtype == BSD)
          fprintf (outfile, "Branch score distances between all pairs of trees in tree file\n\n");
        else {
          if (dtype == SYMMETRIC)
            fprintf (outfile, "Symmetric differences between all pairs of trees in tree file:\n\n");
          else
            fprintf (outfile, "Robinson-Foulds distances between all pairs of trees in tree file:\n\n");
        }
      }
      else if (output_scheme == FULL_MATRIX)
      {
        fprintf(outfile, "\nTree distance program, version %s\n\n", VERSION);
        if (dtype == BSD)
          fprintf (outfile, "Branch score distances between all pairs of trees in tree file:\n\n");
        else {
          if (dtype == SYMMETRIC)
            fprintf (outfile, "Symmetric differences between all pairs of trees in tree file:\n\n");
          else
            fprintf (outfile, "Robinson-Foulds distances between all pairs of trees in tree file:\n\n");
        }
      }
      break;

    case CORR_IN_1_AND_2:

      if (output_scheme == VERBOSE)
      {
        fprintf(outfile, "\nTree distance program, version %s\n\n", VERSION);
        if (dtype == BSD) {
          fprintf (outfile, "Branch score distances between corresponding pairs of trees\n");
          fprintf (outfile, "   from first and second tree files:\n\n");
        }
        else
        {
          if (dtype == SYMMETRIC)
            fprintf (outfile, "Symmetric differences between corresponding pairs of trees\n");
          else
            fprintf (outfile, "Robinson-Foulds distances between corresponding pairs of trees\n");
          fprintf (outfile, "   from first and second tree files:\n\n");
        }
      }
      else if (output_scheme != SPARSE)
        printf ("Error -- cannot output corresponding pairs into a full matrix.\n");
      break;

    case (ALL_IN_1_AND_2) :
      if ( output_scheme == VERBOSE)
      {
        fprintf(outfile, "\nTree distance program, version %s\n\n", VERSION);
        if (dtype == BSD) {
          fprintf (outfile, "Branch score distances between all pairs of trees\n");
          fprintf (outfile, "   from first and second tree files:\n\n");
        }
        else
        {
          if (dtype == SYMMETRIC)
            fprintf(outfile, "Symmetric differences between all pairs of trees\n");
          else
            fprintf(outfile, "Robinson-Foulds distances between all pairs of trees\n");
          fprintf(outfile, "   from first and second tree files:\n\n");
        }
      }
      else if ( output_scheme == FULL_MATRIX)
      {
        fprintf(outfile, "\nTree distance program, version %s\n\n", VERSION);
      }
      break;
  }
  (void)end_tree;                       // RSGnote: Variable set but never read.
} /* print_header */


void output_submenu(void)
{
  /* this allows the user to select a different output of distances scheme. */
  long loopcount;
  boolean done = false;
  Char    ch;

  if (tree_pairing == NO_PAIRING)
    return;

  loopcount = 0;
  while (!done) {
    printf ("\nDistances output options:\n");

    if ((tree_pairing == ALL_IN_1_AND_2) || (tree_pairing == ALL_IN_FIRST)) {
      printf (" F     Full matrix.\n");
      printf (" C     Full matrix with no headings, computer-readable\n");
    }
    printf (" V     One pair per line, verbose.\n");
    printf (" S     One pair per line, sparse.\n");

    if ((tree_pairing == ALL_IN_1_AND_2) || (tree_pairing == ALL_IN_FIRST))
      printf ("\n Choose one: (F,V,S)\n");
    else
      printf ("\n Choose one: (V,S)\n");

    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    uppercase(&ch);

    if (strchr("CFVS", ch) != NULL)
    {
      switch (ch)
      {
        case 'F':
          if ((tree_pairing == ALL_IN_1_AND_2) || (tree_pairing == ALL_IN_FIRST))
            output_scheme = FULL_MATRIX;
          else
            /* If this can't be a full matrix... */
            continue;
          break;

        case 'C':
          if ((tree_pairing == ALL_IN_1_AND_2) || (tree_pairing == ALL_IN_FIRST))
            output_scheme = COMPUTER_READABLE_MATRIX;
          else
            /* If this can't be a full matrix... */
            continue;
          break;

        case 'V':
          output_scheme = VERBOSE;
          break;

        case 'S':
          output_scheme = SPARSE;
          break;
      }
      done = true;
    }
    countup(&loopcount, 10);
  }
}  /* output_submenu */


void pairing_submenu(void)
{
  /* this allows the user to select a different tree pairing scheme. */
  long loopcount;
  boolean done = false;
  Char    ch;

  loopcount = 0;
  while (!done)
  {
    cleerhome();
    printf ("Tree Pairing Submenu:\n");
    printf (" A     Distances between adjacent pairs in tree file.\n");
    printf (" P     Distances between all possible pairs in tree file.\n");
    printf (" C     Distances between corresponding pairs in one tree file and another.\n");
    printf (" L     Distances between all pairs in one tree file and another.\n");

    printf ("\n Choose one: (A,P,C,L)\n");

    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    uppercase(&ch);

    if (strchr("APCL", ch) != NULL)
    {
      switch (ch)
      {
        case 'A':
          tree_pairing = ADJACENT_PAIRS;
          break;

        case 'P':
          tree_pairing = ALL_IN_FIRST;
          break;

        case 'C':
          tree_pairing = CORR_IN_1_AND_2;
          break;

        case 'L':
          tree_pairing = ALL_IN_1_AND_2;
          break;
      }
      output_submenu();
      done = true;
    }
    countup(&loopcount, 10);
  }
}  /* pairing_submenu */


void read_second_file(pattern_elm ***pattern_array, long trees_in_1, long trees_in_2)
{
  boolean firsttree2, haslengths;
  long nextnode, trees_read=0;

  (void)trees_in_2;                     // RSGnote: Parameter never used.

  firsttree2 = false;
  while (!eoff(intree2))
  {
    goteof = false;
    nextnode = 0;
    haslengths = false;
    allocate_nodep(&treeForNodes->nodep, intree2, &spp);
    treeread(treeForNodes, intree2, &treeForNodes->root, treeForNodes->nodep, &goteof, &firsttree2, &nextnode, &haslengths, initconsnode, false, -1);
    missingname(treeForNodes->root);
    reordertips(treeForNodes);
    if (goteof)
      continue;
    ntrees += trweight;
    if (noroot)
    {
      reroot(treeForNodes, treeForNodes->nodep[outgrno - 1], &nextnode);
      didreroot = outgropt;
    }
    accumulate(treeForNodes, treeForNodes->root);
    // BUG.966 -- original code emptied treeForNodes->[nodep] here

    store_pattern (pattern_array, trees_in_1 + trees_read);
    trees_read++;
  }
}  /* read_second_file */


void doinit(void)
{
    if (!javarun)
    {
        getoptions();
    }
}


void getoptions(void)
{
  /* interactively set options */
  long loopcount;
  Char ch;
  boolean done;

  /* Initial settings */
  dtype          = BSD;
  tree_pairing   = ADJACENT_PAIRS;
  output_scheme  = VERBOSE;
  ibmpc          = IBMCRT;
  ansi           = ANSICRT;
  didreroot      = false;
  spp            = 0;
  col            = 0;

  putchar('\n');
  noroot = true;
  numopts = 0;
  outgrno = 1;
  outgropt = false;
  progress = true;

  /* The following are not used by treedist, but may be used
     in functions in cons.c, so we set them here. */
  treeprint = false;
  trout = false;
  prntsets = false;

  loopcount = 0;
  do {
    cleerhome();
    printf("\nTree distance program, version %s\n\n", VERSION);
    printf("Settings for this run:\n");

    printf(" D                         Distance Type:  ");
    switch (dtype)
    {
      case SYMMETRIC:
        printf("Symmetric Difference\n");
        break;
      case BSD:
        printf("Branch Score Distance\n");
        break;
      case RF:
        printf("Robinson-Foulds Distance\n");
        break;
    }
    printf(" R         Trees to be treated as Rooted:");
    if (noroot)
      printf("  No\n");
    else
      printf("  Yes\n");
    printf(" T    Terminal type (IBM PC, ANSI, none):");
    if (ibmpc)
      printf("  IBM PC\n");
    if (ansi)
      printf("  ANSI\n");
    if (!(ibmpc || ansi))
      printf("  (none)\n");
    printf(" 1  Print indications of progress of run:  %s\n", (progress ? "Yes" : "No"));

    printf(" 2                 Tree distance submenu:");
    switch (tree_pairing) {
      case NO_PAIRING:
        printf("\nERROR:  Unallowable option!\n\n");
        exxit(-1);
        break;

      case ADJACENT_PAIRS:
        printf("  Distance between adjacent pairs\n");
        break;

      case CORR_IN_1_AND_2:
        printf("  Distances between corresponding \n");
        printf("                                             pairs in first and second tree files\n");
        break;

      case ALL_IN_FIRST:
        printf("  Distances between all possible\n");
        printf("                                             pairs in tree file.\n");
        break;

      case ALL_IN_1_AND_2:
        printf("  Distances between all pairs in\n");
        printf("                                              first and second tree files\n");
        break;
    }

    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if ((noroot && (ch == 'O')) || strchr("RTD12", ch) != NULL)
      {
        switch (ch)
        {
          case 'D':
            if ( dtype == SYMMETRIC )
              dtype = RF;
            else {if ( dtype == BSD )
                dtype = SYMMETRIC;
              else if ( dtype == RF)
                dtype = BSD;
            }
            break;

          case 'R':
            noroot = !noroot;
            break;

          case 'T':
            initterminal(&ibmpc, &ansi);
            break;

          case '1':
            progress = !progress;
            break;

          case '2':
            pairing_submenu();
            break;
        }
      } else
        printf("Not a possible option!\n");
    }
    countup(&loopcount, 100);
  } while (!done);
}  /* getoptions */


void treedistrun(void)
{
  pattern_elm  ***pattern_array;
  long tip_count = 0, trees_in_1 = 0, trees_in_2 = 0;
  long i, j;

  // debug printout // JRMdebug
  /*
    printf("dtype: %i\n", dtype);
    printf("tree_pairing: %li\n", tree_pairing);
    printf("output_scheme: %li\n", output_scheme);
    printf("didreroot: %i\n", didreroot);
    printf("spp: %li\n", spp);
    printf("col: %li\n", col);
    printf("noroot: %i\n", noroot);
    printf("outgropt: %i\n", outgropt);
    printf("outgrno: %li\n", outgrno);
  */
  // do the work
  ntrees = 0.0;
  lasti  = -1;

  // how many tips?
  countcomma(intree, &tip_count);
  tip_count++; // countcomma does a raw comma count, tips is one greater

  treeForNodes = funcs->tree_new(4*tip_count, 2*tip_count);

  // how many trees do we have?
  trees_in_1 = countsemic(intree);
  if ((tree_pairing == ALL_IN_1_AND_2) ||
      (tree_pairing == CORR_IN_1_AND_2)) {
    // If another intree file should exist,
    trees_in_2 = countsemic(intree2);
    trees2 = trees_in_2;
  }

  // Since "maxgrp" is a limit on the number of items we'll need to put
  // in a hash, we double it to make space for quick hashing
  // limit chosen to make hash arithmetic work
  // we only need enough for spp-1 groups for each tree

  maxgrp = 4*(tip_count-1);

  // Read the (first) tree file and put together grouping, order, and timesseen
  read_groups (&pattern_array, trees_in_1, trees_in_1 + trees_in_2, intree);

  if ((tree_pairing == ADJACENT_PAIRS) || (tree_pairing == ALL_IN_FIRST))
  {
    // Here deal with the adjacent or all-in-first pairing difference computation
    compute_distances (pattern_array, trees_in_1, 0);
  }
  else if ((tree_pairing == CORR_IN_1_AND_2) || (tree_pairing == ALL_IN_1_AND_2))
  {
    // Here, open the other tree file, parse it, and then put together the difference array
    read_second_file (pattern_array, trees_in_1, trees_in_2);
    compute_distances (pattern_array, trees_in_1, trees_in_2);
  }
  else if (tree_pairing == NO_PAIRING)
  {
    // Compute the consensus tree.
    putc('\n', outfile);
    consensus(treeForNodes, pattern_array, trees_in_1+trees_in_2);
  }

  // Free all the buffers needed to compute the differences.
  if (progress)
  {
    sprintf(progbuf, "\nOutput written to file \"%s\".\n\n", outfilename);
    print_progress(progbuf);
  }

  for(i=0; i < setsz; i++)
  {
    for(j=0; j < trees_in_1+trees_in_2; j++)
    {
      // debug   -- not clear we need to free this anyway
      // Patterns need to be freed in a more complex fashion.
      // This removed 'cause it was causing problems:
      // free_patterns (pattern_array, trees_in_1 + trees_in_2);
      // the problem is that pattern_array[i][j][k] has different
      // possible ranges for k depending on i and j
      free(pattern_array[i][j]);
    }
    free(pattern_array[i]);
  }
  free(pattern_array);
}


void treedist(
  char *intreename,
  char *intree2name,
  char *OutfileName,
  char *outfileopt,
  char *DistType,
  char *DistKind,
  int Rooted,
  char *OutputKind,
  int PrintInd)
{
  initdata* funcs;
  int argc;
  Char *argv[1];
  //printf("Hello from TreeDist!\n"); // JRMdebug

  argc = 1;
  argv[0] = "Treedist";

  funcs = Malloc(sizeof(initdata));
  funcs->node_new = cons_node_new;
  phylipinit(argc, argv, funcs, true);
  /*
  //dtype          = BSD;
  //tree_pairing   = ADJACENT_PAIRS;
  //output_scheme  = VERBOSE;
  didreroot      = false;
  spp            = 0;
  col            = 0;
  //noroot = true;
  numopts = 0;
  outgrno = 1;
  outgropt = false;
  //progress = true;
  char *intree,
  char *intree2,
  char *outfile,
  char *outfileopt,
  //char *DistType,
  //char *DistKind,
  //int Rooted,
  //char *OutputType,
  //int PrintInd);
  */
  if (!strcmp(DistType, "BSD"))
  {
    dtype = BSD;
  }
  else if (!strcmp(DistType, "SYMMETRIC"))
  {
    dtype = SYMMETRIC;
  }
  else // RF
  {
    dtype = RF;
  }

  if (!strcmp(DistKind, "ADJACENT_PAIRS"))
  {
    tree_pairing = ADJACENT_PAIRS;
  }
  else if (!strcmp(DistKind, "ALL_IN_FIRST"))
  {
    tree_pairing = ALL_IN_FIRST;
  }
  else if (!strcmp(DistKind, "CORR_IN_1_AND_2"))
  {
    tree_pairing = CORR_IN_1_AND_2;
  }
  else // ALL_IN_1_AND_2
  {
    tree_pairing = ALL_IN_1_AND_2;
  }

  if (Rooted != 0)
  {
    noroot = false;
  }
  else
  {
    noroot = true;
  }

  if (!strcmp(OutputKind, "FULL_MATRIX"))
  {
    output_scheme = FULL_MATRIX;
  }
  else if (!strcmp(OutputKind, "VERBOSE"))
  {
    output_scheme = VERBOSE;
  }
  else // SPARSE
  {
    output_scheme = SPARSE;
  }

  if (PrintInd != 0)
  {
    progress = true;
  }
  else
  {
    progress = false;
  }

  // everything translated, start the run

  // initialize a bunch of stuff
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  didreroot = false;
  spp = 0;
  col = 0;
  numopts = 0;
  outgrno = 1;
  outgropt = false;

  intree = fopen(intreename, "rb");
  outfile = fopen(OutfileName, outfileopt);
  strcpy(outfilename, OutfileName);

  if (progress)
  {
    progfile = fopen("progress.txt", "w");
    fclose(progfile); // make sure it is there for the Java code to detect
    progfile = fopen("progress.txt", "w");
  }

  if ((tree_pairing == CORR_IN_1_AND_2) || (tree_pairing == ALL_IN_1_AND_2))
  {
    intree2 = fopen(intree2name, "rb");
  }

  doinit();

  treedistrun();  // do the actual work

  FClose(intree);
  FClose(outfile);
  if ((tree_pairing == CORR_IN_1_AND_2) || (tree_pairing == ALL_IN_1_AND_2))
  {
    FClose(intree2);
  }
  free(funcs);
  //printf("\ndone\n"); // JRMdebug
}


int main(int argc, Char *argv[])
{
  initdata* funcs;

#ifdef MAC
  argc = 1;                /* macsetup("Treedist", "");        */
  argv[0] = "Treedist";
#endif

  funcs = Malloc(sizeof(initdata));
  funcs->node_new = cons_node_new;
  phylipinit(argc, argv, funcs, false);
  /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
  openfile(&intree, INTREE, "input tree file", "rb", argv[0], intreename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  /* Initialize option-based variables, then ask for changes regarding
     their values. */
  //getoptions();
  doinit();
  if ((tree_pairing == ALL_IN_1_AND_2) ||
      (tree_pairing == CORR_IN_1_AND_2)) {
    // If another intree file should exist,
    // Open in binary: ftell() is broken for UNIX line-endings under WIN32
    openfile(&intree2, INTREE2, "input tree file 2", "rb", argv[0], intree2name);
  }

  treedistrun();

  free(funcs);
  FClose(outtree);
  FClose(intree);
  FClose(outfile);

  if ((tree_pairing == ALL_IN_1_AND_2) || (tree_pairing == CORR_IN_1_AND_2))
    FClose(intree2);

  printf("Done.\n\n");

#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif

  return 0;
}  /* main */


// End.
