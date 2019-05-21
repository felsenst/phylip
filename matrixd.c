/* Version 4.0. (c) Copyright 2012-2013 by the University of Washington.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdlib.h>

#include "phylip.h"
#include "Slist.h"
#include "matrixd.h"


typedef struct Matrix_double_allocator *Matrix_double_allocator;

struct Md_allocator {
  long          dimx, dimy;     /* major, minor dimensions */
  long          nalloc;         /* Number of total matrices created */
  Slist_ptr     free_list;      /* Stack of free matrices */
  Slist_ptr     alloc_list;     /* List of all matrices created */
};


Matrix_double Matrix_double_new(long dimx, long dimy)
{
  Matrix_double md;
  double **xp, **xp_end;
  double *yp;

  md = Malloc(dimx * (sizeof(double *) + dimy * sizeof(double)));
  xp = md;
  xp_end = xp + dimx;
  yp = (double *)xp_end;
  for ( xp = md; xp < xp_end; xp++ ) {
    *xp = yp;
    yp += dimy;
  }
  return md;
}


void Matrix_double_delete(Matrix_double *md_ptr, long dimx)
{
  Matrix_double md;

  (void)dimx;                           // RSGnote: Parameter never used.

  assert( md_ptr != NULL );
  md = *md_ptr;
  assert( md != NULL );

  free(md);
  *md_ptr = NULL;
}


/* IDEA: possible generic pool allocator. (Possibly too generic):
Allocator Allocator_new( (void *)(obj_new*)(void *), const void *new_arg,
    (void)(obj_del*)(void **) );
*/


Md_allocator Md_allocator_new(long dimx, long dimy, long n)
{
  /* Create an allocator for new transmatrices */
  Md_allocator  mda;
  Matrix_double md;
  long          i;

  assert( dimx > 0 && dimy > 0 );
  assert( n >= 0 );

  mda = Malloc(sizeof(struct Md_allocator));
  assert( mda != NULL );

  mda->dimx       = dimx;
  mda->dimy       = dimy;
  mda->nalloc     = 0;
  mda->free_list  = Slist_new();
  mda->alloc_list = Slist_new();

  /* Create n matrices immediately */
  for ( i = 0; i < n; i++ ) {
    md = Matrix_double_new(dimx, dimy);
    Slist_push(mda->free_list, md);
    Slist_push(mda->alloc_list, md);
    mda->nalloc++;
  }

  return mda;
}


void Md_allocator_delete(Md_allocator *mda_ptr)
{
  /* Free all unused matrices and delete allocator */
  Md_allocator  mda;
  Matrix_double md;

  assert( mda_ptr != NULL );
  mda = *mda_ptr;
  assert( mda != NULL );

  if ( mda->nalloc - Slist_get_length(mda->free_list) != 0 ) {
    /* TODO: Issue a warning or assert when not all matrices
     * have been released? */
  }

  while ( !Slist_isempty(mda->free_list) ) {
    md = (Matrix_double)Slist_pop(mda->free_list);
    free(md);
  }

  Slist_delete(mda->free_list);
  free(mda);
  *mda_ptr = NULL;
}


Matrix_double Md_allocator_get(Md_allocator mda)
{
  Matrix_double md;

  assert( mda != NULL );
  /* Get a new transmatrix from the allocator */
  if ( !Slist_isempty(mda->free_list) ) {
    /* Return a free matrix */
    return (Matrix_double)Slist_pop(mda->free_list);
  } else {
    /* Create a new matrix */
    md = Matrix_double_new(mda->dimx, mda->dimy);
    Slist_push(mda->alloc_list, md);
    return md;
  }
  return NULL; /* shouldn't happen */
}


void Md_allocator_release(Md_allocator mda, Matrix_double md)
{
  assert( md != NULL );
  Slist_push(mda->free_list, md);
}


// End.
