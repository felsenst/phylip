/* Version 4.0. (c) Copyright 2012-2013 by the University of Washington.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifndef MATRIXD_H
#define MATRIXD_H

/* A data type for a two-dimensional array of double */

typedef double** Matrix_double;

Matrix_double Matrix_double_new(long dimx, long dimy);
void Matrix_double_delete(Matrix_double *md, long dimx);


/* An allocator data type for reusable identically-sized matrices */

typedef struct Md_allocator *Md_allocator;

/*
 * Create a new allocator for matrices. All matrices will be sized
 * [dimx][dimy]. The integer n specifies how many to allocate initially. May be
 * zero. Additional matrices will be created as needed.
 */
Md_allocator Md_allocator_new(long dimx, long dimy, long n);

/*
 * Get an unused matrix or create a new one if none is available.
 */
Matrix_double Md_allocator_get(Md_allocator mda);

/*
 * Release a matrix when it is no longer needed. The matrix must be one
 * previously returned by a call to Md_allocator_get(mda).
 */
void Md_allocator_release(Md_allocator mda, Matrix_double md);

/*
 * Free all matrices allocated by mda (whether released or not) and
 * delete the allocator. The pointer to the allocator is set to NULL.
 */
void Md_allocator_delete(Md_allocator *mda_ptr);

#endif /* MATRIXD_H */


// End.
