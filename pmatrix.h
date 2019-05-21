/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifndef PMATRIX_H
#define PMATRIX_H

/*** Output options ***/

/* Number of significant figures to display in numeric output */
#define OUTPUT_PRECISION        6

/* Maximum line length of matrix output - 0 for unlimited */
#define OUTPUT_TEXTWIDTH        78

/** output_matrix() flags **/

/* Block output: Matrices are vertically split into blocks that
 * fit within OUTPUT_TEXTWIDTH columns */
#define MAT_BLOCK       0x1
/* Lower triangle: Values on or above the diagonal are not printed */
#define MAT_LOWER       0x2
/* Print a border between headings and data */
#define MAT_BORDER      0x4
/* Do not print the column header */
#define MAT_NOHEAD      0x8
/* Output the number of columns before the matrix */
#define MAT_PCOLS       0x10
/* Do not enforce maximum line width */
#define MAT_NOBREAK     0x20
/* Pad row header with spaces to 10 char */
#define MAT_PADHEAD     0x40
/* Human-readable format. */
#define MAT_HUMAN       MAT_BLOCK
/* Machine-readable format. */
#define MAT_MACHINE     (MAT_PCOLS | MAT_NOHEAD | MAT_PADHEAD)
/* Lower-triangular format. */
#define MAT_LOWERTRI    (MAT_LOWER | MAT_MACHINE)

extern char **stringnames_new(void);
extern void stringnames_delete(char **names);
extern void output_matrix_d(FILE *fp, double **matrix,
    unsigned long rows, unsigned long cols,
    char **row_head, char **col_head, unsigned int flags);

#endif /* PMATRIX_H */


// End.
