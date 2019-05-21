/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


/* pmatrix.c
 *
 * Functions for formatting and printing a matrix.
 */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "pmatrix.h"


/* These functions are temporarily used for translating the fixed-width
 * space-padded nayme array to an array of null-terminated char *. */
char **stringnames_new(void)
{
  /* Copy nayme array to null terminated strings and return array of char *.
   * Spaces are stripped from end of naym's.
   * Returned array size is spp+1; last element is NULL. */

  char **names;
  char *ch;
  long len, i;

  names = (char **)Malloc((spp+1) * sizeof(char *));

  for ( i = 0; i < spp; i++ ) {
    len = strlen(nayme[i]);
    (void)len;                          // RSGnote: Variable set but never referenced.
    names[i] = (char *)Malloc((MAXNCH+1) * sizeof(char));
    strncpy(names[i], nayme[i], MAXNCH);
    names[i][MAXNCH] = '\0';
    /* Strip trailing spaces */
    for ( ch = names[i] + MAXNCH - 1; *ch == ' ' || *ch == '\0'; ch-- )
      *ch = '\0';
  }
  names[spp] = NULL;

  return names;
}


void stringnames_delete(char **names)
{
  /* Free a string array returned by stringnames_new() */
  long i;

  assert( names != NULL );
  for ( i = 0; i < spp; i++ ) {
    assert( names[i] != NULL );
    free(names[i]);
  }
  free(names);
}


static int fieldwidth_double(double val, unsigned int precision)
{
  /* Printf a double to a temporary buffer with specified precision using %g
   * and return its length. Precision must not be greater than 999,999 */

  char format[10];
  char buf[0x200]; /* TODO: What's the largest possible? */

  if (precision > 999999)
    abort();

  sprintf(format, "%%.%uf", precision); /* %.Nf */
  /* snprintf() would be better, but is it avaliable on all systems? */
  return sprintf(buf, format, val);
}


void output_matrix_d(FILE *fp, double **matrix, unsigned long rows, unsigned long cols,
                     char **row_head, char **col_head, unsigned int flags)
{
  /*
   * Print a matrix of double to file. Headings are given in row_head and
   * col_head, either of which may be NULL to indicate that headings should not
   * be printed. Otherwise, they must be null-terminated arrays of pointers to
   * null-terminalted character arrays.
   *
   * The macro OUTPUT_PRECISION defines the number of significant figures to
   * print, and OUTPUT_TEXTWIDTH defines the maximum length of each line.
   *
   * Optional formatting is specified by flags argument, using macros MAT_*
   * defined in pmatrix.h.
   */

  unsigned     *colwidth;               /* [0..spp-1] min width of each column */
  unsigned      headwidth;              /* min width of row header column */
  unsigned long linelen;                /* length of current printed line */
  unsigned      fw;
  unsigned long row, col;
  unsigned long i;
  unsigned long cstart, cend;
  unsigned long textwidth = OUTPUT_TEXTWIDTH;
  const unsigned gutter = 1;            // RSGnote: Changed from SIGNED to UNSIGNED INT.
  boolean       do_block;
  boolean       lower_triangle;
  boolean       border;
  boolean       output_cols;
  boolean       pad_row_head;

  if ( flags & MAT_NOHEAD )
    col_head = NULL;
  if ( flags & MAT_NOBREAK )
    textwidth = 0;
  do_block = (flags & MAT_BLOCK) && (textwidth > 0);
  lower_triangle = flags & MAT_LOWER;
  border = flags & MAT_BORDER;
  output_cols = flags & MAT_PCOLS;
  pad_row_head = flags & MAT_PADHEAD;

  /* Determine minimal width for row headers, if given */
  headwidth = 0;
  if ( row_head != NULL ) {
    for (row = 0; row < rows; row++) {
      fw = strlen(row_head[row]);
      if ( headwidth < fw )
        headwidth = fw;
    }
  }

  /* Enforce minimum of  nmlngth  ch for machine-readable output */
  if ( (pad_row_head) && (headwidth < nmlngth) )
    headwidth = nmlngth;

  /* Determine minimal width for each matrix col */
  colwidth = (unsigned int *)Malloc(spp * sizeof(int));
  for (col = 0; col < cols; col++) {
    if ( col_head != NULL )
      colwidth[col] = strlen(col_head[col]);
    else
      colwidth[col] = 0;
    for (row = 0; row < rows; row++) {
      fw = fieldwidth_double(matrix[row][col], OUTPUT_PRECISION);
      if ( colwidth[col] < fw )
        colwidth[col] = fw;
    }
  }

  /*** Print the matrix ***/
  /* Number of columns if requested */
  if ( output_cols ) {
    fprintf(fp, "%5lu\n", cols);
  }

  /* Omit last column for lower triangle */
  if ( lower_triangle )
    cols--;

  /* Blocks */
  cstart = cend = 0;
  while ( cend != cols ) {
    if ( do_block ) {
      linelen = headwidth;
      for ( col = cstart; col < cols; col++ ) {
        if ( linelen + colwidth[col] + gutter > textwidth ) {
          break;
        }
        linelen += colwidth[col] + gutter;
      }
      cend = col;
      /* Always print at least one, regardless of line len */
      if ( cend == cstart )
        cend++;
    } else {
      cend = cols;
    }

    /* Column headers */
    if ( col_head != NULL ) {
      /* corner space */
      for ( i = 0; i < headwidth; i++ )
        putc(' ', fp);
      if ( border ) {
        for ( i = 0; i < gutter+1; i++ )
          putc(' ', fp);
      }
      /* Names */
      for ( col = cstart; col < cend; col++ ) {
        for ( i = 0; i < gutter; i++ )
          putc(' ', fp);
        /* right justify */
        fw = strlen(col_head[col]);
        for ( i = 0; i < colwidth[col] - fw; i++ )
          putc(' ', fp);
        fputs(col_head[col], fp);
      }
      putc('\n', fp);
    }

    /* Top border */
    if ( border ) {
      for ( i = 0; i < headwidth + gutter; i++ )
        putc(' ', fp);
      putc('\\', fp);
      for ( col = cstart; col < cend; col++ ) {
        for ( i = 0; i < colwidth[col] + gutter; i++ )
          putc('-', fp);
      }
      putc('\n', fp);
    }

    /* Rows */
    for (row = 0; row < rows; row++) {
      /* Row header, if given */
      if ( row_head != NULL ) {
        /* right-justify for non-machine-readable */
        if ( !pad_row_head ) {
          for ( i = strlen(row_head[row]); i < headwidth; i++ )
            putc(' ', fp);
        }
        fputs(row_head[row], fp);
        /* left-justify for machine-readable */
        if ( pad_row_head ) {
          for ( i = strlen(row_head[row]); i < headwidth; i++ )
            putc(' ', fp);
        }
      }
      linelen = headwidth;

      /* Left border */
      if ( border ) {
        for ( i = 0; i < gutter; i++ )
          putc(' ', fp);
        putc('|', fp);
        linelen += 2;
      }

      /* Row data */
      for (col = cstart; col < cend; col++) { /* cols */
        /* Stop after col == row for lower triangle */
        if ( lower_triangle && col >= row )
          break;
        /* Break line if going over max text width */
        if ( !do_block && textwidth > 0 ) {
          if ( linelen + colwidth[col] > textwidth )
          {
            // RSGnote: If printing continuation on next line, print Newline
            // and indent 2 extra spaces before next number starts.
            fprintf(fp, "\n  ");
            linelen = 2;
          }
          linelen += colwidth[col] + gutter;
        }

        for ( i = 0; i < gutter; i++ )
          putc(' ', fp);

        /* Print the datum */
        fprintf(fp, "%*.6f", colwidth[col], matrix[row][col]);
      }
      putc('\n', fp);
    } /* End of row */
//    putc('\n', fp); /* blank line */  debug  (need to put in calling program)
    cstart = cend;
  } /* End of block */
  free(colwidth);
} /* output_matrix_d */


// End.
