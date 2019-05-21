/* Version 4.0. (c) Copyright 2012-2013 by the University of Washington.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>

#include "phylip.h"


/* When this is true, comments can be nested: [a[b]c] */
#define NEWICK_NESTED_COMMENTS  1
/* Buffer size for fgets */
#define LINE_BUFSIZE            0x400

#define SPACE_STR               " "
#define WHITESPACE_CHARS        SPACE_STR "\t\n\r"
#define SINGLE_QUOTE_CHAR       '\''


/* Object for reading trees from a file */
typedef struct Usertrees *Usertrees;

struct Usertrees {
  FILE *fp;
  long ntrees;
  const char *current;
};

/* Validated, but unparsed tree */
typedef struct Usertree *Usertree;

struct Usertree {
  long nspecies;
  long nnodes;
  const char *tree;
};

typedef struct Filereader *Filereader;
typedef struct Ut_node *Ut_node;
typedef enum Ut_node_type Ut_node_type;

enum Ut_node_type { TIP, FORK };

struct Ut_node {
  Ut_node_type  type;
  char *        name;
  double        length;
  int           length_defined;
  Slist         children;
};


/* Object for reading line-oriented files in general */
typedef struct Filereader *Filereader;

struct Filereader {
  FILE *fp;
  const char *buf;
  long bufsize;
  const char *next;
  long line;
};


Filereader Filereader_open(const char *filename)
{
  FILE *fp;
  Filereader reader;

  fp = fopen(filename, "r");
  if ( fp == NULL ) {
    return NULL;
  }

  reader = Malloc(sizeof(struct Filereader));
  if ( reader == NULL )
    return NULL;

  reader->buf = Malloc(sizeof(char) * LINE_BUFSIZE);
  if ( reader.buf == NULL ) {
    free(reader);
    return NULL;
  }

  reader->fp = fp;
  reader->bufsize = sizeof(char) * LINE_BUFSIZE;
  reader->next = reader->buf;
  reader->line = 1;

  return reader;
}


void Filereader_close(Filereader *fr_ptr)
{
  Filereader reader;

  assert( fr_ptr != NULL );
  reader = *fr_ptr;
  assert( reader != NULL );

  fclose(reader->fp);
  free(reader->buf);
  free(reader);
}


char Filereader_getchar(Filereader fr)
{
  char ch;

  assert( fr != NULL );

  ch = getc(fr->fp);
  if ( ch == EOF )
    return EOF;

  /* Convert mac '\r' and dos '\r\n' to C '\n' if reading in binary mode */
  if ( ch == '\r' ) {
    ch = getc(file);
    if ( ch != '\n' )
      ungetc(ch, file);
    ch = '\n';
  }

  /* Keep track of the line number */
  if ( ch == '\n' )
    fr->line++;

  return ch;
}


typedef struct Newick_filereader *Newick_filereader;

struct Newick_filereader {
  Filereader filereader;
  long paren_level;
  long comment_level;
  boolean is_quoted;
  long treenum;
};


Newick_filereader Newick_filereader_open(const char *filename)
{
  Newick_filereader nfr;
  Filereader fr;

  fr = Filereader_open(filename);

  if ( fr == NULL )
    return NULL;

  nfr = Malloc(sizeof(struct Newick_filereader));
  if ( nfr == NULL )
    return NULL;

  nfr->filereader = fr;
  nfr->paren_level = 0;
  nfr->comment_level = 0;
  nfr->is_quoted = false;
  nfr->treenum = 0;

  /* TODO: Scan the tree file for correct syntax and number of trees */

  return nfr;
}


void Newick_filereader_close(Newick_filereader *nfr_ptr)
{
  Newick_filereader nfr;

  assert( nfr_ptr != NULL );
  nfr = *nfr_ptr;
  assert( nfr != NULL );

  Filereader_close(&nfr->filereader);

  free(nfr);
}


Usertree Newick_filereader_next(Newick_filereader nfr)
{
  char ch;
  String treestr;
  long  nnodes;
  long  nspecies;

  treestr = String_new(0x100);
  nspecies = 1;         /* number of ',' + 1 */
  nnodes = 0;           /* number of '(' + nspecies */

  /* Skip to start of tree: first open paren outside comments
   * assuming pre: nfr->paren_level == 0 */
  while ( (ch = Newick_filereader_getchar(nfr)) != EOF ) {
    if ( ch = '(' && nfr->comment_level == 0 ) {
      String_append_char(treestr, ch);
      nnodes = 1;
      break;
    }
  }

  if ( ch == EOF ) {
    /* TODO Handle error */
  }

  while ( (ch = Filereader_getchar(nfr->filereader)) != EOF ) {
    if ( ch == SINGLE_QUOTE_CHAR ) {
      /* Quotes are special inside the outer parens, but not in comments */
      if ( nfr->paren_level > 0 && nfr->comment_level == 0 ) {
        nfr->is_quoted = !nfr->is_quoted;
      }
      else {
        /* TODO: warn or ignore? (or quote?)
         * What are the rules at paren_level 0 */
      }
    }
    else if ( !nfr->is_quoted ) {
      /* Context: not quoted  */
      /* Comments may start and end anywhere except inside quotes */
      if ( ch == '[' ) {
        nfr->comment_level++;
        if ( !NEWICK_NESTED_COMMENTS && nfr->comment_level > 1 )
          nfr->comment_level = 1;
      }
      else if ( ch == ']' ) {
        nfr->comment_level--;
        if ( nfr->comment_level < 0 ) {
          nfr->comment_level = 0;
          if ( NEWICK_NESTED_COMMENTS ) {
            /* TODO: warn about mismatched [] */
            tree_parse_error("mismatched ']'", nfr->filereader->line, nfr->treenum);
          }
        }
        /* Do not include this character in the string, even though commentlevel
         * may be 0 */
        continue;
      }
      else if ( nfr->comment_level == 0 ) {
        /* context: not quoted, not commented */
        if ( ch == '(' ) {
          nfr->paren_level++;
          nnodes++;
        }
        else if ( ch == ')' ) {
          nfr->paren_level--;
          if ( nfr->paren_level < 0 ) {
            nfr->paren_level = 0;
            /* TODO: warn about mismatched () */
            tree_parse_error("mismatched parentheses", nfr->filereader->line, nfr->treenum);
          }
        }
        else if ( ch == ',' ) {
          nspecies++;
          nnodes++;
        }
        if ( *c == ';' ) {
          if ( paren_level == 0 ) {
            nfr->ntrees++;
            nfr->treenum++;
            break;
          }
          else {
            /* warn semicolon not valid */
            tree_parse_error("unexpected ';' (mismatched parentheses?)",
                             nfr->filereader->line, nfr->treenum);
          }
        }
      } /* comment_level == 0 */
    } /* !is_quoted */

    /* Add non-comment characters to string. */
    if ( nfr->comment_level == 0 ) {
      if ( ch != '\0' )
        String_append_char(treestr, ch);
      /* Finish when ';' is reached */
      if ( ch == ';' )
        break;
    }
  }
  else { /* ch == EOF */
    /* Handle unexpected EOF */
  }

  ut = Malloc(sizeof(struct Usertree));
  if ( ut == NULL )
    /* TODO Handle error */;

  ut->tree = treestr;
  ut->nnodes = nnodes;
  ut->nspecies = nspecies;

  return ut;
}


Usertree Newick_filereader_nexttree(Newick_filereader nfr)
{
  /* Read an entire tree from tree file and return a Usertree object */

  char *buf;
  long size;

  assert( nfr != NULL );

  buf = Malloc(LINE_BUFSIZE);
}


void tree_parse_error (const char *msg, long line, long treenum);


Usertrees Usertrees_load_file(const char *filename, Usertrees ret)
{
  FILE *treefile;

  /* Open the file */
  treefile = fopen(filename, "r");

  if ( treefile == NULL ) {
    /* TODO */
  }
}


Usertree Usertree_next(Usertrees)
{
  /* Read next tree from file in Usertrees */
  /* Report parse errors and return Usertree */

  chat buf[LINE_BUFSIZE];
  char *c;

  long comment_level = 0;
  long paren_level = 0;

  assert( Usertrees != NULL );

}


long trees_in_file(FILE *treefile)
{ /* Determine the number of trees in a tree file */

  char buf[LINE_BUFSIZE];
  char *c;

  long n;
  long ntrees_param = 0;
  long ntrees = 0;
  long treenum;
  long comment_level = 0;
  long paren_level = 0;
  long line = 1;

  boolean is_quoted = 0;

  c = buf;

  /* Attempt to read an int at the start */
  if ( fgets(buf, LINE_BUFSIZE, treefile) ) {
    /* Skip whitespace, newlines */
    c += strspn(buf, WHITESPACE_CHARS);
    if ( isdigit(*c) ) {
      n = sscanf(buf, "%ld", &ntrees_param);
      /* Is number valid? (>0) */
      if ( n == 1 && ntrees_param > 0 ) {
        /* If so, just return that number */
        return ntrees_param;
      }
    }
    else {
      /* Count the semicolons outside parens, (possibly-nested)
       * comments [[]] and single-quoted '' names.
       * Issue a warning if parens, quotes, or brackets do not
       * match */

      /* Tree parser states:
       * TPS_PRE  - Before first '(', only comments are special
       * TPS_SUB  - In tree, expecting label or subtree
       * TPS_PRESUB Expecting ')'
       * TPS_LEN  - Expecting ",:)"
       * TPS_LEN2 - Expecting number
       * TPS_ENDSUB Expecting '),'
       * TPS_END  - Expecting ":;"
       * TPS_RLEN - Expecting number
       * TPS_END2 - Expecting ;
       */

      rewind(treefile);

      comment_level = 0;
      paren_level = 0;
      is_quoted = false;
      /* For all lines */
      while ( fgets(buf, LINE_BUFSIZE, treefile) != NULL ) {
        if ( DEBUG ) {
          fprintf(stderr, buf);
        }
        c = buf;
        /* Parse special chars */
        while( c = strpbrk(c, ";[]'()\n,:") ) {
          /* Context: any */
          if ( *c == '\n' ) {
            line++;
          }
          else if ( *c == SINGLE_QUOTE_CHAR ) {
            /* Quotes are special inside the outer parens, but not in comments */
            if ( paren_level > 0 && comment_level == 0 ) {
              is_quoted = !is_quoted;
            }
            else {
              /* TODO: warn or ignore? (or quote?)
               * What are the rules at paren_level 0 */
            }
          }
          else if ( !is_quoted ) {
            /* Context: not quoted  */
            /* Comments may start and end anywhere except inside quotes */
            if ( *c == '[' ) {
              comment_level++;
              if ( !NEWICK_NESTED_COMMENTS && comment_level > 1 )
                comment_level = 1;
            }
            else if ( *c == ']' ) {
              comment_level--;
              if ( comment_level < 0 ) {
                comment_level = 0;
                if ( NEWICK_NESTED_COMMENTS ) {
                  /* TODO: warn about mismatched [] */
                  tree_parse_error("mismatched ']'", line, treenum);
                }
              }
            }
            else if ( comment_level == 0 ) {
              /* context: not quoted, not commented */
              if ( *c == '(' ) {
                paren_level++;
                nnodes++;
              }
              else if ( *c == ')' ) {
                paren_level--;
                if ( paren_level < 0 ) {
                  paren_level = 0;
                  /* TODO: warn about mismatched () */
                  tree_parse_error("mismatched ')'", line, treenum);
                }
              }
              else if ( *c == ',' ) {
                nnodes++;
              }
              if ( *c == ';' ) {
                if ( paren_level == 0 ) {
                  ntrees++;
                  treenum++;
                }
                else {
                  /* TODO: warn semicolon not valid */
                  tree_parse_error("unexpected ';' (mismatched parentheses?)", line, treenum);
                }
              }
            } /* comment_level == 0 */
          } /* !is_quoted */

            /* Move to next char in buf */
          c++;

        } /* while strpbrk */
        if ( DEBUG ) {
          fprintf(stderr, "ntrees_param = %ld\n", ntrees_param);
          fprintf(stderr, "ntrees       = %ld\n", ntrees);
          fprintf(stderr, "paren_level  = %ld\n", paren_level);
          fprintf(stderr, "comment_level= %ld\n", comment_level);
          fprintf(stderr, "is_quoted    = %s\n", is_quoted ? "yes" : "no");
        }
      } /* while fgets */

      if (paren_level > 0 || comment_level > 0 || is_quoted)
      {
        /* TODO */
      }
    } /* !isdigit() */
  } /* fgets */
  else
  {
    /* fgets == NULL */
  }

  /* Rewind the file */

  /* Return Usertrees object */

  /* DEBUG print the results */
  fprintf(stderr, "ntrees_param = %ld\n", ntrees_param);
  fprintf(stderr, "ntrees       = %ld\n", ntrees);
  fprintf(stderr, "paren_level  = %ld\n", paren_level);
  fprintf(stderr, "comment_level= %ld\n", comment_level);

  return ntrees;
} /* tree_parse() */


void tree_parse_error (const char *msg, long line, long treenum)
{
  assert( msg != NULL );
  fprintf(stderr, "Error parsing tree file: %s at line %d.\n", msg, line, treenum);
  fprintf(stderr, "Exiting.");
  exit(-1);
}


int main(void)
{
  FILE *intree;
  intree = fopen("intree", "r");
  if ( intree == NULL )
  {
    perror("ERROR: ");
    abort();
  }

  trees_in_file(intree);
  return 0;
}


#ifdef NOT_DEFINED

void Usertrees_free(Usertrees uts)
{
  /* Free Usertrees object */
}


const char * Usertrees_gettree(void)
{
}


static int Usertree_parse_err = 0;


boolean Usertree_new(const char *treestr)
{
  const char *c;
  long paren_lvl;

  c = treestr;
  paren_lvl = 0;

  /* Find first open paren */
  while ( *c != '(' ) {
    c++;
    if ( *c == '\0' )
      /* fail */;
  }

  /* Parse subtree */
  c = Usertree_parse_subtree(c);

  if ( Usertree_parse_err != 0 )
    return false;

  /* Skip whitespace */
  while ( *c != ';' ) {
    if ( !isspace(*c) )
      return false;
  }

  return true;
}


/* Parse a subtree */

const char * Usertree_parse_subtree(const char *c)
{
  ch = c;
  /* Skip whitespace */

  /* Expect '(' */
  if ( *ch != '(' )
    /* error */;

  c++;
  /* Parse species name or subtree */
  if ( *c == '(' )
    c = Usertree_parse_subtree(c);
  else
  {
    // c = Usertree_parse_label(c);
  }

  /* If optional ':' */
    /* Parse length */
  /* Skip optional '[comment]' */

    /* If ',' */
    /* repeat */

  /* Expect ')' */
}


const char * Usertree_parse_species_name(const char *c)
{
  /* Skip whitespace */

  /* Read name until one of ':,)\0' */

  /* If more than max name length */
  /* fail */
}

#endif /* NOT_DEFINED */


// End.
