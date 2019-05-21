/* Version 4.0. (c) Copyright 2012-2013 by the University of Washington.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


/* restore.c
 *
 * Prototype module for providing rewindable memory operations. Regions of
 * memory can be allocated, changed, and freed, and these operations can then
 * be either rewound or forgotten. Malloc operations result in "newborns,"
 * which are freed if the log is rewound. Free operations result in "zombies,"
 * which are freed if the log is forgotten. Regions of memory can be saved,
 * changed, and subsequently restored to an arbitrary point. Because
 * operations are restored in reverse order and to their original memory
 * locations, this method can safely be used with complex data structures
 * involving pointers.
 *
 * The idea is to provide a versatile tool for local tree rearrangements,
 * where it is desirable to modify parts of a tree temporarily, evaluate, and
 * then possibly restore the tree to a previous state.
 *
 * Performance optimizations may be possible for fixed sized chunks (like
 * nodes) by incorporating a pool allocator to reuse malloc'd memory. This
 * has been included for the saveblocks, and could also be added to the Slist
 * data structure for internal nodes.
 */


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Slist.c"

#ifndef NO_POOL
#  include "pool.c"
#endif


typedef void * (*malloc_t)(size_t);
typedef void (*free_t)(void *);

typedef enum savetype {
  SL_DATA,      /* Stored data */
  SL_NEWBORN,   /* Reversibly allocated space */
  SL_ZOMBIE,    /* Reversibly freed space */
  SL_POINT      /* No-op checkpoint */
} savetype;

typedef struct saveblock {
  savetype type;        /* one of savetype */
  long   id;            /* all: checkpoint ID */
  size_t size;          /* SL_DATA: size of data */
  void *data;           /* SL_DATA: ptr to saved data */
  void *loc;            /* SL_DATA: ptr to restore location
                           SL_NEWBORN, SL_ZOMBIE: ptr to block */
} saveblock;

typedef struct savelog {
  Slist         log;    /* The log itself, a stack of saveblocks */
  long          next_id;/* checkpoint ID of next op */
} savelog;

static Pool sb_pool = NULL;
static long pool_usecount = 0;
static malloc_t sl_malloc = NULL;
static free_t   sl_free   = NULL;

/* Create and return a new savelog */
savelog *savelog_new(void);

/* Forget a savelog and delete */
void    savelog_delete(savelog *log);

/* Make a checkpoint in the log and return its id */
long    remember_checkpoint(savelog *log);

/* Save a region of memory about to be changed. Return a checkpoint. */
long    remember_data(savelog *log, void *loc, size_t size);

/* Remember a recently-allocated block. Return a checkpoint.  */
long    remember_alloc(savelog *log, void *ptr);

/* Remember that a block is no longer used. Return a checkpoint. */
long    remember_free(savelog *log, void *ptr);

/* Reverse the operations until checkpoint id is reached */
void    rewind_to(savelog *log, long id);

/* Forget the logged operations and keep the current state */
void    forget_all(savelog *log);

/* Test the module */
int     main(void);


/* Private block allocator */
static inline saveblock *_new_sb(void)
{
#ifdef NO_POOL
  saveblock *sb;
  sb = (saveblock *)(*sl_malloc)(sizeof(struct saveblock));
  assert(sb != NULL);
  return sb;
#else
  return Pool_alloc(sb_pool);
#endif
}


static inline void _free_sb(saveblock *sb)
{
#ifdef NO_POOL
  free(sb);
#else
  Pool_free(sb_pool, sb);
#endif
}


void savelog_use_malloc(malloc_t malloc_func)
{
  sl_malloc = malloc_func;
}


void savelog_use_free(free_t free_func)
{
  sl_free = free_func;
}


savelog *savelog_new(void)
{
#ifndef NO_POOL
  /* Initialize the saveblock pool if needed */
  if ( sb_pool == NULL ) {
    sb_pool = Pool_new(sizeof(struct saveblock), 10);
    pool_usecount = 1;
  }
  else {
    pool_usecount++;
  }
#endif

  /* Initialize allocator to default */
  if ( sl_malloc == NULL )
    sl_malloc = malloc;
  if ( sl_free == NULL )
    sl_free = free;

  /* Create a new log */
  savelog *sl;
  sl = (savelog *)(*sl_malloc)(sizeof(struct savelog));
  assert( sl != NULL );

  sl->log = Slist_new();
  sl->next_id = 0;

  return sl;
}


void savelog_delete(savelog *log)
{
  /* Forget a log's contents and free its memory */
  forget_all(log);
  Slist_delete(&(log->log));
  (*sl_free)(log);

#ifndef NO_POOL
  /* Free the saveblock pool if unused */
  pool_usecount--;
  if ( pool_usecount == 0 )
    Pool_delete(&sb_pool);
#endif
}


long remember_checkpoint(savelog *log)
{
  /* Do nothing except return a checkpoint to rewind to */
  saveblock *sb;
  sb = _new_sb();

  sb->type = SL_POINT;
  sb->id   = log->next_id;
  log->next_id++;

  /* not strictly necessary */
  sb->size = 0;
  sb->loc  = NULL;
  sb->data = NULL;

  Slist_push(log->log, sb);

  return sb->id;
}


long remember_data(savelog *log, void *loc, size_t size)
{
  /* Store a copy of a region of memory in the log */
  saveblock *sb;

  assert(size > 0);
  assert(loc != NULL);

  sb = _new_sb();

  sb->type = SL_DATA;
  sb->id   = log->next_id;
  log->next_id++;
  sb->size = size;
  sb->loc  = loc;
  sb->data = (*sl_malloc)(size);

  assert( sb->data != NULL );

  memcpy(sb->data, loc, size);

  Slist_push(log->log, sb);

  return sb->id;
}


long remember_alloc(savelog *log, void *ptr)
{
  /* Remember to free a block recently allocated */
  saveblock *sb;

  sb = _new_sb();

  sb->type = SL_NEWBORN;
  sb->size = 0;
  sb->id   = log->next_id;
  log->next_id++;
  sb->data = NULL;
  sb->loc  = ptr;
  assert( sb->loc != NULL );

  Slist_push(log->log, (void *)sb);

  return sb->id;
}


long remember_free(savelog *log, void *ptr)
{
  /* Mark memory to be freed if the log is forgotten */
  saveblock *sb;

  assert(ptr != NULL);

  sb = _new_sb();

  sb->type = SL_ZOMBIE;
  sb->size = 0;
  sb->id   = log->next_id;
  log->next_id++;
  sb->data = NULL;
  sb->loc  = ptr;

  Slist_push(log->log, (void *)sb);

  return sb->id;
}


void rewind_to(savelog *log, long id)
{
  /* Rewind the log before last action returning id. Restore saved data
   * and free newborns. */
  saveblock *sb;
  for (;;) {
    assert( !Slist_isempty(log->log) );
    sb = (saveblock *)Slist_pop(log->log);
    switch (sb->type)
    {
      case SL_DATA:
        memcpy(sb->loc, sb->data, sb->size);
        (*sl_free)(sb->data);
        break;
      case SL_NEWBORN:
        (*sl_free)(sb->loc);
        break;
      case SL_ZOMBIE:
        /* Do nothing */
        break;
      case SL_POINT:
        /* Do nothing */
        break;
      default:
        assert( 0 ); /* Shouldn't happen */
    }
    if ( sb->id == id ) {
      _free_sb(sb);
      break;
    }
    _free_sb(sb);
  }
}


void forget_all(savelog *log)
{
  /* Forget a log. Free log memory and zombies */
  saveblock *sb;
  while ( !Slist_isempty(log->log) ) {
    sb = (saveblock *)Slist_pop(log->log);
    switch (sb->type)
    {
      case SL_DATA:
        (*sl_free)(sb->data);
        break;
      case SL_NEWBORN:
        /* Do nothing */
        break;
      case SL_ZOMBIE:
        (*sl_free)(sb->loc);
        break;
      case SL_POINT:
        /* Do nothing */
        break;
      default:
        assert( 0 ); /* shouldn't happen */
    }
    _free_sb(sb);
  }
}


/* Example alternate malloc/free */
void * the_other_malloc(size_t size)
{
  void *ptr = malloc(size);
  printf("%p: %d bytes.\n", ptr, size);
  return ptr;
}


void the_other_free(void *ptr)
{
  printf("Free %p!\n", ptr);
  return free(ptr);
}


int main(void)
{
  /* Test the module. Valgrind this to check for leaks! */

  double *a, **ap;
  savelog *mylog;
  long  p1, p2, p3;

  savelog_use_malloc(the_other_malloc);
  savelog_use_free(the_other_free);
  mylog = savelog_new();

  /* Reversibly create an array */
  a = malloc(5 * sizeof(double));
  p1 = remember_alloc(mylog, a);
  a[0] = 1.23;
  a[1] = 4.56;
  a[2] = 7.89;
  a[3] = 0.12;
  a[4] = 3.14;
  printf("{%f, %f, %f, %f, %f}\n", a[0], a[1], a[2], a[3], a[4]);
  /* Now an array of pointers into the array */
  ap = the_other_malloc(3 * sizeof(double *));
  remember_alloc(mylog, ap);
  ap[0] = &a[0];
  ap[1] = &a[1];
  ap[2] = &a[2];
  /* Reversibly change the contents */
  p2 = remember_checkpoint(mylog);
  remember_data(mylog, a, 5 * sizeof(double));
  a[0] = a[4];
  a[4] = 6.00;
  a[3] = 2*a[3];
  *(ap[2]) = 495.82;
  *(ap[1]) = a[0] * 4;
  printf("{%f, %f, %f, %f, %f}\n", *ap[0], *ap[1], *ap[2], a[3], a[4]);
  /* Zero the pointers */
  p3 = remember_data(mylog, ap, 3 * sizeof(double *));
  memset(ap, 0, 3 * sizeof(double *));

  /* Reversibly free the array and pointers */
  remember_free(mylog, a);
  remember_free(mylog, ap);

  printf("Keep changes? (y/n) ");
  if ( getchar() == 'y' )
  {
    /* Keep our changes and free log memory */
    forget_all(mylog);
  }
  else
  {
    /* Rewind to before we freed array and zero'd pointers */
    rewind_to(mylog, p3);
    printf("{%f, %f, %f, %f, %f}\n", *ap[0], *ap[1], *ap[2], a[3], a[4]);
    /* Rewind to before we changed data */
    rewind_to(mylog, p2);
    printf("{%f, %f, %f, %f, %f}\n", a[0], a[1], a[2], a[3], a[4]);
    /* Rewind to before we allocated the array */
    rewind_to(mylog, p1);
  }

  savelog_delete(mylog);

  return 0;
}


// End.
