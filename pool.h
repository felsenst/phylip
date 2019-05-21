/* Version 4.0. (c) Copyright 2012-2013 by the University of Washington.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


/* pool.h
 *
 * Pool allocator
 *
 * Rapid re-use of fixed-sized chunks of memory.
 * 2006 Ian Robertson */


#ifndef POOL_H
#define POOL_H

#include <assert.h>
#include <sys/types.h>
#include "Slist.h"


typedef struct Pool * Pool;

struct Pool {
  size_t size;          /* Size of each chunk */
  long   block_size;    /* Number of chunks per malloc'd block*/
  Slist  free;          /* Stack of free chunks */
  Slist  blocks;        /* List of malloc'd blocks */
};


/* Create a new Pool object to allocate chunks of memory of a given size.
 * malloc block_size chunks at a time */
extern Pool Pool_new(size_t size, long block_size);

/* Delete a pool and free all memory allocated. Unless all chunks have been
 * returned to the pool, an assertion is thrown. The pool is set to NULL */
extern void Pool_delete(Pool * pool_ptr);


/* Add a block to to pool. This is called automatically when needed by
 * Pool_alloc() */
extern void Pool_add_block(Pool pool);


/* Get a block from the pool. */
static inline void * Pool_alloc(Pool pool)
{
  if ( Slist_isempty(pool->free) )
  {
    Pool_add_block(pool);
  }

  return Slist_pop(pool->free);
}


/* Private debug function */
extern int _Pool_chunk_owned(Pool, void *);

static inline void Pool_free(Pool pool, void * chunk)
{
  assert( _Pool_chunk_owned(pool, chunk) );
  Slist_push(pool->free, chunk);
}

#endif /* POOL_H */


// End.
