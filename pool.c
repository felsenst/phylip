/* PHYLIP Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   Dan Fineman, Patrick Colacurcio, and Mike Palczewski.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


/* pool.c
 *
 * Pool allocator. Rapid management of fixed-sized chunks of memory.
 *
 * Original author: Ian Robertson <iroberts@u.washington.edu> */


#include <stdlib.h>
#include <sys/types.h>
#include <assert.h>

#include "pool.h"


Pool Pool_new(size_t size, long block_size)
{
  /* Create a new Pool object to allocate chunks of memory of a given size.
   * Allocate at least nreserve blocks */

  Pool pool;

  assert( size > 0 );

  pool = malloc(sizeof(struct Pool));
  assert( pool );

  pool->size       = size;
  pool->block_size = block_size;
  pool->blocks     = Slist_new();
  pool->free       = Slist_new();

  return pool;
}


/* Delete a pool and free all memory allocated */

void Pool_delete(Pool * pool_ptr)
{
  Pool   pool;

  assert( pool_ptr != NULL );
  pool = *pool_ptr;
  assert( pool != NULL );

  /* check that all chunks have been returned */
  assert( Slist_get_length(pool->blocks) * pool->block_size == Slist_get_length(pool->free) );

  while( !Slist_isempty(pool->free) )
    Slist_pop(pool->free);

  Slist_delete(&(pool->free));

  while( !Slist_isempty(pool->blocks) )
    free( Slist_pop(pool->blocks) );

  Slist_delete(&(pool->blocks));

  free(pool);
  *pool_ptr = NULL;
}


void Pool_add_block(Pool pool)
{
  char * block;
  char * block_end;

  block = malloc(pool->block_size * pool->size);
  Slist_push(pool->blocks, block);

  block_end = block + pool->block_size * pool->size;
  for ( ; block < block_end ; block += pool->size )
    Slist_push(pool->free, block);
}


int _Pool_chunk_owned(Pool pool, void * chunk)
{
  /* Return true if a chunk is within one of the Pool's blocks
   * Useful for debugging calls to Pool_free, but imposes a large
   * performance hit */

  Slist_node node;
  Slist      blocks;

  blocks = pool->blocks;
  /* Iterate the list of blocks */
  for( node = blocks->first; node != NULL; node = node->next )
  {
    /* If a chunk address falls within the range, return true */
    if ( chunk >= node->data && chunk < node->data + (pool->block_size * pool->size))
    {
      return 1;
    }
  }
  /* else return false */
  return 0;
}


#ifdef POOL_ADT_TEST

typedef int tenint[10];

int main(void)
{
  Pool pool;
  tenint * my_arrays[14];
  int i;

  pool = Pool_new(sizeof(tenint), 5);

  for ( i = 0; i < 14; i++ )
  {
    my_arrays[i] = (tenint *)Pool_alloc(pool);
    assert( my_arrays[i] != NULL );
    (*my_arrays[i])[9] = 437;
  }

  assert( Slist_get_length(pool->blocks) == 3 );
  assert( Slist_get_length(pool->free) == 1 );

  for ( i = 13; i >=0; i-- )
    (*my_arrays[i])[i % 10] = 7*i;

  for ( i = 0; i < 14; i += 2 )
    Pool_free(pool, my_arrays[i]);

  for ( i = 12; i >= 0; i -= 2 )
    my_arrays[i] = (tenint *)Pool_alloc(pool);

  for ( i = 0; i < 14; i++ )
    Pool_free(pool, my_arrays[i]);

  /* This should throw an assertion */
  tenint *not_from_pool = malloc(sizeof(tenint));
  Pool_free(pool, not_from_pool);

  Pool_delete(&pool);
  assert( pool == NULL );
}

#endif /* POOL_ADT_TEST */


// End.
