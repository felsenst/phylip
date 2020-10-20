/* Version 4.0. (c) Copyright 2013 by the University of Washington.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


/* Slist.c
 *
 * Implementation of singly-linked list
 *
 * by Ian Robertson, 2006
 */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "Slist.h"


/* Define this to include testing function main() */
/* #define LIST_ADT_TEST */

/* Private function prototypes */
void Slist_node_delete(Slist_node_ptr ln);
int  _Slist_checklen(Slist_ptr l);


static Slist_node_ptr Slist_node_new_(Slist_data_ptr data)
{ /* Slist_node constructor which can even accept NULL data,
   * called by Slist_node_new (note difference in names) */
  Slist_node_ptr      node;

  node = (Slist_node_ptr)malloc(sizeof(struct _Slist_node));
  assert( node != NULL );

  node->data = data;     /* have its data pointer point to "data" argument */
  node->next = NULL;                                   /* nothing after it */

  return node;
} /* Slist_node_new_ */


static Slist_node_ptr Slist_node_new(Slist_data_ptr data)
{ /* get a new node from the linked list */
  assert( data != NULL );
  return Slist_node_new_(data);
} /* Slist_node_new */


void Slist_node_delete(Slist_node_ptr ln)
{ /* remove a node from the list, and free it */
  assert( ln != NULL );
  free(ln);
} /* Slist_node_delete */


Slist_ptr Slist_new(void)
{ /* set up a new linked list, initially empty */
  Slist_ptr l;

  l = malloc(sizeof(struct _Slist));
  assert( l != NULL );

  l->length = 0;
  l->first = NULL;
  l->last = NULL;

  return l;
} /* Slist_new */


void Slist_delete(Slist_ptr l)
{ /* delete a linked list */
  assert( l != NULL );
  assert( Slist_isempty(l) );
  free( l );
} /* Slist_delete */


/* DEBUG function  */
int _Slist_checklen(Slist_ptr l)
{
  long i;
  Slist_node_ptr node;

  i = 0;
  for ( node = l->first; node != NULL; node = node->next )
    i++;

  if ( i == l->length )
    return 1;
  else
    printf("List %p->length is %ld, should be %ld!\n", (void *)l, l->length, i);
  return 0;
}


int Slist_isempty(Slist_ptr l)
{  /* check whether the Slist is empty */
  assert( l != NULL );
  return l->length == 0;
} /* Slist_isempty */


void Slist_push(Slist_ptr l, Slist_data_ptr data)
{
  /* make a new list-node and put it on "last"
  Slist_node_ptr node;

  assert( data != NULL );
  node = Slist_node_new(data); /* make new list-node which points to "data" */

  if ( l->first == NULL )                  /* if there's nobody on the list */
  {
    assert(l->last == NULL);
    l->last = node;          /* then have "last" point to the new list-node */
  }

  node->next = l->first;    /* have new list-node point to previous "first" */
  l->first = node;          /* ... and have "first" point to it */

  l->length++;
} /* Slist_push */


Slist_data_ptr Slist_pop(Slist_ptr l)
{ /* pop a node off the linked list */
  Slist_node_ptr    del;
  Slist_data_ptr    data;

  assert( !Slist_isempty(l) );

  del = l->first;                      /* the list-node that was "on first" */
  data = del->data;                        /* the "data" that it pointed to */

  l->first = del->next;  /* make "first" point it the list-node's successor */

  if ( l->first == NULL )
    l->last = NULL;

  Slist_node_delete(del);         /* toss the list-node that was popped off */

  l->length--;

  return data;  /* return the "data" from the list-node that was popped off */
} /* Slist_pop */


void Slist_append(Slist_ptr l, Slist_data_ptr data)
{
  Slist_node_ptr node;

  assert( data != NULL );
  node = Slist_node_new(data);

  if ( l->last == NULL )
  {
    assert( l->first == NULL);
    l->first = node;
  }
  else
  {
    l->last->next = node;
  }
  l->last = node;

  l->length++;
} /* Slist_append */


#ifdef LIST_ADT_TEST
#include <stdio.h>


void Slist_delete_data(Slist_ptr l, Slist_data_delete_t delete_func)
{
  Slist_data_ptr        obj;

  while ( !Slist_isempty(l) )
  {
    obj = Slist_pop(l);
    delete_func(&obj);
  }

} /* Slist_delete_data */


Slist_ptr Slist_new_fromarray(Slist_data_ptr obj[])
{
  Slist_ptr             list;

  Slist_data_ptr        *obj_ptr;

  /* Create new list */
  list = Slist_new();

  /* Insert from obj[] until a NULL Slist_data is encountered */
  for ( obj_ptr = obj; *obj_ptr != NULL; obj_ptr++ )
  {
    Slist_append(list, *obj_ptr);
  }

  return list;
} /* Slist_new_fromarray */


Slist_ptr Slist_copy(Slist_ptr l)
{
  return Slist_copy_deep(l, NULL);
} /* Slist_copy */


Slist_ptr Slist_copy_deep(Slist_ptr l, Slist_data_copy_t copy_func)
{
  Slist_ptr      dst;
  Slist_iter_ptr si;
  Slist_data_ptr data;

  assert( l != NULL );
  dst = Slist_new();

  si = Slist_begin(l);
  while ( !Si_atend(si) )
  {
    data = Si_next(si);
    if ( copy_func )
      data = (*copy_func)(data);
    Slist_append(dst, data);
  }

  Si_delete(&si);

  return dst;
}


Slist_iter_ptr Slist_begin(Slist_ptr l)
{
  Slist_iter_ptr iter;

  iter = (Slist_iter_ptr)malloc(sizeof(struct _Slist_iter));
  assert( iter != NULL );
  iter->next = l->first;
  return iter;
}


void Si_delete(Slist_iter_ptr *iter_ptr)
{
  Slist_iter_ptr    iter;

  assert( iter_ptr != NULL );
  iter = *iter_ptr;
  assert( iter != NULL );
  free(iter);
  *iter_ptr = NULL;
}


void * Si_next(Slist_iter_ptr iter)
{
  void * data;

  assert( iter != NULL );
  if ( iter->next == NULL )
    return NULL;

  data = iter->next->data;
  iter->next = iter->next->next;

  return data;
}


/* Example object destructor */
void nobj_delete(void **nobj_ptr)
{
  printf("Deleting %d.\n", **(int **)nobj_ptr);
  assert( nobj_ptr );
  assert( *nobj_ptr );
  free(*nobj_ptr);
  *nobj_ptr = NULL;
}


/* Example object copy constructor */
void *nobj_copy(void *nobj)
{
  int *dst;

  printf("Copying %d\n", *(int *)nobj);
  dst = (int *)malloc(sizeof(int));
  *dst = *(int *)nobj;

  return (void *)dst;
}


int main(void)
{
  Slist          l, m;
  Slist_iter iter;
  int *data;
  int           n[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  Slist_data_ptr        nobj[11];
  int           i;

  for (i = 0; i < 10; i++)
  {
    nobj[i] = &n[i];
  }
  nobj[10] = NULL;

  printf("Create new Slist...\n");
  l = Slist_new();
  printf("Slist pointer is %p\n", l);
  printf("Slist is %sempty\n", Slist_isempty(l) ? "" : "not ");

  printf("Adding 5, 6 at end...\n");
  Slist_append(l, nobj[5]);
  Slist_append(l, nobj[6]);
  printf("Adding 4, 3 at start...\n");
  Slist_push(l, nobj[4]);
  Slist_push(l, nobj[3]);

  printf("Slist = ");
  while ( !Slist_isempty(l) )
  {
    printf("%d ", *((int *)Slist_pop(l)));
  }
  printf("\n");

  printf("Deleting list...\n");
  Slist_delete(&l);

  printf("Slist pointer is %p\n", l);

  printf("Creating new list from array...\n");
  l = Slist_new_fromarray(nobj);
  printf("Slist pointer is %p\n", l);
  printf("Slist = ");
  while ( !Slist_isempty(l) )
  {
    printf("%d ", *((int *)Slist_pop(l)));
  }
  Slist_delete(&l);
  printf("\n");
  printf("Creating again...\n");
  l = Slist_new_fromarray(nobj);


  printf("Copying l to m...\n");
  m = Slist_copy_deep(l, nobj_copy);

  printf("Deleting l...\n");
  while( !Slist_isempty(l) )
    Slist_pop(l);
  Slist_delete(&l);

  printf("Iterating list m...\n");
  iter = Slist_begin(m);
  while( (data = (int *)Si_next(iter)) != NULL )
  {
    if ( Si_atend(iter) )
      printf("and finally ");
    printf("%d ", *((int *)data));
  }
  printf("\n");
  assert( Si_atend(iter) );
  Si_delete(&iter);


  printf("Deleting m with object destructor...\n");
  Slist_delete_data(&m, nobj_delete);
  printf("Slist pointer is %p\n", l);
  printf("Done.\n");

  return 0;
}

#endif  // LIST_ADT_TEST


// End.
