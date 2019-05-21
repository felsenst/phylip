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
#include <stdlib.h>
#include "Slist.h"

/* Define this to include testing function main() */
/* #define LIST_ADT_TEST */

/* Private function prototypes */
void Slist_node_delete(Slist_node *ln_ptr);
int  _Slist_checklen(Slist l);


static Slist_node Slist_node_new_(Slist_data data)
{ /* Slist_node constructor which accepts NULL data */
  Slist_node      node;

  node = (Slist_node)malloc(sizeof(struct Slist_node));
  assert( node != NULL );

  node->data = data;
  node->next = NULL;

  return node;
} /* Slist_node_new_ */


static Slist_node Slist_node_new(Slist_data data)
{ /* get a new node from the linked list */

  assert( data != NULL );
  return Slist_node_new_(data);
} /* Slist_node_new */


void Slist_node_delete(Slist_node *ln_ptr)
{ /* remove a node from the list, and free it */
  Slist_node      ln;

  assert( ln_ptr != NULL );
  
  ln = *ln_ptr;
  assert( ln != NULL );

  free(ln);

  *ln_ptr = NULL;
} /* Slist_node_delete */


Slist Slist_new()
{ /* set up a new linked list, initially empty */
  Slist l;

  l = malloc(sizeof(struct Slist));
  assert( l != NULL );

  l->length = 0;
  l->first = NULL;
  l->last = NULL;

  return l;
} /* Slist_new */


void Slist_delete(Slist *list_ptr)
{ /* delete a linked list */
  Slist          l = NULL;

  assert( list_ptr != NULL );
  l = *list_ptr;
  assert( l != NULL );
  assert( Slist_isempty(l) );

  free( l );

  *list_ptr = NULL;
}

#include <stdio.h>


/* DEBUG function  */
int _Slist_checklen(Slist l) {

  long i;
  Slist_node node;

  i = 0;
  for ( node = l->first; node != NULL; node = node->next )
    i++;

  if ( i == l->length )
    return 1;
  else
    printf("List %p->length is %ld, should be %ld!\n", (void *)l, l->length, i);
    return 0;
}


int Slist_isempty(Slist l)
{  /* check whether the Slist is empty */
  assert( l != NULL );
  return l->length == 0;
} /* Slist_isempty */


void Slist_delete_data(Slist *list_ptr, Slist_data_delete_t delete_func)
{
  Slist          l = *list_ptr;
  Slist_data        obj;

  while ( !Slist_isempty(l) ) {
    obj = Slist_pop(l);
    delete_func(&obj);
  }

  Slist_delete( list_ptr );
} /* Slist_delete */


Slist Slist_new_fromarray(Slist_data obj[])
{
  Slist          list;

  Slist_data        *obj_ptr;

  /* Create new list */
  list = Slist_new();
  
  /* Insert from obj[] until a NULL Slist_data is encountered */
  for ( obj_ptr = obj; *obj_ptr != NULL; obj_ptr++ ) {
    Slist_append(list, *obj_ptr);
  }

  return list;
} /* Slist_new_fromarray */


Slist Slist_copy(Slist l) {
  return Slist_copy_deep(l, NULL);
} /* Slist_copy */


Slist Slist_copy_deep(Slist l, Slist_data_copy_t copy_func)
{
  Slist      dst;
  Slist_iter si;
  Slist_data data;
  
  assert( l != NULL );
  dst = Slist_new();

  si = Slist_begin(l);
  while ( !Si_atend(si) ) {
    data = Si_next(si);
    if ( copy_func )
      data = (*copy_func)(data);
    Slist_append(dst, data);
  }

  Si_delete(&si);

  return dst;
}


void Slist_append(Slist l, Slist_data data)
{
  Slist_node	node;
  
  assert( data != NULL );

  node = Slist_node_new(data);

  if ( l->last == NULL ) {
    l->first = node;
  }
  else {
    l->last->next = node;
  }
  l->last = node;

  l->length++;
}

void Slist_push(Slist l, Slist_data data)
{
  Slist_node node;

  assert( data != NULL );
  node = Slist_node_new(data);

  if ( l->first == NULL )
    l->last = node;
  node->next = l->first;
  l->first = node;

  l->length++;
}


Slist_data Slist_pop(Slist l)
{ /* pop a node off the linked list */
  Slist_node    del;
  Slist_data    data;

  assert( !Slist_isempty(l) );

  del = l->first;
  data = del->data;

  l->first = del->next;
  
  if ( l->first == NULL )
    l->last = NULL;

  Slist_node_delete(&del);

  l->length--;

  return data;
} /* Slist_pop */


Slist_iter Slist_begin(Slist l)
{
  Slist_iter iter;

  iter = (Slist_iter)malloc(sizeof(struct Slist_iter));
  assert( iter != NULL );
  iter->next = l->first;
  return iter;
}


void Si_delete(Slist_iter *iter_ptr)
{
  Slist_iter    iter;

  assert( iter_ptr != NULL );
  iter = *iter_ptr;
  assert( iter != NULL );
  free(iter);
  *iter_ptr = NULL;
}


void * Si_next(Slist_iter iter)
{
  void * data;
 
  assert( iter != NULL );
  if ( iter->next == NULL )
    return NULL;

  data = iter->next->data;
  iter->next = iter->next->next;

  return data;
}


#ifdef LIST_ADT_TEST
#include <stdio.h>

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


int main()
{
  Slist          l, m;
  Slist_iter iter;
  int *data;
  int           n[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  Slist_data        nobj[11];
  int           i;

  for (i = 0; i < 10; i++) {
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
  while ( !Slist_isempty(l) ) {
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
  while ( !Slist_isempty(l) ) {
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
  while( (data = (int *)Si_next(iter)) != NULL ) {
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

#endif

