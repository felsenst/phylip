/* Version 4.0. (c) Copyright 2020-2023.
   */

/* Slist.h
 *
 * Singly-linked list ADT   debug: what is "ADT"?
 *
 * Keeps an ordered list of objects, each of which contains a list
 * item which iself is a structure which contains a pointer to the 
 * next list item plus a pointer to a free node..
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifndef SLIST_H
#define SLIST_H

/* debug: #include "phylip.h" Slist.h is included in phylip.h!  */

typedef struct _Slist_node _Slist_node;              /* forward declaration */
typedef struct _Slist_node* Slist_node_ptr;
typedef struct node* Slist_data_ptr;


/* Slist node object */
struct _Slist_node {
  Slist_node_ptr next;
  Slist_data_ptr data;
};

/* Slist object */
struct _Slist {
  long length;
  Slist_node_ptr first;
  Slist_node_ptr last;
};

typedef struct _Slist* Slist_ptr;

/* Get the length of the list */
static __inline__ long Slist_get_length(Slist_ptr l)
{
  return l->length;
}

#ifdef LIST_ADT_TEST

/* Slist iterator */
struct _Slist_iter {
  Slist_ptr         l;
  Slist_node_ptr    next;
};

struct _Slist_iter* Slist_iter_ptr;

typedef Slist_data_ptr (*Slist_data_copy_t)(Slist_data_ptr);
typedef void (*Slist_data_delete_t)(Slist_data_ptr *);

/* Create a new list from an array of Data objects. Array must be terminated
 * with a null pointer. */
Slist_ptr Slist_new_fromarray(Slist_data_ptr data[]);

Slist_ptr Slist_copy(Slist_ptr l);
Slist_ptr Slist_copy_deep(Slist_ptr l, Slist_data_copy_t copy_func);

Slist_iter_ptr Slist_begin(Slist_ptr l);

void Si_delete(Slist_iter_ptr *iter_ptr);

void * Si_next(Slist_iter_ptr iter);

static __inline__ int Si_atend(const Slist_iter_ptr iter)
{
  assert( iter != NULL );
  return iter->next == NULL;
}

#endif  // LIST_ADT_TEST

typedef struct _Slist_iter* Slist_iter_ptr;
typedef Slist_data_ptr (*Slist_data_copy_t)(Slist_data_ptr);
typedef void (*Slist_data_delete_t)(Slist_data_ptr *);

#ifndef OLDC  /* if not old original K&R C */
/* prototypes */
Slist_node_ptr Slist_node_new_(Slist_data_ptr data);  /* note diff in names */
Slist_node_ptr Slist_node_new(Slist_data_ptr data);   /* note diff in names */
void Slist_node_delete(Slist_node_ptr ln);
Slist_ptr Slist_new(void);
void Slist_delete(Slist_ptr l);
int _Slist_checklen(Slist_ptr l);
int Slist_isempty(Slist_ptr l);
void Slist_push(Slist_ptr l, Slist_data_ptr data);
Slist_data_ptr Slist_pop(Slist_ptr l);
void Slist_append(Slist_ptr l, Slist_data_ptr data);
void Slist_delete_data(Slist_ptr l, Slist_data_delete_t delete_func);
Slist_ptr Slist_new_fromarray(Slist_data_ptr obj[]);
Slist_ptr Slist_copy(Slist_ptr l);
Slist_ptr Slist_copy_deep(Slist_ptr l, Slist_data_copy_t copy_func);
Slist_iter_ptr Slist_begin(Slist_ptr l);
void * Si_next(Slist_iter_ptr iter);
void nobj_delete(void **nobj_ptr);
void *nobj_copy(void *nobj);
#endif

#endif
/* SLIST_H */


/* End. */
