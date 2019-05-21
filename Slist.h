/* Version 4.0. (c) Copyright 2012 by the University of Washington.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


/* Slist.h
 *
 * Singly-linked list ADT
 *
 * Keeps an ordered list of pointers.
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifndef SLIST_H
#define SLIST_H

typedef void * Slist_data_ptr;

typedef struct _Slist * Slist_ptr;
typedef struct _Slist_node * Slist_node_ptr;

/* Slist object */
struct _Slist {
  long length;
  Slist_node_ptr first;
  Slist_node_ptr last;
};

/* Slist node object */
struct _Slist_node {
  Slist_node_ptr next;
  Slist_data_ptr data;
};

/* Create and return a new empty list */
extern Slist_ptr Slist_new(void);

/* Free a list's memory. List must be empty */
extern void Slist_delete(Slist_ptr list_ptr);

/* Return true if list is empty. */
extern int Slist_isempty(Slist_ptr l);

/* Add obj to start of list */
extern void Slist_push(Slist_ptr l, Slist_data_ptr obj);

/* Remove and return first object */
extern Slist_data_ptr Slist_pop(Slist_ptr l);


/* Get the length of the list */
static __inline__ long Slist_get_length(Slist_ptr l)
{
  return l->length;
}

/* Add obj to end of list */
extern void Slist_append(Slist_ptr l, Slist_data_ptr obj);


#ifdef LIST_ADT_TEST

/* Slist iterator */
struct _Slist_iter {
  Slist_ptr         l;
  Slist_node_ptr    next;
};

typedef struct _Slist_iter * Slist_iter_ptr;
typedef void (*Slist_data_delete_t)(Slist_data_ptr *);
typedef Slist_data_ptr (*Slist_data_copy_t)(Slist_data_ptr);

/* Create a new list from an array of Data objects. Array must be terminated
 * with a null pointer. */
extern Slist_ptr Slist_new_fromarray(Slist_data_ptr data[]);

extern Slist_ptr Slist_copy(Slist_ptr l);
extern Slist_ptr Slist_copy_deep(Slist_ptr l, Slist_data_copy_t copy_func);

/* Delete list and free data objects */
extern void Slist_delete_data(Slist_ptr list_ptr, Slist_data_delete_t delete_func);

Slist_iter_ptr Slist_begin(Slist_ptr l);

void Si_delete(Slist_iter_ptr *iter_ptr);

void * Si_next(Slist_iter_ptr iter);

static __inline__ int Si_atend(const Slist_iter_ptr iter)
{
  assert( iter != NULL );
  return iter->next == NULL;
}

#endif  // LIST_ADT_TEST


#endif
/* SLIST_H */


// End.
