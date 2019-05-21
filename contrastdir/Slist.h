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

typedef void *Slist_data;

typedef Slist_data(*Slist_data_copy_t)(Slist_data);
typedef void (*Slist_data_delete_t)(Slist_data *);

typedef struct Slist *Slist;
typedef struct Slist_node *Slist_node;
typedef struct Slist_iter *Slist_iter;

/* Slist object */
struct Slist {
  long		length;
  Slist_node	first;
  Slist_node    last;
};

/* Slist node object */
struct Slist_node {
  Slist_node	next;
  Slist_data	data;
};

/* Slist iterator */
struct Slist_iter {
  Slist         l;
  Slist_node    next;
};

/* Create and return a new empty list */
extern Slist Slist_new(void);

/* Create a new list from an array of Data objects. Array must be terminated
 * with a null pointer. */
extern Slist Slist_new_fromarray(Slist_data data[]);

/* These are not implemented. Do not call */
extern Slist Slist_copy(Slist l);
extern Slist Slist_copy_deep(Slist l, Slist_data_copy_t copy_func);

/* Add obj to end of list */
extern void Slist_append(Slist l, Slist_data obj);

/* Add obj to start of list */
extern void Slist_push(Slist l, Slist_data obj);

/* Remove and return first object */
extern Slist_data Slist_pop(Slist l);

/* Get the length of the list */
static __inline__ long Slist_get_length(Slist l)
{
  return l->length;
}

/* Free a list's memory. List must be empty */
extern void Slist_delete(Slist *list_ptr);

/* Return true if list is empty. */
extern int Slist_isempty(Slist l);

/* Delete list and free data objects */
extern void Slist_delete_data(Slist *list_ptr, Slist_data_delete_t delete_func);

Slist_iter Slist_begin(Slist l);

void Si_delete(Slist_iter *iter_ptr);

void * Si_next(Slist_iter iter);

static __inline__ int Si_atend(const Slist_iter iter)
{
  assert( iter != NULL );
  return iter->next == NULL;
}

#endif
/* SLIST_H */
