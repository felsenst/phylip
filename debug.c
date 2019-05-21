/* Version 4.0. (c) Copyright 2012-2013 by the University of Washington.
   Written by Bob Giansiracusa for debugging purposes.
   This file contains declarations used only for debugging, not to be included in distributed systems. */


#include "phylip.h"
#include "debug.h"


// If TRUE, every call to TREECHECKER prints diagnostic information, whether a problem is detected or not.
// If FALSE, only the call to TREECHECKER which diagnosed a problem prints, and then it immediately ABORTS.
#define PRINT_ALWAYS false


// Global variable by which COUNT_EM can return an ERRORCODE (in addition to a count) if needed.
static int global_errorcode = 0;


// Prototypes for functions defined and called only in this file.
//
int treeprinter(tree * pt, FILE * stream);
void clearintreebits(Slist_ptr forklist, boolean isringtype);
void tree_traverse(node * pn);
int count_em(node * pn);
int forklistprinter(Slist_ptr forklist, node * prootnode, int nextNodeState, boolean isringtype, FILE * stream);
int nodeloopprinter(const char * label, node * pq, node * prootnode, int nextNodeState, FILE * stream);
int nodegutsprinter(const char * label, node * pq, node * prootnode, int nextNodeState, FILE * stream);
//
// End of prototypes.


// For manual calling in GDB.  Prints to stdout (terminal window).
//
void seenode(node * pq)
{
  if (pq == NULL)
  {
    fprintf(stdout, "\n  Node pointer is NULL.\n\n");
  }
  else
  {
    fprintf(stdout, "\n  Node at %p is ", (void *) pq);

    // Unknown NEXT pointer status when printing an arbitrary node; test on type (TIP vs FORK).
    // Can't check validity of anything because we have no TREE to traverse for INTREE nodes.  Could be printing an INTREE node, a non-INTREE (in nodep[] but not in a TREE), or a garbage node.
    (void)nodegutsprinter("NODE", pq, (node *)0, 0, stdout); // Will always report node as non-ROOT because pointer-to-root argument is NULL here.
    putc('\n', stdout);
  }
}


// For manual calling in GDB.  Prints to stdout (terminal window).
//
void seetree(tree * pt)
{
  (void)treeprinter(pt, stdout);
}


// For calling in source program.  Prints to STDERR ("Error.Output" file).  Prints on every call until an error is detected, and then it aborts.
//
void treechecker(tree * pt, const char * filename, const unsigned int linenum, const char * funcname)
{
  int errorcode;
  static unsigned int callcount = 0;

  ++callcount;                          // Total number of calls to TREECHECKER from any location.

  if ((pt == NULL) || pt->root)         // OK to print or check if pt is NULL or if pt is valid and pt->root is valid, but don't call if pt->root is NULL.
  {
#if PRINT_ALWAYS
    // Label the printout for the upcoming batch of text.  This gets printed before COUNT_EM gets called, which might print an error message, in which case this line gets re-printed with the ErrorCode.
    fprintf(stderr, "Calls: %u, File: \"%s\", Line: %u, Function: \"%s\"\n", callcount, filename, linenum, funcname);
#endif

    // Calling TREEPRINTER with second arg NULL causes it (and all its callees) to CHECK for errors and REPORT (by returning an ERRORCODE > 0) but to not print anything,
    // except that COUNT_EM prints messages on encountering errors, because it is only called in TREEPRINTER's TESTING phase.
    // TREEPRINTER returns an ERRORCODE > 0 (when errors are detected) only in TESTING mode (when STREAM is NULL).
    // Next line expressions evaluated left-to-right, so GLOBAL_ERRORCODE will contain value set by COUNT_EM (if it gets set) after call to TREEPRINTER returns.
    errorcode = treeprinter(pt, NULL) + global_errorcode;

    // Now reprint any non-zero ERRORCODE returns, so that other programs can GREP for these codes.
    if (errorcode > 0)
      fprintf(stderr, "Calls: %u, File: \"%s\", Line: %u, Function: \"%s\", ErrorCode: %d.\n", callcount, filename, linenum, funcname, errorcode);

#if PRINT_ALWAYS
    (void)treeprinter(pt, stderr);      // And call TREEPRINTER again, this time in PRINTING mode, to do the full printing job.
#endif

    if (errorcode > 0)                  // If ANY error is detected, ABORT via ASSERT (so debugger can report BackTrace).
    {
#if ! PRINT_ALWAYS
    (void)treeprinter(pt, stderr);      // And call TREEPRINTER again, this time in PRINTING mode, to do the full printing job.
#endif
      assert(false);
    }
  }
}


int treeprinter(tree * pt, FILE * stream)
{
  node * pq, * pp;
  int nodepidx, tnonodes, tspp, errorcode, itercount = 0;

  if (pt == NULL)
  {
    if (stream)
      fprintf(stream, "The argument is a NULL tree.\n");
  }
  else
  {
    tnonodes = (int)pt->nonodes;
    tspp = (int)pt->spp;

    if (stream)
      fprintf(stream, "Tree at %p, spp: %d, nonodes: %d, t.type: %d, t.nodep: %p, t.root: %p, t.free_forkrings: %p, t.free_forknodes: %p\n",
              (void *) pt, tspp, tnonodes, pt->type, (void *)pt->nodep, (void *)pt->root, (void *)pt->free_forkrings, (void *)pt->free_forknodes);

    // Print the "Nodes in use" first, so we can see if any nodes on the garbage lists point into "in use" node chains.
    if (tnonodes == 0)
    {
      if (stream)
        fprintf(stream, "  NODEP array is size ZERO.\n");
    }
    else
    {
      assert(tnonodes > 0);             // Consistency check; MUST croak if this is not satisfied.
      assert(tspp > 0);                 // Consistency check; MUST croak if this is not satisfied.
      assert(tspp < tnonodes);          // Consistency check; MUST croak if this is not satisfied.

      if (stream == NULL)               // In TESTING mode, clear INTREE bits for all nodes and then set bit for those nodes IN THE TREE.
      {
        for (nodepidx = 0; nodepidx < tnonodes; ++nodepidx) // For each TIPnode or FORKnode (in a ring) in the NODEP array, CLEAR its INTREE bit.
        {
          pq = pt->nodep[nodepidx];
          if (pq)                       // If the array slot is non-NULL,
          {
            pq->type = NODE_T_UNKNOWN;  // Clear its INTREE bit (using MODE_T_UNKNOWN in TYPE slot as "not-INTREE" indication; slot is otherwise unused).
            if (nodepidx >= spp)        // If the node is a FORKnode,
            {
              pp = pq->next;             // Set PP to its successor.
              while (pp && (pp != pq))   // Iterate (as long as node pointer is non-NULL) until we return to starting point.
              {
                pp->type = NODE_T_UNKNOWN; // Clear INTREE bit for rest of FORKnodes in the FORKring.
                pp = pp->next;
                if (++itercount >= 64)  // Protection against infinite loop if the FORKring is invalid (ie, closes not at starting node).
                  break;                // Assume we have cleared all the nodes (multiple times) and break out of this loop.
              }
            }
          }
        }

        // Must also reset INTREE bits for all nodes on FREE_FORKNODES and FREE_FORKRINGS lists, because some of them may not be in the NODEP array.
        clearintreebits(pt->free_forknodes, false);
        clearintreebits(pt->free_forkrings, true);

        // Now we are finished CLEARING all the INTREE bits for all nodes.  Traverse the tree, MARKING all reachable nodes as being INTREE.
        if (pt->root->tip)
        {
          if (pt->root->back)              // If ROOT node (TIP, in this case) has non-NULL BACK, traverse from there; it may have a NULL BACK.
            tree_traverse(pt->root->back); // Then traverse the tree recursively, starting from the FORKnode pointing to this TIPnode, marking all nodes which ARE in the tree as being INTREE nodes.
        }
        else
        {
          tree_traverse(pt->root);      // Or traverse the tree recursively from the ROOT node, starting from this FORKnode, marking all nodes which ARE in the tree as being INTREE nodes.
        }
      }

      // Done with garbage-handling and now back to tree printing.
      if (stream)
        fprintf(stream, "  NODEP array (size %d):\n", tnonodes);

      for (nodepidx = 0; nodepidx < tnonodes; ++nodepidx)
      {
        pq = pt->nodep[nodepidx];

        if (pq == NULL)
        {
          if (stream)
            fprintf(stream, "    NODEP[%d] has a NULL node pointer.\n", nodepidx);
        }
        else
        {
          if ((int)pq->index != (nodepidx + 1))
          {
            if (stream)
              fprintf(stream, "======> DEBUG ERROR (001):  A NODE's INDEX (%d) does not match nodep[] slot index (%d) + 1.\n", (int)pq->index, nodepidx);
            else
              return(1);
          }

          if (nodepidx < tspp)          // TIPnodes.
          {
            if (stream)
              fprintf(stream, "    NODEP[%d] is %p -> ", nodepidx, (void *)pq);
            if (! pq->tip)
            {
              if (stream)
                fprintf(stream, "======> DEBUG ERROR (002):  A TIPnode (via index) is mistakenly marked as a FORKnode (via flag).\n");
              else
                return(2);
            }

            // TIPnodes (INTREE or not) must have NULL NEXT pointers.
            if (pq->type == NODE_T_UNKNOWN) // Non-INTREE nodes (in nodep[] array but not reachable in a tree traversal from ROOT).
            {
              errorcode = nodeloopprinter("non-INTREE", pq, pt->root, 1, stream); // BACK is arbitary (may be partially initialized and not yet on garbage list).
              if (errorcode > 0)
                return(errorcode);      // A non-INTREE TIPnode had a non-NULL NEXT pointer.
            }
            else if (pq == pt->root)    // INTREE and ROOT TIPnode.
            {
              if (pq->back)
              {
                errorcode = nodeloopprinter("yes-INTREE", pq, pt->root, 1, stream);
                if (errorcode > 0)
                  return(errorcode);
              }
              else                      // If BACK is NULL, it must be NULL (obviously; we just need this branch to print the node).
              {
                errorcode = nodeloopprinter("yes-INTREE", pq, pt->root, 1, stream);
                if (errorcode > 0)
                  return(errorcode);
              }
            }
            else                        // INTREE and non-ROOT TIPnode; must have mutually-reflective BACK state.
            {
              errorcode = nodeloopprinter("yes-INTREE", pq, pt->root, 1, stream);
              if (errorcode > 0)
                return(errorcode);
            }
          }
          else                          // FORKnodes in a FORKring.
          {
            if (stream)
              fprintf(stream, "    NODEP[%d] is %p ->:\n", nodepidx, (void *)pq);
            if (pq->tip)
            {
              if (stream)
                fprintf(stream, "======> DEBUG ERROR (003):  A FORKnode (via index) is mistakenly marked as a TIPnode (via flag).\n");
              else
                return(3);
            }
            // FORKnodes (INTREE) must have valid NEXT pointers.
            // FORKnodes (on the FREE_FORKNODES list) must have NULL NEXT pointers.
            // FORKnodes (on the FREE_FORKRINGS list) must have valid NEXT pointers.
            // Assumption is that all Non-INTREE nodes are on the FREE_FORKNODES list or on the FREE_FORKRINGS list and all nodes on those lists are Non-INTREE.
            // This assumption is checked by tests in FORKLISTTESTER and FORKLISTPRINTER.
            if (pq->type == NODE_T_UNKNOWN) // Non-INTREE nodes (in nodep[] array but not reachable in a tree traversal from ROOT).
            {
              errorcode = nodeloopprinter("      non-INTREE", pq, pt->root, 0, stream);
              if (errorcode > 0)
                return(errorcode);      // An non-INTREE FORKnode not reachable in tree traversal from ROOT.
            }
            else if (pq == pt->root)    // INTREE and ROOT FORKnode.
            {
              if (pq->back)             // If BACK is non-NULL, it must be mutually-reflective (point to its BACK's BACK).
              {
                errorcode = nodeloopprinter("      yes-INTREE", pq, pt->root, 2, stream);
                if (errorcode > 0)
                  return(errorcode);
              }
              else                      // If BACK is NULL, it must be NULL (obviously; we just need this branch to print the node).
              {
                errorcode = nodeloopprinter("      yes-INTREE", pq, pt->root, 2, stream);
                if (errorcode > 0)
                  return(errorcode);
              }
            }
            else                        // INTREE and non-ROOT FORKnode; must have mutually-reflective BACK state.
            {
              errorcode = nodeloopprinter("      yes-INTREE", pq, pt->root, 2, stream);
              if (errorcode > 0)
                return(errorcode);
            }
          }
        }
      }

      if (pt->free_forkrings == NULL && pt->free_forknodes == NULL)
      {
        if (stream)
        {
          fprintf(stream, "  FREE_FORKNODES: Empty.\n");
          fprintf(stream, "  FREE_FORKRINGS: Empty.\n");
        }
      }
      else if (pt->free_forknodes == NULL)
      {
        if (stream)
        {
          fprintf(stream, "  FREE_FORKNODES: Empty.\n");
          fprintf(stream, "  FREE_FORKRINGS (length %d):\n", (int)pt->free_forkrings->length);
        }
        errorcode = forklistprinter(pt->free_forkrings, pt->root, 2, true, stream);
        if (errorcode > 0)
          return(errorcode);            // Garbage FORKnodes (on the FREE_FORKRINGS list) must have non-NULL NEXT pointers.
      }
      else if (pt->free_forkrings == NULL)
      {
        if (stream)
          fprintf(stream, "  FREE_FORKNODES (length %d):\n", (int)pt->free_forknodes->length);
        errorcode = forklistprinter(pt->free_forknodes, pt->root, 1, false, stream);
        if (errorcode > 0)
          return(errorcode);            // Garbage FORKnodes (on the FREE_FORKNODES list) must have NULL NEXT pointers.
        if (stream)
          fprintf(stream, "  FREE_FORKRINGS: Empty.\n");
      }
      else
      {
        if (stream)
          fprintf(stream, "  FREE_FORKNODES (length %d):\n", (int)pt->free_forknodes->length);
        errorcode = forklistprinter(pt->free_forknodes, pt->root, 1, false, stream);
        if (errorcode > 0)
          return(errorcode);            // Garbage FORKnodes (on the FREE_FORKNODES list) must have NULL NEXT pointers.
        if (stream)
          fprintf(stream, "  FREE_FORKRINGS (length %d):\n", (int)pt->free_forkrings->length);
        errorcode = forklistprinter(pt->free_forkrings, pt->root, 2, true, stream);
        if (errorcode > 0)
          return(errorcode);            // Garbage FORKnodes (on the FREE_FORKRINGS list) must have non-NULL NEXT pointers.
      }
    }
  }

  if (stream)
    putc('\n', stream);

  return(0);
}


void clearintreebits(Slist_ptr forklist, boolean isringtype)
{
  int listidx, listlen = (int)forklist->length;
  int itercount = 0;
  Slist_node_ptr snp = forklist->first;

  for (listidx = 0; listidx < listlen; ++listidx) // "listidx" is iterating through the Slist here, not through the nodep[] array.
  {
    node * pq = (node *)snp->data;
    assert(pq);

    pq->type = NODE_T_UNKNOWN;          // Clear its INTREE bit (using MODE_T_UNKNOWN in TYPE slot as "not-INTREE" indication; slot is otherwise unused).
    if (isringtype)                     // If the FORKnode is on the FREE_FORKRINGS list (ie, if FORK is of type RING rather than of type NODE) ...
    {
      node * pp = pq->next;             // Set PP to its successor.
      while (pp && (pp != pq))          // Iterate (as long as node pointer is non-NULL) until we return to starting point.
      {
        pp->type = NODE_T_UNKNOWN;      // Clear INTREE bit for rest of FORKnodes in the FORKring.
        pp = pp->next;
        if (++itercount >= 64)          // Protection against infinite loop if the FORKring is invalid (ie, closes not at starting node).
          break;                        // Assume we have cleared all the nodes (multiple times) and break out of this loop.
      }
    }

    snp = snp->next;                    // Go to next data container on list (whose snp->data component is a FORKnode).
  }
}


void tree_traverse(node * pn)
{
  // Traverse tree recursively from root, marking all INTREE nodes (those nodes in nodep[] array that are in the tree) using NODETYPE slot as INTREE bit (it is otherwised unused).
  int sibidx, num_sibs;

  if (! pn)                             // If passed a NULL (tree may not have a root), return.
    return;

  if (pn->type == NODE_T_GENERIC)       // If this node has already been marked as INTREE, return (prevents infinite recursion).
    return;

  pn->type = NODE_T_GENERIC;            // Mark node (TIP or FORK) as being INTREE, not a garbage node or a disconnected chunk.

  if (! pn->tip)                        // If a FORKnode, count its siblings and mark them, recursively traversing each and its BACKee.
  {
    if (pn->back)                       // Recursively mark the FORKnode's BACK (if it is either a TIPnode or another FORKring; could be NULL).
      tree_traverse(pn->back);

    num_sibs = count_em(pn);
    for (sibidx = 0; sibidx < num_sibs; ++sibidx)
    {
      pn = pn->next;                    // Mark starting from startnode's NEXT and ending with STARTNODE.
      if (pn && pn->back)
      {
        pn->type = NODE_T_GENERIC;       // Mark FORKnode as being INTREE, not a garbage node or a disconnected chunk.
        tree_traverse(pn->back);         // Recursively mark all other nodes reachable from this FORKnode.
      }
    }
  }
}


int count_em(node * pq)
{                                       // Count the number of nodes in a ring, and return the total number of nodes excluding the one passed into the function (siblings).
  node * pn, * pp;
  int cnt = 0;

  assert(pq);                           // Make sure argument is a node, not a NULL pointer (should never happen).
  assert(! pq->tip);                    // Make sure argument node is a FORKnode, not a TIPnode (should never happen).

  pn = pq->next;

  while (pn != pq)
  {
    if (! pn)                           // Found a FORKring ending on a NULL NEXT; fatal error; abort immediately.
    {
      fprintf(stderr, "======> DEBUG ERROR (400):  A FORKnode has a NULL NEXT pointer in COUNT_EM.\n");
      global_errorcode = 400;           // Return global indication that COUNT_EM found a problem.
      return(cnt);
    }

    if (pn->next != pq)
    {
      pp = pq;
      do
      {
        if (pp == pn->next)             // Found a FORKring which closes OTHER than at the starting point; fatal error; abort immediately.
        {
          fprintf(stderr, "======> DEBUG ERROR (500):  A FORKnode ring terminates elsewhere than at starting node in COUNT_EM.\n\n");
          global_errorcode = 500;       // Return global indication that COUNT_EM found a problem.
          return(cnt);
        }
        pp = pp->next;
      } while (pp != pn);               // Keep testing until "pp" catches up with "pn" FROM BEHIND.
    }

    // Infinite loop detection; all datasets have less than 64 tips.  This detects case where pointer chains wander all over without ever converging on an end point.
    if (++cnt >= 64)
    {
      fprintf(stderr, "======> DEBUG ERROR (600):  Hit limit without ring closure in COUNT_EM.\n\n");
      global_errorcode = 600;           // Return global indication that COUNT_EM found a problem.
      return(cnt);
    }

    pn = pn->next;
  }

  global_errorcode = 0;                 // Return global indication that COUNT_EM found no problems.
  return(cnt);                          // Return number of sibs (FORKnodes not counting starting one) counted so far.
}


// Called ONLY on GARBAGE FORKnodes (on the FREE_FORKNODES list or on the FREE_FORKRINGS list), not INTREE nodes (in some tree).
// nextNodeState:  1 ==> must be NULL; 2 ==> must be non-NULL.
//
int forklistprinter(Slist_ptr forklist, node * prootnode, int nextNodeState, boolean isringtype, FILE * stream)
{
  int listidx, errorcode, listlen = (int)forklist->length;
  Slist_node_ptr snp = forklist->first;
  node * pq;

  for (listidx = 0; listidx < listlen; ++listidx) // "listidx" is iterating through the Slist here, not through the nodep[] array.
  {
    if (snp->data == NULL)
    {
      if (stream)
        fprintf(stream, "======> DEBUG ERROR (007):  A ListItem.data is NULL.\n");
      return(stream ? 0 : 7);           // Fatal error; cannot continue if printing.
    }

    pq = (node *)snp->data;

    // Garbage nodes (FORKnodes only, no TIPnodes) can have either NULL NEXT pointers (if they are on the FREE_FORKNODES list)
    // or valid NEXT pointers (if they are on the FREE_FORKRINGS list), and their BACK nodes are arbitrary (will be initialized when re-used).
    if (stream)
      fprintf(stream, "    FREELIST[%d] item at %p is:\n", listidx, (void *)pq);

    if (pq->type == NODE_T_GENERIC)     // Found an INTREE node (in nodep[] array and reachable in a tree traversal from ROOT) on a Garbage list.
    {
      if (stream)
      {
        (void)nodeloopprinter("      yes-INTREE", pq, prootnode, nextNodeState, stream);
        fprintf(stream, "======> DEBUG ERROR (008):  A garbage FORKnode on the %s list is INTREE.\n", isringtype ? "FREE_FORKRINGS" : "FREE_FORKNODES");
      }
      else
        return(8);
    }
    else                                // Non-INTREE (on the FREE_FORKNODES list or on the FREE_FORKRINGS list).
    {
      errorcode = nodeloopprinter("      non-INTREE", pq, prootnode, nextNodeState, stream);
      if (errorcode > 0)
        return(errorcode);              // If FORKring is bad, report it.
    }

    if ((listidx == listlen - 1) && (snp != forklist->last))
    {
      if (stream)
        fprintf(stream, "======> DEBUG ERROR (009):  The FORKnode list length is wrong.\n");
      return(stream ? 0 : 9);           // Fatal error; cannot continue if printing.
    }

    snp = snp->next;                    // Go to next data container on list (whose snp->data component is a FORKnode).
  }

  return(0);
}


// nextNodeState:  1 ==> must be NULL; 2 ==> must be non-NULL.
//
int nodeloopprinter(const char * label, node * pq, node * prootnode, int nextNodeState, FILE * stream)
{
  node * pp = pq;                       // Starting node; should also be terminating node for ring.
  int errorcode;

  if (pq == NULL)
  {
    if (stream)
      fprintf(stream, "======> DEBUG ERROR (010):  A NODE is NULL.\n");
    return(stream ? 0 : 10);            // Fatal error; cannot continue if printing.
  }

  do
  {
    if (pq == NULL)                     // NEXT pointer was NULL when it should not have been.
    {
      if (stream)
      {
        // This message is a bit redundant, because on the previous iteration, NODEGUTSPRINTER already printed a message about this FORKnode having a NULL NEXT pointer.
        fprintf(stream, "======> DEBUG ERROR (011):  A FORKnode's NEXT pointer is NULL (in a FORKring).\n");
      }
      return(stream ? 0 : 11);          // Fatal error; cannot continue if printing.
    }

    if (pq == pq->next)                 // Make sure next pointer in loop does not loop immediately to self.
    {
      if (stream)
        fprintf(stream, "======> DEBUG ERROR (012):  A FORKnode's NEXT pointer points to Self.\n");
      return(stream ? 0 : 12);          // Fatal error; cannot continue if printing.
    }

    errorcode = nodegutsprinter(label, pq, prootnode, nextNodeState, stream);
    if (errorcode > 0)
      return(errorcode);

    if (pq->tip)                        // If node is TIPnode, with a single iteration we are done.
    {                                   // If node was a TIPnode, NODEGUTSPRINTER tested its NEXT node state already.
      break;
    }

    if ((nextNodeState == 1) && (pq->next == NULL)) // If NEXT is supposed to be NULL (and it is), the node must be a FORKnode on the FREE_FORKNODES list; break out of the iteration.
    {
      break;
    }

    pq = pq->next;                      // Go on to next node in the loop.
  } while (pq != pp);                   // If node is a FORKnode in a ring, follow ring until loop closes.

  return(0);
}


// nextNodeState:  0 ==> no test (unknown constraint; anything is valid); 1 ==> must be NULL; 2 ==> must be non-NULL.
//
int nodegutsprinter(const char * label, node * pq, node * prootnode, int nextNodeState, FILE * stream)
{
  if (pq == NULL)
  {
    if (stream)
      fprintf(stream, "======> DEBUG ERROR (013):  A NODE-pointer is NULL.\n");
    return(stream ? 0 : 13);            // Fatal error; cannot continue if printing.
  }

  if (stream)
  {
    if (pq->tip)
    {
      fprintf(stream, "%s, %s, TIP at %p", label, (pq == prootnode) ? "yes-ROOT" : "non-ROOT", (void *)pq);
    }
    else
    {
      fprintf(stream, "%s, %s, FORK at %p", label, (pq == prootnode) ? "yes-ROOT" : "non-ROOT", (void *)pq);
    }
    fprintf(stream, ", NODE.next: %p, NODE.index: %d, NODE.init: %d, NODE.back: %p", (void *)pq->next, (int)pq->index, pq->initialized, (void *)pq->back);
    if (pq->back)
    {
      fprintf(stream, ", BACK.back: %p, BACK.next: %p, BACK.index: %d, BACK.init: %d", (void *)pq->back->back, (void *)pq->back->next, (int)pq->back->index, pq->back->initialized);
    }
    putc('\n', stream);
  }

  // Assumption: NEXT pointer either should always be NULL (for FORKnodes on the FREE_FORKNODES list or for TIPnodes)
  // or should always be non-NULL (for FORKnodes on the FREE_FORKRINGS list or for FORKnodes in a NODEP array).
  if ((nextNodeState == 1) && pq->next)  // TIPnode or FORKnode on the FREE_FORKNODES list; NEXT should be NULL.
  {
    if (stream)
      fprintf(stream, "======> DEBUG ERROR (014):  A TIPnode or garbage FORKnode (on the FREE_FORKNODES list) has a Non-NULL NEXT pointer.\n");
    return(stream ? 0 : 14);
  }
  else if ((nextNodeState == 2) && (pq->next == NULL)) // FORKnode NEXT should be non-NULL for FORKnodes in TREE or garbage on the FREE_FORKRINGS list.
  {
    if (stream)
      fprintf(stream, "======> DEBUG ERROR (015):  A FORKnode in a FORKring in the TREE or on the FREE_FORKRINGS list has a NULL NEXT pointer.\n");
    return(stream ? 0 : 15);
  }

  return(0);
}


// End.
