/* Interface
   Version 4.0. (c) Copyright 1992-2013 by the University of Washington.
   Written by Sean T. Lamont and Michal Palczewski.
   For use with the Macintosh version of the Phylogeny Inference Package,
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed.

   interface.h:  access to interface.c, a 2 window text/graphics
   environment, with a scrolling text window and c-like I/O functions.
   This also sets up some defines for the standard c stuff.
*/


#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#ifdef OSX_CARBON
#define MAC_OFFSET 60
#endif


/* function prototypes */
void   macsetup(char *, char *);
void   queryevent(void);
void   eventloop(void);
void    process_window_closure(void);
int    handleevent(void);
void   textmode(void);
void    gfxmode(void);
//pascal void scroll(void);
void scroll(void);
int    process_char(void);
void paint_gfx_window(void);
#ifndef OSX_CARBON
void resize_gfx_window(EventRecord ev);
#endif
void menu_select(long what);
/*debug void fixmacfile(char *);*/


typedef struct {
  char* fn;
  double* xo;
  double* yo;
  double* scale;
  long nt;
  void* root;

} mpreviewparams;

extern mpreviewparams macpreviewparms;

/* function prototypes */


#ifdef __MWERKS__
#define MAC
#endif
#endif


// End.
