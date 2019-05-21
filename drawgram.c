/* Version 4.0.  Copyright (c) 1986-2013 by The University of Washington and
   Written by Joseph Felsenstein and Christopher A. Meacham.  Additional code
   written by Hisashi Horino, Sean Lamont, Andrew Keefe, Daniel Fineman,
   Akiko Fuseki, Doug Buxton and Michal Palczewski. Java interface added
   by Jim McGill.  Permission is granted to copy, distribute, and modify
   this program provided that (1) this copyright message is not removed
   and (2) no fee is charged for this program. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef OSX_CARBON
#include <Carbon/Carbon.h>
#endif

#include "phylip.h"
#include "draw.h"

#ifdef MAC
char* about_message =
  "Drawgram unrooted tree plotting program\r"
  "PHYLIP version 4.0 (c) Copyright 1986-2013\r"
  "by The University of Washington.\r"
  "Written by Joseph Felsenstein and Christopher A. Meacham.\r"
  "Additional code written by Hisashi Horino, Sean Lamont, Andrew Keefe,\r"
  "Daniel Fineman, Akiko Fuseki, Doug Buxton and Michal Palczewski.\r"
  "Permission is granted to copy, distribute and modify this program\r"
  "provided that\r"
  "(1) This copyright message is not removed and\r"
  "(2) no fee is charged for this program.";
#endif

//#define JAVADEBUG

tree * curtree;
long nonodes;


#define gap 0.5    /* distance in character heights between the end
                      of a branch and the start of the name */

char trefilename[FNMLNGTH];
//char *progname;

long nextnode, payge, numlines,  iteration;
boolean canbeplotted, haslengths, uselengths, rescaled, firstscreens,
  multiplot, finished;
double topoflabels, bottomoflabels, rightoflabels,
  leftoflabels, tipspacing,maxheight, scale,
  xoffset, yoffset, nodespace, stemlength, treedepth,
  rootx, rooty, bscale, xx0, yy0, fontheight, maxx, minx, maxy, miny;

double *textlength, *firstlet;
treestyle style;
fonttype font;
Char ch;
double trweight;   /* starting here, needed to make sccs version happy */
boolean goteof;
long *zeros;     /* ... down to here */

static enum {weighted, intermediate, centered, inner, vshaped} nodeposition;

#ifdef PHYLIP_CAN_DRAW_WITH_X
String res[]= {
  "*.input: True",
  "*.menubar.orientation: horizontal",
  "*.menubar.borderWidth: 0",
  "*.drawing_area.background: #CCFFFF",
  "*.drawing_area.foreground: #000000",
  "*.menubar.right: ChainLeft",
  "*.menubar.bottom: ChainTop",
  "*.menubar.top: ChainTop",
  "*.menubar.left: ChainLeft",
  "*.drawing_area.fromVert: menubar",
  "*.drawing_area.top: ChainTop",
  "*.drawing_area.bottom: ChainBottom",
  "*.drawing_area.left: ChainLeft",
  "*.drawing_area.right: ChainRight",
  "*.dialog.label: "
  "Drawgram unrooted tree plotting program\\n"
  "PHYLIP version 4.0 (c) Copyright 1986-2013\\n"
  "by The University of Washington.\\n"
  "Written by Joseph Felsenstein and Christopher A. Meacham.\\n"
  "Additional code written by Hisashi Horino, Sean Lamont, Andrew Keefe,\\n"
  "Daniel Fineman, Akiko Fuseki, Doug Buxton and Michal Palczewski.\\n"
  "Permission is granted to copy, distribute and modify this program\\n"
  "provided that\\n"
  "(1) This copyright message is not removed and\\n"
  "(2) no fee i charged for this program.",
  NULL
};
#endif


#ifndef OLDC
/* function prototypes */
void   initialparms(void);
char   showparms(void);
void   getparms(char);
void   calctraverse(node *, double, double *);
void   calculate(void);
void   rescale(void);
void   setup_environment(int argc, Char *argv[]);
void   user_loop(void);
void   drawgram(char* intreename, char* fontfilename, char* plotfilename, char* plotfileopt,
                char* treegrows, char* stylekind, int usebranchlengths, double labelangle,
                int scalebranchlength, double branchscale, double breadthdepthratio,
                double stemltreedratio, double chhttipspratio, double xmarginratio,
                double ymarginratio, char* ancnodes, int  dofinalplot, char* finalplotkind);
void   makebox(char *, double *, double *, double *, long);
/* function prototypes */
#endif


void initialparms(void)
{
  /* initialize parameters */
  plotter = DEFPLOTTER;
  previewer = DEFPREV;
  paperx = 20.6375;
  pagex = 20.6375;
  papery = 26.9875;
  pagey = 26.9875;
  strcpy(fontname,"Times-Roman");
  plotrparms(spp);   /* initial, possibly bogus, parameters */
  style = phenogram;
  grows = horizontal;
  labelrotation = 90.0;
  nodespace = 3.0;
  stemlength = 0.05;
  treedepth = 0.5 / 0.95;
  rescaled = true;
  bscale = 1.0;
  uselengths = haslengths;
  if (uselengths)
    nodeposition = weighted;
  else
    nodeposition = centered;
  xmargin = 0.08 * xsize;
  ymargin = 0.08 * ysize;
  preview = true;
  hpmargin = 0.02*pagex;
  vpmargin = 0.02*pagey;
}  /* initialparms */


char showparms(void)
{
  char input[32];
  Char ch;
  char cstr[32];

  if (!firstscreens)
    clearit();
  printf("\nRooted tree plotting program version %s\n\n", VERSION);
  printf("Here are the settings: \n");
  printf(" 0  Screen type (IBM PC, ANSI):  %s\n",
         ibmpc ? "IBM PC" : ansi  ? "ANSI"  : "(none)");
  printf(" P       Final plotting device: ");
  switch (plotter) {
    case lw:
      printf(" Postscript printer\n");
      break;
    case svg:
      printf(" SVG scalable vector graphics\n");
      break;
    case pcl:
      printf(" HP Laserjet compatible printer (%d DPI)\n", (int) hpresolution);
      break;
    case epson:
      printf(" Epson dot-matrix printer\n");
      break;
    case pcx:
      printf(" PCX file for PC Paintbrush drawing program (%s)\n",
             (resopts == 1) ? "EGA 640x350" :
             (resopts == 2) ? "VGA 800x600" : "VGA 1024x768");
      break;
    case pict:
      printf(" Macintosh PICT file for drawing program\n");
      break;
    case idraw:
      printf(" Idraw drawing program\n");
      break;
    case fig:
      printf(" Xfig drawing program\n");
      break;
    case hp:
      printf(" HPGL graphics language for HP plotters\n");
      break;
    case xbm:
      printf(" X Bitmap file format (%d by %d resolution)\n", (int)xsize,
             (int)ysize);
      break;
    case bmp:
      printf(" MS-Windows Bitmap (%d by %d resolution)\n", (int)xsize,
             (int)ysize);
      break;
    case gif:
      printf(" Compuserve GIF format (%d by %d)\n",(int)xsize,(int)ysize);
      break;
    case ibm:
      printf(" IBM PC graphics (CGA, EGA, or VGA)\n");
      break;
    case tek:
      printf(" Tektronix graphics screen\n");
      break;
    case decregis:
      printf(" DEC ReGIS graphics (VT240 or DECTerm)\n");
      break;
    case houston:
      printf(" Houston Instruments plotter\n");
      break;
    case toshiba:
      printf(" Toshiba 24-pin dot matrix printer\n");
      break;
    case citoh:
      printf(" Imagewriter or C.Itoh/TEC/NEC 9-pin dot matrix printer\n");
      break;
    case oki:
      printf(" old Okidata 9-pin dot matrix printer\n");
      break;
    case ray:
      printf(" Rayshade ray-tracing program file format\n");
      break;
    case pov:
      printf(" POV ray-tracing program file format\n");
      break;
    case vrml:
      printf(" VRML, Virtual Reality Markup Language\n");
      break;
    case mac:
    case other:
      printf(" (Current output device unannounced)\n");
      break;
    default:        /*case xpreview not handled */
      break;
  }
#ifdef MACOS10
  preview = false;
  printf(" (Preview not available on Macs)\n");
#else
  printf(" V           Previewing device: ");
  if (!preview)
    printf(" (none)\n");
  else {
    switch (previewer) {
      case ibm:
        printf(" IBM PC graphics\n");
        break;
      case xpreview:
        printf(" X Windows display\n");
        break;
      case tek:
        printf(" Tektronix graphics screen\n");
        break;
      case mac:
        printf(" Macintosh graphics screen\n");
        break;
      case decregis:
        printf(" DEC ReGIS graphics (VT240 or DECTerm)\n");
        break;
      case winpreview:
        printf(" MS Windows display\n");
        break;
      case other:
        printf(" (Current previewing device unannounced)\n");
        break;
      default:   /* all other cases */
        break;
    }
  }
#endif
  printf(" H                  Tree grows:  ");
  printf((grows == vertical) ? "Vertically\n" : "Horizontally\n");
  printf(" S                  Tree style:  %s\n",
         (style == cladogram) ? "Cladogram" :
         (style == phenogram)  ? "Phenogram" :
         (style == curvogram) ? "Curvogram" :
         (style == eurogram)  ? "Eurogram"  :
         (style == swoopogram) ? "Swoopogram" : "Circular");
  printf(" B          Use branch lengths:  ");
  if (haslengths) {
    if (uselengths)
      printf("Yes\n");
    else
      printf("No\n");
  } else
    printf("(no branch lengths available)\n");
  if (style != circular) {
    printf(" L             Angle of labels:");
    if (labelrotation < 10.0)
      printf("%5.1f\n", labelrotation);
    else
      printf("%6.1f\n", labelrotation);
  }
  printf(" R      Scale of branch length:");
  if (rescaled)
    printf("  Automatically rescaled\n");
  else
    printf("  Fixed:%6.2f cm per unit branch length\n", bscale);
  printf(" D       Depth/Breadth of tree:%6.2f\n", treedepth);
  printf(" T      Stem-length/tree-depth:%6.2f\n", stemlength);
  printf(" C    Character ht / tip space:%8.4f\n", 1.0 / nodespace);
  printf(" A             Ancestral nodes:  %s\n",
         (nodeposition == weighted)     ? "Weighted"     :
         (nodeposition == intermediate) ? "Intermediate" :
         (nodeposition == centered)     ? "Centered"     :
         (nodeposition == inner)        ? "Inner"        :
         "So tree is V-shaped");
  if (plotter == lw || plotter == idraw || plotter == svg ||
      (plotter == fig ) ||
      (plotter == pict && ((grows == vertical && labelrotation == 0.0) ||
                           (grows == horizontal && labelrotation == 90.0))))
    printf(" F                        Font:  %s\n",fontname);
  if ((plotter == pict && ((grows == vertical && labelrotation == 0.0) ||
                           (grows == horizontal && labelrotation == 90.0)))
      && (strcmp(fontname,"Hershey") != 0))
    printf(" Q        Pict Font Attributes:  %s, %s, %s, %s\n",
           (pictbold   ? "Bold"   : "Medium"),
           (pictitalic ? "Italic" : "Regular"),
           (pictshadow ? "Shadowed": "Unshadowed"),
           (pictoutline ? "Outlined" : "Unoutlined"));
  if (plotter == ray) {
    printf(" M          Horizontal margins:%6.2f pixels\n", xmargin);
    printf(" M            Vertical margins:%6.2f pixels\n", ymargin);
  } else {
    printf(" M          Horizontal margins:%6.2f cm\n", xmargin);
    printf(" M            Vertical margins:%6.2f cm\n", ymargin);
  }
  printf(" #              Pages per tree:  ");
  /* Add 0.5 to clear up truncation problems. */
  if (((int) ((pagex / paperx) + 0.5) == 1) &&
      ((int) ((pagey / papery) + 0.5) == 1))
    /* If we're only using one page per tree, */
    printf ("one page per tree\n") ;
  else
    printf ("%.0f by %.0f pages per tree\n",
            (pagey-vpmargin) / (papery-vpmargin),
            (pagex-hpmargin) / (paperx-hpmargin)) ;

  for (;;) {
    printf("\n Y to accept these or type the letter for one to change\n");
    phyFillScreenColor();
    getstryng(input);
    uppercase(&input[0]);
    ch=input[0];
    if (plotter == idraw || plotter == lw || plotter == svg )
      strcpy(cstr,"#Y0PVHSBLMRDTCAF");
    else if (plotter == pict) {
      if (((grows == vertical && labelrotation == 0.0) ||
           (grows == horizontal && labelrotation == 90.0)))
        strcpy(cstr,"#Y0PVHSBLMRDTCAFQ");
      else
        strcpy(cstr,"#Y0PVHSBLMRDTCA"); }
    else
      strcpy(cstr,"#Y0PVHSBLMRDTCA");

    if  (strchr(cstr,ch))
      break;
    printf(" That letter is not one of the menu choices.  Type\n");
  }
  return ch;
}  /* showparms */


void getparms(char numtochange)
{
  /* get from user the relevant parameters for the plotter and diagram */
  long loopcount;
  Char ch;
  char input[100];
  boolean ok;
  int m, n;

  n = (int)((pagex-hpmargin-0.01)/(paperx-hpmargin)+1.0);
  m = (int)((pagey-vpmargin-0.01)/(papery-vpmargin)+1.0);
  switch (numtochange) {

    case '0':
      initterminal(&ibmpc, &ansi);
      break;

    case 'P':
      getplotter();
      break;

    case 'V':
      getpreview();
      break;

    case 'H':
      if (grows == vertical)
        grows = horizontal;
      else
        grows = vertical;
      break;

    case 'S':
      clearit() ;
      printf("What style tree is this to be (currently set to %s):\n",
             (style == cladogram) ? "Cladogram" :
             (style == phenogram)  ? "Phenogram" :
             (style == curvogram) ? "Curvogram" :
             (style == eurogram)  ? "Eurogram"  :
             (style == swoopogram) ? "Swoopogram" : "Circular") ;
      printf(" C    Cladogram -- v-shaped \n") ;
      printf(" P    Phenogram -- branches are square\n") ;
      printf(" V    Curvogram -- branches are 1/4 of an ellipse\n") ;
      printf(" E    Eurogram -- branches angle outward, then up\n");
      printf(" S    Swoopogram -- branches curve outward then reverse\n") ;
      printf(" O    Circular tree\n");
      do {
        printf("\n Type letter of style to change to (C, P, V, E, S or O):\n");
        phyFillScreenColor();
        if(scanf("%c%*[^\n]", &ch)) {}  // Read char and scan to EOL.
        (void)getchar();
        uppercase(&ch);
      } while (ch != 'C' && ch != 'P' && ch != 'V'
               && ch != 'E' && ch != 'S' && ch != 'O');
      switch (ch) {

        case 'C':
          style = cladogram;
          break;

        case 'P':
          style = phenogram;
          break;

        case 'E':
          style = eurogram;
          break;

        case 'S':
          style = swoopogram;
          break;

        case 'V':
          style = curvogram;
          break;

        case 'O':
          style = circular;
          treedepth = 1.0;
          break;
      }
      break;

    case 'B':
      if (haslengths) {
        uselengths = !uselengths;
        if (!uselengths)
          nodeposition = weighted;
        else
          nodeposition = intermediate;
      } else {
        printf("Cannot use lengths since not all of them exist\n");
        uselengths = false;
      }
      break;

    case 'L':
      clearit();
      printf("\n(Considering the tree as if it \"grew\" vertically:)\n");
      printf("Are the labels to be plotted vertically (90),\n");
      printf(" horizontally (0), or at a 45-degree angle?\n");
      loopcount = 0;
      do {
        printf(" Choose an angle in degrees from 90 to 0:\n");
        phyFillScreenColor();
        if(scanf("%lf%*[^\n]", &labelrotation)) {} // Read number and scan to EOL.
        (void)getchar();
        uppercase(&ch);
        countup(&loopcount, 10);
      } while (labelrotation < 0.0 && labelrotation > 90.0);
      break;

    case 'M':
      clearit();
      printf("\nThe tree will be drawn to fit in a rectangle which has \n");
      printf(" margins in the horizontal and vertical directions of:\n");
      if (plotter == ray) {
        printf(
          "%6.2f pixels (horizontal margin) and%6.2f pixels (vertical margin)\n",
          xmargin, ymargin);
      }
      else {
        printf("%6.2f cm (horizontal margin) and%6.2f cm (vertical margin)\n",
               xmargin, ymargin);
      }
      putchar('\n');
      loopcount = 0;
      do {
        if (plotter == ray)
          printf(" New value (in pixels) of horizontal margin?\n");
        else
          printf(" New value (in cm) of horizontal margin?\n");
        phyFillScreenColor();
        if(scanf("%lf%*[^\n]", &xmargin)) {} // Read number and scan to EOL.
        (void)getchar();
        ok = ((unsigned)xmargin < xsize / 2.0);
        if (!ok)
          printf(" Impossible value.  Please retype it.\n");
        countup(&loopcount, 10);
      } while (!ok);
      loopcount = 0;
      do {
        if (plotter == ray)
          printf(" New value (in pixels) of vertical margin?\n");
        else
          printf(" New value (in cm) of vertical margin?\n");
        phyFillScreenColor();
        if(scanf("%lf%*[^\n]", &ymargin)) {} // Read number and scan to EOL.
        (void)getchar();
        ok = ((unsigned)ymargin < ysize / 2.0);
        if (!ok)
          printf(" Impossible value.  Please retype it.\n");
        countup(&loopcount, 10);
      } while (!ok);

      break;

    case 'R':
      rescaled = !rescaled;
      if (!rescaled) {
        printf("Centimeters per unit branch length?\n");
        phyFillScreenColor();
        if(scanf("%lf%*[^\n]", &bscale)) {} // Read number and scan to EOL.
        (void)getchar();
      }
      break;

    case 'D':
      printf("New value of depth of tree as fraction of its breadth?\n");
      phyFillScreenColor();
      if(scanf("%lf%*[^\n]", &treedepth)) {} // Read number and scan to EOL.
      (void)getchar();
      break;

    case 'T':
      loopcount = 0;
      do {
        printf("New value of stem length as fraction of tree depth?\n");
        phyFillScreenColor();
        if(scanf("%lf%*[^\n]", &stemlength)) {} // Read number and scan to EOL.
        (void)getchar();
        countup(&loopcount, 10);
      } while ((unsigned)stemlength >= 0.9);
      break;

    case 'C':
      printf("New value of character height as fraction of tip spacing?\n");
      phyFillScreenColor();
      if(scanf("%lf%*[^\n]", &nodespace)) {} // Read number and scan to EOL.
      (void)getchar();
      nodespace = 1.0 / nodespace;
      break;

    case '#':
      loopcount = 0;
      for (;;) {
        clearit();
        printf("  Page Specifications Submenu\n\n");
        printf(" L   Output size in pages: %.0f down by %.0f across\n",
               (pagey / papery), (pagex / paperx));
        printf(" P   Physical paper size: %1.5f by %1.5f cm\n",paperx,papery);
        printf(" O   Overlap Region: %1.5f %1.5f cm\n",hpmargin,vpmargin);
        printf(" M   main menu\n");
        getstryng(input);
        ch = input[0];
        uppercase(&ch);
        switch (ch) {
          case 'L':
            printf("Number of pages in height:\n");
            phyFillScreenColor();
            getstryng(input);
            m = atoi(input);
            printf("Number of pages in width:\n");
            getstryng(input);
            n = atoi(input);
            break;
          case 'P':
            printf("Paper Width (in cm):\n");
            phyFillScreenColor();
            getstryng(input);
            paperx = atof(input);
            printf("Paper Height (in cm):\n");
            getstryng(input);
            papery = atof(input);
            break;
          case 'O':
            printf("Horizontal Overlap (in cm):");
            phyFillScreenColor();
            getstryng(input);
            hpmargin = atof(input);
            printf("Vertical Overlap (in cm):");
            phyFillScreenColor();
            getstryng(input);
            vpmargin = atof(input);
          case 'M':
            break;
          default:
            printf("Please enter L, P, O , or M.\n");
            break;
        }
        pagex = ((double)n * (paperx-hpmargin)+hpmargin);
        pagey = ((double)m * (papery-vpmargin)+vpmargin);
        if (ch == 'M')
          break;
        countup(&loopcount, 10);
      }
      break;

    case 'A':
      clearit();
      printf("Should interior node positions:\n");
      printf(" be Intermediate between their immediate descendants,\n");
      printf("    Weighted average of tip positions\n");
      printf("    Centered among their ultimate descendants\n");
      printf("    iNnermost of immediate descendants\n");
      printf(" or so that tree is V-shaped\n");
      loopcount = 0;
      do {
        printf(" (type I, W, C, N or V):\n");
        phyFillScreenColor();
        if(scanf("%c%*[^\n]", &ch)) {}  // Read char and scan to EOL.
        (void)getchar();
        uppercase(&ch);
        countup(&loopcount, 10);
      } while (ch != 'I' && ch != 'W' && ch != 'C' && ch != 'N' && ch != 'V');
      switch (ch) {

        case 'W':
          nodeposition = weighted;
          break;

        case 'I':
          nodeposition = intermediate;
          break;

        case 'C':
          nodeposition = centered;
          break;

        case 'N':
          nodeposition = inner;
          break;

        case 'V':
          nodeposition = vshaped;
          break;
      }
      break;

    case 'F':
      printf("Enter font name or \"Hershey\" for the default font\n");
      phyFillScreenColor();
      getstryng(fontname);
      break;

    case 'Q':
      clearit();
      loopcount = 0;
      do {
        printf("Italic? (Y/N)\n");
        phyFillScreenColor();
        getstryng(input);
        input[0] = toupper(input[0]);
        countup(&loopcount, 10);
      } while (input[0] != 'Y' && input[0] != 'N');
      pictitalic = (input[0] == 'Y');
      loopcount = 0;
      do {
        printf("Bold? (Y/N)\n");
        phyFillScreenColor();
        getstryng(input);
        input[0] = toupper(input[0]);
        countup(&loopcount, 10);
      } while (input[0] != 'Y' && input[0] != 'N');
      pictbold = (input[0] == 'Y');
      loopcount = 0;
      do {
        printf("Shadow? (Y/N)\n");
        phyFillScreenColor();
        getstryng(input);
        input[0] = toupper(input[0]);
        countup(&loopcount, 10);
      } while (input[0] != 'Y' && input[0] != 'N');
      pictshadow = (input[0] == 'Y');
      loopcount = 0;
      do {
        printf("Outline? (Y/N)\n");
        phyFillScreenColor();
        getstryng(input);
        input[0] = toupper(input[0]);
        countup(&loopcount, 10);
      } while (input[0] != 'Y' && input[0] != 'N');
      pictoutline = (input[0] == 'Y');break;
  }
}  /* getparms */


void calctraverse(node *p, double lengthsum, double *tipx)
{
  /* traverse to establish initial node coordinates */
  double x1, y1, x2, y2, x3, x4, x5, w1, w2, sumwx, sumw, nodeheight;
  node *pp, *plast, *panc;

  //printf ("uselengths: %i\n", uselengths);
  if (p == curtree->root)
    nodeheight = 0.0;
  else if (uselengths)
    nodeheight = lengthsum + fabs(p->oldlen);
  else
    nodeheight = 1.0;
  if (nodeheight > maxheight)
    maxheight = nodeheight;
  if (p->tip) {
    p->xcoord = *tipx;
    ((draw_node*)p)->tipsabove = 1;
    if (uselengths)
      p->ycoord = nodeheight;
    else
      p->ycoord = 1.0;
    *tipx += tipspacing;
    return;
  } else {
    sumwx = 0.0;
    sumw = 0.0;
    ((draw_node*)p)->tipsabove = 0;
    pp = p->next;
    x3 = 0.0;
    do {
      calctraverse(pp->back, nodeheight, tipx);
      ((draw_node*)p)->tipsabove += ((draw_node*)pp->back)->tipsabove;
      sumw += ((draw_node*)pp->back)->tipsabove;
      sumwx += ((draw_node*)pp->back)->tipsabove * pp->back->xcoord;
      if (fabs(pp->back->xcoord - 0.5) < fabs(x3 - 0.5))
        x3 = pp->back->xcoord;
      plast = pp;
      pp = pp->next;
    } while (pp != p);
    x1 = p->next->back->xcoord;
    x2 = plast->back->xcoord;
    y1 = p->next->back->ycoord;
    y2 = plast->back->ycoord;

    switch (nodeposition) {

      case weighted:
        w1 = y1 - p->ycoord;
        w2 = y2 - p->ycoord;
        if (w1 + w2 <= 0.0)
          p->xcoord = (x1 + x2) / 2.0;
        else
          p->xcoord = (w2 * x1 + w1 * x2) / (w1 + w2);
        break;

      case intermediate:
        p->xcoord = (x1 + x2) / 2.0;
        break;

      case centered:
        p->xcoord = sumwx / sumw;
        break;

      case inner:
        p->xcoord = x3;
        break;

      case vshaped:
        if (iteration > 1) {
          if (!(p == curtree->root)) {
            panc = curtree->nodep[p->back->index-1];
            w1 = p->ycoord - panc->ycoord;
            w2 = y1 - p->ycoord;
            if (w1+w2 < 0.000001)
              x4 = (x1+panc->xcoord)/2.0;
            else
              x4 = (w1*x1+w2*panc->xcoord)/(w1+w2);
            w2 = y2 - p->ycoord;
            if (w1+w2 < 0.000001)
              x5 = (x2+panc->xcoord)/2.0;
            else
              x5 = (w1*x2+w2*panc->xcoord)/(w1+w2);
            if (panc->xcoord < p->xcoord)
              p->xcoord = x5;
            else
              p->xcoord = x4;
          }
          else {
            if ((y1-2*p->ycoord+y2) < 0.000001)
              p->xcoord = (x1+x2)/2;
            else
              p->xcoord = ((y2-p->ycoord)*x1+(y1-p->ycoord))/(y1-2*p->ycoord+y2);
          }
        }
        break;
    }
    if (uselengths) {
      p->ycoord = nodeheight;
      return;
    }
    if (nodeposition != inner) {
      p->ycoord = (y1 + y2 - sqrt((y1 + y2) * (y1 + y2) - 4 * (y1 * y2 -
                                                               (x2 - p->xcoord) * (p->xcoord - x1)))) / 2.0;
      /* this formula comes from the requirement that the vector from
         (x,y) to (x1,y1) be at right angles to that from (x,y) to (x2,y2) */
      return;
    }
    if (fabs(x1 - 0.5) > fabs(x2 - 0.5)) {
      p->ycoord = y1 + x1 - x2;
      w1 = y2 - p->ycoord;
    } else {
      p->ycoord = y2 + x1 - x2;
      w1 = y1 - p->ycoord;
    }
    if (w1 < epsilon)
      p->ycoord -= fabs(x1 - x2);
  }
}  /* calctraverse */


void calculate(void)
{
  /* compute coordinates for tree */
  double tipx;
  double sum, temp, maxtextlength, maxfirst=0, leftfirst, angle;
  double lef = 0.0, rig = 0.0, top = 0.0, bot = 0.0;
  double *firstlet, *textlength;
  long i;

  firstlet = (double *)Malloc(nextnode*sizeof(double));
  textlength = (double *)Malloc(nextnode*sizeof(double));
  for (i = 0; i < nextnode; i++) {
    curtree->nodep[i]->xcoord = 0.0;
    curtree->nodep[i]->ycoord = 0.0;
    if (curtree->nodep[i]->naymlength > 0)
      firstlet[i] = lengthtext(curtree->nodep[i]->nayme, 1L,fontname,font);
    else
      firstlet[i] = 0.0;
  }
  i = 0;
  do
    i++;
  while (!curtree->nodep[i]->tip);
  leftfirst = firstlet[i];
  maxheight = 0.0;
  maxtextlength = 0.0;
  for (i = 0; i < nextnode; i++) {
    if (curtree->nodep[i]->tip) {
      textlength[i] = lengthtext(curtree->nodep[i]->nayme,
                                 curtree->nodep[i]->naymlength, fontname, font);
      if (textlength[i]-0.5*firstlet[i] > maxtextlength) {
        maxtextlength = textlength[i]-0.5*firstlet[i];
        maxfirst = firstlet[i];
      }
    }
  }
  maxtextlength = maxtextlength + 0.5*maxfirst;
  fontheight = heighttext(font,fontname);
  if (style == circular) {
    if (grows == vertical)
      angle = pi / 2.0;
    else
      angle = 2.0*pi;
  } else
    angle = pi * labelrotation / 180.0;
  maxtextlength /= fontheight;
  maxfirst /= fontheight;
  leftfirst /= fontheight;
  for (i = 0; i < nextnode; i++) {
    if (curtree->nodep[i]->tip) {
      textlength[i] /= fontheight;
      firstlet[i] /= fontheight;
    }
  }
  if (spp > 1)
    labelheight = 1.0 / (nodespace * (spp - 1));
  else
    labelheight = 1.0 / nodespace;
  if (angle < pi / 6.0)
    tipspacing = (nodespace
                  + cos(angle) * (maxtextlength - 0.5*maxfirst)) * labelheight;
  else if (spp > 1) {
    if (style == circular) {
      tipspacing = 1.0 / spp;
    } else
      tipspacing = 1.0 / (spp - 1.0);
  } else
    tipspacing = 1.0;
  finished = false;
  iteration = 1;
  do {
    if (style == circular)
      tipx = 1.0/(2.0*(double)spp);
    else
      tipx = 0.0;
    sum = 0.0;
    calctraverse(curtree->root, sum, &tipx);
    iteration++;
  }
  while ((nodeposition == vshaped) && (iteration < 4*spp));
  rooty = curtree->root->ycoord;
  labelheight *= 1.0 - stemlength;
  //printf("stemlength: %f treedepth: %f rooty: %f maxheight: %f rescaled: %i\n", stemlength, treedepth, rooty, maxheight, rescaled);
  for (i = 0; i < nextnode; i++) {
    //printf("before rescale curtree->nodep[i]->ycoord: %f\n", curtree->nodep[i]->ycoord);
    if (rescaled) {
      if (style != circular)
        curtree->nodep[i]->xcoord *= 1.0 - stemlength;
      curtree->nodep[i]->ycoord = stemlength * treedepth + (1.0 - stemlength) *
        treedepth * (curtree->nodep[i]->ycoord - rooty) / (maxheight - rooty);
      ((draw_node*)curtree->nodep[i])->oldtheta = angle;
    } else {
      curtree->nodep[i]->xcoord = curtree->nodep[i]->xcoord * (maxheight - rooty) / treedepth;
      curtree->nodep[i]->ycoord = stemlength / (1 - stemlength) * (maxheight - rooty) +
        curtree->nodep[i]->ycoord;
      ((draw_node*)curtree->nodep[i])->oldtheta = angle;
    }
    //printf("after rescale curtree->nodep[i]->ycoord: %f\n", curtree->nodep[i]->ycoord);
  }
  topoflabels = 0.0;
  bottomoflabels = 0.0;
  leftoflabels = 0.0;
  rightoflabels = 0.0;
  if  (style == circular) {
    for (i = 0; i < nextnode; i++) {
      temp = curtree->nodep[i]->xcoord;
      if (grows == vertical) {
        curtree->nodep[i]->xcoord = (1.0+curtree->nodep[i]->ycoord
                                     * cos((1.5-2.0*temp)*pi)/treedepth)/2.0;
        curtree->nodep[i]->ycoord = (1.0+curtree->nodep[i]->ycoord
                                     * sin((1.5-2.0*temp)*pi)/treedepth)/2.0;
        ((draw_node*)curtree->nodep[i])->oldtheta = (1.5-2.0*temp)*pi;
      } else {
        curtree->nodep[i]->xcoord = (1.0+curtree->nodep[i]->ycoord
                                     * cos((1.0-2.0*temp)*pi)/treedepth)/2.0;
        curtree->nodep[i]->ycoord = (1.0+curtree->nodep[i]->ycoord
                                     * sin((1.0-2.0*temp)*pi)/treedepth)/2.0;
        ((draw_node*)curtree->nodep[i])->oldtheta = (1.0-2.0*temp)*pi;
      }
    }
    tipspacing *= 2.0*pi;
  }
  maxx = curtree->nodep[0]->xcoord;
  maxy = curtree->nodep[0]->ycoord;
  minx = curtree->nodep[0]->xcoord;
  if (style == circular)
    miny = curtree->nodep[0]->ycoord;
  else
    miny = 0.0;
  //printf("Before loop nmaxy: %f miny: %f maxx: %f minx: %f\n", maxy, miny, maxx, minx);
  for (i = 1; i < nextnode; i++) {
    if (curtree->nodep[i]->xcoord > maxx)
      maxx = curtree->nodep[i]->xcoord;
    if (curtree->nodep[i]->ycoord > maxy)
      maxy = curtree->nodep[i]->ycoord;
    if (curtree->nodep[i]->xcoord < minx)
      minx = curtree->nodep[i]->xcoord;
    if (curtree->nodep[i]->ycoord < miny)
      miny = curtree->nodep[i]->ycoord;
    //printf("end of loop curtree->nodep[i]: (%f, %f) nmaxy: %f miny: %f maxx: %f minx: %f\n", curtree->nodep[i]->xcoord, curtree->nodep[i]->ycoord, maxy, miny, maxx, minx);
  }
  //printf("*****in calculate*****\nmaxy: %f miny: %f maxx: %f minx: %f\n", maxy, miny, maxx, minx);
  if  (style == circular) {
    for (i = 0; i < nextnode; i++) {
      if (curtree->nodep[i]->tip) {
        angle = ((draw_node*)curtree->nodep[i])->oldtheta;
        if (cos(angle) < 0.0)
          angle -= pi;
        top = (curtree->nodep[i]->ycoord - maxy) / labelheight +
          sin(((draw_node*)curtree->nodep[i])->oldtheta);
        rig = (curtree->nodep[i]->xcoord - maxx) / labelheight +
          cos(((draw_node*)curtree->nodep[i])->oldtheta);
        bot = (miny - curtree->nodep[i]->ycoord) / labelheight -
          sin(((draw_node*)curtree->nodep[i])->oldtheta);
        lef = (minx - curtree->nodep[i]->xcoord) / labelheight -
          cos(((draw_node*)curtree->nodep[i])->oldtheta);
        if (cos(((draw_node*)curtree->nodep[i])->oldtheta) > 0) {
          if (sin(angle) > 0.0)
            top += sin(angle) * textlength[i];
          top += sin(angle - 1.25 * pi) * gap * firstlet[i];
          if (sin(angle) < 0.0)
            bot -= sin(angle) * textlength[i];
          bot -= sin(angle - 0.75 * pi) * gap * firstlet[i];
          if (sin(angle) > 0.0)
            rig += cos(angle - 0.75 * pi) * gap * firstlet[i];
          else
            rig += cos(angle - 1.25 * pi) * gap * firstlet[i];
          rig += cos(angle) * textlength[i];
          if (sin(angle) > 0.0)
            lef -= cos(angle - 1.25 * pi) * gap * firstlet[i];
          else
            lef -= cos(angle - 0.75 * pi) * gap * firstlet[i];
        } else {
          if (sin(angle) < 0.0)
            top -= sin(angle) * textlength[i];
          top += sin(angle + 0.25 * pi) * gap * firstlet[i];
          if (sin(angle) > 0.0)
            bot += sin(angle) * textlength[i];
          bot -= sin(angle - 0.25 * pi) * gap * firstlet[i];
          if (sin(angle) > 0.0)
            rig += cos(angle - 0.25 * pi) * gap * firstlet[i];
          else
            rig += cos(angle + 0.25 * pi) * gap * firstlet[i];
          if (sin(angle) < 0.0)
            rig += cos(angle) * textlength[i];
          if (sin(angle) > 0.0)
            lef -= cos(angle + 0.25 * pi) * gap * firstlet[i];
          else
            lef -= cos(angle - 0.25 * pi) * gap * firstlet[i];
          lef += cos(angle) * textlength[i];
        }
      }
      if (top > topoflabels)
        topoflabels = top;
      if (bot > bottomoflabels)
        bottomoflabels = bot;
      if (rig > rightoflabels)
        rightoflabels = rig;
      if (lef > leftoflabels)
        leftoflabels = lef;
    }
    topoflabels *= labelheight;
    bottomoflabels *= labelheight;
    leftoflabels *= labelheight;
    rightoflabels *= labelheight;
  }
  if (style != circular) {
    topoflabels = labelheight *
      (1.0 + sin(angle) * (maxtextlength - 0.5 * maxfirst)
       + cos(angle) * 0.5 * maxfirst);
    rightoflabels = labelheight *
      (cos(angle) * (textlength[nextnode-1] - 0.5 * maxfirst)
       + sin(angle) * 0.5);
    leftoflabels = labelheight * (cos(angle) * leftfirst * 0.5
                                  + sin(angle) * 0.5);
  }
  rooty = miny;
  free(firstlet);
  free(textlength);
}  /* calculate */


void rescale(void)
{
  /* compute coordinates of tree for plot or preview device */
  long i;
  double treeheight, treewidth, extrax, extray, temp;
  //printf("*****in rescale*****\n");
  //printf("maxy: %f miny: %f maxx: %f minx: %f\n", maxy, miny, maxx, minx);
  treeheight = maxy - miny;
  treewidth = maxx - minx;
  if (style == circular) {
    treewidth = 1.0;
    treeheight = 1.0;
    if (!rescaled) {
      if (uselengths) {
        labelheight *= (maxheight - rooty) / treedepth;
        topoflabels *= (maxheight - rooty) / treedepth;
        bottomoflabels *= (maxheight - rooty) / treedepth;
        leftoflabels *= (maxheight - rooty) / treedepth;
        rightoflabels *= (maxheight - rooty) / treedepth;
        treewidth *= (maxheight - rooty) / treedepth;
      }
    }
  }
  treewidth += rightoflabels + leftoflabels;
  treeheight += topoflabels + bottomoflabels;
  //printf("rightoflabels: %f leftoflabels: %f topoflabels: %f bottomoflabels: %f treewidth: %f treeheight: %f\n", rightoflabels, leftoflabels, topoflabels, bottomoflabels, treewidth, treeheight);
  if (grows == vertical) {
    if (!rescaled)
      expand = bscale;
    else {
      expand = (xsize - 2 * xmargin) / treewidth;
      if ((ysize - 2 * ymargin) / treeheight < expand)
        expand = (ysize - 2 * ymargin) / treeheight;
    }
    extrax = (xsize - 2 * xmargin - treewidth * expand) / 2.0;
    extray = (ysize - 2 * ymargin - treeheight * expand) / 2.0;
  } else {
    if (!rescaled)
      expand = bscale;
    else {
      expand = (ysize - 2 * ymargin) / treewidth;
      if ((xsize - 2 * xmargin) / treeheight < expand)
        expand = (xsize - 2 * xmargin) / treeheight;
    }
    extrax = (xsize - 2 * xmargin - treeheight * expand) / 2.0;
    extray = (ysize - 2 * ymargin - treewidth * expand) / 2.0;
  }
  //printf("in rescale rescaled: %i treewidth: %f treeheight: %f expand: %f\n", rescaled, treewidth, treeheight, expand);
  //printf("xsize: %f xmargin: %f ysize: %f ymargin: %f\n", xsize, xmargin, ysize, ymargin);
  for (i = 0; i < nextnode; i++) {
      //printf("before rescale i: %li, curtree->nodep[i]: %f, %f\n", i, curtree->nodep[i]->xcoord, curtree->nodep[i]->ycoord);
    curtree->nodep[i]->xcoord = expand * (curtree->nodep[i]->xcoord + leftoflabels);
    curtree->nodep[i]->ycoord = expand * (curtree->nodep[i]->ycoord + bottomoflabels);
    if ((style != circular) && (grows == horizontal)) {
      temp = curtree->nodep[i]->ycoord;
      curtree->nodep[i]->ycoord = expand * treewidth - curtree->nodep[i]->xcoord;
      curtree->nodep[i]->xcoord = temp;
      //printf("after rescale i: %li, curtree->nodep[i]: %f, %f\n", i, curtree->nodep[i]->xcoord, curtree->nodep[i]->ycoord);
  }
    curtree->nodep[i]->xcoord += xmargin + extrax;
    curtree->nodep[i]->ycoord += ymargin + extray;
  }
  if (style == circular) {
    xx0 = xmargin+extrax+expand*(leftoflabels + 0.5);
    yy0 = ymargin+extray+expand*(bottomoflabels + 0.5);
  }
  else if (grows == vertical)
    rooty = ymargin + extray;
  else
    rootx = xmargin + extrax;
}  /* rescale */


void plottree(node *p, node *q)
{
  /* plot part or all of tree on the plotting device */
  long i;
  double x00=0, y00=0, x1, y1, x2, y2, x3, y3, x4, y4,
    cc, ss, f, g, fract=0, minny, miny;
  node *pp;
    //printf("p: %p q: %p xscale: %f xoffset: %f yscale: %f yoffset %f\n", p, q, xscale, xoffset, yscale, yoffset);
    //fflush(stdout); 

  x2 = xscale * (xoffset + p->xcoord);
  y2 = yscale * (yoffset + p->ycoord);
  if (style == circular) {
    x00 = xscale * (xx0 + xoffset);
    y00 = yscale * (yy0 + yoffset);
  }
  if (p != curtree->root) {
    x1 = xscale * (xoffset + q->xcoord);
    y1 = yscale * (yoffset + q->ycoord);
    //printf("branch p: %f, %f q: %f, %f\n",p->xcoord, p->ycoord, q->xcoord, q->ycoord);
    //printf("line from: %f, %f to: %f, %f\n",x1, y1, x2, y2);
    //fflush(stdout); 

    switch (style) {

      case cladogram:
        plot(penup, x1, y1);
        plot(pendown, x2, y2);
        break;

      case phenogram:
        plot(penup, x1, y1);
        if (grows == vertical)
          plot(pendown, x2, y1);
        else
          plot(pendown, x1, y2);
        plot(pendown, x2, y2);
        break;

      case curvogram:
        plot(penup, x1, y1) ;
        curvespline(x1,y1,x2,y2,(boolean)(grows == vertical),20);
        break;

      case eurogram:
        plot(penup, x1, y1);
        if (grows == vertical)
          plot(pendown, x2, (2 * y1 + y2) / 3);
        else
          plot(pendown, (2 * x1 + x2) / 3, y2);
        plot(pendown, x2, y2);
        break;

      case swoopogram:
        plot(penup, x1, y1);
        if ((grows == vertical && fabs(y1 - y2) >= epsilon) ||
            (grows == horizontal && fabs(x1 - x2) >= epsilon)) {
          if (grows == vertical)
            miny = p->ycoord;
          else
            miny = p->xcoord;
          pp = q->next;
          while (pp != q) {
            if (grows == vertical)
              minny = pp->back->ycoord;
            else
              minny = pp->back->xcoord;
            if (minny < miny)
              miny = minny;
            pp = pp->next;
          }
          if (grows == vertical)
            miny = yscale * (yoffset + miny);
          else
            miny = xscale * (xoffset + miny);
          if (grows == vertical)
            fract = 0.3333 * (miny - y1) / (y2 - y1);
          else
            fract = 0.3333 * (miny - x1) / (x2 - x1);
        } if ((grows == vertical && fabs(y1 - y2) >= epsilon) ||
              (grows == horizontal && fabs(x1 - x2) >= epsilon)) {
          if (grows == vertical)
            miny = p->ycoord;
          else
            miny = p->xcoord;
          pp = q->next;
          while (pp != q) {
            if (grows == vertical)
              minny = pp->back->ycoord;
            else
              minny = pp->back->xcoord;
            if (minny < miny)
              miny = minny;
            pp = pp->next;
          }
          if (grows == vertical)
            miny = yscale * (yoffset + miny);
          else
            miny = xscale * (xoffset + miny);
          if (grows == vertical)
            fract = 0.3333 * (miny - y1) / (y2 - y1);
          else
            fract = 0.3333 * (miny - x1) / (x2 - x1);
        }
        swoopspline(x1,y1,x1+fract*(x2-x1),y1+fract*(y2-y1),
                    x2,y2,(boolean)(grows != vertical),segments);
        break;

      case circular:
        plot(penup, x1, y1);
        if (fabs(x1-x00)+fabs(y1-y00) > 0.00001) {
          g = ((x1-x00)*(x2-x00)+(y1-y00)*(y2-y00))
            /sqrt(((x1-x00)*(x1-x00)+(y1-y00)*(y1-y00))
                  *((x2-x00)*(x2-x00)+(y2-y00)*(y2-y00)));
          if (g > 1.0)
            g = 1.0;
          if (g < -1.0)
            g = -1.0;
          f = acos(g);
          if ((x2-x00)*(y1-y00)>(x1-x00)*(y2-y00))
            f = -f;
          if (fabs(g-1.0) > 0.0001) {
            cc = cos(f/segments);
            ss = sin(f/segments);
            x3 = x1;
            y3 = y1;
            for (i = 1; i <= segments; i++) {
              x4 = x00 + cc*(x3-x00) - ss*(y3-y00);
              y4 = y00 + ss*(x3-x00) + cc*(y3-y00);
              x3 = x4;
              y3 = y4;
              plot(pendown, x3, y3);
            }
          }
        }
        plot(pendown, x2, y2);
        break;

    }
  } else {
    if (style == circular) {
      x1 = x00;
      y1 = y00;
    } else {
      if (grows == vertical) {
        x1 =  xscale *  (xoffset + p->xcoord);
        y1 =  yscale *  (yoffset + rooty);
      } else {
        x1 =  xscale *  (xoffset + rootx);
        y1 =  yscale *  (yoffset + p->ycoord);
      }
    }
    //printf("root from: %f, %f to: %f, %f\n",x1, y1, x2, y2);
    //fflush(stdout); 
    plot(penup, x1, y1);
    plot(pendown, x2, y2);
  }
  if (p->tip)
    return;
  pp = p->next;
  while (pp != p) {
    plottree(pp->back, p);
    pp = pp->next;
  }
}  /* plottree */


void plotlabels(char *fontname)
{
  long i;
  double compr, dx = 0, dy = 0, labangle, cosl, sinl, cosv, sinv, vec;
  boolean left, right;
  node *lp;
  double *firstlet;

  firstlet = (double *)Malloc(nextnode*sizeof(double));
  textlength = (double *)Malloc(nextnode*sizeof(double));
  compr = xunitspercm / yunitspercm;
  if (penchange == yes)
    changepen(labelpen);
  for (i = 0; i < nextnode; i++) {
    if (curtree->nodep[i]->tip) {
      lp = curtree->nodep[i];
      firstlet[i] = lengthtext(curtree->nodep[i]->nayme,1L,fontname,font)
        /fontheight;
      textlength[i] = lengthtext(curtree->nodep[i]->nayme,
                                 curtree->nodep[i]->naymlength, fontname, font)/fontheight;
      labangle = ((draw_node*)curtree->nodep[i])->oldtheta;
      if (cos(labangle) < 0.0)
        labangle += pi;
      cosl = cos(labangle);
      sinl = sin(labangle);
      vec = sqrt(1.0+firstlet[i]*firstlet[i]);
      cosv = firstlet[i]/vec;
      sinv = 1.0/vec;
      if (style == circular) {
        right = cos(((draw_node*)curtree->nodep[i])->oldtheta) > 0.0;
        left = !right;
        if (right) {
          dx = labelheight * expand * cos(((draw_node*)curtree->nodep[i])->oldtheta);
          dy = labelheight * expand * sin(((draw_node*)curtree->nodep[i])->oldtheta);
          dx += labelheight * expand * 0.5 * vec * (-cosl*sinv+sinl*cosv);
          dy += labelheight * expand * 0.5 * vec * (-sinl*sinv-cosl*cosv);
        }
        if (left) {
          dx = labelheight * expand * cos(((draw_node*)curtree->nodep[i])->oldtheta);
          dy = labelheight * expand * sin(((draw_node*)curtree->nodep[i])->oldtheta);
          dx -= labelheight * expand * textlength[i] * cosl;
          dy -= labelheight * expand * textlength[i] * sinl;
          dx += labelheight * expand * 0.5 * vec * (cosl*cosv+sinl*sinv);
          dy += labelheight * expand * 0.5 * vec * (-sinl*cosv-cosl*sinv);
        }
      } else {
        dx = labelheight * expand * cos(((draw_node*)curtree->nodep[i])->oldtheta);
        dy = labelheight * expand * sin(((draw_node*)curtree->nodep[i])->oldtheta);
        dx -= labelheight * expand * 0.5 * firstlet[i] * (cosl-sinl*sinv);
        dy -= labelheight * expand * 0.5 * firstlet[i] * (sinl+cosl*sinv);
      }
      if (style == circular) {
        plottext(lp->nayme, lp->naymlength,
                 labelheight * expand * xscale / compr, compr,
                 xscale * (lp->xcoord + dx + xoffset),
                 yscale * (lp->ycoord + dy + yoffset), 180 * (-labangle) / pi,
                 font,fontname);
      } else {
        if (grows == vertical)
          plottext(lp->nayme, lp->naymlength,
                   labelheight * expand * xscale / compr, compr,
                   xscale * (lp->xcoord + dx + xoffset),
                   yscale * (lp->ycoord + dy + yoffset),
                   -labelrotation, font,fontname);
        else
          plottext(lp->nayme, lp->naymlength, labelheight *
                   expand * yscale, compr, xscale * (lp->xcoord + dy + xoffset),
                   yscale * (lp->ycoord - dx + yoffset), 90.0 - labelrotation,
                   font,fontname);
      }
    }
  }
  if (penchange == yes)
    changepen(treepen);
  free(firstlet);
  free(textlength);
}  /* plotlabels */


void setup_environment(int argc, Char *argv[])
{
  boolean firsttree;
  (void)argc;                           // RSGnote: Parameter never used.

  /* Set up all kinds of fun stuff */
#ifdef MAC
  OSErr retcode;
  FInfo  fndrinfo;
  macsetup("Drawgram","Preview");
#endif
#ifdef TURBOC
  if ((registerbgidriver(EGAVGA_driver) <0) ||
      (registerbgidriver(Herc_driver) <0)   ||
      (registerbgidriver(CGA_driver) <0)) {
    printf("Graphics error: %s ",grapherrormsg(graphresult()));
    exit(-1);
  }
#endif

  if (javarun)
  {
    intree = fopen(trefilename,"r");
  }
  else
  {
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree, INTREE, "input tree file", "rb", argv[0], trefilename);
    printf("DRAWGRAM from PHYLIP version %s\n", VERSION);
    printf("Reading tree ... \n");
  }
  firsttree = true;
  inputNumbersFromTreeFile(intree,&spp,&nonodes);
  curtree = functions.tree_new(nonodes,spp);
  treeread (curtree,intree, &curtree->root, curtree->nodep, &goteof, &firsttree, &nextnode,
            &haslengths, initdrawgramnode, true, -1);
  curtree->root->oldlen = 0.0;

  // if there are no lengths, shut uselengths off
  if (haslengths==0)
    uselengths = false;     
    
  if (uselengths)
    nodeposition = weighted;
  else
    nodeposition = centered;
  
  //printf("in setup haslengths: %i uselengths: %i nodeposition: %i\n", haslengths, uselengths, nodeposition);

  if (!javarun)
  {
    printf("Tree has been read.\n");
    printf("Loading the font .... \n");
    loadfont(font,argv[0]);
    printf("Font loaded.\n");
  }
  previewing = false;
  ansi = ANSICRT;
  ibmpc = IBMCRT;
  firstscreens = true;
  initialparms();
  //printf("after initialparms haslengths: %i uselengths: %i nodeposition: %i\n", haslengths, uselengths, nodeposition);

  canbeplotted = false;
}  /* setup_environment */


void user_loop(void)
{
  char     input_char;
  long stripedepth;

  while (!canbeplotted) {
    do {
      input_char=showparms();
      firstscreens = false;
      if (input_char != 'Y')
        getparms(input_char);
    } while (input_char != 'Y');

#ifdef JAVADEBUG
    // JRM Debug
    printf("\n***internal parameters***\n");
    printf("fontname: %s\n", fontname);
    printf("plotter: %i\n",plotter);
    printf("previewer: %i\n",previewer);
    printf("paperx: %f\n",paperx);
    printf("papery: %f\n",papery);
    printf("hpmargin: %f\n",hpmargin);
    printf("vpmargin: %f\n",vpmargin);
    printf("pagex: %f\n",pagex);
    printf("pagey: %f\n",pagey);
    printf("xmargin: %f\n",xmargin);
    printf("ymargin: %f\n",ymargin);
    printf("previewing: %i\n",previewing);
    printf("firstscreens: %i\n",firstscreens);
    printf("canbeplotted: %i\n",canbeplotted);
    printf("style: %i\n",style);
    printf("grows: %i\n",grows);
    printf("labelrotation: %f\n",labelrotation);
    printf("nodespace: %f\n",nodespace);
    printf("stemlength: %f\n",stemlength);
    printf("treedepth: %f\n",treedepth);
    printf("bscale: %f\n",bscale);
    printf("uselengths: %i\n",uselengths);
    printf("nodeposition: %i\n",nodeposition);
    printf("dotmatrix: %i\n",dotmatrix);
    fflush(stdout);
#endif

    if (dotmatrix) {
      stripedepth = allocstripe(stripe,(strpwide/8),
                                ((long)(yunitspercm * ysize)));
      strpdeep = stripedepth;
      strpdiv  = stripedepth;
    }
    plotrparms(spp);
    numlines = dotmatrix ?
      ((long)floor(yunitspercm * ysize + 0.5) / strpdeep) :1;
    xscale = xunitspercm;
    yscale = yunitspercm;
    calculate();
    rescale();
    canbeplotted = true;
    if (!javarun)
    {
      if (preview) {
        previewing = true;
        canbeplotted = plotpreview(fontname,&xoffset,&yoffset,&scale,spp,curtree->root);
      } else
        canbeplotted = true;

      if ((previewer == winpreview || previewer == xpreview || previewer == mac)
          && (winaction == quitnow)) {
        canbeplotted = true;
      }
    }
  }
} /* user_loop */


void drawgram(
  char   *intreename,
  char   *fontfilename,
  char   *plotfilename,
  char   *plotfileopt,
  char   *treegrows,
  char   *stylekind,
  int     usebranchlengths,
  double  labelangle,
  int     scalebranchlength,
  double  branchscale,
  double  breadthdepthratio,
  double  stemltreedratio,
  double  chhttipspratio,
  double  xmarginratio,
  double  ymarginratio,
  char   *ancnodes,
  int     dofinalplot,
  char*   finalplotkind)

{
    //printf("hello from DrawGram\n");

    initdata* funcs;
    int argc = 1;
    Char *argv[1];
    argv[0] = "Drawgram";
    strcpy(progname, argv[0]);
    funcs = Malloc(sizeof(initdata));
    funcs->node_new = draw_node_new;
    phylipinit(argc, argv, funcs, true);

    // read in the tree and set up things that depend on it
    strcpy(trefilename, intreename); // need this for setup_environment

    setup_environment(argc, argv);

    previewer = lw;  // hardwired to ps
    previewing = false;
    firstscreens = false;
    canbeplotted = true;

    char *previewfilename = "JavaPreview.ps";
    previewer = lw;  // hardwired to ps

    // translate Java doubles to Drawgram local variables
    labelrotation = labelangle;
    stemlength = stemltreedratio;
    treedepth = breadthdepthratio;
    bscale = branchscale;
    nodespace = 1.0 / chhttipspratio;
    xmargin = xmarginratio * paperx;
    ymargin = ymarginratio * papery;

    // set both font strings the same since ps doesn't need Hershey fonts
    strcpy(fontname, fontfilename);

    // translate Java boolean to Phylip	boolean
    boolean doplot;
    if (dofinalplot != 0)
        doplot = true;
    else
        doplot = false;

    //printf("after treeread haslengths: %i\n",haslengths);
    if (haslengths==0){
            uselengths = false; 
        }
    else{
        if (usebranchlengths != 0)
            uselengths = true;
        else
            uselengths = false;
    }    

    if (scalebranchlength != 0)
        rescaled = true;
    else
        rescaled = false;

    // plot style
    style = phenogram;
    if (!strcmp(stylekind,"cladogram"))
        style = cladogram;
    if (!strcmp(stylekind,"phenogram"))
        style = phenogram;
    if (!strcmp(stylekind,"curvogram"))
        style = curvogram;
    if (!strcmp(stylekind,"eurogram"))
        style = eurogram;
    if (!strcmp(stylekind,"swoopogram"))
        style = swoopogram;
    if (!strcmp(stylekind,"circular"))
        style = circular;

    // tree growth direction
    grows = horizontal;
    if (!strcmp(treegrows,"vertical"))
        grows = vertical;

    // node positions
    if (!strcmp(ancnodes,"weighted"))
        nodeposition = weighted;
    if (!strcmp(ancnodes,"intermediate"))
        nodeposition = intermediate;
    if (!strcmp(ancnodes,"centered"))
        nodeposition = centered;
    if (!strcmp(ancnodes,"inner"))
        nodeposition = inner;
    if (!strcmp(ancnodes,"vshaped"))
        nodeposition = vshaped;

  plotter = lw;
  if(dofinalplot)
  {
    if (!strcmp(finalplotkind,"lw"))
      plotter = lw;
    if(!strcmp(finalplotkind,"svg"))
      plotter = svg;
    /*
    // not currently available
    if(!strcmp(finalplotkind,"hp"))
    plotter = hp;
    if(!strcmp(finalplotkind,"tek"))
    plotter = tek;
    if(!strcmp(finalplotkind,"ibm"))
    plotter = ibm;
    if(!strcmp(finalplotkind,"mac"))
    plotter = mac;
    if(!strcmp(finalplotkind,"houston"))
    plotter = houston;
    if(!strcmp(finalplotkind,"decregis"))
    plotter = decregis;
    if(!strcmp(finalplotkind,"epson"))
    plotter = epson;
    if(!strcmp(finalplotkind,"oki"))
    plotter = oki;
    if(!strcmp(finalplotkind,"fig"))
    plotter = fig;
    if(!strcmp(finalplotkind,"citoh"))
    plotter = citoh;
    if(!strcmp(finalplotkind,"toshiba"))
    plotter = toshiba;
    if(!strcmp(finalplotkind,"pcx"))
    plotter = pcx;
    if(!strcmp(finalplotkind,"pcl"))
    plotter = pcl;
    if(!strcmp(finalplotkind,"pict"))
    plotter = pict;
    if(!strcmp(finalplotkind,"ray"))
    plotter = ray;
    if(!strcmp(finalplotkind,"pov"))
    plotter = pov;
    if(!strcmp(finalplotkind,"xbm"))
    plotter = xbm;
    if(!strcmp(finalplotkind,"bmp"))
    plotter = bmp;
    if(!strcmp(finalplotkind,"gif"))
    plotter = gif;
    if(!strcmp(finalplotkind,"idraw"))
    plotter = idraw;
    if(!strcmp(finalplotkind,"vrml"))
    plotter = vrml;
    */
  }
  else
  {
    plotter   = lw;  // hardwired to ps
  }
/*
  printf("***before setup_environment\n");
  printf("firstscreens: %i\n",firstscreens);
  printf("canbeplotted: %i\n",canbeplotted);
  printf("nodeposition: %i\n",nodeposition);

  setup_environment(argc, argv);
  
  previewer = lw;  // hardwired to ps
  previewing = false;
  firstscreens = false;
  canbeplotted = true;
*/

#ifdef JAVADEBUG
  // JRM Debug
  // dump translated parameters from Java to make sure they made it
  printf("Starting parameters: \n");
  printf("intreename: %s\n",intreename);
  printf("fontfilename: %s\n",fontfilename);
  printf("plotfilename: %s\n",plotfilename);
  printf("plotfileopt: %s\n",plotfileopt);
  printf("previewfilename: %s\n",previewfilename);
  printf("plotter: %i\n",plotter);
  printf("grows: %i\n",grows);
  printf("stylekind: %s\n", stylekind);
  printf("style: %i\n",style);
  printf("uselengths: %i\n",uselengths);
  printf("labelrotation: %f\n",labelrotation);
  printf("rescaled: %i\n",rescaled);
  printf("bscale: %f\n",bscale);
  printf("treedepth: %f\n",treedepth);
  printf("stemlength: %f\n",stemlength);
  printf("ancnodes: %s\n", ancnodes);
  printf("nodeposition: %i\n",nodeposition);
  printf("fontname: %s\n", fontname);
  printf("paperx: %f\n",paperx);
  printf("papery: %f\n",papery);
  printf("hpmargin: %f\n",hpmargin);
  printf("vpmargin: %f\n",vpmargin);
  printf("pagex: %f\n",pagex);
  printf("pagey: %f\n",pagey);
  printf("xmarginratio: %f\n",xmarginratio);
  printf("ymarginratio: %f\n",ymarginratio);
  printf("xmargin: %f\n",xmargin);
  printf("ymargin: %f\n",ymargin);
  printf("finalplotkind: %s\n",finalplotkind);

  printf("\n***internal parameters***\n");
  printf("fontname: %s\n", fontname);
  printf("plotter: %i\n",plotter);
  printf("previewer: %i\n",previewer);
  printf("paperx: %f\n",paperx);
  printf("papery: %f\n",papery);
  printf("hpmargin: %f\n",hpmargin);
  printf("vpmargin: %f\n",vpmargin);
  printf("pagex: %f\n",pagex);
  printf("pagey: %f\n",pagey);
  printf("xmargin: %f\n",xmargin);
  printf("ymargin: %f\n",ymargin);
  printf("previewing: %i\n",previewing);
  printf("firstscreens: %i\n",firstscreens);
  printf("canbeplotted: %i\n",canbeplotted);
  printf("style: %i\n",style);
  printf("grows: %i\n",grows);
  printf("labelrotation: %f\n",labelrotation);
  printf("nodespace: %f\n",nodespace);
  printf("stemlength: %f\n",stemlength);
  printf("treedepth: %f\n",treedepth);
  printf("bscale: %f\n",bscale);
  printf("uselengths: %i\n",uselengths);
  printf("nodeposition: %i\n",nodeposition);
  printf("dotmatrix: %i\n",dotmatrix);
  printf("xsize: %f\n", xsize);
  printf("ysize: %f\n", ysize);
  fflush(stdout);
#endif

  //setup_environment(argc, argv);

  // start drawgram code
  boolean wasplotted = false;
  numlines = 1;
  xscale = xunitspercm;
  yscale = yunitspercm;

  //printf("calling calculate\n");
  //fflush(stdout);
  calculate();
  //printf("calling rescale\n");
  //fflush(stdout);
  rescale();
  //printf("back from rescale\n");
  //fflush(stdout);

  int lenname;
  if (doplot)
  {
    lenname = strlen(plotfilename);
  }
  else
  {
    lenname = strlen(previewfilename);
  }

  char* usefilename = (char*) malloc(lenname + 1);

  if (doplot)
  {
    strcpy(usefilename, plotfilename);
  }
  else
  {
    strcpy(usefilename, previewfilename);
  }

  // Open plot files in binary mode.
  //printf("opening plotfile\n");
  //printf("usefilename: %s\n",usefilename);
  //printf("plotfileopt: %s\n",plotfileopt);
  //fflush(stdout);
  plotfile = fopen(usefilename, plotfileopt);

  previewing = false;
  //printf("calling initplotter\n");
  //fflush(stdout);
  initplotter(spp,fontname);
  numlines = dotmatrix ? ((long)floor(yunitspercm * ysize + 0.5)/strpdeep) : 1;
  
  if (!dofinalplot)
  {
      changepen(labelpen);
      makebox(fontname,&xoffset,&yoffset,&scale,spp);
      changepen(treepen);
  }

  //printf("calling drawit\n");
  //fflush(stdout);
  drawit(fontname,&xoffset,&yoffset,numlines,curtree->root);
  finishplotter();
  fclose(plotfile);

  wasplotted = true;
  if (!javarun)
  {
    printf("\nPlot written to file \"%s\".\n\n", pltfilename);
  }
  fclose(intree);
  //printf("Done.\n\n");
  return;
}


int main(int argc, Char *argv[])
{
  boolean wasplotted = false;
  initdata* funcs;
#ifdef MAC
  OSErr retcode;
  FInfo  fndrinfo;
#ifdef OSX_CARBON
  FSRef fileRef;
  FSSpec fileSpec;
#endif
#ifdef __MWERKS__
  SIOUXSetTitle("\pPHYLIP:  Drawtree");
#endif
  argv[0] = "Drawgram";
#endif

  strcpy(progname, argv[0]);

#ifdef PHYLIP_CAN_DRAW_WITH_X
  nargc=argc;
  nargv=argv;
#endif

  funcs = Malloc(sizeof(initdata));
  funcs->node_new = draw_node_new;
  phylipinit(argc, argv, funcs, false);

  setup_environment(argc, argv);

  user_loop();
  if (!((previewer == winpreview || previewer == xpreview || previewer == mac)
        && (winaction == quitnow))) {
    /* Open plot files in binary mode. */
    openfile(&plotfile,PLOTFILE,"plot file", "wb",argv[0],pltfilename);
    previewing = false;
    initplotter(spp,fontname);
    numlines = dotmatrix ? ((long)floor(yunitspercm * ysize + 0.5)/strpdeep) : 1;
    if (plotter != ibm)
      printf("\nWriting plot file ...\n");
    drawit(fontname,&xoffset,&yoffset,numlines,curtree->root);
    finishplotter();
    FClose(plotfile);
    wasplotted = true;
    printf("\nPlot written to file \"%s\".\n\n", pltfilename);
  }
  FClose(intree);
#ifdef MAC
  if (plotter == pict && wasplotted) {
#ifdef OSX_CARBON
    FSPathMakeRef((unsigned char *)pltfilename, &fileRef, NULL);
    FSGetCatalogInfo(&fileRef, kFSCatInfoNone, NULL, NULL, &fileSpec, NULL);
    FSpGetFInfo(&fileSpec, &fndrinfo);
    fndrinfo.fdType='PICT';
    fndrinfo.fdCreator='MDRW';
    FSpSetFInfo(&fileSpec, &fndrinfo);
#else
    retcode=GetFInfo(CtoPstr(PLOTFILE),0,&fndrinfo);
    fndrinfo.fdType='PICT';
    fndrinfo.fdCreator='MDRW';
    retcode=SetFInfo(CtoPstr(PLOTFILE),0,&fndrinfo);
#endif
  }
  if ((plotter == lw || plotter == svg) && wasplotted) {
#ifdef OSX_CARBON
    FSPathMakeRef((unsigned char *)pltfilename, &fileRef, NULL);
    FSGetCatalogInfo(&fileRef, kFSCatInfoNone, NULL, NULL, &fileSpec, NULL);
    FSpGetFInfo(&fileSpec, &fndrinfo);
    fndrinfo.fdType='TEXT';
    FSpSetFInfo(&fileSpec, &fndrinfo);
#else
    retcode=GetFInfo(CtoPstr(PLOTFILE),0,&fndrinfo);
    fndrinfo.fdType='TEXT';
    retcode=SetFInfo(CtoPstr(PLOTFILE),0,&fndrinfo);
#endif
  }
#endif
  printf("Done.\n\n");

  phyRestoreConsoleAttributes();

  return 0;
}


// End.
