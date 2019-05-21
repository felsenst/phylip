/* Version 4.0.  Copyright (c) 1986-2013 by the University of Washington and
   Written by Joseph Felsenstein and Christopher A. Meacham.  Additional code
   written by Sean Lamont, Andrew Keefe, Hisashi Horino, Akiko Fuseki, Doug
   Buxton and Michal Palczewski.  Permission is granted to copy, distribute,
   and modify this program provided that (1) this copyright message is
   not removed and (2) no fee is charged for this program. */


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
  "Drawtree unrooted tree plotting program\r"
  "PHYLIP version 4.0 (c) Copyright 1986-2013\r"
  "by The University of Washington.\r"
  "Written by Joseph Felsenstein and Christopher A. Meacham.\r"
  "Additional code written by Sean Lamont, Andrew Keefe, Hisashi Horino,\r"
  "Akiko Fuseki, Doug Buxton and Michal Palczewski.\r"
  "Permission is granted to copy, distribute and modify this program\r"
  "provided that\r"
  "(1) This copyright message is not removed and\r"
  "(2) no fee is charged for this program.";
#endif

//#define JAVADEBUG

char trefilename[FNMLNGTH];

tree * curtree;
long nonodes;

#define GAP             0.5
#define MAXITERATIONS   100
#define MINIMUMCHANGE   0.0001

/*When 2 Nodes are on top of each other, this is the max.
  force that's allowed.*/
#ifdef INFINITY
#undef INFINITY
#endif
#define INFINITY        (double) 9999999999.0

typedef enum {fixed, radial, along, middle} labelorient;

long        nextnode, payge, numlines;
double      topoflabels, rightoflabels, leftoflabels,
  bottomoflabels, ark, maxx, maxy, minx, miny, scale,
  xoffset, yoffset, charht, treeangle, bscale, maxchange;
boolean     canbeplotted, haslengths,
  uselengths, regular, rotate,  rescaled,
  notfirst, improve, nbody, firstscreens, labelavoid;

labelorient labeldirec;
node *where;
fonttype font;
char ch;
double *textlength, *firstlet;
double trweight;   /* starting here, needed to make sccs version happy */
boolean goteof;

long maxNumOfIter;
struct stackElem
{
  /* This is actually equivalent to a reversed link list; pStackElemBack
     point toward the direction of the bottom of the stack */
  struct stackElem *pStackElemBack;
  node *pNode;
};
typedef struct stackElem stackElemType;

typedef struct drawtree_node {
  draw_node draw_node;
  long depth;
  double r, width, lefttheta, righttheta;
} drawtree_node;

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
          "Drawtree unrooted tree plotting program\\n"
  "PHYLIP version 4.0 (c) Copyright 1986-2013\\n"
  "by The University of Washington.\\n"
  "Written by Joseph Felsenstein and Christopher A. Meacham.\\n"
  "Additional code written by Sean Lamont, Andrew Keefe, Hisashi Horino,\\n"
  "Akiko Fuseki, Doug Buxton and Michal Palczewski.\\n"
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
void   getwidth(node *);
void   plrtrans(node *, double, double, double);
void   coordtrav(node *, double *, double *);
double angleof(double , double );
void   polartrav(node *, double, double, double, double, double *,
        double *, double *, double *);
void   tilttrav(node *, double *, double *, double *, double *);

void   leftrightangle(node *, double, double);
void   improvtrav(node *);
void   force_1to1(node *, node *, double *, double *, double);
void   totalForceOnNode(node *, node *, double *, double *, double);
double dotProduct(double, double, double, double );
double capedAngle(double);
double angleBetVectors(double, double, double, double);
double signOfMoment(double, double, double, double);
double forcePerpendicularOnNode(node *, node *, double);
void   polarizeABranch(node *, double *, double *);

void   pushNodeToStack(stackElemType **, node *);
void   popNodeFromStack(stackElemType **, node **);
double medianOfDistance(node *, boolean);
void   leftRightLimits(node *, double *, double *);
void   branchLRHelper(node *, node *, double *, double *);
void   branchLeftRightAngles(node *, double *,  double *);
void   improveNodeAngle(node *, double);
void   improvtravn(node *);
void   coordimprov(double *, double *);
void   calculate(void);
void   rescale(void);
void   user_loop(void);
void setup_environment(int argc, Char *argv[]);
void polarize(node *p, double *xx, double *yy);
double vCounterClkwiseU(double Xu, double Yu, double Xv, double Yv);
node* drawtree_node_new(node_type type, long index);
void drawtree_node_init(node *n, node_type type, long index);
void drawtree(char* intreename, char* plotfilename, char* plotfileopt, char* fontfilename,
              char* treegrows, int usebranchlengths, char* labeldirecstr, double labelangle,
              double treerotation, double treearc, char* iterationkind, int iterationcount,
              int regularizeangles, int avoidlabeloverlap, int branchrescale,
              double branchscaler, double relcharhgt, double xmarginratio,
              double ymarginratio, int dofinalplot, char* finalplotkind);
void   makebox(char *, double *, double *, double *, long);
/* function prototypes */
#endif


node* drawtree_node_new(node_type type, long index)
{
  node* n = Malloc(sizeof(drawtree_node));
  drawtree_node_init(n, type, index);
  return n;
}


void drawtree_node_init(node *n, node_type type, long index)
{
  drawtree_node *dtn = (drawtree_node *)n;

  draw_node_init(n, type, index);
  n->init = drawtree_node_init;
  dtn->r = 0;
}


void initialparms(void)
{
  /* initialize parameters */
  //printf("in initialparms\n");
  //   fflush(stdout);
  paperx = 20.6375;
  pagex  = 20.6375;
  papery = 26.9875;
  pagey  = 26.9875;
  strcpy(fontname, "Times-Roman");
  plotrparms(spp);
  grows = vertical;
  treeangle = pi / 2.0;
  ark = 2 * pi;
  improve = true;
  nbody = false;
  regular = false;
  rescaled = true;
  bscale = 1.0;
  labeldirec = middle;
  xmargin = 0.08 * xsize;
  ymargin = 0.08 * ysize;
  labelrotation = 0.0;
  charht = 0.3333;
  preview = true;
  hpmargin = 0.02*pagex;
  vpmargin = 0.02*pagey;
  labelavoid = false;
  uselengths = haslengths;
}  /* initialparms */


char showparms(void)
{
  long loopcount;
  char numtochange;
  Char ch, input[64];
  double treea;
  char    options[32];

  strcpy(options, "#YN0OPVBLRIDSMC");
  if (strcmp(fontname, "Hershey") !=0 &&
      (((plotter == pict || plotter == mac)  &&
        (((grows == vertical && labelrotation == 0.0) ||
          (grows == horizontal && labelrotation == 90.0))))))
    strcat(options, "Q");
  if  (plotter == lw || plotter == idraw || plotter == pict || plotter == mac
       || plotter == svg)
    strcat(options, "F");
  if (!improve)
    strcat(options, "GA");
  if (!firstscreens)
    clearit();
  printf("\nUnrooted tree plotting program version %s\n", VERSION);
  putchar('\n');
  printf("Here are the settings: \n\n");
  printf(" 0  Screen type (IBM PC, ANSI)?  %s\n",
         ibmpc ? "IBM PC" : ansi  ? "ANSI"  : "(none)");
  printf(" P       Final plotting device: ");
  switch (plotter) {
    case lw:
      printf(" Postscript printer\n");
      break;
    case svg:
      printf(" SVG Scalable Vector Graphics\n");
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
    case bmp:
      printf(" MS-Windows Bitmap (%d by %d resolution)\n",
             (int)xsize, (int)ysize);
      break;
    case xbm:
      printf(" X Bitmap file format (%d by %d resolution)\n",
             (int)xsize, (int)ysize);
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
    case mac:
    case gif:
    case other:
      break ;
    default:        /* case xpreview not handled        */
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
        printf(" IBM PC graphics (CGA, EGA, or VGA)\n");
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
      case lw:
      case svg:
      case hp:
      case houston:
      case epson:
      case oki:
      case fig:
      case citoh:
      case toshiba:
      case pcx:
      case pcl:
      case pict:
      case ray:
      case pov:
      case bmp:
      case xbm:
      case gif:
      case idraw:
      case other:
        break ;
      default:        /* case vrml not handled        */
        break;
    }
  }
#endif
  printf(" B          Use branch lengths:  ");
  if (haslengths)
    printf("%s\n", uselengths ? "Yes" : "No");
  else
    printf("(no branch lengths available)\n");
  printf(" L             Angle of labels:");
  if (labeldirec == fixed) {
    printf("  Fixed angle of");
    if (labelrotation >= 10.0)
      printf("%6.1f", labelrotation);
    else if (labelrotation <= -10.0)
      printf("%7.1f", labelrotation);
    else if (labelrotation < 0.0)
      printf("%6.1f", labelrotation);
    else
      printf("%5.1f", labelrotation);
    printf(" degrees\n");
  } else if (labeldirec == radial)
    printf("  Radial\n");
  else if (labeldirec == middle)
    printf("  branch points to Middle of label\n");
  else
    printf("  Along branches\n");
  printf(" R            Rotation of tree:");
  treea = treeangle * 180 / pi;
  if (treea >= 100.0)
    printf("%7.1f\n", treea);
  else if (treea >= 10.0)
    printf("%6.1f\n", treea);
  else if (treea <= -100.0)
    printf("%8.1f\n", treea);
  else if (treea <= -10.0)
    printf("%7.1f\n", treea);
  else if (treea < 0.0)
    printf("%6.1f\n", treea);
  else
    printf("%5.1f\n", treea);
  if (!improve) {
    printf(" A       Angle of arc for tree:");
    treea = 180 * ark / pi;
    if (treea >= 100.0)
      printf("%7.1f\n", treea);
    else if (treea >= 10.0)
      printf("%6.1f\n", treea);
    else if (treea <= -100.0)
      printf("%8.1f\n", treea);
    else if (treea <= -10.0)
      printf("%7.1f\n", treea);
    else if (treea < 0.0)
      printf("%6.1f\n", treea);
    else
      printf("%5.1f\n", treea);
  }
  printf(" I     Iterate to improve tree:  ");
  if (improve) {
    if (nbody)
      printf("n-Body algorithm\n");
    else
      printf("Equal-Daylight algorithm\n");
  } else
    printf("No\n");
  if (improve)
    printf(" D  Try to avoid label overlap?  %s\n",
           (labelavoid? "Yes" : "No"));
  printf(" S      Scale of branch length:");
  if (rescaled)
    printf("  Automatically rescaled\n");
  else
    printf("  Fixed:%6.2f cm per unit branch length\n", bscale);
  if (!improve) {
    printf(" G       Regularize the angles:  %s\n",
           (regular ? "Yes" : "No"));
  }
  printf(" C   Relative character height:%8.4f\n", charht);
  if ((((plotter == pict || plotter == mac)  &&
        (((grows == vertical && labelrotation == 0.0) ||
          (grows == horizontal && labelrotation == 90.0))))))
    printf(" F                        Font:  %s\n Q"
           "        Pict Font Attributes:  %s, %s, %s, %s\n",
           fontname, (pictbold   ? "Bold"   : "Medium"),
           (pictitalic ? "Italic" : "Regular"),
           (pictshadow ? "Shadowed": "Unshadowed"),
           (pictoutline ? "Outlined" : "Unoutlined"));
  else if (plotter == lw || plotter == idraw || plotter == svg)
    printf(" F                        Font:  %s\n", fontname);
  if (plotter == ray) {
    printf(" M          Horizontal margins:%6.2f pixels\n", xmargin);
    printf(" M            Vertical margins:%6.2f pixels\n", ymargin);
  } else {
    printf(" M          Horizontal margins:%6.2f cm\n", xmargin);
    printf(" M            Vertical margins:%6.2f cm\n", ymargin);
  }
  printf(" #           Page size submenu:  ");
  /* Add 0.5 to clear up truncation problems. */
  if (((int) ((pagex / paperx) + 0.5) == 1) &&
      ((int) ((pagey / papery) + 0.5) == 1))
    /* If we're only using one page per tree, */
    printf ("one page per tree\n") ;
  else
    printf ("%.0f by %.0f pages per tree\n",
            (pagey-vpmargin) / (papery-vpmargin),
            (pagex-hpmargin) / (paperx-hpmargin)) ;
  loopcount = 0;
  for (;;) {
    printf("\n Y to accept these or type the letter for one to change\n");
    phyFillScreenColor();
    getstryng(input);
    uppercase(&input[0]);
    ch=input[0];
    if (strchr(options, ch))
    {
      numtochange = ch;
      break;
    }
    printf(" That letter is not one of the menu choices.  Type\n");
    countup(&loopcount, 100);
  }
  return numtochange;
}  /* showparms */


void getparms(char numtochange)
{
  /* get from user the relevant parameters for the plotter and diagram */
  long loopcount2;
  Char ch;
  boolean ok;
  char    options[32];
  char    line[32];
  char    input[100];
  int m, n;

  n = (int)((pagex-hpmargin-0.01)/(paperx-hpmargin)+1.0);
  m = (int)((pagey-vpmargin-0.01)/(papery-vpmargin)+1.0);
  strcpy(options, "YNOPVBLRIDSMC");
  if      ((((plotter == pict || plotter == mac)  &&
             (((grows == vertical && labelrotation == 0.0) ||
               (grows == horizontal && labelrotation == 90.0))))))
    strcat(options, "Q");
  if  (plotter == lw || plotter == idraw || plotter == svg)
    strcat(options, "F");
  if (!improve) {
    strcat(options, "GA");
  }
  if (numtochange == '*') {
    do {
      printf(" Type the number of one that you want to change:\n");
      phyFillScreenColor();
      getstryng(line);
      numtochange = line[0];
    } while (strchr(options, numtochange));
  }
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

    case '#':
      loopcount2 = 0;
      for (;;) {
        clearit();
        printf("  Page Specifications Submenu\n\n");
        printf(" L   Output size in pages: %.0f down by %.0f across\n",
               (pagey / papery), (pagex / paperx));
        printf(" P   Physical paper size: %1.5f by %1.5f cm\n", paperx, papery);
        printf(" O   Overlap Region: %1.5f %1.5f cm\n", hpmargin, vpmargin);
        printf(" M   main menu\n");
        phyFillScreenColor();
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
            phyFillScreenColor();
            getstryng(input);
            n = atoi(input);
            break;
          case 'P':
            printf("Paper Width (in cm):\n");
            phyFillScreenColor();
            getstryng(input);
            paperx = atof(input);
            printf("Paper Height (in cm):\n");
            phyFillScreenColor();
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
        countup(&loopcount2, 20);
      }
      break;

    case 'B':
      if (haslengths)
        uselengths = !uselengths;
      else {
        printf("Cannot use lengths since not all of them exist\n");
        uselengths = false;
      }
      break;

    case 'L':
      printf("\nDo you want labels to be Fixed angle, Radial, Along,");
      printf(" or Middle?\n");
      loopcount2 = 0;
      do {
        printf(" Type F, R, A, or M\n");
        phyFillScreenColor();
        if(scanf("%c%*[^\n]", &ch)) {}  // Read char and scan to EOL.
        (void)getchar();
        if (ch == '\n')
          ch = ' ';
        uppercase(&ch);
        countup(&loopcount2, 10);
      } while (ch != 'F' && ch != 'R' && ch != 'A' && ch != 'M');
      switch (ch) {

        case 'A':
          labeldirec = along;
          break;

        case 'F':
          labeldirec = fixed;
          break;

        case 'R':
          labeldirec = radial;
          break;

        case 'M':
          labeldirec = middle;
          break;
      }
      if (labeldirec == fixed) {
        printf("Are the labels to be plotted vertically (90),\n");
        printf(" horizontally (0), or downwards (-90) ?\n");
        loopcount2 = 0;
        do {
          printf(" Choose an angle in degrees from 90 to -90: \n");
          phyFillScreenColor();
          if(scanf("%lf%*[^\n]", &labelrotation)) {} // Read number and scan to EOL.
          (void)getchar();
          countup(&loopcount2, 10);
        } while ((labelrotation < -90.0 || labelrotation > 90.0) &&
                 labelrotation != -99.0);
      }
      break;

    case 'R':
      printf("\n At what angle is the tree to be plotted?\n");
      loopcount2 = 0;
      do {
        printf(" Choose an angle in degrees from 360 to -360: \n");
        phyFillScreenColor();
        if(scanf("%lf%*[^\n]", &treeangle)) {} // Read number and scan to EOL.
        (void)getchar();
        uppercase(&ch);
        countup(&loopcount2, 10);
      } while (treeangle < -360.0 && treeangle > 360.0);
      treeangle = treeangle * pi / 180;
      break;

    case 'A':
      printf(" How many degrees (up to 360) of arc\n");
      printf("  should the tree occupy? (Currently it is %5.1f)\n",
             180 * ark / pi);
      loopcount2 = 0;
      do {
        printf("Enter a number of degrees from 0 up to 360)\n");
        phyFillScreenColor();
        if(scanf("%lf%*[^\n]", &ark)) {} // Read number and scan to EOL.
        (void)getchar();
        countup(&loopcount2, 10);
      } while (ark <= 0.0 || ark > 360.0);
      ark = ark * pi / 180;
      break;

    case 'I':
      if (nbody) {
        improve = false;
        nbody = false;
      } else {
        if (improve)
          nbody = true;
        else improve = true;
      }
      break;

    case 'D':
      labelavoid = !labelavoid;
      break;

    case 'S':
      rescaled = !rescaled;
      if (!rescaled) {
        printf("Centimeters per unit branch length?\n");
        phyFillScreenColor();
        if(scanf("%lf%*[^\n]", &bscale)) {} // Read number and scan to EOL.
        (void)getchar();
      }
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
      loopcount2 = 0;
      do {
        printf(" New value (in cm) of horizontal margin?\n");
        phyFillScreenColor();
        if(scanf("%lf%*[^\n]", &xmargin)) {} // Read number and scan to EOL.
        (void)getchar();
        ok = ((unsigned)xmargin < xsize / 2.0);
        if (!ok)
          printf(" Impossible value.  Please retype it.\n");
        countup(&loopcount2, 10);
      } while (!ok);
      loopcount2 = 0;
      do {
        printf(" New value (in cm) of vertical margin?\n");
        phyFillScreenColor();
        if(scanf("%lf%*[^\n]", &ymargin)) {} // Read number and scan to EOL.
        (void)getchar();
        ok = ((unsigned)ymargin < ysize / 2.0);
        if (!ok)
          printf(" Impossible value.  Please retype it.\n");
        countup(&loopcount2, 10);
      } while (!ok);
      break;

    case 'C':
      printf("New value of character height?\n");
      phyFillScreenColor();
      if(scanf("%lf%*[^\n]", &charht)) {} // Read number and scan to EOL.
      (void)getchar();
      break;

    case 'F':
      printf("Enter font name or \"Hershey\" for default font\n");
      phyFillScreenColor();
      getstryng(fontname);
      break;

    case 'G':
      regular = !regular;
      break;

    case 'Q':
      clearit();
      loopcount2 = 0;
      do {
        printf("Italic? (Y/N)\n");
        phyFillScreenColor();
        getstryng(input);
        input[0] = toupper(input[0]);
        countup(&loopcount2, 10);
      } while (input[0] != 'Y' && input[0] != 'N');
      pictitalic = (input[0] == 'Y');
      loopcount2 = 0;
      do {
        printf("Bold? (Y/N)\n");
        phyFillScreenColor();
        getstryng(input);
        input[0] = toupper(input[0]);
        countup(&loopcount2, 10);
      } while (input[0] != 'Y' && input[0] != 'N');
      pictbold = (input[0] == 'Y');
      loopcount2 = 0;
      do {
        printf("Shadow? (Y/N)\n");
        phyFillScreenColor();
        getstryng(input);
        input[0] = toupper(input[0]);
        countup(&loopcount2, 10);
      } while (input[0] != 'Y' && input[0] != 'N');
      pictshadow = (input[0] == 'Y');
      loopcount2 = 0;
      do {
        printf("Outline? (Y/N)\n");
        phyFillScreenColor();
        getstryng(input);
        input[0] = toupper(input[0]);
        countup(&loopcount2, 10);
      } while (input[0] != 'Y' && input[0] != 'N');
      pictoutline = (input[0] == 'Y');
      break;
  }
}  /* getparms */


void getwidth(node *p)
{
  /* get width and depth beyond each node */
  double nw, nd;
  node *pp, *qq;
  //seetree(p, nodep, nextnode); //JRMDebug
  //dumpnodelinks(p, nodep, nextnode); //JRMDebug

  //printf("  in getwidth p->type: %i nayme: %s index: %li\n", p->type, p->nayme, p->index); //JRMDebug

  nd = 0.0;
  if (p->tip)
    nw = 1.0;
  else {
    nw = 0.0;
    qq = p;
    pp = p->next;
    do {
      //printf("pp->back: %p\n", pp->back); //JRMDebug
      getwidth(pp->back);
      nw += ((drawtree_node*)pp->back)->width;
      if (((drawtree_node*)pp->back)->depth > nd)
        nd = ((drawtree_node*)pp->back)->depth;
      pp = pp->next;
    } while (((p != curtree->root) && (pp != qq)) || ((p == curtree->root) && (pp != p->next)));
  }
  ((drawtree_node*)p)->depth = nd + p->length;
  ((drawtree_node*)p)->width = nw;
}  /* getwidth */


void plrtrans(node *p, double theta, double lower, double upper)
{
  /* polar coordinates of a node relative to start */
  long num;
  double nn, pr, ptheta, angle, angle2, subangle, len;
  node *pp, *qq;

  // RSGnote: Parameter never used, but there is an object slot of same name.
  (void)theta;

  nn = ((drawtree_node*)p)->width;
  subangle = (upper - lower) / nn;
  qq = p;
  pp = p->next;
  if (p->tip)
    return;
  angle = upper;
  do {
    angle -= ((drawtree_node*)pp->back)->width / 2.0 * subangle;
    pr = ((drawtree_node*)p)->r;
    ptheta = ((draw_node*)p)->theta;
    if (regular) {
      num = 1;
      while (num * subangle < 2 * pi)
        num *= 2;
      if (angle >= 0.0)
        angle2 = 2 * pi / num * (long)(num * angle / (2 * pi) + 0.5);
      else
        angle2 = 2 * pi / num * (long)(num * angle / (2 * pi) - 0.5);
    } else
      angle2 = angle;
    if (uselengths)
      len = fabs(pp->back->oldlen);
    else
      len = 1.0;
    ((drawtree_node*)pp->back)->r = sqrt(len * len + pr * pr + 2 * len * pr *
                                         cos(angle2 - ptheta));
    if (fabs(pr * cos(ptheta) + len * cos(angle2)) > epsilon)
      ((draw_node*)pp->back)->theta = atan((pr * sin(ptheta) +
                                            len * sin(angle2)) / (pr *
                                                                  cos(ptheta) + len * cos(angle2)));
    else if (pr * sin(ptheta) + len * sin(angle2) >= 0.0)
      ((draw_node*)pp->back)->theta = pi / 2;
    else
      ((draw_node*)pp->back)->theta = 1.5 * pi;
    if (pr * cos(ptheta) + len * cos(angle2) < -epsilon)
      ((draw_node*)pp->back)->theta += pi;
    if (!pp->back->tip)
      plrtrans(pp->back, ((draw_node*)pp->back)->theta,
               angle - ((drawtree_node*)pp->back)->width * subangle / 2.0,
               angle + ((drawtree_node*)pp->back)->width * subangle / 2.0);
    else
      ((draw_node*)pp->back)->oldtheta = angle2;
    angle -= ((drawtree_node*)pp->back)->width / 2.0 * subangle;
    pp = pp->next;
  } while (((p != curtree->root) && (pp != qq)) || ((p == curtree->root) && (pp != p->next)));
}  /* plrtrans */


void coordtrav(node *p, double *xx, double *yy)
{
  /* compute x and y coordinates */
  node *pp;

  if (!p->tip) {
    pp = p->next;
    while (pp != p) {
      coordtrav(pp->back, xx, yy);
      pp = pp->next;
      if (p == curtree->root)
        coordtrav(p->back, xx, yy);
    }
  }
  (*xx) = ((drawtree_node*)p)->r * cos(((draw_node*)p)->theta);
  (*yy) = ((drawtree_node*)p)->r * sin(((draw_node*)p)->theta);
  if ((*xx) > maxx)
    maxx = (*xx);
  if ((*xx) < minx)
    minx = (*xx);
  if ((*yy) > maxy)
    maxy = (*yy);
  if ((*yy) < miny)
    miny = (*yy);
  p->xcoord = (*xx);
  p->ycoord = (*yy);
}  /* coordtrav */


double angleof(double x, double y)
{
  /* compute the angle of a vector */
  double theta;

  if (fabs(x) > epsilon)
    theta = atan(y / x);
  else if (y >= 0.0)
    theta = pi / 2;
  else
    theta = 1.5 * pi;
  if (x < -epsilon)
    theta = pi + theta;
  while (theta > 2 * pi)
    theta -= 2 * pi;
  while (theta < 0.0)
    theta += 2 * pi;
  return theta;
}  /* angleof */


void polartrav(node *p, double xx, double yy, double firstx,
               double firsty, double *leftx, double *lefty, double *rightx,
               double *righty)
{
  /* go through subtree getting left and right vectors */
  double x, y, xxx, yyy, labangle = 0;
  boolean lookatit;
  node *pp;

  lookatit = true;
  if (!p->tip)
    lookatit = (p->next->next != p || p->index != curtree->root->index);
  if (lookatit) {
    x = curtree->nodep[p->index - 1]->xcoord;
    y = curtree->nodep[p->index - 1]->ycoord;
    if (p->tip) {
      if (labeldirec == fixed) {
        labangle = pi * labelrotation / 180.0;
        if (cos(((draw_node*)p)->oldtheta) < 0.0)
          labangle = labangle - pi;
      }
      if (labeldirec == radial)
        labangle = ((draw_node*)p)->theta;
      else if (labeldirec == along)
        labangle = ((draw_node*)p)->oldtheta;
      else if (labeldirec == middle)
        labangle = 0.0;
      xxx = x;
      yyy = y;
      if (labelavoid) {
        if  (labeldirec == middle) {
          xxx += GAP * labelheight * cos(((draw_node*)p)->oldtheta);
          yyy += GAP * labelheight * sin(((draw_node*)p)->oldtheta);
          xxx += labelheight * cos(labangle) * textlength[p->index - 1];
          if (textlength[p->index - 1] * sin(((draw_node*)p)->oldtheta) < 1.0)
            xxx += labelheight * cos(labangle) * textlength[p->index - 1];
          else
            xxx += 0.5 * labelheight * cos(labangle)
              * textlength[p->index - 1];
          yyy += labelheight * sin(labangle) * textlength[p->index - 1];
        }
        else {
          xxx += GAP * labelheight * cos(((draw_node*)p)->oldtheta);
          yyy += GAP * labelheight * sin(((draw_node*)p)->oldtheta);
          xxx -= labelheight * cos(labangle) * 0.5 * firstlet[p->index - 1];
          yyy -= labelheight * sin(labangle) * 0.5 * firstlet[p->index - 1];
          xxx += labelheight * cos(labangle) * textlength[p->index - 1];
          yyy += labelheight * sin(labangle) * textlength[p->index - 1];
        }
      }
      if ((yyy - yy) * firstx - (xxx - xx) * firsty < 0.0) {
        if ((yyy - yy) * (*rightx) - (xxx - xx) * (*righty) < 0.0) {
          (*rightx) = xxx - xx;
          (*righty) = yyy - yy;
        }
      }
      if ((yyy - yy) * firstx - (xxx - xx) * firsty > 0.0) {
        if ((yyy - yy) * (*leftx) - (xxx - xx) * (*lefty) > 0.0) {
          (*leftx) = xxx - xx;
          (*lefty) = yyy - yy;
        }
      }
    }
    if ((y - yy) * firstx - (x - xx) * firsty < 0.0) {
      if ((y - yy) * (*rightx) - (x - xx) * (*righty) < 0.0) {
        (*rightx) = x - xx;
        (*righty) = y - yy;
      }
    }
    if ((y - yy) * firstx - (x - xx) * firsty > 0.0) {
      if ((y - yy) * (*leftx) - (x - xx) * (*lefty) > 0.0) {
        (*leftx) = x - xx;
        (*lefty) = y - yy;
      }
    }
  }
  if (p->tip)
    return;
  pp = p->next;
  while (pp != p) {
    if (pp != NULL)
      polartrav(pp->back, xx, yy, firstx, firsty, leftx, lefty, rightx, righty);
    pp = pp->next;
  }
}  /* polartrav */


void tilttrav(node *q, double *xx, double *yy, double *sinphi,
              double *cosphi)
{
  /* traverse to move successive nodes */
  double x, y;
  node *pp;

  pp = curtree->nodep[q->index - 1];
  x = pp->xcoord;
  y = pp->ycoord;
  pp->xcoord = (*xx) + (x - (*xx)) * (*cosphi) + ((*yy) - y) * (*sinphi);
  pp->ycoord = (*yy) + (x - (*xx)) * (*sinphi) + (y - (*yy)) * (*cosphi);
  if (q->tip)
    return;
  pp = q->next;
  while (pp != q) {
    /*
      if (pp != root)
    */
    if (pp->back != NULL)
      tilttrav(pp->back, xx, yy, sinphi, cosphi);
    pp = pp->next;
  }
}  /* tilttrav */


void polarize(node *p, double *xx, double *yy)
{
  double TEMP, TEMP1;

  if (fabs(p->xcoord - (*xx)) > epsilon)
    ((draw_node*)p)->oldtheta = atan((p->ycoord - (*yy)) /
                                     (p->xcoord - (*xx)));
  else if (p->ycoord - (*yy) > epsilon)
    ((draw_node*)p)->oldtheta = pi / 2;
  if (p->xcoord - (*xx) < -epsilon)
    ((draw_node*)p)->oldtheta += pi;
  if (fabs(p->xcoord - curtree->root->xcoord) > epsilon)
    ((draw_node*)p)->theta =
      atan((p->ycoord - curtree->root->ycoord) / (p->xcoord - curtree->root->xcoord));
  else if (p->ycoord - curtree->root->ycoord > 0.0)
    ((draw_node*)p)->theta = pi / 2;
  else
    ((draw_node*)p)->theta = 1.5 * pi;
  if (p->xcoord - curtree->root->xcoord < -epsilon)
    ((draw_node*)p)->theta += pi;
  TEMP = p->xcoord - curtree->root->xcoord;
  TEMP1 = p->ycoord - curtree->root->ycoord;
  ((drawtree_node*)p)->r = sqrt(TEMP * TEMP + TEMP1 * TEMP1);
}  /* polarize */


void leftrightangle(node *p, double xx, double yy)
{
  /* get leftmost and rightmost angle of subtree, put them in node p */
  double firstx, firsty, leftx, lefty, rightx, righty;
  double langle, rangle;

  firstx = curtree->nodep[p->back->index-1]->xcoord - xx;
  firsty = curtree->nodep[p->back->index-1]->ycoord - yy;
  leftx = firstx;
  lefty = firsty;
  rightx = firstx;
  righty = firsty;
  if (p->back != NULL)
    polartrav(p->back, xx, yy, firstx, firsty, &leftx, &lefty, &rightx, &righty);
  if ((fabs(leftx) < epsilon) && (fabs(lefty) < epsilon))
    langle = ((draw_node*)p->back)->oldtheta;
  else
    langle = angleof(leftx, lefty);
  if ((fabs(rightx) < epsilon) && (fabs(righty) < epsilon))
    rangle = ((draw_node*)p->back)->oldtheta;
  else
    rangle = angleof(rightx, righty);
  while (langle - rangle > 2*pi)
    langle -= 2 * pi;
  while (rangle > langle) {
    if (rangle > 2*pi)
      rangle -= 2 * pi;
    else
      langle += 2 * pi;
  }
  while (langle > 2*pi) {
    rangle -= 2 * pi;
    langle -= 2 * pi;
  }
  ((drawtree_node*)p)->lefttheta = langle;
  ((drawtree_node*)p)->righttheta = rangle;
} /* leftrightangle */


void improvtrav(node *p)
{
  /* traverse tree trying different tiltings at each node */
  double xx, yy, cosphi, sinphi;
  double langle, rangle, sumrot, olddiff;
  node *pp, *qq, *ppp;;

  if (p->tip)
    return;
  xx = p->xcoord;
  yy = p->ycoord;
  pp = p->next;
  do {
    leftrightangle(pp, xx, yy);
    pp = pp->next;
  } while ((pp != p->next));
  if (p == curtree->root) {
    pp = p->next;
    do {
      qq = pp;
      pp = pp->next;
    } while (pp != curtree->root);
    ((drawtree_node*)p)->righttheta = ((drawtree_node*)qq)->righttheta;
    ((drawtree_node*)p)->lefttheta = ((drawtree_node*)p->next)->lefttheta;
  }
  qq = p;
  pp = p->next;
  ppp = p->next->next;
  do {
    langle = ((drawtree_node*)qq)->righttheta -
      ((drawtree_node*)pp)->lefttheta;
    rangle = ((drawtree_node*)pp)->righttheta -
      ((drawtree_node*)ppp)->lefttheta;
    while (langle > pi)
      langle -= 2*pi;
    while (langle < -pi)
      langle += 2*pi;
    while (rangle > pi)
      rangle -= 2*pi;
    while (rangle < -pi)
      rangle += 2*pi;
    olddiff = fabs(langle-rangle);
    sumrot = (langle - rangle) /2.0;
    if (sumrot > langle)
      sumrot = langle;
    if (sumrot < -rangle)
      sumrot = -rangle;
    cosphi = cos(sumrot);
    sinphi = sin(sumrot);
    if (p != curtree->root) {
      if (fabs(sumrot) > maxchange)
        maxchange = fabs(sumrot);
      ((draw_node*)pp->back)->oldtheta += sumrot;
      tilttrav(pp->back, &xx, &yy, &sinphi, &cosphi);
      polarize(pp->back, &xx, &yy);
      leftrightangle(pp, xx, yy);
      langle = ((drawtree_node*)qq)->righttheta -
        ((drawtree_node*)pp)->lefttheta;
      rangle = ((drawtree_node*)pp)->righttheta -
        ((drawtree_node*)ppp)->lefttheta;
      while (langle > pi)
        langle -= 2*pi;
      while (langle < -pi)
        langle += 2*pi;
      while (rangle > pi)
        rangle -= 2*pi;
      while (rangle < -pi)
        rangle += 2*pi;
      while ((fabs(langle-rangle) > olddiff) && (fabs(sumrot) > 0.01)) {
        sumrot = sumrot /2.0;
        cosphi = cos(-sumrot);
        sinphi = sin(-sumrot);
        ((draw_node*)pp->back)->oldtheta -= sumrot;
        tilttrav(pp->back, &xx, &yy, &sinphi, &cosphi);
        polarize(pp->back, &xx, &yy);
        leftrightangle(pp, xx, yy);
        langle = ((drawtree_node*)qq)->righttheta -
          ((drawtree_node*)pp)->lefttheta;
        rangle = ((drawtree_node*)pp)->righttheta -
          ((drawtree_node*)ppp)->lefttheta;
        if (langle > pi)
          langle -= 2*pi;
        if (langle < -pi)
          langle += 2*pi;
        if (rangle > pi)
          rangle -= 2*pi;
        if (rangle < -pi)
          rangle += 2*pi;
      }
    }
    qq = pp;
    pp = pp->next;
    ppp = ppp->next;
  } while (((p == curtree->root) && (pp != p->next)) || ((p != curtree->root) && (pp != p)));
  pp = p->next;
  do {
    improvtrav(pp->back);
    pp = pp->next;
  }
  while (((p == curtree->root) && (pp != p->next)) || ((p != curtree->root) && (pp != p)));
}  /* improvtrav */


void force_1to1(node *pFromSubNode, node *pToSubNode, double *pForce,
                double *pAngle, double medianDistance)
{
  /* calculate force acting between 2 nodes and return the force in pForce.
     Remember to pass the index subnodes to this function if needed.
     Force should always be positive for repelling.  Angle changes to
     indicate the direction of the force.  The value of INFINITY is the cap
     to the value of Force.
     There might have problem (error msg.) if pFromSubNode and pToSubNode
     are the same node or the coordinates are identical even with double
     precision.  */
  double distanceX, distanceY, distance, norminalDistance;
  distanceX = pFromSubNode->xcoord - pToSubNode->xcoord;
  distanceY = pFromSubNode->ycoord - pToSubNode->ycoord;
  distance = sqrt( distanceX*distanceX + distanceY*distanceY );
  norminalDistance = distance/medianDistance;

  if (norminalDistance < epsilon)
  {
    *pForce = INFINITY;
  }
  else
  {

    *pForce = (double)1 / (norminalDistance * norminalDistance);

    if (*pForce > INFINITY) *pForce = INFINITY;
  }

  *pAngle = computeAngle(pFromSubNode->xcoord, pFromSubNode->ycoord,
                         pToSubNode->xcoord, pToSubNode->ycoord);

  return;
}  /* force_1to1 */


void totalForceOnNode(node *pPivotSubNode, node *pToSubNode,
                      double *pTotalForce, double *pAngle, double medianDistance)
{
  /* pToSubNode is where all the relevent nodes apply forces to.
     All branches are visited except the branch contains pToSubNode.
     pToSubNode must be one of the branch out of the current Node (= out of one
     of the subnode in the current subnodes set.)
     Most likely pPivotSubNode is not the index subNode!  In any case, only
     the leafs are consider in repelling force; so, no worry about index subNode.
     pTotalForce and pAngle must be set to 0 before calling this function for the
     first time, or the result will be invalid.  pPivotSubNode is named for
     external interface.  When calling totalForceOnNode() recursively,
     pPivotSubNode should be thought of as pFromSubNode.
  */
  node *pSubNode;
  double force, angle, forceX, forceY, prevForceX, prevForceY;

  pSubNode = pPivotSubNode;

  /* visit the rest of the branches of current node; the branch attaches to
     the current subNode may be visited in the code down below. */
  while (pSubNode->next != NULL && pSubNode->next != pPivotSubNode)
  {
    pSubNode = pSubNode->next;
    if ( pSubNode->back != NULL && pSubNode->back != pToSubNode)
      totalForceOnNode(pSubNode->back, pToSubNode, pTotalForce, pAngle,
                       medianDistance);
  }

  /* visit this branch; You need to visit it for the first time - at root only!
   *
   * Modified so that all nodes are visited and calculated forces, instead of
   * just the leafs only.
   * use pPivotSubNode instead of pSubNode here because pSubNode stop short
   *  just before pPivotSubNode (the entry node) */
  if ( pPivotSubNode == curtree->root && pPivotSubNode->back != NULL
       && pPivotSubNode->back != pToSubNode)
    totalForceOnNode(pPivotSubNode->back, pToSubNode, pTotalForce, pAngle,
                     medianDistance);

  /* Break down the previous sum of forces to components form */
  prevForceX = *pTotalForce * cos(*pAngle);
  prevForceY = *pTotalForce * sin(*pAngle);

  force_1to1(curtree->nodep[pPivotSubNode->index-1], pToSubNode, &force, &angle,
             medianDistance);
  /* force between 2 nodes */
  forceX = force * cos(angle);
  forceY = force * sin(angle);

  /* Combined force */
  forceX = forceX + prevForceX;
  forceY = forceY + prevForceY;

  /* Write to output parameters */
  *pTotalForce = sqrt( forceX*forceX + forceY*forceY );

  *pAngle = computeAngle((double)0, (double)0, forceX, forceY);

  return;
}  /* totalForceOnNode */


double dotProduct(double Xu, double Yu, double Xv, double Yv)
{
  return Xu * Xv + Yu * Yv;
}  /* dotProduct */


double capedAngle(double angle)
{
/* Return the equivalent value of angle that is within
   0 to 2*pi */
  while (angle < 0 || angle >= 2*pi)
  {
    if(angle < 0)
    {
      angle = angle + 2*pi;
    }
    else if (angle >= 2*pi)
    {
      angle = angle - 2*pi;
    }
  }
  return angle;
}  /* capedAngle */


double angleBetVectors(double Xu, double Yu, double Xv, double Yv)
{
  /* Calculate angle between 2 vectors; value returned is always between
     0 and pi.
     Use vCounterClkwiseU() to get the relative position of the vectors.  */
  double dotProd, cosTheta, theta, lengthsProd;

  dotProd = dotProduct(Xu, Yu, Xv, Yv);
  lengthsProd = sqrt(Xu*Xu+Yu*Yu) * sqrt(Xv*Xv+Yv*Yv);

  if (lengthsProd < epsilon)
  {
    printf("ERROR:  drawtree - division by zero in angleBetVectors()!\n");
    printf("Xu %f Yu %f Xv %f Yv %f\n", Xu, Yu, Xv, Yv);
    exxit(0);
  }
  cosTheta = dotProd / lengthsProd;

  if (cosTheta > 1) /* cosTheta will only be > 1 or < -1 due to rounding errors */
    theta = 0; /* cosTheta = 1 */
  else if (cosTheta < -1)
    theta = pi; /* cosTheta = -1 */
  else
    theta = acos(cosTheta);

  return theta;
}  /* angleBetVectors */


double signOfMoment(double xReferenceVector, double yReferenceVector,
                    double xForce, double yForce)
{
  /* it return the sign of the moment caused by the force, applied
     to the tip of the refereceVector; the root of the refereceVector
     is the pivot. */
  double angleReference, angleForce, sign;

  angleReference = computeAngle((double)0, (double)0, xReferenceVector,
                                yReferenceVector);
  angleForce = computeAngle((double)0, (double)0, xForce, yForce);
  angleForce = capedAngle(angleForce);
  angleReference = capedAngle(angleReference);

  /* reduce angleReference to 0 */
  angleForce = angleForce - angleReference;
  angleForce = capedAngle(angleForce);

  if (angleForce > 0 && angleForce < pi)
  {
    /* positive sign - force pointing toward the left of the reference
       line/vector.  */
    sign = 1;
  }
  else
  {
    /* negative sign */
    sign = -1;
  }
  return sign;
}  /* signOfMoment */


double vCounterClkwiseU(double Xu, double Yu, double Xv, double Yv)
{
/* Return 1 if vector v is counter clockwise from u */
  /* signOfMoment() is doing just that! */
  return signOfMoment(Xu, Yu, Xv, Yv);
}  /* vCounterClkwiseU */


double forcePerpendicularOnNode(node *pPivotSubNode, node *pToSubNode,
                                double medianDistance)
{
/* Read comment for totalForceOnNode */
/* It supposed to return a positive value to indicate that it has a positive
   moment; and negative return value to indicate negative moment.
   force perpendicular at norminal distance 1 is taken to be 1.
   medianDistance is the median of Distances in this graph.  */
/*
                        / Force
                       /
          |  ToNode   o  > alpha
          |              \
   yDelta |               \     theta = pi/2 + alpha
          |                  beta = vector (or angle) from Pivot to ToNode
  Pivot   o-----------
             xDelta

  alpha = theta + beta
 */
  double totalForce, forceAngle, xDelta, yDelta;
  double alpha, theta, forcePerpendicular, sinForceAngle, cosForceAngle;

  totalForce = (double)0;
  forceAngle = (double)0;

  totalForceOnNode(pPivotSubNode, pToSubNode, &totalForce, &forceAngle,
                   medianDistance);

  xDelta = curtree->nodep[pToSubNode->index-1]->xcoord -
    curtree->nodep[pPivotSubNode->index-1]->xcoord;
  yDelta = curtree->nodep[pToSubNode->index-1]->ycoord -
    curtree->nodep[pPivotSubNode->index-1]->ycoord;

  /* Try to avoid the case where 2 nodes are on top of each other. */
  /*
  if (xDelta < 0) tempx = -xDelta;
  else tempx = xDelta;
  if (yDelta < 0) tempy = -yDelta;
  else tempy = yDelta;
  if (tempx < epsilon && tempy < epsilon)
  {
    return;
  }
  */

  sinForceAngle = sin(forceAngle);
  cosForceAngle = cos(forceAngle);
  theta = angleBetVectors(xDelta, yDelta, cosForceAngle, sinForceAngle);

  if (theta > pi/2)
  {
    alpha = theta - pi/2;
  }
  else
  {
    alpha = pi/2 - theta;
  }

  forcePerpendicular = totalForce * cos(alpha);

  if (forcePerpendicular < -epsilon)
  {
    printf("ERROR:  drawtree - forcePerpendicular applied at an angle should"
           " not be less than zero (in forcePerpendicularOnNode()). \n");
    printf("alpha = %f\n", alpha);
    exxit(1);
  }
  /* correct the sign of the moment */
  forcePerpendicular = signOfMoment(xDelta, yDelta, cosForceAngle,
                                    sinForceAngle)
    * forcePerpendicular;
  return forcePerpendicular;
}  /* forcePerpendicularOnNode */


void polarizeABranch(node *pStartingSubNode, double *xx, double *yy)
{
  /* added - danieyek 990128 */
  /* After calling tilttrav(), if you don't polarize all the nodes on the branch
     to convert the x-y coordinates to theta and radius, you won't get result on
     the plot!  This function takes a subnode and branch out of all other subnode
     except the starting subnode (where the parent is), thus converting the x-y
     to polar coordinates for the branch only.  xx and yy are purely "inherited"
     features of polarize().  They should have been passed as values not
     addresses. */
  node *pSubNode;

  pSubNode = pStartingSubNode;

  /* convert the current node (note: not subnode) to polar coordinates. */
  polarize( curtree->nodep[pStartingSubNode->index - 1], xx, yy);

  /* visit the rest of the branches of current node */
  while (pSubNode->next != NULL && pSubNode->next != pStartingSubNode)
  {
    pSubNode = pSubNode->next;
    if ( pSubNode->tip != true )
      polarizeABranch(pSubNode->back, xx, yy);
  }
  return;
}  /* polarizeABranch */


void pushNodeToStack(stackElemType **ppStackTop, node *pNode)
{
  /* added - danieyek 990204 */
  /* pStackTop must be the current top element of the stack, where we add another
     element on top of it.
     ppStackTop must be the location where we can find pStackTop.
     This function "returns" the revised top (element) of the stack through
     the output parameter, ppStackTop.
     The last element on the stack has the "back" (pStackElemBack) pointer
     set to NULL.  So, when the last element is poped, ppStackTop will be
     automatically set to NULL.  If popNodeFromStack() is called with
     ppStackTop = NULL, we assume that it is the error caused by over popping
     the stack.
  */
  stackElemType *pStackElem;

  if (ppStackTop == NULL)
  {
    /* NULL can be stored in the location, but the location itself can't
       be NULL! */
    printf("ERROR:  drawtree - error using pushNodeToStack(); ppStackTop is NULL.\n");
    exxit(1);
  }
  pStackElem = (stackElemType*)Malloc(sizeof(stackElemType) );
  pStackElem->pStackElemBack = *ppStackTop;
  /* push an element onto the stack */
  pStackElem->pNode = pNode;
  *ppStackTop = pStackElem;
  return;
}  /* pushNodeToStack */


void popNodeFromStack(stackElemType **ppStackTop, node **ppNode)
{
  /* added - danieyek 990205 */
  /* pStackTop must be the current top element of the stack, where we pop an
     element from the top of it.
     ppStackTop must be the location where we can find pStackTop.
     This function "returns" the revised top (element) of the stack through
     the output parameter, ppStackTop.
     The last element on the stack has the "back" (pStackElemBack) pointer
     set to NULL.  So, when the last element is poped, ppStackTop will be
     automatically set to NULL.  If popNodeFromStack() is called with
     ppStackTop = NULL, we assume that it is the error caused by over popping
     the stack.
  */
  stackElemType *pStackT;

  if (ppStackTop == NULL)
  {
    printf("ERROR:  drawtree - a call to pop while the stack is empty.\n");
    exxit(1);
  }

  pStackT = *ppStackTop;
  *ppStackTop = pStackT->pStackElemBack;
  *ppNode  = pStackT->pNode;
  free(pStackT);

  return;
}  /* popNodeFromStack */


double medianOfDistance(node *pRootSubNode, boolean firstRecursiveCallP)
{
  /* added - danieyek 990208 */
  /* Find the median of the distance; used to compute the angle to rotate
     in proportion to the size of the graph and forces.
     It is assumed that pRootSubNode is also the pivot (during the first
     call to this function) - the center, with
     respect to which node the distances are calculated.
     If there are only 3 element, element #2 is returned, ie. (2+1)/2.
     This function now finds the median of distances of all nodes, not only
     the leafs!
  */
  node *pSubNode;
  double xDelta, yDelta, distance;
  long i, j;
  struct dblLinkNode
  {
    double value;
    /* Implement reverse Linked List */
    struct dblLinkNode *pBack;
  } *pLink, *pBackElem, *pMidElem, *pFrontElem, junkLink;
  /* must use static to retain values over calls */
  static node *pReferenceNode;
  static long count;
  static struct dblLinkNode *pFrontOfLinkedList;

  /* Remember the reference node so that it doesn't have to be passed
     arround in the function parameter. */
  if (firstRecursiveCallP == true)
  {
    pReferenceNode = pRootSubNode;
    pFrontOfLinkedList = NULL;
    count = 0;
  }

  pSubNode = pRootSubNode;
  /* visit the rest of the branches of current node; the branch attaches to
     the current subNode may be visited in the code further down below. */
  while (pSubNode->next != NULL && pSubNode->next != pRootSubNode)
  {
    pSubNode = pSubNode->next;
    if ( pSubNode->back != NULL)
      medianOfDistance(pSubNode->back, false);
  }

  /* visit this branch; You need to visit it for the first time - at root
     only! use pRootSubNode instead of pSubNode here because pSubNode stop
     short just before pRootSubNode (the entry node) */
  if ( firstRecursiveCallP == true && pRootSubNode->back != NULL)
    medianOfDistance(pRootSubNode->back, false);

  /* Why only leafs count?  Modifying it!  */
  xDelta = curtree->nodep[pSubNode->index-1]->xcoord -
    curtree->nodep[pReferenceNode->index-1]->xcoord;
  yDelta = curtree->nodep[pSubNode->index-1]->ycoord -
    curtree->nodep[pReferenceNode->index-1]->ycoord;
  distance = sqrt( xDelta*xDelta + yDelta*yDelta );

  /* Similar to pushing onto the stack */
  pLink = (struct dblLinkNode*) Malloc(sizeof(struct dblLinkNode) );
  if (pLink == NULL)
  {
    printf("Fatal ERROR:  drawtree - Insufficient Memory in medianOfDistance()!\n");
    exxit(1);
  }
  pLink->value = distance;
  pLink->pBack = pFrontOfLinkedList;
  pFrontOfLinkedList = pLink;
  count = count + 1;

  if (firstRecursiveCallP == true)
  {
    if (count == 0)
    {
      return (double)0;
    }
    else if (count == 1)
    {
      distance = pFrontOfLinkedList->value;
      free(pFrontOfLinkedList);
      return distance;
    }
    else if (count == 2)
    {
      distance = (pFrontOfLinkedList->value +
                  pFrontOfLinkedList->pBack->value)/(double)2;
      free(pFrontOfLinkedList->pBack);
      free(pFrontOfLinkedList);
      return distance;
    }
    else
    {
      junkLink.pBack = pFrontOfLinkedList;

      /* SORT first - use bubble sort; we start with at least 3 elements here. */
      /* We are matching backward when sorting the list and comparing MidElem and
         BackElem along the path;  junkLink is there just to make a symmetric
         operation at the front end. */
      for (j = 0; j < count - 1; j++)
      {
        pFrontElem = &junkLink;
        pMidElem = junkLink.pBack;
        pBackElem = junkLink.pBack->pBack;
        for (i = j; i < count - 1; i++)
        {
          if(pMidElem->value < pBackElem->value)
          {
            /* Swap - carry the smaller value to the root of the linked list. */
            pMidElem->pBack = pBackElem->pBack;
            pBackElem->pBack = pMidElem;
            pFrontElem->pBack = pBackElem;
            /* Correct the order of pFrontElem, pMidElem, pBackElem and match
               one step */
            pFrontElem = pBackElem;
            pBackElem = pMidElem->pBack;
          }
          else
          {
            pFrontElem = pMidElem;
            pMidElem = pBackElem;
            pBackElem = pBackElem->pBack;
          }
        }
        pFrontOfLinkedList = junkLink.pBack;
      }
      /* Sorted; now get the middle element. */
      for (i = 1; i < (count + 1)/(long) 2; i++)
      {
        /* Similar to Poping the stack */
        pLink = pFrontOfLinkedList;
        pFrontOfLinkedList = pLink->pBack;
        free(pLink);
      }

      /* Get the return value!! - only the last return value is the valid one. */
      distance = pFrontOfLinkedList->value;

      /* Continue from the same i value left off by the previous for loop!  */
      for (; i <= count; i++)
      {
        /* Similar to Poping the stack */
        pLink = pFrontOfLinkedList;
        pFrontOfLinkedList = pLink->pBack;
        free(pLink);
      }
    }
  }
  return distance;
}  /* medianOfDistance */


void leftRightLimits(node *pToSubNode, double *pLeftLimit,
                     double *pRightLimit)
/* As usual, pToSubNode->back is the angle
   leftLimit is the max angle you can rotate on the left and
   rightLimit vice versa.
   *pLeftLimit and *pRightLimit must be initialized to 0; without
   initialization, it would introduce bitter bugs into the program;
   they are initialized in this routine.
*/
{
  /* pPivotNode is nodep[pToSubNode->back->index-1], not pPivotSubNode
     which is just pToSubNode->back! */
  node *pLeftSubNode, *pRightSubNode, *pPivotNode, *pSubNode;
  double leftLimit, rightLimit, xToNodeVector, yToNodeVector, xLeftVector,
    yLeftVector, xRightVector, yRightVector, lengthsProd;

  *pLeftLimit = 0;
  *pRightLimit = 0;
  /* Make an assumption first - guess "pToSubNode->back->next" is the right
     and the opposite direction is the left! */
  /* It shouldn't be pivoted at a left, but just checking. */
  if (pToSubNode->back->tip == true)
  {
    /* Logically this should not happen.  But we actually can return pi
       as the limit. */
    printf("ERROR:  In leftRightLimits() - Pivoted at a leaf! Unable to "
           "calculate left and right limit.\n");
    exxit(1);
  }
  else if (pToSubNode->back->next->next == pToSubNode->back)
  {
    *pLeftLimit = 0; /* don't pivot where there is no branching */
    *pRightLimit = 0;
    return;
  }
  /* Else, do this */
  pPivotNode = curtree->nodep[pToSubNode->back->index-1];

  /* 3 or more branches - the regular case. */
  /* First, initialize the pRightSubNode - non-repeative portion of the code */
  pRightSubNode = pToSubNode->back;
  pLeftSubNode = pToSubNode->back;
  xToNodeVector = curtree->nodep[pToSubNode->index-1]->xcoord - pPivotNode->xcoord;
  yToNodeVector = curtree->nodep[pToSubNode->index-1]->ycoord - pPivotNode->ycoord;

  /* If both x and y are 0, then the length must be 0; but this check is not
     enough yet, we need to check the product of length also. */
  if ( fabs(xToNodeVector) < epsilon && fabs(yToNodeVector) < epsilon )
  {
    /* If the branch to rotate is too short, don't rotate it. */
    *pLeftLimit = 0;
    *pRightLimit = 0;
    return;
  }

  while( curtree->nodep[pRightSubNode->index-1]->tip != true )
  {
    /* Repeative code */
    pRightSubNode = pRightSubNode->next->back;

    xRightVector = curtree->nodep[pRightSubNode->index-1]->xcoord - pPivotNode->xcoord;
    yRightVector = curtree->nodep[pRightSubNode->index-1]->ycoord - pPivotNode->ycoord;

    lengthsProd = sqrt(xToNodeVector*xToNodeVector+yToNodeVector*yToNodeVector)
      * sqrt(xRightVector*xRightVector+yRightVector*yRightVector);
    if ( lengthsProd < epsilon )
    {
      continue;
    }
    rightLimit = angleBetVectors(xToNodeVector, yToNodeVector, xRightVector,
                                 yRightVector);

    if ( (*pRightLimit) < rightLimit) *pRightLimit = rightLimit;
  }

  while( curtree->nodep[pLeftSubNode->index-1]->tip != true )
  {
    /* First, let pSubNode be 1 subnode after rightSubNode. */
    pSubNode = pLeftSubNode->next->next;
    /* Then, loop until the last subNode before getting back to the pivot */
    while (pSubNode->next != pLeftSubNode)
    {
      pSubNode = pSubNode->next;
    }
    pLeftSubNode = pSubNode->back;

    xLeftVector = curtree->nodep[pLeftSubNode->index-1]->xcoord - pPivotNode->xcoord;
    yLeftVector = curtree->nodep[pLeftSubNode->index-1]->ycoord - pPivotNode->ycoord;

    lengthsProd = sqrt(xToNodeVector*xToNodeVector+yToNodeVector*yToNodeVector)
      * sqrt(xLeftVector*xLeftVector+yLeftVector*yLeftVector);
    if ( lengthsProd < epsilon )
    {
      continue;
    }
    leftLimit = angleBetVectors(xToNodeVector, yToNodeVector, xLeftVector,
                                yLeftVector);

    if ( (*pLeftLimit) < leftLimit) *pLeftLimit = leftLimit;
  }
  return;
}  /* leftRightLimits */


void branchLRHelper(node *pPivotSubNode, node *pCurSubNode,
                    double *pBranchL, double *pBranchR)
{
  /* added - danieyek 990226 */
  /* Recursive helper function for branchLeftRightAngles().
     pPivotSubNode->back is the pToNode, to which node you apply the forces!
  */
  /* Abandoned as it is similar to day-light algorithm; the first part is
     done implementing but not tested, the second part yet to be
     implemented if necessary. */
  double xCurNodeVector, yCurNodeVector, xPivotVector, yPivotVector;

  /* Base case : a leaf - return 0 & 0.  */
  if ( curtree->nodep[pCurSubNode->index-1]->tip == true )
  {
    xPivotVector = curtree->nodep[pPivotSubNode->back->index-1]->xcoord
      - curtree->nodep[pPivotSubNode->index-1]->xcoord;
    yPivotVector = curtree->nodep[pPivotSubNode->back->index-1]->ycoord
      - curtree->nodep[pPivotSubNode->index-1]->ycoord;
    xCurNodeVector = curtree->nodep[pCurSubNode->index-1]->xcoord
      - curtree->nodep[pPivotSubNode->index-1]->xcoord;
    yCurNodeVector = curtree->nodep[pCurSubNode->index-1]->ycoord
      - curtree->nodep[pPivotSubNode->index-1]->ycoord;

    if ( vCounterClkwiseU(xPivotVector, yPivotVector, xCurNodeVector,
                          yCurNodeVector)
         == 1)
    {
      /* Relevant to Left Angle */
      *pBranchL = angleBetVectors(xPivotVector, yPivotVector,
                                  xCurNodeVector, yCurNodeVector);
      *pBranchR = (double)0;
    }
    else
    {
      /* Relevant to Right Angle */
      *pBranchR = angleBetVectors(xPivotVector, yPivotVector,
                                  xCurNodeVector, yCurNodeVector);
      *pBranchL = (double)0;
    }
    return;
  }
  else
  {
    /* not a leaf */
  }
} /* branchLRHelper */


void improveNodeAngle(node *pToNode, double medianDistance)
{
  double forcePerpendicular, distance,
    xDistance, yDistance, angleRotate, sinAngleRotate, cosAngleRotate,
    norminalDistance, leftLimit, rightLimit, limitFactor;
  node *pPivot;

  /* Limit factor determinte how close the rotation can approach the
     absolute limit before colliding with other branches */
  limitFactor = (double)4 / (double)5;
  pPivot = pToNode->back;

  xDistance = curtree->nodep[pPivot->index-1]->xcoord - curtree->nodep[pToNode->index-1]->xcoord;
  yDistance = curtree->nodep[pPivot->index-1]->ycoord - curtree->nodep[pToNode->index-1]->ycoord;
  distance = sqrt( xDistance*xDistance + yDistance*yDistance  );

  /* convert distance to absolute value and test if it is zero */
  if ( fabs(distance) < epsilon)
  {
    angleRotate = (double)0;
  }
  else
  {
    leftRightLimits(pToNode, &leftLimit, &rightLimit);
    norminalDistance = distance / medianDistance;
    forcePerpendicular = forcePerpendicularOnNode(pPivot, pToNode,
                                                  medianDistance);
    angleRotate = forcePerpendicular / norminalDistance;
    /* Limiting the angle of rotation */
    if ( angleRotate > 0 && angleRotate > limitFactor * leftLimit)
    {
      /* Left */
      angleRotate = limitFactor * leftLimit;
    }
    else if ( -angleRotate > limitFactor * rightLimit )
      /* angleRotate < 0 && */
    {
      /* Right */
      angleRotate = - limitFactor * rightLimit;
    }
  }
  angleRotate = (double).1 * angleRotate;

  sinAngleRotate = sin(angleRotate);
  cosAngleRotate = cos(angleRotate);

  tilttrav(pToNode,
           &(curtree->nodep[pPivot->index - 1]->xcoord),
           &(curtree->nodep[pPivot->index - 1]->ycoord),
           &sinAngleRotate, &cosAngleRotate);

  polarizeABranch(pToNode,
                  &(curtree->nodep[pPivot->index - 1]->xcoord),
                  &(curtree->nodep[pPivot->index - 1]->ycoord));
}  /* improveNodeAngle */


void improvtravn(node *pStartingSubNode)
{
  /* improvtrav for n-body. */
  /* POPStack is the stack that is currently being used (popped);
     PUSHStack is the stack that is for the use of the next round
     (is pushed now) */
  stackElemType *pPUSHStackTop, *pPOPStackTop, *pTempStack;
  node *pSubNode, *pBackStartNode, *pBackSubNode;
  double medianDistance;
  long noOfIteration;

  /* Stack starts with no element on it */
  pPUSHStackTop = NULL;
  pPOPStackTop = NULL;


  /* Get the median to relate force to angle proportionally. */
  medianDistance = medianOfDistance(curtree->root, true);

  /* Set max. number of iteration */
  for ( noOfIteration = (long)0; noOfIteration < maxNumOfIter; noOfIteration++) {

    /* First, push all subNodes in the root node onto the stack-to-be-used
       to kick up the process */
    pSubNode = pStartingSubNode;
    pushNodeToStack(&pPUSHStackTop, pSubNode);
    while(pSubNode->next != pStartingSubNode)
    {
      pSubNode = pSubNode->next;
      pushNodeToStack(&pPUSHStackTop, pSubNode);
    }

    while (true)
    {
      /* Finishes with the current POPStack; swap the function of the stacks if
         PUSHStack is not empty */
      if (pPUSHStackTop == NULL)
      {
        /* Exit infinity loop here if empty. */
        break;
      }
      else
      {
        /* swap */
        pTempStack = pPUSHStackTop;
        pPUSHStackTop = pPOPStackTop;
        pPOPStackTop = pTempStack;
      }

      while (pPOPStackTop != NULL)
      {
        /* We always push the pivot subNode onto the stack!  That's when we
           pop that pivot subNode, subNode.back is the node we apply the force
           to (ToNode). Also, when we pop a pivot subNode, always push all
           pivot subNodes in the same ToNode onto the stack. */
        popNodeFromStack(&pPOPStackTop, &pSubNode);

        pBackStartNode = pSubNode->back;
        if (pBackStartNode->tip == true)
        {
          /* tip indicates if a node is a leaf */
          improveNodeAngle(pSubNode->back, medianDistance);

        }
        else
        {
          /* Push all subNodes in this pSubNode->back onto the stack-to-be-used,
           * stack-to-be-used, after poping a pivot subNode. If
           * pSubNode->back is a leaf, no push on stack. */
          pBackSubNode = pBackStartNode;
          /* Do not push this pBackStartNode onto the stack!  Or the
           * process will never stop. */
          while(pBackSubNode->next != pBackStartNode)
          {
            pBackSubNode = pBackSubNode->next;
            pushNodeToStack(&pPOPStackTop, pBackSubNode);
          }
          /* improve the node even if it is not a leaf */
          improveNodeAngle(pSubNode->back, medianDistance);
        }
      }
    }
  }
} /* improvtravn */


void coordimprov(double *xx, double *yy)
{
  /* use angles calculation to improve node coordinate placement */
  long i;

  (void)xx;                             // RSGnote: Parameter never used.
  (void)yy;                             // RSGnote: Parameter never used.

  if (nbody)       /* n-body algorithm */
    improvtravn(curtree->root);
  else {          /* equal-daylight algorithm */
    i = 0;
    do {
      maxchange = 0.0;
      improvtrav(curtree->root);
      i++;
    } while ((i < MAXITERATIONS) && (maxchange > MINIMUMCHANGE));
  }
}  /* coordimprov */


void calculate(void)
{
  /* compute coordinates for tree */
  double xx, yy;
  long i;
  double nttot, fontheight, labangle=0, top, bot, rig, lef;
  //printf("in calculate\n"); //JRMDebug
  //dumpnodelinks(root, nodep, nextnode); //JRMDebug

  //seetree(root, nodep, nextnode); //JRMDebug

  for (i = 0; i < nextnode; i++)
    ((drawtree_node*)curtree->nodep[i])->width = 1.0;

  //printf("after nodep->width\n"); //JRMDebug
  //dumpnodelinks(root, nodep, nextnode); //JRMDebug

  for (i = 0; i < nextnode; i++)
    curtree->nodep[i]->xcoord = 0.0;

  //printf("after nodep->xcoord\n"); //JRMDebug
  //dumpnodelinks(root, nodep, nextnode); //JRMDebug

  for (i = 0; i < nextnode; i++)
    curtree->nodep[i]->ycoord = 0.0;

  //printf("after nodep->ycoord\n"); //JRMDebug
  //dumpnodelinks(root, nodep, nextnode); //JRMDebug

  if (!uselengths) {
    for (i = 0; i < nextnode; i++)
      curtree->nodep[i]->length = 1.0;
  } else {
    for (i = 0; i < nextnode; i++)
      curtree->nodep[i]->length = fabs(curtree->nodep[i]->oldlen);
  }

  /* //JRMDebug
  for (i=0; i<nextnode; i++)
  {
    printf("i: %li width: %f xcoord: %f ycoord: %f length: %f\n",
           i, ((drawtree_node*)nodep[i])->width, nodep[i]->xcoord, nodep[i]->ycoord, nodep[i]->length);
  }
  */
  //printf("before getwidth\n"); //JRMDebug
  //dumpnodelinks(root, nodep, nextnode); //JRMDebug

  //printf("calling getwidth root: %p root->type: %i nayme: %s index: %li\n", root, root->type, root->nayme, root->index);  //JRMDebug

  getwidth(curtree->root);
  nttot = ((drawtree_node*)curtree->root)->width;
  //printf("nttot: %f\n", nttot); //JRMDebug
  for (i = 0; i < nextnode; i++)
    ((drawtree_node*)curtree->nodep[i])->width = ((drawtree_node*)curtree->nodep[i])->width * spp / nttot;
  if (!improve)
    plrtrans(curtree->root, treeangle, treeangle - ark / 2.0, treeangle + ark / 2.0);
  else plrtrans(curtree->root, treeangle, treeangle - pi, treeangle + pi);
  maxx = 0.0;
  minx = 0.0;
  maxy = 0.0;
  miny = 0.0;
  //printf("calling coordtrav\n"); //JRMDebug
  coordtrav(curtree->root, &xx, &yy);
  fontheight = heighttext(font, fontname);
  if (labeldirec == fixed)
    labangle = pi * labelrotation / 180.0;
  textlength = (double*) Malloc(nextnode * sizeof(double));
  firstlet = (double*) Malloc(nextnode * sizeof(double));
  for (i = 0; i < nextnode; i++) {
    if (curtree->nodep[i]->tip) {
      textlength[i] = lengthtext(curtree->nodep[i]->nayme,
                                curtree->nodep[i]->naymlength, fontname, font);
      textlength[i] /= fontheight;
      firstlet[i] = lengthtext(curtree->nodep[i]->nayme, 1L, fontname, font)
        / fontheight;
    }
  }
  if (spp > 1)
    labelheight = charht * (maxx - minx) / (spp - 1);
  else
    labelheight = charht * (maxx - minx);
  if (improve) {
    coordimprov(&xx, &yy);
    maxx = 0.0;
    minx = 0.0;
    maxy = 0.0;
    miny = 0.0;
    coordtrav(curtree->root, &xx, &yy);
  }
  topoflabels = 0.0;
  bottomoflabels = 0.0;
  rightoflabels = 0.0;
  leftoflabels = 0.0;
  for (i = 0; i < nextnode; i++) {
    if (curtree->nodep[i]->tip) {
      if (labeldirec == radial)
        labangle = ((draw_node*)curtree->nodep[i])->theta;
      else if (labeldirec == along)
        labangle = ((draw_node*)curtree->nodep[i])->oldtheta;
      else if (labeldirec == middle)
        labangle = 0.0;
      if (cos(labangle) < 0.0 && labeldirec != fixed)
        labangle -= pi;
      firstlet[i] = lengthtext(curtree->nodep[i]->nayme, 1L, fontname, font)
        / fontheight;
      top = (curtree->nodep[i]->ycoord - maxy) / labelheight +
        sin(((draw_node*)curtree->nodep[i])->oldtheta);
      rig = (curtree->nodep[i]->xcoord - maxx) / labelheight +
        cos(((draw_node*)curtree->nodep[i])->oldtheta);
      bot = (miny - curtree->nodep[i]->ycoord) / labelheight -
        sin(((draw_node*)curtree->nodep[i])->oldtheta);
      lef = (minx - curtree->nodep[i]->xcoord) / labelheight -
        cos(((draw_node*)curtree->nodep[i])->oldtheta);
      if (cos(labangle) * cos(((draw_node*)curtree->nodep[i])->oldtheta) +
          sin(labangle) * sin(((draw_node*)curtree->nodep[i])->oldtheta) > 0.0) {
        if (sin(labangle) > 0.0)
          top += sin(labangle) * textlength[i];
        top += sin(labangle - 1.25 * pi) * GAP * firstlet[i];
        if (sin(labangle) < 0.0)
          bot -= sin(labangle) * textlength[i];
        bot -= sin(labangle - 0.75 * pi) * GAP * firstlet[i];
        if (sin(labangle) > 0.0)
          rig += cos(labangle - 0.75 * pi) * GAP * firstlet[i];
        else
          rig += cos(labangle - 1.25 * pi) * GAP * firstlet[i];
        rig += cos(labangle) * textlength[i];
        if (sin(labangle) > 0.0)
          lef -= cos(labangle - 1.25 * pi) * GAP * firstlet[i];
        else
          lef -= cos(labangle - 0.75 * pi) * GAP * firstlet[i];
      } else {
        if (sin(labangle) < 0.0)
          top -= sin(labangle) * textlength[i];
        top += sin(labangle + 0.25 * pi) * GAP * firstlet[i];
        if (sin(labangle) > 0.0)
          bot += sin(labangle) * textlength[i];
        bot -= sin(labangle - 0.25 * pi) * GAP * firstlet[i];
        if (sin(labangle) > 0.0)
          rig += cos(labangle - 0.25 * pi) * GAP * firstlet[i];
        else
          rig += cos(labangle + 0.25 * pi) * GAP * firstlet[i];
        if (sin(labangle) < 0.0)
          rig += cos(labangle) * textlength[i];
        if (sin(labangle) > 0.0)
          lef -= cos(labangle + 0.25 * pi) * GAP * firstlet[i];
        else
          lef -= cos(labangle - 0.25 * pi) * GAP * firstlet[i];
        lef += cos(labangle) * textlength[i];
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
  }
  topoflabels *= labelheight;
  bottomoflabels *= labelheight;
  leftoflabels *= labelheight;
  rightoflabels *= labelheight;
}  /* calculate */


void rescale(void)
{
  /* compute coordinates of tree for plot or preview device */
  long i;
  double treeheight, treewidth, extrax, extray, temp;

  treeheight = maxy - miny + topoflabels + bottomoflabels;
  treewidth = maxx - minx + rightoflabels + leftoflabels;
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
  for (i = 0; i < (nextnode); i++) {
    curtree->nodep[i]->xcoord = expand * (curtree->nodep[i]->xcoord - minx + leftoflabels);
    curtree->nodep[i]->ycoord = expand * (curtree->nodep[i]->ycoord - miny + bottomoflabels);
    if (grows == horizontal) {
      temp = curtree->nodep[i]->ycoord;
      curtree->nodep[i]->ycoord = expand * treewidth - curtree->nodep[i]->xcoord;
      curtree->nodep[i]->xcoord = temp;
    }
    curtree->nodep[i]->xcoord += xmargin + extrax;
    curtree->nodep[i]->ycoord += ymargin + extray;
  }
}  /* rescale */


void plottree(node *p, node *q)
{
  /* plot part or all of tree on the plotting device */
  double x1, y1, x2, y2;
  node *pp;

  x2 = xscale * (xoffset + p->xcoord);
  y2 = yscale * (yoffset + p->ycoord);
  if (p != curtree->root) {
    x1 = xscale * (xoffset + q->xcoord);
    y1 = yscale * (yoffset + q->ycoord);
    plot(penup, x1, y1);
    plot(pendown, x2, y2);
  }
  if (p->tip)
    return;
  pp = p->next;
  do {
    plottree(pp->back, p);
    pp = pp->next;
  } while (((p == curtree->root) && (pp != p->next)) || ((p != curtree->root) && (pp != p)));
}  /* plottree */


void plotlabels(char *fontname)
{
  long i;
  double compr, dx = 0, dy = 0, labangle, sino, coso, cosl, sinl,
    cosv, sinv, vec;
  boolean right;
  node *lp;

  compr = xunitspercm / yunitspercm;
  if (penchange == yes)
    changepen(labelpen);
  for (i = 0; i < (nextnode); i++) {
    if (curtree->nodep[i]->tip) {
      lp = curtree->nodep[i];
      labangle = labelrotation * pi / 180.0;
      if (labeldirec == radial)
        labangle = ((draw_node*)curtree->nodep[i])->theta;
      else if (labeldirec == along)
        labangle = ((draw_node*)curtree->nodep[i])->oldtheta;
      else if (labeldirec == middle)
        labangle = 0.0;
      if (cos(labangle) < 0.0)
        labangle -= pi;
      sino = sin(((draw_node*)curtree->nodep[i])->oldtheta);
      coso = cos(((draw_node*)curtree->nodep[i])->oldtheta);
      cosl = cos(labangle);
      sinl = sin(labangle);
      right = ((coso*cosl+sino*sinl) > 0.0) || (labeldirec == middle);
      vec = sqrt(1.0+firstlet[i]*firstlet[i]);
      cosv = firstlet[i]/vec;
      sinv = 1.0/vec;
      if (labeldirec == middle) {
        if ((textlength[i]+1.0)*fabs(tan(((draw_node*)curtree->nodep[i])->oldtheta))
            > 2.0) {
          dx = -0.5 * textlength[i] * labelheight * expand;
          if (sino > 0.0) {
            dy = 0.5 * labelheight * expand;
            if (fabs(((draw_node*)curtree->nodep[i])->oldtheta - pi/2.0) > 1000.0)
              dx += labelheight * expand /
                (2.0*tan(((draw_node*)curtree->nodep[i])->oldtheta));
          } else {
            dy = -1.5 * labelheight * expand;
            if (fabs(((draw_node*)curtree->nodep[i])->oldtheta - pi/2.0) > 1000.0)
              dx += labelheight * expand /
                (2.0*tan(((draw_node*)curtree->nodep[i])->oldtheta));
          }
        }
        else {
          if (coso > 0.0) {
            dx = 0.5 * labelheight * expand;
            dy = (-0.5 + (0.5*textlength[i]+0.5)*tan(((draw_node*)curtree->nodep[i])->oldtheta))
              * labelheight * expand;
          }
          else {
            dx = -(textlength[i]+0.5) * labelheight * expand;
            dy = (-0.5 - (0.5*textlength[i]+0.5)*tan(((draw_node*)curtree->nodep[i])->oldtheta))
              * labelheight * expand;
          }
        }
      } else {
        if (right) {
          dx = labelheight * expand * coso;
          dy = labelheight * expand * sino;
          dx += labelheight * expand * 0.5 * vec * (-cosl*cosv+sinl*sinv);
          dy += labelheight * expand * 0.5 * vec * (-sinl*cosv-cosl*sinv);
        } else {
          dx = labelheight * expand * coso;
          dy = labelheight * expand * sino;
          dx += labelheight * expand * 0.5 * vec * (cosl*cosv+sinl*sinv);
          dy += labelheight * expand * 0.5 * vec * (sinl*cosv-cosl*sinv);
          dx -= textlength[i] * labelheight * expand * cosl;
          dy -= textlength[i] * labelheight * expand * sinl;
        }
      }
      plottext(lp->nayme, lp->naymlength,
               labelheight * expand * xscale / compr, compr,
               xscale * (lp->xcoord + dx + xoffset),
               yscale * (lp->ycoord + dy + yoffset), -180 * labangle / pi,
               font, fontname);
    }
  }
  if (penchange == yes)
    changepen(treepen);
}  /* plotlabels */


void user_loop(void)
{
  /* loop to make preview window and decide what to do with it */

  long loopcount;
  char input_char;

  while (!canbeplotted) {
    loopcount = 0;
    do {
      input_char=showparms();
      firstscreens = false;
      if ( input_char != 'Y')
        getparms(input_char);
      countup(&loopcount, 10);
    } while (input_char != 'Y');

#ifdef JAVADEBUG
    //JRMDebug
    printf("***internal parameters***\n");
    printf("fontname: %s\n", fontname);
    printf("plotter: %i\n", plotter);
    printf("previewer: %i\n", previewer);
    printf("paperx: %f\n", paperx);
    printf("papery: %f\n", papery);
    printf("hpmargin: %f\n", hpmargin);
    printf("vpmargin: %f\n", vpmargin);
    printf("pagex: %f\n", pagex);
    printf("pagey: %f\n", pagey);
    printf("xmargin: %f\n", xmargin);
    printf("ymargin: %f\n", ymargin);
    printf("rescaled: %i\n", rescaled);
    printf("previewing: %i\n", previewing);
    printf("firstscreens: %i\n", firstscreens);
    printf("canbeplotted: %i\n", canbeplotted);
    printf("labelrotation: %f\n", labelrotation);
    printf("treeangle: %f\n", treeangle);
    printf("ark: %f\n", ark);
    printf("bscale: %f\n", bscale);
    printf("charht: %f\n", charht);
    printf("uselengths: %i\n", uselengths);
    printf("regular: %i\n", regular);
    printf("rescaled: %i\n", rescaled);
    printf("labeldirec: %i\n", labeldirec);
    printf("improve: %i\n", improve);
    printf("nbody: %i\n", nbody);
    printf("grows: %i\n", grows);
    printf("dotmatrix: %i\n", dotmatrix);
#endif

    xscale = xunitspercm;
    yscale = yunitspercm;
    plotrparms(spp);
    numlines = dotmatrix
      ? ((long)floor(yunitspercm * ysize + 0.5) / strpdeep):1;

    //printf("calling calculate\n"); //JRMDebug
    calculate();

    //printf("calling rescale\n"); //JRMDebug
    rescale();
    canbeplotted = true;
    if (!javarun)
    {
      if (preview)
        canbeplotted=plotpreview(fontname, &xoffset, &yoffset, &scale, spp, curtree->root);
      else
        canbeplotted = true;
      if ((previewer == winpreview || previewer == xpreview || previewer == mac)
          && (winaction == quitnow))
        canbeplotted = true;
    }
  }
} /* user_loop */


void setup_environment(int argc, Char *argv[])
{
  /* Set up all kinds of fun stuff */
  node *q, *r;

  char *pChar;
  double i;
  boolean firsttree;
#ifdef MAC
  OSErr retcode;
  FInfo  fndrinfo;
  macsetup("Drawtree", "Preview");
#endif

#ifdef TURBOC
  if ((registerbgidriver(EGAVGA_driver) <0) ||
      (registerbgidriver(Herc_driver) <0)   ||
      (registerbgidriver(CGA_driver) <0)) {
    fprintf(stderr, "Graphics error: %s ", grapherrormsg(graphresult()));
    exxit(-1);
  }
#endif

  if(javarun)
  {
    intree = fopen(trefilename, "r");
  }
  else
  {
    printf("DRAWTREE from PHYLIP version %s\n", VERSION);

    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree, INTREE, "input tree file", "rb", argv[0], NULL);
    printf("Reading tree ... \n");
  }
  firsttree = true;
  inputNumbersFromTreeFile(intree, &spp, &nonodes);
  curtree = functions.tree_new(nonodes, spp);
  treeread (curtree, intree, &curtree->root, curtree->nodep, &goteof, &firsttree, &nextnode, &haslengths, initdrawtreenode, true, -1);

  //printf("Raw Tree nextnode: %li root: %p\n", nextnode, root); //JRMDebug
  //dumpnodelinks(root, nodep, nextnode); //JRMDebug

  // if there are no lengths, shut uselengths off
  if (haslengths==0)
    uselengths = false;

  // the following changes the next pointer of the last entry in the tree
  // so it points to the root, which changes it from an open tree to a closed loop
  q = curtree->root;
  r = curtree->root;
  while (!(q->next == curtree->root))
    q = q->next;
  q->next = curtree->root->next;
  curtree->root = q;
  curtree->nodep[spp] = q;
  where = curtree->root;
  rotate = true;
  if(!javarun)
  {
    printf("Tree has been read.\n");
    //printf("Tree has been read. spp: %li\n", spp); //JRMDebug
    //printf("Tree after reroot root: %p  next: %p back: %p\n", root, root->next, root->back); //JRMDebug
    //dumpnodelinks(root, nodep, nextnode); //JRMDebug

    printf("Loading the font ... \n");
    loadfont(font, argv[0]);
    printf("Font loaded.\n");
  }
  previewing = false;
  ansi = ANSICRT;
  ibmpc = IBMCRT;
  firstscreens = true;
  initialparms();
  canbeplotted = false;
  if (argc > 1)
  {
    pChar = argv[1];
    for (i = 0; i < strlen(pChar); i++)
    {
      if ( ! isdigit(*pChar) )
      {
        /* set to default if the 2nd. parameter is not a number */
        maxNumOfIter = 50;
        return;
      }
      else if ( isspace(*pChar) )
      {
        printf("ERROR:  Number of iteration should not contain space!\n");
        exxit(1);
      }
    }
    sscanf(argv[1], "%li", &maxNumOfIter);
  }
  else
  {
    /* 2nd. argument is not entered; use default. */
    maxNumOfIter = 50;
  }
  return;
}  /* setup_environment */


void drawtree(
  char*  intreename,
  char*  plotfilename,
  char*  plotfileopt,
  char*  fontfilename,
  char*  treegrows,
  int    usebranchlengths,
  char*  labeldirecstr,
  double labelangle,
  double treerotation,
  double treearc,
  char*  iterationkind,
  int    iterationcount,
  int    regularizeangles,
  int    avoidlabeloverlap,
  int    branchrescale,
  double branchscaler,
  double relcharhgt,
  double  xmarginratio,
  double  ymarginratio,
  int    dofinalplot,
  char*  finalplotkind)
{
  (void)plotfileopt;                    // RSGnote: Parameter never used.
  //printf("hello from DrawTree\n");    //JRMDebug

  strcpy(trefilename, intreename); // need this for setup_environment

  initdata* funcs;
  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Drawtree";
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = drawtree_node_new;
  phylipinit(argc, argv, funcs, true);
  strcpy(progname, argv[0]);
  ibmpc     = IBMCRT;
  ansi      = ANSICRT;
  //setup_environment(argc, argv);

  char *previewfilename = "JavaPreview.ps";
  previewer = lw;  // hardwired to ps

  // translate Java doubles and ints to Drawtree local variables
  labelrotation = labelangle;
  treeangle     = treerotation * pi / 180;
  ark           = treearc * pi / 180;
  bscale        = branchscaler;
  charht        = relcharhgt;
  maxNumOfIter  = iterationcount;
  
  setup_environment(argc, argv);

  strcpy(fontname, fontfilename); 
  boolean wasplotted = false;

  // translate Java boolean to Phylip boolean
  boolean doplot;
  if (dofinalplot != 0)
    doplot = true;
  else
    doplot = false;

    if (!haslengths) {
        uselengths = false;
    }
    else {
        if (usebranchlengths != 0) {
            uselengths = true;
        }
        else {
            uselengths = false;
        }
    }

  if (regularizeangles != 0)
    regular = true;
  else
    regular = false;

  if (avoidlabeloverlap != 0)
    labelavoid = true;
  else
    labelavoid = false;

  if (branchrescale != 0)
    rescaled = true;
  else
    rescaled = false;

  // label direction
  labeldirec = middle;
  if (!strcmp(labeldirecstr, "fixed"))
    labeldirec = fixed;
  if (!strcmp(labeldirecstr, "middle"))
    labeldirec = middle;
  if (!strcmp(labeldirecstr, "radial"))
    labeldirec = radial;
  if (!strcmp(labeldirecstr, "along"))
    labeldirec = along;

  // iterate to improve tree
  improve = true;
  nbody = false;
  if (!strcmp(iterationkind, "improve"))
    improve = true;
  if (!strcmp(iterationkind, "nbody"))
    nbody = true;

  // tree growth direction
  grows = horizontal;
  if (!strcmp(treegrows, "vertical"))
    grows = vertical;

  // figure out plot type
  plotter = lw;
  if(doplot)
  {
    if (!strcmp(finalplotkind, "lw"))
      plotter = lw;
    if(!strcmp(finalplotkind, "svg"))
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

#ifdef JAVADEBUG
  // dump translated parameters from Java to make sure they made it
  //JRMDebug
  printf("***Java input parameters***\n");
  printf("intreename: %s\n", intreename);
  printf("plotfilename: %s\n", plotfilename);
  printf("plotfileopt: %s\n", plotfileopt);
  printf("fontfilename: %s\n", fontfilename);
  printf("usebranchlengths:: %i\n", usebranchlengths);
  printf("labeldirecstr: %s\n", labeldirecstr);
  printf("labelangle: %f\n", labelangle);
  printf("treerotation: %f\n", treerotation);
  printf("treearc: %f\n", treearc);
  printf("iterationkind: %s\n", iterationkind);
  printf("iterationcount: %i\n", iterationcount);
  printf("regularizeangles: %i\n", regularizeangles);
  printf("avoidlabeloverlap: %i\n", avoidlabeloverlap);
  printf("branchrescale: %i\n", branchrescale);
  printf("branchscaler: %f\n", branchscaler);
  printf("relcharhgt: %f\n", relcharhgt);
  printf("xmarginratio: %f\n", xmarginratio);
  printf("ymarginratio: %f\n", ymarginratio);
  printf("dofinalplot: %i\n", dofinalplot);
  printf("previewfilename: %s\n", previewfilename);
  printf("finalplotkind: %s\n", finalplotkind);
  printf("doplot: %i\n", doplot);
#endif

  // start drawtree code
  numlines = 1;
  xscale = xunitspercm;
  yscale = yunitspercm;
  xmargin = xmarginratio * paperx;
  ymargin = ymarginratio * papery;

#ifdef JAVADEBUG
  printf("\n***internal parameters***\n");
  printf("fontname: %s\n", fontname);
  printf("plotter: %i\n", plotter);
  printf("previewer: %i\n", previewer);
  printf("paperx: %f\n", paperx);
  printf("papery: %f\n", papery);
  printf("hpmargin: %f\n", hpmargin);
  printf("vpmargin: %f\n", vpmargin);
  printf("pagex: %f\n", pagex);
  printf("pagey: %f\n", pagey);
  printf("xmargin: %f\n", xmargin);
  printf("ymargin: %f\n", ymargin);
  printf("rescaled: %i\n", rescaled);
  printf("previewing: %i\n", previewing);
  printf("firstscreens: %i\n", firstscreens);
  printf("canbeplotted: %i\n", canbeplotted);
  printf("labelrotation: %f\n", labelrotation);
  printf("treeangle: %f\n", treeangle);
  printf("maxNumOfIter: %li\n", maxNumOfIter);
  printf("ark: %f\n", ark);
  printf("bscale: %f\n", bscale);
  printf("charht: %f\n", charht);
  printf("uselengths: %i\n", uselengths);
  printf("regular: %i\n", regular);
  printf("rescaled: %i\n", rescaled);
  printf("labeldirec: %i\n", labeldirec);
  printf("improve: %i\n", improve);
  printf("nbody: %i\n", nbody);
  printf("grows: %i\n", grows);
  printf("dotmatrix: %i\n", dotmatrix);
  fflush(stdout);
#endif

  calculate();
  rescale();

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
  plotfile = fopen(usefilename, "wb");

  //printf("calling initplotter plotter = %i\n", plotter); //JRMDebug
  previewing = false;
  initplotter(spp, fontname);
  numlines = dotmatrix ? ((long)floor(yunitspercm * ysize + 0.5)/strpdeep) : 1;

  if (!dofinalplot)
  {
      changepen(labelpen);
      makebox(fontname,&xoffset,&yoffset,&scale,spp);
      changepen(treepen);
  }

  //printf("\nWriting plot file ...\n"); //JRMDebug
  //printf("calling drawit\n"); //JRMDebug
  drawit(fontname, &xoffset, &yoffset, numlines, curtree->root);
  //printf("calling finishplotter\n");
  finishplotter();

  fclose(plotfile);
  wasplotted = true;
  //printf("\nPlot written to file \"%s\".\n\n", usefilename); //JRMDebug

  fclose(intree);

  //printf("Done.\n\n"); //JRMDebug
  return;
}


int main(int argc, Char *argv[])
{
  long stripedepth;
  boolean wasplotted = false;
  initdata* funcs;
#ifdef MAC
  char filename1[FNMLNGTH];
  OSErr retcode;
  FInfo  fndrinfo;
#ifdef OSX_CARBON
  FSRef fileRef;
  FSSpec fileSpec;
#endif
#ifdef __MWERKS__
  SIOUXSetTitle("\pPHYLIP:  Drawtree");
#endif
  argv[0] = "Drawtree";
#endif
#ifdef PHYLIP_CAN_DRAW_WITH_X
  nargc=argc;
  nargv=argv;
#endif
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = drawtree_node_new;
  phylipinit(argc, argv, funcs, false);

  strcpy(progname, argv[0]);
  setup_environment(argc, argv);

  user_loop();

  if (dotmatrix) {
    stripedepth = allocstripe(stripe, (strpwide/8), ((long)(yunitspercm * ysize)));
    strpdeep = stripedepth;
    strpdiv  = stripedepth;
  }

  if (!((previewer == winpreview || previewer == xpreview || previewer == mac)
        && (winaction == quitnow))) {
    /* Open plotfiles in binary mode */
    openfile(&plotfile, PLOTFILE, "plot file", "wb", argv[0], pltfilename);
    previewing = false;
    initplotter(spp, fontname);
    numlines = dotmatrix ? ((long)floor(yunitspercm * ysize + 0.5)/strpdeep) : 1;
    if (plotter != ibm)
      printf("\nWriting plot file ...\n");
    drawit(fontname, &xoffset, &yoffset, numlines, curtree->root);
    finishplotter();
    wasplotted = true;
    FClose(plotfile);
    printf("\nPlot written to file \"%s\".\n\n", pltfilename);
  }

  FClose(intree);
  printf("Done.\n\n");
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
    strcpy(filename1, pltfilename);
    retcode=GetFInfo(CtoPstr(filename1), 0, &fndrinfo);
    fndrinfo.fdType='PICT';
    fndrinfo.fdCreator='MDRW';
    strcpy(filename1, pltfilename);
    retcode=SetFInfo(CtoPstr(PLOTFILE), 0, &fndrinfo);
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
    retcode=GetFInfo(CtoPstr(PLOTFILE), 0, &fndrinfo);
    fndrinfo.fdType='TEXT';
    retcode=SetFInfo(CtoPstr(PLOTFILE), 0, &fndrinfo);
#endif
  }
#endif
  phyRestoreConsoleAttributes();
  exxit(0);
  return 1;
}


// End.
