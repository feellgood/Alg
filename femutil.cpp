#include "fem.h"
#include <set>
#include <limits>       // std::numeric_limits

const double HUGE = std::numeric_limits<double>::max();

void femutil(Fem &fem)
{
const int NOD = fem.NOD;
const int SEG = fem.SEG;
const int TRI = fem.TRI;
pair <string,int> p,p1,p2;
fem.pts= annAllocPts(NOD, 2);

// calcul du diametre et du centrage
double xmin, xmax, ymin, ymax, lx,ly;
xmin = ymin = +HUGE;
xmax = ymax = -HUGE;

for (int i=0; i<NOD; i++){
    double xi,yi;
    xi = fem.node[i].x;      yi = fem.node[i].y;
    fem.pts[i][0]=xi;        fem.pts[i][1]=yi;
    if (xi<xmin) xmin=xi;    
	if (xi>xmax) xmax=xi;
    if (yi<ymin) ymin=yi;    
	if (yi>ymax) ymax=yi;
    }

// allocation de l'arbre de recherche
fem.kdtree = new ANNkd_tree(fem.pts, NOD, 2);
if (!fem.kdtree) exit(1);

lx=xmax-xmin; fem.lx=lx;
ly=ymax-ymin; fem.ly=ly;

fem.diam = lx;
if (fem.diam<ly) fem.diam=ly;

fem.cx = 0.5*(xmax+xmin);
fem.cy = 0.5*(ymax+ymin);

fem.as[0] = lx/fem.diam;
fem.as[1] = ly/fem.diam;

for (int s=0; s<SEG; s++){
   Seg &seg = fem.seg[s];
   int i0,i1;
   i0=seg.ind[0];   i1=seg.ind[1];
   
   double x0,y0, x1,y1;
   x0 = fem.node[i0].x;   y0 = fem.node[i0].y;
   x1 = fem.node[i1].x;   y1 = fem.node[i1].y;
   seg.len=sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
   }

// calcul des surfaces et reorientation des segments si necessaire
double surftot = 0.;

for (int t=0; t<TRI; t++){
   Tri &tri = fem.tri[t];
   int i0,i1,i2;
   i0=tri.ind[0];   i1=tri.ind[1];   i2=tri.ind[2];  
   
   double x0,y0, x1,y1, x2,y2;
   x0 = fem.node[i0].x;   y0 = fem.node[i0].y;
   x1 = fem.node[i1].x;   y1 = fem.node[i1].y;
   x2 = fem.node[i2].x;   y2 = fem.node[i2].y;
 
   double surf;
   surf = 0.5*((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0));
   if (surf<0.) {
      tri.ind[2]=i1; tri.ind[1]=i2;
      surf=-surf;
//      cout << "ill-oriented triangle : " << t << " now corrected!"<< endl;
      }
  
   tri.surf = surf;
   surftot+= surf;
   }
fem.surf = surftot;
}
