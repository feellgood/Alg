#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <list>

#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#include <boost/format.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

#include "alg_core.h"

//#include "alg_sparseMat.h"
//#include "alg_denseMat.h"

const double EPSILON = 1e-12;
const double VACUUM_PERMITTIVITY = 8.854187817e-12;  /* A^2 s^4 / kg m^3 */
const double VACUUM_PERMEABILITY = 1.25663706144e-6; /* kg m / A^2 s^2   */

using namespace std;

struct Node {
    double x, y;
    double sol;                                             /* valeur nodale de la solution */
    };

struct Seg{
    static const int NBN = 2, NPI = 2;
    int reg;                                                /* numero de region */
    double len;                                             /* longueur du segment */
    int ind[NBN];                                           /* noeuds de l'element */
    double x[NPI], y[NPI];                                  /* points de gauss */
    double weight[NPI];                                     /* poids de gauss x detJ */
    double a[NBN][NPI];                                     /* polynomes d'interpolation */
    };

struct Tri{
    static const int NBN = 3, NPI = 3;
    int reg;                                                /* numero de region */
    double surf;                                            /* surface de l'element */
    int ind[NBN];                                           /* noeuds de l'element */
    double x[NPI], y[NPI];                                  /* points de gauss */
    double weight[NPI];                                     /* poids de gauss x detJ */
    double a[NBN][NPI];                                     /* polynomes d'interpolation */
    double dadx[NBN][NPI], dady[NBN][NPI], dadz[NBN][NPI];  /* derivees */
    };
 
struct Fem{
    string  simname, pbname;                                /* nom de la simulation et du fichier probleme */
    int REG, NOD, SEG, TRI;
    double scale, surf;
    double cx, cy, lx, ly, diam, as[2];
    vector <Node> node;                                     /* liste des noeuds du maillage */
    vector <Tri>  tri;                                      /* liste des triangles */
    vector <Seg>  seg;                                      /* liste des segments */
    map <pair<string,int>,double > param;                   /* tableaux associatifs des parametres physiques */
    map <pair<string,int>,double > dir;                     /* tableaux associatifs des conditions de dirichlet */
    map <pair<string,int>,string > per;                     /* tableaux associatifs des conditions periodiques */
    };

void extract_comment(istream &flux);
inline double sq(double x) {return x*x;}

void dialog(Fem& fem);
void lecture(Fem &fem);

void femutil(Fem &fem);
void chapeaux(Fem &fem);
void affichage(Fem &fem);

void plotsol(Fem &fem);
void savevtk(Fem &fem);
void savesol(Fem &fem);

void solve(Fem &fem,const int MAXITER);

void integrales(Fem &fem, Tri &tri, alg::denseMat &AE, vector <double> &BE);
void integrales(Fem &fem, Seg &seg, vector <double> &BE);

template <class T>
void assemblage(T &obj, alg::denseMat    &Ke, vector <double> &Le, alg::sparseMat &K,  vector <double> &L)
{
const size_t NBN = T::NBN;

for (size_t ie=0; ie<NBN; ie++){
    size_t i= obj.ind[ie];             
    for (size_t je=0; je<NBN; je++){
        size_t j= obj.ind[je];
        K.push_back(i, j, Ke(ie, je));  
        }
    L[i]+= Le[ie];
    }
}

