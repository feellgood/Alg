#include <chrono>

#include "fem.h"
#include "alg.h"

inline std::ostream & operator<<(std::ostream & flux, std::vector<size_t> const& v) {
std::for_each(v.begin(),v.end(), [&flux](const size_t& x) { flux << x << " "; });
return flux;
}

inline std::ostream & operator<<(std::ostream & flux, std::vector<double> const& v) {
std::for_each(v.begin(),v.end(), [&flux](const double& x) { flux << x << " "; });
return flux;
}

inline std::ostream & operator<<(std::ostream & flux, alg::sparseMat const& m) {
m.print(flux);
return flux;
}

void solve(Fem &fem, const int MAXITER)
{
typedef map <pair<string,int>,double> mapType;

const size_t NOD = fem.NOD;
const int TRI = fem.TRI;
const int SEG = fem.SEG;

alg::sparseMat K(NOD);
std::vector<double> Lw(NOD, 0.0);
std::vector<double> Xw(NOD);

auto t1 = std::chrono::high_resolution_clock::now();

for (int t=0; t<TRI; t++){
    Tri &tri = fem.tri[t];
    const int NBN = Tri::NBN;
    alg::denseMat K_tri(NBN, NBN, 0.0);
    vector <double> L(NBN);
    integrales(fem, tri, K_tri, L);   
    assemblage<Tri>(tri, K_tri, L, K, Lw);
    }


for (int s=0; s<SEG; s++){
    Seg &seg = fem.seg[s];
    const int NBN = Seg::NBN;
    alg::denseMat K_seg(NBN, NBN, 0.0);
    vector <double> L(NBN);
    integrales(fem, seg, L);    
    assemblage<Seg>(seg, K_seg, L, K, Lw);
    }

auto t2 = std::chrono::high_resolution_clock::now();

std::chrono::duration<double,std::micro> micros = t2-t1;
std::cout << boost::format("%5t assembling %50T. ") << micros.count() << " microsecondes\n";

t1 = std::chrono::high_resolution_clock::now();

//alg::sparseMat Kr(Kw);
K.collect();

std::vector<size_t> ld;
std::vector<double> Vd(NOD, 0.);
vector<bool> cld(NOD);
for (int s=0; s<SEG; s++){
    mapType::iterator it;
    Seg &seg = fem.seg[s];
    const int NBN = Seg::NBN;
    const int reg = seg.reg;
    pair <string,int> p;
    mapType &cond = fem.dir;

    p=make_pair("V",reg);
    it=cond.find(p);
    if (it != cond.end()){ 
       double V=it->second;
       for (int ie=0; ie<NBN; ie++){
           int i= seg.ind[ie];
           Vd[i]=V;
           ld.push_back(i);
           }
       }
    }

std::sort(ld.begin(), ld.end());
auto last = std::unique(ld.begin(), ld.end());
ld.erase(last, ld.end());

/* equilibrage des lignes */
/*
for (int i=0; i<NOD; i++){
    gmm::linalg_traits<gmm::row_matrix<write_vector > >::const_sub_row_type row = mat_const_row(Kw, i);
    double norme=gmm::vect_norminf(row);
    Lw[i]/=norme;
    for (write_vector::const_iterator it=vect_const_begin(row); it!=vect_const_end(row); ++it)
        Kw(i, it->first)/=norme;
    }
*/

std::vector<double> Lr(NOD,0.0);
Lr.assign(Lw.begin(),Lw.end());

alg::mult(K, Xw, Lw);
alg::scaled(Lw, -1.0, Lr); //Lr = -Lw

t2 = std::chrono::high_resolution_clock::now();
micros = t2-t1;
std::cout<< boost::format("%5t conditions %50T. ") << micros.count() << " microsecondes" << std::endl;

alg::iteration iter(1e-6);
iter.set_maxiter(MAXITER);
iter.set_verbosity(false);

t1 = std::chrono::high_resolution_clock::now();

Xw.clear();
Xw.resize(NOD);
double res = alg::cg_dir(K,Xw,Lr,Vd,ld,iter); // Conjugate gradient with dirichlet conditions and diagonal preconditionner

t2 = std::chrono::high_resolution_clock::now();
micros = t2-t1;
std::cout << boost::format("%5t solving (cg dir) %50T. ") << std::fixed << std::setprecision(1) << micros.count() << " microsecondes, in "<< iter.get_iteration() << " iterations, residu = "<< std::scientific << std::setprecision(6) << res << std::endl; 

for (unsigned int i=0; i<NOD; i++) {
    Node &node = fem.node[i];
    node.sol=Xw[i];
    }
}
