#include "fem.h"

void integrales(Fem &fem, Tri &tri, alg::denseMat &AE, vector <double> &BE)
{
const int NBN = Tri::NBN;
const int NPI = Tri::NPI;
const int reg = tri.reg;

pair <string,int> p;
map <pair<string,int>,double> &param = fem.param;

p=make_pair("epsr",reg);    double epsr = param[p];    	
p=make_pair("rho",reg);     double rho  = param[p]; 
//cout << "epsr : " << epsr  << "\trho : " << rho << endl;

/*-------------------------------------------------------*/

for (int npi=0; npi<NPI; npi++){
    double w = tri.weight[npi]; /* := pds[npi]*detJ[npi] */

    for (int ie=0; ie<NBN; ie++){
        double ai = tri.a[ie][npi];
        double dai_dx= tri.dadx[ie][npi];  
	double dai_dy= tri.dady[ie][npi];  

        BE[ie]+= rho* ai *w;

        for (int je=0; je<NBN; je++){
            double daj_dx= tri.dadx[je][npi];  
	    double daj_dy= tri.dady[je][npi];  

	    double Dai_Daj= dai_dx*daj_dx + dai_dy*daj_dy;
            AE(ie,je) += VACUUM_PERMITTIVITY*epsr* Dai_Daj *w;  // OK
	    }
	}
    }
}

void integrales(Fem &fem, Seg &seg, alg::denseMat &AE, vector <double> &BE)
{
const int NBN = Seg::NBN;
const int NPI = Seg::NPI;
const int reg = seg.reg;

pair <string,int> p;
map <pair<string,int>,double> &param = fem.param;

p=make_pair("Dn",reg);      double Dn = param[p];  
p=make_pair("sigma",reg);   double sigma = param[p]; 
//cout << "Dn : " << Dn  << "\tsigma : " << sigma << endl;
/*-------------------------------------------------------*/

for (int npi=0; npi<NPI; npi++){
    double w = seg.weight[npi]; /* := pds[npi]*detJ[npi] */

    for (int ie=0; ie<NBN; ie++){
        double ai = seg.a[ie][npi];

        BE[ie] += (Dn+sigma)*ai *w;
	}
    }
}
