#include "fem.h"
#define DEBUG 0

void lecture(Fem &fem)
{
int NOD, ELEM, tags, reg, TYP;
double scale, value;
string trash, symb;
pair <string,int> p;

//string str = fem.pbname;
ifstream msh(fem.pbname);   // ouverture du fichier probleme en lecture
if (msh.fail()){
    cerr << "Impossible d'ouvrir le fichier " << fem.pbname << endl;
    exit(1);}
 
getline(msh, trash);        // entete 3 lignes
getline(msh, trash);
getline(msh, trash);

/*
fem.REG = REG;    
*/         
//msh >> trash >> scale;      // lecture de l'echelle
//fem.scale = scale;

msh >> trash >> NOD;        // lecture des noeuds
cout << boost::format("%5t nodes %50T. %d\n") % NOD;

fem.NOD = NOD; 
fem.node.resize(NOD);
for (int i=0; i<NOD; i++){
    double x,y,z;
    msh >> trash >> x >> y >> z;
    fem.node[i].x = x;
    fem.node[i].y = y;
    }

msh >> trash;		// End Nodes
msh >> trash >> ELEM;        // lecture des elements

cout << boost::format("%5t elements %50T. %d\n") % ELEM << endl;

for (;;){
    msh >> symb;
    if (symb == "$EndElements")
        break;
    msh >> TYP >> tags >> reg;
    for (int i=1; i<tags; i++)
        msh >> trash;

    switch (TYP){
        case 1:{
            Seg seg;
            seg.reg = reg;	    
            msh >> seg.ind[0] >> seg.ind[1];
	    for (int i=0; i<2; i++)
	        seg.ind[i]--;           // passage convention Matlab/msh a C++
		
            fem.seg.push_back(seg);
            break;
	    }
        case 2:{
            Tri tri;
            tri.reg = reg;
	    
            msh >> tri.ind[0] >> tri.ind[1] >> tri.ind[2];
	    for (int i=0; i<3; i++)
	        tri.ind[i]--;           // passage convention Matlab/msh a C++
		
            fem.tri.push_back(tri);
            break;
	    }

        default:
            getline(msh,trash);
	}
    }

fem.SEG = fem.seg.size();
fem.TRI = fem.tri.size();

msh >> trash >> scale;      // lecture de l'echelle
fem.scale = scale;

for (int i=0; i<NOD; i++){
    fem.node[i].x *= scale;
    fem.node[i].y *= scale;
    }

msh >> trash;			// $Parameters
for (;;){                     // lecture des parametres
    msh >> symb;
    if (symb=="$EndParameters" || symb=="$End")
        break;
    msh >> reg >> value;
    cout << symb << '\t' << reg << '\t' << value << endl;
    p = make_pair(symb,reg);
    fem.param[p] = value;
    }

msh >> trash;			// $Dirichlet
for (;;){                     // lecture des conditions
    msh >> symb;
    if (symb=="$EndDirichlet" || symb=="$End")
        break;
    msh >> reg >> value;
    cout << symb << '\t' << reg << '\t' << value << endl;
    p = make_pair(symb,reg);
    fem.dir[p] = value;
    }

msh.close();
}
