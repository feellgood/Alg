#include "fem.h"

void affichage(Fem &fem)
{
//cout << "process: " << getpid() <<"\n\n";

pair <string,int> p;
map  <pair<string,int>,double> &param = fem.param;

//int REG = fem.REG;
//cout << "\n\t regions\t\t" << REG << endl;

int NOD = fem.NOD;
cout << boost::format("%5t nodes %50T. %d\n") % NOD;

int SEG = fem.SEG; 
cout << boost::format("%5t segments %50T. %d\n") % SEG;

int TRI = fem.TRI;
cout << boost::format("%5t triangles %50T. %d\n") % TRI;

cout << boost::format("%5t surface %50T. %d\n") % fem.surf;

cout << endl;
}
