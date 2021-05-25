#include <iostream>
#include <sstream>

#include "fem.h"

void savevtk(Fem &fem)
{
string str;
ostringstream ostr; 
ostr.str("");
ostr << fem.simname<< ".vtk";
str = ostr.str();

cout << boost::format("%5t saving %50T. ") << str << endl; 

ofstream fout(str.c_str(), ios::out);
if (!fout){
    cerr << "pb ouverture fichier : " << str << "en ecriture" << endl;
    exit(1);}

ostr.str("");
ostr << fem.simname;
str = ostr.str();

const int NOD = fem.NOD;
const int TRI = fem.TRI;
const int NBN = Tri::NBN;

fout << "# vtk DataFile Version 2.0" << endl;
fout << str << endl;
fout << "ASCII" << endl;
fout << "DATASET UNSTRUCTURED_GRID" << endl;
fout << "POINTS "<< NOD << " float" << endl;

for (int i=0; i<NOD; i++){
    Node &node = fem.node[i];
    double x   = node.x;
    double y   = node.y;
    double z   = 0;
    fout << boost::format("%+20.10e %+20.10e %+20.10e") % x % y % z << endl;
    }

fout << boost::format("CELLS %8d %8d") % TRI % (4*TRI) << endl;
for (int t=0; t<TRI; t++){
    Tri &tri = fem.tri[t];
    fout << setw(8) << NBN;
    for (int ie=0; ie<NBN; ie++)
        fout << setw(8) << tri.ind[ie];
    fout << endl;
    }

fout << boost::format("CELL_TYPES %8d") % TRI << endl;
for (int t=0; t<TRI; t++)
    fout << setw(8) << 5 << endl;
   
fout << boost::format("POINT_DATA %8d") % NOD << endl;
fout << boost::format("SCALARS V float %d") % 1 << endl;
fout << "LOOKUP_TABLE my_table" << endl;

for (int i=0; i<NOD; i++){
    Node &node = fem.node[i];
    double u = node.sol;
    fout << boost::format("%+20.10e") % u << endl;
    }
}    
