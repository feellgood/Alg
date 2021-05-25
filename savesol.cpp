#include <iostream>
#include <sstream>

#include "fem.h"

void savesol(Fem &fem)
{
string str;
ostringstream ostr; 
ostr.str("");
ostr << fem.simname  << boost::format(".sol");
str = ostr.str();
cout << boost::format("%5t saving %50T. ") << str << endl; 

ofstream fout(str, ios::out);
if (fout.fail()){
   cerr << "pb ouverture fichier " << str << "en ecriture" << endl;
   exit(1);}

const int    NOD   = fem.NOD;
const double scale = fem.scale;

for (int i=0; i<NOD; i++){
    Node &node = fem.node[i];
    double x = node.x / scale;
    double y = node.y / scale;
    double u = node.sol;
 
    fout << boost::format("%8d %+20.10f %+20.10f %+20.10e") % i % x % y % u << endl;
    }

fout.close();
}
