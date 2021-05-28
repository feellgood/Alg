#include <boost/progress.hpp>
#include <boost/format.hpp>
#include <boost/tokenizer.hpp>

#include "alg.h"
#include "alg_cg_dir.h"

#include "alg_bicg.h"

#include "fem.h"

int main(int argc,char* argv[])
{
Fem fem;
std::string fileMsh;


if (argc == 2)
	{
	fileMsh = argv[1];
	std::cout << "mesh file = " << fileMsh << std::endl;	
	fem.pbname = fileMsh;	
	}
else std::cout << "missing argument : meshFile" << std::endl;
lecture(fem);

cout << "fin de lecture " << endl;

femutil(fem);
chapeaux(fem);
affichage(fem);

solve(fem);  

savesol(fem);
savevtk(fem);
return 0;
}


