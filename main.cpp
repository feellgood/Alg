#include <boost/format.hpp>
#include <boost/tokenizer.hpp>

#include "fem.h"

int main(int argc,char* argv[])
{
Fem fem;
std::string fileMsh;

#if defined __GNUC__
	std::cout << "gnuC detected version: " << __GNUC__ << "." << __GNUC_MINOR__ << std::endl;
	#if ((__GNUC__ > 9)||((__GNUC__ == 9)&&(__GNUC_MINOR__ >= 3)) )
		std::cout << "gnuc version Ok\n";
	#else
		std::cout << "This gnuc is too old to handle C++17.\n";
	#endif
#endif

if (argc == 2)
	{
	fileMsh = argv[1];
	std::cout << "mesh file = " << fileMsh << std::endl;	
	fem.pbname = fileMsh;	
	}
else std::cout << "missing argument : meshFile" << std::endl;
lecture(fem);
femutil(fem);
chapeaux(fem);
affichage(fem);
const int MAXITER = 500;
solve(fem,MAXITER);  

/*
savesol(fem);
savevtk(fem);
*/
return 0;
}


