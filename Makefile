CXX = g++  -c  -Wno-deprecated -ggdb -Wunused-result
INC = -I ../include
LIB = -L -lm

OBJECTS= main.o 
	

.SUFFIXES: .cpp.o

		
all: 		$(OBJECTS)
	g++ -Wno-deprecated $(OBJECTS) -o  exe $(LIB) 

.cc.o:
		$(CXX) $(INC) $<

clean: 
		rm -f *.o *~

cleanall: 
		rm -f *.o *~ exe


main.o : alg.h
