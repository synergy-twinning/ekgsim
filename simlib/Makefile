COMPILEFLAGS = -g -c -std=c++0x
LINKFLAGS = -O3 -L../../lib
CPPFLAGS = -O3 -I../../include/ -I../../../include/  -I../


simlib: lib_main.o simulator.o
	ar rc lib_main.a lib_main.o simulator.o
	mv lib_main.a simlib.a
	ranlib simlib.a

lib_main.o: lib_main.cpp
	g++ ${COMPILEFLAGS} ${CPPFLAGS} -fPIC lib_main.cpp 

simulator.o: simulator.cpp simulator.h support.h Wohlfart.h columnFile.h
	g++ ${COMPILEFLAGS} ${CPPFLAGS} -fPIC simulator.cpp


clean:
	rm -rf *.o
	rm -rf *.a

