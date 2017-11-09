CPPFLAGS = -Wall -O3
BOOST = -I ../../Libraries/boost_1_56_0
SRC1 = operator.cc lattice.cc simulation.cc qecsim.cc colex.cc RandomStuff.cpp RngStream.cpp
SRC2 = parseInput.cc
OBJS = ../blossom5-v2.05.src/misc.o ../blossom5-v2.05.src/PMduals.o ../blossom5-v2.05.src/PMexpand.o ../blossom5-v2.05.src/PMinit.o ../blossom5-v2.05.src/PMinterface.o ../blossom5-v2.05.src/PMmain.o ../blossom5-v2.05.src/PMrepair.o ../blossom5-v2.05.src/PMshrink.o ../blossom5-v2.05.src/GEOM/GPMinit.o ../blossom5-v2.05.src/GEOM/GPMinterface.o ../blossom5-v2.05.src/GEOM/GPMkdtree.o ../blossom5-v2.05.src/GEOM/GPMmain.o ../blossom5-v2.05.src/MinCost/MinCost.o 

all:
	g++ $(CPPFLAGS) $(BOOST) $(SRC1) $(OBJS) -o qecsim 
	g++ $(CPPFLAGS) $(SRC2) -o parseInput 

doc: doc/timestamp.dummy

doc/timestamp.dummy:
	touch doc/timestamp.dummy
	$(DOXYGEN)

clean:
	rm -f *.o qecsim *~
	rm -rf ./doc/*
