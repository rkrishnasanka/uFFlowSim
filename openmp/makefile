CXX = g++
WARNINGS = -Wall -pedantic -Wformat -Wcast-align
CFLAGS = -O3 -std=c++17 -fopenmp $(WARNINGS)

LODEPNGSRC = lodepng/lodepng.cpp
LODEPNGOBJ = lodepng/lodepng.o

SRCF = main.cpp
FUNCTIONS = Functions.cpp
SOLIDBODY = SolidBody.cpp
FLUIDQUANTITY = FluidQuantity.cpp
FLUIDSOLVER = FluidSolver.cpp

OBJF = $(subst .cpp,.o,$(SRCF))
FUNCTIONSF = $(subst .cpp,.o,$(FUNCTIONS))
SOLIDBODYF = $(subst .cpp,.o,$(SOLIDBODY))
FLUIDQUANTITYF = $(subst .cpp,.o,$(FLUIDQUANTITY))
FLUIDSOLVERF = $(subst .cpp,.o,$(FLUIDSOLVER))


all: curved-boundaries

curved-boundaries: $(OBJF) $(FUNCTIONSF) $(SOLIDBODYF) $(FLUIDQUANTITYF) $(FLUIDSOLVERF) $(LODEPNGOBJ)
	$(CXX) -o $@ $^ $(LDFLAGS) 
	
lodepng:
	$(CXX) $(CFLAGS) -c -o $(LODEPNGOBJ) $(LODEPNGSRC)

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $^

clean:
	rm -f curved-boundaries
	rm -f lodepng/lodepng.o
	rm -f *.o
