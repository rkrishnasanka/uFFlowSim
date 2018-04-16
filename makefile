CXX = g++
WARNINGS = -Wall -pedantic -Wformat -Wcast-align
CFLAGS = -O2 $(WARNINGS)

LODEPNGSRC = lodepng/lodepng.cpp
LODEPNGOBJ = lodepng/lodepng.o

SRCF = Fluid.cpp
OBJF = $(subst .cpp,.o,$(SRCF))

all: curved-boundaries

curved-boundaries: $(OBJF) $(LODEPNGOBJ)
	$(CXX) -o $@ $^ $(LDFLAGS) 
	
lodepng:
	$(CXX) $(CFLAGS) -c -o $(LODEPNGOBJ) $(LODEPNGSRC)

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $^

clean:
	rm -f $(DIR_5)curved-boundaries
	rm -f lodepng/lodepng.o
