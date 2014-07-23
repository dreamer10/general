# Jello cube Makefile 
# Jernej Barbic, USC

DEBUG := -g
LIBRARIES = -framework OpenGL -framework GLUT $(DEBUG)

COMPILER = g++ $(DEBUG)
COMPILERFLAGS = -O2 $(DEBUG)
OBJ := basic_math.o camera.o constants.o fps.o grid_data.o mac_grid.o main.o \
	 smoke_sim.o stb_image.o stb_image_write.o vec.o

all: smoke

smoke: $(OBJ)
	$(COMPILER) $(COMPILERFLAGS) -o $@ $^ $(LIBRARIES)

basic_math.o: basic_math.cpp *.h
	$(COMPILER) -c $(COMPILERFLAGS) basic_math.cpp
camera.o: camera.cpp *.h
	$(COMPILER) -c $(COMPILERFLAGS) camera.cpp
constants.o: constants.cpp *.h
	$(COMPILER) -c $(COMPILERFLAGS) constants.cpp
fps.o: fps.cpp *.h
	$(COMPILER) -c $(COMPILERFLAGS) fps.cpp
grid_data.o: grid_data.cpp *.h
	$(COMPILER) -c $(COMPILERFLAGS) grid_data.cpp
mac_grid.o: mac_grid.cpp *.h
	$(COMPILER) -c $(COMPILERFLAGS) mac_grid.cpp
main.o: main.cpp *.h
	$(COMPILER) -c $(COMPILERFLAGS) main.cpp
smoke_sim.o: smoke_sim.cpp *.h
	$(COMPILER) -c $(COMPILERFLAGS) smoke_sim.cpp
stb_image.o: stb_image.c *.h
	$(COMPILER) -c $(COMPILERFLAGS) stb_image.c
stb_image_write.o: stb_image_write.c *.h
	$(COMPILER) -c $(COMPILERFLAGS) stb_image_write.c
vec.o: vec.cpp *.h
	$(COMPILER) -c $(COMPILERFLAGS) vec.cpp

clean:
	-rm -rf core *.o *~ "#"*"#" test


