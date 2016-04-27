#!smake
#
#

MACHINE= $(shell uname -s)

ifeq ($(MACHINE),Darwin)
	OPENGL_INC= -FOpenGL -I/usr/X11R6/include
	OPENGL_LIB= -L./jcUtil -framework OpenGL -framework GLUT -framework Cocoa -ljcUtil-1.0
else
	OPENGL_INC= -I/usr/lib64 -I/usr/include
	OPENGL_LIB= -L./jcUtil -L/usr/lib64 -ljcUtil-1.0 -lglut  -lGLU -lGL -lXmu -lXext -lX11 -lm
endif

CXX=g++
COMPILER_FLAGS= -g

INCLUDE= -I./jcUtil $(OPENGL_INC)
LLDLIBS= $(OPENGL_LIB)

TARGETS = MyViewingPlane

default : $(TARGETS)

SRCS = MyViewingPlane.cpp 

OBJS = $(SRCS:.cpp=.o)

%.o: %.cpp
	$(CXX) -c -o $@ $(COMPILER_FLAGS)  $< $(INCLUDE)

MyViewingPlane: MyViewingPlane.o ivTrackball.o ivview.o svUtil.o
	$(CXX) $(COMPILER_FLAGS) $(LLDLIBS) MyViewingPlane.o ivTrackball.o ivview.o svUtil.o -o MyViewingPlane

clean:
	rm -f *.o $(TARGETS) *~ .*.swp .*.swo
