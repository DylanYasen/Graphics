#!make
# -----------------------------------------------------------------------
#  makefile for vector / matrix lib 
#  Jian Chen (jichen@umbc.edu)
#
#  Tested on mac (Nov 27 2013)
# ----------------------------------------------------------------------
#  compiler options
#  ---------------------------------------------------------------------
#

MACHINE=$(shell uname -s)

ifeq ($(MACHINE),Darwin)
	OPENGL_INC=-FOpenGL
	OPENGL_LIB=-framework OpenGL -framework GLUT -framework Cocoa
else
	OPENGL_INC=-I/use/include -I/usr/X11R6/include -I/usr/local/include 
	OPENGL_LIB=-L/usr/lib64 -L/usr/local/lib 
endif

LLDLIBS= -lGL -lGLU -lglut -lXmu -lXt -lXext -lX11 -lm -lstdc++

CC = g++ -DOPENGL -Wno-deprecated -g 

INCLUDE=$(OPENGL_INC)
LIBS=$(OPENGL_LIB)  $(LLDLIBS)

# ----------------------------------------------------------------------
# the files for cjSveWidget library 
#  ---------------------------------------------------------------------
SOURCES =  MGL.cpp 
TARGET = ./libjcUtil-1.0.a 

# ----------------------------------------------------------------------
# Put all the source together 
#  ---------------------------------------------------------------------
OBJECTS	= $(SOURCES:.cpp=.o)

# ----------------------------------------------------------------------
#  now build the  library
#  ---------------------------------------------------------------------
all: $(TARGET)

MGL.o: MGL.cpp
	$(CC) $(INCLUDE) -c MGL.cpp

$(TARGET): MGL.o
	ar -ruv libjcUtil-1.0.a MGL.o

install: 
	cp ./libjcUtil-1.0.a ../../lib

clean:
	rm -f $(TARGET) .*.swp *.o *~ core
