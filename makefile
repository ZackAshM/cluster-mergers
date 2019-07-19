# makefile for c++ programs
#
# to create binary executable "area" below, type "make area"
#
# If you don't have a "bin" folder for binaries, create it first
#
# the "#" is the comment character for makefile scripts
#
# this section sets some variables/alaises
# The syntax "$(HOME)" gives substitution of the variable HOME
# change "your_username" below to appropriate name
# ** for Mac users, this will probably be /Users/[yourname] instead of /home **
#
HOME = /home/zackashm
#HOME =/Users/gorham   # uncomment this for the mac
BIN  =  $(HOME)/bin
LOCALDIR = $(HOME)/work
HERE = $(LOCALDIR)/refinal
#
# this section for compiler flags and link libraries
# -O3 gives a high degree of code optimization, -W gives all warnings
# and -g turns on the debug options
#
GFLAGS  = -O3 -W -g 
LIBS    = -lm
CC      = g++
#
# Now for the actual compiling instructions:
# make sure that the space before $(CC)
# is a tab whenever you add a new compiler line.
# The syntax below is such that each left-justified name is
# the "target" (or final product executable program) of the make,
# and the tabbed-in lines that follow provide all of the steps needed
# to make that executable binary file and put it into ~/bin

# Note: -o means output name ___, -c means compile only (for example something with no main)
# Structure:
# (target): (dependencies)
# (1 TAB) (instructions)
# Example:
# target: dependencies target.cpp
# 	g++ target.cpp dependency.cpp -o /home/zackashm/bin/target
# dependencies: dependency.cpp
# 	g++ -c dependency.cpp
#
# In terminal:
# >make target
# >(computer stuff)
# >target
# >(target runs)

# Main commands -------------------
main: VecFRK4 cluster
	$(CC) $(GFLAGS) cluster_generator.o VecFRK4.o main.cpp -o $(HERE)/main $(LIBS)
gen_clust: cluster
	$(CC) $(GFLAGS) cluster_generator.o cluster.cpp -o $(HERE)/cluster
	

# Dependencies --------------------
cluster: cluster_generator.cpp
	$(CC) -c cluster_generator.cpp
VecFRK4: VecFRK4.cpp
	$(CC) -c VecFRK4.cpp
	

# Test Programs -------------------
single: VecFRK4 cluster
	$(CC) $(GFLAGS) cluster_generator.o VecFRK4.o single.cpp -o $(HERE)/single $(LIBS)
EM: VecFRK4 EM.cpp
	$(CC) $(GFLAGS) VecFRK4.o EM.cpp -o $(HERE)/EM $(LIBS)
tester: VecFRK4 cluster
	$(CC) $(GFLAGS) cluster_generator.o VecFRK4.o tester.cpp -o $(HERE)/tester
BarnesHut_test: cluster
	$(CC) $(GFLAGS) cluster_generator.o BarnesHut_test.cpp -o $(HERE)/BHT
mergerfirst: 
	$(CC) $(GFLAGS) mergerfirsttry.cpp -o $(HERE)/mergerfirst $(LIBS)
scratch:
	$(CC) $(GFLAGS) scratch.cpp -o $(HERE)/scratch


clean:
	rm *.o;
# Whenever you add a new program to your collection of work,
# you add new lines in the file to govern the compiling of your new program.
# end of makefile

