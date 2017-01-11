SRC = ./src
TARGETDIR = ./build
#TARGETDIR = ./debug
#DEBUG_FLAGS = -g -DDEBUG -DDEBUG_SOLUTION_SEARCH
CPLEXDIR=/home/enigma/opt/ibm/ILOG/CPLEX_Studio127
CC = gcc
CXX = g++
CFLAGS = $(DEBUG_FLAGS) -Wextra -std=c++11 -pedantic -I$(CPLEXDIR)/cplex/include/
CLNFLAGS = -L$(CPLEXDIR)/cplex/lib/x86-64_linux/static_pic/
LIBS=-pthread -lcplex -lboost_program_options

all: executable

executable: $(TARGETDIR)/aira

OBJS = $(TARGETDIR)/symgroup.o $(TARGETDIR)/aira.o $(TARGETDIR)/solutions.o $(TARGETDIR)/result.o $(TARGETDIR)/problem.o

$(TARGETDIR):
	mkdir -p $(TARGETDIR)

clean:
	rm -Rf ./build ./debug

test: executable
	@TARGETDIR=.$(TARGETDIR) make -C tests

$(TARGETDIR)/aira: $(TARGETDIR) $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $(TARGETDIR)/aira $(CLNFLAGS)

$(TARGETDIR)/aira.o: $(SRC)/aira.cpp $(SRC)/env.h $(SRC)/result.h $(SRC)/solutions.h $(SRC)/symgroup.h $(SRC)/symgroup_extern.h $(SRC)/sense.h $(SRC)/problem.h
	$(CXX) -c $(CFLAGS) -o $@ $(SRC)/aira.cpp

$(TARGETDIR)/symgroup.o: $(SRC)/symgroup.h $(SRC)/symgroup.cpp $(SRC)/symgroup_extern.h
	$(CXX) -c $(CFLAGS) -o $@ $(SRC)/symgroup.cpp

$(TARGETDIR)/solutions.o: $(SRC)/solutions.h $(SRC)/solutions.cpp $(SRC)/sense.h $(SRC)/result.h
	$(CXX) -c $(CFLAGS) -o $@ $(SRC)/solutions.cpp

$(TARGETDIR)/result.o: $(SRC)/result.h $(SRC)/result.cpp
	$(CXX) -c $(CFLAGS) -o $@ $(SRC)/result.cpp

$(TARGETDIR)/problem.o: $(SRC)/problem.h $(SRC)/problem.cpp
	$(CXX) -c $(CFLAGS) -o $@ $(SRC)/problem.cpp

$(SRC)/symgroup_extern.h: $(SRC)/mk_symgroup.py
	$(SRC)/mk_symgroup.py

$(SRC)/symgroup.cpp: $(SRC)/mk_symgroup.py
	$(SRC)/mk_symgroup.py
