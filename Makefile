SRC = ./src
BUILD = ./build
CPLEXDIR=/home/enigma/opt/ibm/ILOG/CPLEX_Studio1263
CC = gcc
CXX = g++
CFLAGS = -g -Wextra -std=c++11 -pedantic -I$(CPLEXDIR)/cplex/include/
CLNFLAGS = -L$(CPLEXDIR)/cplex/lib/x86-64_linux/static_pic/ -lcplex -pthread

all: aira

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -R $(BUILD)

aira: $(BUILD) $(BUILD)/symgroup.o $(BUILD)/aira.o
	$(CXX) $(BUILD)/aira.o $(BUILD)/symgroup.o -o $(BUILD)/aira $(CLNFLAGS)

$(BUILD)/%.o: $(SRC)/%.cpp $(SRC)/symgroup.h $(SRC)/symgroup_extern.h
	$(CXX) -c $(CFLAGS) -o $@ $<

$(SRC)/symgroup_extern.h:
	$(SRC)/mk_symgroup.py

