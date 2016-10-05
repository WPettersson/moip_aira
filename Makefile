SRC = ./src
BUILD = ./build
CPLEXDIR=/home/enigma/opt/ibm/ILOG/CPLEX_Studio1263
CC = gcc
CXX = g++
CFLAGS = -g -Wextra -std=c++11 -pedantic -I$(CPLEXDIR)/cplex/include/
CLNFLAGS = -L$(CPLEXDIR)/cplex/lib/x86-64_linux/static_pic/ -lcplex -pthread

all: aira

init:
	mkdir -p $(BUILD)

clean:
	rm -R $(BUILD)

aira: init aira.o
	$(CXX) $(BUILD)/aira.o -o $(BUILD)/aira $(CLNFLAGS)
aira.o: $(SRC)/aira.cpp
	$(CXX) -c $(CFLAGS) $(SRC)/aira.cpp -o $(BUILD)/aira.o
