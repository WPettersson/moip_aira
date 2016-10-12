SRC = ./src
BUILD = ./build
CPLEXDIR=/home/enigma/opt/ibm/ILOG/CPLEX_Studio1263
CC = gcc
CXX = g++
CFLAGS = -g -Wextra -std=c++11 -pedantic -I$(CPLEXDIR)/cplex/include/
CLNFLAGS = -L$(CPLEXDIR)/cplex/lib/x86-64_linux/static_pic/ -lcplex -pthread

OBJS = $(BUILD)/symgroup.o $(BUILD)/aira.o $(BUILD)/solutions.o

all: aira

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -R $(BUILD)

aira: $(BUILD) $(OBJS)
	$(CXX) $(OBJS) -o $(BUILD)/aira $(CLNFLAGS)

$(BUILD)/aira.o: $(SRC)/aira.cpp $(SRC)/env.h $(SRC)/result.h $(SRC)/solutions.h $(SRC)/symgroup.h $(SRC)/symgroup_extern.h $(SRC)/sense.h
	$(CXX) -c $(CFLAGS) -o $@ $(SRC)/aira.cpp

$(BUILD)/symgroup.o: $(SRC)/symgroup.h $(SRC)/symgroup.cpp $(SRC)/symgroup_extern.h
	$(CXX) -c $(CFLAGS) -o $@ $(SRC)/symgroup.cpp

$(BUILD)/solutions.o: $(SRC)/solutions.h $(SRC)/solutions.cpp $(SRC)/sense.h $(SRC)/result.h
	$(CXX) -c $(CFLAGS) -o $@ $(SRC)/solutions.cpp

$(BUILD)/result.o: $(SRC)/result.h $(SRC)/result.cpp
	$(CXX) -c $(CFLAGS) -o $@ $(SRC)/result.cpp

$(SRC)/symgroup_extern.h: $(SRC)/mk_symgroup.py
	$(SRC)/mk_symgroup.py

