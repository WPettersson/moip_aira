SRC = ./src
BUILD = ./build
CPLEXDIR=/home/enigma/opt/ibm/ILOG/CPLEX_Studio1263
CC = gcc
CFLAGS = -Wextra -std=c99 -pedantic -g -DIL_STD -I$(CPLEXDIR)/cplex/include/ `pkg-config --cflags glib-2.0`
CLNFLAGS = -L$(CPLEXDIR)/cplex/lib/x86-64_linux/static_pic/ -lcplex -lm -lpthread `pkg-config --libs glib-2.0`

all: aira

init:
	mkdir -p ./build

clean:
	rm -R ./build

aira: init aira.o
	$(CC) $(BUILD)/aira.o -o $(BUILD)/aira $(CLNFLAGS)
aira.o: $(SRC)/aira.c
	$(CC) -c $(CFLAGS) $(SRC)/aira.c -o $(BUILD)/aira.o
