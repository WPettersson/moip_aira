EXDIR = ../Examples

all: 2AP05 3AP05 3KP10 4AP05 4KP10 moip_2_30_1_knapsack;

%:
	@echo -n "Testing $* "
	@$(TARGETDIR)/aira -p $(EXDIR)/$** -o $*.out
	@diff -w -I'seconds\|solved' $(EXDIR)/$*.out $*.out && echo passed && rm $*.out || echo failed
