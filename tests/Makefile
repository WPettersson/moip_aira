EXDIR = ../Examples

all: 2AP05 3AP05 3KP10 4AP05 4KP10 moip_2_30_1_knapsack;

%:
	@echo -n "Testing $* "
	@$(TARGETDIR)/aira -p $(EXDIR)/$** -o $*.out
	@diff -w -I'seconds\|solved' $(EXDIR)/$*.out $*.out && echo -n "t1-c1 passed " && rm $*.out || echo -n "t1-c1 failed "
	@$(TARGETDIR)/aira -t 2 -p $(EXDIR)/$** -o $*.out
	@diff -w -I'seconds\|solved' $(EXDIR)/$*.out $*.out && echo -n "t2-c1 passed " && rm $*.out || echo -n "t2-c1 failed "
	@$(TARGETDIR)/aira -c 2 -p $(EXDIR)/$** -o $*.out
	@diff -w -I'seconds\|solved' $(EXDIR)/$*.out $*.out && echo -n "t1-c2 passed " && rm $*.out || echo -n "t1-c2 failed "
	@echo

2AP05:
	@echo -n "Testing 2AP05 "
	@$(TARGETDIR)/aira -p $(EXDIR)/2AP05.lp -o 2AP05.out
	@diff -w -I'seconds\|solved' $(EXDIR)/2AP05.out 2AP05.out && echo -n "t1-c1 passed " && rm 2AP05.out || echo -n "t1-c1 failed "
	@$(TARGETDIR)/aira -t 2 -p $(EXDIR)/2AP05* -o 2AP05.out
	@diff -w -I'seconds\|solved' $(EXDIR)/2AP05.out 2AP05.out && echo -n "t2-c1 passed " && rm 2AP05.out || echo -n "t2-c1 failed "
	@$(TARGETDIR)/aira -c 2 -p $(EXDIR)/2AP05* -o 2AP05.out
	@diff -w -I'seconds\|solved' $(EXDIR)/2AP05.out 2AP05.out && echo -n "t1-c2 passed " && rm 2AP05.out || echo -n "t1-c2 failed "
	@$(TARGETDIR)/aira -t 2 -s -p $(EXDIR)/2AP05* -o 2AP05.out
	@diff -w -I'seconds\|solved' $(EXDIR)/2AP05.out 2AP05.out && echo -n "t2-split passed " && rm 2AP05.out || echo -n "t1-c1 failed "
	@echo
