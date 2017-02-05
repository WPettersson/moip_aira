EXDIR = ../Examples

all: 2KP50 2AP05 3AP05 3KP10 4AP05 4KP10 moip_2_30_1_knapsack;

%:
	@echo -n "Testing $* "
	@$(TARGETDIR)/aira -p $(EXDIR)/$** -o $*.out
	@diff -w -I'seconds\|solved' $(EXDIR)/$*.out $*.out && echo -n "t1-c1 passed " && rm $*.out || echo -n "t1-c1 failed "
	@$(TARGETDIR)/aira -t 2 -p $(EXDIR)/$**.lp -o $*.out
	@diff -w -I'seconds\|solved' $(EXDIR)/$*.out $*.out && echo -n "t2-c1 passed " && rm $*.out || echo -n "t2-c1 failed "
	@$(TARGETDIR)/aira -c 2 -p $(EXDIR)/$**.lp -o $*.out
	@diff -w -I'seconds\|solved' $(EXDIR)/$*.out $*.out && echo -n "t1-c2 passed " && rm $*.out || echo -n "t1-c2 failed "
	@echo

moip_2_30_1_knapsack:
	@echo -n "Testing moip_2_30_1_knapsack "
	@$(TARGETDIR)/aira -p $(EXDIR)/moip_2_30_1_knapsack* -o moip_2_30_1_knapsack.out
	@diff -w -I'seconds\|solved' $(EXDIR)/moip_2_30_1_knapsack.out moip_2_30_1_knapsack.out && echo -n "t1-c1 passed " && rm moip_2_30_1_knapsack.out || echo -n "t1-c1 failed "
	@$(TARGETDIR)/aira -t 2 -p $(EXDIR)/moip_2_30_1_knapsack.mop -o moip_2_30_1_knapsack.out
	@diff -w -I'seconds\|solved' $(EXDIR)/moip_2_30_1_knapsack.out moip_2_30_1_knapsack.out && echo -n "t2-c1 passed " && rm moip_2_30_1_knapsack.out || echo -n "t2-c1 failed "
	@$(TARGETDIR)/aira -c 2 -p $(EXDIR)/moip_2_30_1_knapsack.mop -o moip_2_30_1_knapsack.out
	@diff -w -I'seconds\|solved' $(EXDIR)/moip_2_30_1_knapsack.out moip_2_30_1_knapsack.out && echo -n "t1-c2 passed " && rm moip_2_30_1_knapsack.out || echo -n "t1-c2 failed "
	@echo

4AP05:
	@echo -n "Testing 4AP05 "
	@$(TARGETDIR)/aira -p $(EXDIR)/4AP05.lp -o 4AP05.out
	@diff -w -I'seconds\|solved' $(EXDIR)/4AP05.out 4AP05.out && echo -n "t1-c1 passed " && rm 4AP05.out || echo -n "t1-c1 failed "
	@$(TARGETDIR)/aira -t 2 -p $(EXDIR)/4AP05.lp -o 4AP05.out
	@diff -w -I'seconds\|solved' $(EXDIR)/4AP05.out 4AP05.out && echo -n "t2-c1 passed " && rm 4AP05.out || echo -n "t2-c1 failed "
	@$(TARGETDIR)/aira -c 2 -p $(EXDIR)/4AP05.lp -o 4AP05.out
	@diff -w -I'seconds\|solved' $(EXDIR)/4AP05.out 4AP05.out && echo -n "t1-c2 passed " && rm 4AP05.out || echo -n "t1-c2 failed "
	@$(TARGETDIR)/aira -t 4 -p $(EXDIR)/4AP05.lp -o 4AP05.out
	@diff -w -I'seconds\|solved' $(EXDIR)/4AP05.out 4AP05.out && echo -n "t4 passed " && rm 4AP05.out || echo -n "t4 failed "
	@$(TARGETDIR)/aira -t 4 --perms 0,1,6,7 -p $(EXDIR)/4AP05.lp -o 4AP05.out
	@diff -w -I'seconds\|solved' $(EXDIR)/4AP05.out 4AP05.out && echo -n "t4-altperms passed " && rm 4AP05.out || echo -n "t4-altperms failed "
	@echo

4KP10:
	@echo -n "Testing 4KP10 "
	@$(TARGETDIR)/aira -p $(EXDIR)/4KP10.lp -o 4KP10.out
	@diff -w -I'seconds\|solved' $(EXDIR)/4KP10.out 4KP10.out && echo -n "t1-c1 passed " && rm 4KP10.out || echo -n "t1-c1 failed "
	@$(TARGETDIR)/aira -t 2 -p $(EXDIR)/4KP10.lp -o 4KP10.out
	@diff -w -I'seconds\|solved' $(EXDIR)/4KP10.out 4KP10.out && echo -n "t2-c1 passed " && rm 4KP10.out || echo -n "t2-c1 failed "
	@$(TARGETDIR)/aira -c 2 -p $(EXDIR)/4KP10.lp -o 4KP10.out
	@diff -w -I'seconds\|solved' $(EXDIR)/4KP10.out 4KP10.out && echo -n "t1-c2 passed " && rm 4KP10.out || echo -n "t1-c2 failed "
	@$(TARGETDIR)/aira -t 4 -p $(EXDIR)/4KP10.lp -o 4KP10.out
	@diff -w -I'seconds\|solved' $(EXDIR)/4KP10.out 4KP10.out && echo -n "t4 passed " && rm 4KP10.out || echo -n "t4 failed "
	@$(TARGETDIR)/aira -t 4 --perms 0,1,6,7 -p $(EXDIR)/4KP10.lp -o 4KP10.out
	@diff -w -I'seconds\|solved' $(EXDIR)/4KP10.out 4KP10.out && echo -n "t4-altperms passed " && rm 4KP10.out || echo -n "t4-altperms failed "
	@echo

2KP50:
	@echo -n "Testing 2KP50 "
	@$(TARGETDIR)/aira -p $(EXDIR)/2KP50.lp -o 2KP50.out
	@diff -w -I'seconds\|solved' $(EXDIR)/2KP50.out 2KP50.out && echo -n "t1-c1 passed " && rm 2KP50.out || echo -n "t1-c1 failed "
	@$(TARGETDIR)/aira -t 2 -p $(EXDIR)/2KP50.lp -o 2KP50.out
	@diff -w -I'seconds\|solved' $(EXDIR)/2KP50.out 2KP50.out && echo -n "t2-c1 passed " && rm 2KP50.out || echo -n "t2-c1 failed "
	@$(TARGETDIR)/aira -c 2 -p $(EXDIR)/2KP50.lp -o 2KP50.out
	@diff -w -I'seconds\|solved' $(EXDIR)/2KP50.out 2KP50.out && echo -n "t1-c2 passed " && rm 2KP50.out || echo -n "t1-c2 failed "
	@$(TARGETDIR)/aira -t 2 -s -p $(EXDIR)/2KP50.lp -o 2KP50.out
	@diff -w -I'seconds\|solved' $(EXDIR)/2KP50.out 2KP50.out && echo -n "t2-split passed " && rm 2KP50.out || echo -n "t2-split failed "
	@echo

2AP05:
	@echo -n "Testing 2AP05 "
	@$(TARGETDIR)/aira -p $(EXDIR)/2AP05.lp -o 2AP05.out
	@diff -w -I'seconds\|solved' $(EXDIR)/2AP05.out 2AP05.out && echo -n "t1-c1 passed " && rm 2AP05.out || echo -n "t1-c1 failed "
	@$(TARGETDIR)/aira -t 2 -p $(EXDIR)/2AP05.lp -o 2AP05.out
	@diff -w -I'seconds\|solved' $(EXDIR)/2AP05.out 2AP05.out && echo -n "t2-c1 passed " && rm 2AP05.out || echo -n "t2-c1 failed "
	@$(TARGETDIR)/aira -c 2 -p $(EXDIR)/2AP05.lp -o 2AP05.out
	@diff -w -I'seconds\|solved' $(EXDIR)/2AP05.out 2AP05.out && echo -n "t1-c2 passed " && rm 2AP05.out || echo -n "t1-c2 failed "
	@$(TARGETDIR)/aira -t 2 -s -p $(EXDIR)/2AP05.lp -o 2AP05.out
	@diff -w -I'seconds\|solved' $(EXDIR)/2AP05.out 2AP05.out && echo -n "t2-split passed " && rm 2AP05.out || echo -n "t2-split failed "
	@echo
