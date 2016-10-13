#!/usr/bin/env python3

import os
from math import factorial    # size of symmetric group


# Creates the symgroup.{cpp,h} files with static const for symmetric groups up
# to order N
N = 5


if not os.getcwd().endswith("/src"):
    # Assume being run as src/mk_symgroup.py
    os.chdir("src")


# These are ordered in a hopefully optimal way
def next_val(n, done=-1, sofar=[]):
    if done == -1:
        done = n
#    print("n is %d, done is %d"%(n, done))
    if done == 0:
        yield sofar
    else:
        for k in range(0, n):
            if k not in sofar:
                nextlist = list(sofar)
                nextlist.append(k)
                for new in next_val(n, done-1, nextlist):
                    yield new

with open("symgroup.cpp", "w") as hfile:
    hfile.write("#include \"symgroup.h\"\n\n")
    hfile.write("const int S0[1] = { 0 };\n")
    hfile.write("const int S1[1] = { 0 };\n")
    for n in range(2, N + 1):
        count = 0
        total = factorial(n)
        hfile.write("const int S%d[%d] = {\n" % (n, total * n))
        for s in next_val(n):
            hfile.write("\t%s" % (",".join(map(str, s))))
            count += 1
            if count < total:
                hfile.write(",")
            hfile.write("\n")

        hfile.write("};\n")
    hfile.write("\nconst SymGroup S[] = {\n")
    for n in range(0, N + 1):
        count = 0
        hfile.write("\tSymGroup(%d, S%d)" % (n, n))
        count += 1
        if count < N:
            hfile.write(",")
        hfile.write("\n")
    hfile.write("};\n\n")
    hfile.write("""int SymGroup::size() const {
    int res = 1;
    for( int i = 2; i <= N_; ++i)
        res *= i;
    return res;
}
""")

with open("symgroup_extern.h", "w") as hfile:
    hfile.write("extern const SymGroup S[%d];\n" % (N + 1))
