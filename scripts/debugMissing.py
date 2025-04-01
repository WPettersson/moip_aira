#!/usr/bin/env python

import argparse
from collections import defaultdict
import time
import subprocess
import re

FOUND_RE = re.compile(r".*, found ([\d,]+), ")

def run(filename: str, threads: int, output: str):
    start = time.time()
    result = subprocess.run(["build/src/aira", "-t", str(threads), "-p", filename, "-o", output],
                            text=True, capture_output=True)
    time_taken = time.time() - start
    solutions: list[list[int]] = []
    with open(output, "r") as infile:
        for line in infile:
            try:
                solution = [_ for _ in map(int, line.split())]
                if len(solution) >= 2:
                    solutions.append(solution)
            except ValueError:
                if "Solutions found" in line:
                    count = int(line.split()[0])
                pass
    if count != len(solutions):
        print("\n".join(str(s) for s in solutions))
        raise Exception(f"Uh oh: {count=} {len(solutions)=}")
    return result.stdout, solutions, time_taken


def do_diff(stdout: str, solutions: list[list[int]], newout: str, newsoln: list[list[int]]):
    if len(newsoln) > len(solutions):
        goodout = newout
        goodsoln = newsoln
        badout = stdout
        badsoln = solutions
    else:
        goodout = stdout
        goodsoln = solutions
        badout = newout
        badsoln = newsoln
    good_log: dict[int, list[str]] = defaultdict(list)
    good_finds: dict[tuple, int] = {}
    if len(goodout) > 10:
        for line in goodout.split("\n"):
            try:
                thread = int(line.split()[1])
            except IndexError:
                print(f"Couldn't get thread from {line.rstrip()}")
                continue
            good_log[thread].append(line)
            if match_obj := FOUND_RE.match(line):
                solution = tuple(map(int, match_obj.group(1).split(",")))
                good_finds[solution] = thread
    bad_log: dict[int, list[str]] = defaultdict(list)
    bad_finds: dict[tuple, int] = {}
    if len(badout) > 10:
        for line in badout.split("\n"):
            if not line:
                continue
            try:
                thread = int(line.split()[1])
            except IndexError:
                print(f"Couldn't get thread from {line.rstrip()}")
                continue
            bad_log[thread].append(line)
            if match_obj := FOUND_RE.match(line):
                solution = tuple(map(int, match_obj.group(1).split(",")))
                bad_finds[solution] = thread
    missing: list[list[int]] = []
    for soln in goodsoln:
        if soln not in badsoln:
            missing.append(soln)
            print(f"Missing solution: {soln}")
            tsoln = tuple(soln)
            if good_finds:
                print(f"Found by thread {good_finds[tsoln]}")
                with open("good-thread.txt", "w") as outfile:
                    for line in good_log[thread]:
                        outfile.write(f"{line}\n")
                with open("bad-thread.txt", "w") as outfile:
                    for line in bad_log[thread]:
                        outfile.write(f"{line}\n")
                print("Relevant logs written to good-thread.txt and bad-thread.txt, diff them")


def main(filename: str, threads: int=16, iterations=30):
    output_file = "debug-missing.out"
    stdout, solutions, time_taken = run(filename, threads, output=output_file)
    print(f"First run in {time_taken:.2f} seconds, {len(solutions)} solutions found")
    for _ in range(iterations):
        newout, newsoln, time_taken = run(filename, threads, output=output_file)
        print(f"New run in {time_taken:.2f} seconds, {len(newsoln)} solutions found")
        if len(newsoln) != len(solutions):
            print(f"Different results found: {len(newsoln)}")
            do_diff(stdout, solutions, newout, newsoln)
            return
    print(f"No differences found after {iterations} iterations")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--threads", default=16, type=int)
    parser.add_argument("-i", "--iterations", default=30, type=int)
    parser.add_argument("input_file")
    args = parser.parse_args()
    main(args.input_file, threads=args.threads, iterations=args.iterations)
