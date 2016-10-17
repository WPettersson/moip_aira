#!/usr/bin/env python3

import fileinput
import re

lines = []

reg = re.compile("Thread (\d)")

for line in fileinput.input():
    m = reg.match(line)
    if m:
        thread_id = int(m.group(1))
    lines.append({"id": thread_id, "line": line})

lines.sort(key=lambda x: x["id"])

for line in lines:
    print(line["line"].rstrip())
