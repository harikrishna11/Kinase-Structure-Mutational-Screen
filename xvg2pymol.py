#!/usr/bin/env python

import sys

for l in sys.stdin:
    if l.startswith("#"):
        continue
    if l.startswith("@"):
        continue
    bits = l.split()
    sys.stdout.write(bits[1] + "\n")


