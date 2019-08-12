#!/usr/bin/python

import sys

def parsespec(a):

        if a.find("_") != -1 and a.find(">") != -1:
                # Try to parse syntax like 1246_1247SL>RM
                gt_offset = a.find(">")
                newbases = a[gt_offset + 1:].strip()
                oldbases = a[gt_offset - len(newbases) : gt_offset]
                indices = [int(x) for x in a[0:gt_offset - len(newbases)].split("_")]

                if len(indices) != 2 or ((indices[1] - indices[0]) + 1) != len(newbases):
                        raise Exception("Inconsistent lengths in mutation descriptor " + a)

                return oldbases, newbases, indices[0]

        old_base = a[0]
        new_base = a[-1]
        base_pos = int(a[1:-1])
    
        if base_pos < 1:
                raise Exception("base pos < 1 in " + a)

        return old_base, new_base, base_pos
	
def applymutation(in_bases, old_bases, new_bases, base_pos):

        base_pos -= 1

        if len(in_bases) <= (base_pos + (len(old_bases) - 1)):
                raise Exception("Input too short (%d) for mutation position (%d)" % (len(in_bases), base_pos+1))
        
        if in_bases[base_pos : base_pos + len(old_bases)] != old_bases:
                raise Exception("Base at position %d (%s) does not match mutation spec %s" % (base_pos+1, in_bases[base_pos], old_bases))

        return in_bases[:base_pos] + new_bases + in_bases[base_pos + len(new_bases):]


