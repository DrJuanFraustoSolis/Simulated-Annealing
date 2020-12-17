#!/usr/bin/env python
import sys
import shutil
import re

patternstring = r'\![^\n\"\']*$'
patternstring2 = r'\!.*\'.*\'[^\n]*$'
patternstring3 = r'\!.*c\'[^\n]*$'
pattern = re.compile(patternstring, re.M)
pattern2 = re.compile(patternstring2, re.M)
pattern3 = re.compile(patternstring3, re.M)
for fn in sys.argv[1:]:
    shutil.copy2(fn, fn + ".bak")

    f = open(fn, 'r')
    s = f.read()
    t = re.sub(pattern, '', s)
    t = re.sub(pattern2, '', t)
    t = re.sub(pattern3, '', t)
    f.close()
    f = open(fn, 'w')
    print >> f, t
#    print t   
