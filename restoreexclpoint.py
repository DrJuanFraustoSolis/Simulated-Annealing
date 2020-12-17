#!/usr/bin/env python
import sys
import os
for fn in sys.argv[1:]:
  try:
    os.rename(fn + ".bak", fn)
  except:
    pass
