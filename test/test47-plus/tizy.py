import numpy as np
import sys

def skipcomment(f):
    c = '#'
    while c == '#':
        l = f.readline()
        c = l[0]
    return l
def readn(f, n):
    nn = 0
    ret = []
    while nn < n:
        t = skipcomment(f)
        t = t.replace(",", "").split()
        ret.extend(t)
        nn = nn + len(t)
    return ret

def maxerror(d0, d1):
    r = 0.
    for i in range(len(d0)):
        r = max(r, abs(d0[i]-d1[i]))
    return r

f0 = open(sys.argv[1], 'r')

nn0 = skipcomment(f0).split()
ndom0 = int(nn0[0])
nangle0 = int(nn0[1])
nlat0 = int(nn0[2])
# ignore ndom

# ignore ih, ik
ihik = readn(f0, 2*nlat0)

for i in range(nangle0):
    d0 = readn(f0, 1+nlat0)
    d0 = np.array([float(x) for x in d0])
    for j in range(nlat0):
        print(sys.argv[1], sys.argv[2], sys.argv[3], ihik[2*j], ihik[2*j+1], d0[0], '{0:20.15e}'.format(d0[j+1]))
