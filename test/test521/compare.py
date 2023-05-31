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
f1 = open(sys.argv[2], 'r')

nn0 = skipcomment(f0).split()
nn1 = skipcomment(f1).split()
ndom0 = int(nn0[0])
ndom1 = int(nn1[0])
nangle0 = int(nn0[1])
nangle1 = int(nn1[1])
nlat0 = int(nn0[2])
nlat1 = int(nn1[2])
if ndom0 != ndom1 or nangle0 != nangle1 or nlat0 != nlat1:
    print(nn0, nn1)
    print('error')
    exit(1)

# ignore ndom

# ignore ih, ik
readn(f0, 2*nlat0)
readn(f1, 2*nlat0)

mm = 0.
a=[]
for i in range(nangle0):
    d0 = readn(f0, 1+nlat0)
    d0 = np.array([float(x) for x in d0])
    d1 = readn(f1, 1+nlat0)
    d1 = np.array([float(x) for x in d1])
    n0 = np.linalg.norm(d0[1:])
    if n0==0.:
        n0 = 1.
    e0 = np.linalg.norm(d0[1:] - d1[1:])
    e1 = e0/n0
    mm = max(n0, mm)
    a.append([d0[0], e0, e1])

for x in a:
    a0 = x[0]
    a1 = x[1]
    a2 = x[2]
    print(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], a0, '{0:20.15e}'.format(a1/mm), '{0:20.15e}'.format(a2))


