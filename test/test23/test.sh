#! /bin/bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=true
export OMP_PLACES=0:8:2
cp ../../str-lapack/src/surf.exe surf-lapack.exe
t0=`./surf-lapack.exe << xxx
0.01 0
xxx
`
cp surf-bulkP.s test-results/lapack-0.01.txt
python3 compare.py test-results/lapack-0.01.txt surf-bulkP-orig.s lapack 0.01 $t0 0

cp ../../sim-trhepd-rheed/src/surf.exe surf-opt.exe
t0=`./surf-opt.exe << xxx
0.01 7
xxx
`
cp surf-bulkP.s test-results/rk4-0.01.txt
python3 compare.py test-results/rk4-0.01.txt surf-bulkP-orig.s rk4 0.01 $t0 7

t0=`./surf-opt.exe << xxx
0.01 5
xxx
`
cp surf-bulkP.s test-results/sp6-0.01.txt
python3 compare.py test-results/sp6-0.01.txt surf-bulkP-orig.s sp6 0.01 $t0 5

rm *log.txt
