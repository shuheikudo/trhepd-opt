#! /bin/bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=true
export OMP_PLACES=0:8:2
cp ../../str-orig-plus/src/surf.exe surf-orig-plus.exe

exe="99"
dz=0.001
t0=$((time -p ./surf-orig-plus.exe) 2>&1 | grep real | cut -d " " -f 2)
cp surf-bulkP.s surf-bulkP-orig-plus.s
rm *log.txt
echo $exe $dz 0 $t0 > test-results/time-orig-plus.txt

