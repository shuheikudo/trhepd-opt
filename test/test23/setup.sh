#! /bin/bash
cp ../../str-orig/src/bulk.exe bulk-orig.exe
cp ../../str-orig/src/surf.exe surf-orig.exe
./bulk-orig.exe
./surf-orig.exe
cp surf-bulkP.s surf-bulkP-orig.s

cp ../../str-lapack/src/surf.exe surf-lapack.exe
cp ../../sim-trhepd-rheed/src/surf.exe surf-opt.exe
t0=`./surf-opt.exe << xxx
0.001 5
xxx
`
cp surf-bulkP.s surf-bulkP-str-accurate.s
