#! /bin/bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=true
export OMP_PLACES=0:8:2
cp ../../str-lapack/src/surf.exe surf-lapack.exe
cp ../../sim-trhepd-rheed/src/surf.exe surf-opt.exe

rm test-results/perf-comp-orig.txt
rm test-results/perf-comp-acc.txt
rm test-results/time.txt

echo "opt warm-up"
./surf-opt.exe << xxx
0.1 4
xxx
for exe in 4 5 6 7
do
	for dz in "0.68" "0.47" "0.33" "0.22" "0.15" "0.1" "0.068" "0.047" "0.033" "0.022" "0.015" "0.01" \
		#"0.0068" "0.0047" "0.0033" "0.0022" #"0.0015" "0.001"
	do
		for k in 1 2 3 4 5
		do
		t0=`./surf-opt.exe << xxx
$dz $exe
xxx
`
		echo $exe $dz $t0
		if [ $k -eq 1 ]; then
			cp surf-bulkP.s test-results/perf-$exe-$dz.txt
			python3 compare.py test-results/perf-$exe-$dz.txt surf-bulkP-orig.s $exe $dz >> test-results/perf-comp-orig.txt
			python3 compare.py test-results/perf-$exe-$dz.txt surf-bulkP-str-accurate.s $exe $dz >> test-results/perf-comp-acc.txt
		fi
		echo $exe $dz $k $t0 >> test-results/time.txt
		rm *log.txt
		done
	done

done


echo "lapack warm-up"
./surf-lapack.exe << xxx
0.1
xxx

exe=0
echo "lapack run"
for dz in "0.68" "0.47" "0.33" "0.22" "0.15" "0.1" "0.068" "0.047" "0.033" "0.022" "0.015" "0.01" \
	#"0.0068" "0.0047" "0.0033" "0.0022" #"0.0015" "0.001"
do
	for k in 1 2 3 4 5
	do
	t0=`./surf-lapack.exe << xxx
$dz
xxx
`
	echo $exe $dz $t0
	if [ $k -eq 1 ]; then
		cp surf-bulkP.s test-results/perf-$exe-$dz.txt
		python3 compare.py test-results/perf-$exe-$dz.txt surf-bulkP-orig.s $exe $dz >> test-results/perf-comp-orig.txt
		python3 compare.py test-results/perf-$exe-$dz.txt surf-bulkP-str-accurate.s $exe $dz >> test-results/perf-comp-acc.txt
	fi
	echo $exe $dz $k $t0 >> test-results/time.txt
	rm *log.txt
	done
done

echo "orig warm-up"
time ./surf-orig.exe

exe="99"
for k in 1 2 3 4 5
do
	dz="0.01"
	t0=$((time -p ./surf-orig.exe) 2>&1 | grep real | cut -d " " -f 2)
	echo $exe 0.01 $t0
	if [ $k -eq 1 ]; then
		cp surf-bulkP.s test-results/perf-$exe-$dz.txt
		python3 compare.py test-results/perf-$exe-$dz.txt surf-bulkP-orig.s $exe $dz >> test-results/perf-comp-orig.txt
		python3 compare.py test-results/perf-$exe-$dz.txt surf-bulkP-str-accurate.s $exe $dz >> test-results/perf-comp-acc.txt
	fi
	echo $exe $dz $k $t0 >> test-results/time.txt
	rm *log.txt
done


