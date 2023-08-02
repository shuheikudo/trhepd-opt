#! /bin/bash
rm test-results/rock.txt
dirresults=../test47/test-results
for exe in 4 5 6 7
do
	for dz in "0.68" "0.47" "0.33" "0.22" "0.15" "0.1" "0.068" "0.047" "0.033" "0.022" "0.015" "0.01" \
		#"0.0068" "0.0047" "0.0033" "0.0022" #"0.0015" "0.001"
	do
                echo $exe $dz
		python3 tizy.py $dirresults/perf-$exe-$dz.txt $exe $dz >> test-results/rock.txt
	done

done

exe=0
echo "lapack run"
for dz in "0.68" "0.47" "0.33" "0.22" "0.15" "0.1" "0.068" "0.047" "0.033" "0.022" "0.015" "0.01" \
	#"0.0068" "0.0047" "0.0033" "0.0022" #"0.0015" "0.001"
do
        echo $exe $dz
	python3 tizy.py $dirresults/perf-$exe-$dz.txt $exe $dz >> test-results/rock.txt
done

exe="99"
dz="0.01"
echo $exe $dz
python3 tizy.py $dirresults/perf-$exe-$dz.txt $exe $dz >> test-results/rock.txt

exe="999"
dz="0.001"
echo $exe $dz
python3 tizy.py surf-bulkP-orig-plus.s $exe $dz >> test-results/rock.txt


