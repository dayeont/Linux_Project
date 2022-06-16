#!/bin/bash

#SBATCH --job-name bash
#SBATCH --partition cpu
#SBATCH --output="%j.out"
#SBATCH --error="%j.err"

number_of_threads=(1 4 6 8 12)
output='parallel4.txt'
items=('./parallel3_Ofast.exe')
size_of_matrix=(64 128 512)
for item in ${items[*]}
do
	file_exe=$item
	echo $item >> $output
	echo ${number_of_threads[*]} >> $output

	for size in ${size_of_matrix[*]}
	do
		to_file=$size
		for thread in ${number_of_threads[*]}
		do
			to_file=$to_file' '`$file_exe $thread $size`
		done
		echo $to_file >> $output
	done
done
