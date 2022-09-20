#!/bin/bash
wd=/mnt/d/SJ-AAV/final_project
for file in $wd/virus_vs_plasmid/raw_fastq/*.fq.gz
do
	if test -f $file
	then
		arr=(${arr[*]} $file)
	fi
done
arr_length=${#arr[*]}
loop_num=`expr $arr_length / 2`
loop_num=`expr $loop_num - 1`
for ((i=0;i<=${loop_num};i++))
do
	val_1=`expr $i \* 2`
	val_2=`expr $val_1 + 1`
	fastq1=${arr[$val_1]}
	fastq2=${arr[$val_2]} 
	python $wd/analyze_insertion_new.py \
		--read_1=$fastq1  --read_2=$fastq2 \
	       	--backbone_fname=$wd/back_bone.fa --insert_fname=$wd/aa-seq.fa --tmp_file=$wd/tmp \
		--insert_len=37 --kmer_remove=12 --kmer_retain=20 \
		--out=$wd/virus_vs_plasmid/raw_counts
done


