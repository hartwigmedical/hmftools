#!/bin/bash

base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/PCAWG_2020/
gene_ann_dir=$base_dir/gene_ann/

out_txt=$base_dir/scripts/annotate_genes/gene_diplotypes_with_amps_PCAWG.txt.gz

counter=0
for i in $gene_ann_dir/*; do
	counter=$((counter+1))

	sample_name=$(basename $i)
	echo -ne "Processing [$counter]: $sample_name\r"

	in_txt=$i/gene_diplotypes_with_amps.txt.gz
	if [[ ! -f $in_txt ]]; then continue; fi

	## Write header using first file
	if [[ $counter -eq 1 ]]; then
		zcat $in_txt | head -n 1 | 
		awk '{print "sample""\t"$0}' | 
		gzip -c > $out_txt
	fi

	zcat $in_txt | tail -n +2 | 
	#grep BRCA | 
	awk -v sample_name="$sample_name" '{print sample_name"\t"$0}' | 
	gzip -c >> $out_txt

done

echo -e '\n'
