#!/bin/bash

base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/HMF_DR010_DR047/
gene_ann_dir=$base_dir/gene_ann/

out_txt=$base_dir/scripts/annotate_genes/gene_diplotypes_with_amps_HMF_DR010_DR047.txt.gz

counter=0
for i in $gene_ann_dir/*; do
	counter=$((counter+1))

	in_txt=$i/gene_diplotypes_with_amps.txt.gz
	sample_name=$(basename $i)
	echo -ne "Processing [$counter]: $sample_name\r"

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
