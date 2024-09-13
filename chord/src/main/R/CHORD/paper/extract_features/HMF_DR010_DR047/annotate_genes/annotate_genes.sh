#!/bin/bash

base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/HMF_DR010_DR047/

manifest_path=$base_dir/manifest/manifest.txt
parent_out_dir=$base_dir/gene_ann/; mkdir -p $parent_out_dir

RSCRIPT=$base_dir/scripts/annotate_genes/detGeneStatuses_exec.R
JAVA=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/dep/jre1.8.0_191/bin/java
BED_FILE=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/inst/misc/genes_chord_paper.bed

counter=0
cat $manifest_path | tail -n +2 | while read sample_name germ_vcf som_vcf sv_vcf purple_gene_cnv purple_cnv; do
	counter=$((counter+1))

	echo -e "\n######### [$counter] $sample_name #########"

	out_dir=$parent_out_dir/$sample_name/; mkdir -p $out_dir

	job_dir=$out_dir/jobs/; mkdir -p $job_dir
	job_file=$job_dir/${sample_name}_detGeneStatuses.job
	done_file=${job_file}.done

echo "#!/bin/bash
guixr load-profile ~/.guix-profile --<<EOF
Rscript $RSCRIPT $out_dir $germ_vcf $som_vcf $purple_gene_cnv $purple_cnv $sample_name $JAVA $BED_FILE &&
touch $done_file
EOF
" > $job_file
	
	if [[ ! -f $done_file ]]; then
		sbatch --time=01:00:00 --mem=10G --job-name=$(basename $job_file) --output=${job_file}.o $job_file
	else
		echo "Skipping: $(basename $job_file)"
	fi

	#if [[ $counter -eq 1 ]]; then break; fi

done
