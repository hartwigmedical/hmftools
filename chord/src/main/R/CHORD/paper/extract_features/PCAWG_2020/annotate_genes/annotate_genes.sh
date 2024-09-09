#!/bin/bash

base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/PCAWG_2020/

manifest_path=$base_dir/scripts/annotate_genes/manifest_samples_with_all_data.txt
parent_out_dir=$base_dir/gene_ann/; mkdir -p $parent_out_dir

RSCRIPT=$base_dir/scripts/annotate_genes/detGeneStatuses_exec.R
JAVA=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/dep/jre1.8.0_191/bin/java
GENES_BED_FILE=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/inst/misc/cosmic_cancer_gene_census_20200225.bed
EXONS_BED_FILE=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/inst/misc/genes_chord_paper_exons.bed.gz

counter=0
cat $manifest_path | tail -n +2 | while read sample_name germ_vcf som_vcf cnv; do
	counter=$((counter+1))

	echo -e "\n######### [$counter] $sample_name #########"

	out_dir=$parent_out_dir/$sample_name/; mkdir -p $out_dir

	job_dir=$out_dir/jobs/; mkdir -p $job_dir
	job_file=$job_dir/${sample_name}_detGeneStatuses.job
	done_file=${job_file}.done

echo "#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --job-name=$(basename $job_file)
#SBATCH --output=${job_file}.o
guixr load-profile ~/.guix-profile --<<EOF
Rscript $RSCRIPT $out_dir $germ_vcf $som_vcf $cnv $sample_name $JAVA $GENES_BED_FILE $EXONS_BED_FILE &&
touch $done_file
EOF
" > $job_file
	
	if [[ ! -f $done_file ]]; then
		sbatch $job_file
	else
		echo "Skipping: $(basename $job_file)"
	fi

	#if [[ $counter -eq 1 ]]; then break; fi

done
