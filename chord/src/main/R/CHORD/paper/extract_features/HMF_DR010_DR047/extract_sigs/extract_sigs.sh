#!/bin/bash

base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/HMF_DR010_DR047/

manifest_path=$base_dir/manifest/manifest.txt

out_dir=$base_dir/matrices/
job_dir=$base_dir/jobs/extract_sigs/; mkdir -p $job_dir; cd $job_dir

R_SCRIPT=$base_dir/scripts/extract_sigs/extract_sigs.R

counter=0
tail $manifest_path -n +2 | while read sample_name germ som sv gene_cnv cnv; do
	counter=$((counter+1))

	job_file=$job_dir/extract_sigs_${sample_name}.job

echo "#!/bin/bash
guixr load-profile ~/.guix-profile --<<EOF
Rscript $R_SCRIPT $sample_name $som $sv
EOF
" > $job_file

	#qsub -cwd -l h_rt=1:00:00 -l h_vmem=10G $job_file

	sbatch --time=01:00:00 --mem=10G --job-name=$(basename $job_file) --output=${job_file}.o $job_file

	#if [[ $counter -eq 1 ]]; then break; fi
done
