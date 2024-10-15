#!/bin/bash

base_dir=/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/datasets/processed/PCAWG_2020/
#data_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD_data/PCAWG_2020/

manifest_path=$base_dir/manifest/manifest_all_som_data.txt

out_dir=$base_dir/matrices/
done_dir=$out_dir/done_files/
job_dir=$base_dir/jobs/extract_sigs/; mkdir -p $job_dir; cd $job_dir

R_SCRIPT=$base_dir/scripts/extract_sigs/extract_sigs.R

counter=0
tail $manifest_path -n +2 | while read sample_name snv_mnv indel sv cnv; do
	counter=$((counter+1))

	echo "[$counter] $sample_name"

	job_file=$job_dir/extract_sigs_${sample_name}.job
	done_file=$done_dir/${sample_name}.done

	if [[ ! -f $done_file ]]; then

echo "#!/bin/bash
guixr load-profile ~/.guix-profile --<<EOF
Rscript $R_SCRIPT $out_dir $sample_name $snv_mnv $indel $sv
EOF
" > $job_file
			
		sbatch --time=01:00:00 --mem=10G --output=${job_file}.o $job_file
	else
		echo "Skipping $sample_name"
	fi

	#if [[ $counter -eq 1 ]]; then break; fi
done
