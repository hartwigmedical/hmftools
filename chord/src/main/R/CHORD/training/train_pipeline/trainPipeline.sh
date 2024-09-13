#!/bin/bash

script_dir=/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORDv2/processed/training/scripts/train_rf/

## Args
training_data_path=$1
out_dir=$2
master_seed=${3:-1}
do_outer_cv=${4:-1}

mkdir -p $out_dir

trainOneModel () {

	#========= Inputs/ouputs =========#
	training_data_path=$1
	out_dir=$2
	master_seed=${3:-1}
	test_data_path=$4

	r_doCvSelSamples=$script_dir/doCvSelSamples.R
	r_trainFinalOnWhitelist=$script_dir/trainFinalOnWhitelist.R

	cv_dir=$out_dir/cv_sel_samples/; mkdir -p $cv_dir
	job_dir=$out_dir/jobs/; mkdir -p $job_dir

	#========= Repeat CV to select samples =========#
	sel_sample_cv_job=$job_dir/sel_sample_cv.array.job

	if [[ ! -f ${sel_sample_cv_job}.done ]]; then

		echo "  Submitting $sel_sample_cv_job"

		sel_sample_cv_stdout=$job_dir/sel_sample_cv_stdout/; mkdir -p $sel_sample_cv_stdout

echo "#!/bin/bash
#SBATCH --job-name=sel_sample_cv
#SBATCH --output=${sel_sample_cv_stdout}/cv_sel_samples_%a.o
#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --array=1-100

seed=\$SLURM_ARRAY_TASK_ID
out_path=$cv_dir/pred_\${seed}.txt.gz

guixr load-profile ~/.guix-profile --<<EOF
Rscript $r_doCvSelSamples $training_data_path \$seed \$out_path
EOF
" > $sel_sample_cv_job

		sel_sample_cv_ID=$(sbatch --parsable $sel_sample_cv_job)

		## Mark done
		sbatch --depend=afterok:${sel_sample_cv_ID} \
		--job-name=mark_done.sel_sample_cv \
		--wrap="touch ${sel_sample_cv_job}.done" \
		--output=/dev/null
	else
		echo "  Done file exists; Skipping ${sel_sample_cv_job}"
	fi

	#========= Train final model =========#
	train_final_job=$job_dir/train_final.job

	if [[ ! -f ${train_final_job}.done ]]; then

		echo "  Submitting $train_final_job"

echo "#!/bin/bash
#SBATCH --job-name=train_final
#SBATCH --output=${job_dir}/train_final.o
#SBATCH --time=01:00:00
#SBATCH --mem=10G

guixr load-profile ~/.guix-profile --<<EOF
Rscript $r_trainFinalOnWhitelist $training_data_path $master_seed $out_dir $cv_dir $test_data_path
EOF
" > $train_final_job

		train_final_ID=$(sbatch --parsable --depend=afterok:${sel_sample_cv_ID} $train_final_job)

		## Mark done
		sbatch --depend=afterok:${train_final_ID} \
		--job-name=mark_done.train_final \
		--wrap="touch ${train_final_job}.done" \
		--output=/dev/null
	else
		echo "  Done file exists; Skipping ${train_final_job}"
	fi
}

#========= Train final model =========#
echo "## Submitting training final model"

## Copy data with training samples
#cp $training_data_path $out_dir/training_data.txt.gz

# ## Copy data with all samples
# full_data_path=$(dirname $training_data_path)/mFull_$(basename $training_data_path | sed 's/mTrain/mFull/')
# if [[ -f $full_data_path ]]; then
# 	cp $full_data_path $out_dir/full_data.txt.gz
# fi

trainOneModel \
$training_data_path \
$out_dir/ \
$master_seed

#========= Do outer CV =========#
if [[ $do_outer_cv -eq 1 ]]; then
r_spawnOuterCvData=$script_dir/spawnOuterCvData.R
outer_cv_dir=$out_dir/outer_cv/; mkdir -p $outer_cv_dir

echo "Making outer CV train test sets"
guixr load-profile ~/.guix-profile --<<EOF
Rscript $r_spawnOuterCvData $training_data_path $master_seed $outer_cv_dir
EOF

for fold_dir in $outer_cv_dir/*/; do
	echo -e "\n## Submitting training CV fold model @ $fold_dir"
	trainOneModel \
	$fold_dir/training_data.txt.gz \
	$fold_dir/ \
	$master_seed \
	$fold_dir/test_data.txt.gz
done
fi


