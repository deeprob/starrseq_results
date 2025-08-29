# Submit the first job
job0=$(sbatch dataset.sh)
job0_id=$(echo $job0 | awk '{print $4}')
echo "Submitted dataset.sh with Job ID: $job0_id"

# Submit the fourth job with dependency on the third
job1=$(sbatch --dependency=afterok:$job0_id run.sh)
job1_id=$(echo $job1 | awk '{print $4}')
echo "Submitted run.sh with Job ID: $job1_id"

# Submit the second job with dependency on the first
job2=$(sbatch --dependency=afterok:$job1_id evaluate.sh)
job2_id=$(echo $job2 | awk '{print $4}')
echo "Submitted evaluate.sh with Job ID: $job2_id"

# Submit the third job with dependency on the second
job3=$(sbatch --dependency=afterok:$job2_id interpret.sh)
job3_id=$(echo $job3 | awk '{print $4}')
echo "Submitted interpret.sh with Job ID: $job3_id"
