module load tabix
module load ucsc_tools

num_sims=$1

for i in $(seq 1 $num_sims); 
do 
	python3 scripts/simulate_data/simulate_sum_stats.py scripts/simulate_data/config/simulate_diverse_sites.config $i &
	echo $i; 

	while [ `ps -ef | grep "scripts/simulate_data/simulate_sum_stats.py" | wc -l` -ge 16 ]
	do
		sleep 10
	done

done
