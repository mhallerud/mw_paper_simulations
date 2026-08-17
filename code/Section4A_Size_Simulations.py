## First, subsets were made by pulling randomly from each dataset to create subsets for 100, 500, 1000, 1500, and 2000 templates (as dataset size allowed)
## Then, optimizations were run in parallel on Oregon State University's High-Performance Computing cluster
## NOTE on usage (ran from Python scripts instead of full install): 
#python3.9 ../multiplex_wormhole/src/multiplex_wormhole/multiplexWormhole.py = multiplex-wormhole
#python3.9 ../multiplex_wormhole/src/multiplex_wormhole/tabulate_dimers.py = mw-tabulate-dimers
#python3.9 ../multiplex_wormhole/src/multiplex_wormhole/multiple_run_optimizations.py = mw-mult-optimizations



## RUN PANEL DESIGN SIMULATIONS WITH STANDARD ALGORITHM
# set up commands for each simulation in bash
for TEMPLATES in `ls -1 *Templates*.csv`; do
  for N in 20 50 100 150 200 250 300; do
    PREFIX="$(echo ${TEMPLATES} | cut -d '.' -f1 | sed 's/Templates_subset//' | cut -d'_' -f1)"
    echo "python3.9 ../multiplex_wormhole/src/multiplex_wormhole/multiplexWormhole.py -t ${TEMPLATES} -n ${N} -o ${PREFIX}_N${N} -p ${PREFIX}_N${N}_standard --threads 2";
  done
done > run_sims.txt

# submit as batch job on the HPC for parallel processing (2 threads each)
cat run_sims.txt | hqsub -t array -q fw -r mw_sims_log -p 2



## RUN PANEL DESIGN SIMULATIONS WITH DELTA-G ALGORITHM
# set up function for rerunning tabulate dimers, optimization steps (saved in file run_optimOnly.sh):
# able to skip primer design, predict dimers steps because already run above
#---------------------------run_optimOnly.py-------------------------------#
#!/bin/bash
DIR=$1
N_LOCI="$(echo $DIR | cut -d'N' -f2)"

python3.9 ../multiplex_wormhole/src/multiplex_wormhole/tabulate_dimers.py -a ${DIR}/2_PredictedDimers/MFEprimerDimers.txt -e ${DIR}/2_PredictedDimers/MFEprimerDimers_ends.txt \
-o ${DIR}/2_PredictedDimers/PrimerPairDimers --deltaG
python3.9 ../multiplex_wormhole/src/multiplex_wormhole/multiple_run_optimization.py -f ${DIR}/1_PrimerDesign/FilteredPrimers.fa -d ${DIR}/2_PredictedDimers/PrimerPairDimers_mean.csv \
-t ${DIR}/2_PredictedDimers/PrimerPairDimers_wide.csv -o ${DIR}/3_OptimizedMultiplexes/${DIR}_deltaG -n ${N_LOCI} -r 10 --deltaG
#-------------(END)-----------#

# make run_optimOnly.txt file to run on all output dirs from above simulations, e.g.:
#---------------------------run_optimOnly.txt--------------------------------#
run_optimOnly.sh Marten100_N20
run_optimOnly.sh Marten100_N50
run_optimOnly.sh Marten500_N20
run_optimOnly.sh Marten500_N50
run_optimOnly.sh Marten500_N100
run_optimOnly.sh Marten500_N150
run_optimOnly.sh Marten500_N200
....
#-------------(END)-----------#


# then run as batch job:
cat run_optimOnly.txt | hqsub -t array -q fw -r mw_sims_dg_log
