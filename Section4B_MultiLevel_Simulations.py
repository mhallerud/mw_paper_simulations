#!/bin/python3.9
# load dependencies
import sys
import os
import pandas as pd
import traceback

# load multiplex wormhole modules
sys.path.append("../multiplex_wormhole/src/multiplex_wormhole")
from multiplexWormhole import main as mw
from multiple_run_optimization import main as opt
from helpers import CSVtoFasta


# function to subset data, run simulations
def main(TEMPLATES, N1, NAME, CPUs=1):
        # read in CSV
        if "Marten" in TEMPLATES:
                df = pd.read_csv(TEMPLATES, sep="\t")
        else:
             	df = pd.read_csv(TEMPLATES)
        # set up possible orders for running optimization
        files = [[1000,50,20], [1000,20,50], [50,1000,20], [50,20,1000], [20,50,1000], [20,1000,50]]
        N_locis = [[100,20,10], [100,10,20], [20,100,10], [20,10,100], [10,20,100], [10,100,20]]
        i=0
	while i<1: #adjust if wanting to run multiple runs (setting to 1 for now to enable parallelization...)
                print("Running data subset "+str(i+1))
                # create random, non-overlapping subsets of 1000, 50, and 20 templates
                sub1000 = df.sample(n=N1, replace=False)
                df_remove = pd.concat([df, sub1000], axis=0, ignore_index=True).drop_duplicates(subset=["SEQUENCE_ID","SEQUENCE_TEMPLATE","SEQUENCE_TARGET"], keep=False, ignore_index=True)
                sub50 = df_remove.sample(n=50, replace=False)
                df_remove = pd.concat([df_remove, sub50], axis=0, ignore_index=True).drop_duplicates(subset=["SEQUENCE_ID","SEQUENCE_TEMPLATE","SEQUENCE_TARGET"], keep=False, ignore_index=True)
                sub20 = df_remove.sample(n=20, replace=False)
                sub1000.to_csv(NAME+"_Sub1000.csv")
                sub50.to_csv(NAME+"_Sub50.csv")
                sub20.to_csv(NAME+"_Sub20.csv")
                # run all orders on this set of data...
                for j in range(6):
                        # loop through orders 1-6
                        print("...Order "+str(j))
                        outdir1 = NAME+"_N"+str(files[j][0])
                        outdir2 = os.path.join(NAME+"_N"+str(files[j][0]), "N"+str(files[j][1]))
                        outdir3 = os.path.join(outdir2, "N"+str(files[j][2]))
                        keeplist1 = os.path.join(os.path.join(outdir1, "3_OptimizedMultiplexes"),
                                                 "Run01_"+os.path.basename(outdir1)+"_primers.fasta")
                        keeplist2 = os.path.join(os.path.join(outdir2, "3_OptimizedMultiplexes"),
                                                 "Run01_"+os.path.basename(outdir2)+"_primers.fasta")
                        n2 = N_locis[j][0]+N_locis[j][1]
                        n3 = n2 + N_locis[j][2]
                        # run level 1 optimization
                        try:
                            	print("......"+str(files[j][0])+" "+str(N_locis[j][0]))
                              PREFIX = "Run"+str(i)+os.path.basename(outdir1)
                              mw(TEMPLATES=NAME+"_Sub"+str(files[j][0])+".csv", N_LOCI=N_locis[j][0], OUTDIR=outdir1, 
                                 PREFIX=PREFIX, N_RUNS=1, THREADS=CPUs)
                        except Exception as err:
                                print("Error/warning with run " + str(N_locis[j][0]) + " and " + str(files[j][0]))
                                print(traceback.format_exc())
                        # convert level 1 output to FASTA to be used as level 2 keeplist
                        CSVtoFasta.main(os.path.join(outdir1,"3_OptimizedMultiplexes/Final_Primers/")+
                                        "Run"+str(i)+os.path.basename(outdir1)+"_Run01_primers.csv", 
                                        keeplist1)
                        # run level 2 optimization
                        try:
                                PREFIX = "Run"+str(i)+os.path.basename(outdir2)
                                print("......"+str(files[j][1])+" "+str(N_locis[j][1]))
                                mw(TEMPLATES=NAME+"_Sub"+str(files[j][1])+".csv", N_LOCI=n2, OUTDIR=outdir2, 
                                   PREFIX="Run"+str(i)+os.path.basename(outdir2), KEEPLIST_FA=keeplist1, N_RUNS=1, 
                                   THREADS=CPUs)
                        except Exception as err:
                                print("Error/warning with run " + str(N_locis[j][1]) + " and " + str(files[j][1]))
                                print(traceback.format_exc())
                        # convert level 2 output to FASTA to be used as level 3 keeplist
                        CSVtoFasta.main(os.path.join(outdir2, "3_OptimizedMultiplexes/Final_Primers/")+
                                        "Run"+str(i)+os.path.basename(outdir2)+"_Run01_primers.csv", 
                                        keeplist2)
                        try:
                                PREFIX = "Run"+str(i)+os.path.basename(outdir3)
                                print("......"+str(files[j][2])+" "+str(N_locis[j][2]))
                                mw(TEMPLATES=NAME+"_Sub"+str(files[j][2])+".csv", N_LOCI=n3, 
                                   OUTDIR=outdir3, PREFIX="Run"+str(i)+os.path.basename(outdir3), KEEPLIST_FA=keeplist2, 
                                   N_RUNS=1, THREADS=CPUs)
                        except Exception as err:
                                print("Error/warning with run " + str(N_locis[j][2]) + " and " + str(files[j][2]))
                                print(traceback.format_exc())
                i+=1 #move to next set of orders



if __name__=="__main__":
        main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
