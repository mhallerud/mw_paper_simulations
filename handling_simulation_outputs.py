#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:11:59 2026

@author: maggiehallerud
"""

import glob
import pandas as pd
import os
import csv



def process_file(f):
    outdir = os.path.dirname(os.path.dirname(f))
    dataset = outdir.split("_")[0]
    targets = outdir.split("_")[1].replace("N","")
    if 'dG' in f:
        deltaG = "True"
    else:
        deltaG = "False"
    df = pd.read_csv(f)
    df['Filepath'] = f
    df['Templates'] = dataset
    df['N_LOCI'] = targets
    df['DeltaG'] = deltaG
    return df

summs = glob.glob("Marten/*/3_OptimizedMultiplexes/*RunSummary.csv")

normsims = [x for x in summs if 'Run0' not in x]
order1sims = [x for x in summs if 'Run0' in x]
order2sims = glob.glob("*/*/3_OptimizedMultiplexes/*RunSummary.csv")
order3sims = glob.glob("*/*/*/3_OptimizedMultiplexes/*RunSummary.csv")

dfs = [process_file(f) for f in summs]
result = pd.concat(dfs)
result.to_csv("ALL_SUMMARY_INFO.csv")





def process_logs(OUTPATH):
    log_files = glob.glob("*.log")
    out=[["LogFile", "StartTime","EndTime"]]
    log = []
    for l in log_files:
        with open(l, "r") as file:
            for row in file:
                log.append(row.strip())
        start = log[0]
        end = log[len(log)-1]
        start = start.replace("START TIME: ", "")
        end = end.replace("END TIME: ", "")
        out.append([l, start, end])
    with open("LogTimes.csv", "w") as file:
        writer = csv.writer(file)
        for row in out:
            writer.writerow(row)

