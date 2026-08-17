#!/bin/python3.9
import multiplex_wormhole
import os

## run for each panel- input is CSV with "PrimerID" and "Sequence" in fields; SequenceID in format <ID>.1.FWD
## Marten shown as an example
# convert CSV to FASTA
multiplex_wormhole.assessPanel(IN_CSV="Marten/Marten95.csv", OUT_FA="Marten/Marten_primers.fa")

# assess dimer load metrics
multiplex_wormhole.assessPanel(PRIMERS="Marten/Marten95.fa", ALL_DIMERS_dG=-8, END_DIMERS_dG=-5, BAD_DIMERS_dG=-10)
# runs MFEprimer dimer, output: Marten95_MFEprimerDimer.txt, Marten95_MFEprimerDimer_ends.txt
# creates dimer tables: Marten95_PrimerPairDimers_sum.csv, Marten95_PrimerPairDimers_Wide.csv, Marten95_PrimerPairDeltaG_mean.csv, Marten95_PrimerPairDeltaG_wide.csv

# redesign panel for lower dimer load
os.mkdir("Marten/Opt80")
multiplex_wormhole.optimizeMultiplex(PRIMER_FASTA="Marten/Marten95.fa",
                                     DIMER_SUMS="Marten/Marten95_PrimerPairDimers_sum.csv",
                                     DIMER_TABLE="Marten/Marten_PrimerPairDimers_wide.csv,
                                     OUTPATH = "Marten/Opt80/Marten75,
                                     N_LOCI = 75) # set at 80% of panel size for each (calculate by hand)


#-----AFTER RUNNING runPrimerTree STEP in R-----#
# add thermodynamics to specificity check outputs (from runPrimerTree)
# note: this overwrites current output, but should be same just with additional fields
multiplex_wormhole.offtargetThermodynamics(INFILE="Marten/Marten_primerblast.csv", OUTFILE="Marten/Marten_primerblast.csv")
