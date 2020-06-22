#!/usr/bin/env Rscript
rm(list=ls())
args = commandArgs(TRUE)
##############################################################
# -> After defining simulation experiment plan 
# (1_ExperimentalPlan.R)
#------------------------------------------------------------#
# Date: September 2019
# Author: Lucus Pocus
##############################################################
# Note:
##############################################################

exp_name <- "Experiment1"

tx_file <- "~/Evaluation/Experiment1/tx_sim.txt"

fasta_file <-"Experiment1.all_transcripts.fa"
gtf_file <- "Experiment1.all_transcripts.ensembl_97.gtf"

out_folder <- "~/Evaluation/Experiment1"

##############################################################
# R packages

suppressPackageStartupMessages( require(dplyr) )
suppressPackageStartupMessages( require(data.table) )
suppressPackageStartupMessages( require(stringr) )

suppressPackageStartupMessages( require(rtracklayer) )
suppressPackageStartupMessages( require(polyester) )
suppressPackageStartupMessages( require(Biostrings) )

##############################################################

tx_ratios <- fread(file = tx_file)

fasta <- readDNAStringSet(fasta_file)

rescaling_factor <- 2/100
nb_rep <- 3
nb_conditions <- 9

paired <- TRUE
bias <- ""
gc_bias <- ""
strand_specific <- TRUE
lib_sizes <- c()

simulate_experiment(fasta = fasta_file, gtf = gtf_file, 
                    read_len = read_length, 
                    paired = paired,
                    bias = bias, gc_bias = gc_bias, strand_specific = strand_specific,
                    lib_sizes = lib_sizes,
                    reads_per_transcript = round(tx_ratios$mean_count*rescaling_factor), 
                    num_reps = rep(nb_rep, nb_condtions), 
                    fold_changes = fold_changes, outdir = paste(out_folder, exp_name, sep = "/")) 
