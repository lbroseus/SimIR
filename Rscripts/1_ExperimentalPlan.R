#!/usr/bin/env Rscript
rm(list=ls())
args = commandArgs(TRUE)
##############################################################
# For generating the experimental plan 
# IR + PROXIMAL EXON SKIPPING
#------------------------------------------------------------#
# Author: Lucus Pocus
# Date: September 2019
# Updates: 
# octobre 2019
# novembre 2019: simulation d'épissage alternatif (proximal exon)
##############################################################

# Sorties:
# -> description des gènes considérés
# -> description des transcrits simulés (canoniques + IR + Skip)
# -> description des introns retenus

# -> fichier fasta pour les simulations

##############################################################
# Paths and Parameters

outFolder <- "~/Evaluation/Experiment1"
experiment <- "Experiment1"

gtf_file <- "Homo_sapiens.GRCh38.97.gtf"

##############################################################
# R packages

suppressPackageStartupMessages( require(dplyr) )
suppressPackageStartupMessages( require(rtracklayer) )
suppressPackageStartupMessages( require(GenomicFeatures) )

##############################################################
# Calcul des comptages moyens à simuler par couple de Tx (Can, IR)
# Choix des couples Tx? Tx canonique + choix aléatoire de l'intron?

# Scénarios:
# Niveau de rétention:
# 5 scénarios :       0-0.1             0.2 0.4 0.6 0.8
#               unif:       trianglular:
# Niveau d'expression:
# 3 scénarios: LOW MEDIUM HIGH
# Nombre de réplicats par scénario: 3
# 5x3 x3 = 45 mesures pour chaque intron, 45 échantillon de RNA-seq à simuler

# Exclusion des gènes chevauchants
# Annotation des gènes:
# gene_name gene_id, chr start end strand width nbAnnotPseudogenes?
# Annotation des introns:
# [identifiant] chr start end strand width gène nbIntrons intron_rank width_5p width_3p | gc_5p gc_3p gc_deviation mappability? junc_type? 

# import du fichier d'annotation

gtf <- rtracklayer::import(con = gtf_file)

# restriction au chromosome 1 
gtf <- gtf[seqnames(gtf) %in% c(1)] %>% data.frame() 

# restriction aux gènes et transcrits 'protein_coding'

gtf_prot.tx <- gtf %>% filter(type == "transcript" & gene_biotype == "protein_coding" & transcript_biotype == "protein_coding") %>%
  dplyr::select(transcript_name, transcript_id, 
                seqnames, start.tx = start, end.tx = end, strand, transcript_biotype, gene_biotype, gene_name , gene_id)

# compter le nombre d'exon par transcrit
# exclusion des transcrits uni-exoniques / des gènes sans transcrit multi-exonique

gtf_prot.exon <- gtf %>% filter(type == "exon" & gene_biotype == "protein_coding" & transcript_biotype == "protein_coding") %>%
  group_by(transcript_id) %>% summarize(nbExons = n()) %>% data.frame()
gtf_prot.tx <- merge(gtf_prot.tx, gtf_prot.exon, by = "transcript_id")
rm(gtf_prot.exon)

gtf_prot.tx <- gtf_prot.tx %>% filter(nbExons>3)
summary(gtf_prot.tx$nbExons)

gtf_prot.gene <- gtf %>% filter(type == "gene" & gene_biotype == "protein_coding" & gene_name %in% gtf_prot.tx$gene_name) %>% 
  dplyr::select(gene_name , gene_id, seqnames, start, end, strand, gene_biotype)
nrow(gtf_prot.gene)

# compter le nombre d'isoformes annotées par gène

nbIso <- gtf_prot.tx %>% group_by(gene_name) %>% summarize(nbIsoforms = n()) %>% data.frame()
gtf_prot.gene <- merge(gtf_prot.gene, nbIso, by = "gene_name"); rm(nbIso)

summary(gtf_prot.gene$nbIsoforms)

# exclusion des gènes chevauchants

hits <- findOverlaps(query = GRanges(gtf_prot.gene), subject = GRanges(gtf_prot.gene))
hits <- data.frame(hits) %>% group_by(queryHits) %>% summarize(nbOverlap = n()-1) %>% data.frame()
hits <- hits %>% filter(nbOverlap>0)

gtf_prot.gene <- gtf_prot.gene[-hits$queryHits,]
gtf_prot.tx <- gtf_prot.tx %>% filter(gene_name %in% gtf_prot.gene$gene_name)

rm(hits)

# choix aléatoire uniforme d'un transcrit multi-exonique par gène

selected <- gtf_prot.tx %>% group_by(gene_name) %>% summarize(canonical_tx = sample(x = transcript_name, size = 1)) %>% data.frame()

gtf_prot.tx <- gtf_prot.tx %>% filter(transcript_name %in% selected$canonical_tx)
nrow(gtf_prot.tx)

#----> liste des transcrits canoniques sélectionnés

out_gtf_file <- paste(outFolder, "/", experiment, ".canonical_transcripts.ensembl_97.gtf", sep = "")
rtracklayer::export(object = filter(gtf, transcript_name %in% gtf_prot.tx$transcript_name | (type == "gene" & gene_name %in% gtf_prot.tx$gene_name)),
                    con = out_gtf_file)

gtf_prot.gene <- merge(gtf_prot.gene, gtf_prot.tx, by = c("gene_name", "gene_id" , "seqnames", "strand", "gene_biotype"))
rm(gtf_prot.tx)

rm(gtf)

#--------------------------------------------------------------------------#
# Construction des transcrits alternatifs (IR et Skip)
#--------------------------------------------------------------------------#

# II: can vs IR vs proxymal exon skipping 

txdb <-  makeTxDbFromGFF(file = out_gtf_file)

# définition des introns correspondants

transcripts.can <- exonsBy(x = txdb, by = "tx", use.names = TRUE)
intronsByTx <- intronsByTranscript(x = txdb, use.names = TRUE)
length(intronsByTx)

rm(txdb)

#----> annotation des introns considérés transcript_name_IR_numIntron; type = IR
#----> annotation des exons cassettés transcript_name_Skip_numExon; type = Skip

# Choix aléatoire d'un intron pour chaque transcrit canonique

for(i in 1:length(intronsByTx)) intronsByTx[[i]]$intron_rank <- 1:length(intronsByTx[[i]])

selected_introns <- endoapply(X = intronsByTx, FUN = function(x) x[sample(x = 1:length(x), size = 1)])

stopifnot( length(selected_introns) == nrow(gtf_prot.gene) )

# Choix aléatoire de l'exon flanquant qui jouera le rôle d'exon-cassette

exon_index <- unlist( selected_introns )$intron_rank
exon_index <- tapply(exon_index, 1:length(exon_index), function(x) x + sample(x = c(0,1), size = 1))

#----> liste/gtf des transcrits ALT construits

transcripts.skipping <- mapply(FUN = function(x,y) x[-y], transcripts.can, exon_index, USE.NAMES = TRUE, SIMPLIFY = TRUE)
transcripts.skipping <- GRangesList(transcripts.skipping)

names(transcripts.skipping) <- paste(names(transcripts.skipping), "Skip", exon_index, sep = "_")

#----> liste/gtf des transcrits IR construits

transcripts.ir <- union(transcripts.can, selected_introns)
names(transcripts.ir) <- paste(names(transcripts.ir), "IR", unlist( selected_introns )$intron_rank, sep = "_")

selected_introns <- data.frame(selected_introns) %>% dplyr::select(transcript_id = group_name, 
                                                                   start.ir = start, end.ir = end,
                                                                   intron_width = width, 
                                                                   intron_rank)

selected_introns$transcript_id.ir <- paste(selected_introns$transcript_id, "IR", selected_introns$intron_rank, sep = "_")
selected_introns$transcript_id.skip <- paste(selected_introns$transcript_id, "Skip", exon_index, sep = "_")

gtf_prot.gene <- merge(gtf_prot.gene, selected_introns)

head(gtf_prot.gene)

rm(selected_introns, intronsByTx, exon_index)

out_file <- paste(outFolder, "/", experiment, ".by_gene_info.ensembl_97.txt", sep = "")
fwrite(gtf_prot.gene, file = out_file, sep = "\t") 

#-----> Transcrits Canoniques et Transcrits IR
transcripts <- c(transcripts.can, transcripts.ir, transcripts.skipping)

# Shape it as gtf...(tedious)
exons <- data.frame(transcripts) %>% dplyr::select(-group, transcript_id = group_name, start, end, width, strand, exon_rank, -exon_name, -exon_id, -width)
exons$type <- "exon"
transcripts <- exons %>% group_by(transcript_id) %>% summarize(seqnames = seqnames[1],
                                                               start = min(start), end = max(end), strand = strand[1]) %>% data.frame()
exons$root_tx <- str_remove(exons$transcript_id, pattern = "\\_.*")
transcripts$root_tx <- str_remove(transcripts$transcript_id, pattern = "\\_.*")

genes <- gtf_prot.gene %>% dplyr::select(seqnames, start, end, strand, gene_id, gene_name, transcript_id) %>% mutate(exon_rank = NA, type = "gene")

transcripts <- merge(transcripts, dplyr::select(genes, transcript_id, gene_id, gene_name), by.x = "root_tx", by.y = "transcript_id")
exons <- merge(exons, dplyr::select(genes, transcript_id, gene_id, gene_name), by.x = "root_tx", by.y = "transcript_id")

transcripts <- transcripts %>% dplyr::select(-root_tx) %>% mutate(exon_rank = NA, type = "transcript")
exons <- exons %>% dplyr::select(-root_tx)

genes$transcript_id <- NA

new_gtf <- rbind.data.frame(genes, rbind.data.frame(transcripts, exons)) %>% arrange(seqnames, start)
rm(exons, genes)

out_gtf_file <- paste(outFolder, "/", experiment, ".all_transcripts.ensembl_97.gtf", sep = "")
export(con = out_gtf_file, object = new_gtf)

rm(transcripts.can, transcripts.ir, transcripts.skipping)

#---------------------------------------------------#
#Génération du fichier fasta correspondant
#---------------------------------------------------#

suppressPackageStartupMessages( require(Biostrings) )
require(BSgenome.Hsapiens.UCSC.hg38)

genome <- getBSgenome(genome = "hg38", masked=FALSE)

transcripts <- makeTxDbFromGFF(out_gtf_file)
transcripts <- exonsBy(transcripts, use.names = TRUE)
seqlevelsStyle(transcripts) <- seqlevelsStyle(genome)

sequences <- extractTranscriptSeqs(genome, transcripts = transcripts)

out_fasta_file <- paste(outFolder, "/", experiment, ".all_transcripts.fa", sep = "")
writeXStringSet(x = sequences, filepath = out_fasta_file)

rm(genome, sequences, transcripts)

#---------------------------------------------------------------------------------------#
# Sélection des comptages moyens par transcrit pour chaque "condition"
#---------------------------------------------------------------------------------------#

retention_levels <- matrix(c(0.05, 0.2, 0.2, 0.4, 0.4, 0.6, 0.6, 0.8, 0.8, 1), ncol = 2, byrow = TRUE)
perbase_exonique_coverage <- matrix(c(0, 50, 50, 250, 250, 1500), ncol = 2, byrow = TRUE)

colnames(retention_levels) <- c("low_bound", "upper_bound")
colnames(perbase_exonique_coverage) <- c("low_bound", "upper_bound")

head(gtf_prot.gene)

# Gène -> (Can, IR) -> répartition des 9 scénarios 

scenarii <- c(1:(nrow(retention_levels)*nrow(perbase_exonique_coverage)))
scenarii 

plan <- data.frame()
for(g in 1:nrow(gtf_prot.gene)) plan <- rbind.data.frame(plan,sample(x = scenarii, replace = FALSE))

colnames(plan) <- paste("condition", 1:length(scenarii), sep = "")
rownames(plan) <- gtf_prot.gene$gene_id
head(plan)

out_file <- paste(outFolder, "/", experiment, ".sim_plan.txt", sep = "")
fwrite(x = plan, file = out_file)

scenarii <- data.frame(id = 1:length(scenarii), 
                       ratio_class = rep(c(1:nrow(retention_levels)), nrow(perbase_exonique_coverage)), 
                       cov_class = sort( rep(c(1:nrow(perbase_exonique_coverage)), nrow(retention_levels))))
scenarii

transcripts <- new_gtf %>% filter(type == "transcript") %>% dplyr::select(-type, -exon_rank) 
transcripts$class <- "canonical"
transcripts$class[grep(pattern = "IR", x = transcripts$transcript_id)] <- "IR"
transcripts$class[grep(pattern = "Skip", x = transcripts$transcript_id)] <- "Skip"

drawExpressionLevels <- function(scenario_id){
  
  drawnCoverage <- runif(n = 1, min = perbase_exonique_coverage[scenarii$cov_class[scenario_id],1], 
                         max = perbase_exonique_coverage[scenarii$cov_class[scenario_id],2])
  return(drawnCoverage)
  
}

drawIRatio <- function(scenario_id){
  
  drawnRatio <- runif(n = 1, min = retention_levels[scenarii$ratio_class[scenario_id],1], 
                      max = retention_levels[scenarii$ratio_class[scenario_id],2])
  return(drawnRatio)
  
}

drawExonRatio <-  function(min, max){
  
  drawnRatio <- runif(n = 1, min = min, max = max)
  return(drawnRatio)
  
}


tmp <- data.frame()
for(tx in 1:nrow(gtf_prot.gene)){
  tmp <- rbind.data.frame(tmp, data.frame(gtf_prot.gene[tx,], 
                                          condition = scenarii$id, 
                                          sim_scenario = sample(x = scenarii$id, size = nrow(scenarii), replace = FALSE)))
}
head(tmp,25)

#---------------------------#

sim_info.df <- tmp %>% 
  dplyr::select(gene_id, gene_name, transcript_id, seqnames, start=start.tx, end=end.tx, strand, condition, 
                sim_scenario) %>%
  mutate(tx_type = "Can")

selected_introns <- tmp %>% 
  dplyr::select(gene_id, gene_name, transcript_id = transcript_id.ir, seqnames, start=start.ir, end=end.ir, strand, condition, 
                sim_scenario) %>%
  mutate(tx_type = "IR")

selected_skip <- tmp %>% 
  dplyr::select(gene_id, gene_name, transcript_id = transcript_id.skip, seqnames, start=start.tx, end=end.tx, strand, condition, 
                sim_scenario) %>%
  mutate(tx_type = "Skip")

sim_info.df <- rbind.data.frame(sim_info.df, selected_introns, selected_skip)
sim_info.df <- sim_info.df %>% arrange(gene_id, transcript_id)

sim_info.df %>% head(50)

rm(tmp)

# ->  by gene_id, draw can_coverage, draw fractions
sim_info.df <- sim_info.df %>% 
  group_by(gene_id, condition) %>% 
  mutate(can_coverage = drawExpressionLevels(sim_scenario),
         ir_fraction = drawIRatio(sim_scenario),
         skip_fraction = drawExonRatio(min = 0, max = 1)) %>% ungroup() %>% data.frame()

sim_info.df <- sim_info.df %>% arrange(gene_id, transcript_id)

out_file <- paste(outFolder, "/", experiment, ".sim_info.ensembl_97.txt", sep = "")

fwrite(sim_info.df, file = out_file, sep = "\t") 


###########################################
