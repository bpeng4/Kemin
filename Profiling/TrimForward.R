setwd("/work/benson/bpeng4/Kemin/Kemin_All/Kemin_All_16S")
suppressMessages({
library("tidyverse")
library("dada2")
library("gridExtra")
library("devtools")
})
no_of_cores = 36
metadata = "/work/benson/bpeng4/Kemin/Kemin_All/SampleSheet_Final.csv"
metadata_df = read.csv(metadata)
metadata_df = metadata_df[,-1]
fnFs <- metadata_df$fq1

sample.names <- metadata_df$Sample_Name

filt_path <- paste0(getwd(), '/filtered') # don't change this 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))

FORWARD_TRUNC <- 250 # determine from quality plots

out <- filterAndTrim(fnFs, filtFs,
                     truncLen=FORWARD_TRUNC, 
                     trimLeft=20, maxEE=2, 
                     multithread=FALSE,
                     matchIDs=TRUE, compress=TRUE, 
                     verbose=TRUE)

derepFs <- derepFastq(filtFs, n = 1e7, verbose = TRUE)

names(derepFs) <- sample.names

errF <- learnErrors(filtFs, verbose=TRUE, multithread=no_of_cores)

dadaFs <- dada(derepFs, err=errF, pool=TRUE, multithread=no_of_cores, 
               verbose=TRUE)

seqtab <- makeSequenceTable(dadaFs)


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=no_of_cores, verbose=TRUE)

save(metadata_df, seqtab.nochim, file = "/work/benson/bpeng4/Kemin/Kemin_All/intermediate/Forward.rda")
