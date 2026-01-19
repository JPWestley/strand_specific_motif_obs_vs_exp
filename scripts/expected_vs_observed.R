###################################################################
### CALCULATING EXPECTED VS OBSERVED FOR ALL MOTIFS AND GENOMES ###
###################################################################

# 1.0 load dependencies ####

library(Biostrings)
library(data.table)
library(tidyverse)

# 2.0 defining functions ####

## 2.1 expand IUPAC function ####

expand_IUPAC <- function(motif, max_ambiguity = 10) {
  nts <- strsplit(toupper(motif), "")[[1]]
  
  choices <- lapply(nts, function(b) {
    if (!b %in% names(IUPAC_CODE_MAP))
      stop(paste("Unknown nucleotide code:", b))
    strsplit(IUPAC_CODE_MAP[[b]], "")[[1]]
  })
  
  n_ambiguity <- sum(lengths(choices) > 1)
  if (n_ambiguity > max_ambiguity)
    stop(paste0("Too many ambiguous positions (", n_ambiguity, ")"))
  
  apply(expand.grid(choices), 1, paste0, collapse = "")
}

## 2.2 prep genome function ####

prepare_genome <- function(genome_file) { # give this function the path to a fasta file(s)
  genome_seq <- readDNAStringSet(genome_file)[[1]] # read the seq of a fasta file for a given genome
  
  tri_freq <- oligonucleotideFrequency(genome_seq, 3) # count observations of all possible trinucleotides for the genome 
  tri_freq <- tri_freq / sum(tri_freq) # convert counts to a freq for each trinucleotide
  
  di_freq <- oligonucleotideFrequency(genome_seq, 2) # count observations of all possible dinucleotides for the genome 
  di_freq <- di_freq / sum(di_freq) # convert counts to a freq for each dinucleotide
  
  list( # save important info for the genome
    genome_name   = tools::file_path_sans_ext(basename(genome_file)),
    genome_seq    = genome_seq,
    tri_freq      = tri_freq,
    di_freq       = di_freq,
    genome_length = length(genome_seq)
  )
}

## 2.3 calculate expected vs observed function (rm_motif_enrichment_strandwise) ####

rm_motif_enrichment_strandwise <- function(genome_obj, motif_seq) { # give this function a genome_obj as output by prepare_genome function, and a motif sequence
  
  motif_seq <- toupper(motif_seq) # ensure case is correct
  rc_seq    <- as.character(reverseComplement(DNAString(motif_seq))) # get reverse complement of the motif 
  L <- nchar(motif_seq) # get motif length
  
  if (L < 3) stop("Motif must be at least 3 bp")
  
  tri_freq   <- genome_obj$tri_freq # get the relevant tri_freqs for the genome
  di_freq    <- genome_obj$di_freq # get the relevant di_freqs for the genome
  genome_seq <- genome_obj$genome_seq # get the genome seq
  
  expected_prob_markov2 <- function(motif) { # define a function that gets the expected freq of a motif
    nts  <- strsplit(motif, "")[[1]] # split motif into nucleotides
    tris <- vapply(1:(L - 2), function(i) # generates all trinucleotides in the motif, e.g., GATTCT gets GAT, ATT, TTC, TCT
      paste0(nts[i:(i + 2)], collapse = ""), character(1))
    dins <- vapply(2:(L - 2), function(i) # generates the dinucleotides that overlap between each trinucleotide in motif, e.g., GATTCT gets AT, TT, CT
      paste0(nts[i:(i + 1)], collapse = ""), character(1))
    p <- unname(tri_freq[tris[1]]) # gets observed probability of first trinucleotide in motif (GAT)
    if (length(tris) > 1) { # if motif is longer than 3 bases
      for (i in 2:length(tris)) { # for next trinucleotides in the motif (first in example loop is ATT)
        p <- p * (tri_freq[tris[i]] / di_freq[dins[i - 1]]) # so the first loop here gets probability of first trinucleotide (GAT) * (probability of next trinucleotide in sliding window (ATT)/probability of overlapping dinucleotide (AT)) # then loop continues for next trinucleotide (TTC)  
      }
    }
    p # return final p for the motif
  }
  
  compute_expected <- function(motif) { # define function that runs expected_prob_markov2 for all possible motifs a given motif can expand out to e.g., for GATNTCT we expand to GATATCT, GATCTCT, GATGTCT, GATTTCT
    concrete_seqs <- expand_IUPAC(motif) # expand motif
    total_prob <- sum(vapply(concrete_seqs, expected_prob_markov2, numeric(1))) # sum probabilities of all motifs a given motif can expand to
    total_prob * (genome_obj$genome_length - L + 1) # convert probability to a count of expected occurrences for given genom length
  }
  
  obs_fwd <- countPattern(motif_seq, genome_seq, fixed = FALSE) # gets observed count of motifs on forward strand
  obs_rev <- countPattern(rc_seq, genome_seq, fixed = FALSE) # gets observed count of motifs on reverse strand
  
  exp_fwd <- compute_expected(motif_seq) # gets expected counts of motif on forward strand
  exp_rev <- compute_expected(rc_seq) # gets expected counts of motif on reverse strand
  
  data.frame( # builds dataframe of results 
    motif      = motif_seq,
    strand     = c("forward", "reverse"),
    observed   = c(obs_fwd, obs_rev),
    expected   = c(exp_fwd, exp_rev),
    enrichment = c(obs_fwd / exp_fwd, obs_rev / exp_rev),
    genome     = genome_obj$genome_name,
    genome_length = genome_obj$genome_length,
    row.names  = NULL
  )
}

# 3.0 run the fucntions for all genomes and motifs ####

motifs <- read.csv("datafiles/raw/motifs.csv")

motifs <- as.vector(motifs$motif)

genomes <- list.files("datafiles/genomes", full.names = TRUE)

# precompute genomes 
genome_objs <- lapply(genomes, prepare_genome)

# preallocate result list
results <- vector("list", length(genome_objs) * length(motifs))
k <- 1

# counter for genomes
genome_counter <- 0  

startTime <- Sys.time()

for (g in genome_objs) { # run the loop across all given genomes and motifs
  
  for (m in motifs) {
    results[[k]] <- rm_motif_enrichment_strandwise(
      genome_obj = g,
      motif_seq  = m
    )
    k <- k + 1
  }
  
  genome_counter <- genome_counter + 1
  message(genome_counter, " genomes analysed")
}

# bind results
data <- rbindlist(results)

endTime <- Sys.time()
print(endTime - startTime)

# write data
write.csv(data, file = "datafiles/processed/obvs_vs_exp.csv")
