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

prepare_genome <- function(genome_file) {
  genome_seq <- readDNAStringSet(genome_file)[[1]]
  
  tri_freq <- oligonucleotideFrequency(genome_seq, 3)
  tri_freq <- tri_freq / sum(tri_freq)
  
  di_freq <- oligonucleotideFrequency(genome_seq, 2)
  di_freq <- di_freq / sum(di_freq)
  
  list(
    genome_name   = tools::file_path_sans_ext(basename(genome_file)),
    genome_seq    = genome_seq,
    tri_freq      = tri_freq,
    di_freq       = di_freq,
    genome_length = length(genome_seq)
  )
}

## 2.3 calculate expected vs observed function (rm_motif_enrichment_strandwise) ####

rm_motif_enrichment_strandwise <- function(genome_obj, motif_seq) {
  
  motif_seq <- toupper(motif_seq)
  rc_seq    <- as.character(reverseComplement(DNAString(motif_seq)))
  L <- nchar(motif_seq)
  
  if (L < 3) stop("Motif must be at least 3 bp")
  
  tri_freq   <- genome_obj$tri_freq
  di_freq    <- genome_obj$di_freq
  genome_seq <- genome_obj$genome_seq
  
  expected_prob_markov2 <- function(motif) {
    nts  <- strsplit(motif, "")[[1]]
    tris <- vapply(1:(L - 2), function(i)
      paste0(nts[i:(i + 2)], collapse = ""), character(1))
    dins <- vapply(2:(L - 2), function(i)
      paste0(nts[i:(i + 1)], collapse = ""), character(1))
    
    p <- unname(tri_freq[tris[1]])
    if (length(tris) > 1) {
      for (i in 2:length(tris)) {
        p <- p * (tri_freq[tris[i]] / di_freq[dins[i - 1]])
      }
    }
    p
  }
  
  compute_expected <- function(motif) {
    concrete_seqs <- expand_IUPAC(motif)
    total_prob <- sum(vapply(concrete_seqs, expected_prob_markov2, numeric(1)))
    total_prob * (genome_obj$genome_length - L + 1)
  }
  
  obs_fwd <- countPattern(motif_seq, genome_seq, fixed = FALSE)
  obs_rev <- countPattern(rc_seq, genome_seq, fixed = FALSE)
  
  exp_fwd <- compute_expected(motif_seq)
  exp_rev <- compute_expected(rc_seq)
  
  data.frame(
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

for (g in genome_objs) {
  
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
