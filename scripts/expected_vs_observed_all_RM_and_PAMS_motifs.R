##############################################
### EXPECTED VS OBSERVED FOR ALL RM MOTIFS ###
##############################################

# 1.0 load dependencies ####

library(Biostrings)
library(tidyverse)

# 2.0 define functions ####

# expand IUPAC codes such as N into all possible motifs 

expand_IUPAC <- function(motif, max_ambiguity = 10) {
  # Split motif into letters
  nts <- strsplit(toupper(motif), "")[[1]]
  
  # Map each letter to all possible bases using IUPAC_CODE_MAP
  choices <- lapply(nts, function(b) {
    if(!b %in% names(IUPAC_CODE_MAP)) stop(paste("Unknown nucleotide code:", b))
    strsplit(IUPAC_CODE_MAP[[b]], "")[[1]]
  })
  
  # Count total ambiguity positions
  n_ambiguity <- sum(sapply(choices, length) > 1)
  if(n_ambiguity > max_ambiguity) {
    stop(paste0("Too many ambiguous positions (", n_ambiguity,
                "), reduce motif ambiguity or increase max_ambiguity"))
  }
  
  # Expand all combinations
  expanded <- apply(expand.grid(choices), 1, paste0, collapse = "")
  return(expanded)
}

# get expected vs observed for a given genome and motif

rm_motif_enrichment_strandwise <- function(genome, motif_seq) {
  
  # --------------------------------
  # Input validation
  # --------------------------------
  if (!is.character(genome) || length(genome) != 1 || !file.exists(genome)) {
    stop("genome must be a single FASTA filename")
  }
  
  # --------------------------------
  # Parse genome name + sequence
  # --------------------------------
  genome_name <- tools::file_path_sans_ext(basename(genome))
  genome_seq  <- readDNAStringSet(genome)[[1]]
  
  motif_seq <- toupper(motif_seq)
  rc_seq    <- as.character(reverseComplement(DNAString(motif_seq)))
  L <- nchar(motif_seq)
  
  if (L < 3) stop("Motif must be at least 3 bp")
  
  # Background frequencies
  tri_freq <- oligonucleotideFrequency(genome_seq, width = 3)
  tri_freq <- tri_freq / sum(tri_freq)
  
  di_freq <- oligonucleotideFrequency(genome_seq, width = 2)
  di_freq <- di_freq / sum(di_freq)
  
  # Expected probability for a single concrete motif
  expected_prob_markov2 <- function(motif) {
    nts <- strsplit(motif, "")[[1]]
    tris <- sapply(1:(L-2), function(i) paste0(nts[i:(i+2)], collapse=""))
    dins <- sapply(2:(L-2), function(i) paste0(nts[i:(i+1)], collapse=""))
    
    p <- unname(tri_freq[tris[1]])
    if(length(tris) > 1) {
      for(i in 2:length(tris)) {
        p <- p * (tri_freq[tris[i]] / di_freq[dins[i-1]])
      }
    }
    p
  }
  
  # Compute expected counts by expanding Ns and summing
  compute_expected <- function(motif) {
    concrete_seqs <- expand_IUPAC(motif)
    total_prob <- sum(sapply(concrete_seqs, expected_prob_markov2))
    total_prob * (length(genome_seq) - L + 1)
  }
  
  # Observed counts (wildcard N handled)
  obs_fwd <- countPattern(motif_seq, genome_seq, fixed = FALSE)
  obs_rev <- countPattern(rc_seq, genome_seq, fixed = FALSE)
  
  # Expected counts
  exp_fwd <- compute_expected(motif_seq)
  exp_rev <- compute_expected(rc_seq)
  
  # Output
  data.frame(
    motif      = motif_seq,
    strand     = c("forward", "reverse"),
    observed   = c(obs_fwd, obs_rev),
    expected   = c(exp_fwd, exp_rev),
    enrichment = c(obs_fwd / exp_fwd, obs_rev / exp_rev),
    genome = genome_name,
    genome_length = length(genome_seq),
    row.names  = NULL
  )
}

# 3.0 looping over all motifs and genomes ####


motifs <- read.csv("datafiles/motifs.csv")

# quick check to make sure no motifs are duplicated 
# motifs$motif <- as.factor(motifs$motif)
# 
# motifs <- motifs %>% select(motif)%>%
#   unique()

# filter out Type I RM for now as the high counts of Ns makes this extremely slow 
motifs <- motifs %>%
  filter(type != "RM_Type_I")

motifs <- as.vector(motifs$motif)

genomes <- list.files("datafiles/genomes",full.names = T)

data <- data.frame(
  motif = character(),
  strand = character(),
  observed = numeric(),
  expected = numeric(),
  enrichment = numeric(),
  genome = character(),
  genome_length = numeric(),
  stringsAsFactors = FALSE
)

startTime <- Sys.time()

for (i in 1:length(genomes)){
  for( y in 1:length(motifs)){
    d <- rm_motif_enrichment_strandwise(genome = genomes[i], motif_seq = motifs[y])
    data <- rbind(data,d)
  }
} 

endTime <- Sys.time()

print(endTime - startTime) 


