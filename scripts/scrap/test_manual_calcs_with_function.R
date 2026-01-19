library(Biostrings)
library(tidyverse)
#install.packages("tidyverse")

# # 1.0 Using trinucleotide frequencies with overlap correction to predict expected counts ####
# 
# rm_motif_enrichment_strandwise <- function(genome_seq, motif_seq) {
#   
#   genome_seq <- DNAString(genome_seq)
#   motif_seq  <- toupper(motif_seq)
#   rc_seq     <- as.character(reverseComplement(DNAString(motif_seq)))
#   
#   L <- nchar(motif_seq)
#   if (L < 3) stop("Motif must be at least 3 bp")
#   
#   # --------------------------------
#   # Background frequencies
#   # --------------------------------
#   tri_freq <- oligonucleotideFrequency(genome_seq, width = 3)
#   tri_freq <- tri_freq / sum(tri_freq)
#   
#   di_freq <- oligonucleotideFrequency(genome_seq, width = 2)
#   di_freq <- di_freq / sum(di_freq)
#   
#   # --------------------------------
#   # Expected probability (order-2 Markov)
#   # --------------------------------
#   expected_prob_markov2 <- function(motif) {
#     nts <- strsplit(motif, "")[[1]]
#     
#     tris <- sapply(1:(L - 2), function(i)
#       paste0(nts[i:(i + 2)], collapse = "")
#     )
#     
#     dins <- sapply(2:(L - 2), function(i)
#       paste0(nts[i:(i + 1)], collapse = "")
#     )
#     
#     p <- unname(tri_freq[tris[1]])
#     
#     for (i in 2:length(tris)) {
#       p <- p * (tri_freq[tris[i]] / di_freq[dins[i - 1]])
#     }
#     
#     p
#   }
#   
#   # --------------------------------
#   # Observed and expected
#   # --------------------------------
#   obs_fwd <- countPattern(motif_seq, genome_seq)
#   exp_fwd <- expected_prob_markov2(motif_seq) * (length(genome_seq) - L + 1)
#   
#   obs_rev <- countPattern(rc_seq, genome_seq)
#   exp_rev <- expected_prob_markov2(rc_seq) * (length(genome_seq) - L + 1)
#   
#   # --------------------------------
#   # Output (same motif name)
#   # --------------------------------
#   data.frame(
#     motif      = motif_seq,
#     strand     = c("forward", "reverse"),
#     observed   = c(obs_fwd, obs_rev),
#     expected   = c(exp_fwd, exp_rev),
#     enrichment = c(obs_fwd / exp_fwd, obs_rev / exp_rev),
#     genome = genome,
#     row.names = NULL
#   )
# }
# 
# genome <- readDNAStringSet("datafiles/genomes/141.fasta")[[1]]
# 
# data <- rm_motif_enrichment_strandwise(genome, "CAACAC")
# 
# # 2.0 Using only single nucleotide frequencies to get expected counts ####
# 
# library(Biostrings)
# 
# rm_motif_enrichment_strandwise_mono <- function(genome_seq, motif_seq) {
#   
#   genome_seq <- DNAString(genome_seq)
#   motif_seq  <- toupper(motif_seq)
#   rc_seq     <- as.character(reverseComplement(DNAString(motif_seq)))
#   
#   L <- nchar(motif_seq)
#   if (L < 1) stop("Motif must be at least 1 bp")
#   
#   # --------------------------------
#   # Background frequencies (mono)
#   # --------------------------------
#   mono_freq <- letterFrequency(genome_seq, letters = c("A","C","G","T"))
#   mono_freq <- mono_freq / sum(mono_freq)
#   
#   # --------------------------------
#   # Expected probability (iid model)
#   # --------------------------------
#   expected_prob_mono <- function(motif) {
#     nts <- strsplit(motif, "")[[1]]
#     prod(mono_freq[nts])
#   }
#   
#   # --------------------------------
#   # Observed and expected
#   # --------------------------------
#   obs_fwd <- countPattern(motif_seq, genome_seq)
#   exp_fwd <- expected_prob_mono(motif_seq) * (length(genome_seq) - L + 1)
#   
#   obs_rev <- countPattern(rc_seq, genome_seq)
#   exp_rev <- expected_prob_mono(rc_seq) * (length(genome_seq) - L + 1)
#   
#   # --------------------------------
#   # Output
#   # --------------------------------
#   data.frame(
#     motif      = motif_seq,
#     strand     = c("forward", "reverse"),
#     observed   = c(obs_fwd, obs_rev),
#     expected   = c(exp_fwd, exp_rev),
#     enrichment = c(obs_fwd / exp_fwd, obs_rev / exp_rev),
#     row.names  = NULL
#   )
# }
# 
# genome <- readDNAStringSet("datafiles/genomes/141.fasta")[[1]]
# 
# data <- rm_motif_enrichment_strandwise_mono(genome, "CAACAC")
# 
# # 3.0 Now I configure the function to interpret Ns in the sequences ####
# library(Biostrings)
# 
# rm_motif_enrichment_strandwise_window <- function(genome_seq, motif_seq) {
#   
#   genome_seq <- DNAString(genome_seq)
#   motif_seq  <- toupper(motif_seq)
#   rc_seq     <- as.character(reverseComplement(DNAString(motif_seq)))
#   
#   L <- nchar(motif_seq)
#   if (L < 3) stop("Motif must be at least 3 bp")
#   
#   # -----------------------------
#   # Background frequencies
#   # -----------------------------
#   tri_freq <- oligonucleotideFrequency(genome_seq, width = 3)
#   tri_freq <- tri_freq / sum(tri_freq)
#   
#   di_freq <- oligonucleotideFrequency(genome_seq, width = 2)
#   di_freq <- di_freq / sum(di_freq)
#   
#   # -----------------------------
#   # Helper: probability of matching motif for a given window
#   # Handles Ns properly
#   # -----------------------------
#   window_prob <- function(window, motif) {
#     # Split into nucleotides
#     win_nts <- strsplit(as.character(window), "")[[1]]
#     motif_nts <- strsplit(motif, "")[[1]]
#     
#     # Overlapping trinucleotides in motif
#     motif_tris <- sapply(1:(L-2), function(i) paste0(motif_nts[i:(i+2)], collapse=""))
#     motif_dins <- sapply(2:(L-2), function(i) paste0(motif_nts[i:(i+1)], collapse=""))
#     
#     # Probability of this window
#     p <- 1
#     
#     for(i in 1:(L-2)) {
#       # Determine trinucleotide in window
#       tri_win <- paste0(win_nts[i:(i+2)], collapse="")
#       tri_motif <- motif_tris[i]
#       
#       # Compute contribution for each base in the motif trinucleotide
#       # If motif base is N, treat as wildcard (average over all possibilities)
#       prob_tri <- 1
#       for(j in 1:3) {
#         base <- substr(tri_motif, j, j)
#         if(base == "N") {
#           prob_tri <- prob_tri * 1  # keep 1, contribution is handled by averaging later
#         } else if(substr(tri_win, j, j) != base) {
#           prob_tri <- 0  # mismatch
#         }
#       }
#       
#       # If i==1, no division by di
#       if(i == 1) {
#         p <- p * tri_freq[tri_win]
#       } else {
#         di_win <- paste0(win_nts[i:(i+1)], collapse="")
#         p <- p * (tri_freq[tri_win] / di_freq[di_win])
#       }
#       
#       # If tri contains N, we can marginalize over all possibilities
#       # For simplicity, since wildcard treated as 1 in contribution, leave as is
#     }
#     return(p)
#   }
#   
#   # -----------------------------
#   # Expected count by sliding window
#   # -----------------------------
#   expected_count_window <- function(motif) {
#     windows <- DNAStringSet(Views(genome_seq, start=1:(length(genome_seq)-L+1), width=L))
#     probs <- sapply(windows, window_prob, motif=motif)
#     sum(probs)
#   }
#   
#   # -----------------------------
#   # Observed counts
#   # -----------------------------
#   obs_fwd <- countPattern(motif_seq, genome_seq, fixed=FALSE)
#   obs_rev <- countPattern(rc_seq, genome_seq, fixed=FALSE)
#   
#   # -----------------------------
#   # Expected counts
#   # -----------------------------
#   exp_fwd <- expected_count_window(motif_seq)
#   exp_rev <- expected_count_window(rc_seq)
#   
#   # -----------------------------
#   # Output
#   # -----------------------------
#   data.frame(
#     motif      = motif_seq,
#     strand     = c("forward", "reverse"),
#     observed   = c(obs_fwd, obs_rev),
#     expected   = c(exp_fwd, exp_rev),
#     enrichment = c(obs_fwd / exp_fwd, obs_rev / exp_rev),
#     row.names  = NULL
#   )
# }
# 
# genome <- readDNAStringSet("datafiles/genomes/141.fasta")[[1]]
# 
# data2 <- rm_motif_enrichment_strandwise_Ns(genome, "CANNAC")

# 4.0 incorporating Ns full manual expansion THIS WORKS WELL ####
# 
# library(Biostrings)
# 
# expand_IUPAC <- function(motif, max_ambiguity = 3) {
#   # Split motif into letters
#   nts <- strsplit(toupper(motif), "")[[1]]
#   
#   # Map each letter to all possible bases using IUPAC_CODE_MAP
#   choices <- lapply(nts, function(b) {
#     if(!b %in% names(IUPAC_CODE_MAP)) stop(paste("Unknown nucleotide code:", b))
#     strsplit(IUPAC_CODE_MAP[[b]], "")[[1]]
#   })
#   
#   # Count total ambiguity positions
#   n_ambiguity <- sum(sapply(choices, length) > 1)
#   if(n_ambiguity > max_ambiguity) {
#     stop(paste0("Too many ambiguous positions (", n_ambiguity, 
#                 "), reduce motif ambiguity or increase max_ambiguity"))
#   }
#   
#   # Expand all combinations
#   expanded <- apply(expand.grid(choices), 1, paste0, collapse = "")
#   return(expanded)
# }
# 
# # Strand-specific RM motif function, adding expected over all concrete sequences
# rm_motif_enrichment_strandwise_Ns_expand <- function(genome_seq, motif_seq) {
#   
#   genome_seq <- DNAString(genome_seq)
#   motif_seq  <- toupper(motif_seq)
#   rc_seq     <- as.character(reverseComplement(DNAString(motif_seq)))
#   L <- nchar(motif_seq)
#   
#   # Background frequencies
#   tri_freq <- oligonucleotideFrequency(genome_seq, width = 3)
#   tri_freq <- tri_freq / sum(tri_freq)
#   
#   di_freq <- oligonucleotideFrequency(genome_seq, width = 2)
#   di_freq <- di_freq / sum(di_freq)
#   
#   # Expected probability for a single concrete motif
#   expected_prob_markov2 <- function(motif) {
#     nts <- strsplit(motif, "")[[1]]
#     tris <- sapply(1:(L-2), function(i) paste0(nts[i:(i+2)], collapse=""))
#     dins <- sapply(2:(L-2), function(i) paste0(nts[i:(i+1)], collapse=""))
#     
#     p <- unname(tri_freq[tris[1]])
#     if(length(tris) > 1) {
#       for(i in 2:length(tris)) {
#         p <- p * (tri_freq[tris[i]] / di_freq[dins[i-1]])
#       }
#     }
#     p
#   }
#   
#   # Compute expected counts by expanding Ns and summing
#   compute_expected <- function(motif) {
#     concrete_seqs <- expand_IUPAC(motif)
#     total_prob <- sum(sapply(concrete_seqs, expected_prob_markov2))
#     total_prob * (length(genome_seq) - L + 1)
#   }
#   
#   # Observed counts (wildcard N handled)
#   obs_fwd <- countPattern(motif_seq, genome_seq, fixed = FALSE)
#   obs_rev <- countPattern(rc_seq, genome_seq, fixed = FALSE)
#   
#   # Expected counts
#   exp_fwd <- compute_expected(motif_seq)
#   exp_rev <- compute_expected(rc_seq)
#   
#   # Output
#   data.frame(
#     motif      = motif_seq,
#     strand     = c("forward", "reverse"),
#     observed   = c(obs_fwd, obs_rev),
#     expected   = c(exp_fwd, exp_rev),
#     enrichment = c(obs_fwd / exp_fwd, obs_rev / exp_rev),
#     row.names  = NULL
#   )
# }
# 
# genome <- readDNAStringSet("datafiles/genomes/141.fasta")[[1]]
# 
# data3 <- rm_motif_enrichment_strandwise_Ns_expand(genome, "CAACAC")

# 5.0 now adding genome name to dataframe ####

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

motifs <- read.csv("datafiles/motifs.csv")

# motifs$motif <- as.factor(motifs$motif)
# 
# motifs <- motifs %>% select(motif)%>%
#   unique()

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


