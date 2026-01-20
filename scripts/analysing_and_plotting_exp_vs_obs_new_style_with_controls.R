#############################
### ANALYSIS AND PLOTTING ###
#############################

# 1.0 load dependencies ####

library(tidyverse)
library(RColorBrewer)
library(ggpointdensity)
library(Biostrings)
library(viridis)

#install.packages("viridis")

# 2.0 loading and joining dataframes ####

data <- read.csv("datafiles/processed/obvs_vs_exp.csv")

datac <- read.csv("datafiles/processed/obvs_vs_exp_control_motifs.csv")

dict <- read.csv("datafiles/raw/phage_type_dict.csv")
motifs <- read.csv("datafiles/raw/motifs.csv")

data <- left_join(data,dict)
datac <- left_join(datac,dict)

# datac <- datac %>%
#   unique()
# 

data <- left_join(data,motifs)

data <- data %>%
  mutate(
    Palindromic = if_else(
      motif == as.character(reverseComplement(DNAStringSet(motif))),
      "Palindromic",
      "Non-palindromic"
    )
  )

datac$subtype <- NULL

datac <- datac %>%
  mutate(type = "control",
         sub_type = "control",
         motif_abundance = "1",
         Palindromic = if_else(motif == as.character(reverseComplement(DNAStringSet(motif))),"Palindromic","Non-palindromic")
  )

data2 <- rbind(data,datac)

# 3.0 using standardised residuals ####

results <- data2 %>%
  mutate(
    SR = (observed - expected) / sqrt(expected)) %>%
  mutate(signed_log_SR = sign(SR) * log10(abs(SR) + 1))

results_wide <- results %>%
  pivot_wider(id_cols = c(motif,genome,motif_abundance,phage_type,genome_length,type,sub_type,Palindromic),
              names_from = strand,
              values_from = c(observed, expected, enrichment,signed_log_SR,SR))

threshold <- log10(1.96 + 1)  # 95% cutoff on signed log10 SR scale (need cutoff to be 2.5 to aproximate same levels of significance as FDR corrected)

results_wide <- results_wide %>%
  mutate(
    result = case_when(
      signed_log_SR_forward <= -threshold & signed_log_SR_reverse <= -threshold ~ "Double depleted",
      signed_log_SR_forward >=  threshold & signed_log_SR_reverse >=  threshold ~ "Double enriched",
      
      (signed_log_SR_forward <= -threshold & abs(signed_log_SR_reverse) < threshold) |
        (abs(signed_log_SR_forward) < threshold & signed_log_SR_reverse <= -threshold) ~ "Single depleted",
      
      (signed_log_SR_forward >=  threshold & abs(signed_log_SR_reverse) < threshold) |
        (abs(signed_log_SR_forward) < threshold & signed_log_SR_reverse >=  threshold) ~ "Single enriched",
      
      (signed_log_SR_forward <= -threshold & signed_log_SR_reverse >=  threshold) |
        (signed_log_SR_forward >=  threshold & signed_log_SR_reverse <= -threshold) ~ "Depleted and enriched",
      
      abs(signed_log_SR_forward) < threshold & abs(signed_log_SR_reverse) < threshold ~ "Non-significant"
    )
  )


str(results_wide)

results_wide$type <- as.factor(results_wide$type)
# results_wide$result <- as.factor(results_wide$result)

results_wide$motif <- as.character(results_wide$motif)

results_wide$result <- factor(
  results_wide$result,
  levels = c(
    "Non-significant",
    "Single depleted",
    "Double depleted",
    "Single enriched",
    "Double enriched",
    "Depleted and enriched"
  )
)

results_wide_NP <- results_wide %>%
  filter(Palindromic == "Non-palindromic")

# 4.0 Plotting ####

# facet by major types, keep palindromes, use geom_point 

str(results_wide)

# Compute 2D bins per facet
library(ggplot2)
library(dplyr)

# Precompute counts per facet
binned <- results_wide %>%
  group_by(type) %>%
  mutate(
    # Compute bin index along x and y
    x_bin = cut(signed_log_SR_forward, breaks = seq(-1, 1.3, length.out = 41), include.lowest = TRUE),
    y_bin = cut(signed_log_SR_reverse, breaks = seq(-1, 1.3, length.out = 41), include.lowest = TRUE)
  ) %>%
  count(type, x_bin, y_bin) %>%
  group_by(type) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  # Convert bins to numeric midpoints for plotting
  mutate(
    x_mid = as.numeric(sub("\\((.+),.*", "\\1", x_bin)) + diff(range(-1,1.3))/20/2,
    y_mid = as.numeric(sub("\\((.+),.*", "\\1", y_bin)) + diff(range(-1,1.3))/20/2
  )

binned <- binned %>%
  filter(type != "control")

P <- ggplot(binned, aes(x = x_mid, y = y_mid, fill = log10(prop + 1e-6))) +
  geom_tile(width = diff(range(-1,1.3))/20, height = diff(range(-1,1.3))/20) +
  facet_wrap(~type) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = c(- 0.553,  0.553)) +
  geom_vline(xintercept = c(- 0.553,  0.553)) +
  theme_bw() +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(fill = "log10(Proportion)") +
  scale_fill_viridis_c() +
  theme(legend.position = "bottom")

P

jpeg("plots/p1C.jpeg", width = 2400, height = 2400, units = "px", res = 300)

P

dev.off()


p1B <- ggplot(results_wide, aes(x = signed_log_SR_forward, y = signed_log_SR_reverse)) +
  geom_bin2d(
    aes(fill = after_stat(log10((count / sum(count)) + 1e-6))),
    bins = 40
  ) +
  facet_wrap(~type) +
  coord_cartesian(xlim = c(-1, 1.3), ylim = c(-1, 1.3)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = c(-0.471, 0.471)) +
  geom_vline(xintercept = c(-0.471, 0.471)) +
  theme_bw() +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(fill = "log10(Proportion)") +
  scale_fill_viridis_c() +
  theme(legend.position = "bottom")





jpeg("plots/p1B.jpeg", width = 2400, height = 2400, units = "px", res = 300)

p1B

dev.off()

# p1B <- ggplot(results_wide,aes(x = signed_log_SR_forward, y = signed_log_SR_reverse))+
#   geom_pointdensity(adjust= 0.1,size = 0.7) +
#   facet_wrap(~type)+
#   xlim(-1,1.3) +
#   ylim(-1,1.3) +
#   geom_hline(yintercept=0,linewidth = 0.5,linetype = "dashed") +
#   geom_vline(xintercept=0,linewidth = 0.5,linetype = "dashed") +
#   geom_hline(yintercept = c(-0.471, 0.471))+
#   geom_vline(xintercept = c(-0.471, 0.471))+
#   theme_bw() +
#   xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
#   ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
#   labs(colour = "Density")+
#   theme(legend.position = "bottom") +
#   scale_color_viridis(trans = "log10")

jpeg("plots/p1B.jpeg", width = 2400, height = 2400, units = "px", res = 300)

p1B

dev.off()

p2B <- ggplot(results_wide,aes(x = signed_log_SR_forward, y = signed_log_SR_reverse,colour = result))+
  geom_point(size = 1) +
  facet_wrap(~type)+
  xlim(-1,1.3) +
  ylim(-1,1.3) +
  geom_hline(yintercept=0,linewidth = 0.5,linetype = "dashed") +
  geom_vline(xintercept=0,linewidth = 0.5,linetype = "dashed") +
  geom_hline(yintercept = c(-0.471, 0.471))+
  geom_vline(xintercept = c(-0.471, 0.471))+
  theme_bw() +
  scale_color_manual(
    values = c(
      "Non-significant"        = "grey",
      "Depleted and enriched"  = "#1B9E77",
      "Double depleted"        = "#D95F02",
      "Double enriched"        = "#7570B3",
      "Single depleted"        = "#E7298A",
      "Single enriched"        = "#66A61E"
    )
  ) +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom")+
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  )

jpeg("plots/p2B.jpeg", width = 2400, height = 2400, units = "px", res = 300)

p2B

dev.off()

p3B <- ggplot(data=results_wide, aes(x=type, y = ,fill = result)) + 
  geom_bar(position = "fill")+
  scale_fill_manual(
    values = c(
      "Non-significant"        = "grey",
      "Depleted and enriched"  = "#1B9E77",
      "Double depleted"        = "#D95F02",
      "Double enriched"        = "#7570B3",
      "Single depleted"        = "#E7298A",
      "Single enriched"        = "#66A61E"
    )
  ) +
  labs(fill = "Strand specific significance")+
  xlab("Defence system type") +
  ylab("Proportion") +
  theme_bw()

jpeg("plots/p3B.jpeg", width = 1800, height = 1200, units = "px", res = 300)

p3B

dev.off()

results_wide_TII <- results_wide %>%
  filter(type == "RM_Type_II")

p4B <- ggplot(results_wide_TII,aes(x = signed_log_SR_forward, y = signed_log_SR_reverse,colour = result))+
  geom_point(size = 1) +
  facet_wrap(~sub_type)+
  xlim(-1,1.3) +
  ylim(-1,1.3) +
  geom_hline(yintercept=0,linewidth = 0.5,linetype = "dashed") +
  geom_vline(xintercept=0,linewidth = 0.5,linetype = "dashed") +
  geom_hline(yintercept = c(-0.471, 0.471))+
  geom_vline(xintercept = c(-0.471, 0.471))+
  theme_bw() +
  scale_color_manual(
    values = c(
      "Non-significant"        = "grey",
      "Depleted and enriched"  = "#1B9E77",
      "Double depleted"        = "#D95F02",
      "Double enriched"        = "#7570B3",
      "Single depleted"        = "#E7298A",
      "Single enriched"        = "#66A61E"
    )
  ) +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom")+
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  )

jpeg("plots/p4B.jpeg", width = 2400, height = 2400, units = "px", res = 300)

p4B

dev.off()

p5B <- ggplot(data=results_wide_TII, aes(x=sub_type, y = ,fill = result)) + 
  geom_bar(position = "fill")+
  scale_fill_manual(
    values = c(
      "Non-significant"        = "grey",
      "Depleted and enriched"  = "#1B9E77",
      "Double depleted"        = "#D95F02",
      "Double enriched"        = "#7570B3",
      "Single depleted"        = "#E7298A",
      "Single enriched"        = "#66A61E"
    )
  ) +
  labs(fill = "Strand specific significance")+
  xlab("RM Type II subtype") +
  ylab("Proportion") +
  theme_bw()

jpeg("plots/p5B.jpeg", width = 1800, height = 1200, units = "px", res = 300)

p5B

dev.off()


results_wide_TI <- results_wide %>%
  filter(type == "RM_Type_I")

results_wide_TII <- results_wide %>%
  filter(type == "RM_Type_II")

results_wide_TIII <- results_wide %>%
  filter(type == "RM_Type_III")

results_wide_PAMS <- results_wide %>%
  filter(type == "PAMS")


p6B <- ggplot(results_wide_TI,aes(x = signed_log_SR_forward, y = signed_log_SR_reverse,colour = result))+
  geom_point(size = 2) +
  facet_wrap(~motif)+
  xlim(-1,1.3) +
  ylim(-1,1.3) +
  geom_hline(yintercept=0,linewidth = 0.5,linetype = "dashed") +
  geom_vline(xintercept=0,linewidth = 0.5,linetype = "dashed") +
  geom_hline(yintercept = c(-0.471, 0.471))+
  geom_vline(xintercept = c(-0.471, 0.471))+
  theme_bw() +
  scale_color_manual(
    values = c(
      "Non-significant"        = "grey",
      "Depleted and enriched"  = "#1B9E77",
      "Double depleted"        = "#D95F02",
      "Double enriched"        = "#7570B3",
      "Single depleted"        = "#E7298A",
      "Single enriched"        = "#66A61E"
    )
  ) +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom")+
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  )

jpeg("plots/p6B.jpeg", width = 6000, height = 6000, units = "px", res = 300)

p6B

dev.off()


p7B <- ggplot(results_wide_TII,aes(x = signed_log_SR_forward, y = signed_log_SR_reverse,colour = result))+
  geom_point(size = 1) +
  facet_wrap(~motif)+
  xlim(-1,1.3) +
  ylim(-1,1.3) +
  geom_hline(yintercept=0,linewidth = 0.5,linetype = "dashed") +
  geom_vline(xintercept=0,linewidth = 0.5,linetype = "dashed") +
  geom_hline(yintercept = c(-0.471, 0.471))+
  geom_vline(xintercept = c(-0.471, 0.471))+
  theme_bw() +
  scale_color_manual(
    values = c(
      "Non-significant"        = "grey",
      "Depleted and enriched"  = "#1B9E77",
      "Double depleted"        = "#D95F02",
      "Double enriched"        = "#7570B3",
      "Single depleted"        = "#E7298A",
      "Single enriched"        = "#66A61E"
    )
  ) +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom")+
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  )

jpeg("plots/p7B.jpeg", width = 3000, height = 3000, units = "px", res = 300)

p7B

dev.off()


p8B <- ggplot(results_wide_TIII,aes(x = signed_log_SR_forward, y = signed_log_SR_reverse,colour = result))+
  geom_point(size = 1) +
  facet_wrap(~motif)+
  xlim(-1,1.3) +
  ylim(-1,1.3) +
  geom_hline(yintercept=0,linewidth = 0.5,linetype = "dashed") +
  geom_vline(xintercept=0,linewidth = 0.5,linetype = "dashed") +
  geom_hline(yintercept = c(-0.471, 0.471))+
  geom_vline(xintercept = c(-0.471, 0.471))+
  theme_bw() +
  scale_color_manual(
    values = c(
      "Non-significant"        = "grey",
      "Depleted and enriched"  = "#1B9E77",
      "Double depleted"        = "#D95F02",
      "Double enriched"        = "#7570B3",
      "Single depleted"        = "#E7298A",
      "Single enriched"        = "#66A61E"
    )
  ) +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom")+
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  )

jpeg("plots/p8B.jpeg", width = 3000, height = 2500, units = "px", res = 300)

p8B

dev.off()


p9B <- ggplot(results_wide_PAMS,aes(x = signed_log_SR_forward, y = signed_log_SR_reverse,colour = result))+
  geom_point(size = 1) +
  facet_wrap(~motif)+
  xlim(-1,1.3) +
  ylim(-1,1.3) +
  geom_hline(yintercept=0,linewidth = 0.5,linetype = "dashed") +
  geom_vline(xintercept=0,linewidth = 0.5,linetype = "dashed") +
  geom_hline(yintercept = c(-0.471, 0.471))+
  geom_vline(xintercept = c(-0.471, 0.471))+
  theme_bw() +
  scale_color_manual(
    values = c(
      "Non-significant"        = "grey",
      "Depleted and enriched"  = "#1B9E77",
      "Double depleted"        = "#D95F02",
      "Double enriched"        = "#7570B3",
      "Single depleted"        = "#E7298A",
      "Single enriched"        = "#66A61E"
    )
  ) +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom")+
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  )

jpeg("plots/p9B.jpeg", width = 2000, height = 1300, units = "px", res = 300)

p9B

dev.off()

