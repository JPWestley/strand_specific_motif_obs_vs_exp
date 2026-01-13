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
dict <- read.csv("datafiles/raw/phage_type_dict.csv")
motifs <- read.csv("datafiles/raw/motifs.csv")

data <- left_join(data,dict)

data <- left_join(data,motifs)

# 3.0 poisson testing ####

results <- data %>%
  rowwise() %>%
  mutate(
    pval = poisson.test(observed, T = expected, r = 1, alternative = "two.sided")$p.value,
  ) %>%
  ungroup()

results <- results %>%
  mutate(
    pval.adj = p.adjust(pval, method = "fdr"),
  )


results <- results %>%
  mutate(
    sig = case_when(
      pval.adj < 0.001 ~ "***",
      pval.adj < 0.01  ~ "**",
      pval.adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

results <- results %>%
  mutate(
    direction = case_when(
      pval.adj < 0.05 & enrichment < 1 ~ "Significantly lower",
      pval.adj < 0.05 & enrichment > 1 ~ "Significantly higher",
      TRUE ~ "Non-significant"
    )
  )

results <- results %>%
  mutate(normalised_ratio = (observed-expected)/(observed+expected))

results_wide <- results %>%
  pivot_wider(id_cols = c(motif,genome,motif_abundance,phage_type,genome_length,type,sub_type,Palindromic),
              names_from = strand,
              values_from = c(observed, expected, enrichment,pval,pval.adj,sig,direction,normalised_ratio))

results_wide <- results_wide %>%
  mutate(result = case_when(
    direction_forward == "Significantly lower" & direction_reverse == "Significantly lower" ~ "Double depleted",
    
    direction_forward == "Significantly higher" & direction_reverse == "Significantly higher" ~ "Double enriched",
    
    direction_forward == "Significantly lower" & direction_reverse == "Non-significant" ~ "Single depleted",
    direction_forward == "Non-significant" & direction_reverse == "Significantly lower" ~ "Single depleted",
    
    direction_forward == "Significantly higher" & direction_reverse == "Non-significant" ~ "Single enriched",
    direction_forward == "Non-significant" & direction_reverse == "Significantly higher" ~ "Single enriched",
    
    direction_forward == "Significantly lower" & direction_reverse == "Significantly higher" ~ "Depleted and enriched",
    direction_forward == "Significantly higher" & direction_reverse == "Significantly lower" ~ "Depleted and enriched",
    
    direction_forward == "Non-significant" & direction_reverse == "Non-significant" ~ "Non-significant"
    
  ))

str(results_wide)

results_wide$type <- as.factor(results_wide$type)
results_wide$result <- as.factor(results_wide$result)

results_wide$motif <- as.character(results_wide$motif)

results_wide <- results_wide %>%
  mutate(
    Palindromic = if_else(
      motif == as.character(reverseComplement(DNAStringSet(motif))),
      "Palindromic",
      "Non-palindromic"
    )
  )

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

p1 <- ggplot(results_wide,aes(x = normalised_ratio_forward, y = normalised_ratio_reverse,color = result))+
  geom_point(alpha = 0.15) +
  facet_wrap(~type)+
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept=0,linewidth = 0.5) +
  geom_vline(xintercept=0,linewidth = 0.5) +
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
  xlab("Ratio of observed to expected fwd strand (normalised)") +
  ylab("Ratio of observed to expected\n rvs strand (normalised)") +
  labs(colour = "Strand specific significance")+
  theme(legend.position = "bottom") +
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  )


p1

jpeg("plots/p1.jpeg", width = 2400, height = 2400, units = "px", res = 300)

p1

dev.off()


# facet by major types, no palindromes, use geom_point 

p2 <- ggplot(results_wide_NP,aes(x = normalised_ratio_forward, y = normalised_ratio_reverse,color = result))+
  geom_point(alpha = 0.2) +
  facet_wrap(~type)+
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept=0,linewidth = 0.5) +
  geom_vline(xintercept=0,linewidth = 0.5) +
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
  xlab("Ratio of observed to expected fwd strand (normalised)") +
  ylab("Ratio of observed to expected\n rvs strand (normalised)") +
  labs(colour = "Strand specific significance")+
  theme(legend.position = "bottom") +
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  )

jpeg("plots/p2.jpeg", width = 2400, height = 2400, units = "px", res = 300)

p2

dev.off()

# facet by major types, keep palindromes, use geom_pointdensity

p3 <- ggplot(results_wide,aes(x = normalised_ratio_forward, y = normalised_ratio_reverse))+
  geom_pointdensity(adjust= 0.1) +
  facet_wrap(~type)+
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept=0,linewidth = 0.5) +
  geom_vline(xintercept=0,linewidth = 0.5) +
  theme_bw() +
  scale_color_viridis(trans = "log10")+
  xlab("Ratio of observed to expected fwd strand (normalised)") +
  ylab("Ratio of observed to expected\n rvs strand (normalised)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom")

jpeg("plots/p3.jpeg", width = 2400, height = 2400, units = "px", res = 300)

p3

dev.off()

# facet by major types, no palindromes, use geom_pointdensity

p4 <- ggplot(results_wide_NP,aes(x = normalised_ratio_forward, y = normalised_ratio_reverse))+
  geom_pointdensity(adjust= 0.1) +
  facet_wrap(~type)+
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept=0,linewidth = 0.5) +
  geom_vline(xintercept=0,linewidth = 0.5) +
  theme_bw() +
  scale_color_viridis(trans = "log10")+
  xlab("Ratio of observed to expected fwd strand (normalised)") +
  ylab("Ratio of observed to expected\n rvs strand (normalised)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom")

jpeg("plots/p4.jpeg", width = 2400, height = 2400, units = "px", res = 300)

p4

dev.off()

# getting accompanying bar charts

# with palindromes

p5 <- ggplot(data=results_wide, aes(x=type, y = ,fill = result)) + 
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

jpeg("plots/p5.jpeg", width = 1800, height = 1200, units = "px", res = 300)

p5

dev.off()


# without palindromes

p6 <- ggplot(data=results_wide_NP, aes(x=type, y = ,fill = result)) + 
  geom_bar(position = "fill")+
  geom_bar(position = "fill") +
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

jpeg("plots/p6.jpeg", width = 1800, height = 1200, units = "px", res = 300)

p6

dev.off()

# by motif

results_wide$motif <- as.factor(results_wide$motif)

p7 <- ggplot(results_wide,aes(x = normalised_ratio_forward, y = normalised_ratio_reverse))+
  geom_pointdensity(adjust= 0.1) +
  facet_wrap(~motif,ncol = 8)+
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept=0,linewidth = 0.5) +
  geom_vline(xintercept=0,linewidth = 0.5) +
  theme_bw() +
  scale_color_viridis(trans = "log10")+
  xlab("Ratio of observed to expected fwd strand (normalised)") +
  ylab("Ratio of observed to expected\n rvs strand (normalised)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom")

jpeg("plots/p7.jpeg", width = 6000, height = 6000, units = "px", res = 300)

p7

dev.off()

results_wide$motif <- as.factor(results_wide$motif)

results_wide_PONLY <- results_wide %>%
  filter(Palindromic == "Palindromic")

p7.1 <- ggplot(results_wide_PONLY,aes(x = normalised_ratio_forward, y = normalised_ratio_reverse))+
  geom_pointdensity(adjust= 0.1) +
  facet_wrap(~motif,ncol = 8)+
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept=0,linewidth = 0.5) +
  geom_vline(xintercept=0,linewidth = 0.5) +
  theme_bw() +
  scale_color_viridis(trans = "log10")+
  xlab("Ratio of observed to expected fwd strand (normalised)") +
  ylab("Ratio of observed to expected\n rvs strand (normalised)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom")

jpeg("plots/p7.1.jpeg", width = 6000, height = 6000, units = "px", res = 300)

p7.1

dev.off()


# by phage

results_wide$genome <- as.factor(results_wide$genome)
results_wide$type <- as.factor(results_wide$type)

p8 <- ggplot(results_wide,aes(x = normalised_ratio_forward, y = normalised_ratio_reverse,colour = type))+
  geom_point(alpha = 0.3) +
  facet_wrap(~genome,ncol = 15)+
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept=0,linewidth = 0.5) +
  geom_vline(xintercept=0,linewidth = 0.5) +
  theme_bw() +
  scale_color_brewer(palette="Dark2") +
  xlab("Ratio of observed to expected fwd strand (normalised)") +
  ylab("Ratio of observed to expected\n rvs strand (normalised)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom")

jpeg("plots/p8.jpeg", width = 10000, height = 7000, units = "px", res = 300)

p8

dev.off()

# facet by RM sub-types, keep palindromes, use geom_point

results_wide_TII <- results_wide %>%
  filter(type == "RM_Type_II")

p9 <- ggplot(results_wide_TII,aes(x = normalised_ratio_forward, y = normalised_ratio_reverse,color = result))+
  geom_point(alpha = 0.2) +
  facet_wrap(~sub_type)+
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept=0,linewidth = 0.5) +
  geom_vline(xintercept=0,linewidth = 0.5) +
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
  xlab("Ratio of observed to expected fwd strand (normalised)") +
  ylab("Ratio of observed to expected\n rvs strand (normalised)") +
  labs(colour = "Strand specific significance")+
  theme(legend.position = "bottom") +
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  )

jpeg("plots/p9.jpeg", width = 2400, height = 2400, units = "px", res = 300)

p9

dev.off()

p10 <- ggplot(data=results_wide_TII, aes(x=sub_type, y = ,fill = result)) + 
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
  xlab("Type II RM sub-type") +
  ylab("Proportion") +
  theme_bw()

jpeg("plots/p10.jpeg", width = 2200, height = 1200, units = "px", res = 300)

p10

dev.off()


# separating out by phage type 

results_wide_lyt <- results_wide %>%
  filter(phage_type == "lytic")

results_wide_temp <- results_wide %>%
  filter(phage_type == "temperate")

p11 <- ggplot(results_wide_lyt,aes(x = normalised_ratio_forward, y = normalised_ratio_reverse,color = result))+
  geom_point(alpha = 0.15) +
  facet_wrap(~type)+
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept=0,linewidth = 0.5) +
  geom_vline(xintercept=0,linewidth = 0.5) +
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
  xlab("Ratio of observed to expected fwd strand (normalised)") +
  ylab("Ratio of observed to expected\n rvs strand (normalised)") +
  labs(colour = "Strand specific significance")+
  theme(legend.position = "bottom") +
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  )


p11

jpeg("plots/p11.jpeg", width = 2400, height = 2400, units = "px", res = 300)

p11

dev.off()

p12 <- ggplot(results_wide_temp,aes(x = normalised_ratio_forward, y = normalised_ratio_reverse,color = result))+
  geom_point(alpha = 0.15) +
  facet_wrap(~type)+
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept=0,linewidth = 0.5) +
  geom_vline(xintercept=0,linewidth = 0.5) +
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
  xlab("Ratio of observed to expected fwd strand (normalised)") +
  ylab("Ratio of observed to expected\n rvs strand (normalised)") +
  labs(colour = "Strand specific significance")+
  theme(legend.position = "bottom") +
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  )


p12

jpeg("plots/p12.jpeg", width = 2400, height = 2400, units = "px", res = 300)

p12

dev.off()

p13 <- ggplot(data=results_wide_lyt, aes(x=type, y = ,fill = result)) + 
  geom_bar(position = "fill")+
  geom_bar(position = "fill") +
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

jpeg("plots/p13.jpeg", width = 2400, height = 1200, units = "px", res = 300)

p13

dev.off()

p14 <- ggplot(data=results_wide_temp, aes(x=type, y = ,fill = result)) + 
  geom_bar(position = "fill")+
  geom_bar(position = "fill") +
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

jpeg("plots/p14.jpeg", width = 2400, height = 1200, units = "px", res = 300)

p14

dev.off()

results_wide$motif <- as.factor(results_wide$motif)

results_wide_TI <- results_wide %>%
  filter(type == "RM_Type_I")

results_wide_TII <- results_wide %>%
  filter(type == "RM_Type_II")

results_wide_TIII <- results_wide %>%
  filter(type == "RM_Type_III")

results_wide_PAMS <- results_wide %>%
  filter(type == "PAMS")


p15 <- ggplot(results_wide_TI,aes(x = normalised_ratio_forward, y = normalised_ratio_reverse))+
  geom_pointdensity(adjust= 0.1) +
  facet_wrap(~motif,ncol = 8)+
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept=0,linewidth = 0.5) +
  geom_vline(xintercept=0,linewidth = 0.5) +
  theme_bw() +
  scale_color_viridis(trans = "log10")+
  xlab("Ratio of observed to expected fwd strand (normalised)") +
  ylab("Ratio of observed to expected\n rvs strand (normalised)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom")

jpeg("plots/p15.jpeg", width = 4000, height = 3000, units = "px", res = 300)

p15

dev.off()

p16 <- ggplot(results_wide_TII,aes(x = normalised_ratio_forward, y = normalised_ratio_reverse, color = sub_type))+
  geom_pointdensity(adjust= 0.1) +
  facet_wrap(~motif,ncol = 6)+
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept=0,linewidth = 0.5) +
  geom_vline(xintercept=0,linewidth = 0.5) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2")+
  xlab("Ratio of observed to expected fwd strand (normalised)") +
  ylab("Ratio of observed to expected\n rvs strand (normalised)") +
  #labs(colour = "Density")+
  theme(legend.position = "bottom")

jpeg("plots/p16.jpeg", width = 3000, height = 3000, units = "px", res = 300)

p16

dev.off()

p17 <- ggplot(results_wide_TIII,aes(x = normalised_ratio_forward, y = normalised_ratio_reverse))+
  geom_pointdensity(adjust= 0.1) +
  facet_wrap(~motif,ncol = 4)+
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept=0,linewidth = 0.5) +
  geom_vline(xintercept=0,linewidth = 0.5) +
  theme_bw() +
  scale_color_viridis(trans = "log10")+
  xlab("Ratio of observed to expected fwd strand (normalised)") +
  ylab("Ratio of observed to expected\n rvs strand (normalised)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom")

jpeg("plots/p17.jpeg", width = 3000, height = 2500, units = "px", res = 300)

p17

dev.off()

p18 <- ggplot(results_wide_PAMS,aes(x = normalised_ratio_forward, y = normalised_ratio_reverse))+
  geom_pointdensity(adjust= 0.1) +
  facet_wrap(~motif,ncol = 2)+
  xlim(-1,1) +
  ylim(-1,1) +
  geom_hline(yintercept=0,linewidth = 0.5) +
  geom_vline(xintercept=0,linewidth = 0.5) +
  theme_bw() +
  scale_color_viridis(trans = "log10")+
  xlab("Ratio of observed to expected fwd strand (normalised)") +
  ylab("Ratio of observed to expected\n rvs strand (normalised)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom")

jpeg("plots/p18.jpeg", width = 2000, height = 1300, units = "px", res = 300)

p18

dev.off()




