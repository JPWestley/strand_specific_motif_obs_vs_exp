#############################
### ANALYSIS AND PLOTTING ###
#############################

# 1.0 load dependencies ####

library(tidyverse)
library(RColorBrewer)
library(ggpointdensity)
library(Biostrings)
library(viridis)
library(ggpubr)
library(patchwork)
library(cowplot)

# 2.0 loading and joining dataframes ####

data <- read.csv("datafiles/processed/obvs_vs_exp.csv")
dict <- read.csv("datafiles/raw/phage_type_dict.csv")
motifs <- read.csv("datafiles/raw/motifs.csv")

data <- data %>%
  mutate(genome = str_replace(genome, "^CPL0*", "P")) 

data <- left_join(data,dict)

data <- left_join(data,motifs)

# 3.0 using standardised residuals ####

results <- data %>%
  mutate(
   SR = (observed - expected) / sqrt(expected)) %>%
  mutate(signed_log_SR = sign(SR) * log10(abs(SR) + 1))

# results <- results %>%
#   mutate(
#     sig = case_when(
#       pval.adj < 0.001 ~ "***",
#       pval.adj < 0.01  ~ "**",
#       pval.adj < 0.05  ~ "*",
#       TRUE ~ ""
#     )
#   )

# results <- results %>%
#   mutate(
#     direction = case_when(
#       pval.adj < 0.05 & enrichment < 1 ~ "Significantly lower",
#       pval.adj < 0.05 & enrichment > 1 ~ "Significantly higher",
#       TRUE ~ "Non-significant"
#     )
#   )
# 
# results <- results %>%
#   mutate(normalised_ratio = (observed-expected)/(observed+expected))

results_wide <- results %>%
  pivot_wider(id_cols = c(motif,genome,motif_abundance,phage_type,genome_length,type,sub_type,Palindromic),
              names_from = strand,
              values_from = c(observed, expected, enrichment,signed_log_SR,SR))

threshold <- log10(2.576 + 1)  # 95% cutoff on signed log10 SR scale (need cutoff to be 2.5 to aproximate same levels of significance as FDR corrected)

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

# results_wide_NP <- results_wide %>%
#   filter(Palindromic == "Non-palindromic")

results_wide <- results_wide %>%
  mutate(sub_type2 = case_when(
    type == "RM_Type_I" & Palindromic == "Non-palindromic" ~ "Non-palindromic\nType I_RM",
    type == "RM_Type_I" & Palindromic == "Palindromic" ~ "Palindromic\nType I RM",
    sub_type == "RM_Type_IIG" ~ "RM Type IIG & IIS",
    sub_type == "RM_Type_IIS" ~ "RM Type IIG & IIS",
    sub_type == "RM_Type_IIP" ~ "RM Type IIP",
    sub_type == "RM_Type_III" ~ "RM Type III",
    TRUE ~ sub_type)) %>%
  filter(sub_type != "RM_Type_IIM")
  

results_wide$sub_type2 <- as.factor(results_wide$sub_type2)

results_wide$sub_type2 <- factor(
  results_wide$sub_type2,
  levels = c(
    "Non-palindromic\nType I_RM",
    "Palindromic\nType I RM",
    "RM Type IIG & IIS",
    "RM Type IIP",
    "RM Type III",
    "PAMS"
  ))

## 3.1 producing all sub-dataframes for plotting ####

results_wide_lyt <- results_wide %>%
  filter(phage_type == "lytic")

results_wide_temp <- results_wide %>%
  filter(phage_type == "temperate")

results_wide_lyt_NPTI <- results_wide_lyt %>% filter(sub_type2 == "Non-palindromic\nType I_RM")
results_wide_lyt_PTI <- results_wide_lyt %>% filter(sub_type2 == "Palindromic\nType I RM")
results_wide_lyt_NPTII <- results_wide_lyt %>% filter(sub_type2 == "RM Type IIG & IIS")
results_wide_lyt_PTII <- results_wide_lyt %>% filter(sub_type2 == "RM Type IIP")
results_wide_lyt_TIII <- results_wide_lyt %>% filter(sub_type2 == "RM Type III")
results_wide_lyt_PAMS <- results_wide_lyt %>% filter(sub_type2 == "PAMS")

results_wide_temp_NPTI <- results_wide_temp %>% filter(sub_type2 == "Non-palindromic\nType I_RM")
results_wide_temp_PTI <- results_wide_temp %>% filter(sub_type2 == "Palindromic\nType I RM")
results_wide_temp_NPTII <- results_wide_temp %>% filter(sub_type2 == "RM Type IIG & IIS")
results_wide_temp_PTII <- results_wide_temp %>% filter(sub_type2 == "RM Type IIP")
results_wide_temp_TIII <- results_wide_temp %>% filter(sub_type2 == "RM Type III")
results_wide_temp_PAMS <- results_wide_temp %>% filter(sub_type2 == "PAMS")

bins <- 30

binned_lyt <- results_wide_lyt %>%
  group_by(sub_type2) %>%
  mutate(
    x_bin = cut(signed_log_SR_forward, breaks = seq(-1, 1.3, length.out = bins + 1), include.lowest = TRUE, labels = FALSE),
    y_bin = cut(signed_log_SR_reverse, breaks = seq(-1, 1.3, length.out = bins + 1), include.lowest = TRUE, labels = FALSE)
  ) %>%
  count(sub_type2, x_bin, y_bin) %>%
  group_by(sub_type2) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    x_mid = -1 + (x_bin - 0.5) * (1.3 - (-1)) / bins,
    y_mid = -1 + (y_bin - 0.5) * (1.3 - (-1)) / bins
  )

binned_temp <- results_wide_temp %>%
  group_by(sub_type2) %>%
  mutate(
    x_bin = cut(signed_log_SR_forward, breaks = seq(-1, 1.3, length.out = bins + 1), include.lowest = TRUE, labels = FALSE),
    y_bin = cut(signed_log_SR_reverse, breaks = seq(-1, 1.3, length.out = bins + 1), include.lowest = TRUE, labels = FALSE)
  ) %>%
  count(sub_type2, x_bin, y_bin) %>%
  group_by(sub_type2) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    x_mid = -1 + (x_bin - 0.5) * (1.3 - (-1)) / bins,
    y_mid = -1 + (y_bin - 0.5) * (1.3 - (-1)) / bins
  )

binned_lyt_NPTI <- binned_lyt %>% filter(sub_type2 == "Non-palindromic\nType I_RM")
binned_lyt_PTI <- binned_lyt %>% filter(sub_type2 == "Palindromic\nType I RM")
binned_lyt_NPTII <- binned_lyt %>% filter(sub_type2 == "RM Type IIG & IIS")
binned_lyt_PTII <- binned_lyt %>% filter(sub_type2 == "RM Type IIP")
binned_lyt_TIII <- binned_lyt %>% filter(sub_type2 == "RM Type III")
binned_lyt_PAMS <- binned_lyt %>% filter(sub_type2 == "PAMS")

binned_temp_NPTI <- binned_temp %>% filter(sub_type2 == "Non-palindromic\nType I_RM")
binned_temp_PTI <- binned_temp %>% filter(sub_type2 == "Palindromic\nType I RM")
binned_temp_NPTII <- binned_temp %>% filter(sub_type2 == "RM Type IIG & IIS")
binned_temp_PTII <- binned_temp %>% filter(sub_type2 == "RM Type IIP")
binned_temp_TIII <- binned_temp %>% filter(sub_type2 == "RM Type III")
binned_temp_PAMS <- binned_temp %>% filter(sub_type2 == "PAMS")

# 4.0 Plotting ####

## 4.1 new style main fig ####

p1_lyt <- ggplot(binned_lyt, aes(x = x_mid, y = y_mid, fill = log10(prop + 1e-6))) +
  geom_tile(width = diff(range(-1,1.3))/30, height = diff(range(-1,1.3))/30) +
  facet_wrap(~sub_type2) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = c(- 0.553,  0.553),linewidth = 0.75) +
  geom_vline(xintercept = c(- 0.553,  0.553),linewidth = 0.75) +
  theme_bw() +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(fill = "log10(Proportion)") +
  scale_fill_viridis_c() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

#p1_lyt

p1_temp <- ggplot(binned_temp, aes(x = x_mid, y = y_mid, fill = log10(prop + 1e-6))) +
  geom_tile(width = diff(range(-1,1.3))/30, height = diff(range(-1,1.3))/30) +
  facet_wrap(~sub_type2) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = c(- 0.553,  0.553),linewidth = 0.75) +
  geom_vline(xintercept = c(- 0.553,  0.553),linewidth = 0.75) +
  theme_bw() +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(fill = "log10(Proportion)") +
  scale_fill_viridis_c() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

#p1_temp

p2_lyt <- ggplot(data=results_wide_lyt, aes(x=sub_type2, y = ,fill = result)) + 
  geom_bar(position = "fill")+
  scale_fill_manual(
    values = c(
      "Non-significant"        = "#3D3D3D",
      "Depleted and enriched"  = "#12436D",
      "Double depleted"        = "#28A197",
      "Double enriched"        = "#801650",
      "Single depleted"        = "#F46A25",
      "Single enriched"        = "#A285D1"
    ), guide = guide_legend(ncol = 2)
  ) +
  labs(fill = "Strand specific\nsignificance")+
  xlab("Defence system type") +
  ylab("Proportion") +
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "bottom"
  )


p2_lyt

p2_temp <- ggplot(data=results_wide_temp, aes(x=sub_type2, y = ,fill = result)) + 
  geom_bar(position = "fill")+
  scale_fill_manual(
    values = c(
      "Non-significant"        = "#3D3D3D",
      "Depleted and enriched"  = "#12436D",
      "Double depleted"        = "#28A197",
      "Double enriched"        = "#801650",
      "Single depleted"        = "#F46A25",
      "Single enriched"        = "#A285D1"
    ), guide = guide_legend(ncol = 2)
  ) +
  labs(fill = "Strand specific\nsignificance")+
  xlab("Defence system type") +
  ylab("Proportion") +
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "bottom"
  )

#p2_temp

# cowplot getlegend is broken so...
get_legend2 <- function(plot, legend = NULL) {
  if (is.ggplot(plot)) {
    gt <- ggplotGrob(plot)
  } else {
    if (is.grob(plot)) {
      gt <- plot
    } else {
      stop("Plot object is neither a ggplot nor a grob.")
    }
  }
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  indices <- grep(pattern, gt$layout$name)
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}

p1_lyt <- p1_lyt + theme(
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 16),
  legend.key.size = unit(1.3, "cm")
)

p2_lyt <- p2_lyt + theme(
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 16),
  legend.key.size = unit(1.3, "cm")
)


legend1 <- get_legend2(p1_lyt)
legend2 <- get_legend2(p2_lyt)

p1_lyt <- p1_lyt + theme(legend.position = "none")
p2_lyt <- p2_lyt + theme(legend.position = "none")
p1_temp <- p1_temp + theme(legend.position = "none")
p2_temp <- p2_temp + theme(legend.position = "none")



p1_lyt <- p1_lyt + plot_annotation(title = "A) Lytic phages",theme = theme(plot.title = element_text(size = 16)))
p1_temp <- p1_temp + plot_annotation(title = "B) Temperate phages",theme = theme(plot.title = element_text(size = 16)))

p2_lyt <- p2_lyt + plot_annotation(title = " ",theme = theme(plot.title = element_text(size = 16)))
p2_temp <- p2_temp + plot_annotation(title = " ",theme = theme(plot.title = element_text(size = 16)))


figure <- ggarrange(p1_lyt,p2_lyt,p1_temp,p2_temp,legend1,legend2,
                    ncol = 2, nrow = 3,
                    widths =c(2,1),
                    heights =c(3,3,1) )

jpeg("plots/combined_1.jpeg", width = 5500, height = 5700, units = "px", res = 300)

figure

dev.off()

## 4.2 SI fig split by phages ####

p3A_lyt <- ggplot(binned_lyt, aes(x = x_mid, y = y_mid, fill = log10(prop + 1e-6))) +
  geom_tile(width = diff(range(-1,1.3))/30, height = diff(range(-1,1.3))/30) +
  facet_wrap(~genome) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = c(- 0.553,  0.553),linewidth = 0.75) +
  geom_vline(xintercept = c(- 0.553,  0.553),linewidth = 0.75) +
  theme_bw() +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(fill = "log10(Proportion)") +
  scale_fill_viridis_c() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

p3A_lyt

p3A_temp <- ggplot(binned_temp, aes(x = x_mid, y = y_mid, fill = log10(prop + 1e-6))) +
  geom_tile(width = diff(range(-1,1.3))/30, height = diff(range(-1,1.3))/30) +
  facet_wrap(~genome) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = c(- 0.553,  0.553),linewidth = 0.75) +
  geom_vline(xintercept = c(- 0.553,  0.553),linewidth = 0.75) +
  theme_bw() +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(fill = "log10(Proportion)") +
  scale_fill_viridis_c() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

p3A_temp

p3B_lyt <- ggplot(results_wide_lyt, aes(x = signed_log_SR_forward, y = signed_log_SR_reverse)) +
  geom_pointdensity(adjust= 0.1) +
  facet_wrap(~genome,ncol = 8)+
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = c(- 0.553,  0.553),linewidth = 0.5) +
  geom_vline(xintercept = c(- 0.553,  0.553),linewidth = 0.5) +
  theme_bw() +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(colour = "Density")+
  scale_color_viridis(trans = "log10")+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))


p3B_lyt

p3B_temp <- ggplot(results_wide_temp, aes(x = signed_log_SR_forward, y = signed_log_SR_reverse)) +
  geom_pointdensity(adjust= 0.1) +
  facet_wrap(~genome,ncol = 8)+
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = c(- 0.553,  0.553),linewidth = 0.5) +
  geom_vline(xintercept = c(- 0.553,  0.553),linewidth = 0.5) +
  theme_bw() +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(colour = "Density")+
  scale_color_viridis(trans = "log10")+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))


p3B_temp


## 4.3 split by motifs ####

results_wide_lyt <- results_wide %>%
  filter(phage_type == "lytic")

results_wide_temp <- results_wide %>%
  filter(phage_type == "temperate")




binned_lyt <- results_wide_lyt %>%
  group_by(motif,sub_type2) %>%
  mutate(
    x_bin = cut(signed_log_SR_forward, breaks = seq(-1, 1.3, length.out = bins + 1), include.lowest = TRUE, labels = FALSE),
    y_bin = cut(signed_log_SR_reverse, breaks = seq(-1, 1.3, length.out = bins + 1), include.lowest = TRUE, labels = FALSE)
  ) %>%
  count(motif, x_bin, y_bin) %>%
  group_by(motif) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    x_mid = -1 + (x_bin - 0.5) * (1.3 - (-1)) / bins,
    y_mid = -1 + (y_bin - 0.5) * (1.3 - (-1)) / bins
  )

binned_temp <- results_wide_temp %>%
  group_by(motif,sub_type2) %>%
  mutate(
    x_bin = cut(signed_log_SR_forward, breaks = seq(-1, 1.3, length.out = bins + 1), include.lowest = TRUE, labels = FALSE),
    y_bin = cut(signed_log_SR_reverse, breaks = seq(-1, 1.3, length.out = bins + 1), include.lowest = TRUE, labels = FALSE)
  ) %>%
  count(motif, x_bin, y_bin) %>%
  group_by(motif) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    x_mid = -1 + (x_bin - 0.5) * (1.3 - (-1)) / bins,
    y_mid = -1 + (y_bin - 0.5) * (1.3 - (-1)) / bins
  )

p4l1_lyt <- ggplot(binned_lyt_NPTI, aes(x = x_mid, y = y_mid, fill = log10(prop + 1e-6))) +
  geom_tile(width = diff(range(-1,1.3))/30, height = diff(range(-1,1.3))/30) +
  facet_wrap(~motif,ncol = 8)+
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = c(- 0.553,  0.553),linewidth = 0.75) +
  geom_vline(xintercept = c(- 0.553,  0.553),linewidth = 0.75) +
  theme_bw() +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(fill = "log10(Proportion)") +
  scale_fill_viridis_c() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

jpeg("plots/p4l1_lyt.jpeg", width = 4000, height = 3000, units = "px", res = 300)

p4l1_lyt

dev.off()

results_wide_lytNPTI <- results_wide_lyt %>%
  filter(sub_type2 == "Non-palindromic\nType I_RM")

p5l1_lyt <- ggplot(results_wide_lytNPTI, aes(x = signed_log_SR_forward, y = signed_log_SR_reverse)) +
  geom_pointdensity(adjust= 0.1) +
  facet_wrap(~motif,ncol = 8)+
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = c(- 0.553,  0.553),linewidth = 0.5) +
  geom_vline(xintercept = c(- 0.553,  0.553),linewidth = 0.5) +
  theme_bw() +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(colour = "Density")+
  scale_color_viridis(trans = "log10")+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))


p5l1_lyt


jpeg("plots/p5l1_lyt.jpeg", width = 4000, height = 3000, units = "px", res = 300)

p5l1_lyt

dev.off()




# test <- results_wide %>% filter(genome == "3")
# 
# summary(test$sub_type)
# 
# test$sub_type2 <- as.factor(test$sub_type2)
# 
# test <- test %>%
#   mutate(sub_type2 = case_when(
#     type == "RM_Type_I" & Palindromic == "Non-palindromic" ~ "Non-palindromic\nType I_RM",
#     type == "RM_Type_I" & Palindromic == "Palindromic" ~ "Palindromic\nType I RM",
#     sub_type == "RM_Type_IIG" ~ "RM Type IIG & IIS",
#     sub_type == "RM_Type_IIS" ~ "RM Type IIG & IIS",
#     sub_type == "RM_Type_IIP" ~ "RM Type IIP",
#     sub_type == "RM_Type_III" ~ "RM Type III",
#     TRUE ~ sub_type)) %>%
#   filter(sub_type != "RM_Type_IIM")






## 4.3 incorporating abundances ####

results_wide_temp
results_wide_lyt

hist(results_wide_lyt$motif_abundance)

results_wide_lyt <- results_wide_lyt %>%
  mutate(
    abundance_percentile = percent_rank(motif_abundance) * 100,
    abundance_bin = round(abundance_percentile / 20) * 10
  )

results_wide_lyt$abundance_bin <- as.factor(results_wide_lyt$abundance_bin)

p3_lyt <- ggplot(data=results_wide_lyt, aes(x= abundance_bin, y = ,fill = result)) + 
  geom_bar(position = "fill")+
  facet_wrap(~sub_type2) +
  scale_fill_manual(
    values = c(
      "Non-significant"        = "#3D3D3D",
      "Depleted and enriched"  = "#12436D",
      "Double depleted"        = "#28A197",
      "Double enriched"        = "#801650",
      "Single depleted"        = "#F46A25",
      "Single enriched"        = "#A285D1"
    ), guide = guide_legend(ncol = 2)
  ) +
  labs(fill = "Strand specific\nsignificance")+
  xlab("Defence system type") +
  ylab("Proportion") +
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "bottom"
  )

p3_lyt




























## 4.3 SCRAP Older style ####

# facet by major types, keep palindromes, use geom_point 

p1B <- ggplot(results_wide,aes(x = signed_log_SR_forward, y = signed_log_SR_reverse))+
  geom_pointdensity(adjust= 0.1,size = 0.7) +
  facet_wrap(~type)+
  xlim(-1,1.3) +
  ylim(-1,1.3) +
  geom_hline(yintercept=0,linewidth = 0.5,linetype = "dashed") +
  geom_vline(xintercept=0,linewidth = 0.5,linetype = "dashed") +
  geom_hline(yintercept = c(-0.471, 0.471))+
  geom_vline(xintercept = c(-0.471, 0.471))+
  theme_bw() +
  xlab("Signed log10 standardized residual\n(observed − expected, fwd strand)") +
  ylab("Signed log10 standardized residual\n(observed − expected, rvs strand)") +
  labs(colour = "Density")+
  theme(legend.position = "bottom") +
  scale_color_viridis(trans = "log10")

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

