# Supplementary figures for the R1 revision
#  (A) Sequencing depth per sample + ASV rarefaction curves (R1.5, R1.19)
#  (B) Cumulative rarefied reads per genus-level taxon, log scale (R1.5)
# Run from the project root.

library(tidyverse)
library(vegan)
library(ape)
library(ggtree)
library(pals)
library(patchwork)

dir.create("figures/supplementary", showWarnings = FALSE, recursive = TRUE)

# ---- Load + decontaminate (same steps as 3_rarefaction_AMF.Rmd) -----------
data_18S <- "data/18S_AMf_Eukaryome/Robust_MiSeq_18S_AMf_20241120_125_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table2"
lev <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

amf_OTUs <- data_18S %>%
  read_tsv(show_col_types = FALSE) %>%
  rename_with(~sub("\\_.*", "", .)) %>%
  rename_with(~gsub("-", "_", .)) %>%
  rename_with(~sub("_AMF_", "", .)) %>%
  rename_with(~sub("\\_AMF.*", "", .)) %>%
  rename_with(~sub("ROB_T_", "CONT_ROB_", .)) %>%
  rename_with(~sub("ROB_Tex", "CONT_ROB_T_ex", .)) %>%
  mutate(taxonomy = str_replace_all(taxonomy, ":", "_"),
         taxonomy = str_replace_all(taxonomy, "[|]", "; ")) %>%
  dplyr::select(OTU, amplicon, taxonomy, starts_with("CONT", ignore.case = FALSE),
                starts_with("ROB", ignore.case = FALSE)) %>%
  filter(str_detect(taxonomy, "p__Glomeromycota"))

d <- amf_OTUs %>%
  replace(amf_OTUs == 0, NA) %>%
  dplyr::select(amplicon, starts_with("CONT")) %>%
  pivot_longer(-amplicon, names_to = "samples", values_to = "reads") %>%
  filter(!is.na(reads))

amf_OTUs_decont <- amf_OTUs %>%
  replace(amf_OTUs == 0, NA) %>%
  pivot_longer(starts_with("ROB"), names_to = "samples", values_to = "reads") %>%
  filter(!is.na(reads)) %>%
  left_join(count(d, amplicon, wt = reads), by = "amplicon") %>%
  mutate(reads = case_when(is.na(n) ~ reads, n > reads ~ 0, TRUE ~ reads - n)) %>%
  dplyr::select(-n, -starts_with("CONT")) %>%
  pivot_wider(names_from = samples, values_from = reads, values_fill = 0) %>%
  as_tibble() %>%
  column_to_rownames(var = "amplicon")

mat <- amf_OTUs_decont %>% dplyr::select(starts_with("ROB"))
depth <- colSums(mat)

cut_off  <- 2101                       # minimum reads per sample (as in script 3)
kept     <- depth >= cut_off
rar_depth <- min(depth[kept])            # depth passed to rrarefy in script 3
message("Samples: ", length(depth), " | kept: ", sum(kept),
        " | removed: ", sum(!kept))
message("Exclusion threshold = ", round(cut_off),
        " | rarefaction depth (smallest retained) = ", rar_depth)
message("Samples < 2101 reads: ", sum(depth < 2101))

# ---- (A) sequencing depth + rarefaction curves ----------------------------
depth_df <- tibble(sample = names(depth), reads = depth, kept = kept) %>%
  arrange(reads) %>%
  mutate(rank = row_number(),
         status = ifelse(kept, "Retained", "Removed"))

pA <- ggplot(depth_df, aes(rank, reads, colour = status)) +
  geom_point(size = 1.6) +
  geom_hline(yintercept = rar_depth, linetype = "dashed") +
  annotate("text", x = 118, y = rar_depth, vjust = 1.6, hjust = 1, size = 3,
           label = paste0("Rarefaction depth = ", format(rar_depth, big.mark = ","))) +
  scale_y_log10(labels = scales::label_comma(), breaks = 10^(0:5)) +
  scale_colour_manual(values = c(Retained = "black", Removed = "red"), name = NULL) +
  labs(x = "Samples (ranked by depth)", y = "AM fungal reads per sample (log10)") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "top")

rc <- vegan::rarecurve(t(mat[, kept]), step = 50, tidy = TRUE)
pB <- ggplot(rc, aes(Sample, Species, group = Site)) +
  geom_line(alpha = 0.5, linewidth = 0.3) +
  geom_vline(xintercept = rar_depth, linetype = "dashed", colour = "red") +
  coord_cartesian(xlim = c(0, 30000)) +
  labs(x = "Number of reads sampled", y = "ASV richness") +
  theme_bw() + theme(panel.grid = element_blank())

figA <- pA + pB + plot_annotation(tag_levels = "a")
ggsave("figures/supplementary/rarefaction_depth_curves.pdf", figA,
       width = 22, height = 10, units = "cm")

# ---- (B) cumulative reads per genus, log scale ----------------------------
amf_genus_rare <- read.table("outputs/abundance_tables/Robust_MiSeq_18S_AMf_20241120_125_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.afm_rarefy_allgenus.txt",
                             header = TRUE, stringsAsFactors = FALSE)
genus_tree <- read.tree("outputs/phylo_tree/Robust_MiSeq_18S_AMf_20241120_125_samples_Eukaryome_Glomeromycota_genus.tree")
genus_tree <- ape::root(genus_tree, "g__Paraglomus_7")

clean_name <- function(x) x %>%
  str_remove("g__") %>%
  str_replace_all("[a-z]__", "Unk. ") %>%
  str_replace("_(\\d+)$", " (\\1 ASV)") %>%
  str_replace("_", " ")
genus_tree$tip.label <- clean_name(genus_tree$tip.label)

amf_genus_rare <- amf_genus_rare %>%
  mutate(genus = clean_name(genus)) %>%
  remove_rownames() %>% column_to_rownames("genus")

Y_genus <- amf_genus_rare %>%
  filter(rowSums(across(starts_with("ROB")) != 0) > 3) %>% t()
genus_tree <- ape::keep.tip(genus_tree, colnames(Y_genus))
message("Rarefied samples: ", nrow(Y_genus), " | genus-level taxa: ", ncol(Y_genus))

# family of each genus (same logic as script 9)
tax <- read_tsv(data_18S, show_col_types = FALSE) %>%
  mutate(taxonomy = str_replace_all(taxonomy, ":", "_"),
         taxonomy = str_replace_all(taxonomy, "[|]", "; ")) %>%
  dplyr::select(taxonomy) %>%
  filter(str_detect(taxonomy, "p__Glomeromycota")) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "(p__|c__|o__|f__|g__|s__)unclassified", "*")) %>%
  separate_wider_delim(taxonomy, names = lev, delim = "; ", cols_remove = FALSE) %>%
  unite(taxonomy, kingdom, phylum, class, order, family, genus, sep = "; ", remove = FALSE) %>%
  group_by(taxonomy) %>%
  mutate(genus = paste(str_extract(taxonomy,
           "(g_[^;*]+$)|([a-z]_[^;*]+)(?=(;\\s*\\*|\\s*\\*|$))"), n(), sep = "_")) %>%
  group_by(genus) %>% summarise(family = first(family)) %>%
  mutate(genus = clean_name(genus),
         family = str_remove(family, "f__"),
         family = str_replace(family, "[*]", "Unk."))

sum_df <- tibble(name = colnames(Y_genus),
                 reads = colSums(Y_genus),
                 n_samples = colSums(Y_genus > 0)) %>%
  left_join(tax, by = c("name" = "genus"))

fam_cols <- setNames(watlington()[c(3:7, 9:13, 8)][seq_along(unique(sum_df$family))],
                     sort(unique(sum_df$family)))

p_tree <- ggtree(genus_tree, layout = "roundrect") +
  geom_tiplab(size = 3, align = TRUE, linesize = .5) +
  theme_tree() + xlim_tree(2.2)

PA <- "Abundance
(log10 cumulative reads)"
PO <- "Occurrence
(no. of samples, of 88)"
figB <- p_tree +
  geom_facet(panel = PA, data = sum_df, geom = geom_col,
             mapping = aes(x = log10(reads), fill = family),
             orientation = "y", width = .6) +
  geom_facet(panel = PO, data = sum_df, geom = geom_segment,
             mapping = aes(x = 0, xend = n_samples, y = y, yend = y),
             colour = "grey60", linewidth = .4) +
  geom_facet(panel = PO, data = sum_df, geom = geom_point,
             mapping = aes(x = n_samples), shape = 21, fill = "white",
             colour = "black", size = 1.8) +
  scale_fill_manual(values = fam_cols, name = "Family") +
  ggtree::xlim_expand(c(0, 90), panel = PO) +
  ggtree::xlim_expand(c(0, 5), panel = PA) +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"))
figB <- ggtree::facet_widths(figB, widths = c(1.4, 1, 0.6))

ggsave("figures/supplementary/cumulative_reads_per_genus_log.pdf", figB,
       width = 22, height = 20, units = "cm")
write_csv(sum_df, "outputs/abundance_tables/cumulative_reads_per_genus.csv")
