# =============================================================================
# 16S rRNA Amplicon Sequencing — Aquatic Microbiome Analysis
# Baltic Sea Science Centre (BSSC), Skansen, Stockholm
# Course: BB2560 Advanced Microbiology and Metagenomics, KTH
# Author: Tanay Bavchikar
# =============================================================================
# Samples: 75 water + sediment samples from 3 aquarium environments
# Target:  V3-V4 region of 16S rRNA gene (bacteria & archaea)
# Primers: 341F / 805R
# Sequencing: Illumina MiSeq at NGI/SciLifeLab
# Pipeline: DADA2 -> taxonomy (GTDB) -> diversity analyses (vegan, pheatmap)
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────────

library(dada2)
library(edgeR)
library(pheatmap)
library(vegan)

# ── 1. Load sample metadata ───────────────────────────────────────────────────

sample_info        <- read.delim("sample_info.txt")
ngi_sample_id      <- sample_info[, 1]
course_id          <- sample_info[, 2]
sample_label       <- sample_info[, 3]
sample_temperature <- sample_info[, 4]
sample_salinity    <- sample_info[, 5]
sample_pH          <- sample_info[, 6]
sample_PO4_levels  <- sample_info[, 7]
sample_NO3_levels  <- sample_info[, 8]

# ── 2. Locate FASTQ files ─────────────────────────────────────────────────────

r1_files <- list.files(recursive = TRUE, pattern = "_R1.fastq", full.names = TRUE)
r2_files <- list.files(recursive = TRUE, pattern = "_R2.fastq", full.names = TRUE)

ix_fwd <- sapply(ngi_sample_id, function(id) grep(id, r1_files))
ix_rev <- sapply(ngi_sample_id, function(id) grep(id, r2_files))

fnFs <- r1_files[ix_fwd]
fnRs <- r2_files[ix_rev]

# Sanity check — forward and reverse must be in the same order
stopifnot(length(fnFs) == length(fnRs))

# ── 3. Quality inspection ─────────────────────────────────────────────────────

plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

# ── 4. Filter and trim ────────────────────────────────────────────────────────
# truncLen chosen after inspecting quality profiles above
# trimLeft removed (primers already stripped by NGI)

filtFs <- file.path("filtered", paste0(ngi_sample_id, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(ngi_sample_id, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen  = c(220, 220),
  maxN      = 0,
  maxEE     = c(2, 2),
  truncQ    = 2,
  rm.phix   = TRUE,
  compress  = TRUE,
  multithread = TRUE
)

plotQualityProfile(filtFs[1:4])
plotQualityProfile(filtRs[1:4])

# ── 5. Learn error rates ──────────────────────────────────────────────────────

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# ── 6. Dereplicate ────────────────────────────────────────────────────────────

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- course_id
names(derepRs) <- course_id

# ── 7. DADA2 denoising ───────────────────────────────────────────────────────

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# ── 8. Merge paired reads ────────────────────────────────────────────────────

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# ── 9. Sequence table + chimera removal ──────────────────────────────────────

seqtab <- makeSequenceTable(mergers)
seqtab <- removeBimeraDenovo(seqtab, method = "consensus",
                              multithread = TRUE, verbose = TRUE)
seqtab <- t(seqtab)   # ASVs as rows, samples as columns

# Rename rows from sequences to ASV IDs
asv            <- rownames(seqtab)
rownames(seqtab) <- paste0("ASV", seq_along(asv))

# Non-chimeric ASV count per sample
nonchim_df <- data.frame(
  Sample          = colnames(seqtab),
  NonChimeric_ASVs = colSums(seqtab > 0)
)
print(nonchim_df)

# ── 10. Taxonomic annotation (GTDB) ──────────────────────────────────────────

taxa_genus <- assignTaxonomy(
  asv,
  "sbdi-gtdb.assignTaxonomy_no_species_or_kingdom.fna.gz",
  multithread = TRUE,
  taxLevels   = c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
)

taxa <- addSpecies(
  taxa_genus,
  "sbdi-gtdb.20genomes.addSpecies.fna.gz",
  allowMultiple = FALSE,
  tryRC         = FALSE,
  n             = 2000,
  verbose       = FALSE
)
rownames(taxa) <- rownames(seqtab)

# Paste genus + species where species is assigned
ix_sp              <- which(!is.na(taxa[, "Species"]))
taxa[ix_sp, "Species"] <- paste(taxa[ix_sp, "Genus"], taxa[ix_sp, "Species"])

colSums(!is.na(taxa))   # how many ASVs annotated at each rank

# ── 11. Normalise (relative abundance) ───────────────────────────────────────

norm_seqtab <- apply(seqtab, 2, function(x) x / sum(x))

# ── 12. Aggregate counts by taxonomic level ───────────────────────────────────

taxonomic_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

make_clade_abundance <- function(seqtab_asv, taxa, levels = taxonomic_levels) {
  out <- list()
  for (lvl in levels) {
    grp        <- taxa[, lvl]
    grp[is.na(grp) | grp == ""] <- "Unassigned"
    out[[lvl]] <- rowsum(seqtab_asv, group = grp)   # taxon x samples
  }
  out
}

taxa2           <- taxa[rownames(seqtab), , drop = FALSE]
clade_abundance <- make_clade_abundance(seqtab, taxa2)

# Normalised version
norm_clade_abundance <- make_clade_abundance(norm_seqtab, taxa2)

# ── 13. Habitat indices ───────────────────────────────────────────────────────

ord        <- match(colnames(seqtab), course_id)
label_ab   <- as.character(sample_label[ord])
NO3_ab     <- as.numeric(sample_NO3_levels[ord])
PO4_ab     <- as.numeric(sample_PO4_levels[ord])

is_water   <- grepl("\\bwater\\b",   label_ab, ignore.case = TRUE)
is_sediment <- grepl("\\bsediment\\b", label_ab, ignore.case = TRUE)
habitat    <- ifelse(is_water, "Water", ifelse(is_sediment, "Sediment", "Other"))

sorting_wat <- which(is_water)
sorting_sed <- which(is_sediment)

# ── 14. Taxonomic bar graphs ──────────────────────────────────────────────────

make_taxonomic_bargraph <- function(clade_abund, sample_names, tax_level,
                                    cex_names = 0.7, top_x = 10,
                                    sorting = NULL, main = "") {
  if (is.null(sorting)) sorting <- seq_len(ncol(clade_abund[[tax_level]]))
  ok    <- sort(rowMeans(clade_abund[[tax_level]]), index.return = TRUE,
                decreasing = TRUE)$ix[seq_len(top_x)]
  mycols <- colorRampPalette(c(
    "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99",
    "#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A",
    "#FFFF99","#B15928"
  ))
  old_par <- par(mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 10))
  barplot(clade_abund[[tax_level]][ok, sorting],
          col = mycols(length(ok)), las = 2,
          cex.names = cex_names,
          names.arg = sample_names[sorting],
          main = main)
  cex_taxa <- min(1.2, 10 / length(ok))
  par(xpd = NA)
  legend(x = par("usr")[2] + 0.5, y = par("usr")[4],
         legend = rownames(clade_abund[[tax_level]])[ok],
         col = mycols(length(ok)), pch = 19, bty = "n",
         cex = cex_taxa, ncol = 1)
  par(old_par)
}

ab_cols <- colnames(clade_abundance[["Phylum"]])

# Phylum level
make_taxonomic_bargraph(clade_abundance, ab_cols, "Phylum",
                        sorting = sorting_sed, top_x = 10,
                        main = "Sediment samples — Phylum")
make_taxonomic_bargraph(clade_abundance, ab_cols, "Phylum",
                        sorting = sorting_wat, top_x = 10,
                        main = "Water samples — Phylum")

# Family level
make_taxonomic_bargraph(clade_abundance, ab_cols, "Family",
                        sorting = sorting_sed, top_x = 10,
                        main = "Sediment samples — Family")
make_taxonomic_bargraph(clade_abundance, ab_cols, "Family",
                        sorting = sorting_wat, top_x = 10,
                        main = "Water samples — Family")

# ── 15. Alpha diversity — Shannon ─────────────────────────────────────────────

shannon <- diversity(seqtab, MARGIN = 2)

barplot(shannon, las = 2, names.arg = label_ab,
        cex.names = 0.6, main = "Shannon Diversity Index per Sample")

boxplot(shannon[is_water], shannon[is_sediment],
        names = c("Water", "Sediment"),
        col   = c("skyblue", "tan"),
        ylab  = "Shannon Diversity Index",
        ylim  = c(0, 6), yaxt = "n",
        main  = "Alpha Diversity: Water vs Sediment")
axis(2, at = seq(0, 6, 0.2), cex.axis = 1.2)

# Statistical test
wilcox.test(shannon[is_water], shannon[is_sediment])
cat("Median Shannon — Sediment:", median(shannon[is_sediment]), "\n")
cat("Median Shannon — Water:   ", median(shannon[is_water]),    "\n")

# ── 16. Beta diversity — Bray-Curtis ─────────────────────────────────────────

bray_dist <- as.matrix(vegdist(t(norm_seqtab), method = "bray"))

# Heatmap: all samples
pheatmap(bray_dist,
         clustering_distance_rows = as.dist(bray_dist),
         clustering_distance_cols = as.dist(bray_dist),
         labels_row = label_ab, labels_col = label_ab,
         main = "Bray-Curtis: All Samples")

# Heatmap: water only (paired sites)
ix_wat <- which(habitat == "Water")
ix_sed <- which(habitat == "Sediment")

bray_wat <- as.matrix(vegdist(t(norm_seqtab)[ix_wat, , drop = FALSE], method = "bray"))
bray_sed <- as.matrix(vegdist(t(norm_seqtab)[ix_sed, , drop = FALSE], method = "bray"))

pheatmap(bray_wat,
         clustering_distance_rows = as.dist(bray_wat),
         clustering_distance_cols = as.dist(bray_wat),
         labels_row = label_ab[ix_wat], labels_col = label_ab[ix_wat],
         main = "Bray-Curtis: Water Only")

pheatmap(bray_sed,
         clustering_distance_rows = as.dist(bray_sed),
         clustering_distance_cols = as.dist(bray_sed),
         labels_row = label_ab[ix_sed], labels_col = label_ab[ix_sed],
         main = "Bray-Curtis: Sediment Only")

# ── 17. NMDS ordination ───────────────────────────────────────────────────────

set.seed(100)
mds <- metaMDS(as.dist(bray_dist), k = 2, trymax = 100)

col_habitat <- ifelse(habitat == "Water", "steelblue", "sienna")

par(mar = c(5, 4, 4, 8), xpd = TRUE)
plot(mds$points[, 1], mds$points[, 2],
     pch = 21, cex = 2, col = "black", bg = col_habitat,
     xlab = "NMDS1", ylab = "NMDS2",
     main = "NMDS (Bray-Curtis): Water vs Sediment")
legend("topright", legend = c("Water", "Sediment"),
       pt.bg = c("steelblue", "sienna"), pch = 21, bty = "n",
       inset = c(-0.25, 0))

# ── 18. ANOSIM ───────────────────────────────────────────────────────────────
# Tests whether between-group distances > within-group distances

these       <- which(habitat %in% c("Water", "Sediment"))
anosim_result <- anosim(as.dist(bray_dist[these, these]),
                         grouping    = habitat[these],
                         permutations = 9999)
print(anosim_result)
# Result: R = 0.5553, p = 0.0001

# ── 19. PERMANOVA (adonis2) ──────────────────────────────────────────────────
# Tests how much variance in community composition is explained by nutrients

eutro <- NO3_ab + PO4_ab

# All samples: nutrient effect + habitat
adonis2(as.dist(bray_dist) ~ eutro + habitat, permutations = 9999)
# R² = 0.200, F = 9.00, p = 0.0001

# Water only: nutrient effect
bray_wat_dist <- vegdist(t(norm_seqtab)[ix_wat, , drop = FALSE], method = "bray")
adonis2(bray_wat_dist ~ eutro[ix_wat], permutations = 9999)
# R² = 0.158, F = 10.14, p = 0.0001

# Sediment only: nutrient effect
bray_sed_dist <- vegdist(t(norm_seqtab)[ix_sed, , drop = FALSE], method = "bray")
adonis2(bray_sed_dist ~ eutro[ix_sed], permutations = 9999)
# R² = 0.226, F = 4.95, p = 0.0001
