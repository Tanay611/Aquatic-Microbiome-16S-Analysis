# Aquatic Microbiome — 16S rRNA Amplicon Analysis

Comparative metagenomics study of water column vs. sediment microbial communities at the **Baltic Sea Science Centre (BSSC), Skansen, Stockholm**.

**Course:** BB2560 Advanced Microbiology and Metagenomics, KTH Royal Institute of Technology  
**Sequencing:** Illumina MiSeq at NGI/SciLifeLab, Solna

## Overview

- **75 samples** — water column (pelagic) and sediment (benthic) from 3 aquarium environments
- **Target region:** V3–V4 of 16S rRNA gene (primers 341F / 805R)
- **Pipeline:** DADA2 (quality filtering → denoising → chimera removal → GTDB taxonomy)
- **Analyses:** Alpha diversity (Shannon, Wilcoxon), Beta diversity (Bray-Curtis, NMDS), ANOSIM, PERMANOVA

## Key findings

- Sediment showed significantly higher alpha diversity (median Shannon 5.49 vs 3.04, p = 2.1 × 10⁻¹⁰)
- Clear community separation between habitats (ANOSIM R = 0.56, p < 0.001)
- Nutrient levels (NO₃⁻ + PO₄³⁻) explained ~20% of beta diversity variation (PERMANOVA R² = 0.20)
- Dominant phyla: Pseudomonadota, Bacteroidota (water); Actinomycetota, Chloroflexota (sediment)

## R packages

`dada2` · `vegan` · `pheatmap` · `edgeR`

## File

| File | Description |
|------|-------------|
| `metagenomics_analysis.R` | Complete analysis pipeline from raw FASTQ to statistical results |
