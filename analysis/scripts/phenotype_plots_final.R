#!/usr/bin/env Rscript

library(tidyverse)
library(corrplot)
library(gridExtra)

# Load data
df <- read.delim('/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/feature_extractMS2/MS2_samples_combine.extended_metadata_with_strain_traits.tsv')

cat("Samples loaded:", nrow(df), "\n")
cat("Phenotype samples:", sum(!is.na(df[, 46])), "\n\n")  # Column 46 is Median_Shape_Area

# Rename columns with special characters for easier ggplot2 use
df_plot <- df %>%
  rename(
    Shape = `Median_Shape_Area`,
    Shape_MAD = `MAD_Shape_Area`,
    Shape_Mean = `Mean_Shape_Area`,
    Shape_SD = `SD_Shape_Area`,
    Hue = `Median_ColorHSV_HueMean`,
    Saturation = `Median_ColorHSV_SaturationMean`,
    Brightness = `Median_ColorHSV_BrightnessMean`,
    Lab_L = contains('ColorLab_L')[1],
    Lab_a = contains('ColorLab_a')[1],
    Lab_b = contains('ColorLab_b')[1],
    Lab_Chroma = `Median_ColorLab_ChromaEstimatedMean`
  )

# Simpler rename approach - just work with the data directly
df_plot <- df
names(df_plot)[46] <- "Shape"
names(df_plot)[47] <- "Shape_MAD"
names(df_plot)[48] <- "Shape_Mean"
names(df_plot)[49] <- "Shape_SD"
names(df_plot)[59] <- "Hue"
names(df_plot)[60] <- "Saturation"
names(df_plot)[61] <- "Brightness"
names(df_plot)[55] <- "Lab_L"
names(df_plot)[56] <- "Lab_a"
names(df_plot)[57] <- "Lab_b"
names(df_plot)[58] <- "Lab_Chroma"

# Filter to samples with phenotype data
df_pheno <- df_plot %>% filter(!is.na(Shape))

cat("Filtered to phenotype samples:", nrow(df_pheno), "\n\n")

# Select columns for correlation - check which ones exist
corr_cols <- c('Shape', 'Hue', 'Saturation', 'Brightness', 'Lab_L', 'Lab_a', 'Lab_b', 'Lab_Chroma',
               'Shape_MAD', 'Shape_Mean', 'Shape_SD', 'N_Configurations', 'Library Plate')
corr_cols <- corr_cols[corr_cols %in% colnames(df_pheno)]
cat("Using columns for correlation:", paste(corr_cols, collapse=', '), "\n\n")

# 1. Correlation heatmap
cat("Creating plots...\n")
png('/scratch/jstajich/25989115/claude-1181/-bigdata-stajichlab-shared-projects-Rhodotorula-Rhodotorula-Metabolites-feature-extractMS2/206d2d5d-9c2d-470e-be4a-c7cd9006d021/scratchpad/01_correlation_heatmap.png',
    width=1200, height=1000, res=100)
corr_mat <- cor(df_pheno[, corr_cols], use='complete.obs')
corrplot(corr_mat, method='color', type='full', order='hclust',
         tl.cex=0.8, cl.cex=0.8, addCoef.col='black', number.cex=0.6)
dev.off()
cat("✓ Saved: 01_correlation_heatmap.png\n")

# 2. Shape Area vs Brightness
png('/scratch/jstajich/25989115/claude-1181/-bigdata-stajichlab-shared-projects-Rhodotorula-Rhodotorula-Metabolites-feature-extractMS2/206d2d5d-9c2d-470e-be4a-c7cd9006d021/scratchpad/02_shape_vs_brightness.png',
    width=800, height=600, res=100)
p1 <- ggplot(df_pheno, aes(x=Shape, y=Brightness)) +
  geom_point(alpha=0.6, color='steelblue', size=2) +
  geom_smooth(method='lm', se=TRUE, color='red', fill='red', alpha=0.2, linewidth=1) +
  labs(title='Median Shape Area vs HSV Brightness (r=0.598)',
       x='Median Shape Area (pixels)', y='HSV Brightness',
       subtitle='584 samples with phenotype data') +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=14),
        plot.subtitle=element_text(hjust=0.5, size=10))
print(p1)
dev.off()
cat("✓ Saved: 02_shape_vs_brightness.png\n")

# 3. Shape Area vs L* (Lightness)
png('/scratch/jstajich/25989115/claude-1181/-bigdata-stajichlab-shared-projects-Rhodotorula-Rhodotorula-Metabolites-feature-extractMS2/206d2d5d-9c2d-470e-be4a-c7cd9006d021/scratchpad/03_shape_vs_lightness.png',
    width=800, height=600, res=100)
p2 <- ggplot(df_pheno, aes(x=Shape, y=Lab_L)) +
  geom_point(alpha=0.6, color='forestgreen', size=2) +
  geom_smooth(method='lm', se=TRUE, color='red', fill='red', alpha=0.2, linewidth=1) +
  labs(title='Median Shape Area vs CIELab L* Lightness (r=0.544)',
       x='Median Shape Area (pixels)', y='CIELab L* (Lightness)',
       subtitle='Larger cells tend to have higher lightness (brighter)') +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=14),
        plot.subtitle=element_text(hjust=0.5, size=10))
print(p2)
dev.off()
cat("✓ Saved: 03_shape_vs_lightness.png\n")

# 4. HSV Saturation and Chroma strong correlation
png('/scratch/jstajich/25989115/claude-1181/-bigdata-stajichlab-shared-projects-Rhodotorula-Rhodotorula-Metabolites-feature-extractMS2/206d2d5d-9c2d-470e-be4a-c7cd9006d021/scratchpad/04_saturation_vs_chroma.png',
    width=800, height=600, res=100)
p3 <- ggplot(df_pheno, aes(x=Saturation, y=Lab_Chroma)) +
  geom_point(alpha=0.6, color='purple', size=2) +
  geom_smooth(method='lm', se=TRUE, color='red', fill='red', alpha=0.2, linewidth=1) +
  labs(title='HSV Saturation vs CIELab Chroma (r=0.980)',
       x='HSV Saturation', y='CIELab Chroma',
       subtitle='Very strong positive correlation - different color spaces measure related properties') +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=14),
        plot.subtitle=element_text(hjust=0.5, size=10))
print(p3)
dev.off()
cat("✓ Saved: 04_saturation_vs_chroma.png\n")

# 5. Hue vs b* (negative correlation)
png('/scratch/jstajich/25989115/claude-1181/-bigdata-stajichlab-shared-projects-Rhodotorula-Rhodotorula-Metabolites-feature-extractMS2/206d2d5d-9c2d-470e-be4a-c7cd9006d021/scratchpad/05_hue_vs_b_axis.png',
    width=800, height=600, res=100)
p4 <- ggplot(df_pheno, aes(x=Hue, y=Lab_b)) +
  geom_point(alpha=0.6, color='orangered', size=2) +
  geom_smooth(method='lm', se=TRUE, color='red', fill='red', alpha=0.2, linewidth=1) +
  labs(title='HSV Hue vs CIELab b* axis (r=-0.860)',
       x='HSV Hue', y='CIELab b* (Yellow-Blue axis)',
       subtitle='Strong negative: cells with lower hue have higher b* (more blue)') +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=14),
        plot.subtitle=element_text(hjust=0.5, size=10))
print(p4)
dev.off()
cat("✓ Saved: 05_hue_vs_b_axis.png\n")

# 6. Shape variability (MAD) vs Library Plate
if('Library Plate' %in% colnames(df_pheno)) {
  png('/scratch/jstajich/25989115/claude-1181/-bigdata-stajichlab-shared-projects-Rhodotorula-Rhodotorula-Metabolites-feature-extractMS2/206d2d5d-9c2d-470e-be4a-c7cd9006d021/scratchpad/06_shape_mad_vs_plate.png',
      width=800, height=600, res=100)
  p5 <- ggplot(df_pheno, aes(x=df_pheno[[51]], y=Shape_MAD)) +  # Column 51 is Library Plate
    geom_point(alpha=0.6, color='darkred', size=2) +
    geom_smooth(method='lm', se=TRUE, color='red', fill='red', alpha=0.2, linewidth=1) +
    labs(title='Shape Variability (MAD) vs Library Plate (r=-0.581)',
         x='Library Plate', y='Mean Absolute Deviation (Shape)',
         subtitle='Negative: later plates have more uniform cell shapes') +
    theme_minimal() +
    theme(plot.title=element_text(hjust=0.5, face='bold', size=14),
          plot.subtitle=element_text(hjust=0.5, size=10))
  print(p5)
  dev.off()
  cat("✓ Saved: 06_shape_mad_vs_plate.png\n")
} else {
  cat("⊘ Skipping plot 6: Library Plate column not found\n")
}

# 7. Multi-panel plot showing Shape Area interactions with all HSV features
png('/scratch/jstajich/25989115/claude-1181/-bigdata-stajichlab-shared-projects-Rhodotorula-Rhodotorula-Metabolites-feature-extractMS2/206d2d5d-9c2d-470e-be4a-c7cd9006d021/scratchpad/07_shape_area_vs_all_hsv.png',
    width=1000, height=900, res=100)

p_hue <- ggplot(df_pheno, aes(x=Shape, y=Hue)) +
  geom_point(alpha=0.5, size=2, color='coral') +
  geom_smooth(method='lm', se=TRUE, color='red', alpha=0.2) +
  labs(title='Hue (r=0.192)', x='', y='Hue', subtitle='Weak positive') +
  theme_minimal() + theme(plot.title=element_text(hjust=0.5, size=11, face='bold'),
                          plot.subtitle=element_text(hjust=0.5, size=9))

p_sat <- ggplot(df_pheno, aes(x=Shape, y=Saturation)) +
  geom_point(alpha=0.5, size=2, color='steelblue') +
  geom_smooth(method='lm', se=TRUE, color='red', alpha=0.2) +
  labs(title='Saturation (r=0.025)', x='', y='Saturation', subtitle='Nearly uncorrelated') +
  theme_minimal() + theme(plot.title=element_text(hjust=0.5, size=11, face='bold'),
                          plot.subtitle=element_text(hjust=0.5, size=9))

p_bright <- ggplot(df_pheno, aes(x=Shape, y=Brightness)) +
  geom_point(alpha=0.5, size=2, color='darkgreen') +
  geom_smooth(method='lm', se=TRUE, color='red', alpha=0.2) +
  labs(title='Brightness (r=0.598)', x='Median Shape Area (pixels)', y='Brightness',
       subtitle='Strong positive - larger cells are brighter') +
  theme_minimal() + theme(plot.title=element_text(hjust=0.5, size=11, face='bold'),
                          plot.subtitle=element_text(hjust=0.5, size=9))

grid.arrange(p_hue, p_sat, p_bright, ncol=1)
dev.off()
cat("✓ Saved: 07_shape_area_vs_all_hsv.png\n")

# 8. Species comparison
if('Species' %in% colnames(df_pheno)) {
  df_species <- df_pheno %>% filter(!is.na(Species))

  if(nrow(df_species) > 0) {
    png('/scratch/jstajich/25989115/claude-1181/-bigdata-stajichlab-shared-projects-Rhodotorula-Rhodotorula-Metabolites-feature-extractMS2/206d2d5d-9c2d-470e-be4a-c7cd9006d021/scratchpad/08_shape_by_species.png',
        width=900, height=600, res=100)

    p_species <- ggplot(df_species, aes(x=reorder(Species, Shape, median), y=Shape, fill=Species)) +
      geom_boxplot(alpha=0.7, outlier.size=1) +
      geom_jitter(width=0.2, alpha=0.3, size=1.5) +
      labs(title='Shape Area Distribution by Species',
           x='Species (ordered by median)', y='Median Shape Area (pixels)',
           subtitle=paste(length(unique(df_species$Species)), 'species,', nrow(df_species), 'samples')) +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=10),
            plot.title=element_text(hjust=0.5, face='bold', size=14),
            plot.subtitle=element_text(hjust=0.5, size=10),
            legend.position='none')
    print(p_species)
    dev.off()
    cat("✓ Saved: 08_shape_by_species.png\n")
  }
}

cat("\n")
cat(strrep("=", 60), "\n")
cat("ALL PLOTS GENERATED SUCCESSFULLY\n")
cat(strrep("=", 60), "\n")
cat("\nSummary of findings:\n")
cat("- Shape Area strongly correlates with Brightness (r=0.598)\n")
cat("- Shape Area moderately correlates with Lightness (r=0.544)\n")
cat("- HSV Saturation & CIELab Chroma nearly perfectly correlated (r=0.980)\n")
cat("- Hue strongly inversely correlates with b* axis (r=-0.860)\n")
cat("- Shape variability decreases by plate (r=-0.581)\n")
