setwd("/Users/wendergomes/Documents/Christine Plant and Fecal samples/R analysis")

library(tidyverse)
library(mixOmics)
library(ggpubr)
library(CoDaSeq)
library(vegan)

data <- read_csv("Artemisia_project_quant_Blks_Removed_th0.3_Imputed.csv") %>% arrange(Sample)
metadata <- read_csv("Artemisia_Metadata.csv") %>% arrange(Sample)

all(data$Sample == metadata$Sample)

# Generate PCA raw
PCA_whole <- mixOmics::pca(column_to_rownames(data, var = "Sample"), ncomp = 2, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X, metadata)
PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "Attribute_type", alpha = 0.6, 
            title = paste("PCA RAW -", "Attribute_type", sep = " "),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_scores %>% group_by(Attribute_type) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = Attribute_type), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PERMANOVA
dist_metabolites <- vegdist(column_to_rownames(data, var = "Sample"), method = "euclidean")
permanova <- adonis2(dist_metabolites ~ PCA_whole_scores$Attribute_type, PCA_whole_scores, na.action = na.omit)

# Generate PCA relative abundance
data_ra <- (column_to_rownames(data, var = "Sample"))/rowSums(column_to_rownames(data, var = "Sample"))
PCA_whole_ra <- mixOmics::pca(data_ra, ncomp = 2, scale = TRUE)
PCA_whole_ra_scores <- data.frame(PCA_whole_ra$variates$X, metadata)
PCA_ra_plot <- PCA_whole_ra_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "Attribute_type", alpha = 0.6, 
            title = paste("PCA RA -", "Attribute_type", sep = " "),
            xlab = paste("PC1 (", round(PCA_whole_ra$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole_ra$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_ra_scores %>% group_by(Attribute_type) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = Attribute_type), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PERMANOVA
dist_metabolites_ra <- vegdist(data_ra, method = "euclidean")
permanova_ra <- adonis2(dist_metabolites_ra ~ PCA_whole_ra_scores$Attribute_type, PCA_whole_ra_scores, na.action = na.omit)

# Generate PCA CLR
data_clr <- codaSeq.clr(column_to_rownames(data, var = "Sample") + 1) %>% as.data.frame()
PCA_whole_clr <- mixOmics::pca(data_clr, ncomp = 2, scale = FALSE)
PCA_whole_clr_scores <- data.frame(PCA_whole_clr$variates$X, metadata)
PCA_clr_plot <- PCA_whole_clr_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "Attribute_type", alpha = 0.6, 
            title = paste("PCA CLR -", "Attribute_type", sep = " "),
            xlab = paste("PC1 (", round(PCA_whole_clr$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole_clr$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_clr_scores %>% group_by(Attribute_type) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = Attribute_type), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PERMANOVA
dist_metabolites_clr <- vegdist(data_clr, method = "euclidean")
permanova_clr <- adonis2(dist_metabolites_clr ~ PCA_whole_clr_scores$Attribute_type, PCA_whole_clr_scores, na.action = na.omit)

#save and plotting PCA
ggsave(plot = PCA_clr_plot, filename = "PCA_CLR_Artemisia.svg", device = "svg", dpi = "retina")
