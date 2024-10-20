setwd("/Users/wendergomes/Documents/Christine Plant and Fecal samples/DIABLO multiomics/Artemisia data")

# (1) Load Libraries
library(mixOmics)
library(dplyr)
library(tidyr)
library(igraph)
library(BiocParallel)

set.seed(123) # For reproducibility

# (2) Load and prepare the data

# Load metabolomics data
metabolomics_data <- read.csv("Artemisia_Metabolomics_CLR.csv", row.names = 1)

# Load microbiome data
microbiome_data <- read.csv("Artemisia_Microbiome_CLR.csv", row.names = 1)

# Load sample phenotype data
sample_phenotype <- read.csv("Phenotype_samples.csv")

# Load OTU phenotype data
otu_phenotype <- read.csv("Phenotype_OTUs.csv")

# Ensure the sample IDs are consistent across the datasets
stopifnot(all(rownames(metabolomics_data) %in% sample_phenotype$SampleID))
stopifnot(all(rownames(microbiome_data) %in% sample_phenotype$SampleID))

# Set row names for sample phenotype data
row.names(sample_phenotype) <- sample_phenotype$SampleID
sample_phenotype$SampleID <- NULL

# Create the response variable (Y) from the sample phenotype
Y <- as.factor(sample_phenotype$Species)

# Convert dataframes to matrices
metabolomics_matrix <- as.matrix(metabolomics_data)
microbiome_matrix <- as.matrix(microbiome_data)

# Remove zero variance features
metabolomics_matrix <- metabolomics_matrix[, apply(metabolomics_matrix, 2, var) != 0]
microbiome_matrix <- microbiome_matrix[, apply(microbiome_matrix, 2, var) != 0]

# Create data list for integration
data <- list(Metabolomics = metabolomics_matrix, Microbiome = microbiome_matrix)

# Ensure dimensions are correct
print(lapply(data, dim))
print(summary(Y))

# (3) Initial Analysis: Pairwise PLS Comparisons
list.keepX <- c(100, 100) # select arbitrary values of features to keep

# Generate pairwise PLS models
pls1 <- spls(data[["Metabolomics"]], data[["Microbiome"]], keepX = list.keepX)

# Plot features of PLS
plotVar(pls1, cutoff = 0.5, title = "(a) Metabolomics vs Microbiome", 
        legend = c("Metabolomics", "Microbiome"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('#C7D6F0', '#EBB0A6'))

# Calculate correlation
cor(pls1$variates$X, pls1$variates$Y)

#Initial DIABLO Model: Design

# Create the design matrix "values above 0.5 will cause a reduction in predictive ability of the model"
design <- matrix(0.5, ncol = length(data), nrow = length(data), 
                 dimnames = list(names(data), names(data)))
diag(design) <- 0

# Form the basic DIABLO model
basic_diablo_model <- block.splsda(X = data, Y = Y, ncomp = 3, design = design)

# Tuning the number of components
perf_diablo <- perf(basic_diablo_model, validation = 'Mfold', folds = 5, nrepeat = 9)
plot(perf_diablo)

# Set the optimal ncomp value
ncomp <- perf_diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
print(perf_diablo$choice.ncomp$WeightedVote)

# (4) Feature Selection Tuning
# Set a grid of values for each component to test
test_keepX <- list(
  Metabolomics = c(5:9, seq(10, 18, 2), seq(20, 30, 5)), 
  Microbiome = c(5:9, seq(10, 18, 2), seq(20, 30, 5))
)

# Determine the minimum number of samples in any class
min_samples_per_class <- min(table(Y))
num_folds <- min(5, min_samples_per_class)

# Set up parallel processing

# Run the feature selection tuning without parallel processing for debugging
tune_TCGA <- tryCatch({
  tune.block.splsda(
    X = data, Y = Y, ncomp = ncomp, 
    test.keepX = test_keepX, design = design,
    validation = 'Mfold', folds = num_folds, nrepeat = 9,
    dist = "centroids.dist"
  )
}, error = function(e) {
  message("An error occurred during tuning: ", e)
  NULL
})

if (is.null(tune_TCGA)) {
  
  # Run with parallel processing but fewer workers
  BPPARAM <- MulticoreParam(workers = 4)
  tune_TCGA <- tryCatch({
    tune.block.splsda(
      X = data, Y = Y, ncomp = ncomp, 
      test.keepX = test_keepX, design = design,
      validation = 'Mfold', folds = num_folds, nrepeat = 9,
      dist = "centroids.dist", BPPARAM = BPPARAM
    )
  }, error = function(e) {
    message("An error occurred during tuning with fewer workers: ", e)
    NULL
  })
}

if (!is.null(tune_TCGA)) {
  message("Tuning completed successfully.")
} else {
  message("Tuning failed.")
}

# Set the optimal values of features to retain
list_keepX <- tune_TCGA$choice.keepX
print(list_keepX)

# (5) Final Model
# Set the optimized DIABLO model with at least 2 components
final_diablo_model <- block.splsda(X = data, Y = Y, ncomp = max(3, ncomp), 
                                   keepX = list_keepX, design = design)

# Print the design matrix for the final model
print(final_diablo_model$design)

# Plot features selected for the first component
selected_vars <- selectVar(final_diablo_model, block = 'Metabolomics', comp = 1)$Metabolomics$name
print(selected_vars)

# Plotting DIABLO
plotDiablo(final_diablo_model, ncomp = 3)

# Ensure at least 2 components for the individual plot
ncomp_indiv <- max(3, ncomp)

# Plot individual samples
plotIndiv(final_diablo_model, ncomp = ncomp_indiv, ind.names = FALSE, legend = TRUE, title = 'DIABLO Sample Plots')

# Plot arrows with axis values
plotArrow(final_diablo_model, ncomp = ncomp_indiv, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

# Define the colors for the plots
block_dims <- lapply(data, dim)
metabolomics_color <- rep('#C7D6F0', block_dims$Metabolomics[2])
microbiome_color <- rep('#EBB0A6', block_dims$Microbiome[2])
colors <- list(metabolomics_color, microbiome_color)

# Plot the correlation circle plot interactively
plotVar(final_diablo_model, var.names = FALSE, style = 'graphics', legend = TRUE,
        pch = c(16, 17), cex = c(2, 2), col = colors)

# Save the current plot to an EPS file
dev.copy2eps(file = "correlation_circle_plot.eps", width = 12, height = 10, horizontal = FALSE, onefile = FALSE, paper = "special")

# Plot the Circos plot interactively
circosPlot(final_diablo_model, cutoff = 0.7, line = TRUE,
           color.blocks = c('#C7D6F0', '#EBB0A6'),  # Adjust to two colors
           color.cor = c("#d62d20", "#000000"), size.labels = 1.5)

# Save the current plot to an EPS file
dev.copy2eps(file = "circos_plot.eps", width = 8, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")

# Set the cutoff value
cutoff <- 0.5

# Extract the scores for the first component from both metabolites and OTUs
metabolite_scores <- final_diablo_model$variates$Metabolomics[, 1]
otu_scores <- final_diablo_model$variates$Microbiome[, 1]

# Calculate the correlations between metabolites and OTUs
metabolite_correlations <- cor(metabolomics_matrix, metabolite_scores)
otu_correlations <- cor(microbiome_matrix, otu_scores)

# Initialize an empty data frame to store the correlations
correlation_df <- data.frame(Metabolite = character(), OTU = character(), Correlation = numeric(), stringsAsFactors = FALSE)

# Loop through each combination of metabolites and OTUs to calculate pairwise correlations
for (i in 1:ncol(metabolomics_matrix)) {
  for (j in 1:ncol(microbiome_matrix)) {
    corr_value <- metabolite_correlations[i, 1] * otu_correlations[j, 1]
    if (!is.na(corr_value) && abs(corr_value) >= cutoff) {
      correlation_df <- rbind(correlation_df, data.frame(
        Metabolite = colnames(metabolomics_matrix)[i], 
        OTU = colnames(microbiome_matrix)[j], 
        Correlation = corr_value))
    }
  }
}

# Save the connections to a CSV file
write.csv(correlation_df, file = "Metabolite_OTU_Correlations.csv", row.names = FALSE)

print("Correlation data has been exported to Metabolite_OTU_Correlations.csv")

# Reset the graphics device
#graphics.off()

# Save the network plot to a high-resolution PNG file with appropriate dimensions
#png(file = "network_plot.png", width = 4000, height = 3200, res = 300) # Set high resolution (300 DPI)
#par(mar = c(4, 4, 4, 4)) # Adjust margins to provide more space and avoid the "figure margins too large" error
#network(final_diablo_model, blocks = c(1, 2), 
#        color.node = c('#C7D6F0', '#EBB0A6'), cutoff = 0.4)
#dev.off() # Close the PNG device

# Save the network in .gml format
# Generate the network object
#png(tempfile()) # Open a temporary null plotting device
#network_obj <- network(final_diablo_model, blocks = c(1, 2), 
#                       color.node = c('#C7D6F0', '#EBB0A6'), cutoff = 0.4)
#dev.off() # Close the null device

# Extract the igraph object
#g <- network_obj$gR

# Assign colors to nodes based on the data blocks
#V(g)$color <- ifelse(V(g)$name %in% colnames(metabolomics_matrix), '#C7D6F0', '#EBB0A6')

# Save the network in .gml format
#write_graph(g, file = "myNetwork.gml", format = "gml")

# Plot loadings
plotLoadings(final_diablo_model, comp = 1, contrib = 'max', method = 'median')

# Save loading plots as EPS
postscript("loadings_plot.eps", width = 12, height = 10, horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar = c(10, 10, 10, 10), oma = c(5, 5, 5, 5)) # Adjust margins to provide more space
plotLoadings(final_diablo_model, comp = 2, contrib = 'max', method = 'median')
dev.off()

# Save CIM plot as EPS
postscript("cimDiablo_plot.eps", width = 12, height = 10, horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar = c(10, 10, 10, 10), oma = c(5, 5, 5, 5)) # Adjust margins to provide more space
cimDiablo(final_diablo_model)
dev.off()

# Save CIM plot with custom trimming and colors as EPS
postscript("cimDiablo_plot_trimmed.eps", width = 12, height = 10, horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar = c(10, 10, 10, 10), oma = c(5, 5, 5, 5)) # Adjust margins to provide more space
cimDiablo(final_diablo_model, trim = c(-3, 3), color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

# run repeated CV performance evaluation
perf.diablo = perf(final_diablo_model, validation = 'Mfold', 
                   M = 10, nrepeat = 20, 
                   dist = 'centroids.dist') 

perf.diablo$MajorityVote.error.rate
perf.diablo$WeightedVote.error.rate

# Save performance evaluation plot as EPS
postscript("perf_diablo_plot.eps", width = 12, height = 10, horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar = c(10, 10, 10, 10), oma = c(5, 5, 5, 5)) # Adjust margins to provide more space
plot(perf_diablo)
dev.off()

# Calculate the AUROC for the specified DIABLO model
auc_splsda <- auroc(final_diablo_model, roc.block = "Microbiome", roc.comp = 3, print = FALSE)

# Print the AUROC results
print(auc_splsda)

#END CODE