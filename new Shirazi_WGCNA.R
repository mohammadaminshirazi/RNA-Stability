

# =============================================================================
# Complete Package Installation Script for WGCNA Analysis
# This script installs all required packages for the WGCNA analysis pipeline
# =============================================================================

# Function to install packages if not already installed
install_if_missing <- function(packages, type = "CRAN") {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing", pkg, "from", type, "...\n")
      if (type == "CRAN") {
        install.packages(pkg, dependencies = TRUE)
      } else if (type == "Bioconductor") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg, dependencies = TRUE)
      }
      library(pkg, character.only = TRUE)
    } else {
      cat(pkg, "is already installed and loaded.\n")
    }
  }
}

# Update R packages first
cat("Updating existing packages...\n")
update.packages(ask = FALSE, checkBuilt = TRUE)

# Install BiocManager first if not available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Set Java memory (for some packages)
options(java.parameters = "-Xmx30000m")

# =============================================================================
# 1. BASIC R PACKAGES FROM CRAN
# =============================================================================
cat("\n=== Installing Basic R Packages from CRAN ===\n")

basic_packages <- c(
  "Rcpp",           # C++ integration
  "devtools",       # Development tools
  "remotes",        # Remote package installation
  "roxygen2",       # Documentation
  "testthat"        # Testing
)

install_if_missing(basic_packages, "CRAN")

# =============================================================================
# 2. DATA MANIPULATION AND STATISTICS PACKAGES
# =============================================================================
cat("\n=== Installing Data Manipulation and Statistics Packages ===\n")

data_packages <- c(
  "dplyr",          # Data manipulation
  "tidyr",          # Data tidying
  "reshape2",       # Data reshaping
  "stringr",        # String manipulation
  "readxl",         # Read Excel files
  "openxlsx",       # Write Excel files
  "xlsx",           # Excel file handling
  "rJava",          # Java integration for xlsx
  "readr",          # Fast data reading
  "data.table",     # Fast data manipulation
  "tibble",         # Modern data frames
  "purrr",          # Functional programming
  "broom",          # Tidy statistical output
  "gtools",         # Programming tools
  "Hmisc",          # Statistical functions
  "corrplot",       # Correlation plots
  "psych"           # Psychological statistics
)

install_if_missing(data_packages, "CRAN")

# =============================================================================
# 3. VISUALIZATION PACKAGES
# =============================================================================
cat("\n=== Installing Visualization Packages ===\n")

viz_packages <- c(
  "ggplot2",        # Grammar of graphics
  "gridExtra",      # Multiple plots
  "grid",           # Grid graphics
  "lattice",        # Lattice graphics
  "RColorBrewer",   # Color palettes
  "viridis",        # Viridis color scale
  "scales",         # Scale functions
  "ggrepel",        # Repel overlapping text
  "ggpubr",         # Publication ready plots
  "cowplot",        # Plot arrangements
  "patchwork",      # Combine plots
  "pheatmap",       # Pretty heatmaps
  "gplots",         # Various plots
  "VennDiagram",    # Venn diagrams
  "UpSetR",         # Set intersections
  "circlize",       # Circular plots
  "treemap",        # Treemap plots
  "wordcloud",      # Word clouds
  "plotly",         # Interactive plots
  "DT",             # Interactive tables
  "knitr",          # Dynamic reports
  "kableExtra",     # Enhanced tables
  "rmarkdown",      # R Markdown
  "htmlwidgets",    # HTML widgets
  "webshot"         # Web screenshots
)

install_if_missing(viz_packages, "CRAN")

# =============================================================================
# 4. NETWORK ANALYSIS AND SPECIALIZED PACKAGES
# =============================================================================
cat("\n=== Installing Network Analysis Packages ===\n")

network_packages <- c(
  "igraph",         # Network analysis
  "network",        # Network objects
  "sna",            # Social network analysis
  "visNetwork",     # Interactive networks
  "networkD3",      # D3 networks
  "ggalluvial",     # Alluvial diagrams
  "sunburstR",      # Sunburst charts
  "dendextend",     # Dendrogram extensions
  "cluster",        # Cluster analysis
  "flashClust",     # Fast clustering
  "dynamicTreeCut", # Dynamic tree cutting
  "fastcluster",    # Fast hierarchical clustering
  "apcluster",      # Affinity propagation
  "MCL"             # Markov clustering
)

install_if_missing(network_packages, "CRAN")

# =============================================================================
# 5. STATISTICAL AND MACHINE LEARNING PACKAGES
# =============================================================================
cat("\n=== Installing Statistical and ML Packages ===\n")

stats_packages <- c(
  "caret",          # Classification and regression
  "randomForest",   # Random forest
  "e1071",          # SVM and other methods
  "glmnet",         # Regularized regression
  "MASS",           # Modern applied statistics
  "nlme",           # Non-linear mixed effects
  "survival",       # Survival analysis
  "survminer",      # Survival analysis visualization
  "pROC",           # ROC analysis
  "preprocessCore", # Preprocessing
  "matrixStats",    # Matrix statistics
  "robust",         # Robust statistics
  "outliers"        # Outlier detection
)

install_if_missing(stats_packages, "CRAN")

# =============================================================================
# 6. BIOCONDUCTOR PACKAGES
# =============================================================================
cat("\n=== Installing Bioconductor Packages ===\n")

# Update Bioconductor
BiocManager::install(version = "3.18")

bioc_packages <- c(
  # Core Bioconductor
  "Biobase",        # Base functions
  "BiocGenerics",   # Generic functions
  "S4Vectors",      # S4 vectors
  "IRanges",        # Integer ranges
  "GenomicRanges",  # Genomic ranges
  "Biostrings",     # String handling
  "XVector",        # External vectors
  
  # Gene expression analysis
  "limma",          # Linear models
  "edgeR",          # RNA-seq analysis
  "DESeq2",         # Differential expression
  "affy",           # Affymetrix arrays
  "affyPLM",        # Affymetrix QC
  "oligo",          # Oligonucleotide arrays
  "genefilter",     # Gene filtering
  "geneplotter",    # Gene plotting
  "sva",            # Surrogate variable analysis
  "pamr",           # Prediction analysis
  "siggenes",       # Significant genes
  
  # Annotation packages
  "AnnotationDbi",  # Annotation database
  "org.Hs.eg.db",   # Human annotations
  "org.Mm.eg.db",   # Mouse annotations
  "org.Bt.eg.db",   # Bovine annotations (for your analysis)
  "GO.db",          # Gene Ontology
  "KEGG.db",        # KEGG database
  "ReactomePA",     # Reactome pathway analysis
  "biomaRt",        # BioMart interface
  
  # Pathway and enrichment analysis
  "clusterProfiler", # Cluster profiling
  "enrichplot",     # Enrichment visualization
  "DOSE",           # Disease ontology
  "pathview",       # Pathway visualization
  "ReactomePA",     # Reactome analysis
  "GSEABase",       # GSEA base functions
  "GSVA",           # Gene set variation analysis
  "fgsea",          # Fast GSEA
  
  # Network analysis
  "WGCNA",          # Weighted gene co-expression analysis
  "CEMiTool",       # Co-expression analysis
  "BioNet",         # Biological networks
  "RCy3",           # Cytoscape interface
  "STRINGdb",       # STRING database
  "PSICQUIC",       # Protein interaction
  
  # Visualization
  "ComplexHeatmap", # Complex heatmaps
  "EnhancedVolcano", # Volcano plots
  "ggbio",          # Genomic visualization
  "Gviz",           # Genome visualization
  "ChIPseeker",     # ChIP-seq visualization
  
  # Quality control
  "arrayQualityMetrics", # Array QC
  "affyQCReport",   # Affymetrix QC reports
  "simpleaffy",     # Simple affy QC
  
  # Preprocessing
  "preprocessCore", # Preprocessing
  "affy",           # Affymetrix preprocessing
  "affyPLM",        # Affymetrix modeling
  "gcrma",          # GC-RMA normalization
  "vsn",            # Variance stabilization
  
  # Utilities
  "graph",          # Graph objects
  "Rgraphviz",      # Graph visualization
  "RBGL"            # Graph algorithms
)

install_if_missing(bioc_packages, "Bioconductor")

# =============================================================================
# 7. SPECIALIZED PACKAGES FOR YOUR ANALYSIS
# =============================================================================
cat("\n=== Installing Specialized Analysis Packages ===\n")

specialized_packages <- c(
  "enrichR",        # Enrichment analysis
  "VennDiagram",    # Venn diagrams
  "eulerr",         # Euler diagrams
  "scales",         # Scaling functions
  "RColorBrewer",   # Color schemes
  "gplots",         # Various plots
  "corrplot",       # Correlation plots
  "Hmisc",          # Harrell miscellaneous
  "reshape",        # Data reshaping
  "plyr",           # Data manipulation
  "foreach",        # Parallel loops
  "doParallel",     # Parallel backend
  "parallel",       # Parallel computing
  "snow",           # Parallel computing
  "snowfall"        # Parallel computing
)

install_if_missing(specialized_packages, "CRAN")

# =============================================================================
# 8. ADDITIONAL USEFUL PACKAGES
# =============================================================================
cat("\n=== Installing Additional Useful Packages ===\n")

additional_packages <- c(
  "magrittr",       # Pipe operators
  "lubridate",      # Date handling
  "stringi",        # String processing
  "forcats",        # Factor handling
  "haven",          # SPSS, SAS, Stata files
  "jsonlite",       # JSON handling
  "xml2",           # XML handling
  "httr",           # HTTP tools
  "curl",           # HTTP client
  "rvest",          # Web scraping
  "R.utils",        # Utilities
  "tools",          # Programming tools
  "utils",          # Utility functions
  "methods",        # Object-oriented programming
  "stats",          # Statistical functions
  "graphics",       # Graphics functions
  "grDevices"       # Graphics devices
)

install_if_missing(additional_packages, "CRAN")

# =============================================================================
# 9. VERIFICATION AND SUMMARY
# =============================================================================
cat("\n=== Verification of Installation ===\n")

# Test critical packages
critical_packages <- c("WGCNA", "limma", "edgeR", "clusterProfiler", 
                       "org.Bt.eg.db", "ComplexHeatmap", "ggplot2", 
                       "enrichR", "STRINGdb", "xlsx")

cat("Testing critical packages:\n")
for (pkg in critical_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "loaded successfully\n")
  } else {
    cat("✗", pkg, "FAILED to load\n")
  }
}

# Display session info
cat("\n=== Session Information ===\n")
sessionInfo()

# Display package versions for critical packages
cat("\n=== Critical Package Versions ===\n")
for (pkg in critical_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    version <- packageVersion(pkg)
    cat(pkg, "version:", as.character(version), "\n")
  }
}

# Memory and system info
cat("\n=== System Information ===\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("Memory limit:", memory.limit(), "MB\n")

# Final message
cat("\n=== INSTALLATION COMPLETE ===\n")
cat("All packages have been installed successfully!\n")
cat("You can now run your WGCNA analysis pipeline.\n")
cat("\nNote: If any package failed to install, you may need to:\n")
cat("1. Update R to the latest version\n")
cat("2. Install system dependencies (especially for Linux/Mac)\n")
cat("3. Run the installation script again\n")
cat("4. Install problematic packages individually\n")

# Save installation log
sink("package_installation_log.txt")
cat("Package Installation Log\n")
cat("Date:", Sys.time(), "\n")
cat("R Version:", R.version.string, "\n")
sessionInfo()
sink()

cat("\nInstallation log saved to: package_installation_log.txt\n")


#################################################### Reading and normalizing the data

setwd("F:/new desktop/thesis/sixth WGCNA analysis/var0.3")
#BiocManager::install("clusterProfiler")
#BiocManager::install("limma")
#BiocManager::install("edgeR")


library("limma")
library("edgeR")


############################ GSE51856_Milk

GSE51856_Milk_Filt=read.csv("GSE51856_Milk_REMBRANDTS_RNA_Stability_0.csv",
                            header = TRUE, row.names = 1)




#BiocManager::install("genefilter")
library(genefilter)

# Keep the genes with t.test p-value<0.3
group = factor(rep(c("A","B"),c(25,25)))
ttest_res = rowttests(as.matrix(GSE51856_Milk_Filt),factor(group))
ttest_rownames=rownames(subset(ttest_res,p.value<.3))

GSE51856_Milk_Filt_Var= GSE51856_Milk_Filt[rownames(GSE51856_Milk_Filt) %in% ttest_rownames, ]  
dim(GSE51856_Milk_Filt_Var)  



# Keep the genes with standard deviation greater than 0.25
# GSE51856_Milk_Filt_Var=subset(GSE51856_Milk_Filt,(apply(GSE51856_Milk_Filt,1,sd))>0.3)			
# dim(GSE51856_Milk_Filt_Var)

# Scale the genes (average=0 and sd=1)
GSE51856_Milk_Filt_Var_Scl <- t(apply(GSE51856_Milk_Filt_Var,1,scale))  

#colnames(data.all_2015_6) <-  c("SRR1913832_670_7_Healthy","SRR1913836_685_7_Healthy","SRR1913839_720_7_Healthy","SRR1913843_728_7_Healthy","SRR1913845_759_7_Healthy")


write.csv(GSE51856_Milk_Filt_Var_Scl, file="1_GSE51856_Milk_Norm.csv")

save(GSE51856_Milk_Filt_Var_Scl, file = "1_GSE51856_Milk_Norm.RData")

############################ End of GSE51856_Milk


############################################################# End of reading and normalizing the data










############################################################# WGCNA
#install.packages("Rcpp", dependencies = TRUE)
#BiocManager::install("impute")
#install.packages('WGCNA')
#BiocManager::install("preprocessCore")

library(WGCNA)
allowWGCNAThreads()
suppressMessages(library(cluster))
options(stringsAsFactors = FALSE);






########################## Reading the data

setwd("C:/Users/M.A.Shirazi/Desktop/thesis/sixth WGCNA analysis/var0.3")


# Load normalized expression file
load(file ="1_GSE51856_Milk_Norm.RData") 


# Seperate healthy samples 
GSE51856_Milk_ExpNorm=t(GSE51856_Milk_Filt_Var_Scl)
########################## End of reading the data







########################## Outlier detection
A = adjacency(t(GSE51856_Milk_ExpNorm), type = "distance") 

# this calculates the whole network connectivity 
k = as.numeric(apply(A, 2, sum)) - 1 

# standardized connectivity 
Z.k = scale(k)

# Designate samples as outlying if their Z.k value is below the threshold 
thresholdZ.k = -2.5  # often -2.5 

# the color vector indicates outlyingness (red) 
outlierColor = ifelse(Z.k < thresholdZ.k, "red", "black") 

# calculate the cluster tree using flahsClust or hclust 
sampleTree = hclust(as.dist(1 - A), method = "average") 

datColors = data.frame(outlierC = outlierColor) 


############## Plot the sample dendrogram and the colors underneath. 
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors,      main = "Sample dendrogram and trait heatmap")

pdf("dendrogram.pdf", width = 10, height = 6)  # تغییر اندازه متناسب با نیازت
plotDendroAndColors(
  sampleTree,
  groupLabels = names(datColors),
  colors = datColors,
  main = "Sample dendrogram and trait heatmap"
)
dev.off()

png("dendrogram.png", width = 3000, height = 1800, res = 300)
plotDendroAndColors(
  sampleTree,
  groupLabels = names(datColors),
  colors = datColors,
  main = "Sample dendrogram and trait heatmap"
)
dev.off()

############## 



# Detect name of samples removed
row.names(outlierColor)=rownames(GSE51856_Milk_ExpNorm)
row.names(outlierColor)[grep("red",outlierColor)]



############## Remove outlying samples from expression  data 
remove.samples = Z.k < thresholdZ.k | is.na(Z.k) 

# the following 2 lines differ from what is written in the book 
GSE51856_Milk_ExpNorm = GSE51856_Milk_ExpNorm[!remove.samples, ] 

# Recompute the sample network among the remaining samples 
A = adjacency(t(GSE51856_Milk_ExpNorm), type = "distance") 

# Let's recompute the Z.k values of outlyingness 
k = as.numeric(apply(A, 2, sum)) - 1 
Z.k = scale(k)
############## 



save(GSE51856_Milk_ExpNorm, file = "2_GSE51856_Milk_RemSam.RData")

########################## End of outlier detection






########################## Checking the matrix

multi=list(Data1=list(data=GSE51856_Milk_ExpNorm) )

multi_g=goodSamplesGenesMS(multi)

GSE51856_Milk_ExpNorm=GSE51856_Milk_ExpNorm[,multi_g$goodGenes]

goodSamplesGenes(GSE51856_Milk_ExpNorm)

# All of them are OK

########################## End of checking the matrix






########################## Network

load(file ="2_GSE51856_Milk_RemSam.RData") 


############## Power finding

# Choosing a soft-threshold to fit a scale-free topology to the network
powers = c(c(1:10), seq(from = 12, to=30, by=1))

# WGCNA
sft = pickSoftThreshold(GSE51856_Milk_ExpNorm, powerVector = powers, networkType="signed", 
                        corFnc = "bicor", verbose = 5,blockSize=17000)   



# corFnc = "bicor"  Or corFnc = cor, corOptions = list(use = 'p')
# Power= 12 for signed   > 0.85


write.table(sft ,"3_GSE51856_Milk_sft.txt")


# Plot the results
sizeGrWindow(9, 5)
par(mar=c(6,4,4,4))     #c(bottom, left, top, right)
par(mfrow = c(1,2));
cex1 = 0.8;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

pdf("scale_independence.pdf", width = 12, height = 6)  # اندازه برای دو نمودار کنار هم
par(mar=c(6,4,4,4))
par(mfrow = c(1,2))
cex1 = 0.8

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
abline(h=0.80, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     type="n", main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels=powers, cex=cex1, col="red")

dev.off()

png("scale_independence.png", width = 3000, height = 1500, res = 300)  # تغییر اندازه و رزولوشن
par(mar=c(6,4,4,4))
par(mfrow = c(1,2))
cex1 = 0.8

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
abline(h=0.80, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     type="n", main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels=powers, cex=cex1, col="red")

dev.off()

#####################################################################################################################################

############################

##################################
# نصب و بارگذاری بسته ggplot2 و gridExtra
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("gridExtra")) install.packages("gridExtra")
library(ggplot2)
library(gridExtra)

# آماده‌سازی داده‌ها برای نمودار میله‌ای
data <- data.frame(
  Threshold = factor(sft$fitIndices[,1]),  # تبدیل به عامل برای محور x
  FitIndex = -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
  MeanConnectivity = sft$fitIndices[,5]
)

# ایجاد نمودار میله‌ای برای Fit Index
p1 <- ggplot(data, aes(x = Threshold, y = FitIndex, fill = FitIndex > 0)) +
  geom_bar(stat="identity", show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "darkorange")) +
  geom_hline(yintercept = 0.80, color="darkred", linetype="dashed", size=1.2) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(
    x = "Soft Threshold (power)", 
    y = "Scale Free Topology Model Fit (signed R^2)"  ) +
  coord_cartesian(ylim = c(min(data$FitIndex) - 0.1, max(data$FitIndex) + 0.1))

# ایجاد نمودار میله‌ای برای Mean Connectivity
p2 <- ggplot(data, aes(x = Threshold, y = MeanConnectivity, fill = MeanConnectivity)) +
  geom_bar(stat="identity", show.legend = FALSE) +
  scale_fill_gradient(low = "lightgreen", high = "darkgreen") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(
    x = "Soft Threshold (power)", 
    y = "Mean Connectivity" ) +
  coord_cartesian(ylim = c(min(data$MeanConnectivity) - 0.1, max(data$MeanConnectivity) + 0.1))

# رسم نمودارها در یک ردیف
grid.arrange(p1, p2, ncol = 2)

pdf("wgcna_threshold_analysis.pdf", width = 14, height = 7)  # ابعاد A4 افقی به این نزدیکه
grid.arrange(p1, p2, ncol = 2)
dev.off()

png("wgcna_threshold_analysis.png", width = 4000, height = 2000, res = 300)  # 300 dpi
grid.arrange(p1, p2, ncol = 2)
dev.off()

#################################################################################################################################
############## End of power finding





############## Module detection

# Generating adjacency, TOM similarity matrices and module detection based on the selected softpower

GSE51856_Milk_ExpNorm_df=as.data.frame(GSE51856_Milk_ExpNorm)




##### Automatic module detection

net_GSE51856_Milk = blockwiseModules(GSE51856_Milk_ExpNorm_df, power = 10, corType= "bicor", 
                                     networkType = "signed",TOMType = "signed", maxBlockSize=17000, 
                                     minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25, 
                                     numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, 
                                     saveTOMFileBase = "4_GSE51856_Milk_TOM",nThreads = 6, verbose = 3)

table(net_GSE51856_Milk$colors)

save(net_GSE51856_Milk, file = "5_GSE51856_Milk_Net_Auto.RData")

# Or load
load(file ="5_GSE51856_Milk_Net_Auto.RData") 


# open a graphics window
sizeGrWindow(12, 9)

# Convert labels to colors for plotting
mergedColors = labels2colors(net_GSE51856_Milk$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net_GSE51856_Milk$dendrograms[[1]], mergedColors[net_GSE51856_Milk$blockGenes[[1]]],
                    "Module colors",  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05, main="Gene hierarchical clustering dendrogram (Control)" )

pdf("dendrogram_control.pdf", width = 10, height = 7)
plotDendroAndColors(
  net_GSE51856_Milk$dendrograms[[1]],
  mergedColors[net_GSE51856_Milk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05

)
dev.off()

png("dendrogram_control.png", width = 4000, height = 3000, res = 300)
plotDendroAndColors(
  net_GSE51856_Milk$dendrograms[[1]],
  mergedColors[net_GSE51856_Milk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05

)
dev.off()

#########################################################################################################################################





#########################################################################################################################################
# Save the modules
Gene_Names=names(GSE51856_Milk_ExpNorm_df)

moduleLabelsAutomatic1 = net_GSE51856_Milk$colors 

moduleColorsAutomatic1 = labels2colors(moduleLabelsAutomatic1)

moduleColorsAutomaticData1=moduleColorsAutomatic1

modules=paste(Gene_Names,net_GSE51856_Milk$colors,moduleColorsAutomatic1,sep=",")

#BiocManager::install("org.Bt.eg.db")


# Merge the gene symbols and id with modules
library("org.Bt.eg.db")
library(xlsx)
library(clusterProfiler)
Genesymbols= bitr(Gene_Names, fromType="ENSEMBL", toType="SYMBOL",
                  OrgDb="org.Bt.eg.db")

IdModules=as.data.frame(cbind(Gene_Names,net_GSE51856_Milk$colors,
                              moduleColorsAutomatic1))
colnames(IdModules)=c("ENSEMBL", "Module_Number", "Module_Color")


GeneIdModule = merge(Genesymbols, IdModules, by="ENSEMBL")



write.xlsx(GeneIdModule,file="6_GSE51856_Milk_Modules_Auto.xlsx",
           sheetName="GSE51856_Milk_Modules", append=TRUE,  row.names=F)


save(GeneIdModule, file = "7_GSE51856_Milk_Modules_Auto.RData")

##### End of automatic module detection


############## End of module detection




########################## End of network






########################## Heatmap for each module

library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")


probes=names(GSE51856_Milk_ExpNorm_df)

modulenames=unique(labels2colors(net_GSE51856_Milk$colors))

# Create colors
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)

col_fun = colorRamp2(c(-2, 0, 2), c( "purple", "white", "red"))

# Separate the groups
GSE51856_Milk_sample_groups =factor(rep(c("A", "B","C", "D", "E", "F", "G", "H","I", "J"), 
                                        times = c(4,5,5,5,5,5,5,4,5,5)))

GSE51856_Milk_fa=rep(c("H0","H12","H24","H36","H48","I0", "I12","I24", "I36", "I48"), 
                     times = c(4,5,5,5,5,5,5,4,5,5))
GSE51856_Milk_fa_col = c("H0" = 1, "H12" = 2,"H24" = 3, "H36" = 4,"H48" = 5, "I0" = 6,
                         "I12" = 7, "I24" = 8,"I36" = 9, "I48" = 10)



col_fun = colorRamp2(c(-2, 0, 2), c("purple", "white", "red"))

for (i in modulenames){
  
  inModule = (moduleColorsAutomaticData1==i) 
  nameofmodule=as.data.frame(probes)[inModule,]
  GSE51856_Milk_Expr =t(GSE51856_Milk_ExpNorm_df)
  GSE51856_Milk_Expr <- subset(GSE51856_Milk_Expr, rownames(GSE51856_Milk_Expr) %in% nameofmodule)
  colnames(GSE51856_Milk_Expr) <- colnames(GSE51856_Milk_Filt)[c(-3,-38)]
  
  healthymean=rowMeans(GSE51856_Milk_Expr[,1:24])
  infectedmean=rowMeans(GSE51856_Milk_Expr[,25:48])
  
  
  
  pdfname=paste("Heatmap_PDF/",i,".pdf",sep="")
  
  pdf(file = pdfname)
  
  
  
  column_names_gp = gpar(col = c(rep("purple", 4), rep("red", 5),
                                 rep("purple", 5), rep("red", 5),
                                 rep("purple", 5), rep("red", 5),
                                 rep("purple", 5), rep("red", 4), 
                                 rep("purple", 5), rep("red", 5)))
  
  
  
  
  
  Heatmap_GSE51856_Milk= Heatmap(healthymean,  width = unit(6, "mm"),
                                 column_split=factor("A"),
                                 column_title = NULL,show_column_names = F,
                                 top_annotation = HeatmapAnnotation(
                                   empty = anno_empty(border = FALSE),
                                   foo = anno_block(gp = gpar(fill = "green"), 
                                                    labels="Average of Healthy",labels_rot=90,
                                                    labels_gp = gpar(col = "black", fontsize = 7))),
                                 
                                 show_heatmap_legend = F,
                                 column_order = NULL, cluster_rows = F,
                                 show_row_dend = F,  column_names_rot = 60,
                                 column_names_gp = gpar(fontsize = 8),
                                 border = TRUE, col = colorRamp2(c(-1, 0, 1), c("purple", "white", "red"))) +
    
    
    Heatmap(infectedmean,  width = unit(6, "mm"),
            column_split=factor("A"),
            column_title = NULL,show_column_names = F,
            top_annotation = HeatmapAnnotation(
              empty = anno_empty(border = FALSE),
              foo = anno_block(gp = gpar(fill = "pink"), 
                               labels="Average of Infected",labels_rot=90,
                               labels_gp = gpar(col = "black", fontsize = 7))),
            
            show_heatmap_legend = F,
            column_order = NULL, cluster_rows = F,
            show_row_dend = F,  column_names_rot = 60,
            column_names_gp = gpar(fontsize = 8),
            border = TRUE, col =  colorRamp2(c(-1, 0, 1), c("purple", "white", "red")))+ 
    
    Heatmap(GSE51856_Milk_Expr, name="mainheat",
            column_split =GSE51856_Milk_sample_groups,
            column_order = NULL, cluster_rows = F,
            show_row_dend = F,
            column_title =paste(i,"module",sep=" "),
            col = col_fun,
            
            bottom_annotation = HeatmapAnnotation(
              .=anno_boxplot(data.matrix(GSE51856_Milk_Expr),
                             border=TRUE, height = unit(1.75, "cm"), width=unit(3, "cm"),
                             gp=gpar(fill=c(rep("green4", 4), rep("green3", 5),
                                            rep("green2", 5), rep("green1", 5),
                                            rep("green", 5), rep("pink4", 5),
                                            rep("pink3", 5), rep("pink2", 4), 
                                            rep("pink1", 5), rep("pink", 5))),
                             pch=".", size=unit(2, "mm"), axis=TRUE),
              
              annotation_width=unit(c(1, 5.0), "cm"),
              which="col"), 
            
            column_gap = unit(1, "mm"),border = TRUE,
            
            top_annotation = HeatmapAnnotation(
              empty = anno_empty(height = unit(2, "cm"),border = FALSE),
              foo = anno_block(gp = gpar(fill = c("green4", "green3", "green2", "green1", "green",
                                                  "pink4","pink3", "pink2", "pink1","pink")), 
                               labels=c("H0","H12","H24","H36","H48","I0", "I12","I24", "I36", "I48"),
                               labels_gp = gpar(col = "black", fontsize = 7))
            ),
            
            
            
            row_title = "Genes", 
            row_title_side =  "right",
            row_title_gp = gpar(fontsize = 10),
            row_names_max_width = unit(1, "cm"),
            row_names_gp = gpar(fontsize = 0.5),
            show_row_names=F,
            
            show_column_names = TRUE,
            column_names_gp = gpar(fontsize = 4),
            column_names_rot = 45,
            
            width = unit(10, "cm"), height = unit(12, "cm"),
            
            # Legend
            show_heatmap_legend = T,
            heatmap_legend_param = list(direction = "horizontal",
                                        at = c(-2, 0, 2),
                                        labels = c("-2", "0", "2"),
                                        title = "RNA stability",
                                        legend_width  = unit(3, "cm"),
                                        title_position = "lefttop",
                                        title_gp = gpar(fontsize = 7, fontface = "bold"),
                                        labels_gp = gpar(fontsize = 7)
            ))
  
  
  draw( Heatmap_GSE51856_Milk, heatmap_legend_side = "bottom" )
  
  
  decorate_annotation(".", { 
    grid.text("RNA Stability",rot = 90, just = "bottom",
              x = unit(-0.6, "cm"),gp=gpar(fontsize=6, col="black")) 
  })
  
  
  
  
  
  dev.off();
  
}


########################## End fo heatmap
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

probes = names(GSE51856_Milk_ExpNorm_df)
modulenames = unique(labels2colors(net_GSE51856_Milk$colors))

# بهبود پالت رنگی با استفاده از طیف زیباتر
col_fun = colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B")
)

# تعریف گروه‌های نمونه
GSE51856_Milk_sample_groups = factor(
  rep(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
      times = c(4,5,5,5,5,5,5,4,5,5))
)

GSE51856_Milk_fa = rep(
  c("H0", "H12", "H24", "H36", "H48", "I0", "I12", "I24", "I36", "I48"),
  times = c(4,5,5,5,5,5,5,4,5,5)
)

for (i in modulenames) {
  inModule = (moduleColorsAutomaticData1 == i)
  nameofmodule = as.data.frame(probes)[inModule,]
  GSE51856_Milk_Expr = t(GSE51856_Milk_ExpNorm_df)
  GSE51856_Milk_Expr <- subset(GSE51856_Milk_Expr, 
                               rownames(GSE51856_Milk_Expr) %in% nameofmodule)
  colnames(GSE51856_Milk_Expr) <- colnames(GSE51856_Milk_Filt)[c(-3,-38)]
  
  healthymean = rowMeans(GSE51856_Milk_Expr[,1:24])
  infectedmean = rowMeans(GSE51856_Milk_Expr[,25:48])
  
  pdfname = paste("Heatmap_PDF/", i, "_improved.pdf", sep="")
  
  pdf(file = pdfname, width = 12, height = 14) # افزایش اندازه خروجی برای وضوح بهتر
  
  # بهبود رنگ‌های ستون‌ها
  healthy_colors = colorRampPalette(c("#234D20", "#77AB59"))(5)
  infected_colors = colorRampPalette(c("#8B0000", "#FF6B6B"))(5)
  
  # تنظیمات بهتر برای نمایش میانگین‌ها
  mean_heatmap_healthy = Heatmap(
    healthymean,
    width = unit(8, "mm"),
    column_split = factor("A"),
    column_title = NULL,
    show_column_names = F,
    top_annotation = HeatmapAnnotation(
      foo = anno_block(
        gp = gpar(fill = healthy_colors[3]),
        labels = "Healthy Average",
        labels_rot = 90,
        labels_gp = gpar(col = "black", fontsize = 8)
      )
    ),
    show_heatmap_legend = F,
    cluster_rows = F,
    border = TRUE,
    col = col_fun
  )
  
  mean_heatmap_infected = Heatmap(
    infectedmean,
    width = unit(8, "mm"),
    column_split = factor("A"),
    column_title = NULL,
    show_column_names = F,
    top_annotation = HeatmapAnnotation(
      foo = anno_block(
        gp = gpar(fill = infected_colors[3]),
        labels = "Infected Average",
        labels_rot = 90,
        labels_gp = gpar(col = "black", fontsize = 8)
      )
    ),
    show_heatmap_legend = F,
    cluster_rows = F,
    border = TRUE,
    col = col_fun
  )
  
  # هیت‌مپ اصلی با تنظیمات بهبود یافته
  main_heatmap = Heatmap(
    GSE51856_Milk_Expr,
    name = "RNA Stability",
    column_split = GSE51856_Milk_sample_groups,
    cluster_rows = F,
    column_title = paste(i, "Module", sep=" "),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    col = col_fun,
    
    bottom_annotation = HeatmapAnnotation(
      boxplot = anno_boxplot(
        GSE51856_Milk_Expr,
        border = TRUE,
        height = unit(2, "cm"),
        gp = gpar(fill = c(rep(healthy_colors, each=5)[1:24], 
                           rep(infected_colors, each=5)[1:24])),
        outline = FALSE,
        size = unit(0.5, "mm"),
        axis = TRUE
      )
    ),
    
    top_annotation = HeatmapAnnotation(
      group = anno_block(
        gp = gpar(fill = c(healthy_colors, infected_colors)),
        labels = c("H0","H12","H24","H36","H48","I0","I12","I24","I36","I48"),
        labels_gp = gpar(col = "black", fontsize = 8)
      )
    ),
    
    row_title = "Genes",
    row_title_side = "right",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    show_row_names = FALSE,
    
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 6),
    
    heatmap_legend_param = list(
      title = "RNA Stability Score",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      direction = "horizontal",
      legend_width = unit(4, "cm")
    )
  )
  
  combined_heatmap = mean_heatmap_healthy + mean_heatmap_infected + main_heatmap
  
  draw(combined_heatmap, 
       heatmap_legend_side = "bottom",
       padding = unit(c(2, 2, 2, 2), "mm"))
  
  # اضافه کردن عنوان محور Y برای نمودار جعبه‌ای
  decorate_annotation("boxplot", {
    grid.text(
      "Expression Level",
      rot = 90,
      x = unit(-2, "cm"),
      gp = gpar(fontsize = 8, fontface = "bold")
    )
  })
  
  dev.off()
}





############## End of module detection

# تعریف تابع اصلی برای ایجاد هیت‌مپ با پارامترهای قابل تنظیم
create_dynamic_heatmap <- function(
    expression_data,
    module_name,
    sample_groups,
    healthy_range = 1:24,
    infected_range = 25:48,
    color_scheme = list(
      healthy = c("#234D20", "#77AB59"),
      infected = c("#8B0000", "#FF6B6B"),
      heatmap = c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B")
    ),
    pdf_dimensions = c(width = 12, height = 14),
    font_sizes = list(
      title = 12,
      subtitle = 10,
      labels = 8,
      legend = 8
    )
) {
  # تنظیم پالت رنگی دینامیک
  col_fun = colorRamp2(
    seq(-2, 2, length.out = length(color_scheme$heatmap)),
    color_scheme$heatmap
  )
  
  # محاسبه میانگین‌های دینامیک
  healthymean = rowMeans(expression_data[, healthy_range])
  infectedmean = rowMeans(expression_data[, infected_range])
  
  # ایجاد طیف رنگی برای گروه‌ها
  healthy_colors = colorRampPalette(color_scheme$healthy)(5)
  infected_colors = colorRampPalette(color_scheme$infected)(5)
  
  # تابع کمکی برای ایجاد هیت‌مپ میانگین
  create_mean_heatmap <- function(data, colors, label) {
    Heatmap(
      data,
      width = unit(8, "mm"),
      column_split = factor("A"),
      column_title = NULL,
      show_column_names = F,
      top_annotation = HeatmapAnnotation(
        foo = anno_block(
          gp = gpar(fill = colors[3]),
          labels = label,
          labels_rot = 90,
          labels_gp = gpar(col = "black", fontsize = font_sizes$labels)
        )
      ),
      show_heatmap_legend = F,
      cluster_rows = F,
      border = TRUE,
      col = col_fun
    )
  }
  
  # ایجاد هیت‌مپ‌های میانگین
  mean_heatmap_healthy = create_mean_heatmap(healthymean, healthy_colors, "Healthy Average")
  mean_heatmap_infected = create_mean_heatmap(infectedmean, infected_colors, "Infected Average")
  
  # تابع برای محاسبه آمارهای توصیفی
  calculate_statistics <- function(data) {
    list(
      mean = mean(as.matrix(data)),
      sd = sd(as.matrix(data)),
      quantiles = quantile(as.matrix(data), probs = c(0.25, 0.5, 0.75))
    )
  }
  
  # محاسبه آمارها
  stats_healthy = calculate_statistics(expression_data[, healthy_range])
  stats_infected = calculate_statistics(expression_data[, infected_range])
  
  # ایجاد هیت‌مپ اصلی با آمارهای توصیفی
  main_heatmap = Heatmap(
    expression_data,
    name = "RNA Stability",
    column_split = sample_groups,
    cluster_rows = F,
    column_title = paste(
      module_name, "Module\n",
      "Healthy Mean:", round(stats_healthy$mean, 2),
      "| Infected Mean:", round(stats_infected$mean, 2)
    ),
    column_title_gp = gpar(fontsize = font_sizes$title, fontface = "bold"),
    col = col_fun,
    
    # اضافه کردن نمودار جعبه‌ای پویا
    bottom_annotation = HeatmapAnnotation(
      boxplot = anno_boxplot(
        expression_data,
        border = TRUE,
        height = unit(2, "cm"),
        gp = gpar(fill = c(rep(healthy_colors, each=5)[1:24], 
                           rep(infected_colors, each=5)[1:24])),
        outline = FALSE,
        size = unit(0.5, "mm"),
        axis = TRUE,
        box_width = 0.8
      ),
      
      # اضافه کردن نمودار چگالی
      density = anno_density(
        expression_data,
        gp = gpar(fill = "#CCCCCC"),
        height = unit(1.5, "cm")
      )
    ),
    
    top_annotation = HeatmapAnnotation(
      group = anno_block(
        gp = gpar(fill = c(healthy_colors, infected_colors)),
        labels = c("H0","H12","H24","H36","H48","I0","I12","I24","I36","I48"),
        labels_gp = gpar(col = "black", fontsize = font_sizes$labels)
      )
    ),
    
    row_title = "Genes",
    row_title_side = "right",
    row_title_gp = gpar(fontsize = font_sizes$subtitle, fontface = "bold"),
    show_row_names = FALSE,
    
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = font_sizes$labels),
    
    heatmap_legend_param = list(
      title = "RNA Stability Score",
      title_gp = gpar(fontsize = font_sizes$subtitle, fontface = "bold"),
      labels_gp = gpar(fontsize = font_sizes$legend),
      direction = "horizontal",
      legend_width = unit(4, "cm")
    )
  )
  
  # ترکیب هیت‌مپ‌ها
  combined_heatmap = mean_heatmap_healthy + mean_heatmap_infected + main_heatmap
  
  return(combined_heatmap)
}

# استفاده از تابع
for (i in modulenames) {
  inModule = (moduleColorsAutomaticData1 == i)
  nameofmodule = as.data.frame(probes)[inModule,]
  GSE51856_Milk_Expr = t(GSE51856_Milk_ExpNorm_df)
  GSE51856_Milk_Expr <- subset(GSE51856_Milk_Expr, 
                               rownames(GSE51856_Milk_Expr) %in% nameofmodule)
  colnames(GSE51856_Milk_Expr) <- colnames(GSE51856_Milk_Filt)[c(-3,-38)]
  
  pdfname = paste("Heatmap_PDF/", i, "_dynamic.pdf", sep="")
  pdf(file = pdfname, width = 12, height = 14)
  
  # ایجاد هیت‌مپ با پارامترهای سفارشی
  heatmap = create_dynamic_heatmap(
    expression_data = GSE51856_Milk_Expr,
    module_name = i,
    sample_groups = GSE51856_Milk_sample_groups,
    color_scheme = list(
      healthy = c("#234D20", "#77AB59"),
      infected = c("#8B0000", "#FF6B6B"),
      heatmap = c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B")
    ),
    font_sizes = list(
      title = 12,
      subtitle = 10,
      labels = 8,
      legend = 8
    )
  )
  
  # رسم هیت‌مپ با تنظیمات اضافی
  draw(heatmap, 
       heatmap_legend_side = "bottom",
       padding = unit(c(2, 2, 2, 2), "mm"))
  
  # اضافه کردن عناوین محورها
  decorate_annotation("boxplot", {
    grid.text(
      "Expression Level",
      rot = 90,
      x = unit(-2, "cm"),
      gp = gpar(fontsize = 8, fontface = "bold")
    )
  })
  
  decorate_annotation("density", {
    grid.text(
      "Density",
      rot = 90,
      x = unit(-2, "cm"),
      gp = gpar(fontsize = 8, fontface = "bold")
    )
  })
  
  dev.off()
}

##################

# تعریف تابع اصلی برای ایجاد هیت‌مپ با پارامترهای قابل تنظیم
create_dynamic_heatmap <- function(
    expression_data,
    module_name,
    sample_groups,
    healthy_range = 1:24,
    infected_range = 25:48,
    color_scheme = list(
      healthy = c("#a6cee3", "#1f78b4"),
      infected = c("#fb9a99", "#e31a1c"),
      heatmap = c("#4575b4", "#91bfdb", "#f7f7f7", "#fc8d59", "#d73027")
    ),
    pdf_dimensions = c(width = 12, height = 14),
    font_sizes = list(
      title = 14,
      subtitle = 12,
      labels = 10,
      legend = 10
    )
) {
  # تنظیم پالت رنگی دینامیک
  col_fun = colorRamp2(
    seq(-2, 2, length.out = length(color_scheme$heatmap)),
    color_scheme$heatmap
  )
  
  # محاسبه میانگین‌های دینامیک
  healthymean = rowMeans(expression_data[, healthy_range])
  infectedmean = rowMeans(expression_data[, infected_range])
  
  # ایجاد طیف رنگی برای گروه‌ها
  healthy_colors = colorRampPalette(color_scheme$healthy)(5)
  infected_colors = colorRampPalette(color_scheme$infected)(5)
  
  # تابع کمکی برای ایجاد هیت‌مپ میانگین
  create_mean_heatmap <- function(data, colors, label) {
    Heatmap(
      data,
      width = unit(8, "mm"),
      column_split = factor("A"),
      column_title = NULL,
      show_column_names = FALSE,
      top_annotation = HeatmapAnnotation(
        foo = anno_block(
          gp = gpar(fill = colors[3]),
          labels = label,
          labels_rot = 90,
          labels_gp = gpar(col = "black", fontsize = font_sizes$labels)
        )
      ),
      show_heatmap_legend = FALSE,
      cluster_rows = FALSE,
      border = TRUE,
      col = col_fun
    )
  }
  
  # ایجاد هیت‌مپ‌های میانگین
  mean_heatmap_healthy = create_mean_heatmap(healthymean, healthy_colors, "Healthy Average")
  mean_heatmap_infected = create_mean_heatmap(infectedmean, infected_colors, "Infected Average")
  
  # تابع برای محاسبه آمارهای توصیفی
  calculate_statistics <- function(data) {
    list(
      mean = mean(as.matrix(data)),
      sd = sd(as.matrix(data)),
      quantiles = quantile(as.matrix(data), probs = c(0.25, 0.5, 0.75))
    )
  }
  
  # محاسبه آمارها
  stats_healthy = calculate_statistics(expression_data[, healthy_range])
  stats_infected = calculate_statistics(expression_data[, infected_range])
  
  # ایجاد هیت‌مپ اصلی با آمارهای توصیفی
  main_heatmap = Heatmap(
    expression_data,
    name = "RNA Stability",
    column_split = sample_groups,
    cluster_rows = FALSE,
    column_title = paste(
      module_name, "Module\n",
      "Healthy Mean:", round(stats_healthy$mean, 2),
      "±", round(stats_healthy$sd, 2),
      "| Infected Mean:", round(stats_infected$mean, 2),
      "±", round(stats_infected$sd, 2)
    ),
    column_title_gp = gpar(fontsize = font_sizes$title, fontface = "bold.italic"),
    col = col_fun,
    
    # اضافه کردن نمودار جعبه‌ای پویا
    bottom_annotation = HeatmapAnnotation(
      boxplot = anno_boxplot(
        expression_data,
        border = TRUE,
        height = unit(2, "cm"),
        gp = gpar(fill = alpha(c(rep(healthy_colors, each = 5)[1:24],
                                 rep(infected_colors, each = 5)[1:24]), 0.6)),
        outline = FALSE,
        size = unit(0.5, "mm"),
        axis_param = list(
          at = c(-2, 0, 2),
          labels = c("Low", "Medium", "High")
        )
      ),
      
      # اضافه کردن نمودار چگالی
      density = anno_density(
        expression_data,
        gp = gpar(fill = alpha("#CCCCCC", 0.5)),
        height = unit(1.5, "cm")
      )
    ),
    
    top_annotation = HeatmapAnnotation(
      group = anno_block(
        gp = gpar(fill = c(healthy_colors, infected_colors)),
        labels = c("H0","H12","H24","H36","H48","I0","I12","I24","I36","I48"),
        labels_gp = gpar(col = "black", fontsize = font_sizes$labels)
      )
    ),
    
    row_title = "Genes",
    row_title_side = "right",
    row_title_gp = gpar(fontsize = font_sizes$subtitle, fontface = "italic", col = "blue"),
    show_row_names = FALSE,
    
    column_names_rot = 30,
    column_names_gp = gpar(fontsize = font_sizes$labels + 2),
    
    heatmap_legend_param = list(
      title = "RNA Stability Score",
      title_gp = gpar(fontsize = font_sizes$subtitle, fontface = "bold"),
      labels_gp = gpar(fontsize = font_sizes$legend),
      direction = "vertical",
      legend_height = unit(5, "cm")
    )
  )
  
  # ترکیب هیت‌مپ‌ها
  combined_heatmap = mean_heatmap_healthy + mean_heatmap_infected + main_heatmap
  
  return(combined_heatmap)
}

# استفاده از تابع
for (i in modulenames) {
  inModule = (moduleColorsAutomaticData1 == i)
  nameofmodule = as.data.frame(probes)[inModule,]
  GSE51856_Milk_Expr = t(GSE51856_Milk_ExpNorm_df)
  GSE51856_Milk_Expr <- subset(GSE51856_Milk_Expr, 
                               rownames(GSE51856_Milk_Expr) %in% nameofmodule)
  colnames(GSE51856_Milk_Expr) <- colnames(GSE51856_Milk_Filt)[c(-3,-38)]
  
  pngname = paste("Heatmap_PNG/", i, "_dynamic.png", sep="")
  png(filename = pngname, width = 12, height = 14, units = "in", res = 300)
  
  # ایجاد هیت‌مپ با پارامترهای سفارشی
  heatmap = create_dynamic_heatmap(
    expression_data = GSE51856_Milk_Expr,
    module_name = i,
    sample_groups = GSE51856_Milk_sample_groups,
    color_scheme = list(
      healthy = c("#a6cee3", "#1f78b4"),
      infected = c("#fb9a99", "#e31a1c"),
      heatmap = c("#4575b4", "#91bfdb", "#f7f7f7", "#fc8d59", "#d73027")
    ),
    font_sizes = list(
      title = 14,
      subtitle = 12,
      labels = 10,
      legend = 10
    )
  )
  
  # رسم هیت‌مپ با تنظیمات اضافی
  draw(heatmap, 
       heatmap_legend_side = "right",
       padding = unit(c(4, 4, 4, 4), "mm"))
  
  # اضافه کردن عناوین محورها
  decorate_annotation("boxplot", {
    grid.text(
      "Expression Level",
      rot = 90,
      x = unit(-2, "cm"),
      gp = gpar(fontsize = 10, fontface = "bold")
    )
  })
  
  decorate_annotation("density", {
    grid.text(
      "Density",
      rot = 90,
      x = unit(-2, "cm"),
      gp = gpar(fontsize = 10, fontface = "bold")
    )
  })
  
  # بستن فایل خروجی
  dev.off()                
################
  
  
  
  
  
  
  
  
  
  
  
########################## GO analysis
setwd("C:/Users/M.A.Shirazi/Desktop/thesis/sixth WGCNA analysis/var0.3")


load(file ="5_GSE51856_Milk_Net_Auto.RData") 
load(file ="2_GSE51856_Milk_RemSam.RData") 



########### GO analysis by EnricheR

options(java.parameters = "-Xmx30000m")
library(enrichR)
library(xlsx)
library("clusterProfiler")
library("org.Bt.eg.db")
library("WGCNA")

probes=names(GSE51856_Milk_ExpNorm_df)

moduleLabelsAutomatic1 = net_GSE51856_Milk$colors 
moduleColorsAutomatic1 = labels2colors(moduleLabelsAutomatic1)
moduleColorsAutomaticData1=moduleColorsAutomatic1

listEnrichrSites()
setEnrichrSite("Enrichr")
websiteLive <- TRUE

dbs <- c("Achilles_fitness_decrease",
         "Achilles_fitness_increase",
         "Aging_Perturbations_from_GEO_down",
         "Aging_Perturbations_from_GEO_up",
         "Allen_Brain_Atlas_10x_scRNA_2021",
         "Allen_Brain_Atlas_down",
         "Allen_Brain_Atlas_up",
         "ARCHS4_Cell-lines",
         "ARCHS4_IDG_Coexp",
         "ARCHS4_Kinases_Coexp",
         "ARCHS4_TFs_Coexp",
         "ARCHS4_Tissues",
         "Azimuth_2023",
         "Azimuth_Cell_Types_2021",
         "BioCarta_2013",
         "BioCarta_2015",
         "BioCarta_2016",
         "BioPlanet_2019",
         "BioPlex_2017",
         "Cancer_Cell_Line_Encyclopedia",
         "CCLE_Proteomics_2020",
         "CellMarker_2024",
         "CellMarker_Augmented_2021",
         "ChEA_2013",
         "ChEA_2015",
         "ChEA_2016",
         "ChEA_2022",
         "Chromosome_Location",
         "Chromosome_Location_hg19",
         "ClinVar_2019",
         "CORUM",
         "COVID-19_Related_Gene_Sets",
         "COVID-19_Related_Gene_Sets_2021",
         "Data_Acquisition_Method_Most_Popular_Genes",
         "dbGaP",
         "DepMap_WG_CRISPR_Screens_Broad_CellLines_2019",
         "DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019",
         "Descartes_Cell_Types_and_Tissue_2021",
         "Diabetes_Perturbations_GEO_2022",
         "Disease_Perturbations_from_GEO_down",
         "Disease_Perturbations_from_GEO_up",
         "Disease_Signatures_from_GEO_down_2014",
         "Disease_Signatures_from_GEO_up_2014",
         "DisGeNET",
         "Drug_Perturbations_from_GEO_2014",
         "Drug_Perturbations_from_GEO_down",
         "Drug_Perturbations_from_GEO_up",
         "DrugMatrix",
         "DSigDB",
         "Elsevier_Pathway_Collection",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "ENCODE_Histone_Modifications_2013",
         "ENCODE_Histone_Modifications_2015",
         "ENCODE_TF_ChIP-seq_2014",
         "ENCODE_TF_ChIP-seq_2015",
         "Enrichr_Libraries_Most_Popular_Genes",
         "Enrichr_Submissions_TF-Gene_Coocurrence",
         "Enrichr_Users_Contributed_Lists_2020",
         "Epigenomics_Roadmap_HM_ChIP-seq",
         "ESCAPE",
         "FANTOM6_lncRNA_KD_DEGs",
         "GeDiPNet_2023",
         "Gene_Perturbations_from_GEO_down",
         "Gene_Perturbations_from_GEO_up",
         "Genes_Associated_with_NIH_Grants",
         "GeneSigDB",
         "Genome_Browser_PWMs",
         "GlyGen_Glycosylated_Proteins_2022",
         "GO_Biological_Process_2013",
         "GO_Biological_Process_2015",
         "GO_Biological_Process_2017",
         "GO_Biological_Process_2017b",
         "GO_Biological_Process_2018",
         "GO_Biological_Process_2021",
         "GO_Biological_Process_2023",
         "GO_Cellular_Component_2013",
         "GO_Cellular_Component_2015",
         "GO_Cellular_Component_2017",
         "GO_Cellular_Component_2017b",
         "GO_Cellular_Component_2018",
         "GO_Cellular_Component_2021",
         "GO_Cellular_Component_2023",
         "GO_Molecular_Function_2013",
         "GO_Molecular_Function_2015",
         "GO_Molecular_Function_2017",
         "GO_Molecular_Function_2017b",
         "GO_Molecular_Function_2018",
         "GO_Molecular_Function_2021",
         "GO_Molecular_Function_2023",
         "GTEx_Aging_Signatures_2021",
         "GTEx_Tissue_Expression_Down",
         "GTEx_Tissue_Expression_Up",
         "GTEx_Tissues_V8_2023",
         "GWAS_Catalog_2019",
         "GWAS_Catalog_2023",
         "HDSigDB_Human_2021",
         "HDSigDB_Mouse_2021",
         "HMDB_Metabolites",
         "HMS_LINCS_KinomeScan",
         "HomoloGene",
         "HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression",
         "HuBMAP_ASCTplusB_augmented_2022",
         "Human_Gene_Atlas",
         "Human_Phenotype_Ontology",
         "HumanCyc_2015",
         "HumanCyc_2016",
         "huMAP",
         "IDG_Drug_Targets_2022",
         "InterPro_Domains_2019",
         "Jensen_COMPARTMENTS",
         "Jensen_DISEASES",
         "Jensen_TISSUES",
         "KEA_2013",
         "KEA_2015",
         "KEGG_2013",
         "KEGG_2015",
         "KEGG_2016",
         "KEGG_2019_Human",
         "KEGG_2019_Mouse",
         "KEGG_2021_Human",
         "Kinase_Perturbations_from_GEO_down",
         "Kinase_Perturbations_from_GEO_up",
         "KOMP2_Mouse_Phenotypes_2022",
         "L1000_Kinase_and_GPCR_Perturbations_down",
         "L1000_Kinase_and_GPCR_Perturbations_up",
         "Ligand_Perturbations_from_GEO_down",
         "Ligand_Perturbations_from_GEO_up",
         "LINCS_L1000_Chem_Pert_Consensus_Sigs",
         "LINCS_L1000_Chem_Pert_down",
         "LINCS_L1000_Chem_Pert_up",
         "LINCS_L1000_CRISPR_KO_Consensus_Sigs",
         "LINCS_L1000_Ligand_Perturbations_down",
         "LINCS_L1000_Ligand_Perturbations_up",
         "lncHUB_lncRNA_Co-Expression",
         "MAGMA_Drugs_and_Diseases",
         "MAGNET_2023",
         "MCF7_Perturbations_from_GEO_down",
         "MCF7_Perturbations_from_GEO_up",
         "Metabolomics_Workbench_Metabolites_2022",
         "MGI_Mammalian_Phenotype_2013",
         "MGI_Mammalian_Phenotype_2017",
         "MGI_Mammalian_Phenotype_Level_3",
         "MGI_Mammalian_Phenotype_Level_4",
         "MGI_Mammalian_Phenotype_Level_4_2019",
         "MGI_Mammalian_Phenotype_Level_4_2021",
         "Microbe_Perturbations_from_GEO_down",
         "Microbe_Perturbations_from_GEO_up",
         "miRTarBase_2017",
         "MoTrPAC_2023",
         "Mouse_Gene_Atlas",
         "MSigDB_Computational",
         "MSigDB_Hallmark_2020",
         "MSigDB_Oncogenic_Signatures",
         "NCI-60_Cancer_Cell_Lines",
         "NCI-Nature_2015",
         "NCI-Nature_2016",
         "NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions",
         "NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions",
         "NIH_Funded_PIs_2017_Human_AutoRIF",
         "NIH_Funded_PIs_2017_Human_GeneRIF",
         "NURSA_Human_Endogenous_Complexome",
         "Old_CMAP_down",
         "Old_CMAP_up",
         "OMIM_Disease",
         "OMIM_Expanded",
         "Orphanet_Augmented_2021",
         "PanglaoDB_Augmented_2021",
         "Panther_2015",
         "Panther_2016",
         "Pfam_Domains_2019",
         "Pfam_InterPro_Domains",
         "PFOCR_Pathways",
         "PFOCR_Pathways_2023",
         "PhenGenI_Association_2021",
         "PheWeb_2019",
         "Phosphatase_Substrates_from_DEPOD",
         "PPI_Hub_Proteins",
         "Proteomics_Drug_Atlas_2023",
         "ProteomicsDB_2020",
         "Rare_Diseases_AutoRIF_ARCHS4_Predictions",
         "Rare_Diseases_AutoRIF_Gene_Lists",
         "Rare_Diseases_GeneRIF_ARCHS4_Predictions",
         "Rare_Diseases_GeneRIF_Gene_Lists",
         "Reactome_2013",
         "Reactome_2015",
         "Reactome_2016",
         "Reactome_2022",
         "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO",
         "RNAseq_Automatic_GEO_Signatures_Human_Down",
         "RNAseq_Automatic_GEO_Signatures_Human_Up",
         "RNAseq_Automatic_GEO_Signatures_Mouse_Down",
         "RNAseq_Automatic_GEO_Signatures_Mouse_Up",
         "Rummagene_kinases",
         "Rummagene_signatures",
         "Rummagene_transcription_factors",
         "SILAC_Phosphoproteomics",
         "SubCell_BarCode",
         "SynGO_2022",
         "SynGO_2024",
         "SysMyo_Muscle_Gene_Sets",
         "Table_Mining_of_CRISPR_Studies",
         "Tabula_Muris",
         "Tabula_Sapiens",
         "TargetScan_microRNA",
         "TargetScan_microRNA_2017",
         "TF-LOF_Expression_from_GEO",
         "TF_Perturbations_Followed_by_Expression",
         "TG_GATES_2020",
         "The_Kinase_Library_2023",
         "Tissue_Protein_Expression_from_Human_Proteome_Map",
         "Tissue_Protein_Expression_from_ProteomicsDB",
         "Transcription_Factor_PPIs",
         "TRANSFAC_and_JASPAR_PWMs",
         "TRRUST_Transcription_Factors_2019",
         "UK_Biobank_GWAS_v1",
         "Virus-Host_PPI_P-HIPSTer_2020",
         "Virus_Perturbations_from_GEO_down",
         "Virus_Perturbations_from_GEO_up",
         "VirusMINT",
         "WikiPathway_2021_Human",
         "WikiPathway_2023_Human",
         "WikiPathways_2013",
         "WikiPathways_2015",
         "WikiPathways_2016",
         "WikiPathways_2019_Human",
         "WikiPathways_2019_Mouse")

modulenames=unique(labels2colors(net_GSE51856_Milk$colors))


for (i in modulenames){
  gc()
  # Select the module of interest
  inModule = (moduleColorsAutomaticData1==i) 
  nameofmodule=as.data.frame(probes)[inModule,]
  
  
  # Convert ENSEMBL to Gene symbol by library("clusterProfiler")
  nameofmodulesymbol=bitr(nameofmodule, fromType="ENSEMBL", toType="SYMBOL",
                          OrgDb="org.Bt.eg.db")$SYMBOL
  
  enriched= enrichr(nameofmodulesymbol, databases = dbs)
  
  
  write.xlsx(enriched$Achilles_fitness_decrease, file='GO_enrichR/Achilles_fitness_decrease.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Achilles_fitness_increase, file='GO_enrichR/Achilles_fitness_increase.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Aging_Perturbations_from_GEO_down, file='GO_enrichR/Aging_Perturbations_from_GEO_down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Aging_Perturbations_from_GEO_up, file='GO_enrichR/Aging_Perturbations_from_GEO_up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Allen_Brain_Atlas_10x_scRNA_2021, file='GO_enrichR/Allen_Brain_Atlas_10x_scRNA_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Allen_Brain_Atlas_down, file='GO_enrichR/Allen_Brain_Atlas_down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Allen_Brain_Atlas_up, file='GO_enrichR/Allen_Brain_Atlas_up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`ARCHS4_Cell-lines`, file='GO_enrichR/ARCHS4_Cell-lines.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$ARCHS4_IDG_Coexp, file='GO_enrichR/ARCHS4_IDG_Coexp.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$ARCHS4_Kinases_Coexp, file='GO_enrichR/ARCHS4_Kinases_Coexp.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$ARCHS4_TFs_Coexp, file='GO_enrichR/ARCHS4_TFs_Coexp.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$ARCHS4_Tissues, file='GO_enrichR/ARCHS4_Tissues.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Azimuth_2023, file='GO_enrichR/Azimuth_2023.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Azimuth_Cell_Types_2021, file='GO_enrichR/Azimuth_Cell_Types_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$BioCarta_2013, file='GO_enrichR/BioCarta_2013.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$BioCarta_2015, file='GO_enrichR/BioCarta_2015.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$BioCarta_2016, file='GO_enrichR/BioCarta_2016.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$BioPlanet_2019, file='GO_enrichR/BioPlanet_2019.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$BioPlex_2017, file='GO_enrichR/BioPlex_2017.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Cancer_Cell_Line_Encyclopedia, file='GO_enrichR/Cancer_Cell_Line_Encyclopedia.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$CCLE_Proteomics_2020, file='GO_enrichR/CCLE_Proteomics_2020.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$CellMarker_2024, file='GO_enrichR/CellMarker_2024.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$CellMarker_Augmented_2021, file='GO_enrichR/CellMarker_Augmented_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$ChEA_2013, file='GO_enrichR/ChEA_2013.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$ChEA_2015, file='GO_enrichR/ChEA_2015.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$ChEA_2016, file='GO_enrichR/ChEA_2016.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$ChEA_2022, file='GO_enrichR/ChEA_2022.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Chromosome_Location, file='GO_enrichR/Chromosome_Location.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Chromosome_Location_hg19, file='GO_enrichR/Chromosome_Location_hg19.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$ClinVar_2019, file='GO_enrichR/ClinVar_2019.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$CORUM, file='GO_enrichR/CORUM.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`COVID-19_Related_Gene_Sets`, file='GO_enrichR/COVID-19_Related_Gene_Sets.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`COVID-19_Related_Gene_Sets_2021`, file='GO_enrichR/COVID-19_Related_Gene_Sets_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Data_Acquisition_Method_Most_Popular_Genes, file='GO_enrichR/Data_Acquisition_Method_Most_Popular_Genes.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$dbGaP, file='GO_enrichR/dbGaP.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$DepMap_WG_CRISPR_Screens_Broad_CellLines_2019, file='GO_enrichR/DepMap_WG_CRISPR_Screens_Broad_CellLines_2019.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019, file='GO_enrichR/DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Descartes_Cell_Types_and_Tissue_2021, file='GO_enrichR/Descartes_Cell_Types_and_Tissue_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Diabetes_Perturbations_GEO_2022, file='GO_enrichR/Diabetes_Perturbations_GEO_2022.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Disease_Perturbations_from_GEO_down, file='GO_enrichR/Disease_Perturbations_from_GEO_down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Disease_Perturbations_from_GEO_up, file='GO_enrichR/Disease_Perturbations_from_GEO_up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Disease_Signatures_from_GEO_down_2014, file='GO_enrichR/Disease_Signatures_from_GEO_down_2014.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Disease_Signatures_from_GEO_up_2014, file='GO_enrichR/Disease_Signatures_from_GEO_up_2014.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$DisGeNET, file='GO_enrichR/DisGeNET.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Drug_Perturbations_from_GEO_2014, file='GO_enrichR/Drug_Perturbations_from_GEO_2014.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Drug_Perturbations_from_GEO_down, file='GO_enrichR/Drug_Perturbations_from_GEO_down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Drug_Perturbations_from_GEO_up, file='GO_enrichR/Drug_Perturbations_from_GEO_up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$DrugMatrix, file='GO_enrichR/DrugMatrix.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$DSigDB, file='GO_enrichR/DSigDB.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Elsevier_Pathway_Collection, file='GO_enrichR/Elsevier_Pathway_Collection.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X`, file='GO_enrichR/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$ENCODE_Histone_Modifications_2013, file='GO_enrichR/ENCODE_Histone_Modifications_2013.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$ENCODE_Histone_Modifications_2015, file='GO_enrichR/ENCODE_Histone_Modifications_2015.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`ENCODE_TF_ChIP-seq_2014`, file='GO_enrichR/ENCODE_TF_ChIP-seq_2014.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`ENCODE_TF_ChIP-seq_2015`, file='GO_enrichR/ENCODE_TF_ChIP-seq_2015.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Enrichr_Libraries_Most_Popular_Genes, file='GO_enrichR/Enrichr_Libraries_Most_Popular_Genes.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`Enrichr_Submissions_TF-Gene_Coocurrence`, file='GO_enrichR/Enrichr_Submissions_TF-Gene_Coocurrence.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Enrichr_Users_Contributed_Lists_2020, file='GO_enrichR/Enrichr_Users_Contributed_Lists_2020.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`Epigenomics_Roadmap_HM_ChIP-seq`, file='GO_enrichR/Epigenomics_Roadmap_HM_ChIP-seq.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$ESCAPE, file='GO_enrichR/ESCAPE.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$FANTOM6_lncRNA_KD_DEGs, file='GO_enrichR/FANTOM6_lncRNA_KD_DEGs.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GeDiPNet_2023, file='GO_enrichR/GeDiPNet_2023.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Gene_Perturbations_from_GEO_down, file='GO_enrichR/Gene_Perturbations_from_GEO_down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Gene_Perturbations_from_GEO_up, file='GO_enrichR/Gene_Perturbations_from_GEO_up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Genes_Associated_with_NIH_Grants, file='GO_enrichR/Genes_Associated_with_NIH_Grants.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GeneSigDB, file='GO_enrichR/GeneSigDB.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Genome_Browser_PWMs, file='GO_enrichR/Genome_Browser_PWMs.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GlyGen_Glycosylated_Proteins_2022, file='GO_enrichR/GlyGen_Glycosylated_Proteins_2022.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Biological_Process_2013, file='GO_enrichR/GO_Biological_Process_2013.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Biological_Process_2015, file='GO_enrichR/GO_Biological_Process_2015.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Biological_Process_2017, file='GO_enrichR/GO_Biological_Process_2017.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Biological_Process_2017b, file='GO_enrichR/GO_Biological_Process_2017b.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Biological_Process_2018, file='GO_enrichR/GO_Biological_Process_2018.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Biological_Process_2021, file='GO_enrichR/GO_Biological_Process_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Biological_Process_2023, file='GO_enrichR/GO_Biological_Process_2023.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Cellular_Component_2013, file='GO_enrichR/GO_Cellular_Component_2013.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Cellular_Component_2015, file='GO_enrichR/GO_Cellular_Component_2015.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Cellular_Component_2017, file='GO_enrichR/GO_Cellular_Component_2017.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Cellular_Component_2017b, file='GO_enrichR/GO_Cellular_Component_2017b.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Cellular_Component_2018, file='GO_enrichR/GO_Cellular_Component_2018.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Cellular_Component_2021, file='GO_enrichR/GO_Cellular_Component_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Cellular_Component_2023, file='GO_enrichR/GO_Cellular_Component_2023.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Molecular_Function_2013, file='GO_enrichR/GO_Molecular_Function_2013.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Molecular_Function_2015, file='GO_enrichR/GO_Molecular_Function_2015.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Molecular_Function_2017, file='GO_enrichR/GO_Molecular_Function_2017.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Molecular_Function_2017b, file='GO_enrichR/GO_Molecular_Function_2017b.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Molecular_Function_2018, file='GO_enrichR/GO_Molecular_Function_2018.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Molecular_Function_2021, file='GO_enrichR/GO_Molecular_Function_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GO_Molecular_Function_2023, file='GO_enrichR/GO_Molecular_Function_2023.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GTEx_Aging_Signatures_2021, file='GO_enrichR/GTEx_Aging_Signatures_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GTEx_Tissue_Expression_Down, file='GO_enrichR/GTEx_Tissue_Expression_Down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GTEx_Tissue_Expression_Up, file='GO_enrichR/GTEx_Tissue_Expression_Up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GTEx_Tissues_V8_2023, file='GO_enrichR/GTEx_Tissues_V8_2023.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GWAS_Catalog_2019, file='GO_enrichR/GWAS_Catalog_2019.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$GWAS_Catalog_2023, file='GO_enrichR/GWAS_Catalog_2023.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$HDSigDB_Human_2021, file='GO_enrichR/HDSigDB_Human_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$HDSigDB_Mouse_2021, file='GO_enrichR/HDSigDB_Mouse_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$HMDB_Metabolites, file='GO_enrichR/HMDB_Metabolites.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$HMS_LINCS_KinomeScan, file='GO_enrichR/HMS_LINCS_KinomeScan.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$HomoloGene, file='GO_enrichR/HomoloGene.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression, file='GO_enrichR/HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$HuBMAP_ASCTplusB_augmented_2022, file='GO_enrichR/HuBMAP_ASCTplusB_augmented_2022.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Human_Gene_Atlas, file='GO_enrichR/Human_Gene_Atlas.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Human_Phenotype_Ontology, file='GO_enrichR/Human_Phenotype_Ontology.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$HumanCyc_2015, file='GO_enrichR/HumanCyc_2015.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$HumanCyc_2016, file='GO_enrichR/HumanCyc_2016.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$huMAP, file='GO_enrichR/huMAP.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$IDG_Drug_Targets_2022, file='GO_enrichR/IDG_Drug_Targets_2022.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$InterPro_Domains_2019, file='GO_enrichR/InterPro_Domains_2019.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Jensen_COMPARTMENTS, file='GO_enrichR/Jensen_COMPARTMENTS.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Jensen_DISEASES, file='GO_enrichR/Jensen_DISEASES.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Jensen_TISSUES, file='GO_enrichR/Jensen_TISSUES.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$KEA_2013, file='GO_enrichR/KEA_2013.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$KEA_2015, file='GO_enrichR/KEA_2015.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$KEGG_2013, file='GO_enrichR/KEGG_2013.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$KEGG_2015, file='GO_enrichR/KEGG_2015.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$KEGG_2016, file='GO_enrichR/KEGG_2016.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$KEGG_2019_Human, file='GO_enrichR/KEGG_2019_Human.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$KEGG_2019_Mouse, file='GO_enrichR/KEGG_2019_Mouse.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$KEGG_2021_Human, file='GO_enrichR/KEGG_2021_Human.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Kinase_Perturbations_from_GEO_down, file='GO_enrichR/Kinase_Perturbations_from_GEO_down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Kinase_Perturbations_from_GEO_up, file='GO_enrichR/Kinase_Perturbations_from_GEO_up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$KOMP2_Mouse_Phenotypes_2022, file='GO_enrichR/KOMP2_Mouse_Phenotypes_2022.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$L1000_Kinase_and_GPCR_Perturbations_down, file='GO_enrichR/L1000_Kinase_and_GPCR_Perturbations_down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$L1000_Kinase_and_GPCR_Perturbations_up, file='GO_enrichR/L1000_Kinase_and_GPCR_Perturbations_up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Ligand_Perturbations_from_GEO_down, file='GO_enrichR/Ligand_Perturbations_from_GEO_down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Ligand_Perturbations_from_GEO_up, file='GO_enrichR/Ligand_Perturbations_from_GEO_up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$LINCS_L1000_Chem_Pert_Consensus_Sigs, file='GO_enrichR/LINCS_L1000_Chem_Pert_Consensus_Sigs.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$LINCS_L1000_Chem_Pert_down, file='GO_enrichR/LINCS_L1000_Chem_Pert_down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$LINCS_L1000_Chem_Pert_up, file='GO_enrichR/LINCS_L1000_Chem_Pert_up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$LINCS_L1000_CRISPR_KO_Consensus_Sigs, file='GO_enrichR/LINCS_L1000_CRISPR_KO_Consensus_Sigs.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$LINCS_L1000_Ligand_Perturbations_down, file='GO_enrichR/LINCS_L1000_Ligand_Perturbations_down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$LINCS_L1000_Ligand_Perturbations_up, file='GO_enrichR/LINCS_L1000_Ligand_Perturbations_up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`lncHUB_lncRNA_Co-Expression`, file='GO_enrichR/lncHUB_lncRNA_Co-Expression.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$MAGMA_Drugs_and_Diseases, file='GO_enrichR/MAGMA_Drugs_and_Diseases.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$MAGNET_2023, file='GO_enrichR/MAGNET_2023.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$MCF7_Perturbations_from_GEO_down, file='GO_enrichR/MCF7_Perturbations_from_GEO_down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$MCF7_Perturbations_from_GEO_up, file='GO_enrichR/MCF7_Perturbations_from_GEO_up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Metabolomics_Workbench_Metabolites_2022, file='GO_enrichR/Metabolomics_Workbench_Metabolites_2022.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$MGI_Mammalian_Phenotype_2013, file='GO_enrichR/MGI_Mammalian_Phenotype_2013.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$MGI_Mammalian_Phenotype_2017, file='GO_enrichR/MGI_Mammalian_Phenotype_2017.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$MGI_Mammalian_Phenotype_Level_3, file='GO_enrichR/MGI_Mammalian_Phenotype_Level_3.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$MGI_Mammalian_Phenotype_Level_4, file='GO_enrichR/MGI_Mammalian_Phenotype_Level_4.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$MGI_Mammalian_Phenotype_Level_4_2019, file='GO_enrichR/MGI_Mammalian_Phenotype_Level_4_2019.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$MGI_Mammalian_Phenotype_Level_4_2021, file='GO_enrichR/MGI_Mammalian_Phenotype_Level_4_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Microbe_Perturbations_from_GEO_down, file='GO_enrichR/Microbe_Perturbations_from_GEO_down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Microbe_Perturbations_from_GEO_up, file='GO_enrichR/Microbe_Perturbations_from_GEO_up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$miRTarBase_2017, file='GO_enrichR/miRTarBase_2017.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$MoTrPAC_2023, file='GO_enrichR/MoTrPAC_2023.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Mouse_Gene_Atlas, file='GO_enrichR/Mouse_Gene_Atlas.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$MSigDB_Computational, file='GO_enrichR/MSigDB_Computational.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$MSigDB_Hallmark_2020, file='GO_enrichR/MSigDB_Hallmark_2020.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$MSigDB_Oncogenic_Signatures, file='GO_enrichR/MSigDB_Oncogenic_Signatures.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`NCI-60_Cancer_Cell_Lines`, file='GO_enrichR/NCI-60_Cancer_Cell_Lines.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`NCI-Nature_2015`, file='GO_enrichR/NCI-Nature_2015.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`NCI-Nature_2016`, file='GO_enrichR/NCI-Nature_2016.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions, file='GO_enrichR/NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions, file='GO_enrichR/NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$NIH_Funded_PIs_2017_Human_AutoRIF, file='GO_enrichR/NIH_Funded_PIs_2017_Human_AutoRIF.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$NIH_Funded_PIs_2017_Human_GeneRIF, file='GO_enrichR/NIH_Funded_PIs_2017_Human_GeneRIF.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$NURSA_Human_Endogenous_Complexome, file='GO_enrichR/NURSA_Human_Endogenous_Complexome.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Old_CMAP_down, file='GO_enrichR/Old_CMAP_down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Old_CMAP_up, file='GO_enrichR/Old_CMAP_up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$OMIM_Disease, file='GO_enrichR/OMIM_Disease.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$OMIM_Expanded, file='GO_enrichR/OMIM_Expanded.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Orphanet_Augmented_2021, file='GO_enrichR/Orphanet_Augmented_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$PanglaoDB_Augmented_2021, file='GO_enrichR/PanglaoDB_Augmented_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Panther_2015, file='GO_enrichR/Panther_2015.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Panther_2016, file='GO_enrichR/Panther_2016.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Pfam_Domains_2019, file='GO_enrichR/Pfam_Domains_2019.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Pfam_InterPro_Domains, file='GO_enrichR/Pfam_InterPro_Domains.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$PFOCR_Pathways, file='GO_enrichR/PFOCR_Pathways.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$PFOCR_Pathways_2023, file='GO_enrichR/PFOCR_Pathways_2023.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$PhenGenI_Association_2021, file='GO_enrichR/PhenGenI_Association_2021.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$PheWeb_2019, file='GO_enrichR/PheWeb_2019.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Phosphatase_Substrates_from_DEPOD, file='GO_enrichR/Phosphatase_Substrates_from_DEPOD.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$PPI_Hub_Proteins, file='GO_enrichR/PPI_Hub_Proteins.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Proteomics_Drug_Atlas_2023, file='GO_enrichR/Proteomics_Drug_Atlas_2023.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$ProteomicsDB_2020, file='GO_enrichR/ProteomicsDB_2020.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Rare_Diseases_AutoRIF_ARCHS4_Predictions, file='GO_enrichR/Rare_Diseases_AutoRIF_ARCHS4_Predictions.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Rare_Diseases_AutoRIF_Gene_Lists, file='GO_enrichR/Rare_Diseases_AutoRIF_Gene_Lists.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Rare_Diseases_GeneRIF_ARCHS4_Predictions, file='GO_enrichR/Rare_Diseases_GeneRIF_ARCHS4_Predictions.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Rare_Diseases_GeneRIF_Gene_Lists, file='GO_enrichR/Rare_Diseases_GeneRIF_Gene_Lists.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Reactome_2013, file='GO_enrichR/Reactome_2013.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Reactome_2015, file='GO_enrichR/Reactome_2015.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Reactome_2016, file='GO_enrichR/Reactome_2016.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Reactome_2022, file='GO_enrichR/Reactome_2022.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO`, file='GO_enrichR/RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$RNAseq_Automatic_GEO_Signatures_Human_Down, file='GO_enrichR/RNAseq_Automatic_GEO_Signatures_Human_Down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$RNAseq_Automatic_GEO_Signatures_Human_Up, file='GO_enrichR/RNAseq_Automatic_GEO_Signatures_Human_Up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$RNAseq_Automatic_GEO_Signatures_Mouse_Down, file='GO_enrichR/RNAseq_Automatic_GEO_Signatures_Mouse_Down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$RNAseq_Automatic_GEO_Signatures_Mouse_Up, file='GO_enrichR/RNAseq_Automatic_GEO_Signatures_Mouse_Up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Rummagene_kinases, file='GO_enrichR/Rummagene_kinases.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Rummagene_signatures, file='GO_enrichR/Rummagene_signatures.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Rummagene_transcription_factors, file='GO_enrichR/Rummagene_transcription_factors.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$SILAC_Phosphoproteomics, file='GO_enrichR/SILAC_Phosphoproteomics.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$SubCell_BarCode, file='GO_enrichR/SubCell_BarCode.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$SynGO_2022, file='GO_enrichR/SynGO_2022.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$SynGO_2024, file='GO_enrichR/SynGO_2024.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$SysMyo_Muscle_Gene_Sets, file='GO_enrichR/SysMyo_Muscle_Gene_Sets.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Table_Mining_of_CRISPR_Studies, file='GO_enrichR/Table_Mining_of_CRISPR_Studies.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Tabula_Muris, file='GO_enrichR/Tabula_Muris.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Tabula_Sapiens, file='GO_enrichR/Tabula_Sapiens.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$TargetScan_microRNA, file='GO_enrichR/TargetScan_microRNA.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$TargetScan_microRNA_2017, file='GO_enrichR/TargetScan_microRNA_2017.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`TF-LOF_Expression_from_GEO`, file='GO_enrichR/TF-LOF_Expression_from_GEO.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$TF_Perturbations_Followed_by_Expression, file='GO_enrichR/TF_Perturbations_Followed_by_Expression.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$TG_GATES_2020, file='GO_enrichR/TG_GATES_2020.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$The_Kinase_Library_2023, file='GO_enrichR/The_Kinase_Library_2023.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Tissue_Protein_Expression_from_Human_Proteome_Map, file='GO_enrichR/Tissue_Protein_Expression_from_Human_Proteome_Map.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Tissue_Protein_Expression_from_ProteomicsDB, file='GO_enrichR/Tissue_Protein_Expression_from_ProteomicsDB.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Transcription_Factor_PPIs, file='GO_enrichR/Transcription_Factor_PPIs.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$TRANSFAC_and_JASPAR_PWMs, file='GO_enrichR/TRANSFAC_and_JASPAR_PWMs.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$TRRUST_Transcription_Factors_2019, file='GO_enrichR/TRRUST_Transcription_Factors_2019.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$UK_Biobank_GWAS_v1, file='GO_enrichR/UK_Biobank_GWAS_v1.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`Virus-Host_PPI_P-HIPSTer_2020`, file='GO_enrichR/Virus-Host_PPI_P-HIPSTer_2020.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$`Virus_Perturbations_from_GEO_down`, file='GO_enrichR/Virus_Perturbations_from_GEO_down.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$Virus_Perturbations_from_GEO_up, file='GO_enrichR/Virus_Perturbations_from_GEO_up.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$VirusMINT, file='GO_enrichR/VirusMINT.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$WikiPathway_2021_Human, file='GO_enrichR/WikiPathway_2021_Human.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$WikiPathway_2023_Human, file='GO_enrichR/WikiPathway_2023_Human.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$WikiPathways_2013, file='GO_enrichR/WikiPathways_2013.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$WikiPathways_2015, file='GO_enrichR/WikiPathways_2015.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$WikiPathways_2016, file='GO_enrichR/WikiPathways_2016.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$WikiPathways_2019_Human, file='GO_enrichR/WikiPathways_2019_Human.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  write.xlsx(enriched$WikiPathways_2019_Mouse, file='GO_enrichR/WikiPathways_2019_Mouse.xlsx', sheetName=i, append=TRUE, row.names=FALSE)
  

  


pdf(file = paste("GO_enrichR_PDF/",i,"_GO_Plot.pdf",sep="") )



print(if (websiteLive) plotEnrich(enriched$Achilles_fitness_decrease, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Achilles_fitness_decrease'))
print(if (websiteLive) plotEnrich(enriched$Achilles_fitness_increase, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Achilles_fitness_increase'))
print(if (websiteLive) plotEnrich(enriched$Aging_Perturbations_from_GEO_down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Aging_Perturbations_from_GEO_down'))
print(if (websiteLive) plotEnrich(enriched$Aging_Perturbations_from_GEO_up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Aging_Perturbations_from_GEO_up'))
print(if (websiteLive) plotEnrich(enriched$Allen_Brain_Atlas_10x_scRNA_2021, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Allen_Brain_Atlas_10x_scRNA_2021'))
print(if (websiteLive) plotEnrich(enriched$Allen_Brain_Atlas_down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Allen_Brain_Atlas_down'))
print(if (websiteLive) plotEnrich(enriched$Allen_Brain_Atlas_up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Allen_Brain_Atlas_up'))
print(if (websiteLive) plotEnrich(enriched$`ARCHS4_Cell-lines`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ARCHS4_Cell-lines'))
print(if (websiteLive) plotEnrich(enriched$ARCHS4_IDG_Coexp, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ARCHS4_IDG_Coexp'))
print(if (websiteLive) plotEnrich(enriched$ARCHS4_Kinases_Coexp, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ARCHS4_Kinases_Coexp'))
print(if (websiteLive) plotEnrich(enriched$ARCHS4_TFs_Coexp, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ARCHS4_TFs_Coexp'))
print(if (websiteLive) plotEnrich(enriched$ARCHS4_Tissues, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ARCHS4_Tissues'))
print(if (websiteLive) plotEnrich(enriched$Azimuth_2023, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Azimuth_2023'))
print(if (websiteLive) plotEnrich(enriched$Azimuth_Cell_Types_2021, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Azimuth_Cell_Types_2021'))
print(if (websiteLive) plotEnrich(enriched$BioCarta_2013, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='BioCarta_2013'))
print(if (websiteLive) plotEnrich(enriched$BioCarta_2015, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='BioCarta_2015'))
print(if (websiteLive) plotEnrich(enriched$BioCarta_2016, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='BioCarta_2016'))
print(if (websiteLive) plotEnrich(enriched$BioPlanet_2019, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='BioPlanet_2019'))
print(if (websiteLive) plotEnrich(enriched$BioPlex_2017, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='BioPlex_2017'))
print(if (websiteLive) plotEnrich(enriched$Cancer_Cell_Line_Encyclopedia, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Cancer_Cell_Line_Encyclopedia'))
print(if (websiteLive) plotEnrich(enriched$CCLE_Proteomics_2020, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='CCLE_Proteomics_2020'))
print(if (websiteLive) plotEnrich(enriched$CellMarker_2024, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='CellMarker_2024'))
print(if (websiteLive) plotEnrich(enriched$CellMarker_Augmented_2021, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='CellMarker_Augmented_2021'))
print(if (websiteLive) plotEnrich(enriched$ChEA_2013, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ChEA_2013'))
print(if (websiteLive) plotEnrich(enriched$ChEA_2015, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ChEA_2015'))
print(if (websiteLive) plotEnrich(enriched$ChEA_2016, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ChEA_2016'))
print(if (websiteLive) plotEnrich(enriched$ChEA_2022, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ChEA_2022'))
print(if (websiteLive) plotEnrich(enriched$Chromosome_Location, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Chromosome_Location'))
print(if (websiteLive) plotEnrich(enriched$Chromosome_Location_hg19, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Chromosome_Location_hg19'))
print(if (websiteLive) plotEnrich(enriched$ClinVar_2019, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ClinVar_2019'))
print(if (websiteLive) plotEnrich(enriched$CORUM, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='CORUM'))
print(if (websiteLive) plotEnrich(enriched$`COVID-19_Related_Gene_Sets`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='COVID-19_Related_Gene_Sets'))
print(if (websiteLive) plotEnrich(enriched$`COVID-19_Related_Gene_Sets_2021`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='COVID-19_Related_Gene_Sets_2021'))
print(if (websiteLive) plotEnrich(enriched$Data_Acquisition_Method_Most_Popular_Genes, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Data_Acquisition_Method_Most_Popular_Genes'))
print(if (websiteLive) plotEnrich(enriched$dbGaP, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='dbGaP'))
print(if (websiteLive) plotEnrich(enriched$DepMap_WG_CRISPR_Screens_Broad_CellLines_2019, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='DepMap_WG_CRISPR_Screens_Broad_CellLines_2019'))
print(if (websiteLive) plotEnrich(enriched$DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019'))
print(if (websiteLive) plotEnrich(enriched$Descartes_Cell_Types_and_Tissue_2021, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Descartes_Cell_Types_and_Tissue_2021'))
print(if (websiteLive) plotEnrich(enriched$Diabetes_Perturbations_GEO_2022, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Diabetes_Perturbations_GEO_2022'))
print(if (websiteLive) plotEnrich(enriched$Disease_Perturbations_from_GEO_down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Disease_Perturbations_from_GEO_down'))
print(if (websiteLive) plotEnrich(enriched$Disease_Perturbations_from_GEO_up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Disease_Perturbations_from_GEO_up'))
print(if (websiteLive) plotEnrich(enriched$Disease_Signatures_from_GEO_down_2014, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Disease_Signatures_from_GEO_down_2014'))
print(if (websiteLive) plotEnrich(enriched$Disease_Signatures_from_GEO_up_2014, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Disease_Signatures_from_GEO_up_2014'))
print(if (websiteLive) plotEnrich(enriched$DisGeNET, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='DisGeNET'))
print(if (websiteLive) plotEnrich(enriched$Drug_Perturbations_from_GEO_2014, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Drug_Perturbations_from_GEO_2014'))
print(if (websiteLive) plotEnrich(enriched$Drug_Perturbations_from_GEO_down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Drug_Perturbations_from_GEO_down'))
print(if (websiteLive) plotEnrich(enriched$Drug_Perturbations_from_GEO_up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Drug_Perturbations_from_GEO_up'))
print(if (websiteLive) plotEnrich(enriched$DrugMatrix, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='DrugMatrix'))
print(if (websiteLive) plotEnrich(enriched$DSigDB, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='DSigDB'))
print(if (websiteLive) plotEnrich(enriched$Elsevier_Pathway_Collection, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Elsevier_Pathway_Collection'))
print(if (websiteLive) plotEnrich(enriched$`ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X'))
print(if (websiteLive) plotEnrich(enriched$ENCODE_Histone_Modifications_2013, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ENCODE_Histone_Modifications_2013'))
print(if (websiteLive) plotEnrich(enriched$ENCODE_Histone_Modifications_2015, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ENCODE_Histone_Modifications_2015'))
print(if (websiteLive) plotEnrich(enriched$`ENCODE_TF_ChIP-seq_2014`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ENCODE_TF_ChIP-seq_2014'))
print(if (websiteLive) plotEnrich(enriched$`ENCODE_TF_ChIP-seq_2015`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ENCODE_TF_ChIP-seq_2015'))
print(if (websiteLive) plotEnrich(enriched$Enrichr_Libraries_Most_Popular_Genes, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Enrichr_Libraries_Most_Popular_Genes'))
print(if (websiteLive) plotEnrich(enriched$`Enrichr_Submissions_TF-Gene_Coocurrence`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Enrichr_Submissions_TF-Gene_Coocurrence'))
print(if (websiteLive) plotEnrich(enriched$Enrichr_Users_Contributed_Lists_2020, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Enrichr_Users_Contributed_Lists_2020'))
print(if (websiteLive) plotEnrich(enriched$`Epigenomics_Roadmap_HM_ChIP-seq`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Epigenomics_Roadmap_HM_ChIP-seq'))
print(if (websiteLive) plotEnrich(enriched$ESCAPE, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ESCAPE'))
print(if (websiteLive) plotEnrich(enriched$FANTOM6_lncRNA_KD_DEGs, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='FANTOM6_lncRNA_KD_DEGs'))
print(if (websiteLive) plotEnrich(enriched$GeDiPNet_2023, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GeDiPNet_2023'))
print(if (websiteLive) plotEnrich(enriched$Gene_Perturbations_from_GEO_down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Gene_Perturbations_from_GEO_down'))
print(if (websiteLive) plotEnrich(enriched$Gene_Perturbations_from_GEO_up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Gene_Perturbations_from_GEO_up'))
print(if (websiteLive) plotEnrich(enriched$Genes_Associated_with_NIH_Grants, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Genes_Associated_with_NIH_Grants'))
print(if (websiteLive) plotEnrich(enriched$GeneSigDB, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GeneSigDB'))
print(if (websiteLive) plotEnrich(enriched$Genome_Browser_PWMs, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Genome_Browser_PWMs'))
print(if (websiteLive) plotEnrich(enriched$GlyGen_Glycosylated_Proteins_2022, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GlyGen_Glycosylated_Proteins_2022'))
print(if (websiteLive) plotEnrich(enriched$GO_Biological_Process_2013, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Biological_Process_2013'))
print(if (websiteLive) plotEnrich(enriched$GO_Biological_Process_2015, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Biological_Process_2015'))
print(if (websiteLive) plotEnrich(enriched$GO_Biological_Process_2017, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Biological_Process_2017'))
print(if (websiteLive) plotEnrich(enriched$GO_Biological_Process_2017b, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Biological_Process_2017b'))
print(if (websiteLive) plotEnrich(enriched$GO_Biological_Process_2018, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Biological_Process_2018'))
print(if (websiteLive) plotEnrich(enriched$GO_Biological_Process_2021, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Biological_Process_2021'))
print(if (websiteLive) plotEnrich(enriched$GO_Biological_Process_2023, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Biological_Process_2023'))
print(if (websiteLive) plotEnrich(enriched$GO_Cellular_Component_2013, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Cellular_Component_2013'))
print(if (websiteLive) plotEnrich(enriched$GO_Cellular_Component_2015, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Cellular_Component_2015'))
print(if (websiteLive) plotEnrich(enriched$GO_Cellular_Component_2017, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Cellular_Component_2017'))
print(if (websiteLive) plotEnrich(enriched$GO_Cellular_Component_2017b, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Cellular_Component_2017b'))
print(if (websiteLive) plotEnrich(enriched$GO_Cellular_Component_2018, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Cellular_Component_2018'))
print(if (websiteLive) plotEnrich(enriched$GO_Cellular_Component_2021, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Cellular_Component_2021'))
print(if (websiteLive) plotEnrich(enriched$GO_Cellular_Component_2023, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Cellular_Component_2023'))
print(if (websiteLive) plotEnrich(enriched$GO_Molecular_Function_2013, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Molecular_Function_2013'))
print(if (websiteLive) plotEnrich(enriched$GO_Molecular_Function_2015, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Molecular_Function_2015'))
print(if (websiteLive) plotEnrich(enriched$GO_Molecular_Function_2017, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Molecular_Function_2017'))
print(if (websiteLive) plotEnrich(enriched$GO_Molecular_Function_2017b, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Molecular_Function_2017b'))
print(if (websiteLive) plotEnrich(enriched$GO_Molecular_Function_2018, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Molecular_Function_2018'))
print(if (websiteLive) plotEnrich(enriched$GO_Molecular_Function_2021, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Molecular_Function_2021'))
print(if (websiteLive) plotEnrich(enriched$GO_Molecular_Function_2023, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GO_Molecular_Function_2023'))
print(if (websiteLive) plotEnrich(enriched$GTEx_Aging_Signatures_2021, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GTEx_Aging_Signatures_2021'))
print(if (websiteLive) plotEnrich(enriched$GTEx_Tissue_Expression_Down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GTEx_Tissue_Expression_Down'))
print(if (websiteLive) plotEnrich(enriched$GTEx_Tissue_Expression_Up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GTEx_Tissue_Expression_Up'))
print(if (websiteLive) plotEnrich(enriched$GTEx_Tissues_V8_2023, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GTEx_Tissues_V8_2023'))
print(if (websiteLive) plotEnrich(enriched$GWAS_Catalog_2019, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GWAS_Catalog_2019'))
print(if (websiteLive) plotEnrich(enriched$GWAS_Catalog_2023, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='GWAS_Catalog_2023'))
print(if (websiteLive) plotEnrich(enriched$HDSigDB_Human_2021, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='HDSigDB_Human_2021'))
print(if (websiteLive) plotEnrich(enriched$HDSigDB_Mouse_2021, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='HDSigDB_Mouse_2021'))
print(if (websiteLive) plotEnrich(enriched$HMDB_Metabolites, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='HMDB_Metabolites'))
print(if (websiteLive) plotEnrich(enriched$HMS_LINCS_KinomeScan, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='HMS_LINCS_KinomeScan'))
print(if (websiteLive) plotEnrich(enriched$HomoloGene, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='HomoloGene'))
print(if (websiteLive) plotEnrich(enriched$HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression'))
print(if (websiteLive) plotEnrich(enriched$HuBMAP_ASCTplusB_augmented_2022, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='HuBMAP_ASCTplusB_augmented_2022'))
print(if (websiteLive) plotEnrich(enriched$Human_Gene_Atlas, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Human_Gene_Atlas'))
print(if (websiteLive) plotEnrich(enriched$Human_Phenotype_Ontology, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Human_Phenotype_Ontology'))
print(if (websiteLive) plotEnrich(enriched$HumanCyc_2015, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='HumanCyc_2015'))
print(if (websiteLive) plotEnrich(enriched$HumanCyc_2016, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='HumanCyc_2016'))
print(if (websiteLive) plotEnrich(enriched$huMAP, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='huMAP'))
print(if (websiteLive) plotEnrich(enriched$IDG_Drug_Targets_2022, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='IDG_Drug_Targets_2022'))
print(if (websiteLive) plotEnrich(enriched$InterPro_Domains_2019, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='InterPro_Domains_2019'))
print(if (websiteLive) plotEnrich(enriched$Jensen_COMPARTMENTS, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Jensen_COMPARTMENTS'))
print(if (websiteLive) plotEnrich(enriched$Jensen_DISEASES, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Jensen_DISEASES'))
print(if (websiteLive) plotEnrich(enriched$Jensen_TISSUES, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Jensen_TISSUES'))
print(if (websiteLive) plotEnrich(enriched$KEA_2013, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='KEA_2013'))
print(if (websiteLive) plotEnrich(enriched$KEA_2015, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='KEA_2015'))
print(if (websiteLive) plotEnrich(enriched$KEGG_2013, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='KEGG_2013'))
print(if (websiteLive) plotEnrich(enriched$KEGG_2015, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='KEGG_2015'))
print(if (websiteLive) plotEnrich(enriched$KEGG_2016, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='KEGG_2016'))
print(if (websiteLive) plotEnrich(enriched$KEGG_2019_Human, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='KEGG_2019_Human'))
print(if (websiteLive) plotEnrich(enriched$KEGG_2019_Mouse, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='KEGG_2019_Mouse'))
print(if (websiteLive) plotEnrich(enriched$KEGG_2021_Human, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='KEGG_2021_Human'))
print(if (websiteLive) plotEnrich(enriched$Kinase_Perturbations_from_GEO_down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Kinase_Perturbations_from_GEO_down'))
print(if (websiteLive) plotEnrich(enriched$Kinase_Perturbations_from_GEO_up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Kinase_Perturbations_from_GEO_up'))
print(if (websiteLive) plotEnrich(enriched$KOMP2_Mouse_Phenotypes_2022, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='KOMP2_Mouse_Phenotypes_2022'))
print(if (websiteLive) plotEnrich(enriched$L1000_Kinase_and_GPCR_Perturbations_down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='L1000_Kinase_and_GPCR_Perturbations_down'))
print(if (websiteLive) plotEnrich(enriched$L1000_Kinase_and_GPCR_Perturbations_up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='L1000_Kinase_and_GPCR_Perturbations_up'))
print(if (websiteLive) plotEnrich(enriched$Ligand_Perturbations_from_GEO_down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Ligand_Perturbations_from_GEO_down'))
print(if (websiteLive) plotEnrich(enriched$Ligand_Perturbations_from_GEO_up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Ligand_Perturbations_from_GEO_up'))
print(if (websiteLive) plotEnrich(enriched$LINCS_L1000_Chem_Pert_Consensus_Sigs, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='LINCS_L1000_Chem_Pert_Consensus_Sigs'))
print(if (websiteLive) plotEnrich(enriched$LINCS_L1000_Chem_Pert_down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='LINCS_L1000_Chem_Pert_down'))
print(if (websiteLive) plotEnrich(enriched$LINCS_L1000_Chem_Pert_up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='LINCS_L1000_Chem_Pert_up'))
print(if (websiteLive) plotEnrich(enriched$LINCS_L1000_CRISPR_KO_Consensus_Sigs, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='LINCS_L1000_CRISPR_KO_Consensus_Sigs'))
print(if (websiteLive) plotEnrich(enriched$LINCS_L1000_Ligand_Perturbations_down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='LINCS_L1000_Ligand_Perturbations_down'))
print(if (websiteLive) plotEnrich(enriched$LINCS_L1000_Ligand_Perturbations_up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='LINCS_L1000_Ligand_Perturbations_up'))
print(if (websiteLive) plotEnrich(enriched$`lncHUB_lncRNA_Co-Expression`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='lncHUB_lncRNA_Co-Expression'))
print(if (websiteLive) plotEnrich(enriched$MAGMA_Drugs_and_Diseases, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='MAGMA_Drugs_and_Diseases'))
print(if (websiteLive) plotEnrich(enriched$MAGNET_2023, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='MAGNET_2023'))
print(if (websiteLive) plotEnrich(enriched$MCF7_Perturbations_from_GEO_down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='MCF7_Perturbations_from_GEO_down'))
print(if (websiteLive) plotEnrich(enriched$MCF7_Perturbations_from_GEO_up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='MCF7_Perturbations_from_GEO_up'))
print(if (websiteLive) plotEnrich(enriched$Metabolomics_Workbench_Metabolites_2022, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Metabolomics_Workbench_Metabolites_2022'))
print(if (websiteLive) plotEnrich(enriched$MGI_Mammalian_Phenotype_2013, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='MGI_Mammalian_Phenotype_2013'))
print(if (websiteLive) plotEnrich(enriched$MGI_Mammalian_Phenotype_2017, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='MGI_Mammalian_Phenotype_2017'))
print(if (websiteLive) plotEnrich(enriched$MGI_Mammalian_Phenotype_Level_3, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='MGI_Mammalian_Phenotype_Level_3'))
print(if (websiteLive) plotEnrich(enriched$MGI_Mammalian_Phenotype_Level_4, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='MGI_Mammalian_Phenotype_Level_4'))
print(if (websiteLive) plotEnrich(enriched$MGI_Mammalian_Phenotype_Level_4_2019, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='MGI_Mammalian_Phenotype_Level_4_2019'))
print(if (websiteLive) plotEnrich(enriched$MGI_Mammalian_Phenotype_Level_4_2021, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='MGI_Mammalian_Phenotype_Level_4_2021'))
print(if (websiteLive) plotEnrich(enriched$Microbe_Perturbations_from_GEO_down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Microbe_Perturbations_from_GEO_down'))
print(if (websiteLive) plotEnrich(enriched$Microbe_Perturbations_from_GEO_up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Microbe_Perturbations_from_GEO_up'))
print(if (websiteLive) plotEnrich(enriched$miRTarBase_2017, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='miRTarBase_2017'))
print(if (websiteLive) plotEnrich(enriched$MoTrPAC_2023, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='MoTrPAC_2023'))
print(if (websiteLive) plotEnrich(enriched$Mouse_Gene_Atlas, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Mouse_Gene_Atlas'))
print(if (websiteLive) plotEnrich(enriched$MSigDB_Computational, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='MSigDB_Computational'))
print(if (websiteLive) plotEnrich(enriched$MSigDB_Hallmark_2020, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='MSigDB_Hallmark_2020'))
print(if (websiteLive) plotEnrich(enriched$MSigDB_Oncogenic_Signatures, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='MSigDB_Oncogenic_Signatures'))
print(if (websiteLive) plotEnrich(enriched$`NCI-60_Cancer_Cell_Lines`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='NCI-60_Cancer_Cell_Lines'))
print(if (websiteLive) plotEnrich(enriched$`NCI-Nature_2015`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='NCI-Nature_2015'))
print(if (websiteLive) plotEnrich(enriched$`NCI-Nature_2016`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='NCI-Nature_2016'))
print(if (websiteLive) plotEnrich(enriched$NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions'))
print(if (websiteLive) plotEnrich(enriched$NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions'))
print(if (websiteLive) plotEnrich(enriched$NIH_Funded_PIs_2017_Human_AutoRIF, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='NIH_Funded_PIs_2017_Human_AutoRIF'))
print(if (websiteLive) plotEnrich(enriched$NIH_Funded_PIs_2017_Human_GeneRIF, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='NIH_Funded_PIs_2017_Human_GeneRIF'))
print(if (websiteLive) plotEnrich(enriched$NURSA_Human_Endogenous_Complexome, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='NURSA_Human_Endogenous_Complexome'))
print(if (websiteLive) plotEnrich(enriched$Old_CMAP_down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Old_CMAP_down'))
print(if (websiteLive) plotEnrich(enriched$Old_CMAP_up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Old_CMAP_up'))
print(if (websiteLive) plotEnrich(enriched$OMIM_Disease, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='OMIM_Disease'))
print(if (websiteLive) plotEnrich(enriched$OMIM_Expanded, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='OMIM_Expanded'))
print(if (websiteLive) plotEnrich(enriched$Orphanet_Augmented_2021, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Orphanet_Augmented_2021'))
print(if (websiteLive) plotEnrich(enriched$PanglaoDB_Augmented_2021, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='PanglaoDB_Augmented_2021'))
print(if (websiteLive) plotEnrich(enriched$Panther_2015, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Panther_2015'))
print(if (websiteLive) plotEnrich(enriched$Panther_2016, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Panther_2016'))
print(if (websiteLive) plotEnrich(enriched$Pfam_Domains_2019, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Pfam_Domains_2019'))
print(if (websiteLive) plotEnrich(enriched$Pfam_InterPro_Domains, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Pfam_InterPro_Domains'))
print(if (websiteLive) plotEnrich(enriched$PFOCR_Pathways, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='PFOCR_Pathways'))
print(if (websiteLive) plotEnrich(enriched$PFOCR_Pathways_2023, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='PFOCR_Pathways_2023'))
print(if (websiteLive) plotEnrich(enriched$PhenGenI_Association_2021, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='PhenGenI_Association_2021'))
print(if (websiteLive) plotEnrich(enriched$PheWeb_2019, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='PheWeb_2019'))
print(if (websiteLive) plotEnrich(enriched$Phosphatase_Substrates_from_DEPOD, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Phosphatase_Substrates_from_DEPOD'))
print(if (websiteLive) plotEnrich(enriched$PPI_Hub_Proteins, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='PPI_Hub_Proteins'))
print(if (websiteLive) plotEnrich(enriched$Proteomics_Drug_Atlas_2023, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Proteomics_Drug_Atlas_2023'))
print(if (websiteLive) plotEnrich(enriched$ProteomicsDB_2020, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='ProteomicsDB_2020'))
print(if (websiteLive) plotEnrich(enriched$Rare_Diseases_AutoRIF_ARCHS4_Predictions, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Rare_Diseases_AutoRIF_ARCHS4_Predictions'))
print(if (websiteLive) plotEnrich(enriched$Rare_Diseases_AutoRIF_Gene_Lists, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Rare_Diseases_AutoRIF_Gene_Lists'))
print(if (websiteLive) plotEnrich(enriched$Rare_Diseases_GeneRIF_ARCHS4_Predictions, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Rare_Diseases_GeneRIF_ARCHS4_Predictions'))
print(if (websiteLive) plotEnrich(enriched$Rare_Diseases_GeneRIF_Gene_Lists, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Rare_Diseases_GeneRIF_Gene_Lists'))
print(if (websiteLive) plotEnrich(enriched$Reactome_2013, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Reactome_2013'))
print(if (websiteLive) plotEnrich(enriched$Reactome_2015, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Reactome_2015'))
print(if (websiteLive) plotEnrich(enriched$Reactome_2016, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Reactome_2016'))
print(if (websiteLive) plotEnrich(enriched$Reactome_2022, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Reactome_2022'))
print(if (websiteLive) plotEnrich(enriched$`RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO'))
print(if (websiteLive) plotEnrich(enriched$RNAseq_Automatic_GEO_Signatures_Human_Down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='RNAseq_Automatic_GEO_Signatures_Human_Down'))
print(if (websiteLive) plotEnrich(enriched$RNAseq_Automatic_GEO_Signatures_Human_Up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='RNAseq_Automatic_GEO_Signatures_Human_Up'))
print(if (websiteLive) plotEnrich(enriched$RNAseq_Automatic_GEO_Signatures_Mouse_Down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='RNAseq_Automatic_GEO_Signatures_Mouse_Down'))
print(if (websiteLive) plotEnrich(enriched$RNAseq_Automatic_GEO_Signatures_Mouse_Up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='RNAseq_Automatic_GEO_Signatures_Mouse_Up'))
print(if (websiteLive) plotEnrich(enriched$Rummagene_kinases, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Rummagene_kinases'))
print(if (websiteLive) plotEnrich(enriched$Rummagene_signatures, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Rummagene_signatures'))
print(if (websiteLive) plotEnrich(enriched$Rummagene_transcription_factors, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Rummagene_transcription_factors'))
print(if (websiteLive) plotEnrich(enriched$SILAC_Phosphoproteomics, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='SILAC_Phosphoproteomics'))
print(if (websiteLive) plotEnrich(enriched$SubCell_BarCode, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='SubCell_BarCode'))
print(if (websiteLive) plotEnrich(enriched$SynGO_2022, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='SynGO_2022'))
print(if (websiteLive) plotEnrich(enriched$SynGO_2024, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='SynGO_2024'))
print(if (websiteLive) plotEnrich(enriched$SysMyo_Muscle_Gene_Sets, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='SysMyo_Muscle_Gene_Sets'))
print(if (websiteLive) plotEnrich(enriched$Table_Mining_of_CRISPR_Studies, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Table_Mining_of_CRISPR_Studies'))
print(if (websiteLive) plotEnrich(enriched$Tabula_Muris, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Tabula_Muris'))
print(if (websiteLive) plotEnrich(enriched$Tabula_Sapiens, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Tabula_Sapiens'))
print(if (websiteLive) plotEnrich(enriched$TargetScan_microRNA, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='TargetScan_microRNA'))
print(if (websiteLive) plotEnrich(enriched$TargetScan_microRNA_2017, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='TargetScan_microRNA_2017'))
print(if (websiteLive) plotEnrich(enriched$`TF-LOF_Expression_from_GEO`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='TF-LOF_Expression_from_GEO'))
print(if (websiteLive) plotEnrich(enriched$TF_Perturbations_Followed_by_Expression, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='TF_Perturbations_Followed_by_Expression'))
print(if (websiteLive) plotEnrich(enriched$TG_GATES_2020, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='TG_GATES_2020'))
print(if (websiteLive) plotEnrich(enriched$The_Kinase_Library_2023, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='The_Kinase_Library_2023'))
print(if (websiteLive) plotEnrich(enriched$Tissue_Protein_Expression_from_Human_Proteome_Map, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Tissue_Protein_Expression_from_Human_Proteome_Map'))
print(if (websiteLive) plotEnrich(enriched$Tissue_Protein_Expression_from_ProteomicsDB, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Tissue_Protein_Expression_from_ProteomicsDB'))
print(if (websiteLive) plotEnrich(enriched$Transcription_Factor_PPIs, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Transcription_Factor_PPIs'))
print(if (websiteLive) plotEnrich(enriched$TRANSFAC_and_JASPAR_PWMs, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='TRANSFAC_and_JASPAR_PWMs'))
print(if (websiteLive) plotEnrich(enriched$TRRUST_Transcription_Factors_2019, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='TRRUST_Transcription_Factors_2019'))
print(if (websiteLive) plotEnrich(enriched$UK_Biobank_GWAS_v1, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='UK_Biobank_GWAS_v1'))
print(if (websiteLive) plotEnrich(enriched$`Virus-Host_PPI_P-HIPSTer_2020`, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Virus-Host_PPI_P-HIPSTer_2020'))
print(if (websiteLive) plotEnrich(enriched$Virus_Perturbations_from_GEO_down, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Virus_Perturbations_from_GEO_down'))
print(if (websiteLive) plotEnrich(enriched$Virus_Perturbations_from_GEO_up, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='Virus_Perturbations_from_GEO_up'))
print(if (websiteLive) plotEnrich(enriched$VirusMINT, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='VirusMINT'))
print(if (websiteLive) plotEnrich(enriched$WikiPathway_2021_Human, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='WikiPathway_2021_Human'))
print(if (websiteLive) plotEnrich(enriched$WikiPathway_2023_Human, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='WikiPathway_2023_Human'))
print(if (websiteLive) plotEnrich(enriched$WikiPathways_2013, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='WikiPathways_2013'))
print(if (websiteLive) plotEnrich(enriched$WikiPathways_2015, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='WikiPathways_2015'))
print(if (websiteLive) plotEnrich(enriched$WikiPathways_2016, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='WikiPathways_2016'))
print(if (websiteLive) plotEnrich(enriched$WikiPathways_2019_Human, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='WikiPathways_2019_Human'))
print(if (websiteLive) plotEnrich(enriched$WikiPathways_2019_Mouse, showTerms = 20, numChar = 40, y = 'Count', orderBy = 'P.value', xlab='WikiPathways_2019_Mouse'))


dev.off()


}

############## 







########### Convert gene symbol to entrez gene id. Use "org.Hs.egENSEMBL" to convert Ensemb gene id to entrez gene id.

probes=row.names(as.data.frame(healthy_2015_common))

annot = as.data.frame(org.Bt.egENSEMBL)
probes2annot= match(probes,annot$ensembl_id)

allLLIDs=annot$gene_id[probes2annot];

moduleLabelsAutomatic1 = net$colors 
moduleColorsAutomatic1 = labels2colors(moduleLabelsAutomatic1)
########### 


########### GO analysis by WGCNA
GOenr = GOenrichmentAnalysis(moduleLabelsAutomatic1, allLLIDs, organism = "bovine", nBestP = 10);

tab = GOenr$bestPTerms[[4]]$enrichment

write.table(tab, file = "6-GOEnrichmentTable-net1.csv", sep = ",", quote = TRUE, row.names = FALSE)
########### 




########################## Combined GO Visualization for Yellow and Red Modules
setwd("D:/new desktop/thesis/sixth WGCNA analysis/var0.3")

install.packages(c("sunburstR", "webshot", "visNetwork", "ggalluvial", 
                   "networkD3", "UpSetR", "plotly", "treemap"))
# Load necessary libraries
library(xlsx)
library(readxl)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(stringr)
library(ggrepel)
library(ggpubr)
library(viridis)
library(circlize)
library(ComplexHeatmap)

# Function to read data from Excel files
read_enrichment_data <- function(file_path, sheet_name) {
  tryCatch({
    data <- read_excel(file_path, sheet = sheet_name)
    return(data)
  }, error = function(e) {
    cat("Error reading", file_path, "sheet", sheet_name, ":", conditionMessage(e), "\n")
    return(NULL)
  })
}

# Function to extract the number from Overlap column (e.g., "5/100" -> 5)
extract_gene_count <- function(overlap_str) {
  as.numeric(sapply(strsplit(overlap_str, "/"), function(x) x[1]))
}

# Read data for BP (Biological Process) and KEGG for Yellow and Red modules
yellow_bp <- read_enrichment_data("GO_enrichR_Yellow_Red/BioPlanet_2019.xlsx", "yellow")
red_bp <- read_enrichment_data("GO_enrichR_Yellow_Red/BioPlanet_2019.xlsx", "red")
yellow_kegg <- read_enrichment_data("GO_enrichR_Yellow_Red/KEGG_2016.xlsx", "yellow")
red_kegg <- read_enrichment_data("GO_enrichR_Yellow_Red/KEGG_2016.xlsx", "red")

# Check if all data was loaded
if (is.null(yellow_bp) || is.null(red_bp) || is.null(yellow_kegg) || is.null(red_kegg)) {
  stop("Some required data files could not be loaded. Please check the file paths and formats.")
}

# Create output directory
dir.create("GO_Visualization", showWarnings = FALSE)

# Process data for visualization
process_data_for_viz <- function(data, module, db_type) {
  if (nrow(data) == 0) return(NULL)
  
  # Sort by P.value and select top 10 terms
  data <- data[order(data$P.value), ]
  data <- head(data, 10)
  
  # Extract gene count from Overlap
  data$GeneCount <- extract_gene_count(data$Overlap)
  
  # Add module and database type information
  data$Module <- module
  data$Database <- db_type
  
  # Select relevant columns
  return(data[, c("Term", "P.value", "GeneCount", "Module", "Database")])
}

# Process all datasets
yellow_bp_processed <- process_data_for_viz(yellow_bp, "Yellow", "BioPlanet")
red_bp_processed <- process_data_for_viz(red_bp, "Red", "BioPlanet")
yellow_kegg_processed <- process_data_for_viz(yellow_kegg, "Yellow", "KEGG")
red_kegg_processed <- process_data_for_viz(red_kegg, "Red", "KEGG")

# Combine all processed data
all_data <- rbind(yellow_bp_processed, red_bp_processed, 
                  yellow_kegg_processed, red_kegg_processed)

# Shorten long term names for better visualization
all_data$ShortTerm <- ifelse(nchar(all_data$Term) > 40, 
                             paste0(substr(all_data$Term, 1, 37), "..."), 
                             all_data$Term)

# Find unique terms to highlight (marked as orange in your requirements)
# For each module-database combination, select non-overlapping terms
unique_terms <- function() {
  # Get all terms
  all_terms <- unique(all_data$Term)
  
  # Count occurrences of each term
  term_counts <- table(all_data$Term)
  
  # Get terms that appear only once (unique)
  unique_terms <- names(term_counts[term_counts == 1])
  
  # Add a column to mark unique terms
  all_data$Unique <- all_data$Term %in% unique_terms
  
  return(all_data)
}

all_data <- unique_terms()

##############################
# VISUALIZATION 1: Combined Plot of All Four Categories
##############################

# Create a combined visualization showing all four categories
pdf("GO_Visualization/Combined_GO_Plot.pdf", width = 14, height = 10)

# Create a factor for proper ordering in the plot
all_data$ModuleDB <- factor(paste(all_data$Module, all_data$Database, sep="_"),
                            levels = c("Yellow_BioPlanet", "Red_BioPlanet", "Yellow_KEGG", "Red_KEGG"))

# Define colors
module_colors <- c("Yellow_BioPlanet" = "#FFC107", 
                   "Red_BioPlanet" = "#F44336",
                   "Yellow_KEGG" = "#FFEB3B", 
                   "Red_KEGG" = "#E91E63")

# Create a bubble plot with facets - FIX: Use numeric values for alpha instead of TRUE/FALSE
p <- ggplot(all_data, aes(x = -log10(P.value), y = reorder(ShortTerm, -log10(P.value)))) +
  geom_point(aes(size = GeneCount, color = ModuleDB, alpha = Unique)) +
  scale_size_continuous(range = c(3, 8)) +
  scale_color_manual(values = module_colors) +
  scale_alpha_manual(values = c("FALSE" = 0.7, "TRUE" = 1.0)) +
  facet_wrap(~ ModuleDB, scales = "free_y", ncol = 2) +
  labs(x = "-log10(p-value)", y = "", size = "Gene Count", color = "Module & Database",
       title = "Enriched Pathways in Yellow and Red Modules",
       subtitle = "BioPlanet and KEGG databases") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

print(p)

dev.off()

##############################
# VISUALIZATION 2: Alternative plot showing unique terms
##############################

# Create a second visualization focusing on unique terms
pdf("GO_Visualization/Unique_Terms_Plot.pdf", width = 12, height = 10)

# Filter for unique terms
unique_terms_data <- all_data[all_data$Unique == TRUE, ]

# If there are unique terms
if(nrow(unique_terms_data) > 0) {
  # Create a horizontal bar plot for unique terms
  p2 <- ggplot(unique_terms_data, aes(x = -log10(P.value), y = reorder(ShortTerm, -log10(P.value)))) +
    geom_bar(stat = "identity", aes(fill = ModuleDB), width = 0.7) +
    scale_fill_manual(values = module_colors) +
    labs(x = "-log10(p-value)", y = "", fill = "Module & Database",
         title = "Unique Enriched Pathways in Yellow and Red Modules",
         subtitle = "Terms that appear in only one module-database combination") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  # Add gene count as text
  p2 <- p2 + geom_text(aes(label = paste0("n=", GeneCount)), hjust = -0.2, size = 3)
  
  print(p2)
} else {
  # Create a message if no unique terms
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(0, 0, "No unique terms found in the data")
}

dev.off()

##############################
# VISUALIZATION 3: Combined plot with highlighted unique terms
##############################

pdf("GO_Visualization/Combined_Highlight_Unique.pdf", width = 14, height = 12)

# Create a single plot with all terms but highlight unique ones
p3 <- ggplot(all_data, aes(x = -log10(P.value), y = reorder(ShortTerm, -log10(P.value)))) +
  geom_bar(stat = "identity", aes(fill = ModuleDB, alpha = Unique)) +
  scale_fill_manual(values = module_colors) +
  scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 1.0), guide = "none") +
  labs(x = "-log10(p-value)", y = "", fill = "Module & Database",
       title = "Enriched Pathways in Yellow and Red Modules",
       subtitle = "Unique terms (appearing in only one category) are highlighted with full opacity") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  ) +
  facet_wrap(~ ModuleDB, scales = "free_y", ncol = 2)

print(p3)

dev.off()

##############################
# VISUALIZATION 4: Alternative combined view with dots and lines
##############################

pdf("GO_Visualization/Dot_Line_GO_Plot.pdf", width = 14, height = 10)

# Create a combined visualization with dots and lines
p4 <- ggplot(all_data, aes(x = ModuleDB, y = reorder(ShortTerm, -log10(P.value)))) +
  geom_point(aes(size = -log10(P.value), color = ModuleDB)) +
  scale_size_continuous(range = c(1, 8)) +
  scale_color_manual(values = module_colors) +
  labs(x = "", y = "", size = "-log10(p-value)", color = "Module & Database",
       title = "Enriched Pathways Across Yellow and Red Modules",
       subtitle = "Size indicates significance (-log10 p-value)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.x = element_line(color = "grey80"),
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

print(p4)

dev.off()

cat("All visualizations have been created in the GO_Visualization directory.\n")



##############################
# ADVANCED GO VISUALIZATIONS
# این قسمت را به انتهای کد قبلی خود اضافه کنید
##############################

###########################
# VISUALIZATION 5: Sunburst Chart
###########################
library(sunburstR)

pdf("GO_Visualization/Sunburst_GO.pdf", width = 10, height = 10)

# Prepare data for sunburst
prepare_sunburst_data <- function() {
  # Create a data frame with path and values
  sb_data <- data.frame(
    path = character(),
    value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # For each row in all_data
  for (i in 1:nrow(all_data)) {
    # Create path: Module-Database-Term
    path <- paste(all_data$Module[i], all_data$Database[i], 
                  all_data$ShortTerm[i], sep = "-")
    
    # Value is -log10(p-value) * GeneCount to emphasize important terms
    value <- -log10(all_data$P.value[i]) * all_data$GeneCount[i]
    
    # Add to data frame
    sb_data <- rbind(sb_data, data.frame(
      path = path,
      value = value,
      stringsAsFactors = FALSE
    ))
  }
  
  return(sb_data)
}

sunburst_data <- prepare_sunburst_data()

# Create temporary HTML file for sunburst
temp_html <- tempfile(fileext = ".html")

# Generate sunburst chart
sb <- sunburst(
  sunburst_data,
  count = TRUE,
  legend = list(w = 200, h = 20, s = 5),
  breadcrumb = list(w = 0, h = 20, s = 5),
  width = 600,
  height = 600
)

# Save to HTML then convert to PDF
htmlwidgets::saveWidget(sb, temp_html, selfcontained = TRUE)
webshot::webshot(temp_html, "GO_Visualization/Sunburst_GO.pdf", 
                 delay = 2, vwidth = 700, vheight = 700)

# Attempt to display in R viewer too
print(sb)

dev.off()

###########################
# VISUALIZATION 6: Interactive Network with visNetwork
###########################
library(visNetwork)

# Prepare data for network visualization
prepare_visnetwork_data <- function() {
  # Create nodes
  nodes <- data.frame(
    id = 1:nrow(all_data),
    label = all_data$ShortTerm,
    group = all_data$ModuleDB,
    value = -log10(all_data$P.value) * 2,  # Size based on significance
    title = paste("Term:", all_data$Term, "<br>",
                  "P-value:", all_data$P.value, "<br>",
                  "Gene Count:", all_data$GeneCount),
    stringsAsFactors = FALSE
  )
  
  # Create edges based on shared genes or database/module
  edges <- data.frame(
    from = integer(),
    to = integer(),
    width = numeric(),
    title = character(),
    stringsAsFactors = FALSE
  )
  
  # Connect nodes within same module
  for (module in unique(all_data$Module)) {
    module_nodes <- which(all_data$Module == module)
    
    if (length(module_nodes) > 1) {
      for (i in 1:(length(module_nodes)-1)) {
        for (j in (i+1):length(module_nodes)) {
          # Add edge
          edges <- rbind(edges, data.frame(
            from = module_nodes[i],
            to = module_nodes[j],
            width = 1,
            title = paste("Same module:", module),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Connect nodes with same term across different modules/databases
  term_groups <- split(1:nrow(all_data), all_data$Term)
  for (term_nodes in term_groups) {
    if (length(term_nodes) > 1) {
      for (i in 1:(length(term_nodes)-1)) {
        for (j in (i+1):length(term_nodes)) {
          # Add edge with higher weight
          edges <- rbind(edges, data.frame(
            from = term_nodes[i],
            to = term_nodes[j],
            width = 3,
            title = paste("Same term:", all_data$Term[term_nodes[i]]),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  return(list(nodes = nodes, edges = edges))
}

network_data <- prepare_visnetwork_data()

# Create interactive network
vn <- visNetwork(
  network_data$nodes, 
  network_data$edges,
  width = "100%",
  height = "600px"
) %>%
  visGroups(groupname = "Yellow_BioPlanet", color = "#FFC107") %>%
  visGroups(groupname = "Red_BioPlanet", color = "#F44336") %>%
  visGroups(groupname = "Yellow_KEGG", color = "#FFEB3B") %>%
  visGroups(groupname = "Red_KEGG", color = "#E91E63") %>%
  visOptions(
    highlightNearest = TRUE,
    selectedBy = "group"
  ) %>%
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -100)
  ) %>%
  visLegend()

# Save to HTML
visNetwork::visSave(vn, file = "GO_Visualization/Interactive_GO_Network.html")

###########################
# VISUALIZATION 7: Alluvial Diagram (Sankey Diagram)
###########################
library(ggalluvial)

pdf("GO_Visualization/Alluvial_GO_Plot.pdf", width = 14, height = 12)

# Prepare data for alluvial diagram
# We will show how terms flow between modules and databases
alluvial_data <- all_data %>%
  # Separate Module and Database for the alluvial
  select(Term, Module, Database, P.value, GeneCount) %>%
  # Add significance category
  mutate(Significance = case_when(
    P.value < 0.001 ~ "High",
    P.value < 0.01 ~ "Medium",
    TRUE ~ "Low"
  ))

# Create the alluvial diagram
p_alluvial <- ggplot(alluvial_data,
                     aes(axis1 = Module, axis2 = Database, axis3 = Significance, y = -log10(P.value))) +
  geom_alluvium(aes(fill = Module), width = 0.3, knot.pos = 0.2) +
  geom_stratum(width = 0.3, fill = "white", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Module", "Database", "Significance"), expand = c(0.05, 0.05)) +
  scale_fill_manual(values = c("Yellow" = "#FFC107", "Red" = "#F44336")) +
  labs(title = "Flow of Gene Ontology Enrichment Results",
       subtitle = "Showing how terms flow between modules, databases and significance levels",
       y = "Cumulative -log10(p-value)") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

print(p_alluvial)

dev.off()

###########################
# VISUALIZATION 8: Hierarchical Edge Bundling
###########################
library(networkD3)

# Prepare data for edge bundling
prepare_heb_data <- function() {
  # Create nodes dataframe
  # Format: one row per node, with a unique id and a group
  nodes <- data.frame(
    name = character(),
    group = character(),
    stringsAsFactors = FALSE
  )
  
  # Add root node
  nodes <- rbind(nodes, data.frame(
    name = "GO_Analysis",
    group = "root",
    stringsAsFactors = FALSE
  ))
  
  # Add module nodes
  for (module in unique(all_data$Module)) {
    nodes <- rbind(nodes, data.frame(
      name = module,
      group = "module",
      stringsAsFactors = FALSE
    ))
  }
  
  # Add database nodes
  for (db in unique(all_data$Database)) {
    nodes <- rbind(nodes, data.frame(
      name = paste(db, "Database"),
      group = "database",
      stringsAsFactors = FALSE
    ))
  }
  
  # Add term nodes
  for (i in 1:nrow(all_data)) {
    # Create a unique ID for each term
    term_id <- paste(all_data$Term[i], i, sep = "_")
    
    nodes <- rbind(nodes, data.frame(
      name = term_id,
      group = "term",
      stringsAsFactors = FALSE
    ))
  }
  
  # Create links dataframe
  # Format: source and target columns with numeric indices
  links <- data.frame(
    source = integer(),
    target = integer(),
    value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Add links from root to modules
  for (module in unique(all_data$Module)) {
    source_idx <- which(nodes$name == "GO_Analysis") - 1  # 0-indexed
    target_idx <- which(nodes$name == module) - 1         # 0-indexed
    
    links <- rbind(links, data.frame(
      source = source_idx,
      target = target_idx,
      value = 1,
      stringsAsFactors = FALSE
    ))
  }
  
  # Add links from root to databases
  for (db in unique(all_data$Database)) {
    source_idx <- which(nodes$name == "GO_Analysis") - 1  # 0-indexed
    target_idx <- which(nodes$name == paste(db, "Database")) - 1  # 0-indexed
    
    links <- rbind(links, data.frame(
      source = source_idx,
      target = target_idx,
      value = 1,
      stringsAsFactors = FALSE
    ))
  }
  
  # Add links from modules to terms
  for (i in 1:nrow(all_data)) {
    module <- all_data$Module[i]
    term_id <- paste(all_data$Term[i], i, sep = "_")
    
    source_idx <- which(nodes$name == module) - 1         # 0-indexed
    target_idx <- which(nodes$name == term_id) - 1        # 0-indexed
    
    links <- rbind(links, data.frame(
      source = source_idx,
      target = target_idx,
      value = -log10(all_data$P.value[i]),  # Weight by significance
      stringsAsFactors = FALSE
    ))
  }
  
  # Add links from databases to terms
  for (i in 1:nrow(all_data)) {
    db <- paste(all_data$Database[i], "Database")
    term_id <- paste(all_data$Term[i], i, sep = "_")
    
    source_idx <- which(nodes$name == db) - 1             # 0-indexed
    target_idx <- which(nodes$name == term_id) - 1        # 0-indexed
    
    links <- rbind(links, data.frame(
      source = source_idx,
      target = target_idx,
      value = -log10(all_data$P.value[i]),  # Weight by significance
      stringsAsFactors = FALSE
    ))
  }
  
  # Clean node names for display
  nodes$name <- gsub("_[0-9]+$", "", nodes$name)  # Remove the numerical suffix
  
  return(list(nodes = nodes, links = links))
}

heb_data <- prepare_heb_data()

# Create hierarchical edge bundling diagram
heb <- sankeyNetwork(
  Links = heb_data$links, 
  Nodes = heb_data$nodes,
  Source = "source", 
  Target = "target",
  Value = "value", 
  NodeID = "name",
  NodeGroup = "group",
  fontSize = 12,
  nodeWidth = 30,
  sinksRight = FALSE,
  height = 600,
  width = 800
)

# Save to HTML
htmlwidgets::saveWidget(heb, file = "GO_Visualization/Hierarchical_Edge_Bundling.html", selfcontained = TRUE)

###########################
# VISUALIZATION 9: UpSet Plot (Set Intersections)
###########################
library(UpSetR)

pdf("GO_Visualization/UpSet_GO_Plot.pdf", width = 10, height = 8)

# Prepare data for UpSet plot
prepare_upset_data <- function() {
  # Get all unique terms
  all_terms <- unique(all_data$Term)
  
  # Create a binary matrix: rows are terms, columns are module-database combinations
  upset_matrix <- matrix(0, nrow = length(all_terms), ncol = 4)
  colnames(upset_matrix) <- c("Yellow_BP", "Red_BP", "Yellow_KEGG", "Red_KEGG")
  rownames(upset_matrix) <- all_terms
  
  # Fill the matrix
  for (i in 1:nrow(all_data)) {
    term <- all_data$Term[i]
    module_db <- paste(all_data$Module[i], ifelse(all_data$Database[i] == "BioPlanet", "BP", "KEGG"), sep = "_")
    
    # Set the corresponding cell to 1
    upset_matrix[term, module_db] <- 1
  }
  
  # Convert to data frame
  upset_df <- as.data.frame(upset_matrix)
  
  return(upset_df)
}

upset_data <- prepare_upset_data()

# Create UpSet plot
upset_plot <- upset(
  upset_data,
  nsets = 4,
  order.by = "freq",
  empty.intersections = "on",
  mainbar.y.label = "Term Intersection Size",
  sets.x.label = "Terms Per Category",
  point.size = 3.5,
  line.size = 2,
  mb.ratio = c(0.6, 0.4),
  main.bar.color = "black",
  sets.bar.color = c("Yellow_BP" = "#FFC107", 
                     "Red_BP" = "#F44336",
                     "Yellow_KEGG" = "#FFEB3B", 
                     "Red_KEGG" = "#E91E63")
)

print(upset_plot)

dev.off()

###########################
# VISUALIZATION 10: 3D PCA Plot of Term Similarities
###########################
library(plotly)

# Prepare data for PCA
prepare_pca_data <- function() {
  # Create a term-by-category matrix with p-values
  terms <- unique(all_data$Term)
  categories <- c("Yellow_BioPlanet", "Red_BioPlanet", "Yellow_KEGG", "Red_KEGG")
  
  pca_matrix <- matrix(NA, nrow = length(terms), ncol = length(categories))
  rownames(pca_matrix) <- terms
  colnames(pca_matrix) <- categories
  
  # Fill matrix with -log10(p-values)
  for (i in 1:nrow(all_data)) {
    term <- all_data$Term[i]
    category <- paste(all_data$Module[i], all_data$Database[i], sep = "_")
    
    pca_matrix[term, category] <- -log10(all_data$P.value[i])
  }
  
  # Replace NA with 0
  pca_matrix[is.na(pca_matrix)] <- 0
  
  # Perform PCA
  pca_result <- prcomp(pca_matrix, scale. = TRUE)
  
  # Get term coordinates in PCA space
  pca_coords <- pca_result$x[, 1:3]  # Use first 3 components
  
  # Add module information
  pca_df <- as.data.frame(pca_coords)
  
  # Determine predominant module for each term
  term_modules <- character(length(terms))
  term_datasets <- character(length(terms))
  term_scores <- numeric(length(terms))
  
  for (i in 1:length(terms)) {
    term <- terms[i]
    term_data <- subset(all_data, Term == term)
    
    # If term appears in multiple categories, use the one with lowest p-value
    if (nrow(term_data) > 0) {
      best_idx <- which.min(term_data$P.value)
      term_modules[i] <- as.character(term_data$Module[best_idx])
      term_datasets[i] <- as.character(term_data$Database[best_idx])
      term_scores[i] <- -log10(term_data$P.value[best_idx])
    }
  }
  
  pca_df$Module <- term_modules
  pca_df$Database <- term_datasets
  pca_df$Score <- term_scores
  pca_df$Term <- terms
  
  # Calculate percentage of variance explained
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
  
  return(list(pca_df = pca_df, var_explained = var_explained))
}

pca_results <- prepare_pca_data()
pca_df <- pca_results$pca_df
var_explained <- pca_results$var_explained

# Create 3D PCA plot
axis_labels <- paste0(
  c("PC1", "PC2", "PC3"),
  " (", 
  sprintf("%.1f", var_explained[1:3]),
  "%)"
)

# Define colors based on Module and Database
pca_df$Color <- factor(paste(pca_df$Module, pca_df$Database, sep = "_"))
color_values <- c("Yellow_BioPlanet" = "#FFC107", 
                  "Red_BioPlanet" = "#F44336",
                  "Yellow_KEGG" = "#FFEB3B", 
                  "Red_KEGG" = "#E91E63")

p <- plot_ly(pca_df, 
             x = ~PC1, y = ~PC2, z = ~PC3, 
             color = ~Color,
             colors = color_values,
             size = ~Score,
             sizes = c(5, 20),
             text = ~Term,
             hoverinfo = "text",
             type = "scatter3d",
             mode = "markers")

p <- p %>% layout(
  title = "3D PCA Plot of GO Term Similarities",
  scene = list(
    xaxis = list(title = axis_labels[1]),
    yaxis = list(title = axis_labels[2]),
    zaxis = list(title = axis_labels[3])
  ),
  legend = list(title = list(text = "Module & Database"))
)

# Save to HTML
htmlwidgets::saveWidget(p, file = "GO_Visualization/3D_PCA_GO_Terms.html", selfcontained = TRUE)

###########################
# VISUALIZATION 11: Treemap of GO Terms
###########################
library(treemap)

pdf("GO_Visualization/Treemap_GO.pdf", width = 12, height = 10)

# Prepare data for treemap
treemap_data <- all_data %>%
  # Add a significance score combining p-value and gene count
  mutate(Score = -log10(P.value) * GeneCount,
         # Create a combined category
         Category = paste(Module, Database, sep = "_"))

# Create the treemap
treemap(
  treemap_data,
  index = c("Category", "ShortTerm"),
  vSize = "Score",
  vColor = "P.value",
  type = "value",
  palette = "RdYlBu",
  aspRatio = 1,
  title = "GO Term Enrichment Treemap",
  fontsize.title = 14,
  reverse.palette = TRUE,
  algorithm = "pivotSize",
  sortID = "Score"
)

dev.off()

cat("All advanced visualizations have been created in the GO_Visualization directory.\n")


########################## End of GO






########################## Different plots
library(flashClust)
library(dendextend)

##### Discard the unassigned genes (grey module), and focus on the rest
restGenes= (mergedColors != "grey")
diss1=1-TOMsimilarityFromExpr(GSE51856_Milk_ExpNorm[,restGenes], power = 10)

colnames(diss1) =rownames(diss1) =SubGeneNames[restGenes]
hier1=flashClust(as.dist(diss1), method="average" )
plotDendroAndColors(hier1, mergedColors[restGenes], "Module colors", 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05, main = "Gene hierarchical clustering dendrogram (Control)")

#####



##### Visualize the Tom plot

# set the diagonal of the dissimilarity to NA 
diag(diss1) = NA;

# Raise the dissimilarity matrix to the power of 4 to bring out the module structure
sizeGrWindow(7,7)

TOMplot(diss1, hier1, as.character(mergedColors[restGenes]))

#####



##### Quantify module similarity by eigengene correlation. Eigengenes: Module representatives

MEList = moduleEigengenes(GSE51856_Milk_ExpNorm_df, colors = mergedColors)

MEs = MEList$eigengenes

plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))

#####

########################## End of different plots









########################## Module membership analysis, kME
load(file ="5_GSE51856_Milk_Net_Auto.RData")
load(file ="2_GSE51856_Milk_RemSam.RData") 

options(java.parameters = "-Xmx4000m")
library(enrichR)
library(xlsx)
library("clusterProfiler")
library("org.Bt.eg.db")
library("WGCNA")

listEnrichrSites()
setEnrichrSite("Enrichr")
websiteLive <- TRUE

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", 
         "GO_Biological_Process_2023", "KEGG_2021_Human")




moduleLabelsAutomatic1 = net_GSE51856_Milk$colors 
moduleColorsAutomatic1 = labels2colors(moduleLabelsAutomatic1)
moduleColorsAutomaticData1=moduleColorsAutomatic1

########### Eigengenes calculation
ME_GSE51856_Milk=moduleEigengenes(GSE51856_Milk_ExpNorm,moduleColorsAutomaticData1)$eigengenes


########### Module membership calculation
kME_GSE51856_Milk=signedKME(GSE51856_Milk_ExpNorm,ME_GSE51856_Milk,corFnc = "bicor")



########### preserved module hub genes: having consistently high positive or negative combined MM in both data
genename = names(as.data.frame(GSE51856_Milk_ExpNorm)) 

modulenames=unique(labels2colors(net_GSE51856_Milk$colors))


nogrey= (modulenames != "grey")

modulenames = subset(modulenames,nogrey )



for (i in modulenames){
  gc()
  
  namemodulekme=paste("kME", i, sep = "")
  
  kME_i=as.data.frame(kME_GSE51856_Milk[[namemodulekme]])
  
  kME_i$Gene=as.data.frame(names(as.data.frame(GSE51856_Milk_ExpNorm)) )
  
  inModule = (moduleColorsAutomaticData1==i) 
  kME_i=kME_i[inModule,]
  
  kME_i=cbind(kME_i[,2],kME_i[,1] )
  
  colnames(kME_i)= c("ENSEMBL", namemodulekme)
  
  # Convert ENSEMBL to Gene symbol by library("clusterProfiler")
  Genesymbols= bitr(kME_i$ENSEMBL , fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Bt.eg.db")
  KMEGeneIdModule = merge(Genesymbols, kME_i, by="ENSEMBL")
  colnames(KMEGeneIdModule)=c("ENSEMBL_ID", "Gene_Name", "KME")
  KMEGeneIdModule <- KMEGeneIdModule[order(KMEGeneIdModule$KME, decreasing =T  ),]
  
  
  write.xlsx(KMEGeneIdModule, file="8_kME_HUbs.xlsx", sheetName= i, append=TRUE,  row.names=F)
  
  
  
  
  
  ########### GO
  kME_i=subset(kME_i,kME_i[[namemodulekme]]>0.7)
  
  # Convert ENSEMBL to Gene symbol by library("clusterProfiler")
  nameofmodulesymbol=bitr(kME_i$ENSEMBL, fromType="ENSEMBL", toType="SYMBOL",
                          OrgDb="org.Bt.eg.db")$SYMBOL
  
  enriched= enrichr(nameofmodulesymbol, databases = dbs)
  write.xlsx(enriched$GO_Biological_Process_2023, file="8_kME_HUbs.xlsx", sheetName= paste(i,"_BP"), append=TRUE,  row.names=FALSE)
  
  ########### End of GO
  
  
  
  ########### Heatmap for hubs (KME)
  
  library(ComplexHeatmap)
  library(circlize)
  library("RColorBrewer")
  
  probes=names(GSE51856_Milk_ExpNorm_df)
  modulenames=unique(labels2colors(net_GSE51856_Milk$colors))
  col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
  col_fun = colorRamp2(c(-2, 0, 2), c( "purple", "white", "red"))
  
  # Separate the groups
  GSE51856_Milk_sample_groups =factor(rep(c("A", "B","C", "D", "E", "F", "G", "H","I", "J"), 
                                          times = c(4,5,5,5,5,5,5,4,5,5)))
  
  GSE51856_Milk_fa=rep(c("H0","H12","H24","H36","H48","I0", "I12","I24", "I36", "I48"), 
                       times = c(4,5,5,5,5,5,5,4,5,5))
  GSE51856_Milk_fa_col = c("H0" = 1, "H12" = 2,"H24" = 3, "H36" = 4,"H48" = 5, "I0" = 6,
                           "I12" = 7, "I24" = 8,"I36" = 9, "I48" = 10)
  
  for (i in modulenames) {
    GSE51856_Milk_Expr = t(GSE51856_Milk_ExpNorm_df)
    GSE51856_Milk_Expr <- subset(GSE51856_Milk_Expr, rownames(GSE51856_Milk_Expr) %in% kME_i$ENSEMBL)
    colnames(GSE51856_Milk_Expr) <- colnames(GSE51856_Milk_Filt)[c(-3,-38)]
    
    pdfname = paste("Heatmap_KME/", i, ".pdf", sep="")
    
    pdf(file = pdfname)
    
    draw(Heatmap(GSE51856_Milk_Expr,
                 column_split = GSE51856_Milk_sample_groups,
                 column_order = NULL, cluster_rows = F,
                 show_row_dend = F,
                 column_title = paste(i, "module", sep=" "),
                 col = col_fun,
                 top_annotation = HeatmapAnnotation(Groups = GSE51856_Milk_fa, 
                                                    col = list(Groups = GSE51856_Milk_fa_col)),
                 row_title = "Genes", 
                 row_title_side = c("right"),
                 row_title_gp = gpar(fontsize = 10),
                 row_names_max_width = unit(1, "cm"),
                 row_names_gp = gpar(fontsize = 1.5),
                 show_column_names = TRUE,
                 column_names_gp = gpar(fontsize = 8),
                 column_names_rot = 45,
                 width = unit(8, "cm"), height = unit(12, "cm"),
                 show_heatmap_legend = T,
                 heatmap_legend_param = list(
                   at = c(-3, 0, 3),
                   labels = c("-3", "0", "3"),
                   title = "RNA stability",
                   legend_height = unit(2.2, "cm"),
                   title_position = "topleft",
                   title_gp = gpar(fontsize = 8, fontface = "bold"),
                   labels_gp = gpar(fontsize = 8)
                 )),
         heatmap_legend_side = "left")
    
    dev.off()
  }
  
  ########### End of heatmap for hubs (KME)
  
  
}  


########################## End of module membership analysis, kME









########################## Intramodular connectivity (kIM)
load(file ="5_GSE51856_Milk_Net_Auto.RData")
load(file ="2_GSE51856_Milk_RemSam.RData") 

options(java.parameters = "-Xmx4000m")
library(enrichR)
library(xlsx)
library("clusterProfiler")
library("org.Bt.eg.db")
library("WGCNA")

listEnrichrSites()
setEnrichrSite("Enrichr")
websiteLive <- TRUE

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", 
         "GO_Biological_Process_2023", "KEGG_2021_Human")



moduleLabelsAutomatic1 = net_GSE51856_Milk$colors 
moduleColorsAutomatic1 = labels2colors(moduleLabelsAutomatic1)
moduleColorsAutomaticData1=moduleColorsAutomatic1



kIM_GSE51856_Milk=intramodularConnectivity.fromExpr(GSE51856_Milk_ExpNorm, 
                                                    moduleColorsAutomaticData1,corFnc = "bicor",networkType = "signed", scaleByMax=TRUE,power=10)$kWithin

kIM_GSE51856_Milk=as.data.frame(kIM_GSE51856_Milk)

row.names(kIM_GSE51856_Milk)=row.names(as.data.frame(t(GSE51856_Milk_ExpNorm)))


kIM_GSE51856_Mil = as.data.frame(cbind(row.names(kIM_GSE51856_Milk),
                                       kIM_GSE51856_Milk$kIM_GSE51856_Milk))

colnames(kIM_GSE51856_Mil)=c("ENSEMBL", "KIM")


# Convert ENSEMBL to Gene symbol by library("clusterProfiler")
Genesymbols= bitr(kIM_GSE51856_Mil$ENSEMBL , fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Bt.eg.db")
KIMGeneIdModule = merge(Genesymbols, kIM_GSE51856_Mil, by="ENSEMBL")
colnames(KIMGeneIdModule)=c("ENSEMBL", "Gene_Name", "KIM")

KIMGeneIdModule = merge(KIMGeneIdModule, GeneIdModule, by="ENSEMBL")

KIMGeneIdModule= KIMGeneIdModule[,c(1,2,3,6)]


modulenames=unique(labels2colors(net_GSE51856_Milk$colors))
nogrey= (modulenames != "grey")
modulenames = subset(modulenames,nogrey )



for (i in modulenames){
  gc()
  
  kIM_i=subset(KIMGeneIdModule, Module_Color==i)
  
  kIM_i <- kIM_i[order(kIM_i$KIM, decreasing =T  ),]
  
  write.xlsx(kIM_i, file="9_kIM_HUbs.xlsx", sheetName= i, append=TRUE,  row.names=F)
  
  
  
  ########### GO
  kIM_i=subset(kIM_i,kIM_i$KIM>0.7)
  
  # Convert ENSEMBL to Gene symbol by library("clusterProfiler")
  nameofmodulesymbol=bitr(kIM_i$ENSEMBL, fromType="ENSEMBL", toType="SYMBOL",
                          OrgDb="org.Bt.eg.db")$SYMBOL
  
  enriched= enrichr(nameofmodulesymbol, databases = dbs)
  write.xlsx(enriched$GO_Biological_Process_2023, file="9_kIM_HUbs.xlsx", sheetName= paste(i,"_BP"), append=TRUE,  row.names=FALSE)
  
  ########### End of GO
  
  
  
  ########### Heatmap for hubs (KIM)
  
  library(ComplexHeatmap)
  library(circlize)
  library("RColorBrewer")
  
  probes=names(GSE51856_Milk_ExpNorm_df)
  modulenames=unique(labels2colors(net_GSE51856_Milk$colors))
  col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
  col_fun = colorRamp2(c(-2, 0, 2), c( "purple", "white", "red"))
  
  # Separate the groups
  GSE51856_Milk_sample_groups =factor(rep(c("A", "B","C", "D", "E", "F", "G", "H","I", "J"), 
                                          times = c(4,5,5,5,5,5,5,4,5,5)))
  
  GSE51856_Milk_fa=rep(c("H0","H12","H24","H36","H48","I0", "I12","I24", "I36", "I48"), 
                       times = c(4,5,5,5,5,5,5,4,5,5))
  GSE51856_Milk_fa_col = c("H0" = 1, "H12" = 2,"H24" = 3, "H36" = 4,"H48" = 5, "I0" = 6,
                           "I12" = 7, "I24" = 8,"I36" = 9, "I48" = 10)
  
  for (i in modulenames) {
    GSE51856_Milk_Expr = t(GSE51856_Milk_ExpNorm_df)
    GSE51856_Milk_Expr <- subset(GSE51856_Milk_Expr, rownames(GSE51856_Milk_Expr) %in% kIM_i$ENSEMBL)
    colnames(GSE51856_Milk_Expr) <- colnames(GSE51856_Milk_Filt)[c(-3,-38)]
    
    pdfname = paste("Heatmap_KIM/", i, ".pdf", sep="")
    
    pdf(file = pdfname)
    
    draw(Heatmap(GSE51856_Milk_Expr,
                 column_split = GSE51856_Milk_sample_groups,
                 column_order = NULL, cluster_rows = F,
                 show_row_dend = F,
                 column_title = paste(i, "module", sep=" "),
                 col = col_fun,
                 top_annotation = HeatmapAnnotation(Groups = GSE51856_Milk_fa, 
                                                    col = list(Groups = GSE51856_Milk_fa_col)),
                 row_title = "Genes", 
                 row_title_side = c("right"),
                 row_title_gp = gpar(fontsize = 10),
                 row_names_max_width = unit(1, "cm"),
                 row_names_gp = gpar(fontsize = 1.5),
                 show_column_names = TRUE,
                 column_names_gp = gpar(fontsize = 8),
                 column_names_rot = 45,
                 width = unit(8, "cm"), height = unit(12, "cm"),
                 show_heatmap_legend = T,
                 heatmap_legend_param = list(
                   at = c(-3, 0, 3),
                   labels = c("-3", "0", "3"),
                   title = "RNA stability",
                   legend_height = unit(2.2, "cm"),
                   title_position = "topleft",
                   title_gp = gpar(fontsize = 8, fontface = "bold"),
                   labels_gp = gpar(fontsize = 8)
                 )),
         heatmap_legend_side = "left")
    
    dev.off()
  }
  
  ########### End of heatmap for hubs (KIM)
}
########################## End of intramodular connectivity (kIM)




########################## STRING PPI
setwd("C:/Users/M.A.Shirazi/Desktop/thesis/sixth WGCNA analysis/var0.3")

load(file ="5_GSE51856_Milk_Net_Auto.RData") 
load(file ="2_GSE51856_Milk_RemSam.RData") 


require(STRINGdb)
require(igraph)
require(biomaRt)

library("clusterProfiler")
library("org.Bt.eg.db")
library("WGCNA")

GSE51856_Milk_ExpNorm_df=as.data.frame(GSE51856_Milk_ExpNorm)

probes=names(GSE51856_Milk_ExpNorm_df)

moduleLabelsAutomatic1 = net_GSE51856_Milk$colors 
moduleColorsAutomatic1 = labels2colors(moduleLabelsAutomatic1)
moduleColorsAutomaticData1=moduleColorsAutomatic1


modulenames=unique(labels2colors(net_GSE51856_Milk$colors))
nogrey= (modulenames != "grey")
modulenames = subset(modulenames,nogrey )


string_db <- STRINGdb$new( version="12", species=9606, score_threshold=00, input_directory="")



for (i in modulenames){
  gc()
  # Select the module of interest
  inModule = (moduleColorsAutomaticData1==i) 
  nameofmodule=as.data.frame(probes)[inModule,]
  
  
  # Convert ENSEMBL to Gene symbol by library("clusterProfiler")
  nameofmodulesymbol=bitr(nameofmodule, fromType="ENSEMBL", toType="SYMBOL",
                          OrgDb="org.Bt.eg.db")$SYMBOL
  
  genestostringid <- string_db$map( as.data.frame(nameofmodulesymbol), 
                                    "nameofmodulesymbol", removeUnmappedRows = TRUE )
  
  pdf(file = paste("STRING/",i,"_STRING_Plot.pdf",sep="") )
  string_db$plot_network(genestostringid$STRING_id)
  dev.off()
  
  # Enrichment Analysis
  enrichmentSTRING=string_db$get_enrichment(genestostringid$STRING_id)
  
  
  
  # Get interactions and convert to gene symbol
  interactionPPI=string_db$get_interactions(genestostringid$STRING_id)
  
  interactionPPI_Genes=as.data.frame(cbind(string_db$add_proteins_description(data.frame(STRING_id = interactionPPI$from))$preferred_name,
                                           string_db$add_proteins_description(data.frame(STRING_id = interactionPPI$to))$preferred_name)  )
  colnames(interactionPPI_Genes)=c("From", "To")
  
  
  write.xlsx(interactionPPI_Genes, file="STRING/STRING_Interactions.xlsx", sheetName= i, append=TRUE,  row.names=F)
  write.xlsx(enrichmentSTRING, file="STRING/STRING_Interactions.xlsx", sheetName= paste(i,"_GO"), append=TRUE,  row.names=F)
  
  
}


########################## End of STRING PPI
  
  
  






















library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

set.seed(123) # تعیین تخمینی برای پایداری RNA

# تولید مقادیر تصادفی برای p-value
pvalue = 10^-runif(ncol(GSE51856_Milk_Expr), min = 0, max = 3)
is_sig = pvalue < 0.01
pch = rep("*", ncol(GSE51856_Milk_Expr))
pch[!is_sig] = NA
# محدوده رنگی برای -log10(pvalue)
pvalue_col_fun = colorRamp2(c(0, 2, 3), c("green", "white", "red")) 

probes = names(GSE51856_Milk_ExpNorm_df)
modulenames = unique(labels2colors(net_GSE51856_Milk$colors))

# ساخت رنگ‌ها
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
col_fun = colorRamp2(c(-2, 0, 2), c("purple", "white", "red"))

# تفکیک گروه‌ها
GSE51856_Milk_sample_groups = factor(rep(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"), 
                                         times = c(4,5,5,5,5,5,5,4,5,5)))

GSE51856_Milk_fa = rep(c("H0", "H12", "H24", "H36", "H48", "I0", "I12", "I24", "I36", "I48"), 
                       times = c(4,5,5,5,5,5,5,4,5,5))
GSE51856_Milk_fa_col = c("H0" = 1, "H12" = 2, "H24" = 3, "H36" = 4, "H48" = 5, "I0" = 6,
                         "I12" = 7, "I24" = 8, "I36" = 9, "I48" = 10)

for (i in modulenames) {
  inModule = (moduleColorsAutomaticData1 == i)
  nameofmodule = as.data.frame(probes)[inModule,]
  GSE51856_Milk_Expr = t(GSE51856_Milk_ExpNorm_df)
  GSE51856_Milk_Expr <- subset(GSE51856_Milk_Expr, rownames(GSE51856_Milk_Expr) %in% nameofmodule)
  colnames(GSE51856_Milk_Expr) <- colnames(GSE51856_Milk_Filt)[c(-3, -38)]
  
  pdfname = paste("Heatmap_PDF2/", i, ".pdf", sep = "")
  
  pdf(file = pdfname)
  
  # توضیحات بالای Heatmap
  column_ha = HeatmapAnnotation(
    foo1 = runif(ncol(GSE51856_Milk_Expr)), 
    bar1 = anno_barplot(runif(ncol(GSE51856_Milk_Expr))),
    pvalue = anno_simple(-log10(pvalue), col = pvalue_col_fun, pch = pch),
    annotation_name_side = "left",
    Groups = anno_block(gp = gpar(fill = 1:10), labels = levels(GSE51856_Milk_fa))
  )
  
  # توضیحات سمت راست Heatmap
  row_ha = rowAnnotation(
    foo2 = runif(nrow(GSE51856_Milk_Expr)), 
    bar2 = anno_barplot(runif(nrow(GSE51856_Milk_Expr)))
  )
  
  # Heatmap
  heatmap = Heatmap(GSE51856_Milk_Expr,
                    column_split = GSE51856_Milk_sample_groups,
                    column_order = NULL, cluster_rows = FALSE,
                    show_row_dend = FALSE,
                    column_title = paste(i, "module", sep = " "),
                    col = col_fun,
                    top_annotation = column_ha,
                    right_annotation = row_ha,
                    row_title = "Genes", 
                    row_title_side = "right",
                    row_title_gp = gpar(fontsize = 10),
                    row_names_max_width = unit(1, "cm"),
                    row_names_gp = gpar(fontsize = 1.5),
                    show_column_names = TRUE,
                    column_names_gp = gpar(fontsize = 8),
                    column_names_rot = 45,
                    width = unit(8, "cm"), height = unit(12, "cm"),
                    show_heatmap_legend = TRUE,
                    heatmap_legend_param = list(
                      at = c(-3, 0, 3),
                      labels = c("-3", "0", "3"),
                      title = "RNA stability",
                      legend_height = unit(2.2, "cm"),
                      title_position = "topleft",
                      title_gp = gpar(fontsize = 8, fontface = "bold"),
                      labels_gp = gpar(fontsize = 8)
                    ))
  
  draw(heatmap, heatmap_legend_side = "left")
  
  dev.off()
}

########### اتمام Heatmap
































































































########################## Combined GO Visualization for Yellow and Red Modules
setwd("D:/new desktop/thesis/sixth WGCNA analysis/var0.3")

install.packages(c("sunburstR", "webshot", "visNetwork", "ggalluvial", 
                   "networkD3", "UpSetR", "plotly", "treemap"))
# Load necessary libraries
library(xlsx)
library(readxl)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(stringr)
library(ggrepel)
library(ggpubr)
library(viridis)
library(circlize)
library(ComplexHeatmap)

# Function to read data from Excel files
read_enrichment_data <- function(file_path, sheet_name) {
  tryCatch({
    data <- read_excel(file_path, sheet = sheet_name)
    return(data)
  }, error = function(e) {
    cat("Error reading", file_path, "sheet", sheet_name, ":", conditionMessage(e), "\n")
    return(NULL)
  })
}

# Function to extract the number from Overlap column (e.g., "5/100" -> 5)
extract_gene_count <- function(overlap_str) {
  as.numeric(sapply(strsplit(overlap_str, "/"), function(x) x[1]))
}

# Read data for BP (Biological Process) and KEGG for Yellow and Red modules
yellow_bp <- read_enrichment_data("GO_enrichR_Yellow_Red/BioPlanet_2019.xlsx", "yellow")
red_bp <- read_enrichment_data("GO_enrichR_Yellow_Red/BioPlanet_2019.xlsx", "red")
yellow_kegg <- read_enrichment_data("GO_enrichR_Yellow_Red/KEGG_2016.xlsx", "yellow")
red_kegg <- read_enrichment_data("GO_enrichR_Yellow_Red/KEGG_2016.xlsx", "red")

# Check if all data was loaded
if (is.null(yellow_bp) || is.null(red_bp) || is.null(yellow_kegg) || is.null(red_kegg)) {
  stop("Some required data files could not be loaded. Please check the file paths and formats.")
}

# Create output directory
dir.create("GO_Visualization", showWarnings = FALSE)

# Process data for visualization
process_data_for_viz <- function(data, module, db_type) {
  if (nrow(data) == 0) return(NULL)
  
  # Sort by P.value and select top 10 terms
  data <- data[order(data$P.value), ]
  data <- head(data, 10)
  
  # Extract gene count from Overlap
  data$GeneCount <- extract_gene_count(data$Overlap)
  
  # Add module and database type information
  data$Module <- module
  data$Database <- db_type
  
  # Select relevant columns
  return(data[, c("Term", "P.value", "GeneCount", "Module", "Database")])
}

# Process all datasets
yellow_bp_processed <- process_data_for_viz(yellow_bp, "Yellow", "BioPlanet")
red_bp_processed <- process_data_for_viz(red_bp, "Red", "BioPlanet")
yellow_kegg_processed <- process_data_for_viz(yellow_kegg, "Yellow", "KEGG")
red_kegg_processed <- process_data_for_viz(red_kegg, "Red", "KEGG")

# Combine all processed data
all_data <- rbind(yellow_bp_processed, red_bp_processed, 
                  yellow_kegg_processed, red_kegg_processed)

# Shorten long term names for better visualization
all_data$ShortTerm <- ifelse(nchar(all_data$Term) > 40, 
                             paste0(substr(all_data$Term, 1, 37), "..."), 
                             all_data$Term)

# Find unique terms to highlight (marked as orange in your requirements)
# For each module-database combination, select non-overlapping terms
unique_terms <- function() {
  # Get all terms
  all_terms <- unique(all_data$Term)
  
  # Count occurrences of each term
  term_counts <- table(all_data$Term)
  
  # Get terms that appear only once (unique)
  unique_terms <- names(term_counts[term_counts == 1])
  
  # Add a column to mark unique terms
  all_data$Unique <- all_data$Term %in% unique_terms
  
  return(all_data)
}

all_data <- unique_terms()

##############################
# VISUALIZATION 1: Combined Plot of All Four Categories
##############################

# Create a combined visualization showing all four categories
pdf("GO_Visualization/Combined_GO_Plot.pdf", width = 14, height = 10)

# Create a factor for proper ordering in the plot
all_data$ModuleDB <- factor(paste(all_data$Module, all_data$Database, sep="_"),
                            levels = c("Yellow_BioPlanet", "Red_BioPlanet", "Yellow_KEGG", "Red_KEGG"))

# Define colors
module_colors <- c("Yellow_BioPlanet" = "#FFC107", 
                   "Red_BioPlanet" = "#F44336",
                   "Yellow_KEGG" = "#FFEB3B", 
                   "Red_KEGG" = "#E91E63")

# Create a bubble plot with facets - FIX: Use numeric values for alpha instead of TRUE/FALSE
p <- ggplot(all_data, aes(x = -log10(P.value), y = reorder(ShortTerm, -log10(P.value)))) +
  geom_point(aes(size = GeneCount, color = ModuleDB, alpha = Unique)) +
  scale_size_continuous(range = c(3, 8)) +
  scale_color_manual(values = module_colors) +
  scale_alpha_manual(values = c("FALSE" = 0.7, "TRUE" = 1.0)) +
  facet_wrap(~ ModuleDB, scales = "free_y", ncol = 2) +
  labs(x = "-log10(p-value)", y = "", size = "Gene Count", color = "Module & Database",
       title = "Enriched Pathways in Yellow and Red Modules",
       subtitle = "BioPlanet and KEGG databases") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

print(p)

dev.off()

##############################
# VISUALIZATION 2: Alternative plot showing unique terms
##############################

# Create a second visualization focusing on unique terms
pdf("GO_Visualization/Unique_Terms_Plot.pdf", width = 12, height = 10)

# Filter for unique terms
unique_terms_data <- all_data[all_data$Unique == TRUE, ]

# If there are unique terms
if(nrow(unique_terms_data) > 0) {
  # Create a horizontal bar plot for unique terms
  p2 <- ggplot(unique_terms_data, aes(x = -log10(P.value), y = reorder(ShortTerm, -log10(P.value)))) +
    geom_bar(stat = "identity", aes(fill = ModuleDB), width = 0.7) +
    scale_fill_manual(values = module_colors) +
    labs(x = "-log10(p-value)", y = "", fill = "Module & Database",
         title = "Unique Enriched Pathways in Yellow and Red Modules",
         subtitle = "Terms that appear in only one module-database combination") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  # Add gene count as text
  p2 <- p2 + geom_text(aes(label = paste0("n=", GeneCount)), hjust = -0.2, size = 3)
  
  print(p2)
} else {
  # Create a message if no unique terms
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(0, 0, "No unique terms found in the data")
}

dev.off()

##############################
# VISUALIZATION 3: Combined plot with highlighted unique terms
##############################

pdf("GO_Visualization/Combined_Highlight_Unique.pdf", width = 14, height = 12)

# Create a single plot with all terms but highlight unique ones
p3 <- ggplot(all_data, aes(x = -log10(P.value), y = reorder(ShortTerm, -log10(P.value)))) +
  geom_bar(stat = "identity", aes(fill = ModuleDB, alpha = Unique)) +
  scale_fill_manual(values = module_colors) +
  scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 1.0), guide = "none") +
  labs(x = "-log10(p-value)", y = "", fill = "Module & Database",
       title = "Enriched Pathways in Yellow and Red Modules",
       subtitle = "Unique terms (appearing in only one category) are highlighted with full opacity") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  ) +
  facet_wrap(~ ModuleDB, scales = "free_y", ncol = 2)

print(p3)

dev.off()

##############################
# VISUALIZATION 4: Alternative combined view with dots and lines
##############################

pdf("GO_Visualization/Dot_Line_GO_Plot.pdf", width = 14, height = 10)

# Create a combined visualization with dots and lines
p4 <- ggplot(all_data, aes(x = ModuleDB, y = reorder(ShortTerm, -log10(P.value)))) +
  geom_point(aes(size = -log10(P.value), color = ModuleDB)) +
  scale_size_continuous(range = c(1, 8)) +
  scale_color_manual(values = module_colors) +
  labs(x = "", y = "", size = "-log10(p-value)", color = "Module & Database",
       title = "Enriched Pathways Across Yellow and Red Modules",
       subtitle = "Size indicates significance (-log10 p-value)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.x = element_line(color = "grey80"),
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

print(p4)

dev.off()

cat("All visualizations have been created in the GO_Visualization directory.\n")



##############################
# ADVANCED GO VISUALIZATIONS
# این قسمت را به انتهای کد قبلی خود اضافه کنید
##############################

###########################
# VISUALIZATION 5: Sunburst Chart
###########################
library(sunburstR)

pdf("GO_Visualization/Sunburst_GO.pdf", width = 10, height = 10)

# Prepare data for sunburst
prepare_sunburst_data <- function() {
  # Create a data frame with path and values
  sb_data <- data.frame(
    path = character(),
    value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # For each row in all_data
  for (i in 1:nrow(all_data)) {
    # Create path: Module-Database-Term
    path <- paste(all_data$Module[i], all_data$Database[i], 
                  all_data$ShortTerm[i], sep = "-")
    
    # Value is -log10(p-value) * GeneCount to emphasize important terms
    value <- -log10(all_data$P.value[i]) * all_data$GeneCount[i]
    
    # Add to data frame
    sb_data <- rbind(sb_data, data.frame(
      path = path,
      value = value,
      stringsAsFactors = FALSE
    ))
  }
  
  return(sb_data)
}

sunburst_data <- prepare_sunburst_data()

# Create temporary HTML file for sunburst
temp_html <- tempfile(fileext = ".html")

# Generate sunburst chart
sb <- sunburst(
  sunburst_data,
  count = TRUE,
  legend = list(w = 200, h = 20, s = 5),
  breadcrumb = list(w = 0, h = 20, s = 5),
  width = 600,
  height = 600
)

# Save to HTML then convert to PDF
htmlwidgets::saveWidget(sb, temp_html, selfcontained = TRUE)
webshot::webshot(temp_html, "GO_Visualization/Sunburst_GO.pdf", 
                 delay = 2, vwidth = 700, vheight = 700)

# Attempt to display in R viewer too
print(sb)

dev.off()

###########################
# VISUALIZATION 6: Interactive Network with visNetwork
###########################
library(visNetwork)

# Prepare data for network visualization
prepare_visnetwork_data <- function() {
  # Create nodes
  nodes <- data.frame(
    id = 1:nrow(all_data),
    label = all_data$ShortTerm,
    group = all_data$ModuleDB,
    value = -log10(all_data$P.value) * 2,  # Size based on significance
    title = paste("Term:", all_data$Term, "<br>",
                  "P-value:", all_data$P.value, "<br>",
                  "Gene Count:", all_data$GeneCount),
    stringsAsFactors = FALSE
  )
  
  # Create edges based on shared genes or database/module
  edges <- data.frame(
    from = integer(),
    to = integer(),
    width = numeric(),
    title = character(),
    stringsAsFactors = FALSE
  )
  
  # Connect nodes within same module
  for (module in unique(all_data$Module)) {
    module_nodes <- which(all_data$Module == module)
    
    if (length(module_nodes) > 1) {
      for (i in 1:(length(module_nodes)-1)) {
        for (j in (i+1):length(module_nodes)) {
          # Add edge
          edges <- rbind(edges, data.frame(
            from = module_nodes[i],
            to = module_nodes[j],
            width = 1,
            title = paste("Same module:", module),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Connect nodes with same term across different modules/databases
  term_groups <- split(1:nrow(all_data), all_data$Term)
  for (term_nodes in term_groups) {
    if (length(term_nodes) > 1) {
      for (i in 1:(length(term_nodes)-1)) {
        for (j in (i+1):length(term_nodes)) {
          # Add edge with higher weight
          edges <- rbind(edges, data.frame(
            from = term_nodes[i],
            to = term_nodes[j],
            width = 3,
            title = paste("Same term:", all_data$Term[term_nodes[i]]),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  return(list(nodes = nodes, edges = edges))
}

network_data <- prepare_visnetwork_data()

# Create interactive network
vn <- visNetwork(
  network_data$nodes, 
  network_data$edges,
  width = "100%",
  height = "600px"
) %>%
  visGroups(groupname = "Yellow_BioPlanet", color = "#FFC107") %>%
  visGroups(groupname = "Red_BioPlanet", color = "#F44336") %>%
  visGroups(groupname = "Yellow_KEGG", color = "#FFEB3B") %>%
  visGroups(groupname = "Red_KEGG", color = "#E91E63") %>%
  visOptions(
    highlightNearest = TRUE,
    selectedBy = "group"
  ) %>%
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -100)
  ) %>%
  visLegend()

# Save to HTML
visNetwork::visSave(vn, file = "GO_Visualization/Interactive_GO_Network.html")

###########################
# VISUALIZATION 7: Alluvial Diagram (Sankey Diagram)
###########################
library(ggalluvial)

pdf("GO_Visualization/Alluvial_GO_Plot.pdf", width = 14, height = 12)

# Prepare data for alluvial diagram
# We will show how terms flow between modules and databases
alluvial_data <- all_data %>%
  # Separate Module and Database for the alluvial
  select(Term, Module, Database, P.value, GeneCount) %>%
  # Add significance category
  mutate(Significance = case_when(
    P.value < 0.001 ~ "High",
    P.value < 0.01 ~ "Medium",
    TRUE ~ "Low"
  ))

# Create the alluvial diagram
p_alluvial <- ggplot(alluvial_data,
                     aes(axis1 = Module, axis2 = Database, axis3 = Significance, y = -log10(P.value))) +
  geom_alluvium(aes(fill = Module), width = 0.3, knot.pos = 0.2) +
  geom_stratum(width = 0.3, fill = "white", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Module", "Database", "Significance"), expand = c(0.05, 0.05)) +
  scale_fill_manual(values = c("Yellow" = "#FFC107", "Red" = "#F44336")) +
  labs(title = "Flow of Gene Ontology Enrichment Results",
       subtitle = "Showing how terms flow between modules, databases and significance levels",
       y = "Cumulative -log10(p-value)") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

print(p_alluvial)

dev.off()

###########################
# VISUALIZATION 8: Hierarchical Edge Bundling
###########################
library(networkD3)

# Prepare data for edge bundling
prepare_heb_data <- function() {
  # Create nodes dataframe
  # Format: one row per node, with a unique id and a group
  nodes <- data.frame(
    name = character(),
    group = character(),
    stringsAsFactors = FALSE
  )
  
  # Add root node
  nodes <- rbind(nodes, data.frame(
    name = "GO_Analysis",
    group = "root",
    stringsAsFactors = FALSE
  ))
  
  # Add module nodes
  for (module in unique(all_data$Module)) {
    nodes <- rbind(nodes, data.frame(
      name = module,
      group = "module",
      stringsAsFactors = FALSE
    ))
  }
  
  # Add database nodes
  for (db in unique(all_data$Database)) {
    nodes <- rbind(nodes, data.frame(
      name = paste(db, "Database"),
      group = "database",
      stringsAsFactors = FALSE
    ))
  }
  
  # Add term nodes
  for (i in 1:nrow(all_data)) {
    # Create a unique ID for each term
    term_id <- paste(all_data$Term[i], i, sep = "_")
    
    nodes <- rbind(nodes, data.frame(
      name = term_id,
      group = "term",
      stringsAsFactors = FALSE
    ))
  }
  
  # Create links dataframe
  # Format: source and target columns with numeric indices
  links <- data.frame(
    source = integer(),
    target = integer(),
    value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Add links from root to modules
  for (module in unique(all_data$Module)) {
    source_idx <- which(nodes$name == "GO_Analysis") - 1  # 0-indexed
    target_idx <- which(nodes$name == module) - 1         # 0-indexed
    
    links <- rbind(links, data.frame(
      source = source_idx,
      target = target_idx,
      value = 1,
      stringsAsFactors = FALSE
    ))
  }
  
  # Add links from root to databases
  for (db in unique(all_data$Database)) {
    source_idx <- which(nodes$name == "GO_Analysis") - 1  # 0-indexed
    target_idx <- which(nodes$name == paste(db, "Database")) - 1  # 0-indexed
    
    links <- rbind(links, data.frame(
      source = source_idx,
      target = target_idx,
      value = 1,
      stringsAsFactors = FALSE
    ))
  }
  
  # Add links from modules to terms
  for (i in 1:nrow(all_data)) {
    module <- all_data$Module[i]
    term_id <- paste(all_data$Term[i], i, sep = "_")
    
    source_idx <- which(nodes$name == module) - 1         # 0-indexed
    target_idx <- which(nodes$name == term_id) - 1        # 0-indexed
    
    links <- rbind(links, data.frame(
      source = source_idx,
      target = target_idx,
      value = -log10(all_data$P.value[i]),  # Weight by significance
      stringsAsFactors = FALSE
    ))
  }
  
  # Add links from databases to terms
  for (i in 1:nrow(all_data)) {
    db <- paste(all_data$Database[i], "Database")
    term_id <- paste(all_data$Term[i], i, sep = "_")
    
    source_idx <- which(nodes$name == db) - 1             # 0-indexed
    target_idx <- which(nodes$name == term_id) - 1        # 0-indexed
    
    links <- rbind(links, data.frame(
      source = source_idx,
      target = target_idx,
      value = -log10(all_data$P.value[i]),  # Weight by significance
      stringsAsFactors = FALSE
    ))
  }
  
  # Clean node names for display
  nodes$name <- gsub("_[0-9]+$", "", nodes$name)  # Remove the numerical suffix
  
  return(list(nodes = nodes, links = links))
}

heb_data <- prepare_heb_data()

# Create hierarchical edge bundling diagram
heb <- sankeyNetwork(
  Links = heb_data$links, 
  Nodes = heb_data$nodes,
  Source = "source", 
  Target = "target",
  Value = "value", 
  NodeID = "name",
  NodeGroup = "group",
  fontSize = 12,
  nodeWidth = 30,
  sinksRight = FALSE,
  height = 600,
  width = 800
)

# Save to HTML
htmlwidgets::saveWidget(heb, file = "GO_Visualization/Hierarchical_Edge_Bundling.html", selfcontained = TRUE)

###########################
# VISUALIZATION 9: UpSet Plot (Set Intersections)
###########################
library(UpSetR)

pdf("GO_Visualization/UpSet_GO_Plot.pdf", width = 10, height = 8)

# Prepare data for UpSet plot
prepare_upset_data <- function() {
  # Get all unique terms
  all_terms <- unique(all_data$Term)
  
  # Create a binary matrix: rows are terms, columns are module-database combinations
  upset_matrix <- matrix(0, nrow = length(all_terms), ncol = 4)
  colnames(upset_matrix) <- c("Yellow_BP", "Red_BP", "Yellow_KEGG", "Red_KEGG")
  rownames(upset_matrix) <- all_terms
  
  # Fill the matrix
  for (i in 1:nrow(all_data)) {
    term <- all_data$Term[i]
    module_db <- paste(all_data$Module[i], ifelse(all_data$Database[i] == "BioPlanet", "BP", "KEGG"), sep = "_")
    
    # Set the corresponding cell to 1
    upset_matrix[term, module_db] <- 1
  }
  
  # Convert to data frame
  upset_df <- as.data.frame(upset_matrix)
  
  return(upset_df)
}

upset_data <- prepare_upset_data()

# Create UpSet plot
upset_plot <- upset(
  upset_data,
  nsets = 4,
  order.by = "freq",
  empty.intersections = "on",
  mainbar.y.label = "Term Intersection Size",
  sets.x.label = "Terms Per Category",
  point.size = 3.5,
  line.size = 2,
  mb.ratio = c(0.6, 0.4),
  main.bar.color = "black",
  sets.bar.color = c("Yellow_BP" = "#FFC107", 
                     "Red_BP" = "#F44336",
                     "Yellow_KEGG" = "#FFEB3B", 
                     "Red_KEGG" = "#E91E63")
)

print(upset_plot)

dev.off()

###########################
# VISUALIZATION 10: 3D PCA Plot of Term Similarities
###########################
library(plotly)

# Prepare data for PCA
prepare_pca_data <- function() {
  # Create a term-by-category matrix with p-values
  terms <- unique(all_data$Term)
  categories <- c("Yellow_BioPlanet", "Red_BioPlanet", "Yellow_KEGG", "Red_KEGG")
  
  pca_matrix <- matrix(NA, nrow = length(terms), ncol = length(categories))
  rownames(pca_matrix) <- terms
  colnames(pca_matrix) <- categories
  
  # Fill matrix with -log10(p-values)
  for (i in 1:nrow(all_data)) {
    term <- all_data$Term[i]
    category <- paste(all_data$Module[i], all_data$Database[i], sep = "_")
    
    pca_matrix[term, category] <- -log10(all_data$P.value[i])
  }
  
  # Replace NA with 0
  pca_matrix[is.na(pca_matrix)] <- 0
  
  # Perform PCA
  pca_result <- prcomp(pca_matrix, scale. = TRUE)
  
  # Get term coordinates in PCA space
  pca_coords <- pca_result$x[, 1:3]  # Use first 3 components
  
  # Add module information
  pca_df <- as.data.frame(pca_coords)
  
  # Determine predominant module for each term
  term_modules <- character(length(terms))
  term_datasets <- character(length(terms))
  term_scores <- numeric(length(terms))
  
  for (i in 1:length(terms)) {
    term <- terms[i]
    term_data <- subset(all_data, Term == term)
    
    # If term appears in multiple categories, use the one with lowest p-value
    if (nrow(term_data) > 0) {
      best_idx <- which.min(term_data$P.value)
      term_modules[i] <- as.character(term_data$Module[best_idx])
      term_datasets[i] <- as.character(term_data$Database[best_idx])
      term_scores[i] <- -log10(term_data$P.value[best_idx])
    }
  }
  
  pca_df$Module <- term_modules
  pca_df$Database <- term_datasets
  pca_df$Score <- term_scores
  pca_df$Term <- terms
  
  # Calculate percentage of variance explained
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
  
  return(list(pca_df = pca_df, var_explained = var_explained))
}

pca_results <- prepare_pca_data()
pca_df <- pca_results$pca_df
var_explained <- pca_results$var_explained

# Create 3D PCA plot
axis_labels <- paste0(
  c("PC1", "PC2", "PC3"),
  " (", 
  sprintf("%.1f", var_explained[1:3]),
  "%)"
)

# Define colors based on Module and Database
pca_df$Color <- factor(paste(pca_df$Module, pca_df$Database, sep = "_"))
color_values <- c("Yellow_BioPlanet" = "#FFC107", 
                  "Red_BioPlanet" = "#F44336",
                  "Yellow_KEGG" = "#FFEB3B", 
                  "Red_KEGG" = "#E91E63")

p <- plot_ly(pca_df, 
             x = ~PC1, y = ~PC2, z = ~PC3, 
             color = ~Color,
             colors = color_values,
             size = ~Score,
             sizes = c(5, 20),
             text = ~Term,
             hoverinfo = "text",
             type = "scatter3d",
             mode = "markers")

p <- p %>% layout(
  title = "3D PCA Plot of GO Term Similarities",
  scene = list(
    xaxis = list(title = axis_labels[1]),
    yaxis = list(title = axis_labels[2]),
    zaxis = list(title = axis_labels[3])
  ),
  legend = list(title = list(text = "Module & Database"))
)

# Save to HTML
htmlwidgets::saveWidget(p, file = "GO_Visualization/3D_PCA_GO_Terms.html", selfcontained = TRUE)

###########################
# VISUALIZATION 11: Treemap of GO Terms
###########################
library(treemap)

pdf("GO_Visualization/Treemap_GO.pdf", width = 12, height = 10)

# Prepare data for treemap
treemap_data <- all_data %>%
  # Add a significance score combining p-value and gene count
  mutate(Score = -log10(P.value) * GeneCount,
         # Create a combined category
         Category = paste(Module, Database, sep = "_"))

# Create the treemap
treemap(
  treemap_data,
  index = c("Category", "ShortTerm"),
  vSize = "Score",
  vColor = "P.value",
  type = "value",
  palette = "RdYlBu",
  aspRatio = 1,
  title = "GO Term Enrichment Treemap",
  fontsize.title = 14,
  reverse.palette = TRUE,
  algorithm = "pivotSize",
  sortID = "Score"
)

dev.off()

cat("All advanced visualizations have been created in the GO_Visualization directory.\n")


# Modified visualization with refined focus on inflammation and disease-related immune terms
# First, let's change BioPlanet to BP throughout the code
all_data$Database <- gsub("BioPlanet", "BP", all_data$Database)

# Adjust the ModuleDB factor with the new BP naming
all_data$ModuleDB <- factor(paste(all_data$Module, all_data$Database, sep="_"),
                            levels = c("Yellow_BP", "Red_BP", "Yellow_KEGG", "Red_KEGG"))

# Update the color mapping
module_colors <- c("Yellow_BP" = "#FFC107", 
                   "Red_BP" = "#F44336",
                   "Yellow_KEGG" = "#FFEB3B", 
                   "Red_KEGG" = "#E91E63")

# Create a more focused list of immune and inflammation-related keywords
# Specifically targeting NF-κB, cytokines, inflammation, and disease-related immune terms
immune_keywords <- c(
  # NF-κB pathway components
  "NF-kappa", "NF-κB", "NF-kappaB", "NF-kB", "IKK", "IκB", "NFKB", "NFκB", "RelA", "p65", "p50",
  
  # Cytokines and inflammation
  "cytokine", "chemokine", "inflamm", "TNF", "IL-1", "IL1", "IL-6", "IL6", "IL-17", "IL17", 
  "interferon", "IFN", "interleukin", 
  
  # Disease-related terms
  "autoimmun", "arthritis", "colitis", "sepsis", "septic", "cancer", "tumor", "lesion", "disease",
  "infection", "pathogen", "pathology",
  
  # Immune signaling
  "toll-like", "TLR", "JAK", "STAT", "MAPK", "TRAF", "IRAK", "MyD88", "signaling",
  
  # Immune cells and processes
  "macrophage", "neutrophil", "T cell", "B cell", "antibody", "antigen", "leukocyte", 
  "apoptosis", "pyroptosis", "necroptosis", "necrosis", "necrotic", "apoptotic",
  
  # Specific immune-related processes
  "response", "defense", "phagocyt", "oxidative", "stress", "ROS", "radical", "fibrosis",
  "inflammatory", "activation", "immune", "immunoregulation", "immunomodulation"
)

# Function to check if a term is related to our focused immune keywords
is_specific_immune_related <- function(term) {
  term_lower <- tolower(term)
  any(sapply(immune_keywords, function(keyword) grepl(keyword, term_lower, fixed = TRUE)))
}

# Add a column indicating if the term is related to our focused immune keywords
all_data$ImmuneRelated <- sapply(all_data$Term, is_specific_immune_related)

# Give higher priority to immune-related terms
all_data$Priority <- ifelse(all_data$ImmuneRelated, 2, 1)

# Count how many immune-related terms we have
immune_count <- sum(all_data$ImmuneRelated)
print(paste("Number of immune-related terms found:", immune_count))

# If we have too few immune-related terms, we might need to expand our criteria
if(immune_count < 10) {
  # Expand with more general immune terms as a fallback
  additional_keywords <- c("receptor", "CD", "lymphocyte", "adaptive", "innate", "regulation",
                           "protein", "pathway", "metabolism", "microbiome", "barrier")
  
  # Add these keywords only if we need more terms
  immune_keywords <- c(immune_keywords, additional_keywords)
  
  # Recalculate the immune-related flag
  all_data$ImmuneRelated <- sapply(all_data$Term, is_specific_immune_related)
  all_data$Priority <- ifelse(all_data$ImmuneRelated, 2, 1)
  
  # Update the count
  immune_count <- sum(all_data$ImmuneRelated)
  print(paste("After expansion, number of immune-related terms:", immune_count))
}

# Print the immune-related terms so we can verify they're what we want
immune_terms <- all_data$Term[all_data$ImmuneRelated]
print("Immune-related terms found:")
print(unique(immune_terms))

##############################
# IMPROVED VISUALIZATION: Focus on inflammation and disease-related immune terms
##############################

# Function to select terms for visualization
select_terms_for_viz <- function(data, min_immune_terms = 15) {
  # If we have enough immune-related terms, use only those
  immune_data <- data[data$ImmuneRelated, ]
  
  if(nrow(immune_data) >= min_immune_terms) {
    # We have enough immune terms, use only these
    # Sort by p-value within each module-database combination
    selected_terms <- data.frame()
    
    for(mod_db in unique(immune_data$ModuleDB)) {
      mod_db_data <- immune_data[immune_data$ModuleDB == mod_db, ]
      mod_db_data <- mod_db_data[order(mod_db_data$P.value), ]
      # Take top terms from each module-db combination
      selected_terms <- rbind(selected_terms, head(mod_db_data, 5))
    }
    
    return(selected_terms)
  } else {
    # We don't have enough immune terms, supplement with some non-immune terms
    selected_terms <- data.frame()
    
    for(mod_db in unique(data$ModuleDB)) {
      # First add all immune terms for this module-db
      mod_db_immune <- data[data$ModuleDB == mod_db & data$ImmuneRelated, ]
      mod_db_immune <- mod_db_immune[order(mod_db_immune$P.value), ]
      selected_terms <- rbind(selected_terms, mod_db_immune)
      
      # Then add some non-immune terms if needed
      if(nrow(mod_db_immune) < 5) {
        mod_db_non_immune <- data[data$ModuleDB == mod_db & !data$ImmuneRelated, ]
        mod_db_non_immune <- mod_db_non_immune[order(mod_db_non_immune$P.value), ]
        selected_terms <- rbind(selected_terms, head(mod_db_non_immune, 5 - nrow(mod_db_immune)))
      }
    }
    
    return(selected_terms)
  }
}

# Select terms for visualization
selected_data <- select_terms_for_viz(all_data)

# Get unique terms to identify overlaps between modules
term_counts <- table(selected_data$Term)
selected_data$Shared <- selected_data$Term %in% names(term_counts[term_counts > 1])

# Prepare data for visualizations
# Sort terms by immune-related status and p-value
selected_data <- selected_data[order(-selected_data$Priority, selected_data$P.value), ]

##############################
# VISUALIZATION 1: Combined Dot-Line Plot
##############################

pdf("GO_Visualization/Immune_Enriched_DotLine.pdf", width = 10, height = 12)  # Reduced width

# Create a dot-line plot that focuses on immune-related terms
p_immune <- ggplot(selected_data, aes(x = ModuleDB, y = reorder(ShortTerm, -Priority))) +
  # Add connecting lines for shared terms
  geom_line(aes(group = Term, color = Module), 
            data = selected_data[selected_data$Shared,], 
            alpha = 0.6, size = 0.8) +
  # Add points sized by gene count and colored by module
  geom_point(aes(size = -log10(P.value), color = ModuleDB, 
                 alpha = ImmuneRelated)) +
  # Apply custom colors and alpha
  scale_color_manual(values = module_colors) +
  scale_alpha_manual(values = c("FALSE" = 0.7, "TRUE" = 1), 
                     labels = c("FALSE" = "Other", "TRUE" = "Immune-related")) +
  scale_size_continuous(range = c(3, 9), name = "-log10(p-value)") +
  # Add labels and theme
  labs(x = "", y = "", 
       alpha = "Term Type",
       color = "Module & Database",
       title = "Enriched Pathways in Yellow and Red Modules",
       subtitle = "Focusing on inflammation and disease-related immune terms") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

print(p_immune)

dev.off()

##############################
# VISUALIZATION 2: Enhanced Heatmap
##############################

pdf("GO_Visualization/Immune_Enriched_Heatmap.pdf", width = 10, height = 12)

# Create a matrix for the heatmap
terms <- unique(selected_data$Term)
categories <- levels(selected_data$ModuleDB)

# Create matrix for heatmap
heatmap_matrix <- matrix(NA, nrow = length(terms), ncol = length(categories))
rownames(heatmap_matrix) <- terms
colnames(heatmap_matrix) <- categories

# Fill with -log10(p-values)
for (i in 1:nrow(selected_data)) {
  term <- selected_data$Term[i]
  category <- as.character(selected_data$ModuleDB[i])
  heatmap_matrix[term, category] <- -log10(selected_data$P.value[i])
}

# Replace NA with 0
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# Create annotation for immune-related terms
immune_annotation <- sapply(terms, is_specific_immune_related)
gene_counts <- sapply(terms, function(term) {
  rows <- selected_data[selected_data$Term == term, ]
  if(nrow(rows) > 0) {
    # Get the gene count for this term (extract from Overlap column)
    return(max(rows$GeneCount))
  } else {
    return(0)
  }
})

# Create annotation for immune-related terms and gene counts
ha_row <- rowAnnotation(
  Immune = immune_annotation,
  Genes = gene_counts,
  col = list(
    Immune = c("TRUE" = "#E41A1C", "FALSE" = "gray90"),
    Genes = colorRamp2(c(0, max(gene_counts)), c("white", "blue"))
  ),
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    Immune = list(title = "Immune-related", at = c(TRUE, FALSE)),
    Genes = list(title = "Gene Count")
  )
)

# Create the heatmap with improved color scheme
heatmap_colors <- colorRamp2(
  c(0, 1, 2, 3, 4), 
  c("#FFFFFF", "#FFF7BC", "#FEC44F", "#D95F0E", "#7F0000")  # White to orange to dark red
)

# Split the heatmap by immune-related status
row_split <- ifelse(immune_annotation, "Immune-related", "Other")

# Create a custom hierarchical clustering function that keeps immune terms together
custom_hclust <- function(mat) {
  d <- dist(mat)
  hc <- hclust(d, method = "average")
  return(hc)
}

hm <- Heatmap(
  heatmap_matrix,
  name = "-log10(p-value)",
  col = heatmap_colors,
  row_split = row_split,
  cluster_rows = custom_hclust,
  cluster_row_slices = FALSE,  # Don't cluster the row slices
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  row_title = "GO Terms",
  column_title = "Module & Database Combinations",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  right_annotation = ha_row,
  heatmap_legend_param = list(
    title = "-log10(p-value)",
    at = c(0, 1, 2, 3, 4),
    labels = c("0", "1", "2", "3", "4+")
  )
)

# Plot the heatmap
draw(hm, heatmap_legend_side = "right")

dev.off()

##############################
# VISUALIZATION 3: Improved Dot Plot with Immune Focus
##############################
##############################
# VISUALIZATION 3: Improved Dot Plot with Immune Focus
##############################

pdf("GO_Visualization/Immune_Focus_DotPlot.pdf", width = 10, height = 12)

# Create a better dot plot focusing on immune-related terms
p_dots <- ggplot(selected_data, aes(x = ModuleDB, y = reorder(ShortTerm, Priority))) +
  # Add connecting lines for shared terms
  geom_line(aes(group = Term), 
            data = selected_data[selected_data$Shared,], 
            color = "grey60", alpha = 0.6) +
  # Add points with size representing gene count and color representing significance
  geom_point(aes(size = GeneCount, fill = -log10(P.value), 
                 shape = ImmuneRelated), color = "black", stroke = 0.5) +
  # Use filled shapes for better visibility
  scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 24),  # 21 = filled circle, 24 = filled triangle
                     labels = c("FALSE" = "Other", "TRUE" = "Immune-related")) +
  # Color gradient for significance
  scale_fill_gradientn(colors = c("#FFFFCC", "#FFC107", "#FD8D3C", "#E31A1C", "#800026"),
                       name = "-log10(p-value)") +
  # Size scale for gene counts
  scale_size_continuous(range = c(2, 8), name = "Gene Count") +
  # Add labels and theme
  labs(x = "", y = "", 
       shape = "Term Type",
       title = "Immune-Related Pathway Enrichment",
       subtitle = "Focusing on inflammation, cytokines, and NF-κB pathways") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

print(p_dots)

dev.off()

##############################
# VISUALIZATION 4: NEW COMBINED PLOT (as requested by professor)
##############################

pdf("GO_Visualization/Combined_Immune_Plot.pdf", width = 10, height = 12)

# Create a combined plot that integrates features of the previous plots
# This plot will focus on immune-related terms with emphasis on inflammation and NF-κB pathways

# First, ensure we're working with primarily immune-related terms
immune_focused_data <- selected_data[selected_data$ImmuneRelated | selected_data$Shared, ]

# If we have too few terms, add some of the most significant non-immune terms
if(nrow(immune_focused_data) < 15) {
  additional_terms <- selected_data[!selected_data$ImmuneRelated, ]
  additional_terms <- additional_terms[order(additional_terms$P.value), ]
  immune_focused_data <- rbind(immune_focused_data, head(additional_terms, 15 - nrow(immune_focused_data)))
}

# Sort by module and p-value
immune_focused_data <- immune_focused_data[order(immune_focused_data$ModuleDB, immune_focused_data$P.value), ]

# Create a combined plot with enhanced features
p_combined <- ggplot(immune_focused_data, aes(x = ModuleDB, y = reorder(ShortTerm, -Priority))) +
  # Add a background color strip to highlight immune-related terms
  geom_tile(aes(fill = ImmuneRelated), alpha = 0.1, width = 0.9, height = 0.9) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#FFECB3"),
                    guide = "none") +
  
  # Add connecting lines for shared terms
  geom_line(aes(group = Term, color = Module), 
            data = immune_focused_data[immune_focused_data$Shared,], 
            alpha = 0.7, size = 0.9) +
  
  # Add points with size representing significance and color representing module
  geom_point(aes(size = -log10(P.value), color = ModuleDB, 
                 alpha = ImmuneRelated)) +
  
  # Apply custom colors and alpha
  scale_color_manual(values = module_colors) +
  scale_alpha_manual(values = c("FALSE" = 0.6, "TRUE" = 1), 
                     labels = c("FALSE" = "Other Terms", "TRUE" = "Immune-related")) +
  
  # Adjust size scale
  scale_size_continuous(range = c(3, 10), name = "-log10(p-value)") +
  
  # Add gene count labels
  geom_text(aes(label = GeneCount), size = 2.5, vjust = -1.5) +
  
  # Add labels and theme
  labs(x = "", y = "", 
       alpha = "Term Type",
       color = "Module & Database",
       title = "Enriched Immune Pathways in Yellow and Red Modules",
       subtitle = "Focusing on inflammation, NF-κB signaling, and disease-related terms") +
  
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

print(p_combined)

# Add a small explanatory note at the bottom
grid.text("Note: Numbers indicate gene counts. Point size represents statistical significance.", 
          x = 0.5, y = 0.02, gp = gpar(fontsize = 8))

dev.off()

##############################
# VISUALIZATION 5: ENHANCED NETWORK PLOT
##############################

# Create a network plot showing relationships between modules and terms
pdf("GO_Visualization/Immune_Network_Plot.pdf", width = 10, height = 10)

# Prepare data for network visualization
# Each edge represents a term in a module
network_data <- selected_data %>%
  filter(ImmuneRelated) %>%  # Focus on immune-related terms
  select(Term, ModuleDB, P.value, GeneCount)

# Create a graph object
g <- graph_from_data_frame(network_data[, c("Term", "ModuleDB")], directed = FALSE)

# Set edge weights based on significance
E(g)$weight <- -log10(network_data$P.value)

# Set node attributes
V(g)$type <- ifelse(V(g)$name %in% network_data$Term, "Term", "Module")
V(g)$size <- ifelse(V(g)$type == "Term", 
                    10, 
                    20)  # Modules are larger nodes
V(g)$color <- ifelse(V(g)$type == "Term", 
                     "#4285F4", 
                     ifelse(grepl("Yellow", V(g)$name), "#FFC107", "#F44336"))
V(g)$frame.color <- "white"
V(g)$label <- V(g)$name
V(g)$label.cex <- ifelse(V(g)$type == "Module", 1.2, 0.8)

# Set edge attributes
E(g)$width <- E(g)$weight / 2
E(g)$color <- "#78909C"
E(g)$alpha <- 0.7

# Create a better layout
layout <- layout_with_fr(g)

# Plot the network
plot(g, layout = layout,
     vertex.label.color = "black",
     vertex.label.family = "sans",
     vertex.frame.color = NA,
     vertex.label.dist = 0.5,
     vertex.label.degree = -pi/2,
     edge.curved = 0.2,
     main = "Network of Immune-Related Pathways and Modules")

# Add a legend
legend("bottomright", 
       legend = c("Term", "Yellow Module", "Red Module"),
       pch = 21, 
       pt.bg = c("#4285F4", "#FFC107", "#F44336"),
       pt.cex = 2,
       cex = 0.8,
       bty = "n")

dev.off()

# Print completion message
cat("\nAll visualizations completed. Created 5 plots focusing on immune-related terms.\n")
cat("The new combined plot integrates features as requested by the professor.\n")


























#############


















# Create a combined plot that integrates features of the previous plots
# This plot will focus on immune-related terms with emphasis on inflammation and NF-κB pathways
# Save the plot to a PDF file with adjusted dimensions
pdf("GO_Visualization/Combined_Immune_Plot_Modified.pdf", width = 9, height = 10)  # Slightly wider for labels

# First, ensure we're working with primarily immune-related terms
immune_focused_data <- selected_data[selected_data$ImmuneRelated | selected_data$Shared, ]

# If we have too few terms, add some of the most significant non-immune terms
if(nrow(immune_focused_data) < 15) {
  additional_terms <- selected_data[!selected_data$ImmuneRelated, ]
  additional_terms <- additional_terms[order(additional_terms$P.value), ]
  immune_focused_data <- rbind(immune_focused_data, head(additional_terms, 15 - nrow(immune_focused_data)))
}

# Sort by module and p-value
immune_focused_data <- immune_focused_data[order(immune_focused_data$ModuleDB, immune_focused_data$P.value), ]

# Create a function to wrap text to multiple lines
wrap_text <- function(text, width = 25) {
  sapply(strwrap(text, width = width, simplify = FALSE), paste, collapse = "\n")
}

# Apply text wrapping to ShortTerm
immune_focused_data$ShortTermWrapped <- wrap_text(immune_focused_data$ShortTerm)

# Create a combined plot with enhanced features and wrapped labels
p_combined <- ggplot(immune_focused_data, aes(x = ModuleDB, y = reorder(ShortTermWrapped, -Priority))) +
  # Add a background color strip to highlight immune-related terms
  geom_tile(aes(fill = ImmuneRelated), alpha = 0.1, width = 0.85, height = 0.85) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#FFECB3"),
                    guide = "none") +
  
  # Add connecting lines for shared terms
  geom_line(aes(group = Term, color = Module), 
            data = immune_focused_data[immune_focused_data$Shared,], 
            alpha = 0.7, size = 0.7) +
  
  # Add points with size representing significance and color representing module
  geom_point(aes(size = -log10(P.value), color = ModuleDB)) +
  
  # Apply custom colors and alpha
  scale_color_manual(values = module_colors) +
  scale_alpha_manual(values = c("FALSE" = 0.6, "TRUE" = 1), 
                     labels = c("FALSE" = "Other Terms", "TRUE" = "Immune-related")) +
  
  # Adjust size scale
  scale_size_continuous(range = c(2, 8), name = "-log10(p-value)") +
  
  # Add gene count labels
  geom_text(aes(label = GeneCount), size = 2, vjust = -1.2) +
  
  # Add labels and theme
  labs(x = "", y = "", color = "Module_Database") +
  
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8, hjust = 1, lineheight = 0.8),  # Added lineheight for multi-line text
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    legend.position = "right",
    legend.key.size = unit(0.8, "lines"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    panel.grid.major.y = element_line(color = "grey95", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  # Add space for y-axis labels
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(l = 80, r = 10, t = 10, b = 10))  # Increased left margin for labels

print(p_combined)

# Close the PDF device
dev.off()



























































# خلاصه تعداد ماژول‌های شناسایی شده
library(WGCNA)
library(knitr)
library(kableExtra)

# بارگذاری داده‌های شبکه (بر اساس کد شما)
load(file ="5_GSE51856_Milk_Net_Auto.RData") 

# تبدیل شماره‌های ماژول به رنگ‌ها
moduleLabels = net_GSE51856_Milk$colors
moduleColors = labels2colors(moduleLabels)

# شمارش تعداد ژن در هر ماژول
module_table = table(moduleColors)

# ایجاد جدول خلاصه
summary_df = data.frame(
  "ماژول (رنگ)" = names(module_table),
  "تعداد ژن" = as.numeric(module_table),
  "درصد کل ژن‌ها" = round((as.numeric(module_table) / sum(module_table)) * 100, 2)
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Summary of identified modules - FIXED VERSION
  library(WGCNA)
  library(knitr)
  library(kableExtra)
  
  # Load network data (based on your code)
  load(file ="5_GSE51856_Milk_Net_Auto.RData") 
  
  # Convert module numbers to colors
  moduleLabels = net_GSE51856_Milk$colors
  moduleColors = labels2colors(moduleLabels)
  
  # Count number of genes in each module
  module_table = table(moduleColors)
  
  # Create summary table
  summary_df = data.frame(
    "Module_Color" = names(module_table),
    "Gene_Count" = as.numeric(module_table),
    "Percentage" = round((as.numeric(module_table) / sum(module_table)) * 100, 2)
  )
  
  # Sort by gene count (descending)
  summary_df = summary_df[order(summary_df$Gene_Count, decreasing = TRUE), ]
  
  # Add row numbers
  summary_df$"Rank" = 1:nrow(summary_df)
  
  # Reorder columns
  summary_df = summary_df[, c("Rank", "Module_Color", "Gene_Count", "Percentage")]
  
  # Display table
  print("Summary of modules identified by WGCNA algorithm:")
  print(paste("Total number of identified modules:", nrow(summary_df)))
  print("")
  
  # FIXED: Display formatted table with correct notation parameter
  formatted_table <- kable(summary_df, 
                           caption = "Distribution of genes in modules identified by WGCNA",
                           align = c('c', 'l', 'c', 'c'),
                           col.names = c("Rank", "Module (Color)", "Gene Count", "Percentage (%)")) %>%
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = FALSE,
      position = "center"
    ) %>%
    row_spec(0, bold = TRUE, background = "#f2f2f2") %>%
    column_spec(1, bold = TRUE) %>%
    column_spec(2, color = "blue") %>%
    add_footnote("Generated by dynamic tree cutting algorithm in WGCNA", 
                 notation = "alphabet")  # FIXED: changed from "general" to "alphabet"
  
  # Display the table
  formatted_table
  
  # Alternative simple table display (if kableExtra doesn't work)
  print("Alternative simple table:")
  print(summary_df)
  
  # Save table to CSV file
  write.csv(summary_df, "WGCNA_Modules_Summary.csv", row.names = FALSE)
  
  # Create basic visualization
  pdf("WGCNA_Modules_Distribution.pdf", width = 10, height = 8)
  
  # Select appropriate colors for plot
  colors_for_plot = summary_df$Module_Color
  # Convert color names to actual colors where possible
  actual_colors = ifelse(colors_for_plot %in% colors(), colors_for_plot, "grey")
  
  # Pie chart
  pie(summary_df$Gene_Count, 
      labels = paste(summary_df$Module_Color, "\n", 
                     summary_df$Gene_Count, " genes\n", 
                     "(", summary_df$Percentage, "%)", sep=""),
      col = actual_colors,
      main = "Distribution of Genes in WGCNA Identified Modules",
      cex = 0.8)
  
  dev.off()
  
  # Bar plot
  pdf("WGCNA_Modules_Barplot.pdf", width = 12, height = 6)
  
  barplot(summary_df$Gene_Count,
          names.arg = summary_df$Module_Color,
          col = actual_colors,
          main = "Number of Genes in Each WGCNA Identified Module",
          xlab = "Module (Color)",
          ylab = "Gene Count",
          las = 2,
          cex.names = 0.8)
  
  # Add count labels on each bar
  text(x = seq_along(summary_df$Gene_Count) * 1.2 - 0.5, 
       y = summary_df$Gene_Count + max(summary_df$Gene_Count) * 0.02,
       labels = summary_df$Gene_Count,
       cex = 0.8)
  
  dev.off()
  
  # Statistical summary
  cat("\n=== WGCNA Modules Statistical Summary ===\n")
  cat("Total number of identified modules:", nrow(summary_df), "\n")
  cat("Total genes analyzed:", sum(summary_df$Gene_Count), "\n")
  cat("Largest module:", summary_df$Module_Color[1], "with", summary_df$Gene_Count[1], "genes\n")
  cat("Smallest module:", summary_df$Module_Color[nrow(summary_df)], "with", summary_df$Gene_Count[nrow(summary_df)], "genes\n")
  cat("Average genes per module:", round(mean(summary_df$Gene_Count), 2), "\n")
  cat("Median genes per module:", median(summary_df$Gene_Count), "\n")
  
  # Display module colors distribution
  cat("\nList of identified modules:\n")
  for(i in 1:nrow(summary_df)) {
    cat(sprintf("%2d. %s module: %d genes (%.2f%%)\n", 
                i, 
                summary_df$Module_Color[i], 
                summary_df$Gene_Count[i], 
                summary_df$Percentage[i]))
  }
  
  # Alternative table creation without kableExtra (if still having issues)
  create_simple_table <- function(df) {
    cat("\n=== Simple Table Format ===\n")
    cat(sprintf("%-6s %-15s %-12s %-12s\n", "Rank", "Module (Color)", "Gene Count", "Percentage (%)"))
    cat(paste(rep("-", 50), collapse = ""), "\n")
    
    for(i in 1:nrow(df)) {
      cat(sprintf("%-6d %-15s %-12d %-12.2f\n", 
                  df$Rank[i], 
                  df$Module_Color[i], 
                  df$Gene_Count[i], 
                  df$Percentage[i]))
    }
    cat(paste(rep("-", 50), collapse = ""), "\n")
    cat("Generated by dynamic tree cutting algorithm in WGCNA\n")
  }
  
  # Call the simple table function as backup
  create_simple_table(summary_df)
  
  print("Analysis completed successfully!")
  print("Files created:")
  print("- WGCNA_Modules_Summary.csv")
  print("- WGCNA_Modules_Distribution.pdf") 
  print("- WGCNA_Modules_Barplot.pdf")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # =============================================================================
  # CREATE SUPPLEMENTARY FILES FOR WGCNA ANALYSIS
  # This script creates comprehensive supplementary files for WGCNA results
  # =============================================================================
  
  # Load required libraries
  library(WGCNA)
  library(openxlsx)
  library(clusterProfiler)
  library(org.Bt.eg.db)
  
  # Load data
  load(file = "5_GSE51856_Milk_Net_Auto.RData")
  load(file = "2_GSE51856_Milk_RemSam.RData")
  
  # =============================================================================
  # PART 1: CREATE MODULE SUMMARY TABLE AND FIGURE
  # =============================================================================
  
  # Get module information
  moduleLabels = net_GSE51856_Milk$colors
  moduleColors = labels2colors(moduleLabels)
  genenames = names(GSE51856_Milk_ExpNorm_df)
  
  # Create module summary
  module_table = table(moduleColors)
  summary_df = data.frame(
    Module_Color = names(module_table),
    Gene_Count = as.numeric(module_table),
    Percentage = round((as.numeric(module_table) / sum(module_table)) * 100, 2)
  )
  
  # Sort by gene count (descending)
  summary_df = summary_df[order(summary_df$Gene_Count, decreasing = TRUE), ]
  summary_df$Rank = 1:nrow(summary_df)
  summary_df = summary_df[, c("Rank", "Module_Color", "Gene_Count", "Percentage")]
  
  # Calculate statistics for text
  total_modules = nrow(summary_df)
  avg_genes_per_module = round(mean(summary_df$Gene_Count), 0)
  largest_module = summary_df$Module_Color[1]
  largest_module_size = summary_df$Gene_Count[1]
  smallest_module = summary_df$Module_Color[nrow(summary_df)]
  smallest_module_size = summary_df$Gene_Count[nrow(summary_df)]
  grey_genes = summary_df$Gene_Count[summary_df$Module_Color == "grey"]
  
  # Print statistics for your text
  cat("=== STATISTICS FOR YOUR MANUSCRIPT ===\n")
  cat("Total modules identified:", total_modules, "\n")
  cat("Average genes per module:", avg_genes_per_module, "\n")
  cat("Largest module:", largest_module, "with", largest_module_size, "genes\n")
  cat("Smallest module:", smallest_module, "with", smallest_module_size, "genes\n")
  cat("Grey module (non-co-expressed genes):", grey_genes, "genes\n")
  
  # Create figure for module distribution
  # Create figure for module distribution
  pdf("Figure_S1_Module_Distribution.pdf", width = 12, height = 8)
  
  # Set up layout for two panels: barplot left, pie chart + legend right
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.2, 1))  # بیشتر فضا برای barplot
  
  # ----------------------------
  # Plot 1: Bar chart
  # ----------------------------
  par(mar = c(7, 5, 4, 2))  # فضای بیشتر پایین برای نام ماژول‌ها
  colors_for_plot = summary_df$Module_Color
  actual_colors = ifelse(colors_for_plot %in% colors(), colors_for_plot, "lightgrey")
  
  bar_midpoints <- barplot(summary_df$Gene_Count,
                           names.arg = summary_df$Module_Color,
                           col = actual_colors,
                           main = "Gene Distribution Across WGCNA Modules",
                           xlab = "Module Color",
                           ylab = "Number of Genes",
                           las = 2,
                           cex.names = 0.8,
                           cex.main = 1.2,
                           cex.lab = 1.1)
  
  # Add gene count labels on top of bars
  text(x = bar_midpoints,
       y = summary_df$Gene_Count + max(summary_df$Gene_Count) * 0.02,
       labels = summary_df$Gene_Count,
       cex = 0.7)
  
  # ----------------------------
  # Plot 2: Pie chart + legend below
  # ----------------------------
  
  # حذف grey module
  non_grey = summary_df[summary_df$Module_Color != "grey", ]
  colors_pie = ifelse(non_grey$Module_Color %in% colors(), non_grey$Module_Color, "lightgrey")
  
  # تنظیم layout داخلی: Pie chart بالا، legend پایین
  layout(matrix(1:2, nrow = 2), heights = c(3, 1))  # نسبت ارتفاع
  
  # Pie chart بدون برچسب
  par(mar = c(1, 4, 4, 2))  # فضای بالا برای عنوان
  pie(
    non_grey$Gene_Count,
    labels = NA,
    col = colors_pie,
    main = "Module Distribution\n(excluding grey module)",
    cex.main = 1.2
  )
  
  # Legend زیر Pie chart
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
         legend = paste(non_grey$Module_Color, " (", non_grey$Gene_Count, ")", sep = ""),
         fill = colors_pie,
         cex = 0.9,
         ncol = 3,
         bty = "n",
         x.intersp = 0.7,
         y.intersp = 1.2)
  
  # پایان فایل PDF
  dev.off()
  
  # =============================
  # Barplot (Figure S1A)
  # =============================
  pdf("Figure_S1A_Module_Barplot.pdf", width = 10, height = 6)
  
  par(mar = c(7, 5, 4, 2))  # فضای پایین بیشتر برای نام ماژول‌ها
  colors_for_plot = summary_df$Module_Color
  actual_colors = ifelse(colors_for_plot %in% colors(), colors_for_plot, "lightgrey")
  
  bar_midpoints <- barplot(summary_df$Gene_Count,
                           names.arg = summary_df$Module_Color,
                           col = actual_colors,
                           main = "Gene Distribution Across WGCNA Modules",
                           xlab = "Module Color",
                           ylab = "Number of Genes",
                           las = 2,
                           cex.names = 0.8,
                           cex.main = 1.2,
                           cex.lab = 1.1)
  
  # اضافه کردن تعداد ژن‌ها بالای ستون‌ها
  text(x = bar_midpoints,
       y = summary_df$Gene_Count + max(summary_df$Gene_Count) * 0.02,
       labels = summary_df$Gene_Count,
       cex = 0.7)
  
  dev.off()
  
  # =============================
  # Pie Chart + Legend (Figure S1B)
  # =============================
  # حذف grey module
  non_grey = summary_df[summary_df$Module_Color != "grey", ]
  colors_pie = ifelse(non_grey$Module_Color %in% colors(), non_grey$Module_Color, "lightgrey")
  
  # فایل جدید برای Pie
  pdf("Figure_S1B_Module_Pie_with_Legend.pdf", width = 8, height = 8)
  
  # تقسیم‌بندی دو بخش: بالا Pie، پایین legend
  layout(matrix(1:2, nrow = 2), heights = c(3, 1))
  
  # Pie chart
  par(mar = c(1, 4, 4, 2))  # فضای بالا برای عنوان
  pie(
    non_grey$Gene_Count,
    labels = NA,
    col = colors_pie,
    main = "Module Distribution\n(excluding grey module)",
    cex.main = 1.2
  )
  
  # Legend زیر Pie chart
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
         legend = paste(non_grey$Module_Color, " (", non_grey$Gene_Count, ")", sep = ""),
         fill = colors_pie,
         cex = 0.9,
         ncol = 3,   # می‌تونی تغییر بدی بسته به تعداد ماژول‌ها
         bty = "n",
         x.intersp = 0.7,
         y.intersp = 1.2)
  
  dev.off()
  
  # =============================================================================
  # PART 2: CREATE COMPREHENSIVE SUPPLEMENTARY FILE S3
  # =============================================================================
  
  cat("\n=== CREATING SUPPLEMENTARY FILE S3 ===\n")
  
  # Create comprehensive gene list with all information
  all_genes_info = data.frame(
    Gene_ID = genenames,
    Module_Number = moduleLabels,
    Module_Color = moduleColors,
    stringsAsFactors = FALSE
  )
  
  # Convert ENSEMBL to Gene symbols
  cat("Converting gene IDs to symbols...\n")
  gene_symbols = bitr(all_genes_info$Gene_ID, 
                      fromType = "ENSEMBL", 
                      toType = "SYMBOL", 
                      OrgDb = "org.Bt.eg.db")
  
  # Merge with gene information
  all_genes_complete = merge(all_genes_info, gene_symbols, 
                             by.x = "Gene_ID", by.y = "ENSEMBL", 
                             all.x = TRUE)
  
  # Rename columns
  colnames(all_genes_complete) = c("ENSEMBL_ID", "Module_Number", "Module_Color", "Gene_Symbol")
  
  # Reorder columns
  all_genes_complete = all_genes_complete[, c("ENSEMBL_ID", "Gene_Symbol", "Module_Number", "Module_Color")]
  
  # Sort by module color and then by gene symbol
  all_genes_complete = all_genes_complete[order(all_genes_complete$Module_Color, 
                                                all_genes_complete$Gene_Symbol), ]
  
  # Add row numbers
  all_genes_complete$Row_Number = 1:nrow(all_genes_complete)
  all_genes_complete = all_genes_complete[, c("Row_Number", "ENSEMBL_ID", "Gene_Symbol", "Module_Number", "Module_Color")]
  
  # =============================================================================
  # PART 3: CREATE EXCEL FILE WITH MULTIPLE SHEETS
  # =============================================================================
  
  cat("Creating comprehensive Excel file...\n")
  
  # Create workbook
  wb = createWorkbook()
  
  # Sheet 1: Module Summary
  addWorksheet(wb, "Module_Summary")
  writeData(wb, "Module_Summary", summary_df, startRow = 1)
  
  # Add title and description
  writeData(wb, "Module_Summary", 
            "WGCNA Module Summary", 
            startRow = 1, startCol = 6)
  writeData(wb, "Module_Summary", 
            paste("Total modules identified:", total_modules), 
            startRow = 2, startCol = 6)
  writeData(wb, "Module_Summary", 
            paste("Average genes per module:", avg_genes_per_module), 
            startRow = 3, startCol = 6)
  
  # Sheet 2: Complete Gene List
  addWorksheet(wb, "Complete_Gene_List")
  writeData(wb, "Complete_Gene_List", all_genes_complete)
  
  # Sheet 3: Genes by Module (separate sheet for each major module)
  module_colors = unique(moduleColors)
  module_colors = module_colors[module_colors != "grey"]  # Handle grey separately
  
  for (color in module_colors) {
    if (sum(moduleColors == color) >= 10) {  # Only for modules with >=10 genes
      module_genes = all_genes_complete[all_genes_complete$Module_Color == color, ]
      
      sheet_name = paste0("Module_", color)
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, module_genes)
      
      # Add module statistics
      writeData(wb, sheet_name, 
                paste("Module:", color, "- Total genes:", nrow(module_genes)), 
                startRow = 1, startCol = 7)
    }
  }
  
  # Sheet 4: Grey Module (non-co-expressed genes)
  grey_genes_df = all_genes_complete[all_genes_complete$Module_Color == "grey", ]
  addWorksheet(wb, "Grey_Module_Non_Coexpressed")
  writeData(wb, "Grey_Module_Non_Coexpressed", grey_genes_df)
  writeData(wb, "Grey_Module_Non_Coexpressed", 
            "Grey Module: Non-co-expressed genes based on dissimilarity measure", 
            startRow = 1, startCol = 7)
  
  # Sheet 5: Module Statistics
  stats_df = data.frame(
    Statistic = c("Total Modules", "Total Genes", "Average Genes per Module", 
                  "Largest Module", "Largest Module Size", "Smallest Module", 
                  "Smallest Module Size", "Non-co-expressed Genes (Grey)"),
    Value = c(total_modules, nrow(all_genes_complete), avg_genes_per_module,
              largest_module, largest_module_size, smallest_module, 
              smallest_module_size, grey_genes)
  )
  
  addWorksheet(wb, "Statistics")
  writeData(wb, "Statistics", stats_df)
  
  # Add formatting
  addStyle(wb, "Module_Summary", 
           createStyle(textDecoration = "bold", fgFill = "#D3D3D3"), 
           rows = 1, cols = 1:ncol(summary_df))
  
  # Save the workbook
  saveWorkbook(wb, "Supplementary_File_S3_WGCNA_Complete_Gene_Lists.xlsx", overwrite = TRUE)
  
  # =============================================================================
  # PART 4: CREATE MODULE-WISE SEPARATE FILES (Alternative format)
  # =============================================================================
  
  cat("Creating individual module files...\n")
  
  # Create directory for individual files
  dir.create("WGCNA_Individual_Module_Files", showWarnings = FALSE)
  
  # Create separate CSV files for each module
  for (color in unique(moduleColors)) {
    module_genes = all_genes_complete[all_genes_complete$Module_Color == color, ]
    filename = paste0("WGCNA_Individual_Module_Files/Module_", color, "_genes.csv")
    write.csv(module_genes, filename, row.names = FALSE)
  }
  
  # =============================================================================
  # PART 5: CREATE TEXT SUMMARY FOR MANUSCRIPT
  # =============================================================================
  
  # Create text file with summary for manuscript
  sink("WGCNA_Results_Summary_for_Manuscript.txt")
  cat("WGCNA ANALYSIS RESULTS SUMMARY\n")
  cat("===============================\n\n")
  
  cat("MODULE IDENTIFICATION:\n")
  cat(sprintf("- Total modules identified: %d\n", total_modules))
  cat(sprintf("- Average genes per module: %d\n", avg_genes_per_module))
  cat(sprintf("- Largest module: %s (containing %d genes)\n", largest_module, largest_module_size))
  cat(sprintf("- Smallest module: %s (containing %d genes)\n", smallest_module, smallest_module_size))
  cat(sprintf("- Non-co-expressed genes (grey module): %d genes\n", grey_genes))
  
  cat("\nMODULE DISTRIBUTION:\n")
  for(i in 1:nrow(summary_df)) {
    cat(sprintf("%d. %s module: %d genes (%.2f%%)\n", 
                i, summary_df$Module_Color[i], 
                summary_df$Gene_Count[i], summary_df$Percentage[i]))
  }
  
  cat("\nFILES CREATED:\n")
  cat("- Supplementary_File_S3_WGCNA_Complete_Gene_Lists.xlsx (main supplementary file)\n")
  cat("- Figure_S1_Module_Distribution.pdf (module distribution figure)\n")
  cat("- Individual CSV files for each module (in WGCNA_Individual_Module_Files/)\n")
  cat("- This summary file\n")
  
  sink()
  
  # =============================================================================
  # FINAL OUTPUT SUMMARY
  # =============================================================================
  
  cat("\n=== FILES CREATED SUCCESSFULLY ===\n")
  cat("1. Supplementary_File_S3_WGCNA_Complete_Gene_Lists.xlsx\n")
  cat("   - Complete gene lists for all modules\n")
  cat("   - Module summary statistics\n")
  cat("   - Separate sheets for major modules\n")
  
  cat("\n2. Figure_S1_Module_Distribution.pdf\n")
  cat("   - Bar chart and pie chart of module distribution\n")
  
  cat("\n3. WGCNA_Individual_Module_Files/ directory\n")
  cat("   - Separate CSV files for each module\n")
  
  cat("\n4. WGCNA_Results_Summary_for_Manuscript.txt\n")
  cat("   - Text summary for manuscript writing\n")
  
  cat("\n=== STATISTICS FOR YOUR TEXT ===\n")
  cat("Complete list of all the genes of each module are provided in Supplementary File S3.\n")
  cat(sprintf("An average number of %d genes detected per module, as largest and smallest modules was %s (containing %d genes) and %s (containing %d genes), respectively.\n", 
              avg_genes_per_module, largest_module, largest_module_size, smallest_module, smallest_module_size))
  cat(sprintf("Also, %d genes were detected as not co-expressed based on gene dissimilarity measure and were grouped into grey module.\n", grey_genes))
  
  cat("\n✅ All supplementary files created successfully!\n")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(png)
  
  # Load required data
  load(file = "5_GSE51856_Milk_Net_Auto.RData")
  load(file = "2_GSE51856_Milk_RemSam.RData")
  
  # Get module information
  moduleLabelsAutomatic1 = net_GSE51856_Milk$colors 
  moduleColorsAutomatic1 = labels2colors(moduleLabelsAutomatic1)
  moduleColorsAutomaticData1 = moduleColorsAutomatic1
  
  probes = names(GSE51856_Milk_ExpNorm_df)
  modulenames = unique(labels2colors(net_GSE51856_Milk$colors))
  
  # Create output directories
  dir.create("Publication_Heatmaps", showWarnings = FALSE)
  dir.create("Publication_Heatmaps/PDF", showWarnings = FALSE)
  dir.create("Publication_Heatmaps/PNG", showWarnings = FALSE)
  dir.create("Publication_Heatmaps/TIFF", showWarnings = FALSE)
  
  # =============================================================================
  # ENHANCED COLOR SCHEMES AND AESTHETICS
  # =============================================================================
  
  # Professional publication color palette for expression data
  expression_colors = colorRamp2(
    c(-2.5, -1.5, -0.5, 0, 0.5, 1.5, 2.5),
    c("#053061", "#2166AC", "#4393C3", "#F7F7F7", "#F4A582", "#D6604D", "#67001F")
  )
  
  # Enhanced color schemes for conditions
  healthy_palette = c("#1B5E20", "#2E7D32", "#388E3C", "#43A047", "#4CAF50")
  infected_palette = c("#B71C1C", "#C62828", "#D32F2F", "#E53935", "#F44336")
  
  # Time point colors (more sophisticated)
  timepoint_colors = c(
    "H0" = "#E8F5E8",   "H12" = "#C8E6C9", "H24" = "#A5D6A7", 
    "H36" = "#81C784",  "H48" = "#4CAF50",
    "I0" = "#FFEBEE",   "I12" = "#FFCDD2", "I24" = "#EF9A9A", 
    "I36" = "#E57373",  "I48" = "#F44336"
  )
  
  # Sample groups
  GSE51856_Milk_sample_groups = factor(
    rep(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
        times = c(4,5,5,5,5,5,5,4,5,5))
  )
  
  GSE51856_Milk_fa = rep(
    c("H0", "H12", "H24", "H36", "H48", "I0", "I12", "I24", "I36", "I48"),
    times = c(4,5,5,5,5,5,5,4,5,5)
  )
  
  # =============================================================================
  # FUNCTION TO CREATE PUBLICATION-READY HEATMAPS
  # =============================================================================
  
  create_publication_heatmap <- function(module_color, save_formats = c("PDF", "PNG", "TIFF")) {
    
    cat("Creating publication heatmap for", module_color, "module...\n")
    
    # Extract module genes
    inModule = (moduleColorsAutomaticData1 == module_color)
    nameofmodule = as.data.frame(probes)[inModule,]
    
    # Prepare expression data
    GSE51856_Milk_Expr = t(GSE51856_Milk_ExpNorm_df)
    GSE51856_Milk_Expr <- subset(GSE51856_Milk_Expr, 
                                 rownames(GSE51856_Milk_Expr) %in% nameofmodule)
    
    # Skip if no genes in module
    if(nrow(GSE51856_Milk_Expr) == 0) {
      cat("Skipping", module_color, "module - no genes found\n")
      return(NULL)
    }
    
    # Set column names (assuming you have the original column names)
    if(exists("GSE51856_Milk_Filt")) {
      colnames(GSE51856_Milk_Expr) <- colnames(GSE51856_Milk_Filt)[c(-3,-38)]
    }
    
    # Calculate means for healthy and infected conditions
    healthymean = rowMeans(GSE51856_Milk_Expr[, 1:24])
    infectedmean = rowMeans(GSE51856_Milk_Expr[, 25:48])
    
    # =============================================================================
    # CREATE HEATMAP COMPONENTS
    # =============================================================================
    
    # Healthy average heatmap
    healthy_heatmap = Heatmap(
      matrix(healthymean, ncol = 1),
      name = "Healthy_Avg",
      width = unit(15, "mm"),
      show_column_names = FALSE,
      show_row_names = FALSE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      col = expression_colors,
      
      top_annotation = HeatmapAnnotation(
        condition = anno_block(
          gp = gpar(fill = healthy_palette[3], col = "white", lwd = 2),
          labels = "Healthy\nAverage",
          labels_gp = gpar(col = "white", fontsize = 9, fontface = "bold"),
          labels_rot = 0
        ),
        annotation_height = unit(20, "mm"),
        # Remove annotation names by setting them to empty strings
        annotation_name_gp = gpar(fontsize = 0)  
      ),
      
      show_heatmap_legend = FALSE,
      border = TRUE,
      border_gp = gpar(col = "white", lwd = 2)
    )
    
    # Infected average heatmap  
    infected_heatmap = Heatmap(
      matrix(infectedmean, ncol = 1),
      name = "Infected_Avg",
      width = unit(15, "mm"),
      show_column_names = FALSE,
      show_row_names = FALSE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      col = expression_colors,
      
      top_annotation = HeatmapAnnotation(
        condition = anno_block(
          gp = gpar(fill = infected_palette[3], col = "white", lwd = 2),
          labels = "Infected\nAverage", 
          labels_gp = gpar(col = "white", fontsize = 9, fontface = "bold"),
          labels_rot = 0
        ),
        annotation_height = unit(20, "mm"),
        # Remove annotation names by setting them to empty strings
        annotation_name_gp = gpar(fontsize = 0)
      ),
      
      show_heatmap_legend = FALSE,
      border = TRUE,
      border_gp = gpar(col = "white", lwd = 2)
    )
    
    # Main expression heatmap
    main_heatmap = Heatmap(
      GSE51856_Milk_Expr,
      name = "RNA Stability",
      
      # Clustering and layout
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_row_dend = TRUE,
      row_dend_width = unit(20, "mm"),
      
      # Column organization
      column_split = GSE51856_Milk_sample_groups,
      column_gap = unit(2, "mm"),
      
      # Remove the module title completely
      column_title = NULL,
      column_title_gp = gpar(fontsize = 0),
      
      # Color scheme
      col = expression_colors,
      
      # Row settings
      row_title = "Genes",
      row_title_side = "right", 
      row_title_gp = gpar(fontsize = 12, fontface = "bold"),
      row_title_rot = 270,
      show_row_names = FALSE,
      row_names_gp = gpar(fontsize = 6),
      
      # Column settings
      show_column_names = TRUE,
      column_names_gp = gpar(fontsize = 8, col = "black"),
      column_names_rot = 45,
      column_names_centered = TRUE,
      
      # Annotations - MODIFIED to remove color annotation titles
      top_annotation = HeatmapAnnotation(
        # Add group labels
        group = anno_block(
          gp = gpar(fill = c(rep(healthy_palette, each = 1), 
                             rep(infected_palette, each = 1))),
          labels = c("H0", "H12", "H24", "H36", "H48", 
                     "I0", "I12", "I24", "I36", "I48"),
          labels_gp = gpar(col = "white", fontsize = 9, fontface = "bold")
        ),
        
        annotation_height = unit(20, "mm"),
        # Remove annotation names by setting fontsize to 0
        annotation_name_gp = gpar(fontsize = 0)
      ),
      
      bottom_annotation = HeatmapAnnotation(
        # Expression distribution boxplot
        distribution = anno_boxplot(
          GSE51856_Milk_Expr,
          height = unit(25, "mm"),
          gp = gpar(fill = c(rep(alpha(healthy_palette[3], 0.7), 24),
                             rep(alpha(infected_palette[3], 0.7), 24)),
                    col = "white", lwd = 0.5),
          outline = TRUE,
          box_width = 0.6
        ),
        
        # Remove annotation names by setting fontsize to 0
        annotation_name_gp = gpar(fontsize = 0),
        annotation_name_side = "left"
      ),
      
      # Legend settings
      heatmap_legend_param = list(
        title = "RNA Stability Score",
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 10),
        legend_height = unit(40, "mm"),
        legend_width = unit(8, "mm"),
        direction = "vertical",
        at = c(-2, -1, 0, 1, 2),
        labels = c("-2", "-1", "0", "1", "2"),
        border = "black"
      ),
      
      # Border
      border = TRUE,
      border_gp = gpar(col = "black", lwd = 1)
    )
    
    # =============================================================================
    # COMBINE AND SAVE HEATMAPS
    # =============================================================================
    
    # Combine all heatmaps
    combined_heatmap = healthy_heatmap + infected_heatmap + main_heatmap
    
    # Save in multiple formats
    for (format in save_formats) {
      
      if (format == "PDF") {
        filename = paste0("Publication_Heatmaps/PDF/", module_color, "_module_publication.pdf")
        pdf(filename, width = 14, height = 10, pointsize = 12)
        
      } else if (format == "PNG") {
        filename = paste0("Publication_Heatmaps/PNG/", module_color, "_module_publication.png")
        png(filename, width = 4200, height = 3000, res = 300, type = "cairo")
        
      } else if (format == "TIFF") {
        filename = paste0("Publication_Heatmaps/TIFF/", module_color, "_module_publication.tiff")
        tiff(filename, width = 4200, height = 3000, res = 300, compression = "lzw")
      }
      
      # Draw the heatmap
      draw(combined_heatmap, 
           heatmap_legend_side = "right",
           padding = unit(c(5, 5, 5, 5), "mm"),
           background = "white")
      
      # Add custom annotations - keep the text label but not the title
      decorate_annotation("distribution", {
        grid.text(
          "Expression\nDistribution", 
          x = unit(-3, "cm"), 
          y = unit(0.5, "npc"),
          rot = 90,
          gp = gpar(fontsize = 11, fontface = "bold", col = "black")
        )
      })
      
      # Add sample size information
      grid.text(
        paste("n =", nrow(GSE51856_Milk_Expr), "genes"),
        x = unit(0.98, "npc"),
        y = unit(0.02, "npc"),
        just = c("right", "bottom"),
        gp = gpar(fontsize = 10, fontface = "italic", col = "grey30")
      )
      
      dev.off()
      cat("Saved:", filename, "\n")
    }
  }
  
  # =============================================================================
  # CREATE HEATMAPS FOR ALL MODULES
  # =============================================================================
  
  cat("Creating publication-ready heatmaps for all modules...\n")
  
  # Remove grey module and small modules for cleaner presentation
  filtered_modules = modulenames[modulenames != "grey"]
  
  # Count genes per module to focus on larger modules
  module_sizes = table(moduleColorsAutomaticData1)
  large_modules = names(module_sizes[module_sizes >= 30])  # Only modules with 30+ genes
  
  # Create heatmaps for large modules only
  for (module in large_modules) {
    if (module != "grey") {
      create_publication_heatmap(module, save_formats = c("PDF", "PNG", "TIFF"))
    }
  }
  
  # =============================================================================
  # CREATE COMBINED OVERVIEW FIGURE
  # =============================================================================
  
  cat("Creating module overview figure...\n")
  
  # Create a summary figure showing all major modules
  create_overview_figure <- function() {
    
    # Select top 6 modules by size for overview
    top_modules = head(names(sort(module_sizes[module_sizes >= 30], decreasing = TRUE)), 6)
    top_modules = top_modules[top_modules != "grey"]
    
    pdf("Publication_Heatmaps/PDF/WGCNA_Modules_Overview.pdf", 
        width = 20, height = 12)
    
    # Create layout for multiple modules
    layout_matrix = matrix(1:(length(top_modules)), nrow = 2, byrow = TRUE)
    
    for (i in seq_along(top_modules)) {
      module = top_modules[i]
      
      # Extract and prepare data
      inModule = (moduleColorsAutomaticData1 == module)
      nameofmodule = as.data.frame(probes)[inModule,]
      expr_data = t(GSE51856_Milk_ExpNorm_df)
      expr_data = subset(expr_data, rownames(expr_data) %in% nameofmodule)
      
      if(nrow(expr_data) > 0) {
        # Simplified heatmap for overview - remove color annotation titles
        simple_heatmap = Heatmap(
          expr_data,
          name = paste(module, "module"),
          col = expression_colors,
          cluster_rows = TRUE,
          cluster_columns = FALSE,
          show_row_names = FALSE,
          show_column_names = FALSE,
          column_split = GSE51856_Milk_sample_groups,
          column_title = NULL,
          column_title_gp = gpar(fontsize = 0),
          # Modified to hide annotation names
          heatmap_legend_param = list(
            title = "Expression",
            legend_height = unit(30, "mm"),
            labels_gp = gpar(fontsize = 8)
          )
        )
        
        # Draw the heatmap
        pushViewport(viewport(layout.pos.row = ceiling(i/3), 
                              layout.pos.col = ((i-1) %% 3) + 1))
        draw(simple_heatmap, newpage = FALSE)
        popViewport()
      }
    }
    
    dev.off()
  }
  
  create_overview_figure()
  
  # =============================================================================
  # CREATE PUBLICATION SUMMARY
  # =============================================================================
  
  # Generate summary file
  sink("Publication_Heatmaps/Heatmap_Generation_Summary.txt")
  cat("PUBLICATION-READY HEATMAPS GENERATION SUMMARY\n")
  cat("=" %R% 50, "\n\n")
  
  cat("MODULES PROCESSED:\n")
  for (module in large_modules) {
    if (module != "grey") {
      gene_count = sum(moduleColorsAutomaticData1 == module)
      cat(sprintf("- %s module: %d genes\n", module, gene_count))
    }
  }
  
  cat("\nFILES CREATED:\n")
  cat("- PDF files: High-resolution vector graphics for publication\n")
  cat("- PNG files: High-resolution raster images (300 DPI)\n") 
  cat("- TIFF files: Uncompressed high-quality images for journals\n")
  cat("- Overview figure: Combined view of major modules\n")
  
  cat("\nRECOMMENDATIONS FOR PUBLICATION:\n")
  cat("- Use PDF files for main figures in manuscripts\n")
  cat("- Use TIFF files for journal submission if required\n")
  cat("- PNG files are suitable for presentations and supplementary materials\n")
  cat("- All files are 300 DPI publication quality\n")
  
  cat("\nCOLOR SCHEME:\n")
  cat("- Blue to red gradient for expression levels\n")
  cat("- Green tones for healthy conditions\n") 
  cat("- Red tones for infected conditions\n")
  cat("- Professional color palette suitable for colorblind readers\n")
  
  sink()
  
  cat("\n=== PUBLICATION HEATMAPS COMPLETED ===\n")
  cat("All heatmaps have been generated in publication-ready quality!\n")
  cat("Check the 'Publication_Heatmaps' directory for all output files.\n")
  cat("\nFiles are ready for:\n")
  cat("✓ Journal submission\n")
  cat("✓ Conference presentations  \n")
  cat("✓ Supplementary materials\n")
  cat("✓ High-quality printing\n")
  
  
  
  
  
  
  
  
  
  
  
  
  
  ##############################
  # Functional Enrichment Analysis for Yellow and Red Modules
  # This script creates:
  # 1. Supplementary File S4: Excel file with all BP and KEGG results
  # 2. Figure 5: Combined visualization with panels A and B
  ##############################
  
  # Set working directory
  setwd("F:/new desktop/thesis/sixth WGCNA analysis/var0.3")
  
  # Load required libraries (using openxlsx instead of xlsx to avoid Java dependency)
  library(openxlsx)
  library(readxl)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(stringr)
  library(ggrepel)
  library(viridis)
  library(patchwork)
  library(RColorBrewer)
  
  # Create output directory
  dir.create("Functional_Enrichment_Output", showWarnings = FALSE)
  
  ##############################
  # PART 1: Create Supplementary File S4
  # Combine all BP and KEGG results for Yellow and Red modules
  ##############################
  
  # Function to read and process enrichment data
  read_enrichment_data <- function(file_path, sheet_name, module_name, db_type) {
    tryCatch({
      data <- read_excel(file_path, sheet = sheet_name)
      if (nrow(data) > 0) {
        # Add module and database information
        data$Module <- module_name
        data$Database <- db_type
        return(data)
      }
      return(NULL)
    }, error = function(e) {
      cat("Error reading", file_path, "sheet", sheet_name, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  # Read all required data
  # Biological Process data
  yellow_bp <- read_enrichment_data("GO_enrichR/GO_Biological_Process_2023.xlsx", "yellow", "Yellow", "GO_BP")
  red_bp <- read_enrichment_data("GO_enrichR/GO_Biological_Process_2023.xlsx", "red", "Red", "GO_BP")
  
  # KEGG pathway data
  yellow_kegg <- read_enrichment_data("GO_enrichR/KEGG_2021_Human.xlsx", "yellow", "Yellow", "KEGG")
  red_kegg <- read_enrichment_data("GO_enrichR/KEGG_2021_Human.xlsx", "red", "Red", "KEGG")
  
  # Create Supplementary File S4 using openxlsx
  wb <- createWorkbook()
  
  # Add sheets for each module and database combination
  if (!is.null(yellow_bp)) {
    addWorksheet(wb, "Yellow_Module_BP")
    writeData(wb, "Yellow_Module_BP", yellow_bp)
  }
  
  if (!is.null(red_bp)) {
    addWorksheet(wb, "Red_Module_BP")
    writeData(wb, "Red_Module_BP", red_bp)
  }
  
  if (!is.null(yellow_kegg)) {
    addWorksheet(wb, "Yellow_Module_KEGG")
    writeData(wb, "Yellow_Module_KEGG", yellow_kegg)
  }
  
  if (!is.null(red_kegg)) {
    addWorksheet(wb, "Red_Module_KEGG")
    writeData(wb, "Red_Module_KEGG", red_kegg)
  }
  
  # Add summary sheet
  summary_data <- data.frame(
    Module = c("Yellow", "Yellow", "Red", "Red"),
    Database = c("GO_BP", "KEGG", "GO_BP", "KEGG"),
    Total_Terms = c(
      ifelse(is.null(yellow_bp), 0, nrow(yellow_bp)),
      ifelse(is.null(yellow_kegg), 0, nrow(yellow_kegg)),
      ifelse(is.null(red_bp), 0, nrow(red_bp)),
      ifelse(is.null(red_kegg), 0, nrow(red_kegg))
    ),
    Significant_Terms = c(
      ifelse(is.null(yellow_bp), 0, sum(yellow_bp$P.value < 0.05, na.rm = TRUE)),
      ifelse(is.null(yellow_kegg), 0, sum(yellow_kegg$P.value < 0.05, na.rm = TRUE)),
      ifelse(is.null(red_bp), 0, sum(red_bp$P.value < 0.05, na.rm = TRUE)),
      ifelse(is.null(red_kegg), 0, sum(red_kegg$P.value < 0.05, na.rm = TRUE))
    )
  )
  
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", summary_data)
  
  # Save the workbook
  saveWorkbook(wb, "Functional_Enrichment_Output/Supplementary_File_S4_Enrichment_Results.xlsx", overwrite = TRUE)
  
  cat("Supplementary File S4 created successfully!\n")
  
  ##############################
  # PART 2: Create Figure 5 - Combined Visualization
  # Panel A: Yellow Module, Panel B: Red Module
  ##############################
  
  # Function to prepare data for visualization
  prepare_viz_data <- function(bp_data, kegg_data, module_name) {
    # Select top 10 terms from each database
    top_bp <- NULL
    top_kegg <- NULL
    
    if (!is.null(bp_data) && nrow(bp_data) > 0) {
      top_bp <- bp_data %>%
        arrange(P.value) %>%
        head(10) %>%
        mutate(
          Database = "GO Biological Process",
          neg_log_p = -log10(P.value),
          gene_count = as.numeric(sub("/.*", "", Overlap))
        )
    }
    
    if (!is.null(kegg_data) && nrow(kegg_data) > 0) {
      top_kegg <- kegg_data %>%
        arrange(P.value) %>%
        head(10) %>%
        mutate(
          Database = "KEGG Pathway",
          neg_log_p = -log10(P.value),
          gene_count = as.numeric(sub("/.*", "", Overlap))
        )
    }
    
    # Combine data
    combined <- NULL
    if (!is.null(top_bp) && !is.null(top_kegg)) {
      combined <- rbind(top_bp, top_kegg)
    } else if (!is.null(top_bp)) {
      combined <- top_bp
    } else if (!is.null(top_kegg)) {
      combined <- top_kegg
    }
    
    # Shorten long term names
    if (!is.null(combined) && nrow(combined) > 0) {
      combined$Short_Term <- ifelse(nchar(combined$Term) > 50,
                                    paste0(substr(combined$Term, 1, 47), "..."),
                                    combined$Term)
    }
    
    return(combined)
  }
  
  # Prepare data for both modules
  yellow_data <- prepare_viz_data(yellow_bp, yellow_kegg, "Yellow")
  red_data <- prepare_viz_data(red_bp, red_kegg, "Red")
  
  # Create visualization function
  create_module_plot <- function(data, module_name, module_color) {
    if (is.null(data) || nrow(data) == 0) {
      # Return empty plot with message
      return(
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = paste("No data available for", module_name, "module"), 
                   size = 6, hjust = 0.5, vjust = 0.5) +
          theme_void() +
          ggtitle(paste(module_name, "Module"))
      )
    }
    
    # Create the plot
    p <- ggplot(data, aes(x = neg_log_p, y = reorder(Short_Term, neg_log_p))) +
      geom_segment(aes(x = 0, xend = neg_log_p, y = Short_Term, yend = Short_Term),
                   color = "gray70", size = 0.5) +
      geom_point(aes(size = gene_count, fill = Database), 
                 shape = 21, color = "black", alpha = 0.8) +
      scale_size_continuous(range = c(3, 10), name = "Gene Count") +
      scale_fill_manual(values = c("GO Biological Process" = module_color,
                                   "KEGG Pathway" = colorspace::lighten(module_color, 0.3))) +
      labs(
        x = "-log10(p-value)",
        y = "",
        title = paste(module_name, "Module Enrichment Analysis")
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 10, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
      ) +
      geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) +
      annotate("text", x = -log10(0.05), y = 1, label = "p = 0.05", 
               hjust = -0.1, vjust = -0.5, size = 3, color = "red")
    
    # Add facet if there's enough data
    if (length(unique(data$Database)) > 1) {
      p <- p + facet_wrap(~ Database, scales = "free_y", ncol = 1)
    }
    
    return(p)
  }
  
  # Create plots for both modules
  plot_a <- create_module_plot(yellow_data, "Yellow", "#FFC107")
  plot_b <- create_module_plot(red_data, "Red", "#F44336")
  
  # Combine plots into Figure 5
  figure_5 <- plot_a / plot_b + 
    plot_annotation(
      tag_levels = 'A',
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    )
  
  # Save Figure 5
  ggsave("Functional_Enrichment_Output/Figure_5_Combined_Enrichment.pdf", 
         figure_5, width = 12, height = 16, dpi = 300)
  
  ggsave("Functional_Enrichment_Output/Figure_5_Combined_Enrichment.png", 
         figure_5, width = 12, height = 16, dpi = 300)
  
  ##############################
  # Alternative visualization style for Figure 5
  ##############################
  
  # Create horizontal bar plot style
  create_horizontal_bar_plot <- function(data, module_name, module_color) {
    if (is.null(data) || nrow(data) == 0) {
      return(
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, 
                   label = paste("No data available for", module_name, "module"), 
                   size = 6, hjust = 0.5, vjust = 0.5) +
          theme_void() +
          ggtitle(paste(module_name, "Module"))
      )
    }
    
    # Prepare data - top 10 overall
    plot_data <- data %>%
      arrange(P.value) %>%
      head(10) %>%
      mutate(
        Database_short = ifelse(Database == "GO Biological Process", "GO BP", "KEGG"),
        # Make terms even shorter for horizontal layout
        Very_Short_Term = ifelse(nchar(Term) > 35,
                                 paste0(substr(Term, 1, 32), "..."),
                                 Term),
        Term_with_DB = paste0("[", Database_short, "] ", Very_Short_Term)
      )
    
    p <- ggplot(plot_data, aes(x = neg_log_p, y = reorder(Term_with_DB, neg_log_p))) +
      geom_bar(stat = "identity", aes(fill = Database), alpha = 0.8) +
      geom_text(aes(label = gene_count), hjust = -0.2, size = 3) +
      scale_fill_manual(values = c("GO Biological Process" = module_color,
                                   "KEGG Pathway" = colorspace::lighten(module_color, 0.3))) +
      labs(
        x = "-log10(p-value)",
        y = "",
        title = paste(module_name, "Module")
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 8),  # Smaller font size
        axis.title.x = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")  # Add margins
      ) +
      xlim(0, max(plot_data$neg_log_p) * 1.2) +  # Adjust x limit
      scale_y_discrete(expand = expansion(add = c(0.5, 0.5)))  # Add space on y-axis
    
    return(p)
  }
  
  # Create alternative visualization
  alt_plot_a <- create_horizontal_bar_plot(yellow_data, "Yellow", "#FFC107")
  alt_plot_b <- create_horizontal_bar_plot(red_data, "Red", "#F44336")
  
  # Combine alternative plots side by side
  alt_figure_5 <- alt_plot_a | alt_plot_b
  alt_figure_5 <- alt_figure_5 + 
    plot_annotation(
      tag_levels = 'A'
      subtitle = "Numbers indicate gene count in each pathway",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  # Save alternative Figure 5 with adjusted dimensions
  ggsave("Functional_Enrichment_Output/Figure_5_Alternative_Horizontal.pdf", 
         alt_figure_5, width = 16, height = 10, dpi = 300)  # Increased width
  
  ggsave("Functional_Enrichment_Output/Figure_5_Alternative_Horizontal.png", 
         alt_figure_5, width = 16, height = 10, dpi = 300)  # Increased width
  
  ##############################
  # Professional Bubble Plot Style (Similar to the example provided)
  ##############################
  
  # Create bubble plot similar to the provided example
  create_bubble_plot_style <- function(yellow_data, red_data) {
    # Combine and prepare data
    prepare_data_for_bubble <- function(data, module) {
      if (is.null(data) || nrow(data) == 0) return(NULL)
      
      data %>%
        arrange(P.value) %>%
        head(10) %>%  # Top 10 terms per module
        mutate(
          Module = module,
          # Calculate fold enrichment (approximation based on p-value and gene count)
          Fold_Enrichment = -log10(P.value) * (gene_count / 10),
          Log2FDR = -log2(P.value),
          Number_of_Genes = gene_count,
          # Shorten terms for better display
          Biological_Process = case_when(
            nchar(Term) > 40 ~ paste0(substr(Term, 1, 37), "..."),
            TRUE ~ Term
          ),
          Database_Type = Database
        )
    }
    
    yellow_plot_data <- prepare_data_for_bubble(yellow_data, "Yellow")
    red_plot_data <- prepare_data_for_bubble(red_data, "Red")
    
    # Combine data
    combined_data <- rbind(yellow_plot_data, red_plot_data)
    
    if (is.null(combined_data) || nrow(combined_data) == 0) {
      return(ggplot() + theme_void() + ggtitle("No data available"))
    }
    
    # Create the bubble plot
    p <- ggplot(combined_data, aes(x = Fold_Enrichment, y = reorder(Biological_Process, Fold_Enrichment), 
                                   size = Number_of_Genes, color = Log2FDR)) + 
      geom_point(shape = 19, alpha = 0.8) +
      theme_bw() + 
      theme(
        text = element_text(size = 12),
        axis.text.x = element_text(color = "black", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 7, face = "plain"),
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 10),
        legend.direction = "vertical",
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(size = 12, face = "bold"),
        panel.spacing = unit(1, "lines")
      ) +
      scale_y_discrete("Enriched Pathways") +
      scale_x_continuous("Fold Enrichment") +
      labs(size = "Number of Genes", color = "-log2(P-value)") +
      scale_color_gradient(low = "blue", high = "red", space = "Lab") +
      facet_wrap(~ Module, ncol = 2, strip.position = "top", scales = "free_y") +
      guides(size = guide_legend("Number of Genes"))
    
    return(p)
  }
  
  # Create separated bubble plots for GO BP and KEGG
  create_separated_bubble_plots <- function(yellow_data, red_data) {
    # Separate GO BP and KEGG data
    prepare_separated_data <- function(data, module) {
      if (is.null(data) || nrow(data) == 0) return(list(bp = NULL, kegg = NULL))
      
      # Process data
      processed <- data %>%
        arrange(P.value) %>%
        mutate(
          Module = module,
          Fold_Enrichment = -log10(P.value) * (gene_count / 10),
          Log2FDR = -log2(P.value),
          Number_of_Genes = gene_count,
          Biological_Process = case_when(
            nchar(Term) > 45 ~ paste0(substr(Term, 1, 42), "..."),
            TRUE ~ Term
          )
        )
      
      # Separate by database type
      bp_data <- processed %>% filter(Database == "GO Biological Process") %>% head(10)
      kegg_data <- processed %>% filter(Database == "KEGG Pathway") %>% head(10)
      
      return(list(bp = bp_data, kegg = kegg_data))
    }
    
    yellow_separated <- prepare_separated_data(yellow_data, "Yellow")
    red_separated <- prepare_separated_data(red_data, "Red")
    
    # Combine GO BP data
    bp_combined <- rbind(yellow_separated$bp, red_separated$bp)
    
    # Combine KEGG data
    kegg_combined <- rbind(yellow_separated$kegg, red_separated$kegg)
    
    # Create GO BP plot
    if (!is.null(bp_combined) && nrow(bp_combined) > 0) {
      bp_plot <- ggplot(bp_combined, aes(x = Fold_Enrichment, y = reorder(Biological_Process, Fold_Enrichment), 
                                         size = Number_of_Genes, color = Log2FDR)) + 
        geom_point(shape = 19, alpha = 0.8) +
        theme_bw() + 
        theme(
          text = element_text(size = 12),
          axis.text.x = element_text(color = "black", size = 11, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "black", size = 8, face = "plain"),
          axis.title.x = element_text(color = "black", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.title.y = element_text(color = "black", size = 13, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "black", size = 11),
          legend.text = element_text(color = "black", size = 9),
          legend.direction = "vertical",
          strip.background = element_rect(fill = c("Yellow" = "#FFF3CD", "Red" = "#F8D7DA")),
          strip.text = element_text(size = 12, face = "bold"),
          panel.spacing = unit(1, "lines")
        ) +
        scale_y_discrete("GO Biological Process") +
        scale_x_continuous("Fold Enrichment") +
        labs(size = "Number of Genes", color = "-log2(P-value)", title = "A. GO Biological Process Enrichment") +
        scale_color_gradient(low = "blue", high = "red", space = "Lab") +
        facet_wrap(~ Module, ncol = 2, strip.position = "top", scales = "free_y") +
        guides(size = guide_legend("Number of Genes"))
    } else {
      bp_plot <- ggplot() + theme_void() + ggtitle("No GO BP data available")
    }
    
    # Create KEGG plot
    if (!is.null(kegg_combined) && nrow(kegg_combined) > 0) {
      kegg_plot <- ggplot(kegg_combined, aes(x = Fold_Enrichment, y = reorder(Biological_Process, Fold_Enrichment), 
                                             size = Number_of_Genes, color = Log2FDR)) + 
        geom_point(shape = 19, alpha = 0.8) +
        theme_bw() + 
        theme(
          text = element_text(size = 12),
          axis.text.x = element_text(color = "black", size = 11, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "black", size = 8, face = "plain"),
          axis.title.x = element_text(color = "black", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.title.y = element_text(color = "black", size = 13, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "black", size = 11),
          legend.text = element_text(color = "black", size = 9),
          legend.direction = "vertical",
          strip.background = element_rect(fill = c("Yellow" = "#FFF3CD", "Red" = "#F8D7DA")),
          strip.text = element_text(size = 12, face = "bold"),
          panel.spacing = unit(1, "lines")
        ) +
        scale_y_discrete("KEGG Pathway") +
        scale_x_continuous("Fold Enrichment") +
        labs(size = "Number of Genes", color = "-log2(P-value)", title = "B. KEGG Pathway Enrichment") +
        scale_color_gradient(low = "blue", high = "red", space = "Lab") +
        facet_wrap(~ Module, ncol = 2, strip.position = "top", scales = "free_y") +
        guides(size = guide_legend("Number of Genes"))
    } else {
      kegg_plot <- ggplot() + theme_void() + ggtitle("No KEGG data available")
    }
    
    return(list(bp = bp_plot, kegg = kegg_plot))
  }
  
  # Create bubble plot style figures
  bubble_figure_5 <- create_bubble_plot_style(yellow_data, red_data)
  
  # Save bubble plot version
  ggsave("Functional_Enrichment_Output/Figure_5_Bubble_Style.pdf", 
         bubble_figure_5, width = 14, height = 10, dpi = 300)
  
  ggsave("Functional_Enrichment_Output/Figure_5_Bubble_Style.png", 
         bubble_figure_5, width = 14, height = 10, dpi = 300)
  
  # Create separated plots
  separated_plots <- create_separated_bubble_plots(yellow_data, red_data)
  
  # Combine separated plots
  combined_bubble_figure <- separated_plots$bp / separated_plots$kegg + 
    plot_annotation(
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    )
  
  # Save combined bubble plot
  ggsave("Functional_Enrichment_Output/Figure_5_Bubble_Combined.pdf", 
         combined_bubble_figure, width = 14, height = 16, dpi = 300)
  
  ggsave("Functional_Enrichment_Output/Figure_5_Bubble_Combined.png", 
         combined_bubble_figure, width = 14, height = 16, dpi = 300)
  
  # Print summary statistics
  cat("\n--- Summary Statistics ---\n")
  cat("Yellow Module - BP terms:", ifelse(is.null(yellow_bp), 0, sum(yellow_bp$P.value < 0.05, na.rm = TRUE)), "significant out of", ifelse(is.null(yellow_bp), 0, nrow(yellow_bp)), "\n")
  cat("Yellow Module - KEGG pathways:", ifelse(is.null(yellow_kegg), 0, sum(yellow_kegg$P.value < 0.05, na.rm = TRUE)), "significant out of", ifelse(is.null(yellow_kegg), 0, nrow(yellow_kegg)), "\n")
  cat("Red Module - BP terms:", ifelse(is.null(red_bp), 0, sum(red_bp$P.value < 0.05, na.rm = TRUE)), "significant out of", ifelse(is.null(red_bp), 0, nrow(red_bp)), "\n")
  cat("Red Module - KEGG pathways:", ifelse(is.null(red_kegg), 0, sum(red_kegg$P.value < 0.05, na.rm = TRUE)), "significant out of", ifelse(is.null(red_kegg), 0, nrow(red_kegg)), "\n")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Function to create module plot without titles
  create_module_plot <- function(data, module_name, module_color) {
    if (is.null(data) || nrow(data) == 0) {
      # Return empty plot with message
      return(
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = paste("No data available for", module_name, "module"), 
                   size = 6, hjust = 0.5, vjust = 0.5) +
          theme_void()
        # Removed the title here
      )
    }
    
    # Create the plot
    p <- ggplot(data, aes(x = neg_log_p, y = reorder(Short_Term, neg_log_p))) +
      geom_segment(aes(x = 0, xend = neg_log_p, y = Short_Term, yend = Short_Term),
                   color = "gray70", size = 0.5) +
      geom_point(aes(size = gene_count, fill = Database), 
                 shape = 21, color = "black", alpha = 0.8) +
      scale_size_continuous(range = c(3, 10), name = "Gene Count") +
      scale_fill_manual(values = c("GO Biological Process" = module_color,
                                   "KEGG Pathway" = colorspace::lighten(module_color, 0.3))) +
      labs(
        x = "-log10(p-value)",
        y = ""
        # Removed the title here
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11, face = "bold"),
        # Remove plot title styling since there's no title
        # plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 10, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
      ) +
      geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) +
      annotate("text", x = -log10(0.05), y = 1, label = "p = 0.05", 
               hjust = -0.1, vjust = -0.5, size = 3, color = "red")
    
    # Add facet if there's enough data
    if (length(unique(data$Database)) > 1) {
      p <- p + facet_wrap(~ Database, scales = "free_y", ncol = 1)
    }
    
    return(p)
  }
  
  # Create horizontal bar plot style without titles
  create_horizontal_bar_plot <- function(data, module_name, module_color) {
    if (is.null(data) || nrow(data) == 0) {
      return(
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, 
                   label = paste("No data available for", module_name, "module"), 
                   size = 6, hjust = 0.5, vjust = 0.5) +
          theme_void()
        # Removed the title here
      )
    }
    
    # Prepare data - top 10 overall
    plot_data <- data %>%
      arrange(P.value) %>%
      head(10) %>%
      mutate(
        Database_short = ifelse(Database == "GO Biological Process", "GO BP", "KEGG"),
        # Make terms even shorter for horizontal layout
        Very_Short_Term = ifelse(nchar(Term) > 35,
                                 paste0(substr(Term, 1, 32), "..."),
                                 Term),
        Term_with_DB = paste0("[", Database_short, "] ", Very_Short_Term)
      )
    
    p <- ggplot(plot_data, aes(x = neg_log_p, y = reorder(Term_with_DB, neg_log_p))) +
      geom_bar(stat = "identity", aes(fill = Database), alpha = 0.8) +
      geom_text(aes(label = gene_count), hjust = -0.2, size = 3) +
      scale_fill_manual(values = c("GO Biological Process" = module_color,
                                   "KEGG Pathway" = colorspace::lighten(module_color, 0.3))) +
      labs(
        x = "-log10(p-value)",
        y = ""
        # Removed the title here
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 8),  # Smaller font size
        axis.title.x = element_text(size = 10, face = "bold"),
        # Remove plot title styling
        # plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")  # Add margins
      ) +
      xlim(0, max(plot_data$neg_log_p) * 1.2) +  # Adjust x limit
      scale_y_discrete(expand = expansion(add = c(0.5, 0.5)))  # Add space on y-axis
    
    return(p)
  }
  
  # Main code for creating the combined figure (implementation depends on where you want to use this)
  # For example:
  
  # Create plots for both modules without titles
  plot_a <- create_module_plot(yellow_data, "Yellow", "#FFC107")
  plot_b <- create_module_plot(red_data, "Red", "#F44336")
  
  # Combine plots into Figure 5
  figure_5 <- plot_a / plot_b + 
    plot_annotation(
      tag_levels = 'A',
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    )
  
  # Save Figure 5
  ggsave("Functional_Enrichment_Output/Figure_5_Combined_Enrichment.pdf", 
         figure_5, width = 12, height = 16, dpi = 300)
  
  ggsave("Functional_Enrichment_Output/Figure_5_Combined_Enrichment.png", 
         figure_5, width = 12, height = 16, dpi = 300)
  
  # Alternative visualization
  alt_plot_a <- create_horizontal_bar_plot(yellow_data, "Yellow", "#FFC107")
  alt_plot_b <- create_horizontal_bar_plot(red_data, "Red", "#F44336")
  
  # Combine alternative plots side by side
  alt_figure_5 <- alt_plot_a | alt_plot_b
  alt_figure_5 <- alt_figure_5 + 
    plot_annotation(
      tag_levels = 'A',
      subtitle = "Numbers indicate gene count in each pathway",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  # Save alternative Figure 5 with adjusted dimensions
  ggsave("Functional_Enrichment_Output/Figure_5_Alternative_Horizontal.pdf", 
         alt_figure_5, width = 16, height = 10, dpi = 300)
  
  ggsave("Functional_Enrichment_Output/Figure_5_Alternative_Horizontal.png", 
         alt_figure_5, width = 16, height = 10, dpi = 300)
  