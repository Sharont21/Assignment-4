##*************************
## BINF 6210
##
## Assignment 4
##
## Sharon Tsukernik
##
## Friday, December 5, 2025

# ---- Packages used -----

library(tidyverse)
library(ggplot2)
#install.packages("randomForest")
library(randomForest)
#install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("Biostrings")
library(Biostrings)
#install.packages("rentrez")
library(rentrez)
#install.packages("pROC")
library(pROC)
#install.packages("e1071")
library(e1071)

# ---- Loading data from BOLD ----
## ---- Load butterfly data ----

dfDanainae <- read_tsv("../data/Danainae.tsv")
dfHeliconiinae <- read_tsv("../data/Heliconiinae.tsv")

# ---- Exploratory analysis on BOLD data ----

# Getting summary of datasets
summary(dfDanainae)
summary(dfHeliconiinae)

# Looking at column names
names(dfDanainae)
names(dfHeliconiinae)

# Getting the counts for each gene
table(dfDanainae$marker_code)
#6112 COI-5P and 8 COII
table(dfHeliconiinae$marker_code)
#10127 COI-5P and 31 COII

# less COII in comparison to COI, so sampling bias is involved

# I will be using COI-5P as my COI gene

# ---- Downloading data from NCBI ----

## ---- Download COI sequences for each subfamily of butterfly ----

# Creating a reusable function to download and save to a file
download_gene_fasta <- function(organism, gene, outfile, retmax = 1000) {
  
  # Creating NCBI search term
  search_term <- paste0(organism, "[Organism] AND ", gene, "[Gene]")
  
  # Searching NCBI nucleotides
  search_res <- entrez_search(
    db = "nucleotide",
    term = search_term,
    retmax = retmax,
    use_history = TRUE
  )
  
  # Getting fasta sequences
  fasta_res <- entrez_fetch(
    db = "nucleotide",
    web_history = search_res$web_history,
    rettype = "fasta"
  )
  
  # Writing to a file
  write(fasta_res, file = outfile)
}

# Using function to download COI and COII genes for both butterfly subfamily (Danainae and Heliconiinae)
download_gene_fasta(
  organism = "Danainae",
  gene     = "COI",
  outfile  = "../data/Danainae_COI_NCBI.fasta"
)

download_gene_fasta(
  organism = "Heliconiinae",
  gene     = "COI",
  outfile  = "../data/Heliconiinae_COI_NCBI.fasta"
)

download_gene_fasta(
  organism = "Danainae",
  gene     = "COII",
  outfile  = "../data/Danainae_COII_NCBI.fasta"
)

download_gene_fasta(
  organism = "Heliconiinae",
  gene     = "COII",
  outfile  = "../data/Heliconiinae_COII_NCBI.fasta"
)

# ---- Load FASTA files into R ----

Danainae_COI <- readDNAStringSet("../data/Danainae_COI_NCBI.fasta")
Heliconiinae_COI <- readDNAStringSet("../data/Heliconiinae_COI_NCBI.fasta")

Danainae_COII <- readDNAStringSet("../data/Danainae_COII_NCBI.fasta")
Heliconiinae_COII <- readDNAStringSet("../data/Heliconiinae_COII_NCBI.fasta")

# ---- Converting each FASTA to a tibble ----

fasta_to_df <- function(x, subfamily, gene) {
  tibble(
    id = names(x),
    nuc = as.character(x),
    subfamily = subfamily,
    gene = gene
  )
}

df_Dan_COI <- fasta_to_df(Danainae_COI, "Danainae", "COI")
df_Hel_COI <- fasta_to_df(Heliconiinae_COI, "Heliconiinae", "COI")

df_Dan_COII <- fasta_to_df(Danainae_COII, "Danainae", "COII")
df_Hel_COII <- fasta_to_df(Heliconiinae_COII, "Heliconiinae", "COII")

# Removing unused variables
rm(Danainae_COI, Danainae_COII, Heliconiinae_COI, Heliconiinae_COII)

# ---- Standardizing the BOLD data ----

dfDan_BOLD <- dfDanainae %>%
  filter(!is.na(nuc)) %>%
  transmute(
    id = processid,
    nuc = nuc,
    subfamily = "Danainae",
    gene = marker_code
  )

dfHel_BOLD <- dfHeliconiinae %>%
  filter(!is.na(nuc)) %>%
  transmute(
    id = processid,
    nuc = nuc,
    subfamily = "Heliconiinae", 
    gene = marker_code
  )

# Standardizing gene name across data sources by renaming BOLD's "COI-5P" to "COI" so it matches NCBI gene annotation
dfDan_BOLD$gene[dfDan_BOLD$gene == "COI-5P"] <- "COI"
dfHel_BOLD$gene[dfHel_BOLD$gene == "COI-5P"] <- "COI"

# Combining all BOLD and NCBI sequences into one unified dataframe
df_all <- bind_rows(
  dfDan_BOLD,
  dfHel_BOLD,
  df_Dan_COI,
  df_Hel_COI,
  df_Dan_COII,
  df_Hel_COII
)

# Removing unused variables
rm(dfDan_BOLD, dfHel_BOLD, df_Dan_COI, df_Hel_COI, df_Dan_COII, df_Hel_COII)

# ---- Exploratory analysis on combined data ----

summary(df_all)
table(df_all$gene)
# 25796 COI and 2832 COII

table(df_all$subfamily)
# 11217 Danainae and 20067 Heliconiinae


# Histogram of sequence lengths for COI
Len_Dan_COI <- df_all %>%
  filter(gene == "COI",
         subfamily == "Danainae",
         nchar(nuc) >= 600,
         nchar(nuc) <= 800) %>%
  pull(nuc) %>%
  nchar()

Len_Hel_COI <- df_all %>%
  filter(gene == "COI",
         subfamily == "Heliconiinae",
         nchar(nuc) >= 600,
         nchar(nuc) <= 800) %>%
  pull(nuc) %>%
  nchar()

# Restricting to 600–800 bp to match the expected insect COI barcode length and remove outliers (partial or whole-mitochondrial sequences)

png("../figs/Histogram_COI.png", width = 800, height = 600, res = 120)

hist(Len_Dan_COI,
     xlab = "Sequence Length",
     ylab = "Frequency",
     main = "COI Sequence Lengths: Danainae vs Heliconiinae",
     col = rgb(0, 0, 1, 0.4),   
     border = "blue"
)

# Adding Heliconiinae on top
hist(Len_Hel_COI,
     col = rgb(1, 0, 0, 0.4),   
     border = "red",
     add = TRUE
)

legend(
  "topright",
  legend = c("Danainae", "Heliconiinae"),
  fill = c(rgb(0,0,1,0.4), rgb(1,0,0,0.4)),
  border = c("blue", "red"),
  bty = "n", 
  cex = 0.9
  )

dev.off()

# Interpretation: Most COI sequences for both Danainae and Heliconiinae cluster tightly around ~650 bp, which is the expected insect barcode length. Heliconiinae shows slightly greater length variability, while Danainae is more uniform. 

# Histogram of sequence lengths for COII
Len_Dan_COII <- nchar(df_all$nuc[df_all$gene == "COII" & df_all$subfamily == "Danainae" &
                                    nchar(df_all$nuc) <= 800])
Len_Hel_COII <- nchar(df_all$nuc[df_all$gene == "COII" & df_all$subfamily == "Heliconiinae" &
                                   nchar(df_all$nuc) <= 800])

# Restricting to 600–800 bp to match the expected insect COII barcode length and remove outliers (partial or whole-mitochondrial sequences)

# Adding a big right margin for legend
png("../figs/Histogram_COII.png", width = 800, height = 600, res = 120)

par(mar = c(5, 5, 4, 8))

hist(Len_Dan_COII,
     xlab = "Sequence Length",
     ylab = "Frequency",
     main = "COII Sequence Lengths: Danainae vs Heliconiinae",
     col = rgb(0, 0, 1, 0.4),
     border = "blue"
)

hist(Len_Hel_COII,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     add = TRUE
)

# Allowing drawing outside plot for the legend
par(xpd = NA)

legend("topright",
       inset = c(-0.40, 0), 
       legend = c("Danainae", "Heliconiinae"),
       fill = c(rgb(0,0,1,0.4), rgb(1,0,0,0.4)),
       border = c("blue", "red"),
       bty = "n",
       cex = 0.9)

dev.off()

# Interpretation: COII sequence lengths in both Danainae and Heliconiinae are tightly clustered, with substantial overlap between groups, indicating highly conserved barcode lengths and no strong subfamily-specific length differences

# Removing unused variables
rm(Len_Dan_COI, Len_Dan_COII, Len_Hel_COI, Len_Hel_COII)

# ---- Clean and filter sequences ----

# Removing sequences with gaps or N's at the ends
df_all <- df_all %>%
  mutate(nuc2 = str_remove(nuc, "^[-N]+")) %>%
  mutate(nuc2 = str_remove(nuc2, "[-N]+$")) %>%
  mutate(nuc2 = str_remove_all(nuc2, "-+"))

# Removing sequences with >5% Ns
df_all <- df_all %>%
  filter(str_count(nuc2, "N") <= 0.05 * nchar(nuc))

# Filtering to standard barcode-length sequences for insects (600–800 bp) to ensure comparability across subfamilies and databases
df_all <- df_all %>%
  filter(
    (gene == "COI"  & nchar(nuc) >= 600 & nchar(nuc) <= 800) |
      (gene == "COII" & nchar(nuc) >= 600 & nchar(nuc) <= 800)
  )

# Comparing old and new sequence columns
df_seq_compare <- cbind(df_all$nuc, df_all$nuc2)
view(df_seq_compare)
# Successfully removed sequences with gaps or N's

rm(df_seq_compare)

# Building quartiles to restrict COI sequence lengths
q1 <- quantile(nchar(df_all$nuc2[df_all$gene == "COI"]), 0.25)
q3 <- quantile(nchar(df_all$nuc2[df_all$gene == "COI"]), 0.75)

# Filtering COI sequences to the interquartile length range to remove unusually short or long outliers, and retaining all COII sequences
df_clean <- df_all %>%
  filter(
    (gene == "COI" & nchar(nuc2) >= q1 & nchar(nuc2) <= q3 | gene == "COII")
  )

# Checking the maximum sequence length after filtering to confirm that outliers were removed
max(nchar(df_clean$nuc))
#727 - falls within the expected 600–800 bp barcode range, confirming that long outlier sequences (full mitochondrial genomes) were successfully removed

# Converting sequences to DNAStringSet
df_clean <- as.data.frame(df_clean)
df_clean$nuc2 <- DNAStringSet(df_clean$nuc2)

# Extracting nucleotides and k-mer features
df_clean <- cbind(
  df_clean,
  as.data.frame(letterFrequency(df_clean$nuc2, letters = c("A", "C", "G", "T")))
)

# Building proportions
df_clean$Aprop <- df_clean$A / (df_clean$A + df_clean$C + df_clean$G + df_clean$T)
df_clean$Tprop <- df_clean$T / (df_clean$A + df_clean$C + df_clean$G + df_clean$T)
df_clean$Gprop <- df_clean$G / (df_clean$A + df_clean$C + df_clean$G + df_clean$T)

# Adding dinucleotide frequency (k-mers length of 2)
df_clean <- cbind(
  df_clean,
  as.data.frame(dinucleotideFrequency(df_clean$nuc2, as.prob = TRUE)))

# Converting DNAStringSet back to character format and removing DNAStringSet column
df_clean$nuc2_char <- as.character(df_clean$nuc2)
df_clean <- df_clean %>%
  select(-nuc2)

# ---- Build Random Forest Classifier ----

## ---- Classifier 1: predict gene (COI vs COII) ----

# Checking the number of sequences per gene after filtering
table(df_clean$gene)
# Maximum sample size is 62 for the COII gene, so will sample 15 genes (about 25% of the total) for the validation set

# Creating validation set
set.seed(123)

df_gene_val <- df_clean %>%
  group_by(gene) %>%
  sample_n(15)

# Confirming validation set
table(df_gene_val$gene)
# 15 samples from each

# Remaining set after 155 samples taken for validation set
df_remaining <- df_clean %>% filter(!id %in% df_gene_val$id)
df_remaining %>% count(gene)
# 12310 COI and 47 COII so will choose sample size of 40

# Creating training set
df_gene_train <- df_remaining %>%
  group_by(gene) %>%
  sample_n(40)

# Confirming training set
table(df_gene_train$gene)
# 40 samples from each

# Building a classifier
gene_classifier <- randomForest(
  x = df_gene_train[, c("Aprop", "Tprop", "Gprop")],
  y = as.factor(df_gene_train$gene),
  ntree = 200,
  importance = TRUE
)

gene_classifier
# COI vs COII performance was excellent, 100% accuracy with 0% OOB error and  no misclassification, indicating that these two genes are easily distinguishable using simple sequence composition features.

# Validation

pred_gene <- predict(gene_classifier, df_gene_val[, c("Aprop", "Tprop", "Gprop")])
table(observed = df_gene_val$gene, predicted = pred_gene)
# Works well

## ---- Classifier 2: predict subfamily (Danainae vs Heliconiinae) ----

# Checking the number of sequences per subfamily after filtering
table(df_clean$subfamily)
# Maximum sample size is 3714 for Danainae, so will sample 928 genes (about 25% of the total) for the validation set

# Creating validation set

set.seed(456)

df_subfamily_val <- df_clean %>%
  group_by(subfamily) %>%
  sample_n(928)

table(df_subfamily_val$subfamily)
#928 samples from each


df_subfamily_remaining <- df_clean %>% filter(!id %in% df_subfamily_val$id)
df_subfamily_remaining %>% count(subfamily)
# 2786 Danainae and 7745 Heliconiinae so will choose sample size of 2500

# Creating training set

df_subfamily_train <- df_subfamily_remaining %>%
  group_by(subfamily) %>%
  sample_n(2500)

table(df_subfamily_train$subfamily)
#2500 samples from each

# building a classifier

subfamily_classifier <- randomForest(
  x = df_subfamily_train[, c("Aprop", "Tprop", "Gprop")],
  y = as.factor(df_subfamily_train$subfamily),
  ntree = 200,
  importance = TRUE
)

subfamily_classifier
# The subfamily classifier performs well overall, with a low OOB error rate (7.36%). Most Danainae and Heliconiinae sequences are correctly classified, but subfamily separation is more difficult than gene classification due to greater overlap in nucleotide composition

# Validation
pred_subfamily <- predict(subfamily_classifier, df_subfamily_val[, c("Aprop", "Tprop", "Gprop")])

table(observed = df_subfamily_val$subfamily, predicted = pred_subfamily)
# Works pretty well

# ---- ROC Curves ----
# Testing to see how good both gene_classifier and subfamily_classifier are at telling 2 classes apart

## ---- ROC for gene classifier (COI vs COII) ----

pred_gene_prob <- predict(
  gene_classifier,
  df_gene_val[, c("Aprop", "Tprop", "Gprop")],
  type = "prob"
)

# Probability of one class - choosing COI as the control
roc_gene <- roc(
  response = df_gene_val$gene,
  predictor = pred_gene_prob[, "COI"]
)

# Plotting it
png("../figs/ROC_genes.png", width = 800, height = 600, res = 120)

plot(
  roc_gene,
  col = "blue",
  lwd = 2,
  main = "ROC Curve: Gene Classifier (COI vs COII)"
)

dev.off()
# Interpretation: this curve reaches the top-left corner with an AUC of 1.0, indicating perfect classification performance. This means the model separates COI and COII with 100% sensitivity and specificity across all thresholds

auc(roc_gene)
# Area under the curve = 1, meaning there is perfect classification
# The model can separate COI and COII with 100% accuracy across all thresholds

## ---- ROC for subfamily classifier (Danainae vs Heliconiinae) ----

pred_subfamily_prob <- predict(
  subfamily_classifier,
  df_subfamily_val[, c("Aprop", "Tprop", "Gprop")],
  type = "prob"
)

# Probability of one class - choosing Danainae as the control
roc_subfamily <- roc(
  response = df_subfamily_val$subfamily,
  predictor = pred_subfamily_prob[, "Danainae"]
)

# Plotting it
png("../figs/ROC_subfamilies.png", width = 800, height = 600, res = 120)

plot(
  roc_subfamily,
  col = "red",
  lwd = 2,
  main = "ROC Curve: Subfamily Classifier (Danainae vs Heliconiinae)"
)

dev.off()
# Interpretation: this curve shows strong classification performance for separating Danainae and Heliconiinae, with high sensitivity achieved at very high specificity. This indicates the subfamily classifier is highly accurate, though not perfectly separable like the gene classifier

auc(roc_subfamily)
# Area under the curve = 0.9781, meaning very strong classification
# The classifier can almost always distinguish between the two butterfly subfamilies

# Overlaying both ROC curves
png("../figs/ROC_comparison.png", width = 800, height = 600, res = 120)

plot(roc_gene, col = "blue", lwd = 2, main = "ROC Comparison")
plot(roc_subfamily, col = "red", lwd = 2, add = TRUE)

legend(
  "bottomright",
  legend = c("Gene (COI vs COII)", "Subfamily (Danainae vs Heliconiinae)"),
  col = c("blue", "red"),
  lwd = 2, 
  bty = "n"
)

dev.off()
# Interpretation: both classifiers perform very well, but the gene classifier (COI vs COII) shows perfect discrimination with an AUC of 1. The subfamily classifier (Danainae vs Heliconiinae) also performs strongly but with slightly lower accuracy, reflecting greater biological overlap between subfamilies than between genes

# ---- SVM Classifiers ----

# adding a second classifier to see which one works better

## ---- Gene Classifier (COI vs COII) ----

# Training SVM

gene_feat <- c("Aprop", "Tprop", "Gprop")

# Extracting feature matrix for gene training set and converting gene labels to factor for classification
x_gene_train <- df_gene_train[, gene_feat]
y_gene_train <- as.factor(df_gene_train$gene)

# Training linear SVM model to classify COI vs COII
svm_gene <- svm(
  x = x_gene_train,
  y = y_gene_train, 
  kernel = "linear",
  probability = TRUE,
  scale = TRUE
)

svm_gene
# Interpretation: the model required only 6 support vectors, indicating that the two gene classes are very well separated using the nucleotide proportion features. This confirms that the classification problem is very easy to solve with these features

# Validate SVM

x_gene_val <- df_gene_val[, gene_feat]
y_gene_val <- as.factor(df_gene_val$gene)

svm_gene_pred <- predict(svm_gene, x_gene_val)

# confusion matrix
table(observed = y_gene_val, predicted = svm_gene_pred)
# The model works well

# ROC and AUC for SVM for genes

# Predicting SVM class probabilities for the validation dataset (COI vs non-COI)
svm_gene_pred_prob <- attr(
  predict(svm_gene, x_gene_val, probability = TRUE),
  "probabilities"
)[, "COI"]

# Generating ROC curve for the SVM classifier using predicted probabilities
roc_gene_svm <- roc(response = y_gene_val, predictor = svm_gene_pred_prob)
auc(roc_gene_svm)
# Area under the curve = 1, meaning there is perfect classification
# The model can separate COI and COII with 100% accuracy across all thresholds

# Plotting the Random Forest ROC curve as a reference
png("../figs/ROC_vs_SVM_gene_comparison.png", width = 800, height = 600, res = 120)

plot(roc_gene, col = "blue", lwd = 2, main = "Gene ROC: RF vs SVM")

# Overlaying the SVM ROC curve for direct model comparison
plot(roc_gene_svm, col = "red", lwd = 2, add = TRUE)
legend("bottomright",
       legend = c("Random Forest", "SVM"),
       col = c("blue", "red"),
       lwd = 2, bty = "n")

dev.off()

# Interpretation: both the Random Forest and SVM models show identical ROC curves, indicating that they perform equally well at distinguishing between gene classes

## ---- Subfamily Classifier (Danainae vs Heliconiinae) ----

# Training SVM

subfamily_feat <- c("Aprop", "Tprop", "Gprop")

# Extracting feature matrix for subfamily training set and converting subfamily labels to factor for classification
x_subfamily_train <- df_subfamily_train[, subfamily_feat]
y_subfamily_train <- as.factor(df_subfamily_train$subfamily)

# Training linear SVM model to classify Danainae vs Heliconiinae
svm_subfamily <- svm(
  x = x_subfamily_train,
  y = y_subfamily_train, 
  kernel = "linear",
  probability = TRUE,
  scale = TRUE
)

svm_subfamily
# Interpretation: the model required a large number of support vectors (3669), indicating that the subfamily classes are more complex and less cleanly separated using the available features. This suggests that the classification task is more difficult and likely involves substantial overlap among subfamilies.

# Validate SVM

x_subfamily_val <- df_subfamily_val[, subfamily_feat]
y_subfamily_val <- as.factor(df_subfamily_val$subfamily)

svm_subfamily_pred <- predict(svm_subfamily, x_subfamily_val)

# confusion matrix
table(observed = y_subfamily_val, predicted = svm_subfamily_pred)
# The model performs with an accuracy of 68%, with substantial misclassification between Danainae and Heliconiinae indicating overlap between the two subfamilies

# ROC and AUC for SVM for subfamilies

# Predicting SVM class probabilities for the validation dataset (Dinainae vs Heliconiinae)
svm_subfamily_pred_prob <- attr(
  predict(svm_subfamily, x_subfamily_val, probability = TRUE),
  "probabilities"
)[, "Danainae"]

# Generating ROC curve for the SVM classifier using predicted probabilities
roc_subfamily_svm <- roc(response = y_subfamily_val, predictor = svm_subfamily_pred_prob)
auc(roc_subfamily_svm)
# Area under the curve = 0.737, indicating moderate classification performance
# The classifier can distinguish between the two butterfly subfamilies better than random, but with noticeable overlap
png("../figs/ROC_vs_SVM_sunfamily_comparison.png", width = 800, height = 600, res = 120)

# Plotting the Random Forest ROC curve as a reference
plot(roc_subfamily, col = "blue", lwd = 2, main = "Subfamily ROC: RF vs SVM")

# Overlaying the SVM ROC curve for direct model comparison
plot(roc_subfamily_svm, col = "red", lwd = 2, add = TRUE)
legend("bottomright",
       legend = c("Random Forest", "SVM"),
       col = c("blue", "red"),
       lwd = 2, bty = "n")

dev.off()

# Interpretation: The Random Forest model substantially outperforms the SVM at the subfamily level, showing much higher sensitivity across most specificity values, while the SVM displays only moderate discriminatory ability. This indicates that Random Forest captures subfamily-level sequence patterns far more effectively than the linear SVM
