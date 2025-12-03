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
library(rentrez)
#install.packages("pROC")
library(pROC)
#install.packages("e1071")
library(e1071)

# ---- Load butterfly data from BOLD ----

dfDanainae <- read_tsv("../data/Danainae.tsv")
dfHeliconiinae <- read_tsv("../data/Heliconiinae.tsv")

# ---- Checking and filtering data ----

summary(dfDanainae)
summary(dfHeliconiinae)

names(dfDanainae)
names(dfHeliconiinae)

table(dfDanainae$marker_code)
#6112 COI-5P and 8 COII
table(dfHeliconiinae$marker_code)
#10127 COI-5P and 31 COII

#histogram of sequence lengths for COI-5P
Len_Dan_COI <- nchar(dfDanainae$nuc[dfDanainae$marker_code == "COI-5P"])
Len_Hel_COI <- nchar(dfHeliconiinae$nuc[dfHeliconiinae$marker_code == "COI-5P"])

hist(Len_Dan_COI,
     xlab = "Sequence Length",
     ylab = "Frequency",
     main = "COI-5P Sequence Lengths: Danainae vs Heliconiinae",
     col = rgb(0, 0, 1, 0.4),   
     border = "blue"
)

# add Heliconiinae on top
hist(Len_Hel_COI,
     col = rgb(1, 0, 0, 0.4),   
     border = "red",
     add = TRUE
)

legend(
  x = 950,
  y = 4000, 
  legend = c("Danainae", "Heliconiinae"),
  fill = c(rgb(0,0,1,0.4), rgb(1,0,0,0.4)),
  border = c("blue", "red"),
  bty = "n", 
  cex = 0.9
  )


#histogram of sequence lengths for COII
Len_Dan_COII <- nchar(dfDanainae$nuc[dfDanainae$marker_code == "COII"])
Len_Hel_COII <- nchar(dfHeliconiinae$nuc[dfHeliconiinae$marker_code == "COII"])

par(mar = c(5,4,10,2))

hist(Len_Dan_COII,
     xlab = "Sequence Length",
     ylab = "Frequency",
     main = "COII Sequence Lengths: Danainae vs Heliconiinae",
     col = rgb(0, 0, 1, 0.4),   
     border = "blue"
)

# add Heliconiinae on top
hist(Len_Hel_COII,
     col = rgb(1, 0, 0, 0.4),   
     border = "red",
     add = TRUE
)

legend(
  x = 671.1,
  y = 6.5,
  legend = c("Danainae", "Heliconiinae"),
  fill = c(rgb(0,0,1,0.4), rgb(1,0,0,0.4)),
  border = c("blue", "red"),
  bty = "n",      
  cex = 0.9       
)

rm(Len_Dan_COI, Len_Dan_COII, Len_Hel_COI, Len_Hel_COII)

# ---- Loading data from NCBI ----
# ---- Download COI sequences for each species of butterfly ----

Danainae_COI_search <- entrez_search(
  db = "nucleotide", 
  term = "Danaus plexippus[Organism] AND COI[Gene]", 
  retmax = 1000, 
  use_history = TRUE
)

Danainae_COI_fasta <- entrez_fetch(
  db = "nucleotide", 
  web_history = Danainae_COI_search$web_history,
  rettype = "fasta"
)

write(Danainae_COI_fasta, file = "../data/Danainae_COI_NCBI.fasta")

Heliconiinae_COI_search <- entrez_search(
  db = "nucleotide",
  term = "Heliconius erato[Organism] AND COI[Gene]", 
  retmax = 1000, 
  use_history = TRUE
)

Heliconiinae_COI_fasta <- entrez_fetch(
  db = "nucleotide", 
  web_history = Heliconiinae_COI_search$web_history, 
  rettype = "fasta"
)

write(Heliconiinae_COI_fasta, file = "../data/Heliconiinae_COI_NCBI.fasta")

# ---- Download COII sequences for each species of butterfly

Danainae_COII_search <- entrez_search(
  db = "nucleotide", 
  term = "Danaus plexippus[Organism] AND COII[Gene]", 
  retmax = 1000,
  use_history = TRUE
)

Danainae_COII_fasta <- entrez_fetch(
  db = "nucleotide", 
  web_history = Danainae_COII_search$web_history,
  rettype = "fasta"
)

write(Danainae_COII_fasta, file = "../data/Danainae_COII_NCBI.fasta")

Heliconiinae_COII_search <- entrez_search(
  db = "nucleotide",
  term = "Heliconius erato[Organism] AND COII[Gene]",
  retmax = 1000, 
  use_history = TRUE
)


Heliconiinae_COII_fasta <- entrez_fetch(
  db = "nucleotide",
  web_history = Heliconiinae_COII_search$web_history,
  rettype = "fasta"
)

write(Heliconiinae_COII_fasta, file = "../data/Heliconiinae_COII_NCBI.fasta")

# ---- Load FASTA files into R ----

Danainae_COI <- readDNAStringSet("../data/Danainae_COI_NCBI.fasta")
Heliconiinae_COI <- readDNAStringSet("../data/Heliconiinae_COI_NCBI.fasta")

Danainae_COII <- readDNAStringSet("../data/Danainae_COII_NCBI.fasta")
Heliconiinae_COII <- readDNAStringSet("../data/Heliconiinae_COII_NCBI.fasta")

rm(Danainae_COI_search, Danainae_COII_search, Heliconiinae_COI_search, Heliconiinae_COII_search)
rm(Danainae_COI_fasta, Danainae_COII_fasta, Heliconiinae_COI_fasta, Heliconiinae_COII_fasta)

# ---- Converting each FASTA to a tibble ----

fasta_to_df <- function(x, species, gene) {
  tibble(
    id = names(x),
    nuc = as.character(x),
    species = species,
    gene = gene
  )
}

df_Dan_COI <- fasta_to_df(Danainae_COI, "Danainae", "COI")
df_Hel_COI <- fasta_to_df(Heliconiinae_COI, "Heliconiinae", "COI")

df_Dan_COII <- fasta_to_df(Danainae_COII, "Danainae", "COII")
df_Hel_COII <- fasta_to_df(Heliconiinae_COII, "Heliconiinae", "COII")

rm(Danainae_COI, Danainae_COII, Heliconiinae_COI, Heliconiinae_COII)

# ---- Standardize the dataframes ----

dfDan_BOLD <- dfDanainae %>%
  filter(!is.na(nuc)) %>%
  transmute(
    id = processid,
    nuc = nuc,
    species = "Danainae",
    gene = marker_code
  )

dfHel_BOLD <- dfHeliconiinae %>%
  filter(!is.na(nuc)) %>%
  transmute(
    id = processid,
    nuc = nuc,
    species = "Heliconiinae", 
    gene = marker_code
  )

# rename COI-5P to COI since that is what it's called in BOLD

dfDan_BOLD$gene[dfDan_BOLD$gene == "COI-5P"] <- "COI"
dfHel_BOLD$gene[dfHel_BOLD$gene == "COI-5P"] <- "COI"

# combine all the datasets into one dataframe

df_all <- bind_rows(
  dfDan_BOLD,
  dfHel_BOLD,
  df_Dan_COI,
  df_Hel_COI,
  df_Dan_COII,
  df_Hel_COII
)

rm(dfDan_BOLD, dfHel_BOLD, df_Dan_COI, df_Hel_COI, df_Dan_COII, df_Hel_COII)

# ---- Clean and filter sequences ----
# remove sequences with gaps or N's at the ends

df_all <- df_all %>%
  mutate(nuc2 = str_remove(nuc, "^[-N]+")) %>%
  mutate(nuc2 = str_remove(nuc2, "[-N]+$")) %>%
  mutate(nuc2 = str_remove_all(nuc2, "-+"))
# remove sequences with >5% Ns

df_all <- df_all %>%
  filter(str_count(nuc2, "N") <= 0.05 * nchar(nuc))

# compare old and new sequence columns
df_seq_compare <- cbind(df_all$nuc, df_all$nuc2)

rm(df_seq_compare)

# build quartiles to restrict COI sequence lengths

q1 <- quantile(nchar(df_all$nuc2[df_all$gene == "COI"]), 0.25)
q3 <- quantile(nchar(df_all$nuc2[df_all$gene == "COI"]), 0.75)

df_clean <- df_all %>%
  filter(
    (gene == "COI" & nchar(nuc2) >= q1 & nchar(nuc2) <= q3 | gene == "COII")
  )

# convert sequences to DNAStringSet

df_clean <- as.data.frame(df_clean)
df_clean$nuc2 <- DNAStringSet(df_clean$nuc2)

# extract nucleotides and k-mer features

df_clean <- cbind(
  df_clean,
  as.data.frame(letterFrequency(df_clean$nuc2, letters = c("A", "C", "G", "T")))
)

#proportions
df_clean$Aprop <- df_clean$A / (df_clean$A + df_clean$C + df_clean$G + df_clean$T)
df_clean$Tprop <- df_clean$T / (df_clean$A + df_clean$C + df_clean$G + df_clean$T)
df_clean$Gprop <- df_clean$G / (df_clean$A + df_clean$C + df_clean$G + df_clean$T)

# adding dinucleotide frequency (k-mers length of 2)
df_clean <- cbind(
  df_clean,
  as.data.frame(dinucleotideFrequency(df_clean$nuc2, as.prob = TRUE)))

df_clean$nuc2_char <- as.character(df_clean$nuc2)
df_clean <- df_clean %>%
  select(-nuc2)

# ---- Build Random Forest Classifier ----

table(df_clean$gene)
# maximum sample size is 620 for the COII gene, so will sample 155 genes (about 25% of the total) for the validation set

# Classifier 1: predict gene (COI vs COII)

# Creating validation set

set.seed(123)

df_gene_val <- df_clean %>%
  group_by(gene) %>%
  sample_n(155)

table(df_gene_val$gene)
#155 samples from each


df_remaining <- df_clean %>% filter(!id %in% df_gene_val$id)
df_remaining %>% count(gene)
#465 COII and 9003 COI so will choose sample size of 400

# Creating training set

df_gene_train <- df_remaining %>%
  group_by(gene) %>%
  sample_n(400)

table(df_gene_train$gene)
#400 samples from each

# building a classifier

gene_classifier <- randomForest(
  x = df_gene_train[, c("Aprop", "Tprop", "Gprop")],
  y = as.factor(df_gene_train$gene),
  ntree = 200,
  importance = TRUE
)

gene_classifier
# good performance, error rate of 0.12%, misclassified 1

# validate

pred_gene <- predict(gene_classifier, df_gene_val[, c("Aprop", "Tprop", "Gprop")])
table(observed = df_gene_val$gene, predicted = pred_gene)
# works well

# Classifier 2: predict species (Danainae vs Heliconiinae)

unique(df_clean$species)

table(df_clean$species)
# maximum sample size is 3025 for Danainae, so will sample 755 genes (about 25% of the total) for the validation set

# Classifier 1: predict species (Danainae vs Heliconiinae)

# Creating validation set

set.seed(456)

df_species_val <- df_clean %>%
  group_by(species) %>%
  sample_n(755)

table(df_species_val$species)
#755 samples from each


df_species_remaining <- df_clean %>% filter(!id %in% df_species_val$id)
df_species_remaining %>% count(species)
#2270 Danainae and 5998 Heliconiinae so will choose sample size of 2000

# Creating training set

df_species_train <- df_species_remaining %>%
  group_by(species) %>%
  sample_n(2000)

table(df_species_train$species)
#2000 samples from each

# building a classifier

species_classifier <- randomForest(
  x = df_species_train[, c("Aprop", "Tprop", "Gprop")],
  y = as.factor(df_species_train$species),
  ntree = 200,
  importance = TRUE
)

species_classifier
# good performance, error rate of 7.15%

# validate

pred_species <- predict(species_classifier, df_species_val[, c("Aprop", "Tprop", "Gprop")])

table(observed = df_species_val$species, predicted = pred_species)
# works pretty good

# ---- ROC Curves ----
# test to see how good both gene_classifier and species_classifier are at telling 2 classes apart

# ROC for gene classifier (COI vs COII)

pred_gene_prob <- predict(
  gene_classifier,
  df_gene_val[, c("Aprop", "Tprop", "Gprop")],
  type = "prob"
)

# probability of one class - choosing COI as the control

roc_gene <- roc(
  response = df_gene_val$gene,
  predictor = pred_gene_prob[, "COI"]
)

#plot it
plot(
  roc_gene,
  col = "blue",
  lwd = 2,
  main = "ROC Curve: Gene Classifier (COI vs COII)"
)

auc(roc_gene)
#Area under the curve = 1, meaning there is perfect classification
# The model can separate COI and COII with 100% accuracy across all thresholds

# ROC for species classifier (Danainae vs Heliconiinae)

pred_species_prob <- predict(
  species_classifier,
  df_species_val[, c("Aprop", "Tprop", "Gprop")],
  type = "prob"
)

# probability of one class - choosing Danainae as the control

roc_species <- roc(
  response = df_species_val$species,
  predictor = pred_species_prob[, "Danainae"]
)

#plot it
plot(
  roc_species,
  col = "red",
  lwd = 2,
  main = "ROC Curve: Species Classifier (Danainae vs Heliconiinae)"
)

auc(roc_species)
#Area under the curve = 0.9699, meaning very strong classification
# The classifier can almost always distinguish between the two butterfly subfamilies

# Overlaying both ROC curves

plot(roc_gene, col = "blue", lwd = 2, main = "ROC Comparison")
plot(roc_species, col = "red", lwd = 2, add = TRUE)

legend(
  "bottomright",
  legend = c("Gene (COI vs COII)", "Species (Danainae vs Heliconiinae)"),
  col = c("blue", "red"),
  lwd = 2, 
  bty = "n"
)

# ---- SVM Classifiers ----

# adding a second classifier to see which one works better

# ---- Gene Classifier (COI vs COII) ----

# Training SVM

gene_feat <- c("Aprop", "Tprop", "Gprop")

x_gene_train <- df_gene_train[, gene_feat]
y_gene_train <- as.factor(df_gene_train$gene)

svm_gene <- svm(
  x = x_gene_train,
  y = y_gene_train, 
  kernel = "linear",
  probability = TRUE,
  scale = TRUE
)

svm_gene

# Validate SVM

x_gene_val <- df_gene_val[, gene_feat]
y_gene_val <- as.factor(df_gene_val$gene)

svm_gene_pred <- predict(svm_gene, x_gene_val)

# confusion matrix
table(observed = y_gene_val, predicted = svm_gene_pred)

# ROC and AUC for SVM 

svm_gene_pred_prob <- attr(
  predict(svm_gene, x_gene_val, probability = TRUE),
  "probabilities"
)[, "COI"]

roc_gene_svm <- roc(response = y_gene_val, predictor = svm_gene_pred_prob)
auc(roc_gene_svm)

plot(roc_gene, col = "blue", lwd = 2, main = "Gene ROC: RF vs SVM")
plot(roc_gene_svm, col = "red", lwd = 2, add = TRUE)
legend("bottomright",
       legend = c("Random Forest", "SVM"),
       col = c("blue", "red"),
       lwd = 2, bty = "n")

# ---- Species Classifier (Danainae vs Heliconiinae) ----

# Training SVM

species_feat <- c("Aprop", "Tprop", "Gprop")

x_species_train <- df_species_train[, species_feat]
y_species_train <- as.factor(df_species_train$species)

svm_species <- svm(
  x = x_species_train,
  y = y_species_train, 
  kernel = "linear",
  probability = TRUE,
  scale = TRUE
)

svm_species

# Validate SVM

x_species_val <- df_species_val[, species_feat]
y_species_val <- as.factor(df_species_val$species)

svm_species_pred <- predict(svm_species, x_species_val)

# confusion matrix
table(observed = y_species_val, predicted = svm_species_pred)

# ROC and AUC for SVM 

svm_species_pred_prob <- attr(
  predict(svm_species, x_species_val, probability = TRUE),
  "probabilities"
)[, "Danainae"]

roc_species_svm <- roc(response = y_species_val, predictor = svm_species_pred_prob)
auc(roc_species_svm)

plot(roc_species, col = "blue", lwd = 2, main = "Species ROC: RF vs SVM")
plot(roc_species_svm, col = "red", lwd = 2, add = TRUE)
legend("bottomright",
       legend = c("Random Forest", "SVM"),
       col = c("blue", "red"),
       lwd = 2, bty = "n")
