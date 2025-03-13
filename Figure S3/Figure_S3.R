# Essential libraries
library(readxl)        # Reading Excel files
library(ggplot2)       # For creating plots
library(dplyr)         # Data manipulation
library(ggpubr)        # Arranging multiple ggplot objects
library(multcompView)  # Extracting Tukey's test letters
library(vegan)         # PERMANOVA analysis and distance calculations


# Set working directory and read data
setwd("/")

##
##
##qPCR boxplots
##
##

qPCR <- read_excel("./Users/guillaumechesneau/Documents/Post_doc/Postdoc/Scripts/Chesneau_et_al_2025/Figure S3/FigS3.xlsx", sheet = "All")  # Read qPCR data

# Define color palettes
Calaaas3 <- c("soil"="black", "COL0"="#F3766E", "Medium"="grey", "Water"="grey", "Inoculum"="grey")
calas3 <- c("soil"="black", "COL0"="black", "Medium"="black", "Water"="black", "Inoculum"="black")

# Generate boxplots for qPCR results
a <- ggplot(qPCR, aes(x=factor(Type, level=c("soil", "COL0", "Medium", "Water", "Inoculum")), 
                      y=Cq_16S, fill=Type, color=Type)) +
  geom_boxplot(position=position_dodge(width=0.75), width=0.6, outlier.shape=NA) +
  geom_point(size=1, alpha=0.4, color="black") +
  ylim(10, 35) +
  scale_fill_manual(values=Calaaas3) +
  scale_color_manual(values=calas3) +
  theme_bw()

b <- ggplot(qPCR, aes(x=factor(Type, level=c("soil", "COL0", "Medium", "Water", "Inoculum")), 
                      y=Cq_ITS, fill=Type, color=Type)) +
  geom_boxplot(position=position_dodge(width=0.75), width=0.6, outlier.shape=NA) +
  geom_point(size=1, alpha=0.4, color="black") +
  ylim(10, 35) +
  scale_fill_manual(values=Calaaas3) +
  scale_color_manual(values=calas3) +
  theme_bw()

ggarrange(a, b, ncol=2)  # Arrange plots side by side

# ANOVA and Tukey's HSD test
anova_result <- aov(Cq_ITS ~ Type, data=qPCR)
summary(anova_result)

tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Tukey letters for grouping
tukey_letters <- multcompLetters4(anova_result, tukey_result)
print(tukey_letters)


##
##
##PCAs of exudates/soil wash/medium
##
##


# PCA analysis
log10_pub2 <- read_excel("./FigS3.xlsx", sheet="log10_pub")  # Read PCA input data
log10_pub2[, 5:41] <- data.frame(lapply(log10_pub2[, 5:41], as.numeric))  # Ensure numeric values

# Replace NA with column means
log10_pub2[, 5:41] <- log10_pub2[, 5:41] %>% mutate(across(everything(), ~ replace_na(., mean(., na.rm=TRUE))))

# Remove redundant columns (adjust column index as needed)
pub2 <- data.frame(log10_pub2[, -40])

# PCA computation
pca <- prcomp(pub2[, 5:40], center=TRUE, scale.=TRUE)
rownames(pub2) <- pub2$SAMPLE_ID
summary(pca)

# PCA visualization
autoplot(pca, data=pub2, colour='Genotype', frame=TRUE, size=5,
         label=FALSE, label.size=3,
         loadings=FALSE, loadings.colour='grey',
         loadings.label=FALSE, loadings.label.size=3) +
  theme_bw()

# PERMANOVA analysis
distance_matrix <- dist(pub2[, 5:40], method="euclidean")
permanova_result <- adonis2(distance_matrix ~ Genotype, data=pub2, method="euclidean", permutations=999)
print(permanova_result)
