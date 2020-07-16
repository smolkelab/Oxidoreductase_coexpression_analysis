# Import A. belladonna RNA sequencing dataset, downloaded from MSU Medicinal Plant Genomics Resource
logFPKMs <- read.csv(file="/path/to/dataset/aba.matrix.FPKM.vf.082511.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
logFPKMs$X <- NULL
colnames(logFPKMs)[3:13] <- unlist(lapply(logFPKMs[1,c(3:13)], as.character), use.names=FALSE)
logFPKMs <- logFPKMs[-1,]
logFPKMvals <- data.matrix(logFPKMs[,3:13])
rownames(logFPKMvals) <- paste(logFPKMs$Locus.ID, logFPKMs$Unified.Functional.Annotation, logFPKMs$PFAM, sep=">")
logFPKMvals <- 2^logFPKMvals #Using raw reads instead of log-vals

# Identify initial list of oxidoreductase candidates based on PFAM and functional annotations
candidates <- logFPKMvals[grep("PF00106|PF13561|PF08659|PF08240|PF00107|PF00248|PF00465|PF13685|PF13823|PF13602|PF16884|PF00248|alcohol dehydrogenase|aldehyde reductase|short chain|aldo/keto|littorine|hyoscyamine|putrescine|tropinone|tropine", rownames(logFPKMvals), ignore.case=TRUE),]

# Using CYP80F1 littorine monooxygenase as bait gene
baitdata <- colMeans(candidates[grep("littorine", rownames(candidates), ignore.case = TRUE),]) 
model <- apply(candidates, 1, function(x) summary(lm(baitdata~x))$coefficients[,4])
CYP80F1_pvals <- data.frame(p=sapply(model, function(x) x[2]))
CYP80F1_pvals <- na.omit(CYP80F1_pvals[order(CYP80F1_pvals$p), , drop=FALSE])

# Using hyoscyamine 6b-hydroxylase/oxygenase as bait gene
baitdata <- colMeans(candidates[grep("hyoscyamine", rownames(candidates), ignore.case = TRUE),]) 
model <- apply(candidates, 1, function(x) summary(lm(baitdata~x))$coefficients[,4])
H6H_pvals <- data.frame(p=sapply(model, function(x) x[2]))
H6H_pvals <- na.omit(H6H_pvals[order(H6H_pvals$p), , drop=FALSE])

# Take hits for each of the bait genes and trim duplicates, then compute product of the p-values
top_hits <- data.frame(ID=unique(c(rownames(CYP80F1_pvals), rownames(H6H_pvals)))) #compile lists
top_hits$ID <- sub('>.*', '', top_hits$ID) #remove everything but locus IDs from rownames
top_hits$combined_p <- sapply(seq(1:length(top_hits$ID)), function(x) 
  log10(CYP80F1_pvals$p[grep(top_hits$ID[x], rownames(CYP80F1_pvals), ignore.case = TRUE)])
  + log10(H6H_pvals$p[grep(top_hits$ID[x], rownames(H6H_pvals), ignore.case = TRUE)]))
top_hits <- top_hits[order(top_hits$combined_p), , drop=FALSE] #order by combined p value
top_hits <- top_hits[top_hits$combined_p < -1.3,] #drop any with combined P < 0.05 (log10 < -1.3)

# Generate subset of original log2FPKM values with the top_hits candidates
top_hits_FPKMs <- data.matrix(logFPKMs[sapply(top_hits$ID, function(x)
  grep(x, logFPKMs$Locus.ID, ignore.case = TRUE)), 3:13, drop=FALSE])
rownames(top_hits_FPKMs) <- top_hits$ID
top_hits_FPKMs <- 2^top_hits_FPKMs

# OPTIONAL: Normalize top_hits_FPKMs by highest expression level for that gene
norm_FPKMs <- t(apply(top_hits_FPKMs, 1, function(x)(x-min(x))/(max(x)-min(x))))

# Generate heatmap for top hits, non-normalized
library(gplots)
library(RColorBrewer)
heatmap.2(top_hits_FPKMs, key=TRUE, col=colorRampPalette(c('red', 'black', 'green')), Colv=FALSE, dendrogram="row",
          scale="row", margins=c(10, 40), trace="none", sepwidth=c(0,0), density.info = 'none', 
          key.title = NA, colsep=1:ncol(top_hits_FPKMs), rowsep=1:nrow(top_hits_FPKMs), keysize = 1)

# Generate heatmap for top hits, normalized
heatmap.2(norm_FPKMs, key=TRUE, col=colorRampPalette(c('red', 'black', 'green')), Colv=FALSE, dendrogram="row",
          margins=c(10, 30), tracecol="white", trace="none", sepwidth=c(0,0), density.info = 'none', key.title = NA,
          colsep=1:ncol(norm_FPKMs), rowsep=1:nrow(norm_FPKMs), keysize = 1)