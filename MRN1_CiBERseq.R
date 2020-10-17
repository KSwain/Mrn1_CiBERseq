## Load in gene and guide annotation files
if (!file.exists("SGD_features.tab")) {
  sgd <- download.file('https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab', destfile="SGD_features.tab")
}

sgd <- read.delim("SGD_features.tab", header=FALSE, quote="",
                  col.names=c("sgdid", "type", "qual", "name", "gene", "alias",
                              "parent", "sgdid2", "chrom", "start", "end",
                              "strand", "genpos", "cver", "sver", "desc"))
head(sgd)

## load in barcode to guide dataframe
options(stringsAsFactors=FALSE)
if (!exists("grna.assign.barcode.grna.good")) {
  grna.assign.barcode.grna.good <- read.delim("/mnt/ingolialab/kswain/NIPD003/grna-assign-barcode-grna-good.txt",
                                              stringsAsFactors=FALSE)
}

grna.assign.barcode.grna.good$guide <- as.character(grna.assign.barcode.grna.good$guide)
grna.assign.barcode.grna.good$guide[grna.assign.barcode.grna.good$guide == "No_gRNA"] <- paste("No_gRNA", seq(1:788), sep=" ")
head(grna.assign.barcode.grna.good)

## load in guide target df
options(stringsAsFactors=FALSE)
if (!exists("guide.good.targets")) {
  guide.good.targets <- read.delim("/mnt/ingolialab/kswain/NIPD003/guide.good.targets_plus_empty.txt",
                                   stringsAsFactors=FALSE)
}

head(guide.good.targets)
head(grna.assign.barcode.grna.good)

library(DESeq2)

#load in count data 
# Load in raw count data and relevant gRNA barcode assignment and gRNA assignment of target files
options(stringsAsFactors=FALSE)
if (!exists("mhalo")) {
  mhalo <- read.delim("~/NIKS022/longname-txt/mhalo_seq1_counts.txt",
                     stringsAsFactors=FALSE)
}
head(mhalo)
dim(mhalo)
halo_complete <- mhalo[complete.cases(mhalo[,14]),]
dim(halo_complete)
write.csv(halo_complete, "~/MRN1/complete_guidebc_tcaga.csv")


options(stringsAsFactors=FALSE)
if (!exists("nterm")) {
  nterm <- read.delim("~/NIKS022/longname-txt/mfrag_seq1_counts.txt",
                      stringsAsFactors=FALSE)
}
head(nterm)  
dim(nterm)
nterm_complete <- nterm[complete.cases(nterm[,14]),]
dim(nterm_complete)
write.csv(nterm_complete, "~/MRN1/complete_guidebc_ccttg.csv")
postn_GO <- postn_GO[complete.cases(postn_GO[,8]),]

## make colnames more legible: 
colnames(nterm) <- sub("X.mnt.ingolialab.kswain.NIKS022.", "", colnames(nterm))
colnames(nterm) <- sub(".count.txt", "", colnames(nterm))
colnames(nterm) <- sub("mfrag", "", colnames(nterm))
colnames(nterm) <- sub("L001", "", colnames(nterm))
colnames(nterm) <- sub("S[0-9]*", "", colnames(nterm))
colnames(nterm) <- sub("___R1_001", "", colnames(nterm))
head(nterm)

colnames(mhalo) <- sub("X.mnt.ingolialab.kswain.NIKS022.", "", colnames(mhalo))
colnames(mhalo) <- sub(".count.txt", "", colnames(mhalo))
colnames(mhalo) <- sub("mhalo", "", colnames(mhalo))
colnames(mhalo) <- sub("L001", "", colnames(mhalo))
colnames(mhalo) <- sub("S[0-9]*", "", colnames(mhalo))
colnames(mhalo) <- sub("___R1_001", "", colnames(mhalo))
head(mhalo)


## Aggregate by guide: match barcodes to guides and then add together
mhalo$guide <- grna.assign.barcode.grna.good[match(mhalo$barcode, grna.assign.barcode.grna.good$barcode), "guide"]
head(mhalo)
mhalo$Undetermined = NULL
nterm$guide <- grna.assign.barcode.grna.good[match(nterm$barcode, grna.assign.barcode.grna.good$barcode), "guide"]
nterm$Undetermined = NULL
head(nterm)

## aggregate barcodes by guide 

halo_agg <- aggregate(mhalo[, c("LePoD", "LePoR", "LePrD", "LePrR", "LeReD", "LeReR",
                               "RiPoD", "RiPoR", "RiPrD", "RiPrR", "RiReD", "RiReR")],
                      by=list(mhalo$guide), FUN=sum)
head(halo_agg)
dim(halo_agg)
names(halo_agg)[1] <- "guide"



nterm_agg <- aggregate(nterm[, c("LePoD", "LePoR", "LePrD", "LePrR", "LeReD", "LeReR",
                                "RiPoD", "RiPoR", "RiPrD", "RiPrR", "RiReD", "RiReR")],
                      by=list(nterm$guide), FUN=sum)
head(nterm_agg)
dim(nterm_agg)
names(nterm_agg)[1] <- "guide"


## take just the columns we want --> RNA only 

halo_rna <- halo_agg[, c("guide", "LePoR", "LePrR", "LeReR", "RiPoR", "RiPrR", "RiReR")]
head(halo_rna)

nterm_rna <- nterm_agg[, c("guide", "LePoR", "LePrR", "LeReR", "RiPoR", "RiPrR", "RiReR")]
head(nterm_rna)

## put RNA counts into one df

rna_counts <- halo_rna
rna_counts$LePoRh <- rna_counts$LePoR
rna_counts$LePrRh <- rna_counts$LePrR
rna_counts$RiPoRh <- rna_counts$RiPoR
rna_counts$RiPrRh <- rna_counts$RiPrR
head(rna_counts)

rna_counts$LePoR = NULL
rna_counts$LePrR = NULL
rna_counts$LeReR = NULL
rna_counts$RiPoR = NULL
rna_counts$RiPrR = NULL
rna_counts$RiReR = NULL
head(rna_counts)

## match in Mrn1 Nterm column data
head(nterm_rna)
rna_counts$LePoNt <- nterm_rna[match(rna_counts$guide, nterm_rna$guide), "LePoR"]
names(rna_counts) <- c("guide", "LePoHalo", "LePrHalo", "RiPoHalo", "RiPrHalo", "LePoNt")
rna_counts$LePrNt <- nterm_rna[match(rna_counts$guide, nterm_rna$guide), "LePrR"]
rna_counts$RiPoNt <- nterm_rna[match(rna_counts$guide, nterm_rna$guide), "RiPoR"]
rna_counts$RiPrNt <- nterm_rna[match(rna_counts$guide, nterm_rna$guide), "RiPrR"]
head(rna_counts)

## set NAs to 0
rna_counts[is.na(rna_counts)] = 0
cor(rna_counts$LePrHalo, rna_counts$RiPrHalo)
cor(rna_counts$LePrNt, rna_counts$RiPrNt)
cor(rna_counts$LePrHalo, rna_counts$LePrNt)

plot(log10(rna_counts$LePrHalo), log10(rna_counts$RiPrHalo))
plot(log10(rna_counts$LePrNt), log10(rna_counts$RiPrNt))
plot(log10(rna_counts$LePrHalo), log10(rna_counts$LePrNt))

## want only genes with 32 or more counts for PreRNA samples

dim(rna_counts)
sig_counts <- filter(rna_counts,
                     rna_counts$LePrHalo > 31 & rna_counts$RiPrHalo > 31 & rna_counts$LePrNt > 31 & rna_counts$RiPrNt > 31)
dim(sig_counts)
head(sig_counts)
cor(sig_counts$LePrHalo, sig_counts$RiPrHalo)
cor(sig_counts$RiPrNt, sig_counts$LePrNt)

## plot count data:
png(file="~/NIKS022/DESeq2_analysis/Halo_precounts.png", width=1000,height=1000,res=144)
ggplot(sig_counts, aes(log10(LePrHalo), log10(RiPrHalo))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("Pre-induction counts: Halo control replicates") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(Halo rep1 counts)", y="Log10(Halo rep2 counts)")
dev.off()

png(file="~/NIKS022/DESeq2_analysis/Halo_postcounts.png", width=1000,height=1000,res=144)
ggplot(sig_counts, aes(log10(LePoHalo), log10(RiPoHalo))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("Post-induction counts: Halo control replicates") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(Halo rep1 counts)", y="Log10(Halo rep2 counts)")
dev.off()

png(file="~/NIKS022/DESeq2_analysis/Mrn1_precounts.png", width=1000,height=1000,res=144)
ggplot(sig_counts, aes(log10(LePrNt), log10(RiPrNt))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("Pre-induction counts: Mrn1(1-200) replicates") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(Mrn1(1-200) rep1 counts)", y="Log10(Mrn1(1-200) rep2 counts)")
dev.off()

png(file="~/NIKS022/DESeq2_analysis/Mrn1_postcounts.png", width=1000,height=1000,res=144)
ggplot(sig_counts, aes(log10(LePoNt), log10(RiPoNt))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("Post-induction counts: Mrn1(1-200) replicates") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(Mrn1(1-200) rep1 counts)", y="Log10(Mrn1(1-200) rep2 counts)")
dev.off()


## set up DESeq analysis
counts = sig_counts
head(counts)
row.names(counts) = counts$guide
head(counts)
counts$guide = NULL
head(counts)

conditions = data.frame(row.names=c("LePoHalo",  "LePrHalo", "RiPoHalo", "RiPrHalo", "LePoNt", "LePrNt", "RiPoNt", "RiPrNt"),
                        Atc=factor(c("Post", "Pre", "Post", "Pre", "Post", "Pre", "Post", "Pre"), levels=c("Pre", "Post")),
                        geno=factor(c("Halo", "Halo", "Halo", "Halo", "Nterm", "Nterm", "Nterm", "Nterm"), levels=c("Halo", "Nterm")),
                        turb=c("Left", "Left", "Right", "Right", "Left", "Left", "Right", "Right"))
conditions

dds <- DESeqDataSetFromMatrix(countData = counts, colData = conditions, design = ~ Atc * geno + turb)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

## set up pre-only matrix to calculate dispersions 
head(counts)
pre_counts <- counts 
pre_counts$LePoHalo = NULL
pre_counts$RiPoHalo = NULL
pre_counts$LePoNt = NULL
pre_counts$RiPoNt = NULL
head(pre_counts)

cond.disp = data.frame(row.names=c("LePrHalo", "RiPrHalo", "LePrNt", "RiPrNt"),
                       geno=c("Halo", "Halo", "Nterm", "Nterm"))
cond.disp
dds.disp <- DESeqDataSetFromMatrix(countData = pre_counts, colData = cond.disp, design = ~ geno)
dds.disp <- estimateSizeFactors(dds.disp)
sizeFactors(dds.disp)
dds.disp <- estimateDispersions(dds.disp)
plotDispEsts(dds.disp, ylim=c(1e-4,1))

# Transfer dispersions from "dds.disp" over to full "dds"
dispersions(dds) <- dispersions(dds.disp)
dds <- nbinomWaldTest(dds, betaPrior=FALSE)

resultsNames(dds)
## run one or the other of the following:
#res <- results(dds, name="AtcPost.genoNterm")
#res <- results(dds, name="Atc_Post_vs_Pre")
res= as.data.frame(res)
head(res)

# save df based on which results ran
#post_nterm <- res
#Atc_prevpost <- res

post_nterm <- cbind(rownames(post_nterm), data.frame(post_nterm, row.names=NULL))
names(post_nterm)[1] <- "guide"
head(post_nterm)

post_nterm$Yorf1 <- guide.good.targets[match(post_nterm$guide, guide.good.targets$Guide), "Yorf1"]
post_nterm$gene <- sgd[match(post_nterm$Yorf1, sgd$name), "gene"]
post_nterm$desc <- sgd[match(post_nterm$Yorf1, sgd$name), "desc"]
View(post_nterm)

write.csv(post_nterm, "~/NIKS022/DESeq2_analysis/post_Mrn1_1-200_DEseq.csv")

Atc_prevpost <- cbind(rownames(Atc_prevpost), data.frame(Atc_prevpost, row.names=NULL))
names(Atc_prevpost)[1] <- "guide"
head(Atc_prevpost)

Atc_prevpost$Yorf1 <- guide.good.targets[match(Atc_prevpost$guide, guide.good.targets$Guide), "Yorf1"]
Atc_prevpost$gene <- sgd[match(Atc_prevpost$Yorf1, sgd$name), "gene"]
Atc_prevpost$desc <- sgd[match(Atc_prevpost$Yorf1, sgd$name), "desc"]
View(Atc_prevpost)

write.csv(Atc_prevpost, "~/NIKS022/DESeq2_analysis/prevpost_DEseq.csv")



library(ggplot2)
library(dplyr)
library(ggrepel)

##Volcano plot halo post v pre data
plot_prevpost <- mutate(Atc_prevpost, sig=ifelse(Atc_prevpost$padj<0.05, "padj < 0.05", "padj > 0.05")) #Will have different colors depending on significance
head(plot_prevpost)

Atcplot <- ggplot(plot_prevpost, aes(log2FoldChange, -log10(pvalue)))  +
  geom_point(size=0.5, aes(col=sig)) + #ylim(0,15) +
  scale_color_manual(values=c("goldenrod1", "black")) + 
  ggtitle("Pre vs. post Atc induction")
Atcplot  
Atcplot+geom_text_repel(data=filter(plot_prevpost, padj<0.005 & log2FoldChange < -4.5 |
                                      padj<0.005 & log2FoldChange > 2.5 ), aes(label=gene), size=2.5, segment.alpha=0.15)

##
post_nterm <- mutate(post_nterm, sig=ifelse(post_nterm$padj<0.05, "padj < 0.05", "padj > 0.05")) #Will have different colors depending on significance
head(post_nterm)

plot_post_nterm <- ggplot(post_nterm, aes(log2FoldChange, -log10(pvalue)))  +
  geom_point(size=0.5, aes(col=sig)) +
  scale_color_manual(values=c("aquamarine1", "black")) + 
  ggtitle("Mrn1(1-200) vs. Halo: post-induction")
plot_post_nterm 
plot_post_nterm+geom_text_repel(data=filter(post_nterm, padj<0.005 & log2FoldChange < -2 |
                                      padj<0.005 & log2FoldChange > 2 ), aes(label=gene), size=2.5, segment.alpha=0.15)
