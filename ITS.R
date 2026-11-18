library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
path <- "C:/Users/xuzy0/Desktop/xzy/split"
list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fq", full.names = TRUE))
FWD <- "GGAAGTAAAAGTCGTAACAAGG"
REV <- "GCTGCGTTCTTCATCGATGC"
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
cutadapt <- "D:/software/anaconda/envs/cutadaptenv/Scripts/cutadapt.exe"
system2(cutadapt, args = "--version")
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC)
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i],
                             fnFs.filtN[i], fnRs.filtN[i]))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
cutFs <- sort(list.files(path.cut, pattern = "_R1.fq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fq", full.names = TRUE))
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(3.5, 3.5),
                     truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = F)
head(out)
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
plotErrors(errF, nominalQ = TRUE)
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
table(nchar(getSequences(seqtab.nochim)))
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,
                                                                       getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
                     "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track,"C:/Users/xuzy0/Desktop/xzy/split/quality_control.txt",row.names = T,col.names = T,sep = "\t",quote = F)
library(vegan)
asv_Flattening = as.matrix(rrarefy(seqtab.nochim, min(rowSums(seqtab.nochim))))
unite.ref <- "C:/Users/xuzy0/Desktop/xzy/split/unite_its/sh_general_release_dynamic_04.04.2024_dev.fasta"
taxa <- assignTaxonomy(asv_Flattening, unite.ref, multithread = TRUE, tryRC = TRUE)
feature_table <- data.frame(cbind(taxa,t(asv_Flattening)))
rownames(feature_table) <- paste("ASV",seq(1:nrow(feature_table)),sep = "")
write.table(feature_table,"C:/Users/xuzy0/Desktop/xzy/split/feature_table.txt",row.names = T,col.names = T,sep = "\t",quote = F)
feature_table[,8:107] <- apply(feature_table[,8:107],2,as.numeric)
feature_table_relative <- sweep(feature_table[,8:107],2,colSums(feature_table[,8:107]),"/")
feature_table_relative <- cbind(feature_table[,1:7],feature_table_relative)
write.table(feature_table_relative,"C:/Users/xuzy0/Desktop/xzy/split/feature_table_relative.txt",row.names = T,col.names = T,sep = "\t",quote = F)
feature_table_unflatted <- data.frame(cbind(taxa,t(seqtab.nochim)))
rownames(feature_table_unflatted) <- paste("ASV",seq(1:nrow(feature_table_unflatted)),sep = "")
write.table(feature_table_unflatted,"C:/Users/xuzy0/Desktop/xzy/split/feature_table_unflatted.txt",row.names = T,col.names = T,sep = "\t",quote = F)
feature_table <- data.frame(lapply(feature_table, function(x) ifelse(is.na(x), "Unclassified", x)))
feature_table_unflatted <- data.frame(lapply(feature_table_unflatted, function(x) ifelse(is.na(x), "Unclassified", x)))
feature_table[,8:107] <- apply(feature_table[,8:107], 2, as.numeric)
feature_table_unflatted[,8:107] <- apply(feature_table_unflatted[,8:107], 2, as.numeric)
phylum_table <- aggregate(feature_table[,8:107],by = list(feature_table$Phylum),FUN = sum)
family_table <- aggregate(feature_table[,8:107],by = list(feature_table$Family),FUN = sum)
genus_table <- aggregate(feature_table[,8:107],by = list(feature_table$Genus),FUN = sum)
species_table <- aggregate(feature_table[,8:107],by = list(feature_table$Species),FUN = sum)
colnames(phylum_table)[1] <- "phylum"
colnames(family_table)[1] <- "family"
colnames(genus_table)[1] <- "genus"
colnames(species_table)[1] <- "species"
write.table(phylum_table,"C:/Users/xuzy0/Desktop/xzy/split/chouping/phylum_table.txt",row.names = F,col.names = T,sep = "\t",quote = F)
write.table(family_table,"C:/Users/xuzy0/Desktop/xzy/split/chouping/family_table.txt",row.names = F,col.names = T,sep = "\t",quote = F)
write.table(genus_table,"C:/Users/xuzy0/Desktop/xzy/split/chouping/genus_table.txt",row.names = F,col.names = T,sep = "\t",quote = F)
write.table(species_table,"C:/Users/xuzy0/Desktop/xzy/split/chouping/species_table.txt",row.names = F,col.names = T,sep = "\t",quote = F)
feature_table_unflatted <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/feature_table_unflatted.txt", row.names = 1)
phylum_table_unflatted <- aggregate(feature_table_unflatted[,8:107],by = list(feature_table_unflatted$Phylum),FUN = sum)
family_table_unflatted <- aggregate(feature_table_unflatted[,8:107],by = list(feature_table_unflatted$Family),FUN = sum)
genus_table_unflatted <- aggregate(feature_table_unflatted[,8:107],by = list(feature_table_unflatted$Genus),FUN = sum)
species_table_unflatted <- aggregate(feature_table_unflatted[,8:107],by = list(feature_table_unflatted$Species),FUN = sum)
colnames(phylum_table_unflatted)[1] <- "phylum"
colnames(family_table_unflatted)[1] <- "family"
colnames(genus_table_unflatted)[1] <- "genus"
colnames(species_table_unflatted)[1] <- "species"
write.table(phylum_table_unflatted,"C:/Users/xuzy0/Desktop/xzy/split/weichouping/phylum_table_unflatted.txt",row.names = F,col.names = T,sep = "\t",quote = F)
write.table(family_table_unflatted,"C:/Users/xuzy0/Desktop/xzy/split/weichouping/family_table_unflatted.txt",row.names = F,col.names = T,sep = "\t",quote = F)
write.table(genus_table_unflatted,"C:/Users/xuzy0/Desktop/xzy/split/weichouping/genus_table_unflatted.txt",row.names = F,col.names = T,sep = "\t",quote = F)
write.table(species_table_unflatted,"C:/Users/xuzy0/Desktop/xzy/split/weichouping/species_table_unflatted.txt",row.names = F,col.names = T,sep = "\t",quote = F)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(FSA)
library(pairwiseAdonis)
library(cowplot)
library(ggvenn)
file_path <-"C:/Users/xuzy0/Desktop/xzy/split/feature_table.txt"
setwd("C:/Users/xuzy0/Desktop/ITS\u5206\u6790\u56fe")
feature_table <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/feature_table.txt", row.names = 1)
asv <- feature_table[,8:107]
rarecurve(t(asv), step = 5000,col = "blue", cex = 0.6)
sequence_depth <- c()
asv_number <- c()
sample <- c()
for (i in seq(0,sum(feature_table$A859.Jul),5000)) {
  test <- rrarefy(t(asv),i)
  sequence_depth <- c(sequence_depth,rep(i,100))
  sample <- c(sample,rownames(test))
  for (i in colnames(asv)) {
    asv_number <- c(asv_number,sum(t(test[i,]) != 0))
  }
}
rarefaction_curve <- data.frame(sequence_depth,asv_number,sample)
rarefaction_curve$period <- ""
rarefaction_curve$period[grep("D",rarefaction_curve$sample)] <- "Spring"
rarefaction_curve$period[grep("A",rarefaction_curve$sample)] <- "Summer"
rarefaction_curve$period[grep("B",rarefaction_curve$sample)] <- "Autumn"
rarefaction_curve$period[grep("C",rarefaction_curve$sample)] <- "Winter"
p1 <- ggplot(rarefaction_curve,aes(sequence_depth,asv_number,fill=sample,color=period)) + geom_line()+
  scale_color_manual(values = c("#80cba4","#fbda83","#e69a9b","#4965b0"))+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA,size = .5),
        legend.background = element_blank())
p1
local({
library(FSA)
library(vegan)
input_file <- "D:/split/feature_table.txt"
raw <- read.delim(input_file, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, quote="", comment.char="")
sample_cols <- grep("^[ABCD][0-9]+[.]", names(raw), value=TRUE)
if (length(sample_cols) != 100L) stop("Expected 100 sample columns. Check the header.")
if (anyDuplicated(sample_cols)) stop("Duplicate sample names.")
counts <- as.matrix(raw[,sample_cols,drop=FALSE])
suppressWarnings(storage.mode(counts) <- "numeric")
if (anyNA(counts) || any(!is.finite(counts)) || any(counts<0) || any(abs(counts-round(counts))>1e-7)) {
stop("ASV abundances must be non-negative integer counts, not relative abundances.") }
if (any(colSums(counts)==0)) stop("Samples with zero total counts were found.")
counts <- counts[rowSums(counts)>0,,drop=FALSE]
X <- t(counts)
richness <- vegan::estimateR(X)
alpha_data <- data.frame(
Shannon=vegan::diversity(X,index="shannon"), Simpson=vegan::diversity(X,index="simpson"), Chao1=richness["S.chao1",], ACE=richness["S.ACE",], row.names=rownames(X) )
indices <- c("Shannon", "Simpson", "Chao1", "ACE")
for (index in indices) {
original <- alpha_data[[index]]
alpha_data[[index]] <- suppressWarnings(
as.numeric(as.character(original)) )
if (any(!is.na(original) & is.na(alpha_data[[index]]))) {
stop(paste(index, "contains values that cannot be converted to numeric.")) } }
sample_ids <- rownames(alpha_data)
if (anyDuplicated(sample_ids)) stop("Duplicate sample names.")
sample_month <- sub(".*\\.", "", sample_ids)
season_by_month <- c(
Apr="Spring", May="Spring", Jun="Spring", Jul="Summer", Aug="Summer", Sep="Summer", Oct="Autumn", Nov="Autumn", Dec="Winter", Jan="Winter", Feb="Winter", Mar="Winter" )
alpha_data$Seasons <- factor(unname(season_by_month[sample_month]), levels=c("Spring","Summer","Autumn","Winter"))
if (anyNA(alpha_data$Seasons)) stop("Unrecognized month suffix in sample names.")
if (any(table(alpha_data$Seasons)!=25L)) {
print(table(alpha_data$Seasons))
stop("Expected 25 samples per season. Check the data.") }
prefix_season <- c(A="Summer",B="Autumn",C="Winter",D="Spring")
expected_season <- unname(prefix_season[substr(sample_ids,1,1)])
if (anyNA(expected_season) || any(as.character(alpha_data$Seasons)!=expected_season)) {
stop("Sample prefixes and month-based seasons do not match.") }
alpha_data$Months <- factor(
sub(".*\\.", "", rownames(alpha_data)), levels = c(
"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" ) )
if (anyNA(alpha_data$Seasons)) {
stop("Missing season assignments. Check the season metadata.") }
if (anyNA(alpha_data$Months)) {
stop("Cannot identify months from sample names. Expected names such as A859.Jul.") }
print(table(alpha_data$Seasons))
print(table(alpha_data$Months))
output_dir <- "Alpha_Dunn_BH_results"
dir.create(output_dir, showWarnings = FALSE)
write.csv(data.frame(Sample=rownames(alpha_data),alpha_data), file.path(output_dir,'Sample_groups_and_alpha_indices.csv'),row.names=FALSE)
kw_list <- list()
dunn_list <- list()
for (group_name in c("Seasons", "Months")) {
for (index in indices) {
test_data <- data.frame(
value = alpha_data[[index]], group = alpha_data[[group_name]] )
valid <- is.finite(test_data$value) &
!is.na(test_data$group)
test_data <- test_data[valid, , drop = FALSE]
test_data$group <- droplevels(test_data$group)
if (nlevels(test_data$group) < 2) {
stop(paste(group_name, index, "has fewer than two valid groups.")) }
kw <- kruskal.test(
value ~ group, data = test_data )
key <- paste(group_name, index, sep = "_")
kw_list[[key]] <- data.frame(
grouping = group_name, index = index, N = nrow(test_data), excluded_N = sum(!valid), H = unname(kw$statistic), Df = unname(kw$parameter), p = kw$p.value )
dt <- FSA::dunnTest(
value ~ group, data = test_data, method = "bh", two.sided = TRUE
)$res
result <- data.frame(
grouping = group_name, index = index, comparison = dt$Comparison, Z = dt$Z, p = dt$P.unadj, q_BH = dt$P.adj )
result$significance <- ifelse(
result$q_BH < 0.001, "***", ifelse(
result$q_BH < 0.01, "**", ifelse(result$q_BH < 0.05, "*", "ns") ) )
result$KW_p <- kw$p.value
result$KW_p_less_0.05 <- kw$p.value < 0.05
dunn_list[[key]] <- result
write.csv(
result, file.path(
output_dir, paste0(key, "_Dunn_BH.csv")
), row.names = FALSE ) } }
kw_result <- do.call(rbind, kw_list)
rownames(kw_result) <- NULL
kw_result$q_BH <- NA_real_
for (group_name in c("Seasons", "Months")) {
select <- kw_result$grouping == group_name
kw_result$q_BH[select] <- p.adjust(
kw_result$p[select], method = "BH" ) }
write.csv(
kw_result, file.path(output_dir, "Overall_Kruskal_Wallis_BH.csv"), row.names = FALSE )
dunn_result <- do.call(rbind, dunn_list)
rownames(dunn_result) <- NULL
write.csv(
dunn_result, file.path(output_dir, "All_Dunn_BH_results.csv"), row.names = FALSE )
write.csv(
dunn_result[dunn_result$q_BH < 0.05, ], file.path(output_dir, "Dunn_BH_q_below_0.05.csv"), row.names = FALSE )
print(kw_result)
print(dunn_result)
cat(
"\nResults saved to: ", normalizePath(output_dir, winslash = "/"), "\n" )
library(ggplot2)
library(multcompView)
season_colors <- c(
Spring = "#d1e6c1", Summer = "#faf0b6", Autumn = "#e69a9b", Winter = "#c0dbe5" )
month_colors <- c(
Jan = "#5ab1ef", Feb = "#ca93b5", Mar = "#a3cc63", Apr = "#f39a41", May = "#fde8d3", Jun = "#8bd2dc", Jul = "#faf0b6", Aug = "#ebbdbd", Sep = "#d1e6c1", Oct = "#f0f8ff", Nov = "#eabaa1", Dec = "#33adad" )
make_alpha_plot <- function(group_name,index) {
z <- data.frame(Group=alpha_data[[group_name]],Value=alpha_data[[index]])
z <- z[is.finite(z$Value)&!is.na(z$Group),,drop=FALSE]
z$Group <- droplevels(z$Group)
kw <- kw_result[kw_result$grouping==group_name & kw_result$index==index,,drop=FALSE]
dt <- dunn_result[dunn_result$grouping==group_name & dunn_result$index==index,,drop=FALSE]
groups <- levels(z$Group)
med <- tapply(z$Value,z$Group,median)
ordered <- names(sort(med,decreasing=TRUE))
qmat <- matrix(NA_real_,length(groups),length(groups),dimnames=list(groups,groups))
diag(qmat) <- 1
for(i in seq_len(nrow(dt))) {
pair <- strsplit(as.character(dt$comparison[i]),' - ',fixed=TRUE)[[1]]
if(length(pair)!=2 || !all(pair%in%groups)) stop('Cannot match pairwise comparison group names.')
qmat[pair[1],pair[2]] <- qmat[pair[2],pair[1]] <- dt$q_BH[i] }
if(anyNA(qmat)) stop('Pairwise comparisons are incomplete or contain missing q-values; letters cannot be generated.')
letters_df <- data.frame(Group=factor(groups,levels=groups),Letters=NA_character_)
show_letters <- is.finite(kw$q_BH[1]) && kw$q_BH[1]<.05
if(show_letters) {
cld <- multcompView::multcompLetters(qmat[ordered,ordered],threshold=.05)
lm <- cld$LetterMatrix
for(i in seq_len(length(groups)-1)) for(j in (i+1):length(groups)) {
shared <- any(lm[groups[i],] & lm[groups[j],])
if(shared != (qmat[groups[i],groups[j]]>=.05)) stop('Letter assignments do not match Dunn-BH results.') }
letters_df$Letters <- unname(cld$Letters[groups]) }
span <- diff(range(z$Value)); if(span==0) span <- max(abs(z$Value),1)
letters_df$y <- max(z$Value) + .07*span
letters_df$grouping <- group_name; letters_df$index <- index
letters_df$overall_KW_q <- kw$q_BH[1]
write.csv(letters_df,file.path(output_dir,paste0(group_name,'_',index,'_Letter_groups.csv')),row.names=FALSE)
colors <- if(group_name=='Seasons') season_colors else month_colors
note <- if(show_letters) 'Dunn-BH q<0.05: groups without a shared letter differ significantly.' else
'Overall KW BH q >= 0.05; post hoc letters omitted.'
p <- ggplot(z,aes(Group,Value,fill=Group))+ geom_boxplot(width=.65,outlier.shape=16,outlier.size=1.5,outlier.colour='black',outlier.alpha=.6)+ scale_fill_manual(values=colors)+scale_y_continuous(expand=expansion(mult=c(.05,.15)))+ labs(x=NULL,y=index,title=paste(index,group_name,sep=' - '), subtitle=bquote('Kruskal-Wallis: H = '*.(sprintf('%.3f',kw$H[1]))*
', '*italic(p)*' = '*.(sprintf('%.4g',kw$p[1]))*
', '*italic(p)[BH]*' = '*.(sprintf('%.4g',kw$q_BH[1]))),caption=note)+ theme_classic(base_size=13)+theme(legend.position='none',panel.border=element_rect(fill=NA,linewidth=.5), axis.text.x=element_text(angle=if(group_name=='Months') 45 else 0,hjust=if(group_name=='Months') 1 else .5), plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0))
if(show_letters) p <- p+geom_text(data=letters_df,aes(x=Group,y=y,label=Letters), inherit.aes=FALSE,size=5,vjust=0)
filename <- paste0(group_name,'_',index,'_Dunn_BH_Letters')
w <- if(group_name=='Months') 11 else 7
ggsave(file.path(output_dir,paste0(filename,'.pdf')),p,width=w,height=6)
p }
alpha_plots <- list()
for(group_name in c('Seasons','Months')) for(index in indices) {
key <- paste(group_name,index,sep='_')
alpha_plots[[key]] <- make_alpha_plot(group_name,index)
print(alpha_plots[[key]]) }
})
local({
library(vegan)
library(ggplot2)
set.seed(123)
nperm <- 999
feature_table <- read.delim( "D:/split/feature_table.txt", header = TRUE, row.names = 1, check.names = FALSE )
asv <- feature_table[, -(1:7), drop = FALSE]
data <- as.data.frame(t(asv))
month_names <- sub(".*\\.", "", rownames(data))
group <- data.frame( row.names = rownames(data), Month = factor( month_names, levels = c( "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" ) ) )
if (anyNA(group$Month)) {
stop("Unrecognized month in sample names. Check month abbreviations after the final dot.") }
print(table(group$Month))
adonis <- adonis2( data ~ Month, data = group, permutations = nperm, method = "bray" )
permanova_result <- data.frame( grouping = "month", N = nrow(data), Df = adonis$Df[1], R2 = adonis$R2[1], pseudo_F = adonis$F[1], permutations = nperm, p = adonis$`Pr(>F)`[1] )
permanova_result$q_BH <- p.adjust( permanova_result$p, method = "BH" )
print(permanova_result)
write.csv( permanova_result, "Month_PERMANOVA_overall.csv", row.names = FALSE )
adonis_lab <- sprintf( paste0( "PERMANOVA: R² = %.3f, pseudo-F = %.3f", "\nPermutations = %d, p = %.4g, BH q = %.4g"
), permanova_result$R2, permanova_result$pseudo_F, permanova_result$permutations, permanova_result$p, permanova_result$q_BH )
distance <- vegdist(data, method = "bray")
month_pairs <- combn( levels(droplevels(group$Month)), 2, simplify = FALSE )
pairwise_result <- do.call( rbind, lapply(month_pairs, function(z) {
index <- group$Month %in% z
sub_distance <- as.dist( as.matrix(distance)[index, index] )
sub_group <- droplevels( group[index, , drop = FALSE] )
fit <- adonis2( sub_distance ~ Month, data = sub_group, permutations = nperm )
data.frame( Group1 = z[1], Group2 = z[2], N = sum(index), Df = fit$Df[1], R2 = fit$R2[1], pseudo_F = fit$F[1], permutations = nperm, p = fit$`Pr(>F)`[1] ) }) )
pairwise_result$q_BH <- p.adjust( pairwise_result$p, method = "BH" )
print(pairwise_result)
write.csv( pairwise_result, "Month_PERMANOVA_pairwise_BH.csv", row.names = FALSE )
dispersion <- betadisper( distance, group = group$Month, bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE )
x <- permutest( dispersion, permutations = nperm )
permdisp_result <- data.frame( grouping = "month", F = x$tab$F[1], permutations = nperm, p = x$tab$`Pr(>F)`[1] )
permdisp_result$q_BH <- p.adjust( permdisp_result$p, method = "BH" )
print(permdisp_result)
write.csv( permdisp_result, "Month_PERMDISP_overall.csv", row.names = FALSE )
betadisper_lab <- sprintf( "PERMDISP: F = %.3f, p = %.4g, BH q = %.4g", permdisp_result$F, permdisp_result$p, permdisp_result$q_BH )
all_lab <- paste( adonis_lab, betadisper_lab, sep = "\n" )
pcoa <- cmdscale( distance, k = nrow(data) - 1, eig = TRUE )
plot_data <- as.data.frame( pcoa$points[, 1:2, drop = FALSE] )
names(plot_data) <- c("PCoA1", "PCoA2")
plot_data$Months <- group[
match(rownames(plot_data), rownames(group)), "Month"
]
eig <- pcoa$eig
month_colors <- c( Jan = "#5ab1ef", Feb = "#ca93b5", Mar = "#a3cc63", Apr = "#f39a41", May = "#b39b42", Jun = "#8bd2dc", Jul = "#faf0b6", Aug = "#ebbdbd", Sep = "#d1e6c1", Oct = "#9986bd", Nov = "#eabaa1", Dec = "#33adad" )
p_month <- ggplot( plot_data, aes(PCoA1, PCoA2, color = Months)
) + geom_point(alpha = 1, size = 3) + scale_color_manual(values = month_colors) + labs( x = paste0( "PCoA 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)"
), y = paste0( "PCoA 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)"
), title = all_lab
) + stat_ellipse( level = 0.95, show.legend = FALSE, linewidth = 0.5
) + theme( panel.background = element_blank(), panel.border = element_rect( fill = NA, linewidth = 0.5
), plot.title = element_text( hjust = 0, size = 10, margin = margin(b = 10)
), legend.position = "bottom", legend.title = element_text(size = 14), legend.text = element_text(size = 12), legend.key = element_blank(), axis.text = element_text(size = 14), axis.title = element_text(size = 14), aspect.ratio = 1.1
) + guides( color = guide_legend(nrow = 3, byrow = TRUE) )
print(p_month)
ggsave( filename = "D:/split/Month_PCoA.pdf", plot = p_month, device = "pdf", width = 8, height = 7, units = "in" )
library(vegan)
library(ggplot2)
set.seed(123)
nperm <- 999
data <- as.data.frame(t(asv))
season_map <- c( A = "Summer", B = "Autumn", C = "Winter", D = "Spring" )
group <- data.frame( row.names = rownames(data), V2 = factor( unname(season_map[substr(rownames(data), 1, 1)]), levels = c("Spring", "Summer", "Autumn", "Winter") ) )
if (anyNA(group$V2)) {
stop("Unrecognized season. Sample names must start with A, B, C, or D.") }
print(table(group$V2))
adonis <- adonis2( data ~ V2, data = group, permutations = nperm, method = "bray" )
permanova_result <- data.frame( grouping = "season", N = nrow(data), Df = adonis$Df[1], R2 = adonis$R2[1], pseudo_F = adonis$F[1], permutations = nperm, p = adonis$`Pr(>F)`[1] )
permanova_result$q_BH <- p.adjust( permanova_result$p, method = "BH" )
print(permanova_result)
write.csv( permanova_result, "Season_PERMANOVA_overall.csv", row.names = FALSE )
adonis_lab <- sprintf( paste0( "PERMANOVA: R² = %.3f, pseudo-F = %.3f", "\nPermutations = %d, p = %.4g, BH q = %.4g"
), permanova_result$R2, permanova_result$pseudo_F, permanova_result$permutations, permanova_result$p, permanova_result$q_BH )
pairwise_distance <- vegdist(data, method = "bray")
season_pairs <- combn( levels(droplevels(group$V2)), 2, simplify = FALSE )
pairwise_result <- do.call( rbind, lapply(season_pairs, function(z) {
index <- group$V2 %in% z
sub_distance <- as.dist( as.matrix(pairwise_distance)[index, index] )
sub_group <- droplevels( group[index, , drop = FALSE] )
fit <- adonis2( sub_distance ~ V2, data = sub_group, permutations = nperm )
data.frame( Group1 = z[1], Group2 = z[2], N = sum(index), Df = fit$Df[1], R2 = fit$R2[1], pseudo_F = fit$F[1], permutations = nperm, p = fit$`Pr(>F)`[1] ) }) )
pairwise_result$q_BH <- p.adjust( pairwise_result$p, method = "BH" )
print(pairwise_result)
write.csv( pairwise_result, "Season_PERMANOVA_pairwise_BH.csv", row.names = FALSE )
distance <- vegdist( data, method = "bray" )
dispersion <- betadisper( distance, group = group$V2, bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE )
x <- permutest( dispersion, permutations = nperm )
permdisp_result <- data.frame( grouping = "season", F = x$tab$F[1], permutations = nperm, p = x$tab$`Pr(>F)`[1] )
permdisp_result$q_BH <- p.adjust( permdisp_result$p, method = "BH" )
print(permdisp_result)
write.csv( permdisp_result, "Season_PERMDISP_overall.csv", row.names = FALSE )
betadisper_lab <- sprintf( "PERMDISP: F = %.3f, p = %.4g, BH q = %.4g", permdisp_result$F, permdisp_result$p, permdisp_result$q_BH )
all_lab <- paste( adonis_lab, betadisper_lab, sep = "\n" )
pcoa <- cmdscale( distance, k = nrow(data) - 1, eig = TRUE )
plot_data <- as.data.frame(pcoa$points[, 1:2, drop = FALSE])
names(plot_data) <- c("PCoA1", "PCoA2")
plot_data$Seasons <- group[
match(rownames(plot_data), rownames(group)), "V2"
]
eig <- pcoa$eig
p6 <- ggplot( data = plot_data, aes( x = PCoA1, y = PCoA2, color = Seasons )
) + geom_point( alpha = 1, size = 3
) + scale_color_manual( values = c( "Spring" = "#80cba4", "Summer" = "#fbda83", "Autumn" = "#e69a9b", "Winter" = "#4965b0" )
) + labs( x = paste0( "PCoA 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)"
), y = paste0( "PCoA 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)"
), title = all_lab
) + stat_ellipse( level = 0.95, show.legend = FALSE, linewidth = 0.5, alpha = 1
) + theme( panel.background = element_blank(), panel.border = element_rect( fill = NA, linewidth = 0.5
), panel.spacing.x = grid::unit(0.3, "cm"), plot.title = element_text( hjust = 0, size = 10, margin = margin(b = 10)
), legend.position = "bottom", legend.title = element_text(size = 14), legend.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title = element_text(size = 14), legend.key = element_blank(), strip.background = element_rect( color = "black", linewidth = 1
), strip.text = element_text(size = 14), strip.placement = "inside", aspect.ratio = 1.1
) + guides(fill = "none")
print(p6)
ggsave( filename = "Season_PCoA.pdf", plot = p6, width = 8, height = 7 )
ggsave( filename = "Season_PCoA.png", plot = p6, width = 8, height = 7, dpi = 300 )
ggsave( filename = "D:/split/Season_PCoA.pdf", plot = p6, device = "pdf", width = 8, height = 7, units = "in" )
})
local({
library(vegan)
set.seed(20261110)
nperm <- 999L
output_dir <- 'Aitchison_results'
dir.create(output_dir,showWarnings=FALSE)
save_table <- function(x,name) write.csv(x,file.path(output_dir,paste0(name,'.csv')),row.names=FALSE,fileEncoding='UTF-8')
input_file <- 'D:/split/feature_table.txt'
if (!file.exists(input_file)) stop('Input file not found: ',input_file)
feature_table <- read.delim(input_file, header=TRUE, row.names=1, sep='\t', check.names=FALSE, quote='', comment.char='')
colnames(feature_table) <- trimws(sub('^\ufeff','',colnames(feature_table)))
if(anyDuplicated(colnames(feature_table))) stop('Duplicate column names after trimming whitespace.')
sample_cols <- grepl('^[ABCD][0-9]+\\.(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$', colnames(feature_table))
if (!any(sample_cols)) stop('No sample columns matching A859.Jul were found. Check delimiters and column names.')
annotation_names <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
unknown <- setdiff(colnames(feature_table)[!sample_cols],annotation_names)
if (length(unknown)) stop(paste('Unrecognized columns: ',paste(sprintf('[%s]',unknown),collapse=', '), '; check column names.'))
asv <- feature_table[,sample_cols,drop=FALSE]
if (anyDuplicated(rownames(asv))) stop('Duplicate ASV identifiers.')
message('Loaded ',nrow(asv),' ASVs and ',ncol(asv),' samples.')
if(is.data.frame(asv) && !all(vapply(asv,is.numeric,logical(1)))) stop('The ASV table contains non-numeric columns.')
X <- t(as.matrix(asv))
stopifnot(is.numeric(X),all(is.finite(X)),all(X>=0),all(rowSums(X)>0), !is.null(rownames(X)),!anyDuplicated(rownames(X)))
X <- X[,colSums(X)>0,drop=FALSE]
if(file.exists('metadata.tsv')) {
meta <- read.delim('metadata.tsv',stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(all(c('sample_id','season','month')%in%names(meta)),!anyDuplicated(meta$sample_id), setequal(meta$sample_id,rownames(X)))
meta <- meta[match(rownames(X),meta$sample_id),]
} else {
ids <- rownames(X)
meta <- data.frame(sample_id=ids, season=unname(c(A='Summer',B='Autumn',C='Winter',D='Spring')[substr(ids,1,1)]), month=sub('.*\\.','',ids)) }
month_season <- c(Jan='Winter',Feb='Winter',Mar='Winter',Apr='Spring',May='Spring',Jun='Spring', Jul='Summer',Aug='Summer',Sep='Summer',Oct='Autumn',Nov='Autumn',Dec='Winter')
if(anyNA(meta[,c('season','month')]) || !all(meta$month%in%names(month_season)))
stop('Cannot identify seasons or months. Supply metadata.tsv rather than assigning groups by column position.')
if(!all(meta$season==unname(month_season[meta$month]))) stop('Months and study season definitions do not match. Check metadata.')
meta$season <- factor(meta$season,levels=c('Spring','Summer','Autumn','Winter'))
meta$month <- factor(meta$month,levels=c('Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'))
save_table(meta,'Sample_group_check')
print(table(meta$season,meta$month))
distance <- vegdist(X,method='robust.aitchison')
global <- list(); disp <- list(); pairwise <- list()
for(term in c('season','month')) {
g <- droplevels(meta[[term]])
if(nlevels(g)<2 || any(table(g)<2)) stop('At least two samples per group are required.')
fit <- adonis2(distance~g,permutations=nperm)
global[[term]] <- data.frame(grouping=term,N=nrow(X),Df=fit$Df[1],pseudo_F=fit$F[1], R2=fit$R2[1],p=fit$`Pr(>F)`[1],requested_permutations=nperm)
capture.output(fit,file=file.path(output_dir,paste0(term,'_PERMANOVA_full_output.txt')))
bd <- betadisper(distance,g,type='median',bias.adjust=FALSE,add=FALSE)
bt <- permutest(bd,permutations=nperm)
disp[[term]] <- data.frame(grouping=term,F=bt$tab$F[1],p=bt$tab$`Pr(>F)`[1],requested_permutations=nperm)
capture.output(bt,file=file.path(output_dir,paste0(term,'_PERMDISP_full_output.txt')))
save_table(data.frame(sample_id=meta$sample_id,group=g,distance_to_median=bd$distances),paste0(term,'_Dispersion'))
pdf(file.path(output_dir,paste0(term,'_PERMDISP.pdf')),width=9,height=5)
boxplot(bd,ylab='Distance to group spatial median',xlab=term)
stripchart(bd$distances~g,vertical=TRUE,method='jitter',add=TRUE,pch=16,col='#4477AA80')
dev.off()
ptab <- lapply(combn(levels(g),2,simplify=FALSE),function(z) {
keep <- g%in%z; dsub <- as.dist(as.matrix(distance)[keep,keep]); subgroup <- droplevels(g[keep])
f <- adonis2(dsub~subgroup,permutations=nperm)
data.frame(group1=z[1],group2=z[2],N=sum(keep),pseudo_F=f$F[1],R2=f$R2[1], p=f$`Pr(>F)`[1],requested_permutations=nperm) })
ptab <- do.call(rbind,ptab); ptab$q_BH <- p.adjust(ptab$p,method='BH')
pairwise[[term]] <- ptab
save_table(ptab,paste0(term,'_Pairwise_PERMANOVA_BH')) }
global <- do.call(rbind,global)
global$q_BH <- p.adjust(global$p,'BH')
disp <- do.call(rbind,disp)
disp$q_BH <- p.adjust(disp$p,'BH')
save_table(global,'Overall_PERMANOVA')
save_table(disp,'Overall_PERMDISP')
print(global); print(disp)
sentences <- vapply(seq_len(nrow(global)),function(i) {
a <- global[i,]
sprintf('%s: Robust Aitchison PERMANOVA, pseudo-F = %.4f, R2 = %.4f, requested permutations = %d, p = %.5g, BH q = %.5g.', a$grouping,a$pseudo_F,a$R2,a$requested_permutations,a$p,a$q_BH)
},character(1))
writeLines(sentences,file.path(output_dir,'PERMANOVA_report.txt'),useBytes=TRUE)
capture.output(sessionInfo(),file=file.path(output_dir,'sessionInfo.txt'))
library(ggplot2)
pc <- cmdscale(distance,k=2,eig=TRUE)
pct <- 100*pc$eig[1:2]/sum(pc$eig)
month_colors <- c(Jan='#5ab1ef',Feb='#ca93b5',Mar='#a3cc63',Apr='#f39a41',May='#fde8d3',Jun='#8bd2dc',Jul='#faf0b6',Aug='#ebbdbd',Sep='#d1e6c1',Oct="#d87a80",Nov='#eabaa1',Dec='#33adad')
season_colors <- c(Spring='#80cba4',Summer='#fbda83',Autumn='#e69a9b',Winter='#4965b0')
for(term in c('season','month')) {
g <- meta[[term]]
if(term=='month') g <- factor(g,levels=names(month_colors))
pd <- data.frame(PCoA1=pc$points[,1],PCoA2=pc$points[,2],Group=g)
a <- global[global$grouping==term,]; b <- disp[disp$grouping==term,]
label <- bquote(atop('Robust Aitchison', atop('PERMANOVA: '*italic(R)^2*' = '*.(sprintf('%.3f',a$R2))*
', '*italic('pseudo-F')*' = '*.(sprintf('%.3f',a$pseudo_F))*
', permutations = '*.(nperm), atop(italic(p)*' = '*.(sprintf('%.4g',a$p))*
', BH '*italic(q)*' = '*.(sprintf('%.4g',a$q_BH)), 'PERMDISP: '*italic(F)*' = '*.(sprintf('%.3f',b$F))*
', '*italic(p)*' = '*.(sprintf('%.4g',b$p))*
', BH '*italic(q)*' = '*.(sprintf('%.4g',b$q_BH))))))
colors <- if(term=='month') month_colors else season_colors
p <- ggplot(pd,aes(PCoA1,PCoA2,color=Group))+geom_point(size=3)+ scale_color_manual(values=colors)+stat_ellipse(level=0.95,show.legend=FALSE,linewidth=0.5)+ labs(x=sprintf('PCoA 1 (%.3f%%)',pct[1]),y=sprintf('PCoA 2 (%.3f%%)',pct[2]),title=label,color=if(term=='month') 'Months' else 'Seasons')+ theme_classic(base_size=14)+theme(panel.border=element_rect(fill=NA,linewidth=0.5),plot.title=element_text(size=10,hjust=0),legend.position='bottom',aspect.ratio=1.1)
if(term=='month') p <- p+guides(color=guide_legend(nrow=3,ncol=4,byrow=FALSE))
print(p)
ggsave(file.path(output_dir,paste0(term,'_Robust_Aitchison.pdf')),p,width=8,height=7)
ggsave(file.path(output_dir,paste0(term,'_Robust_Aitchison.png')),p,width=8,height=7,dpi=300) }
message('Results directory: ',normalizePath(output_dir,winslash='/'))
})
local({
library(vegan)
library(ggplot2)
library(permute)
input_file <- 'D:/split/chouping/core_genus.txt'
output_dir <- file.path(getwd(),'106_shared_genera_Beta_results')
requested <- 999L
seed <- 123L
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
save_table <- function(x,n) write.csv(x,file.path(output_dir,paste0(n,'.csv')),row.names=FALSE)
con <- file(input_file,'rb'); bom <- readBin(con,'raw',n=3); close(con)
enc <- if(length(bom)>=2 && identical(bom[1:2],as.raw(c(255,254)))) 'UTF-16LE' else
       if(length(bom)>=2 && identical(bom[1:2],as.raw(c(254,255)))) 'UTF-16BE' else 'UTF-8'
G <- read.delim(input_file,header=TRUE,row.names=1,check.names=FALSE,fileEncoding=enc,
                quote='',comment.char='')
if(!all(vapply(G,is.numeric,logical(1)))) stop('Non-numeric sample columns were found. Check the input format.')
X <- t(as.matrix(G))
stopifnot(all(is.finite(X)),all(X>=0),all(rowSums(X)>0),
          !anyDuplicated(rownames(X)),!anyDuplicated(colnames(X)))
if(ncol(X)!=106||nrow(X)!=100) warning('Input dimensions differ from 106 genera and 100 samples. Check the genus list.')
X <- X[,colSums(X)>0,drop=FALSE]
months <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
meta <- data.frame(sample_id=rownames(X),Month=factor(sub('.*\\.','',rownames(X)),levels=months))
if(anyNA(meta)||any(table(meta$Month)<2)) stop('Unrecognized months or fewer than two samples per month.')
save_table(data.frame(meta,abundance_sum=rowSums(X)),'Input_sample_check')
d <- vegdist(X,method='bray')
saveRDS(d,file.path(output_dir,'Bray_distance.rds'))
make_perm <- function(n,s) {
  set.seed(s)
  shuffleSet(n,nset=requested,control=how(within=Within(type='free')))
}
P <- make_perm(nrow(X),seed)
extract_fit <- function(fit,n,actual) {
  j <- match('Residual',rownames(fit))
  if(is.na(j)) stop('The PERMANOVA output has no Residual row.')
  f <- fit$F[1]; r <- fit$R2[1]; df1 <- fit$Df[1]; df2 <- fit$Df[j]
  rr <- f*df1/(f*df1+df2)
  if(!isTRUE(all.equal(as.numeric(r),as.numeric(rr),tolerance=1e-7))) stop('Pseudo-F and R-squared are inconsistent. Check the output structure.')
  data.frame(N=n,Df=df1,Df_residual=df2,pseudo_F=f,R2=r,
    requested_permutations=requested,actual_permutations=actual,p=fit$`Pr(>F)`[1])
}
fit <- adonis2(d~Month,data=meta,permutations=P,add=FALSE)
overall <- extract_fit(fit,nrow(X),nrow(P))
overall$q_BH <- p.adjust(overall$p,'BH')
save_table(overall,'PERMANOVA_overall_original_Bray')
capture.output(fit,file=file.path(output_dir,'PERMANOVA_full_output.txt'))
bd <- betadisper(d,meta$Month,type='median',bias.adjust=FALSE,add='lingoes')
bt <- permutest(bd,permutations=P)
disp <- data.frame(method='Lingoes; spatial median; bias.adjust=FALSE',
 F=bt$tab$F[1],requested_permutations=requested,actual_permutations=nrow(P),p=bt$tab$`Pr(>F)`[1])
disp$q_BH <- p.adjust(disp$p,'BH')
save_table(disp,'PERMDISP_overall_Lingoes')
save_table(data.frame(meta,distance_to_median=bd$distances),'PERMDISP_sample_distances')
capture.output(bt,file=file.path(output_dir,'PERMDISP_full_output.txt'))
pairs <- combn(months,2,simplify=FALSE)
pair_results <- lapply(seq_along(pairs),function(i) {
  z <- pairs[[i]]; ix <- meta$Month%in%z
  dd <- as.dist(as.matrix(d)[ix,ix,drop=FALSE]); mm <- droplevels(meta[ix,,drop=FALSE])
  pm <- make_perm(sum(ix),seed+i)
  ff <- adonis2(dd~Month,data=mm,permutations=pm,add=FALSE)
  data.frame(group1=z[1],group2=z[2],extract_fit(ff,sum(ix),nrow(pm)))
})
pair_results <- do.call(rbind,pair_results)
pair_results$q_BH <- p.adjust(pair_results$p,'BH')
save_table(pair_results,'Month_Pairwise_PERMANOVA_BH')
colors <- c(Jan='#5ab1ef',Feb='#ca93b5',Mar='#a3cc63',Apr='#f39a41',May='#fde8d3',Jun='#8bd2dc',
 Jul='#faf0b6',Aug='#ebbdbd',Sep='#d1e6c1',Oct='#f0f8ff',Nov='#eabaa1',Dec='#33adad')
for(mode in c('Uncorrected','Lingoes')) {
  pc <- if(mode=='Uncorrected') wcmdscale(d,k=2,eig=TRUE,add=FALSE) else
    wcmdscale(d,k=2,eig=TRUE,add='lingoes')
  if(ncol(pc$points)<2) stop('Fewer than two valid PCoA axes.')
  pct <- 100*pc$eig[1:2]/sum(pc$eig[pc$eig>0])
  pd <- data.frame(meta,PCoA1=pc$points[,1],PCoA2=pc$points[,2])
  save_table(pd,paste0('PCoA_coordinates_',mode))
  save_table(data.frame(axis=seq_along(pc$eig),eigenvalue=pc$eig),paste0('Eigenvalues_',mode))
  txt <- sprintf("atop('PERMANOVA: '*italic('pseudo-F')*' = %.3f, '*italic(R)^2*' = %.3f',atop('Permutations = %d, '*italic(p)*' = %.4g, '*italic(q)*' = %.4g','PERMDISP (Lingoes): '*italic(F)*' = %.3f, '*italic(p)*' = %.4g'))",
    overall$pseudo_F,overall$R2,nrow(P),overall$p,overall$q_BH,disp$F,disp$p)
  p <- ggplot(pd,aes(PCoA1,PCoA2,color=Month))+geom_point(size=2.5)+
    stat_ellipse(level=.95,show.legend=FALSE,linewidth=.6)+scale_color_manual(values=colors,breaks=months)+
    annotate('text',x=-Inf,y=Inf,label=txt,parse=TRUE,hjust=-.01,vjust=1.1,size=3,color='black')+
    labs(x=sprintf('PCoA 1 (%.3f%%)',pct[1]),y=sprintf('PCoA 2 (%.3f%%)',pct[2]),
      color='Months',tag='f',caption=paste('PCoA:',mode,'; PERMANOVA: original Bray-Curtis; PERMDISP: Lingoes.'))+
    scale_y_continuous(expand=expansion(mult=c(.1,.4)))+theme_classic(base_size=12)+
    theme(panel.border=element_rect(fill=NA,color='grey50',linewidth=.5),axis.line=element_blank(),
      legend.position='bottom',aspect.ratio=1.1,plot.caption=element_text(size=8))+
    guides(color=guide_legend(nrow=3,ncol=4,byrow=FALSE))
  print(p)
  ggsave(file.path(output_dir,paste0('Month_PCoA_',mode,'.pdf')),p,width=8,height=7)
  ggsave(file.path(output_dir,paste0('Month_PCoA_',mode,'.png')),p,width=8,height=7,dpi=300)
}
saveRDS(list(input=X,metadata=meta,permutations=P,adonis=fit,betadisper=bd,permdisp=bt),file.path(output_dir,'Full_analysis_objects.rds'))
writeLines(c('Input relative abundances are used unchanged.',
 'PERMANOVA: original Bray-Curtis, no Lingoes correction; pairwise BH across all 66 tests.',
 'PERMDISP: spatial median, Lingoes correction, bias.adjust=FALSE.',
 'PCoA figures: uncorrected and Lingoes alternatives; axis percentages use positive eigenvalues.',
 'Unrestricted permutations assume exchangeable samples; known repeated individuals need a design-aware scheme.',
 'No manual reconstruction of PERMDISP p-values; reported values come directly from vegan.',
 'The single monthly overall test has BH q=p. Do not copy old values into new outputs.'),file.path(output_dir,'Analysis_settings.txt'))
capture.output(sessionInfo(),file=file.path(output_dir,'sessionInfo.txt'))
print(overall); print(disp); message('Results path: ',normalizePath(output_dir,winslash='/'))
})
phylum_table <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/chouping/phylum_table.txt",row.names = 1)
phylum <- sweep(phylum_table,2,colSums(phylum_table),"/")
phylum <- phylum[order(rowSums(phylum),decreasing = T),]
phylum <- phylum[1:3,]
stat <- data.frame(row.names = rownames(phylum),mean=apply(phylum,1,mean),sd=apply(phylum,1,sd))
phylum["Other",] <- (1-colSums(phylum))
x <- rownames(phylum)
phylum <- cbind(rownames(phylum),phylum)
phylum <- melt(phylum)
colnames(phylum) <- c("phylum","seasons","abundance")
phylum$phylum <- factor(phylum$phylum,levels = rev(x))
phylum$seasons <- ifelse(grepl("^A", phylum$seasons), "Summer",
                         ifelse(grepl("^B", phylum$seasons), "Autumn",
                                ifelse(grepl("^C", phylum$seasons), "Winter",
                                       ifelse(grepl("^D", phylum$seasons), "Spring", phylum$seasons))))
phylum$seasons <- factor(phylum$seasons, levels = c( "Spring","Summer", "Autumn", "Winter"))
phylum$abundance <- phylum$abundance / 25
p8 <- ggplot(phylum,aes(seasons,abundance,fill=phylum))+geom_bar(stat = "identity",width = .75)+
  scale_fill_manual(values = c("#e5cf0d","#077535","#fee3ce", "#a8dde1"))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
p8
phylum_table <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/chouping/phylum_table.txt")
phylum_group <- data.frame(row.names=phylum_table$phylum,
                           Jul=rowSums(phylum_table[,2:11],na.rm = TRUE),
                           Aug=rowSums(phylum_table[,12:21],na.rm = TRUE),
                           Sep=rowSums(phylum_table[,22:26],na.rm = TRUE),
                           Oct=rowSums(phylum_table[,27:40],na.rm = TRUE),
                           Nov=rowSums(phylum_table[,41:51],na.rm = TRUE),
                           Dec=rowSums(phylum_table[,52:56],na.rm = TRUE),
                           Jan=rowSums(phylum_table[,57:66],na.rm = TRUE),
                           Feb=rowSums(phylum_table[,67:71],na.rm = TRUE),
                           Mar=rowSums(phylum_table[,72:76],na.rm = TRUE),
                           Apr=rowSums(phylum_table[,77:86],na.rm = TRUE),
                           May=rowSums(phylum_table[,87:96],na.rm = TRUE),
                           Jun=rowSums(phylum_table[,97:101],na.rm = TRUE))
x <- rownames(phylum_group)
phylum_group <- sweep(phylum_group,2,colSums(phylum_group),"/")
rownames(phylum_group) <- x
phylum_group <- phylum_group[order(rowSums(phylum_group),decreasing = T),]
phylum_group <- phylum_group[c(1:3),]
phylum_group["Other",] <- (1-colSums(phylum_group))
x <- rownames(phylum_group)
phylum_group <- cbind(rownames(phylum_group),phylum_group)
phylum_group <- melt(phylum_group)
colnames(phylum_group) <- c("phylum","months","abundance")
phylum_group$phylum <- factor(phylum_group$phylum,levels = rev(x))
phylum_group$months <- factor(phylum_group$months,levels = c("Jan","Feb","Mar","Apr","May","Jun",
  "Jul","Aug","Sep","Oct","Nov","Dec"))
P9 <- ggplot(phylum_group,aes(months,abundance,fill=phylum))+geom_bar(stat = "identity")+
  scale_fill_manual(values = c( "#e5cf0d","#077535","#fee3ce", "#a8dde1"))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45))
P9
family_table <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/chouping/family_table.txt",row.names = 1)
family <- sweep(family_table,2,colSums(family_table),"/")
family <- family[order(rowSums(family),decreasing = T),]
family <- family[1:10,]
stat <- data.frame(row.names = rownames(family),mean=apply(family,1,mean),sd=apply(family,1,sd))
family["Other",] <- (1-colSums(family))
x <- rownames(family)
family <- cbind(rownames(family),family)
family <- melt(family)
colnames(family) <- c("family","samples","abundance")
family$family <- factor(family$family,levels = rev(x))
family$samples <- ifelse(grepl("^A", family$samples), "Summer",
                         ifelse(grepl("^B", family$samples), "Autumn",
                                ifelse(grepl("^C", family$samples), "Winter",
                                       ifelse(grepl("^D", family$samples), "Spring", family$samples))))
family$samples <- factor(family$samples, levels = c("Summer", "Autumn", "Winter", "Spring"))
family$abundance <- family$abundance / 25
p10 <- ggplot(family,aes(samples,abundance,fill=family))+geom_bar(stat = "identity",width = .75)+
  scale_fill_manual(values = c("#B6b3d6","#F2F3C7","#d87a80","#B2E2DC","#CFBACF","#FACE37","#E9CDDF",
                    "#9DC3E7","#AED96A","#D097BA","#5ab1ef","#f5e8f1","#aecce4"))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
p10
family_table <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/chouping/family_table.txt")
family_group <- data.frame(row.names=family_table$family,
                           Jul=rowSums(family_table[,2:11],na.rm = TRUE),
                           Aug=rowSums(family_table[,12:21],na.rm = TRUE),
                           Sep=rowSums(family_table[,22:26],na.rm = TRUE),
                           Oct=rowSums(family_table[,27:40],na.rm = TRUE),
                           Nov=rowSums(family_table[,41:51],na.rm = TRUE),
                           Dec=rowSums(family_table[,52:56],na.rm = TRUE),
                           Jan=rowSums(family_table[,57:66],na.rm = TRUE),
                           Feb=rowSums(family_table[,67:71],na.rm = TRUE),
                           Mar=rowSums(family_table[,72:76],na.rm = TRUE),
                           Apr=rowSums(family_table[,77:86],na.rm = TRUE),
                           May=rowSums(family_table[,87:96],na.rm = TRUE),
                           Jun=rowSums(family_table[,97:101],na.rm = TRUE))
x <- rownames(family_group)
family_group <- sweep(family_group,2,colSums(family_group),"/")
rownames(family_group) <- x
family_group <- family_group[order(rowSums(family_group),decreasing = T),]
family_group <- family_group[c(1:10),]
family_group["Other",] <- (1-colSums(family_group))
x <- rownames(family_group)
family_group <- cbind(rownames(family_group),family_group)
family_group <- melt(family_group)
colnames(family_group) <- c("family","months","abundance")
family_group$family <- factor(family_group$family,levels = rev(x))
family_group$months <- factor(family_group$months,levels = c("Jul","Aug","Sep","Oct","Nov","Dec",
  "Jan","Feb","Mar","Apr","May","Jun"))
p11 <- ggplot(family_group,aes(months,abundance,fill=family))+geom_bar(stat = "identity")+
  scale_fill_manual(values = c( "#B6b3d6","#F2F3C7","#F0ACAD","#B2E2DC","#CFBACF","#FACE37","#E9CDDF",
                                "#9DC3E7","#AED96A","#D097BA","#5ab1ef","#f5e8f1","#aecce4"))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45))
p11
genus_table <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/chouping/genus_table.txt", row.names = 1)
genus <- sweep(genus_table,2,colSums(genus_table),"/")
genus <- genus[order(rowSums(genus),decreasing = T),]
genus <- genus[c(1:20),]
stat <- data.frame(row.names = rownames(genus),mean=apply(genus,1,mean),sd=apply(genus,1,sd))
genus["Other",] <- (1-colSums(genus))
x <- rownames(genus)
genus <- cbind(rownames(genus),genus)
genus <- melt(genus)
colnames(genus) <- c("genus","samples","abundance")
genus$genus <- factor(genus$genus,levels = rev(x))
genus$samples <- ifelse(grepl("^A", genus$samples), "Summer",
                        ifelse(grepl("^B", genus$samples), "Autumn",
                               ifelse(grepl("^C", genus$samples), "Winter",
                                      ifelse(grepl("^D", genus$samples), "Spring", genus$samples))))
genus$samples <- factor(genus$samples, levels = c("Summer", "Autumn", "Winter", "Spring"))
genus$abundance <- genus$abundance / 25
p12 <- ggplot(genus,aes(samples,abundance,fill=genus))+geom_bar(stat = "identity", width=0.8)+
  scale_fill_manual(values = c( "#B6b3d6","#F0F8FF","#F2F3C7","#F0ACAD","#B2E2DC",
                                "#CFBACF","#FACE37","#E9CDDF","#99bcac","#fde8d3",
                                "#c49892","#aeD96A","#D097BA","#ffa288","#cbe4b1",
                                "#d87a80","#5ab1ef","#E5B2BD","#f5e8f1","#ffd47d","#7bc0cd"))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))+
  guides(fill = guide_legend(ncol = 2))
p12
genus_table <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/chouping/genus_table.txt")
genus_group <- data.frame(row.names=genus_table$genus,
                          Jul=rowSums(genus_table[,2:11],na.rm = TRUE),
                          Aug=rowSums(genus_table[,12:21],na.rm = TRUE),
                          Sep=rowSums(genus_table[,22:26],na.rm = TRUE),
                          Oct=rowSums(genus_table[,27:40],na.rm = TRUE),
                          Nov=rowSums(genus_table[,41:51],na.rm = TRUE),
                          Dec=rowSums(genus_table[,52:56],na.rm = TRUE),
                          Jan=rowSums(genus_table[,57:66],na.rm = TRUE),
                          Feb=rowSums(genus_table[,67:71],na.rm = TRUE),
                          Mar=rowSums(genus_table[,72:76],na.rm = TRUE),
                          Apr=rowSums(genus_table[,77:86],na.rm = TRUE),
                          May=rowSums(genus_table[,87:96],na.rm = TRUE),
                          Jun=rowSums(genus_table[,97:101],na.rm = TRUE))
x <- rownames(genus_group)
genus_group <- sweep(genus_group,2,colSums(genus_group),"/")
rownames(genus_group) <- x
genus_group <- genus_group[order(rowSums(genus_group),decreasing = T),]
genus_group <- genus_group[c(1:20),]
genus_group["Other",] <- (1-colSums(genus_group))
x <- rownames(genus_group)
genus_group <- cbind(rownames(genus_group),genus_group)
genus_group <- melt(genus_group)
colnames(genus_group) <- c("genus","months","abundance")
genus_group$genus <- factor(genus_group$genus,levels = rev(x))
genus_group$months <- factor(genus_group$months,levels = c("Jan","Feb","Mar","Apr","May","Jun",
  "Jul","Aug","Sep","Oct","Nov","Dec"))
p13 <- ggplot(genus_group,aes(months,abundance,fill=genus))+geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#d1e6c1", "#Face37",  "#E15759", "#76B7B2", "#666666", "#F28E2B",
                               "#B07AA1", "#FF9DA7", "#8c564b",  "#BAB0AC", "#66A61E","#B6992D",
                               "#D37295",  "#B2DF8A", "#9467bd", "#5ab1ef", "#F2F3C7","#4E79A7",
                               "#C4A38a", "#86BCB6", "#f5e8f1" ))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45))
p13
species_table <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/chouping/species_table.txt", row.names = 1)
species <- sweep(species_table,2,colSums(species_table),"/")
species <- species[order(rowSums(species),decreasing = T),]
species <- species[c(1:20),]
stat <- data.frame(row.names = rownames(species),mean=apply(species,1,mean),sd=apply(species,1,sd))
species["Other",] <- (1-colSums(species))
x <- rownames(species)
species <- cbind(rownames(species),species)
species <- melt(species)
colnames(species) <- c("species","samples","abundance")
species$species <- factor(species$species,levels = rev(x))
species$samples <- ifelse(grepl("^A", species$samples), "Summer",
                          ifelse(grepl("^B", species$samples), "Autumn",
                                 ifelse(grepl("^C", species$samples), "Winter",
                                       ifelse(grepl("^D", species$samples), "Spring", species$samples))))
species$samples <- factor(species$samples, levels = c("Summer", "Autumn", "Winter", "Spring"))
species$abundance <-species$abundance / 25
p14 <- ggplot(species,aes(samples,abundance,fill=species))+geom_bar(stat = "identity", width=0.8)+
  scale_fill_manual(values = c( "#B6b3d6","#F0F8FF","#F2F3C7","#F0ACAD","#B2E2DC","#CFBACF","#FACE37","#E9CDDF","#99bcac","#fde8d3",
                                "#c49892","#aeD96A","#D097BA","#ffa288","#cbe4b1","#d87a80","#5ab1ef","#E5B2BD","#f5e8f1","#ffd47d","#7bc0cd"))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))+
  guides(fill = guide_legend(ncol = 2))
p14
plot_grid(P9, p13)
species_table <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/chouping/species_table1.txt")
species_group <- data.frame(row.names=species_table$Species,
                            Jul=rowSums(species_table[,2:11],na.rm = TRUE),
                            Aug=rowSums(species_table[,12:21],na.rm = TRUE),
                            Sep=rowSums(species_table[,22:26],na.rm = TRUE),
                            Oct=rowSums(species_table[,27:40],na.rm = TRUE),
                            Nov=rowSums(species_table[,41:51],na.rm = TRUE),
                            Dec=rowSums(species_table[,52:56],na.rm = TRUE),
                            Jan=rowSums(species_table[,57:66],na.rm = TRUE),
                            Feb=rowSums(species_table[,67:71],na.rm = TRUE),
                            Mar=rowSums(species_table[,72:76],na.rm = TRUE),
                            Apr=rowSums(species_table[,77:86],na.rm = TRUE),
                            May=rowSums(species_table[,87:96],na.rm = TRUE),
                            Jun=rowSums(species_table[,97:101],na.rm = TRUE))
x <- rownames(species_group)
species_group <- sweep(species_group,2,colSums(species_group),"/")
rownames(species_group) <- x
species_group <- species_group[order(rowSums(species_group),decreasing = T),]
species_group <- species_group[c(1:20),]
species_group["Other",] <- (1-colSums(species_group))
x <- rownames(species_group)
species_group <- cbind(rownames(species_group),species_group)
species_group <- melt(species_group)
colnames(species_group) <- c("species","months","abundance")
species_group$species <- factor(species_group$species,levels = rev(x))
species_group$months <- factor(species_group$months,levels = c("Jan","Feb","Mar","Apr","May","Jun",
                                                               "Jul","Aug","Sep","Oct","Nov","Dec" ))
p15 <- ggplot(species_group,aes(months,abundance,fill=species))+geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#B2E2DC","#F0F8FF","#F2F3C7","#F0ACAD","#b72a01","#CFBACF","#FACE37","#E9CDDF","#99bcac","#fde8d3",
                               "#8c564b","#AED96A","#D097BA","#e5cf0d","#b6a2de","#d87a80","#5ab1ef","#E5B2BD","#66A61E","#DAD2FC","#aecce4"))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12, face = "italic"))
p15
core_genus <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/beiyong/core_genus.txt",row.names = 1)
core_genus <- sweep(core_genus,2,colSums(core_genus),"/")
core_genus <- core_genus[order(rowSums(core_genus),decreasing = T),]
core_genus <- core_genus[1:20,]
stat <- data.frame(row.names = rownames(core_genus),mean=apply(core_genus,1,mean),sd=apply(core_genus,1,sd))
core_genus["Other",] <- (1-colSums(core_genus))
x <- rownames(core_genus)
core_genus <- cbind(rownames(core_genus),core_genus)
core_genus <- melt(core_genus)
colnames(core_genus) <- c("core_genus","samples","abundance")
core_genus$core_genus <- factor(core_genus$core_genus,levels = rev(x))
core_genus$samples <- ifelse(grepl("^A", core_genus$samples), "Summer",
                             ifelse(grepl("^B", core_genus$samples), "Autumn",
                                    ifelse(grepl("^C", core_genus$samples), "Winter",
                                           ifelse(grepl("^D", core_genus$samples), "Spring", core_genus$samples))))
core_genus$samples <- factor(core_genus$samples, levels = c("Summer", "Autumn", "Winter", "Spring"))
core_genus$abundance <- core_genus$abundance / 25
p16 <- ggplot(core_genus,aes(samples,abundance,fill=core_genus))+geom_bar(stat = "identity",width = .75)+
  scale_fill_manual(values = c("#B2E2DC","#F0F8FF","#F2F3C7","#F0ACAD","#b72a10","#CFBACF","#FACE37","#E9CDDF","#99bcac","#fde8d3",
                               "#9DC3E7","#AED96A","#D097BA","#e5cf0d","#b6a2de","#d87a80","#5ab1ef","#E5B2BD","#f5e8f1","#DAD2FC","#aecce4"))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45))
p16
plot_grid(p8,p9,p10,p11,p12,p13,align = "hv",ncol = 2)
ggsave("tax_bar.pdf",width = 10,height = 10)
data <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/feature_table.txt", row.names = 1)
data <- data[,8:107]
data <- data.frame(Summer=rowSums(data[,1:25]),Autumn=rowSums(data[,26:50]),
                   Winter=rowSums(data[,51:75]),Spring=rowSums(data[,76:100]))
Spring <- data %>% filter(Spring != 0) %>% rownames()
Summer <- data %>% filter(Summer != 0) %>% rownames()
Autumn <- data %>% filter(Autumn != 0) %>% rownames()
Winter <- data %>% filter(Winter != 0) %>% rownames()
a <- list(`Spring` = c(Spring),
          `Summer` = c(Summer),
          `Autumn` = c(Autumn),
          `Winter` = c(Winter))
p1 <- ggvenn(a,show_elements = F,show_percentage = F,
       fill_color = c("#d1e6c1","#faf0b6","#e69a9b","#c0dbe5"),
       fill_alpha = 1,
       stroke_color = "white",stroke_alpha = 0.5,stroke_linetype = 'dashed',
       text_size = 6,set_name_size = 8)
p1
data <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/feature_table.txt",row.names = 1)
data <- sweep(data[,8:107],2,colSums(data[,8:107]),"/")
pathgen <- data[grep("Animal Pathogen",rownames(data)),]
pathgen <- data.frame(colSums(pathgen),c(rep("Spring",25),rep("Summer",25),rep("Autumn",25),rep("Winter",25)))
colnames(pathgen) <- c("Relative_abundace","Seasons")
pathgen$Seasons <- factor(pathgen$Seasons,levels = c("Spring","Summer","Autumn","Winter"))
p2 <- ggplot(pathgen,aes(Seasons,Relative_abundace,fill=Seasons))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=1,color="black",size=0.5)+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",width = 0.15,position = position_dodge( .9))+
  stat_compare_means(method = "kruskal")+
  scale_fill_manual(values = c("#d1e6c1","#faf0b6","#e69a9b","#c0dbe5"))+ylab("Animal_Pathogen_Relative_abundace")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA,linewidth = .5))+guides(fill="none")
p2
dunnTest(as.numeric(Relative_abundace)~Seasons,pathgen,method="bh")
plot_grid(p1,p2)
ggsave("D:/Desktop/seq/ITS/\u62bd\u5e73/venn+guild.pdf",width = 8,height = 6)
data <- read.csv("D:/Desktop/seq/ITS/\u62bd\u5e73/MEAN/Mar_node_attribute_for_Gephi].csv",sep = "\t")
phylum <- read.delim("D:/Desktop/seq/ITS/\u62bd\u5e73/feature_table.txt")
data <- cbind(data,phylum[data$ID,2])
colnames(data)[12] <- "phylum"
data$phylum[is.na(data$phylum)] <- "Others"
for (i in 1:nrow(data)) {
  if (data$phylum[i] != "p__Ascomycota" & data$phylum[i] != "p__Basidiomycota"){
    data$phylum[i] <- "Others"
  }
}
write.csv(data,"D:/Desktop/seq/ITS/\u62bd\u5e73/MEAN/Mar_node.csv",row.names = F,quote = F)
install.packages("microeco")
library(microeco)
library(tidytree)
library(pheatmap)
library(magrittr)
library(ggtree)
library(igraph)
data <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/feature_table.txt",row.names = 1)
data <- sweep(data,2,colSums(data),"/")
data <- sweep(data[,8:107],2,colSums(data[,8:107]),"/")
data <- data[c("g__Mrakia","g__Wickerhamiella","g__Beauveria","g__Cladophialophora",
               "g__Phaeosphaeria","g__Ramularia","g__Saccharomycopsis","g__Dimorphospora",
               "g__Hymenoscyphus","g__Tausonia","g__Cryptococcus","g__Barnettozyma","g__Apiotrichum",
               "g__Rachicladosporium","g__Pseudopeziza","g__Alternaria","g__Malassezia","g__Sporormiella",
               "g__Ceratobasidiaceae_gen_Incertae_sedis","g__Papiliotrema"),]
annotation <- data.frame(row.names=colnames(data),months=c(rep("Dec",10),rep("Jan",10),rep("Feb",10),rep("Mar",10)))
pheatmap(log(data+0.0001),fondsize=10,scale ="row",legend_breaks = c(-1, 0, 1),cluster_rows = F,cluster_cols = F,border_color=NA,
         color = colorRampPalette(colors = c("#5ab1ef","white","#ffb980"))(100),
         angle_col = 45,cellwidth = 15,cellheight = 8,gaps_col = c(10,20,30),gaps_row = c(4,10,13),
         annotation_col = annotation)
asv <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/asv.txt", sep="\t", header=T, check.names = FALSE, row.names = 1)
tax <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/asv.txt", sep="\t", header=T, check.names = FALSE, row.names = 1)
sample <-read.delim("C:/Users/xuzy0/Desktop/xzy/split/sample.txt", sep="\t", header=T, check.names = FALSE, row.names = 1)
tax %<>% tidy_taxonomy
df <- microtable$new(sample_table = sample,
                     otu_table = asv,
                     tax_table = tax,
                     auto_tidy = F)
df
lefse<-trans_diff$new(dataset = df,
                      method = "lefse",
                      group = "group")
lefse$plot_diff_bar(threshold=4,
                    width=0.8,
                    group_order = c("spring", "summer", "autumn", "winter"))
data <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/feature_table_relative.txt",row.names = 1)
data <- data[,8:107]
for (i in nrow(data):1) {
  if (sum(data[i,c(1:25)]) == 0 | sum(data[i,c(26:50)]) == 0 |
      sum(data[i,c(51:75)]) == 0 |sum(data[i,c(76:100)]) == 0){
    data <- data[-i,]
  }
}
annotation <- data.frame(row.names=colnames(data),seasons=c(rep("Summer",25),rep("Autumn",25),rep("Winter",25),rep("Spring",25)))
pheatmap(log(data+0.0001),fondsize=10,scale ="row",legend_breaks = c(-1, 0, 1),cluster_rows = T,cluster_cols = F,border_color=NA,
         color = colorRampPalette(colors = c("navy", "white", "red"))(200),
         angle_col = 45,cellwidth = 5,cellheight = 2.5,gaps_col = c(25,50,75),show_rownames = F,show_colnames = F,
         annotation_col = annotation)
data_plot <- cbind(annotation,abundance=colSums(data))
ggplot(data_plot,aes(seasons,abundance,fill=seasons))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=1,color="black",size=0.3)+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",width = 0.15,position = position_dodge( .9))+
  stat_compare_means(method = "kruskal")+
  scale_fill_manual(values = c("#d1e6c1","#faf0b6","#e69a9b","#c0dbe5"))+ylab("Core_Relative_abundace")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA,linewidth = .5))+guides(fill="none")
mean(data_plot$abundance[1:25]);sd((data_plot$abundance[1:25]))
mean(data_plot$abundance[26:50]);sd((data_plot$abundance[26:50]))
mean(data_plot$abundance[51:75]);sd((data_plot$abundance[51:75]))
mean(data_plot$abundance[76:100]);sd((data_plot$abundance[76:100]))
dunnTest(as.numeric(abundance)~seasons,data_plot,method="bh")
library(FSA)
data <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/chouping/genus_table.txt",row.names = 1)
data <- sweep(data,2,colSums(data),"/")
data <- data[order(rowSums(data),decreasing = T),]
data <- data.frame(t(data[c(2:22),]),seasons=c(rep("Summer",25),rep("Autumn",25),rep("Winter",25),rep("Spring",25)) )
data$seasons <- factor(data$seasons,levels = c("Summer","Autumn","Winter","Spring"))
colnames(data)
p1 <- ggplot(data, aes(x=seasons, y=as.numeric(g__Gautieria),fill=seasons)) +
  geom_boxplot() + stat_compare_means(method = "kruskal")+xlab("")+ylab("abundace")+
  scale_fill_manual(values = c("#d1e6c1","#faf0b6","#e69a9b","#c0dbe5"))+
  ggtitle("g__Gautieria")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA,size = .5))+guides(fill="none")
p2 <- ggplot(data, aes(x=seasons, y=as.numeric(g__Cortinarius),fill=seasons)) +
  geom_boxplot() + stat_compare_means(method = "kruskal")+xlab("")+ylab("abundace")+
  scale_fill_manual(values = c("#d1e6c1","#faf0b6","#e69a9b","#c0dbe5"))+
  ggtitle("g__Cortinarius")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA,size = .5))+guides(fill="none")
p3<- ggplot(data, aes(x=seasons, y=as.numeric(g__Pholiota),fill=seasons)) +
  geom_boxplot() + stat_compare_means(method = "kruskal")+xlab("")+ylab("abundace")+
  scale_fill_manual(values = c("#d1e6c1","#faf0b6","#e69a9b","#c0dbe5"))+
  ggtitle("g__Pholiota")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA,size = .5))+guides(fill="none")
p4<- ggplot(data, aes(x=seasons, y=as.numeric(g__Tausonia),fill=seasons)) +
  geom_boxplot() + stat_compare_means(method = "kruskal")+xlab("")+ylab("abundace")+
  scale_fill_manual(values = c("#d1e6c1","#faf0b6","#e69a9b","#c0dbe5"))+
  ggtitle("g__Tausonia")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA,size = .5))+guides(fill="none")
plot_grid(p1,p2,p3,p4,rows = 2,align = "hv")
dunnTest(as.numeric(p__Ascomycota)~months,data,method="bh")
dunnTest(as.numeric(p__Basidiomycota)~months,data,method="bh")
dunnTest(as.numeric(f__Mrakiaceae)~months,data,method="bh")
dunnTest(as.numeric(f__Bulleribasidiaceae)~months,data,method="bh")
dunnTest(as.numeric(f__Nectriaceae)~months,data,method="bh")
plot_grid(p1,p2,p3,p4,p5,rows = 1,align = "hv")
data <- read.delim("D:/Desktop/seq/ITS/\u62bd\u5e73/feature_table.txt",row.names = 1)
data[,8:47] <- sweep(data[,8:47],2,colSums(data[,8:47]),"/")
data <- data[order(rowSums(data[,8:47]),decreasing = T),]
write.csv(data,"D:/Desktop/feature_relative.csv",row.names = T)




