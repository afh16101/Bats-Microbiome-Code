# Bats-Microbiome-Code

setwd("C:/Users/ghost/Desktop/Ph.D.Program/Dr. Graf Lab/BAT/DADA2/")

# Filtering and inferring sequence variants must be performed per sequence run
library(dada2); packageVersion("dada2")

# File parsing
pathF <- "path/to/FWD" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "path/to/REV" # CHANGE ME ...
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(240,200), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

# File parsing
filtpathF <- "path/to/FWD/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "path/to/FWD/filtered" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# Merge multiple runs (if necessary)
st1 <- readRDS("C:/Users/ghost/Desktop/Ph.D.Program/Dr. Graf Lab/BAT/DADA2/Run1_2017_11_13/seqtab.rds")
st2 <- readRDS("C:/Users/ghost/Desktop/Ph.D.Program/Dr. Graf Lab/BAT/DADA2/Run2_2018_02_12/seqtab.rds")
st3 <- readRDS("C:/Users/ghost/Desktop/Ph.D.Program/Dr. Graf Lab/BAT/DADA2/Run3_4_2018_03_16/seqtab.rds")
st4 <- readRDS("C:/Users/ghost/Desktop/Ph.D.Program/Dr. Graf Lab/BAT/DADA2/Run5_2018_06_08/seqtab.rds")
st5 <- readRDS("C:/Users/ghost/Desktop/Ph.D.Program/Dr. Graf Lab/BAT/DADA2/Run6_2018_06_15/seqtab.rds")
st6 <- readRDS("C:/Users/ghost/Desktop/Ph.D.Program/Dr. Graf Lab/BAT/DADA2/Run7_2018_08_13/seqtab.rds")
st7 <- readRDS("C:/Users/ghost/Desktop/Ph.D.Program/Dr. Graf Lab/BAT/DADA2/Run8_2018_09_25/seqtab.rds")
st8 <- readRDS("C:/Users/ghost/Desktop/Ph.D.Program/Dr. Graf Lab/BAT/DADA2/Run9_2018_11_07/seqtab.rds")
st.all <- mergeSequenceTables(st1, st2, st3, st4, st5, st6, st7, st8)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/ghost/Desktop/Ph.D.Program/Dr. Graf Lab/BAT/DADA2/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "C:/Users/ghost/Desktop/Ph.D. Program/Dr. Graf Lab/BAT/DADA2/silva_species_assignment_v128.fa.gz")

# Making phylogenetic tree #
library(dada2)
library(phangorn)
library(DECIPHER)

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align) ## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0)) # this step takes a while (minutes to 24hr based on quantity of data)
detach("package:phangorn", unload=TRUE)


# Uploading into phyloseq
library(phyloseq); packageVersion("phyloseq")

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa), 
               phy_tree(tree))

# Remove samples/taxa below certain number of reads #
ps <- prune_samples(sample_sums(ps)>=10, ps)

# decontam #
library(phyloseq)
library(ggplot2)
library(decontam)

# Freqyency method
contam.freq <- isContaminant(ps, method = "frequency", conc = "Qubit")
table(contam.freq$contaminant)
head(which(contam.freq$contaminant))
plot_frequency(ps, taxa_names(ps)[c(1,93)], conc="Qubit") + 
  xlab("DNA Concentration")
set.seed(101)
plot_frequency(ps, taxa_names(ps)[sample(which(contam.freq$contaminant),10)], conc="Qubit") +
  xlab("DNA Concentration")
ps.noncontam <- prune_taxa(!contam.freq$contaminant, ps)


# Create table, number of features for each phyla
table(tax_table(ps.noncontam)[, "Phylum"], exclude = NULL)

# Removing NA and uncharacterized features
ps0 <- subset_taxa(ps.noncontam, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Rarefying phyloseq object
library(phyloseq)
set.seed(5)

ps0.1k <- rarefy_even_depth(ps0, sample.size = 1000,
                                 rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps0.5k <- rarefy_even_depth(ps0, sample.size = 5000,
                          rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps0.10k <- rarefy_even_depth(ps0, sample.size = 10000,
                             rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Subsetting phyloseq object
library(phyloseq)

ps0.1k.o <- subset_samples(ps0.1k, Body_Site=="oral")
ps0.1k.r <- subset_samples(ps0.1k, Body_Site=="rectal")
Aj.oral <- subset_samples(ps0.1k.o, Sample_Species=="Artibeus_jamaicensis")
Aj.rectal <- subset_samples(ps0.1k.r, Sample_Species=="Artibeus_jamaicensis")
Bc.oral <- subset_samples(ps0.1k.o, Sample_Species=="Brachyphylla_cavernarum")
Bc.rectal <- subset_samples(ps0.1k.r, Sample_Species=="Brachyphylla_cavernarum")
Es.oral <- subset_samples(ps0.1k.o, Sample_Species=="Erophylla_sezekorni")
Es.rectal <- subset_samples(ps0.1k.r, Sample_Species=="Erophylla_sezekorni")
Mr.oral <- subset_samples(ps0.1k.o, Sample_Species=="Monophyllus_redmani")
Mr.rectal <- subset_samples(ps0.1k.r, Sample_Species=="Monophyllus_redmani")
Nl.oral <- subset_samples(ps0.1k.o, Sample_Species=="Noctilio_leporinus")
Nl.rectal <- subset_samples(ps0.1k.r, Sample_Species=="Noctilio_leporinus")
Ef.oral <- subset_samples(ps0.1k.o, Sample_Species=="Eptesicus_fuscus")
Ef.rectal <- subset_samples(ps0.1k.r, Sample_Species=="Eptesicus_fuscus")
Mb.oral <- subset_samples(ps0.1k.o, Sample_Species=="Mormoops_blainvillii")
Mb.rectal <- subset_samples(ps0.1k.r, Sample_Species=="Mormoops_blainvillii")
Pq.oral <- subset_samples(ps0.1k.o, Sample_Species=="Pteronotus_quadridens")
Pq.rectal <- subset_samples(ps0.1k.r, Sample_Species=="Pteronotus_quadridens")
Pp.oral <- subset_samples(ps0.1k.o, Sample_Species=="Pteronotus_portoricensis")
Pp.rectal <- subset_samples(ps0.1k.r, Sample_Species=="Pteronotus_portoricensis")


#### Standard OTU format table from dada2/phyloseq ####
library(phyloseq)

ASV.table <- Mb.rectal # save original phyloseq object because you are going to overwrite the original and change its format

# This changes the header from the actual sequence to Seq_001, Seq_002 etc
taxa_names(ASV.table)
n_seqs <- seq(ntaxa(ASV.table))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(ASV.table) <- paste("Seq", formatC(n_seqs, 
                                       width = len_n_seqs, 
                                       flag = "0"), sep = "_")

# A possible way to get taxonomy included into the header is the following (continuing from above):
  
# Generate a vector containing the full taxonomy path for all ASVs
wholetax <- do.call(paste, c(as.data.frame(tax_table(ASV.table))
                               [c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")], 
                               sep = ","))  # to distinguish from "_" within tax ranks

# Turn the otu_table into a data.frame
otu_export <- as.data.frame(otu_table(ASV.table))
tmp <- names(otu_export)

# Paste wholetax and OTU_ids together
for(i in 1:length(tmp)){
  names(tmp)[i] = paste(wholetax[i], tmp[i], sep = ",")
}

# Overwrite old names
names(otu_export) <- names(tmp)

head(otu_export)[5]
otu_export <- t(otu_export)
write.csv(otu_export, file = "ASV_tables/Mb_rectal.csv")


##### Rarefaction Curve ####

library(phyloseq)
library(vegan)
library(ggplot2)

col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
pars <- expand.grid(col = col, stringsAsFactors = FALSE)
head(pars)

ps0.rarecurve <- with(pars,rarecurve(otu_table(ps0.5k), sample = 5000 ,step = 10, col = col,label = FALSE))


########## Bat NMDS ##########
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
theme_set(theme_bw())

bray <-"bray" # define metric for analysis

distOrd = phyloseq::distance(ps0.1k.r, method = c(bray)) # calculate distances

ord = ordinate(ps0.1k.r, method = "NMDS", ddistance = distOrd) # calculate ordination

# NMDS. May want to change color= , shape= , alpha= , and size=

nmds.Geo <-plot_ordination(ps0.1k.r, ord, color = "Sample_Species") +
  
  geom_point(size=6,mapping = aes(color=Sample_Species)) +
  
  stat_ellipse(type = "t", level = 0.95, linetype = 2)+
  
  ggtitle("Bray Curtis")+
  
  theme(plot.title = element_text(hjust = 0.5))

nmds.Geo


### vegan::adonis ###

library(phyloseq)
library(vegan)

distOrd = phyloseq::distance(Aj.oral, "wunifrac") # metrics used "bray", "Jaccard", "unifrac", "wunifrac"

distOrd2 = phyloseq::distance(Aj.rectal, "wunifrac") # metrics used "bray", "Jaccard", "unifrac", "wunifrac"

adonis(distOrd ~ Sex, as(sample_data(Aj.oral), "data.frame"))
adonis(distOrd2 ~ Sex, as(sample_data(Aj.rectal), "data.frame"))
adonis(distOrd ~ Cave, as(sample_data(Aj.oral), "data.frame"))
adonis(distOrd2 ~ Cave, as(sample_data(Aj.rectal), "data.frame"))

adonis(distOrd ~ Sex, as(sample_data(Bc.oral), "data.frame"))
adonis(distOrd2 ~ Sex, as(sample_data(Bc.rectal), "data.frame"))
adonis(distOrd ~ Cave, as(sample_data(Bc.oral), "data.frame"))
adonis(distOrd2 ~ Cave, as(sample_data(Bc.rectal), "data.frame"))

adonis(distOrd ~ Diet, as(sample_data(ps0.1k.o), "data.frame"))
adonis(distOrd2 ~ Diet, as(sample_data(ps0.1k.r), "data.frame"))

adonis(distOrd ~ Sample_Species, as(sample_data(ps0.1k.o), "data.frame"))
adonis(distOrd2 ~ Sample_Species, as(sample_data(ps0.1k.r), "data.frame"))


# Differential Abundance
library(DESeq2)
library(dplyr)
library(phyloseq)
library("ggplot2")
theme_set(theme_bw())

Aj.Bc.o <- subset_samples(ps0.o, Sample_Species %in% c("Artibeus_jamaicensis","Brachyphylla_cavernarum"))
Aj.Bc.r <- subset_samples(ps0.r, Sample_Species %in% c("Artibeus_jamaicensis","Brachyphylla_cavernarum"))

Aj.Es.o <- subset_samples(ps0.o, Sample_Species %in% c("Artibeus_jamaicensis","Erophylla_sezekorni"))
Aj.Es.r <- subset_samples(ps0.r, Sample_Species %in% c("Artibeus_jamaicensis","Erophylla_sezekorni"))

Aj.Mr.o <- subset_samples(ps0.o, Sample_Species %in% c("Artibeus_jamaicensis","Monophyllus_redmani"))
Aj.Mr.r <- subset_samples(ps0.r, Sample_Species %in% c("Artibeus_jamaicensis","Monophyllus_redmani"))

Aj.Nl.o <- subset_samples(ps0.o, Sample_Species %in% c("Artibeus_jamaicensis","Noctilio_leporinus"))
Aj.Nl.r <- subset_samples(ps0.r, Sample_Species %in% c("Artibeus_jamaicensis","Noctilio_leporinus"))

Aj.Ef.o <- subset_samples(ps0.o, Sample_Species %in% c("Artibeus_jamaicensis","Erophylla_sezekorni"))
Aj.Ef.r <- subset_samples(ps0.r, Sample_Species %in% c("Artibeus_jamaicensis","Erophylla_sezekorni"))

Aj.Mb.o <- subset_samples(ps0.o, Sample_Species %in% c("Artibeus_jamaicensis","Mormoops_blainvillii"))
Aj.Mb.r <- subset_samples(ps0.r, Sample_Species %in% c("Artibeus_jamaicensis","Mormoops_blainvillii"))

Aj.Pp.o <- subset_samples(ps0.o, Sample_Species %in% c("Artibeus_jamaicensis","Pteronotus_portoricensis"))
Aj.Pp.r <- subset_samples(ps0.r, Sample_Species %in% c("Artibeus_jamaicensis","Pteronotus_portoricensis"))

Aj.Pq.o <- subset_samples(ps0.o, Sample_Species %in% c("Artibeus_jamaicensis","Pteronotus_quadridens"))
Aj.Pq.r <- subset_samples(ps0.r, Sample_Species %in% c("Artibeus_jamaicensis","Pteronotus_quadridens"))


diagdds = phyloseq_to_deseq2(Aj.Pq.r, ~ Sample_Species)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric", sfType = c("poscounts"))
# Results of DESeq2
res = results(diagdds, cooksCutoff = FALSE)
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Aj.Pq.r)[rownames(sigtab), ], "matrix"))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

p.r.8 <- ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  ggtitle("P. quadridens vs A. jamaicensis")

# Combine plots
library(ggpubr)

ggarrange(p.o.1,p.o.2,p.o.3,p.o.4,p.o.5,p.o.6,p.o.7,p.o.8,
          ncol = 3, nrow = 3, common.legend = F, legend = "right")
ggarrange(p.o.7,p.o.8,
          ncol = 2, nrow = 1, common.legend = F, legend = "right")
rectal.diff <-ggarrange(p.r.1,p.r.2,p.r.3,p.r.4,p.r.5,p.r.6,p.r.7,p.r.8,
          ncol = 3, nrow = 3, common.legend = F, legend = "right")
ggarrange(p.r.5,p.r.6,p.r.7,p.r.8,
          ncol = 2, nrow = 1, common.legend = F, legend = "right")

ggsave(rectal.diff, filename="figures/rectal_diff.tiff", dpi="retina",width=16,height=12,units="in") # Save FIGURE A to .eps file

### Core Microbiome ###
library(microbiome)
library(phyloseq)

Aguas.core <- core(Aguas1k.o, detection = 1/100, prevalence = 68/100) # oral 68%, rectal 71%
Aguas.core
plot_bar(Aguas.core, fill = "Genus")
write.csv(psmelt(Aguas1k.r), "Aguas1k.r.csv")

Mata.core <- core(Mata1k.o, detection = 1/100, prevalence = 62/100) # oral 62%, rectal 54%
Mata.core
plot_bar(Mata.core, fill = "Genus")
write.csv(psmelt(Mata1k.r), "Mata1k.r.csv")

Rio.core <- core(Rio1k.o, detection = 1/100, prevalence = 73/100) # oral 73%, rectal 56%
Rio.core
plot_bar(Rio.core, fill = "Genus")
write.csv(psmelt(Rio1k.r), "Rio1k.r.csv")

oral.core <- core(ps0.5k.o, detection = 1/100, prevalence = 50/100)
plot_bar(oral.core, fill = "Genus")
oral.core.abun <- core_abundance(oral.NMDS, detection = 1/100, prevalence = 64/100)
rectal.core <- core(ps0.1k.r, detection = 1/100, prevalence = 50/100)
plot_bar(rectal.core, fill = "Family")
rectal.core.abun <- core_abundance(rectal.NMDS, detection = 1/100, prevalence = 51/100)

Aj.core.o <- core(Aj.oral.5k, detection = 1/100, prevalence = 99.9/100)
plot_bar(Aj.core.o, fill = "Genus")
Aj.core.r <- core(Aj.rectal, detection = 1/100, prevalence = 88/100)
plot_bar(Aj.core.r, fill = "Family")

Bc.core.o <- core(Bc.oral, detection = 1/100, prevalence = 99.9/100)
plot_bar(Bc.core.o, fill = "Genus")
Bc.core.r <- core(Bc.rectal, detection = 1/100, prevalence = 50/100)
plot_bar(Bc.core.r, fill = "Family")

Es.core.o <- core(Es.oral, detection = 1/100, prevalence = 50/100)
plot_bar(Es.core.o, fill = "Genus")
Es.core.r <- core(Es.rectal, detection = 1/100, prevalence = 50/100)
plot_bar(Es.core.r, fill = "Genus")

Ef.core.o <- core(Ef.oral, detection = 1/100, prevalence = 99.9/100)
plot_bar(Ef.core.o, fill = "Genus")
Ef.core.r <- core(Ef.rectal, detection = 1/100, prevalence = 99.9/100)
plot_bar(Ef.core.r, fill = "Genus")

Nl.core.o <- core(Nl.oral, detection = 1/100, prevalence = 99.9/100)
plot_bar(Nl.core.o, fill = "Genus")
Nl.core.r <- core(Nl.rectal, detection = 1/100, prevalence = 72/100)
plot_bar(Nl.core.r, fill = "Genus")

Mr.core.o <- core(Mr.oral, detection = 1/100, prevalence = 40/100)
plot_bar(Mr.core.o, fill = "Genus")
Mr.core.r <- core(Mr.rectal, detection = 1/100, prevalence = 60/100)
plot_bar(Mr.core.r, fill = "Genus")

Mb.core.o <- core(Mb.oral, detection = 1/100, prevalence = 50/100)
plot_bar(Mb.core.o, fill = "Genus")
Mb.core.r <- core(Mb.rectal, detection = 1/100, prevalence = 20/100)
plot_bar(Mb.core.r, fill = "Genus")

Pp.core.o <- core(Pp.oral, detection = 1/100, prevalence = 40/100)
plot_bar(Pp.core.o, fill = "Genus")
Pp.core.r <- core(Pp.rectal, detection = 1/100, prevalence = 60/100)
plot_bar(Pp.core.r, fill = "Genus")

Pq.core.o <- core(Pq.oral, detection = 1/100, prevalence = 40/100)
plot_bar(Pq.core.o, fill = "Genus")
Pq.core.r <- core(Pq.rectal, detection = 1/100, prevalence = 60/100)
plot_bar(Pq.core.r, fill = "Genus")

Neg.core <- core(Neg, detection = 1/100, prevalence = 10/100)
plot_bar(Neg.core, fill = "Genus")


### Heatmap ###
library(phyloseq)
library(dplyr)
library(pheatmap)

# Make sure to use glomerated phyloseq object and have both otu and taxonomy table rows match
ps0.1k.g <- tax_glom(ps0.1k, taxrank = "Family", NArm = TRUE)
otu_table(ps0.1k.g) <- t(otu_table(ps0.1k.g))

top50 <- names(sort(taxa_sums(ps0.1k.g), decreasing=TRUE))[1:50]
ps.top50 <- transform_sample_counts(ps0.1k.g, function(OTU) OTU/sum(OTU))
ps.top50 <- prune_taxa(top50, ps.top50)

# Pull out tax labels at your designated taxonomic level,
# grab abundance data from phyloseq object using genfac as index. 
# This will be your matrix to plot the heatmap from with taxonomy as rows, abundance data as values, and sample name as columns
genfac = factor(tax_table(ps.top50)[, "Family"])
gentab = apply(otu_table(ps.top50), MARGIN = 2, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})

# These reordering steps are only if you want to order your samples in a specific fashion (vs hierarchal), they can be skipped if not 
melted_phy_table = psmelt(ps.top50)
reordered_melted_phyla = melted_phy_table[order(melted_phy_table$Phylum),]

# Creating row annotations. Creating a dataframe to map annotations, this just ends up with phylum as rownames paired with a genus in a single column. 
melted_unique = reordered_melted %>% distinct(Phylum, .keep_all = T)
phyla_frame = subset(melted_unique, select = c(Phylum, Family))
rownames(phyla_frame) = phyla_frame$Phylum
phyla_frame = subset(phyla_frame, select = -c(Family))

# Creating column annotations. Same story as above but for samples this time. I'm putting sample ID as the rowname with the disease state in the column. 
reordered_melted_unique_sample = reordered_melted %>% distinct(SampleID, .keep_all = T)
NAFLD_frame = subset(reordered_melted_unique_sample, select = c(Sample_Species,Diet,Sex,Cave,Body_Site ,SampleID))
rownames(NAFLD_frame) = NAFLD_frame$SampleID
NAFLD_frame = subset(NAFLD_frame, select = -c(SampleID))


reordered_melted_phyla_no_genus_dups = reordered_melted_phyla %>% distinct(Family, .keep_all = TRUE)
gentab_reorder = gentab[reordered_melted_phyla_no_genus_dups$Family, order(colnames(gentab))]

phyla_frame1 = data.frame(Phylum = reordered_melted_phyla_no_genus_dups$Phylum)
rownames(phyla_frame1) = reordered_melted_phyla_no_genus_dups$Family
phyla_frame

# Assigning colors to variables
ann_colors = list(
  Sample_Species = c("Artibeus_jamaicensis"="red","Brachyphylla_cavernarum"="turquoise","Erophylla_sezekorni"="#229954","Monophyllus_redmani"="#2E86C1","Noctilio_leporinus"="#884EA0","Mormoops_blainvillii"="#839192","Pteronotus_portoricensis"="#FC8D62","Pteronotus_quadridens"="black","Eptesicus_fuscus"="yellow"),
  Diet = c("Frugivore"="#f0e442","Euryphagous"="#0072b2","Nectarivore"="#D55e00","Insectivore"="#56B4E9","Piscivore"="#66C2A5"),
  Cave = c("AguasBuenas"="red","MataDePlatano"="blue","RioEncantado"="yellow")) # assign colors and values. If you do not want to assign by hand, remove 'annotation_colors = ann_colors' in plotting below

# Heatmap command with some added customization. I log transform the data with a pseudocount as this highlights differences better and scales the legend in a more useful way. 
pheatmap_final_sample_hierarchal = pheatmap(log10(gentab_reorder+1), cluster_rows = FALSE, 
                                            annotation_row = phyla_frame1,
                                            annotation_col = NAFLD_frame,
                                            angle_col = 315, fontsize_row = 6,
                                            fontsize_col = 6, border_color = NA,
                                            show_colnames = FALSE, annotation_colors = ann_colors)

ggsave(pheatmap_final_sample_hierarchal, filename="figures/heatmap.tiff", dpi="retina",width=13,height=8.5,units="in") # save figure as a tiff pic
