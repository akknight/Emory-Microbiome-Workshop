# Emory-Microbiome-Workshop


AWS Workshop
8/16/19


Qiime2 Code
​
source activate qiime2-2019.7
​
qiime tools import --type SampleData[PairedEndSequencesWithQuality] --input-path Manifest.csv --output-path paired-end.qza —input-format PairedEndFastqManifestPhred33
​
qiime demux summarize --i-data paired-end.qza --o-visualization paired-end.qzv
​
qiime dada2 denoise-paired --i-demultiplexed-seqs paired-end.qza --o-table table.qza --o-representative-sequences rep-seqs.qza --p-trim-left-f 8 --p-trim-left-r 9 --p-trunc-len-f 150 --p-trunc-len-r 150 --p-n-threads 0 --o-denoising-stats denoisingstats.qza
​
qiime metadata tabulate --m-input-file denoisingstats.qza  --o-visualization stats-dada2.qzv
​
qiime feature-table summarize --i-table table.qza --o-visualization table.qzv  --m-sample-metadata-file Map.tsv

qiime feature-table tabulate-seqs  --i-data rep-seqs.qza  --o-visualization rep-seqs.qzv
​
​
qiime feature-classifier classify-sklearn  --i-classifier classifier.qza  --i-reads rep-seqs.qza --o-classification taxonomyV3V5.qza --p-n-jobs -1 --p-confidence .8

qiime metadata tabulate  --m-input-file taxonomyV3V5.qza  --o-visualization taxonomyV3V5.qzv
​
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza
​
qiime diversity core-metrics-phylogenetic  --i-phylogeny rooted-tree.qza  --i-table table.qza  --p-sampling-depth 723  --m-metadata-file Map.tsv  --output-dir core-metrics-results
​
qiime diversity alpha-group-significance  --i-alpha-diversity ./core-metrics-results/evenness_vector.qza  --m-metadata-file ./Map.tsv  --o-visualization ./core-metrics-results/evenness_statistics.qzv
​
​
qiime tools export --input-path taxonomy.qza   --output-path exported
qiime tools export --input-path table.qza   --output-path exported
biom convert –i feature-table.biom -o table.tsv --to-tsv
​
​
Phyloseq Code

install.packages("devtools")
install.packages("ggplot2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("phyloseq")
BiocManager::install("Biostrings")
BiocManager::install("plyr")

library("devtools")
library("phyloseq")
library("ggplot2")
library("Biostrings")
library("plyr")


#Phyloseq works in Objects

#Get sample datat:
data(GlobalPatterns)

#Look at the Object:
GlobalPatterns

#Examples to see data stored in object

head(sample_data(GlobalPatterns))
otu_table(GlobalPatterns)[1:5,1:5]
tax_table(GlobalPatterns)[1:5,]

#Plot alpha diversity

plot_richness(GlobalPatterns)

#Show only the ones you want:
GP = GlobalPatterns

plot_richness(GlobalPatterns, x="SampleType", measures=c("Chao1", "Shannon")) + geom_boxplot()
plot_richness(GP, x="SampleType", measures=c("Chao1", "Shannon"))

#Make a human category:

sample_data(GP)$human <- get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
plot_richness(GP, x="human", color="SampleType", measures=c("Chao1", "Shannon"))

#Bar Charts (make me sad...)

#for speed, plot top twenty

TopNOTUs = names(sort(taxa_sums(GP), TRUE)[1:20])
t20 = prune_taxa(TopNOTUs, GP)
plot_bar(t20, fill="Family")

p = plot_bar(t20, "Genus", fill="Genus", facet_grid=sample_data(t20)$human )
p + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
#Ordination Plots

#Remove OTUs that do not show appear more than 5 times in more than half the samples

wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
GP1 = prune_taxa(wh0, GP)

#Transform to even sampling depth (we can argue about this later).

GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))

#Keep only the most abundant five phyla.

phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
GP1 = prune_taxa((tax_table(GP1)[, "Phylum"] %in% top5phyla), GP1)

GP.ord <- ordinate(GP1, "NMDS", "bray")
p1 <- plot_ordination(GP1, GP.ord, type="taxa", color="Phylum", title="taxa")
print (p1)

#Break out into Taxa

p1 + facet_wrap(~Phylum, 3)

#Let's look at the samples by type:

p2 = plot_ordination(GP1, GP.ord, type="samples", color="SampleType") 
p2 + geom_polygon(aes(fill=SampleType)) + geom_point(size=5) + ggtitle("samples")

#MDS plot

ordu = ordinate(GP1, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(GP1, ordu, color="SampleType", )
