library("dada2")

for(i in c(
    5,
    7,
    8,
    10,
    11,
    13,
    14,
    15,
    16,
    19,
    21,
    22,
    23,
    24,
    49,
    50,
    51,
    52,
    54,
    55,
    56,
    57,
    58,
    59,
    60,
    61,
    62,
    63
    )) {
  barcode = i
  if (barcode < 10) {
    barcode = paste("0", i, sep="")
  }
  print(paste("Analyzing barcode #:", barcode))
path <- paste("/Users/g/Desktop/Classes_24-25/Biotech/BIOL_28_Bioinformatics/eDNA/Python-Programs/BioTech/Bioinformatics/Chaparral/wildfirePlants-DNA-nanopore-sequence/fastq_pass/barcode", barcode, sep="")
fnFs <- sort(list.files(path, pattern = ".fastq", full.names = TRUE, recursive=TRUE))
FpairFITS<-"GGAAGTAAAAGTCGTAACAAGG...AACTTTCAACAACGGATCTCTTG"
cutadapt <- "/opt/anaconda3/envs/cutadaptenv/bin/cutadapt"
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-a", FpairFITS, "-o", fnFs.cut[i], fnFs[i], "--revcomp", "--discard-untrimmed"))}
cutFs <- sort(list.files(path.cut, pattern = ".fastq", full.names = TRUE))
filtFs <- file.path(path.cut, "filtered", basename(cutFs))

out <- filterAndTrim(cutFs, filtFs, minLen = 50, maxN = 0, maxEE = Inf, truncQ = 2, compress = TRUE, multithread = TRUE)
out
errF <- learnErrors(filtFs, multithread = TRUE)
get.sample.name <- function(fname) strsplit(basename(fname), "\\.")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
setDadaOpt(DETECT_SINGLETONS = TRUE)
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
names(dadaFs) <- sample.names
seqtab <- makeSequenceTable(dadaFs)
unite.ref <- "~/Desktop/Classes_24-25/Biotech/BIOL_28_Bioinformatics/eDNA/TAXONOMY_REFS/ITS_sh_general_release_dynamic_19.02.2025.fasta"
taxa <- assignTaxonomy(seqtab, unite.ref, multithread = TRUE, verbose = TRUE)
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
taxa.print
write.csv(taxa.print, file = paste("~/Desktop/Classes_24-25/Biotech/BIOL_28_Bioinformatics/eDNA/dada2_taxa_output/barcode", barcode,"_FITS_taxaPrint.csv", sep=""), row.names = TRUE)
}
# room primers: 16s, 18s, fits
# matt primers: 18s, fits, pits

# Pb_meta<-read.csv("lead_levels.csv", header = TRUE, sep = ",")
# rownames(Pb_meta)<-sample.names
# Pb_meta$environment<-"Safe"
# Pb_meta$environment[Pb_meta$lead_level_PPM>400]<-"Unhealthy"
# library("phyloseq")
# ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
#                sample_data(Pb_meta), 
#                tax_table(taxa))
# dna <- Biostrings::DNAStringSet(taxa_names(ps))
# names(dna) <- taxa_names(ps)
# ps <- merge_phyloseq(ps, dna)
# taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
# ps
# plot_richness(ps, x="lead_level_PPM", measures=c("Shannon"), color="environment")
# top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
# ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
# ps.top20 <- prune_taxa(top20, ps.top20)
# plot_bar(ps.top20, x="environment", fill="Phylum")
