#-----MAKING ASV AND TAXA TABLES-----#

#using terminal to move fastq files into one folder
## do this from inside the folder where the folders holding the fastq files are
#(base) [lgschaer@smtechtm-pc TA_TPA_Enrichment_02052021]$ for FILE in */*.fastq.gz; do mv $FILE /home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_02052021/Fastqs; done

library(dada2)
library(csv)

#path to fastq files
path <- "/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/Fastqs/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)
head(sample.names)

#visualize quality profiles
plotQualityProfile(fnFs[2:5])           #forward reads
plotQualityProfile(fnRs[2:5])           #reverse reads

#place filtered files in "filtered", a subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

#standard filtering parameters: 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(250, 250),
                     maxN=0, maxEE=Inf, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft = 12) # On Windows set multithread=FALSE

head(out)
#View(out)

under.1k.reads <- out %>%
  as.data.frame() %>%
  rownames_to_column(var = "FileName") %>%
  filter(reads.in < 1000)
View(under.1k.reads)

no.reads <- out %>%
  as.data.frame() %>%
  rownames_to_column(var = "FileName") %>%
  filter(reads.out == 0)
head(no.reads)

read.stats <- out %>%
  as_tibble() %>%
  summarise(
    max.in = max(reads.in),
    min.in = min(reads.in),
    mean.in = mean(reads.in),
    median.in = median(reads.in),
    sum.in = sum(reads.in),
    count.in = sum(ifelse(reads.in > 0, 1, 0)),
    under1000.in = sum(ifelse(reads.in < 1000, 1, 0)),
    max.out = max(reads.out),
    min.out = min(reads.out),
    mean.out = mean(reads.out),
    median.out = median(reads.out),
    sum.out = sum(reads.out),
    count.out = sum(ifelse(reads.out > 0, 1, 0)),
    under1000.out = sum(ifelse(reads.out < 1000, 1, 0))
  )
View(t(read.stats))

#save csvs
write.csv(out, "/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/out.csv")
write.csv(read.stats, "/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/read.stats.csv")

#learn about error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plot errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

#name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#sample inferance
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#inspecting dada-class object
dadaFs[[1]]
dadaRs[[1]]

#merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#inspect the merger data.frame from the first sample
head(mergers[[1]])

#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/home/lgschaer/old/Silva_Versions/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#save sequence table and taxa table
write_rds(seqtab.nochim, "/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/seqtab.rds")     #sequence table
write_rds(taxa.print, "/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/taxa.rds")          #taxa table

