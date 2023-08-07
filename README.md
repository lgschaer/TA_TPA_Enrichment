# Code and output files for the following paper:

# Microorganisms able to degrade chemically deconstructed plastic model compounds are widespread in the environment

Polyethylene terephthalate (PET) is one of the most common plastics worldwide, with millions of tons produced each year. Only ~30% of PET is currently recycled in the United States, leaving much room for improvement in current PET recycling practices. Hydrolysis is a well-known chemical reaction which can depolymerize PET into the aromatic molecule terephthalate, which can be biodegraded by many known microorganisms. However, hydrolysis produces a product with high osmolarity, which may be inhibitory to optimal microbial growth. Aminolysis is another chemical reaction which depolymerizes PET resulting in a product with lower osmolarity, however, this reaction produces terephthalamide, an aromatic molecule which is thought to be antimicrobial.
In this study we co-cultured sediments from five unique environments with either terephthalate or terephthalamide and performed biweekly transfers to fresh media and substrate. We use 16S rRNA sequencing to identify the dominant taxa in the enrichment cultures which may have terephthalate or terephthalamide-degrading metabolisms and compare them to the control enrichments. The goals of this study are to evaluate (1) whether terephthalamide is biodegradable and identify microorganisms able to degrade it, and (2) determine how widespread terephthalate and terephthalamide degrading metabolisms are in natural environments.


## Scripts


#### Raw sequencing file processing: dada2_TPATAEnrich_09062021.R

#### Diversity & community composition analysis: phyloseq_TPATAENrich_06242022.R

#### Differential abundance analysis: deseq2_10012021.R


## Data


### Dada2

#### Sequence table: seqtab.rds

#### Taxa table: taxa.rds

#### Metadata: TA_TPA_Metadata_09072021.csv


### Phyloseq


#### Mitochondria, chloroplasts and inoculum removed, subset to only include T4: ps_nomito_nochloro_noinocula_T4.rds

#### Unrarefied ps_unrarefied_09072021

#### Rarefied, chloroplasts removed tatpa_rarefied_nochloroplasts_rds


### DESeq2

#### folder dsq_out

#### enriched ASVs with taxonomy: deseq_enriched_w_tax_01062023.csv
