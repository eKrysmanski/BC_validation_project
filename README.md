# Project: BAIT-CAPTURE ENRICHEMENT VALIDATION (BC_validation)

**The Data:** The dataset we will be using consists of both conventional RNA-sequencing data and bait-capture target-enriched RNA-sequencing data. Conventional RNA-sequencing data were collected from gill, liver, anterior intestine, and trunk kidney samples, consisting of 3 biological replicates (with 3 tissues pooled per sample) from both male and female zebrafish. The target-enriched RNA-sequencing data were for male and female liver libraries for the validation of the approach, followed by the second round of enrichment for the male and female gill, intestine, and kidney libraries, along with the blank. 

Count data is in a new data directory (main/data/…). 

## Overall Objective:

Compare enriched transcriptomic data to unenriched data to validate that the approach does not significantly distort expression data.

## Initial Questions

Does target enrichment uniformly increase transcript counts for low, moderate, and highly expressed genes, or does it disproportionally enrich low or high count transcripts from the unenriched dataset? 

If there are disproportionate increases in enriched read counts, can they be explained by differences in the probe design, such as length, GC content, or melting temperature?

Is the proportional representation of transcripts distorted in the enriched data?

Is the relative expression of genes preserved in the enriched data set?

## Initial Plans

### Equivalence Tests

Compared to previous analysis of the enrichment approach which used Spearman's coefficients to determine the amount distortion in the datasets we will aim to use equivalence tests. Spearmen's coefficients are a correlation of the rank of data points rather which means that certain genes may be disportionally detected but their rank in expression may be conserved between enriched and unenriched libraries. For example if higher expressed genes are detected disproportionally less than lower expressed genes and the number of reads doesn't surpass those of higher expressed genes the ranking will be similar but expression ratios will be inaccurate. To determine if the representation of expression levels between the two datasets is equal (i.e., within tolerances we determine as acceptable), we will first normalize the counts to counts per million (CPM) and perform a linear regression. We will then use TOST (two one-sided t-tests) on the coefficients, using the TOSTER package. This will inform us if globally, across the ~600 bait-captured transcripts, the proportion of counts by gene are practically equivalent within tolerances. On the individual transcript level, we will similarly apply TOST to determine if the representation of individual transcripts in the unenriched and enriched datasets is equivalent, to identify any counts that may be biased. Depending on the results, we can perform an enrichment analysis using probe characteristics to determine if any particular probe characteristics may be causing distortion of the counts. 

Before we can perform equivalence tests, we need to determine what to set our tolerances to. To inform our choice of equivalence bounds, we will need to determine how much transcriptomics data naturally varies. To do this, we will first try empirically by comparing experimental or biological replicates of RNA-seq data to itself to determine how much noise is typical. Alternatively, we will review literature and consult with others to determine a defensible threshold to use. 
