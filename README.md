## Extract meta-analysis results

This is a combination of several Nextflow pipeline for extracting subsets of data from eQTLGen phase II genome-wide eQTL meta-analyses results. Because of the sheer size of the genome-wide association summary statistics, this task utilizes the Nextflow parallel processing.

The typical usage scenarious are:
- Extracting the subset of effects which reach to certain significance threshold (e.g. P<5e-8 or the one determined by stricter multiple testing threshold). This more manageable subset of results can be further summarised, visualised and analysed. (ExtractResults.nf)
- Extracting the global association summary statistics for certain subset of genes of interest. These global association summary statistics can be visualised on Manhattan plot or used in downstream analyses. (ExtractResults.nf)
- Defining a set of uncorrelated variants (ExtractUncorrelatedVariants.nf)
- Calculating gene-gene correlations (ExtractGeneCorrelations.nf)
- Preparing in-sample LD using permuted summary statistics (PrepareLdPanel.nf)
- Fine-mapping (FineMapping.nf)

### Usage instructions

#### Requirements for the system

- Have access to HPC with multiple cores.
- Have Bash >=3.2 installed.
- Have Java >=8 installed.
- Have Slurm scheduler managing the jobs in the HPC.
- HPC has Singularity installed and running.

#### Setup of the pipeline
You can either clone it by using git (if available in HPC):

`git clone https://github.com/eQTLGen/ExtractMetaAnalysisResults.git`

Or just download this from the gitlab/github download link and unzip.

#### Input files

- Folder with the output `.parquet` files from [MetaAnalysis pipeline](https://github.com/eQTLGen/MetaAnalysis).