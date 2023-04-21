## Extract meta-analysis results

This is a small Nextflow pipeline for extracting the subsets of data from eQTLGen phase II genome-wide eQTL meta-analyses results. Because of the sheer size of the genome-wide association summary statistics, this task utilizes the Nextflow parallel processing.

Currently the typical usage scenarious are:

- Extracting the subset of effects which reach to certain significance threshold (e.g. P<5e-8 or the one determined by stricter multiple testing threshold). This more manageable subset of results can be further summarised, visualised and analysed.
- Extracting the global association summary statistics for certain subset of genes of interest. These global association summary statistics can be visualised on Manhattan plot or used in downstream analyses.

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

#### Usage instructions

You need to put 1000G-30x reference files to the folder `ExtractMetaAnalysisResults/bin/hase/data`. These files are provided in Dropbox, together with imputation reference (~30GB).

- Download the reference folder: `wget https://www.dropbox.com/s/6g58ygjg9d2fvbi/eQTLGenReferenceFiles.tar.gz?dl=1`
- Download md5sum: `wget https://www.dropbox.com/s/ekfciajzevn6o1l/eQTLGenReferenceFiles.tar.gz.md5?dl=1`
- Check the md5sum: `md5sum --check eQTLGenReferenceFiles.tar.gz.md5`
- Unzip: `tar -xfvz eQTLGenReferenceFiles.tar.gz`
- Copy the required files to correct place: `cp hg38/hase_reference/* bin/hase/data/.`

##### Arguments of the pipeline

###### Mandatory arguments:

`--inputfolder `    Path to the folder with HASE result `.parquet` files.

`--outputfile  `    Path to output file where results are written as tab-delimited text file.

`--NumberOfCovariates`  Number of covariates in the model. This is used for calculating the degrees of freedom and subsequently P-value. Currently it assumes that the same set of covariates (e.g. 24) was used in each cohort. More precise settings are **TBA**.

`--SnpRef`  Reference file for variant mapper `1000G-30x.ref.gz`. Has to be gzipped and space-delimited. This is used to include information about alleles.

###### Optional arguments:

`--PThresh` P-value threshold for filtering the results. Defaults `None`.

`--PhenFilter`  Filter to include phenotypes into the association results. One of two: first, file with the list of phenotypes to include. In that case, no header expected. Second, individual phenotype IDs specified and separated by space. Defaults to `None`.

`--SnpFilter`   Filter to include SNPs/variants into the association results.  One of two: first, file with the list of SNPs/variants to include. In that case, no header expected. Second, individual SNP IDs specified and separated by space. Defaults to `None`.

`--Chunks`    Number of chunks which was used for original meta-analysis. Defaults to 100.

#### Running the meta-analysis command

Go to folder `ExtractMetaAnalysisResults` and modify the Slurm script template `submit_ExtractMetaAnalysisResults_template.sh` with your input paths. Below is a example template for Slurm scheduler. Here we assume that meta-analysis was ran on 100 chunks and we want to get all effects reaching P<5e-8.

```bash
#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --job-name="ExtractHaseResults"

# Here load needed system tools (Java 1.8 is strictly required, one of two: singularity or conda for python 2.7 are needed,
# depending on the method which is used for dependency management)
module load any/jdk/1.8.0_265
module load any/singularity/3.5.3
module load squashfs/4.4

# NB! Following will extract all the genes and SNPs! In case of full mapping, 
# it is advisable to use gene filter or P-value filter to get manageable subset!

nextflow_path=[path to folder where Nextflow executable is]

input_folder=[folder with input .parquet files]
output_file=[path to tab-separated output file]
nr_covariates=[Number of covariates in the model]
snp_ref=bin/hase/data/1000G-30x.ref.gz

#Optional arguments:
#--PThresh               P-value threshold for filtering the results. Defaults None.
#--PhenFilter            Filter to include phenotypes into the association results. One of two: first, file with the list of phenotypes to include. In that case, no header expected. Second, individual phenotype IDs specified and separated by space. Defaults to None.
#--SnpFilter             Filter to include SNPs/variants into the association results.  One of two: first, file with the list of SNPs/variants to include. In that case, no header expected. Second, individual SNP IDs specified and separated by space. Defaults to None.
#--Chunks                Number of chunks which was used for original meta-analysis. Defaults to 100.

NXF_VER=20.10.0 ${nextflow_path}/nextflow run analysis.nf \
--inputfolder ${input_folder} \
--outputfile ${output_file} \
--NumberOfCovariates ${nr_covariates} \
--SnpRef ${snp_ref} \
--PThresh 5e-8 \
-resume \
-profile slurm,singularity

# Compress the output
gzip ${output_file}
```

You can save the modified script version with different name, e.g. `submit_ExtractMetaAnalysisResults_p5e-8.sh`.

Then submit the job `sbatch submit_ExtractMetaAnalysisResults_p5e-8.sh`. This initiates pipeline, makes analysis environment (using singularity) and automatically submits the steps in correct order and parallel way. Separate work directory is made to the folder and contains all interim files.

##### Monitoring and debugging

- Monitoring:
  - Monitor the `slurm-***.out` log file and check if all the steps finish without error. Trick: command `watch tail -n 20 slurm-***.out` helps you to interactively monitor the status of the jobs.
  - Use `squeue -u [YourUserName]` to see if individual tasks are in the queue.
- If the pipeline crashes (e.g. due to walltime), you can just resubmit the same script after the fixes. Nextflow does not rerun completed steps and continues only from the steps which had not completed.
- When the work has finished, download and check the job report. This file is automatically written to your output folder pipeline_info subfolder, for potential errors or warnings. E.g. `output/pipeline_info/MetaAnalysis_report.html.`
- When you need to do some debugging, then you can use the last section of aforementioned report to figure out in which subfolder from work folder the actual step was run. You can then navigate to this folder and investigate the following hidden files:
  - `.command.sh`: script which was submitted
  - `.command.log`: log file for seeing the analysis outputs/errors.
  - `.command.err`: file which lists the errors, if any.

**NB!** Pipeline does not overwrite the files in the output folder. Therefore, **before re-running the pipeline, delete the output files in the output folder!**

**NB!** Be warned, if you run the current pipeline without P-value filter/gene filter/SNP filter, then the pipeline writes the full eQTL summary statistics into one enormous unzipped tab-delimited `.txt` file, many TBs large. **You probably want to use some filter here.** 

##### Output

Tab-delimited `.txt` file with the following information:

- phenotype: gene
- SNP: variant
- ref_all: reference allele
- alt_all: alternative allele
- beta: beta, corresponding to alt_all
- se: standard error of the beta
- P: P-value, calculated by the pipeline
- i_squared: heterogeneity I<sup>2</sup> for given effect, useful for QC.
- N: sample size, in how many samples this effect was tested, useful for QC.

### TODO

- Implement functionality of writing out individual global association profile for each gene. Files should be compressed as good as possible, to save storage.