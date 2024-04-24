#!/usr/bin/env Rscript

## ----
## Author:  T. van Lieshout
## Email:   t.van.lieshout@umcg.nl
##
## Copyright (c) T. van Lieshout, 2024
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## A copy of the GNU General Public License can be found in the LICENSE file in the
## root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
## ----

# Load libraries
library(data.table)
library(tidyverse)
library(arrow)
library(susieR)

# Declare constants

# Declare function definitions

# Main
print = T
direct = T
formatted = T

if(formatted){
  outMat = NULL
}

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv=NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  # get path of parquet from arguments, but TMP for testing
  #sumstats_path = "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/fine_mapping/ExtractMetaAnalysisResults/work/b5/a52ebf395ab7b62c9b9904b29dffe2/tmp_empirical/phenotype=ENSG00000154719/af25b5073d3c4a92a7356b0ff37e276c.parquet"
  sumstats_path = argv[1]
  # get path of LD from arguments, but TMP for testing
  #ld_path = "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/fine_mapping/ExtractMetaAnalysisResults/work/b5/a52ebf395ab7b62c9b9904b29dffe2/ld.21_27778421_31314679.csv.gz"
  ld_path = argv[2]
  # load in parquet with Robert's strategy
  eQTL <- open_dataset(sumstats_path)
  eQTL_rel <- as_tibble(eQTL)

  # load in LD matrix as data.table
  ldInfo <- fread(ld_path,data.table=F)
  ldInfo = ldInfo[,-1]
  rownames(ldInfo) = colnames(ldInfo)
  relLdInfo <- ldInfo 
  
  
  eQTL_rel = eQTL_rel[which(eQTL_rel$variant %in% rownames(relLdInfo)),]
  relLdInfo = relLdInfo[which(rownames(relLdInfo) %in% eQTL_rel$variant),which(colnames(relLdInfo) %in% eQTL_rel$variant)]
  relLdInfo = relLdInfo[match(eQTL_rel$variant,rownames(relLdInfo)),match(eQTL_rel$variant,colnames(relLdInfo))]

  ##For debug.
#  png("/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/fine_mapping/ExtractMetaAnalysisResults/work/00/cc7f0912b1e3ea04ba834618529d98/TMP_debug.png")
#  plot(as.numeric(eQTL_rel$snp_position),-log10(as.numeric(eQTL_rel$p_value)))
#  dev.off()

  if(all(eQTL_rel$variant == rownames(relLdInfo))){
#    fitted_rss2 = tryCatch({
#      fitted_rss2 = susie_rss(bhat = eQTL_rel$beta, shat = eQTL_rel$standard_error, R = as.matrix(relLdInfo), n = eQTL_rel$sample_size[1], L = 5, estimate_residual_variance = TRUE)
#    },error = function(e){
#      fitted_rss2 = susie_rss(bhat = eQTL_rel$beta, shat = eQTL_rel$standard_error, R = as.matrix(relLdInfo), n = eQTL_rel$sample_size[1], L = 5, estimate_residual_variance = FALSE)
#      estimated_res_var = F
#      return(fitted_rss2)
#    })
    fitted_rss2 = susie_rss(bhat = eQTL_rel$beta, shat = eQTL_rel$standard_error, R = as.matrix(relLdInfo), n = eQTL_rel$sample_size[1], L = 10, estimate_residual_variance = FALSE, verbose=T)

    print("Finished!")
    print(fitted_rss2$converged)
    #susie_plot(fitted_rss2, y= "PIP", eQTL_rel$beta)

    if(fitted_rss2$converged){
      if(direct){
        print(summary(fitted_rss2))
      }

      if(formatted){
        eQTL_rel["SusieRss_pip"] = fitted_rss2$pip
        eQTL_rel["SusieRss_CS"] = NA
        if(length(fitted_rss2$sets[[1]])!=0){
          for(l in 1:length(fitted_rss2$sets[[1]])){
            eQTL_rel$SusieRss_CS[fitted_rss2$sets[[1]][[l]]]=l
          }
        }
        lbfOut = t(fitted_rss2$lbf_variable)
        colnames(lbfOut) = paste("lbf_cs",1:10,sep="_")
        eQTL_rel = cbind(eQTL_rel,lbfOut)
        outMat = rbind(outMat,eQTL_rel)
      }
    } else {
      print(paste("problem with",feat," not converged."))
    }
    
  }
  if(formatted & print){
    #write.table(outMat,"/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/fine_mapping/ExtractMetaAnalysisResults/work/b5/a52ebf395ab7b62c9b9904b29dffe2/finemap_results.txt",quote=F,sep="\t",row.names=F)
    write.table(outMat,argv[3],quote=F,sep="\t",row.names=F)
  }
}
if(direct & print){
  sink()
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
