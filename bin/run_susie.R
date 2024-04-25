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

  # get path of parquet from arguments
  sumstats_path = argv[1]
  # get path of LD from arguments
  ld_path = argv[2]
  # get path of variant reference from arguments
  varRef_path = argv[3] 

  varRef <- fread(varRef_path,data.table=F)

  chr = as.numeric(argv[4])
  start = as.numeric(argv[5])
  end = as.numeric(argv[6])
  varRef_rel = varRef %>% filter(CHR == chr & bp >= start & bp <= end)
  print(as_tibble(varRef))

  # load in parquet with Robert's strategy
  eQTL_ds <- open_dataset(sumstats_path)

  # load in LD matrix as data.frame
  ldInfo <- fread(ld_path,data.table=F)
  ldInfo = ldInfo[,-1]
  rownames(ldInfo) = colnames(ldInfo)
  relLdInfo <- ldInfo 

  # filter on relevant sub loci
  relLdInfo = relLdInfo[which(rownames(relLdInfo) %in% varRef_rel$ID),which(colnames(relLdInfo) %in% varRef_rel$ID)]
  
  # filter and order on same information in sumstats & LD matrix
  eQTL_rel <- eQTL_ds %>% filter(variant %in% rownames(relLdInfo)) %>% collect()
  relLdInfo = relLdInfo[which(rownames(relLdInfo) %in% eQTL_rel$variant),which(colnames(relLdInfo) %in% eQTL_rel$variant)]
  relLdInfo = relLdInfo[match(eQTL_rel$variant,rownames(relLdInfo)),match(eQTL_rel$variant,colnames(relLdInfo))]

  print(as_tibble(eQTL_rel))
  print(as_tibble(relLdInfo))

  if(all(eQTL_rel$variant == rownames(relLdInfo))){
    L = 5
    fitted_rss2 = susie_rss(bhat = eQTL_rel$beta, shat = eQTL_rel$standard_error, R = as.matrix(relLdInfo), n = eQTL_rel$sample_size[1], L = L, estimate_residual_variance = FALSE, verbose=T)

    print("Finished!")
    print(fitted_rss2$converged)

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
        colnames(lbfOut) = paste("lbf_cs",1:L,sep="_")
        eQTL_rel = cbind(eQTL_rel,lbfOut)
        outMat = rbind(outMat,eQTL_rel)
      }
    } else {
      print(paste("problem with",feat," not converged."))
    }
    
  }
  if(formatted & print){
    write.table(outMat,argv[7],quote=F,sep="\t",row.names=F)
  }
}
if(direct & print){
  sink()
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
