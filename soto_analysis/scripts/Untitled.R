#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
SFARI_AD_AA_coords <- read.csv(args[1]

library(ensembldb)
library(EnsDb.Hsapiens.v86)
edbx <- EnsDb.Hsapiens.v86

uniprotID_to_genomic_coord_bed <- function(uniprotID, start, end, folder, idType) { 
  uni_rng <- IRanges(start = start, end = end,
                     names = uniprotID)
  uni_gnm <- proteinToGenome(uni_rng, edbx, idType = idType)
  
  ## Choosing one transcript if there are > 1 transcript for the uniprotID
  # if(length(uni_gnm[[1]]) > 1) {
  #if (typeof(uni_gnm[[1]]) == "list") {
  if (length(GRangesList(uni_gnm[[1]])) > 1) {
    uni_gnm[[1]] = unlist(GRangesList(uni_gnm[[1]][1]))
  }
  
  gr <- unlist(GRangesList(uni_gnm))
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=elementMetadata(gr)[,1],
                   ensts=elementMetadata(gr)[,2],
                   strands=strand(gr))
  
  savefilepath <- paste("~/Desktop/Staller_Lab/SFARI/soto_analysis/outputs/mutations/",  
                        folder, "/", uniprotID, ".bed", sep = "")
  write.table(df, savefilepath, quote=F, sep="\t", row.names=F, col.names=F, append = TRUE)
}

# Looping through all uniprotIDs and Ends in the AD AA coords
# SFARI_AD_AA_coords <- read.csv("~/Desktop/Staller_Lab/SFARI/data/SFARI_ADs_AA_coords.csv")
for (i in 1:nrow(SFARI_AD_AA_coords)){
  input_uniprotID <- SFARI_AD_AA_coords$uniprotID[i]
  input_start <- SFARI_AD_AA_coords$Start[i]
  input_end <- SFARI_AD_AA_coords$End[i]
  print(input_uniprotID)
  print(input_start)
  print(input_end)
  uniprotID_to_genomic_coord_bed(input_uniprotID, input_start, input_end, "domains_bed_format", "uniprot_id")
}

