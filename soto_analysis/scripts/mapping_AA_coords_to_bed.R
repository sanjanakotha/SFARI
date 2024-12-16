library(ensembldb)
library(EnsDb.Hsapiens.v86)
edbx <- EnsDb.Hsapiens.v86

uniprotID_to_genomic_coord_bed <- function(uniprotID, start, end, folder, base = "~/Desktop/Staller_Lab/SFARI/soto_analysis/outputs/mutations/") { 
  uni_rng <- IRanges(start = start, end = end,
                     names = uniprotID)
  uni_gnm <- proteinToGenome(uni_rng, edbx)

  ## Choosing one transcript if there are > 1 transcript for the uniprotID
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
  
  dir.create(file.path(base, folder), recursive = TRUE)
  savefilepath <- paste(base,  
                        folder, "/", uniprotID, ".bed", sep = "")
  write.table(df, savefilepath, quote=F, sep="\t", row.names=F, col.names=F, append = TRUE)
}

# Looping through all uniprotIDs and Ends in the TF AA coords
SFARI_TF_AA_coords <- read.csv("~/Desktop/Staller_Lab/SFARI/data/SFARI_TFs_with_knownADs_coords_ENSP.csv")
for (i in 1:nrow(SFARI_TF_AA_coords)){
  input_uniprotID <- SFARI_TF_AA_coords$uniprotID[i]
  input_start <- SFARI_TF_AA_coords$Start[i]
  input_end <- SFARI_TF_AA_coords$End[i]
  print(input_uniprotID)
  print(input_start)
  print(input_end)
  uniprotID_to_genomic_coord_bed(input_uniprotID, input_start, input_end, "cds_bed_format")
}

# Looping through all uniprotIDs and Ends in the AD AA coords
SFARI_AD_AA_coords <- read.csv("~/Desktop/Staller_Lab/SFARI/data/SFARI_ADs_AA_coords_redone_ENSP.csv")
for (i in 1:nrow(SFARI_AD_AA_coords)){
  input_uniprotID <- SFARI_AD_AA_coords$uniprotID[i]
  input_start <- SFARI_AD_AA_coords$Start[i]
  input_end <- SFARI_AD_AA_coords$End[i]
  print(input_uniprotID)
  print(input_start)
  print(input_end)
  uniprotID_to_genomic_coord_bed(input_uniprotID, input_start, input_end, "domains_bed_format/AD")
}

# Looping through all uniprotIDs and Ends in the DBD coords of SFARI TFs with known ADs
DBD_coords <- read.csv("~/Desktop/Staller_Lab/SFARI/data/SFARI_TFs_with_known_ADs_DBD_coords_ENSP.csv")
for (i in 1:nrow(DBD_coords)){
  input_uniprotID <- DBD_coords$uniprotID[i]
  input_start <- DBD_coords$Start[i]
  input_end <- DBD_coords$End[i]
  print(input_uniprotID)
  print(input_start)
  print(input_end)
  uniprotID_to_genomic_coord_bed(input_uniprotID, input_start, input_end, "domains_bed_format/DBD")
}

# # All known ADs
# AD_coords <- read.csv("~/Desktop/Staller_Lab/SFARI/data/known_ADs_on_TFs_coords.csv")
# for (i in 1:nrow(AD_coords)){
#   input_uniprotID <- AD_coords$uniprotID[i]
#   input_start <- AD_coords$Start[i]
#   input_end <- AD_coords$End[i]
#   print(input_uniprotID)
#   print(input_start)
#   print(input_end)
#   uniprotID_to_genomic_coord_bed(input_uniprotID, input_start, input_end, "ADs_bed_format", base = "~/Desktop/Staller_Lab/SFARI/data/")
# }



