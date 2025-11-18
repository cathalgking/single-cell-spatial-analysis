#!/usr/bin/env Rscript

# convert_seurat_greyscale.R
# Usage:
#   Rscript convert_seurat_greyscale.R file1.rds file2.rds ...
# For each input file, a new file "<original>_gray.rds" is written.

suppressPackageStartupMessages({
  library(Seurat)
  library(tools)  # for file_path_sans_ext
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No input files provided.\nUsage: Rscript convert_seurat_greyscale.R file1.rds [file2.rds ...]", call. = FALSE)
}

############################################################
## TRUE GREYSCALE CONVERSION FOR ALL IMAGES IN A SEURAT OBJECT
############################################################

## Function to convert an RGB (or RGBA) numeric array to greyscale
greyscale_image_array <- function(img_array) {
  if (length(dim(img_array)) != 3) {
    stop("Expected a 3D array: height x width x channels.")
  }
  
  h  <- dim(img_array)[1]
  w  <- dim(img_array)[2]
  ch <- dim(img_array)[3]
  
  if (ch < 3) stop("Image must have at least 3 channels (RGB).")
  
  # Standard luminance conversion
  gray <- 0.299 * img_array[, , 1] +
    0.587 * img_array[, , 2] +
    0.114 * img_array[, , 3]
  
  # Copy grayscale to R/G/B channels
  img_gray <- img_array
  img_gray[, , 1] <- gray
  img_gray[, , 2] <- gray
  img_gray[, , 3] <- gray
  
  # If a 4th channel exists (alpha), leave it unchanged
  return(img_gray)
}

## Function to apply to all images inside a Seurat object
convert_seurat_images_to_grayscale <- function(seu) {
  if (length(seu@images) == 0) {
    message("  No images found in this Seurat object.")
    return(seu)
  }
  
  img_names <- names(seu@images)
  
  for (img in img_names) {
    message("  Converting to greyscale: ", img)
    
    arr <- seu@images[[img]]@image
    
    if (!is.array(arr)) {
      stop(paste("Image", img, "is not stored as an array."))
    }
    
    arr_gray <- greyscale_image_array(arr)
    
    seu@images[[img]]@image <- arr_gray
  }
  
  return(seu)
}

############################################################
## PROCESS ALL FILES PASSED ON THE COMMAND LINE
############################################################

for (f in args) {
  message("Processing: ", f)
  
  if (!file.exists(f)) {
    warning("  File not found, skipping: ", f)
    next
  }
  
  seu <- readRDS(f)
  
  seu_gray <- convert_seurat_images_to_grayscale(seu)
  
  out_file <- paste0(file_path_sans_ext(f), "_gray.rds")
  saveRDS(seu_gray, out_file)
  
  message("  Saved greyscale Seurat object to: ", out_file, "\n")
}
