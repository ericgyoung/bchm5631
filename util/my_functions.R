###############################################
### functions!
###############################################


###############################################
#' A function to multiply two numbers
#'
#' @description 
#' This function will multiply the input values of X and Y
#' 
#' @param x one number you'd like to multiply
#' y the other number you'd like to multiply
fun <- function(x, y) {
  ans <- x * y
  return(ans)
}
###############################################
#' A function to find import peaks from list of dbp's
#'
#' @description 
#' establishing the parameter consensus_file_path is needed
#' extracting TF name as a variable along the way (check env)
#' 
#' @param consensus_file_path one number you'd like to multiply
#' y the other number you'd like to multiply
import_peaks <- function(consensus_file_path = "/scratch/Shares/rinnclass/CLASS_2022/data/peaks") {
  
  # Setting some variables needed in main part of function (same as above -- peak_files & tf_name)
  peak_files <- list.files(consensus_file_path, full.names = T)
  
  # Make an object with each TF name for indexing and merging later
  tf_name <- sapply(peak_files, function(x){
    y <-  str_extract(x, "([^\\/]+$)")
    unlist(strsplit(y, "_"))[[1]]
  })
  # Here is the heart of the function that will import each file as GRanges (we can use for overlaps)
  # takes 
  
  # peak_list <- c()
  # for(i in 1:length(peak_files)) {
  #   # Import peak files
  #   peaks <- rtracklayer::import(peak_files[i])
  #   # Append this GRanges object to the of the list.
  #   peak_list <- c(peak_list, peaks)
  #   # Name the list elements by their TF name (we just made above)
  #   names(peak_list)[length(peak_list)] <- tf_name[i]
  # }
  peak_list = lapply(peak_files,rtracklayer::import)
  names(peak_list) = tf_name
  return(peak_list)
}
###############################################
#' @description 
#' take in broad peaks and outputs canonical chromosome GRanges object
#' 
#' @param broad_peak_file 
read_peaks <- function(broad_peak_file, filter_to_canonical_chr = TRUE) {
  dat <- read.table(broad_peak_file, sep = "\t")
  if(filter_to_canonical_chr == TRUE) {
    dat <- dat[dat$V1 %in% c(paste0("chr", 1:22), "chrM", "chrX", "chrY"),]
  }
  gr <- GRanges(seqnames = dat$V1,
                ranges = IRanges(start=dat$V2,end=dat$V3))
  return(gr)
}
###############################################
#' intersect replicates into a "consensus peak list" 
#' 
#' @description 
#' this function will take the  union of peak widths across replicates for a given
#' DNA binding protein. the function that will take a list of granges objects and return 
#  one granges object with merged peaks that are in all replicates
#' 
#' @param 
#'  the path to consensus peak files
#' # We're going to iterate over all the files to make it work. 
intersect_peaks <- function(peak_list) {

    combined_peaks <- peak_list[[1]]
    for(i in 2:length(peak_list)) {
      suppressWarnings(pl_ov <- findOverlaps(combined_peaks, peak_list[[i]]))
      pl1 <- combined_peaks[unique(pl_ov@from)]
      pl2 <- peak_list[[i]][unique(pl_ov@to)]
      suppressWarnings(combined_peaks <- GenomicRanges::reduce(union(pl1, pl2)))

    }
    return(combined_peaks)
  }



###############################################
#' intersect replicates into a "consensus peak list" 
#' 
#' @description 
#' this function will take the  union of peak widths across replicates for a given
#' DNA binding protein. the function that will take a list of granges objects and return 
#  one granges object with merged peaks that are in all replicates
#' 
#' @param 
#'  the path to consensus peak files
#' # We're going to iterate over all the files to make it work. 
create_consensus_peaks <- function(broadpeakfilepath = "/scratch/Shares/rinnclass/CLASS_2022/data/peaks/") {
  
  # For now we can set broadpeakfilepath
  
  # broadpeakfilepath <- "/Shares/rinn_class/data/CLASS_2022/class_exeRcises/analysis/11_consensus_peak_exercise"
  # making a list of file paths to the (similar to import_peaks function)
  fl <- list.files(broadpeakfilepath, 
                   full.names=TRUE)
  fl <- fl[grep("peaks.broadPeak", fl)]
  # getting a DBP name for same index as each file path
  tf_name <- sapply(fl, function(x){
    y <-  str_extract(x, "([^\\/]+$)")
    unlist(strsplit(y, "_"))[[1]]
  })
  
  
  # making sure there is a replicate and creating "unique_tf" index
  # This will be used in a forloop
  tf_df <- data.frame(table(tf_name)) %>%  # data.frame(table(tf_name))
    filter(Freq > 1)
  unique_tf <- as.character(tf_df$tf_name) # unique_tf
  
  # Now a nested for loop (2 for loops) to make GRanges of peak files.
  # This is similar to read_peaks
  consensus_peaks <- list()
  for(i in 1:length(unique_tf)) {
    
    # load all the peak files corresponding to this DBP[i] in unique_tf.
    # tf <- unique_tf[1] -- allows us to look at output
    tf <- unique_tf[i]
    print(tf)
    # indexing unique DBP name to file path (e.g., first 8 are CTCF files)
    tf_index <- grep(tf, tf_name)
    # takes the TF name and grabs the index in fl for those replicates
    tf_files <- fl[tf_index]
    
    # now make a list of GRanges in a peak_list using another for loop
    # READ_PEAKS being used 
    peak_list <- c()
    for(j in 1:length(tf_files)) {
      # See the read peaks function to know what subfunctions are called.
      peak_list <- c(peak_list, read_peaks(tf_files[j]))
      # same read peaks function and we now have each DBP indexed in tf_files
    }
    # READ_PEAKS now being used
    # filtering chromosomes -- redundant since read peaks does this too -- oh well.
    canonical_chr <- c(paste0("chr", 1:22), "chrM", "chrX", "chrY")
    for(i in 1:length(peak_list)) {
      peak_list[[i]] <-peak_list[[i]][which(seqnames(peak_list[[i]]) %in% canonical_chr)]
    }
    # Now we use intersect_peaks functino to find overlaps 
    # INTERSECT_PEAKS now being used
    final_peakset <- intersect_peaks(peak_list = peak_list)
    if(length(final_peakset) > 0) {
      final_peakset$name <- paste0(tf, "_", 1:length(final_peakset))
    }
    
    consensus_peaks <- c(consensus_peaks, list(final_peakset))
    names(consensus_peaks)[length(consensus_peaks)] <- tf
  }
  return(consensus_peaks)
}
###############################################
#' function to subset features for promomters. 
#' 
#' @description feature_subset
#' Take a gencode gtf to subset the biotype of promoters we want as a set of GRanges
#' 
#' @param gencode_gr
#'  set of genomic features as a GRanges object
#'  
#' @param biotype
#' this takes "lncRNA" or "protein-coding" as input for promoter type
#'
#' @param upstream
#'To add upstream sequence to feature subset
#'
#' @param downstream
#'To add downstream sequence to feature subset

get_promoter_regions <- function(gencode_gr, biotype, upstream = 3e3, downstream = 3e3) {
  
  genes <- gencode_gr[gencode_gr$type == "gene"]
  genes <- genes[genes$gene_type %in% biotype]
  
  proms <- GenomicRanges::promoters(genes, upstream = upstream, downstream = downstream)
  
  return(proms)
  
}
























