###############################################
### functions!
###############################################


###############################################
# fun
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
# import peaks
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
# read_peaks
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
# intersect_peaks
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
# create consensus peaks
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
# get promoter regions
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






###############################################
# profile tss
###############################################
#' do the whole tss thing
#' 
#' @description 
#' 
#' @param gencode_gr
#'  
#'  
#' @param biotype
#' 
#'
#' @param upstream
#'
#'
#' @param downstream
#'

profile_tss <- function(peaks, 
                        promoters_gr,
                        upstream = 3e3,
                        downstream = 3e3) {
  
  # performing coverage function 
  peak_coverage <- coverage(peaks)
  # keeping track of overlaps in RLE
  coverage_length <- elementNROWS(peak_coverage)
  # Defining a GRanges of the promter window
  coverage_gr <- GRanges(seqnames = names(coverage_length),
                         IRanges(start = rep(1, length(coverage_length)), 
                                 end = coverage_length))
  
  # defining the promoters 
  promoters_gr <- subsetByOverlaps(promoters_gr, 
                                   coverage_gr, 
                                   type="within", 
                                   ignore.strand=TRUE)
  # making sure the chromosomes represented are used (prevent error if chr is missing)
  chromosomes <- intersect(names(peak_coverage), 
                           unique(as.character(seqnames(promoters_gr))))
  # arranging chromosomes in the same order
  peak_coverage <- peak_coverage[chromosomes]
  # converting to InterRangesList
  promoters_ir <- as(promoters_gr, "IntegerRangesList")[chromosomes]
  # creating a views object for promoter coverage (because in RLE)
  promoter_peak_view <- Views(peak_coverage, promoters_ir)
  # turning into a vector with ViewApply (because in RLE keeping track of where overlaps are)
  promoter_peak_view <- lapply(promoter_peak_view, function(x) t(viewApply(x, as.vector)))
  # binding each of the view vectors
  promoter_peak_matrix <- do.call("rbind", promoter_peak_view)
  # grabing and reversing promoters on the - strand
  minus_idx <- which(as.character(strand(promoters_gr)) == "-")
  
  # reversing the order from 6,000 - 1 to 1- 6000
  promoter_peak_matrix[minus_idx,] <- promoter_peak_matrix[minus_idx,
                                                           ncol(promoter_peak_matrix):1]
  # eliminating promoters with no binding 
  promoter_peak_matrix <- promoter_peak_matrix[rowSums(promoter_peak_matrix) > 1,]
  # summing all the vectors of a given DBP to the promoter window
  peak_sums <- colSums(promoter_peak_matrix)
  # calculating the density at each position in the promoter window
  peak_dens <- peak_sums/sum(peak_sums)
  # making it go from -3K to + 3K and creating a df
  metaplot_df <- data.frame(x = -upstream:(downstream-1),
                            dens = peak_dens)
  
  return(metaplot_df)
}
















