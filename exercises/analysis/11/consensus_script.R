library(tidyverse)
library(GenomicRanges)

# source
source("/scratch/Shares/rinnclass/CLASS_2022/EricY/bchm/util/my_functions.R")

consensus_peaks = create_consensus_peaks("/scratch/Shares/rinnclass/CLASS_2022/data/peaks")

# export
for (i in 1:length(consensus_peaks)) {
  rtracklayer::export(consensus_peaks[[i]],
                      paste0("/scratch/Shares/rinnclass/CLASS_2022/EricY/bchm/exercises/analysis/11/consensus_peaks",
                             names(consensus_peaks)[i],"_consensus_peaks.bed"))
}