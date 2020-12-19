#' simple function to emit a trimming script for STORM-seq preps
#' 
#' @param stub  the shared part of the fastq filenames 
#' 
#' @export 
make_storm_trimmer <- function(stub) { 

  message("# launch trim_galore to remove illumina adapters")
  message("trim_galore --illumina \\")
  message("  --trim-n \\")
  message("  --length 36 \\")
  message("  --paired \\")
  message("  --clip_R2 3 \\")
  message("  --fastqc \\")
  message("  ", stub, "_L000_R1_001.fastq.gz \\")
  message("  ", stub, "_L000_R2_001.fastq.gz")

}
