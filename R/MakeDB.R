MakeDB <- function(dbtype = "nucl", fsa_file){
  syscommand <- paste("makeblastdb -in ", fsa_file, " -dbtype ", dbtype, " -parse_seqids", sep = "")
  system(syscommand)
}
