library(Rcpp)
#sourceCpp("src/FastaSlice_ver04.cpp")

BlastN <- function(query, db, slicesize=100, leaps=900){
  if (file.exists(query) && file.exists(db)){
  x <- FastaSlice(query, db, slicesize, leaps)
  }

  else{
    print("File does not exist. ")
  }
  #write(query, file = "query/query.fsa")
  #syscommand <- paste("blastn -db ", target, " -query query/query.fsa -outfmt '6 qseqid sseqid length pident score' -out rtest.out", sep = "")
  #system(syscommand)
}
