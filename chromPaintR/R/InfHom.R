InfHom <- function(chromdat,
                   scaff.lim = F,
                   dbchrom.lim = F,
                   max.scaffs = 10,
                   max.dbchroms = 10,
                   resolution =1000){
  require(dplyr)
  ordered.scaffs.indices <- order(Reduce(c,chromdat$summary$`Scaffold Lengths`), decreasing = T)
  query.hits.list <- chromdat$summary$`Database hits`[ordered.scaffs.indices]
  totalhits <- lapply(query.hits.list, length)
  alldbhits <- Reduce(c, chromdat$summary$`Database hits`)

  db.chroms <- unique(alldbhits)
  ordered.dbchroms.indices <- order(as.numeric(chromdat$lengthdata[[3]]), decreasing = T)
  ordered.dbchroms <- chromdat$lengthdata[[4]][ordered.dbchroms.indices]

  if(dbchrom.lim){
    ordered.dbchroms.pared <- sapply(db.chroms, x=ordered.dbchroms, grep)[1:max.dbchroms]
  }else{
    ordered.dbchroms.pared <- sapply(db.chroms, x=ordered.dbchroms, grep)
  }

  db.hits.list <- list()
  db.chroms <- names(sort(ordered.dbchroms.pared))

  for(chrom in db.chroms){
    db.hits.list <- c(db.hits.list, sum(as.numeric(alldbhits==chrom)))
  }

  if(scaff.lim & (length(query.hits.list)>max.scaffs)){
    actual.genome <- as.data.frame(matrix(0, max.scaffs, length(db.chroms)+1))
    colnames(actual.genome) <- c("scaff.names",db.chroms)
    actual.genome$scaff.names <- names(query.hits.list)[1:max.scaffs]
  }else{
    actual.genome <- as.data.frame(matrix(0,length(query.hits.list), length(db.chroms)+1))
    colnames(actual.genome) <- c("scaff.names",db.chroms)
    actual.genome$scaff.names <- names(query.hits.list)
  }

  for(chrom in db.chroms){
    actual.query.col <- c()
    if(scaff.lim & (length(query.hits.list)>max.scaffs)){
      for(i in 1:max.scaffs)
      {
        actual.query.col <- c(actual.query.col, sum(as.numeric(query.hits.list[[i]] == chrom)))
      }
      actual.genome[,chrom] <- actual.query.col
    }else{
      for(i in 1:length(query.hits.list))
      {
        actual.query.col <- c(actual.query.col, sum(as.numeric(query.hits.list[[i]] == chrom)))
      }
      actual.genome[,chrom] <- actual.query.col
    }
  }
  actual.genome.counts <- actual.genome[,-1]

  actual.genome.counts <- as.data.frame(actual.genome.counts)
  if(scaff.lim&(length(query.hits.list)>max.scaffs))
  {
    rownames(actual.genome.counts) <- names(query.hits.list)[1:max.scaffs]
  }else{
    rownames(actual.genome.counts) <- names(query.hits.list)
  }

  if(scaff.lim&(length(query.hits.list)>max.scaffs)){
    hyp.test.results <- matrix(0, max.scaffs, length(db.chroms))
    rownames(hyp.test.results) <- names(query.hits.list)[1:max.scaffs]
  }else{
    hyp.test.results <- matrix(0, length(query.hits.list), length(db.chroms))
    rownames(hyp.test.results) <- names(query.hits.list)
  }
  colnames(hyp.test.results) <- db.chroms
  for(scaff in actual.genome$scaff.names){
    scaff.size <- sum(actual.genome[which(actual.genome$scaff.names == scaff),-1])
    for(chrom in db.chroms){
      scaff.chromhits <- actual.genome[which(actual.genome$scaff.names == scaff), chrom]
      null.dist <- c()
      for(i in 1:resolution){
        null.dist[i] <- sum(sample(rep(db.chroms, times = Reduce(c,db.hits.list)),
                                   size = scaff.size,
                                   replace=F)==chrom)
      }
      scaff.chrom.pval <- sum(null.dist>=scaff.chromhits)/resolution
      hyp.test.results[scaff,chrom] <- scaff.chrom.pval
    }
  }

  actual.genome.counts.sig <- actual.genome.counts
  for(row in 1:nrow(hyp.test.results)){
    for (col in 1:ncol(hyp.test.results)){
      if(hyp.test.results[row,col]<.001){
        actual.genome.counts.sig[row,col] <- paste(signif(actual.genome.counts[row,col],4), "***", sep = "")
      }
      else if(hyp.test.results[row,col]<.01){
        actual.genome.counts.sig[row,col] <- paste(signif(actual.genome.counts[row,col],4), "**", sep = "")
      }
      else if(hyp.test.results[row, col]<.05){
        actual.genome.counts.sig[row,col] <- paste(signif(actual.genome.counts[row,col],4), "*", sep = "")
      }
      else{
        actual.genome.counts.sig[row,col] <- signif(actual.genome.counts[row,col],4)
      }
    }
  }
  return(actual.genome.counts.sig)
}






# library(dplyr)
#
# query.hits.list <- chrom.test$`Database hits`
# totalhits <- lapply(query.hits.list, length)
# alldbhits <- Reduce(c, chrom.test$`Database hits`)
# db.chroms <- unique(alldbhits)
# db.hits.list <- list()
#
# for(chrom in db.chroms){
#   db.hits.list <- c(db.hits.list, sum(as.numeric(alldbhits==chrom)))
# }
#
# actual.genome <- as.data.frame(matrix(0,length(query.hits.list), length(db.chroms)+1))
# colnames(actual.genome) <- c("scaff.names",db.chroms)
# actual.genome$scaff.names <- names(query.hits.list)
#
# for(chrom in db.chroms){
#   actual.query.col <- c()
#   for(i in 1:length(query.hits.list))
#   {
#     actual.query.col <- c(actual.query.col, sum(as.numeric(query.hits.list[[i]] == chrom)))
#   }
#   actual.genome[,chrom] <- actual.query.col
# }
#
# hyp.test.results <- matrix(0, length(query.hits.list), length(db.chroms))
# rownames(hyp.test.results) <- names(query.hits.list)
# colnames(hyp.test.results) <- db.chroms
# for(scaff in actual.genome$scaff.names){
#   scaff.size <- sum(actual.genome[which(actual.genome$scaff.names == scaff),-1])
#   for(chrom in db.chroms){
#     scaff.chromhits <- actual.genome[which(actual.genome$scaff.names == scaff), chrom]
#     null.dist <- c()
#     for(i in 1:10000){
#       null.dist[i] <- sum(sample(rep(db.chroms, times = Reduce(c,db.hits.list)),
#                                  size = scaff.size,
#                                  replace=F)==chrom)
#     }
#     scaff.chrom.pval <- sum(null.dist>=scaff.chromhits)/10000
#     hyp.test.results[scaff,chrom] <- scaff.chrom.pval
#   }
# }
#
# size.scaff1 <- sum(actual.genome[1,-1])
# size.scaff2 <- sum(actual.genome[2,-1])
#
# pvals<-c()
# for(j in 1:1000){
#   scaff1 <- sample(rep(db.chroms, times=Reduce(c,db.hits.list)),
#                    size=size.scaff1,
#                    replace=F)
#   null.dist <- c()
#   for(i in 1:1000){
#     null.dist[i] <- sum(sample(rep(db.chroms, times=Reduce(c,db.hits.list)),
#                                size=size.scaff1,
#                                replace=F)=="3L")
#   }
#   pvals[j] <- sum(null.dist>=sum(scaff1=="3L"))/1000
# }
#
# sum(pvals<=.05)

#run on chicken/anolis/rattlesnake
