AnalyzeBlastR<-function(data = "results_v03.out",
                        chromlengthfile = "chromlength.out",
                        min.chain = 5,
                        allhits = F
                        ){
  if(file.exists(data) && file.exists(chromlengthfile)&&min.chain>0){
    data<-read.table(data, stringsAsFactors = F)
    colnames(data) <- c("qloc", "scaffname", "dbloc", "dbchrom", "length", "pident", "score", "evalue")
    chained.hits <- 0
    scaffnames <- c()
    dbseqs <- c()
    qlocs <- c()
    dblocs <- c()
    divs <- c()
    numscaffs <- 0
    totalhits <- 0
    flag <- 0
  
    if(allhits){
      for(row in 1:(nrow(data)-1)){
        totalhits <- totalhits +1
        dbseqs <- c(dbseqs, data[row, "dbchrom"])
        qlocs <- c(qlocs, data[row, "qloc"])
        dblocs <- c(dblocs, data[row,"dbloc"])
        if(data[row, "scaffname"] != data[row+1, "scaffname"]){
          scaffnames <- c(scaffnames, data[row,"scaffname"])
          divs <- c(divs, totalhits)
          numscaffs <- numscaffs +1
        }else if(row == (nrow(data)-1)){
            scaffnames <- c(scaffnames, data[row, "scaffname"])
            divs <- c(divs, totalhits)
            numscaffs <- numscaffs +1
          }
      }
    }else{
      for(row in 1:(nrow(data)-1)){
        if(data[row,"scaffname"] == data[row+1,"scaffname"] &&
           data[row,"dbchrom"] == data[row+1, "dbchrom"] &&
           data[row, "qloc"]!=data[row+1, "qloc"]){
  
          chained.hits <- chained.hits +1
  
          if(chained.hits >= min.chain){
            flag <- 1
            totalhits <- totalhits +1
            dbseqs <- c(dbseqs, data[row, "dbchrom"])
            qlocs <- c(qlocs, data[row, "qloc"])
            dblocs <- c(dblocs, data[row,"dbloc"])
          }
  
        }
        else if(data[row,"scaffname"] != data[row+1, "scaffname"] & flag == 1){
          scaffnames <- c(scaffnames, data[row,"scaffname"])
          divs <- c(divs, totalhits)
          numscaffs <- numscaffs +1
          flag <- 0
          chained.hits <- 0
        }
        else if(row == (nrow(data)-1) && flag == 1){
          scaffnames <- c(scaffnames, data[row, "scaffname"])
          divs <- c(divs, totalhits)
          numscaffs <- numscaffs +1
          flag <- 0
          chained.hits <- 0
        }
        else{
          chained.hits<-0
        }
      }
    }
  
    lengthdata <- scan(file = chromlengthfile, what = character(),sep = "\n")
    lengthdata <- sapply(lengthdata, split = "\t", strsplit)
    scaffstrs <- sapply(lengthdata[[1]], split = " ", strsplit)
    for(i in 1:length(lengthdata[[1]])){
      lengthdata[[1]][i] <- scaffstrs[[i]][1]
    }
    scaff_lengths <- data.frame(names = lengthdata[[1]], lengths = as.numeric(lengthdata[[2]]))
    data.list <- list(`Number of scaffolds` = numscaffs,
                      `Scaffold Names` = scaffnames,
                      `DB seq names` = dbseqs,
                      `Query seq nums` = qlocs,
                      `DB seq nums` = dblocs,
                      `Divisions` = divs,
                      `Query scaffold lengths` = as.numeric(lengthdata[[2]]))
  
    dbseqnames = c()
    qseqnums = c()
    dbseqnums = c()
    j = 1
    for(i in data.list$Divisions){
      dbseqnames <- c(dbseqnames, list(data.list$`DB seq names`[j:i]))
      qseqnums <- c(qseqnums, list(data.list$`Query seq nums`[j:i]))
      dbseqnums <- c(dbseqnums, list(data.list$`DB seq nums`[j:i]))
      j = i+1
    }
    i = length(data.list$`DB seq names`)
    dbseqnames <- c(dbseqnames, list(data.list$`DB seq names`[j:i]))
    qseqnums <- c(qseqnums, list(data.list$`Query seq nums`[j:i]))
    dbseqnums <- c(dbseqnums, list(data.list$`DB seq nums`[j:i]))
  
  
    data.list.summary <-cbind(vector(mode = "list", length=data.list$`Number of scaffolds`),
                      vector(mode = "list", length=data.list$`Number of scaffolds`),
                      vector(mode = "list", length=data.list$`Number of scaffolds`),
                      vector(mode = "numeric", length = data.list$`Number of scaffolds`),
                      vector(mode = "character", length = data.list$`Number of scaffolds`),
                      vector(mode = "numeric", length = data.list$`Number of scaffolds`))
    colnames(data.list.summary) <- c("Database hits", "Locations on query scaffold", "Locations on database", "Scaffold Lengths", "Chromosome match", "Match probability")
    rownames(data.list.summary) <- data.list$`Scaffold Names`
    data.list.summary <- as.data.frame(data.list.summary)
  
  
  
    for(i in seq(1, data.list$`Number of scaffolds`, by = 1)){
      data.list.summary$`Database hits`[i] <- dbseqnames[i]
      data.list.summary$`Locations on query scaffold`[i] <- qseqnums[i]
      data.list.summary$`Locations on database`[i] <- dbseqnums[i]
      data.list.summary$`Scaffold Lengths`[i] <- as.numeric(scaff_lengths[scaff_lengths$names == rownames(data.list.summary)[i],2])
    }
    data.complete <- list(summary = data.list.summary, lengthdata = lengthdata)
    return(data.complete)
  }else{
    print("One or more parameters are faulty.")
  }
}
