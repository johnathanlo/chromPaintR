GraphChroms <- function(chromdat,
                        scaff.lim = T,
                        dbchrom.lim = T,
                        sample = 1,
                        max.scaffs = 10,
                        max.dbchroms = 10){###have user input colors
  require(ggplot2);require(ggpubr);require(RColorBrewer);require(dplyr)
  alldbhits <- Reduce(c, chromdat$summary$`Database hits`)
  allscaffpos <- Reduce(c, chromdat$summary$`Locations on query scaffold`)
  alldbpos <- Reduce(c, chromdat$summary$`Locations on database`)
  hits_per_scaff <- lapply(chromdat$summary$`Database hits`, length)
  allscaffs <- Reduce(c,
                      mapply(rep, rownames(chromdat$summary),
                             each = hits_per_scaff))
  graph_scaffs <- data.frame(scaffolds = allscaffs,
                             chromosomes = alldbhits,
                             positions = allscaffpos,
                             origins = alldbpos)
  ordered.chroms <- chromdat$lengthdata[[4]][order(as.numeric(chromdat$lengthdata[[3]]), decreasing = T)]
  if(dbchrom.lim){
    graph_scaffs<-graph_scaffs[graph_scaffs$chromosomes %in% ordered.chroms[1:max.dbchroms],]
  }
  if(scaff.lim){
    largest.scaffs <- chromdat$lengthdata[[1]][order(as.numeric(chromdat$lengthdata[[2]]), decreasing = T)[1:max.scaffs]]
    graph_scaffs <- graph_scaffs[graph_scaffs$scaffolds %in% largest.scaffs,]
  }
  graph_scaffs$chromosomes_num <-as.numeric(factor(graph_scaffs$chromosomes))
  graph_scaffs$scaffolds_num <- as.numeric(factor(graph_scaffs$scaffolds))
  scaff_lengths <- c()
  for(scaff in graph_scaffs$scaffolds){
    scaff_lengths <- c(scaff_lengths, chromdat$summary$`Scaffold Lengths`[[scaff]])
  }
  graph_scaffs$scafflengths <- scaff_lengths


  dbchrom_locarray <- Reduce(c,mapply(seq, from = 0, to = as.numeric(chromdat$lengthdata[[3]][1:max.dbchroms]), by=10000))
  dbchrom_namearray <- Reduce(c,mapply(rep, chromdat$lengthdata[[4]][1:max.dbchroms], each = ceiling(as.numeric(chromdat$lengthdata[[3]][1:max.dbchroms])/10000)))
  graph_dbchroms <- data.frame(chromosomes = dbchrom_namearray, positions = dbchrom_locarray)
  graph_dbchroms$chromosomes_num <- as.numeric(factor(graph_dbchroms$chromosomes))

  if(sample < 1){
    graph_scaffs <- sample_frac(graph_scaffs, size = sample)
  }
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  colors <- sample(max.dbchroms, col_vector)
  #colors <- brewer.pal(max.dbchroms, "Paired")
  names(colors) <- levels(graph_dbchroms$Chromosomes)
  colorScale <- scale_color_manual(name = "Chromosomes", values = colors)
  if(length(unique(graph_scaffs$scaffolds))<max.scaffs){
    cp2.breaks <- length(unique(graph_scaffs$scaffolds))
  }else{cp2.breaks=max.scaffs}

  chromplot <- ggplot(data=graph_dbchroms, aes(x= chromosomes_num, y =positions, color = chromosomes )) +
    geom_segment(aes(x= chromosomes_num-.1, y = positions, xend = chromosomes_num +.1, yend = positions), size = .5) +
    scale_x_continuous(name = "Chromosomes", breaks = 1:max.dbchroms, labels = levels(factor(graph_dbchroms$chromosomes)))+
    scale_y_continuous(name = "Length", labels = fancy_scientific)+
    colorScale +
    coord_flip()+
    guides(colour = guide_legend(override.aes = list(size = max.dbchroms)))

  chromplot2 <- ggplot(data = graph_scaffs, aes(x=scaffolds_num, y = positions, color = chromosomes)) +
    geom_rect(data = graph_scaffs, aes(xmin = scaffolds_num-.1, xmax = scaffolds_num+.1, ymin = 0, ymax = scafflengths, fill = I("white")), color = NA)+
      geom_segment(aes(x=scaffolds_num-.1, y = positions, xend = scaffolds_num+.1, yend = positions), size = .5) +
      scale_x_continuous(name = "Scaffolds", breaks = 1:cp2.breaks, labels = levels(factor(graph_scaffs$scaffolds)))+
      scale_y_continuous(name = "Positions", labels = fancy_scientific)+
      colorScale +
      coord_flip() +
      theme(legend.position = "none")

  ggarrange(chromplot2, chromplot, labels = c("Query Scaffolds", "Database chromosomes"))
}

