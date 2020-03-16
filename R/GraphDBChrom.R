GraphScaff<-function(chromdat, dbchrom){
  require(ggplot2);require(ggpubr);require(RColorBrewer);require(dplyr)
  scaffs <- c()
  dblocs <- c()
  chrom.length <- as.numeric(chromdat$lengthdata[[3]][grep(dbchrom, chromdat$lengthdata[[4]])])
  numscaffs <- length(chromdat$summary$`Database hits`)
  hitsperscaff <- lapply(chromdat$summary$`Database hits`, length)
  for(i in 1:numscaffs){
    for(j in 1:hitsperscaff[[i]])
    if(chromdat$summary$`Database hits`[[i]][j] == dbchrom){
      scaffs <- c(scaffs, names(chromdat$summary$`Database hits`)[i])
      dblocs <- c(dblocs, chromdat$summary$`Locations on database`[[i]][j])
    }
  }
  
  chrom <- data.frame(scaffs = scaffs, dblocs = dblocs)
  
  uniquescaffs <- unique(chrom$scaffs)
  colorScale <- scale_color_manual(name = "Chromosomes", values = 1:length(uniquescaffs))
  
  chromplot <- ggplot(data=chrom, aes(x = 1, y=dblocs, color = scaffs )) +
    geom_rect(data = chrom, aes(xmin = .9, xmax = 1.1, ymin = 0, ymax = chrom.length, fill = I("white")), color = NA) +
    geom_segment(aes(x= .9, y = dblocs, xend = 1.1, yend = dblocs), size = .5) +
    scale_x_continuous(name = "Scaffold", limits = c(0,2))+
    scale_y_continuous(name = "Position", labels = fancy_scientific)+
    guides(colour = guide_legend(override.aes = list(size = 10)))
  chromplot
}