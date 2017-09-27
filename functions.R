## Caris

plot.heatmap <- function(df, x, y, response,
                         category.x=NULL, category.y=NULL,
                         title='', xlab='Patient ID', ylab='Marker',
                         title.size=10, axis.size=10, strip.text.size=9, geom.size=4.5,facetangle=90,
                         show.all=TRUE, color="green_red",
                         legends=c('Over Expressed','Under Expressed'),extra.layers=NULL, ...){
  require(ggplot2)
  
  if (!(nrow(df) > 0)){
    stop('Data.frame is empty.')
  }
  
  if (!all(df[,response] %in% 0:1)){
    warning(paste("Response variable should be coded as numeric values",
                  "'1' for relevant and '0' for not relevant."))
  }
  
  col.select <- c(x, y, category.x, category.y)
  plot.data <- aggregate(list(Positive=df[,response]),
                         by=df[,col.select],
                         FUN=max, simplify=TRUE)
  
  ## Sort the biomarkers with the most mutated times
  ### Sort by individual categories
  # plot.data<-cbind(plot.data,Number=0)
  # plot.data.split<-split(plot.data,plot.data[,category.x])
  #   sum.gene<-by(plot.data,plot.data[,category.x],function(x) unlist(by(x,x[,y], function(x) sum(x[,response]))))
  #   #sum.gene<-sapply(sum.gene, function(x) names(sort(x,decreasing=TRUE)))
  #   
  #     for(j in 1:length(plot.data.split)){
  #      for(i in 1:nrow(plot.data.split[[j]])) {
  #        plot.data.split[[j]][,"Number"][i]<-sum.gene[[j]][plot.data.split[[j]][i,][,y]]
  #       }}
  #   #plot.data.split<-lapply(plot.data.split,function(x) x[order(x[,"Number"],decreasing=T),])
  # plot.data<-do.call(rbind.data.frame, plot.data.split)
  # plot.data<-plot.data[order(plot.data[,"Number"],decreasing=T),]
  
  ### Sort by combining all categories
  plot.data<-cbind(plot.data,Number=0)
  sum.gene<-by(plot.data,plot.data[,y],function(x) sum(x[,response]))
  for (i in 1:nrow(plot.data)){
    
    plot.data[i,][,"Number"]<-sum.gene[which(names(sum.gene)==plot.data[i,][,y])]
    
  }
  plot.data<-plot.data[order(plot.data[,"Number"],decreasing=T),]
  
  
  if (!show.all){
    show <- aggregate(list(Populated=plot.data$Positive),
                      by=data.frame(plot.data[,c(y, category.y)]),
                      FUN=max)
    show <- show[show$Populated > 0,]
    names(show)[length(c(y, category.y))] <- c(y, category.y)
    plot.data <- merge(plot.data, show[,c(y, category.y), drop=FALSE])
  }
  
  if (!is.null(category.y)){
    
    multi <- aggregate(list(Multiple=plot.data$Positive),
                       by=plot.data[,c(x, y, category.x)],
                       FUN=sum, simplify=TRUE)
    
    plot.data <- merge(plot.data, multi, all.x=TRUE)
    plot.data <- within(plot.data, Positive <- ifelse(Positive>0, Multiple, Positive))
    
  }
  #c("green","red","mediumorchid1","darkmagenta")
  plot.data$Ranges <- cut(plot.data$Positive, breaks=c(-Inf,0,1))
  plot.data$Ranges <-factor(plot.data$Ranges,levels=c("(0,1]",'(-Inf,0]'))
  plot.data[,x] <- factor(plot.data[,x])
  plot.data[,y] <- factor(plot.data[,y],levels=rev(unique(plot.data[,y])))
  
  ## Create legends 
  if (sum(plot.data$Positive)>0){
    legends=legends
  }
  else{
     legends=legends[2]
  }
  g <- ggplot(data=plot.data,aes_string(x=x,y=y)
  ) +
    
    geom_tile(colour='white', fill='white')+
    geom_point(aes(fill=Ranges), shape=22, col='white', size=geom.size) +
    # Change borders "gray90" to "white"
    
    scale_fill_manual(" ",values=c('(0,1]'=strsplit(color,"_")[[1]][1], '(-Inf,0]'=strsplit(color,"_")[[1]][2]), 
                       labels=legends) +
    
    
    theme(panel.background=element_blank(),panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    
    theme(axis.text.y=element_text(size=axis.size, color='black')) +
    theme(axis.text.x=element_text(size=axis.size, angle=45, hjust=1, color='black')) +
    theme(axis.title=element_text(size=title.size)) +
    
    theme(title=element_text(size=title.size)) +
    
    theme(legend.position='right') + 
    theme(legend.text=element_text(colour='black', size=axis.size)) +
    
    ggtitle(title) +
    xlab(xlab) +
    ylab(ylab) +
    
    theme(strip.text.x=element_text(size=strip.text.size, angle=90)) +
    guides(shape=guide_legend(override.aes=list(size=6)))
    
  
  
  ## Create a Labeller
  plot.data[,category.x]<-as.character(plot.data[,category.x])
  plot.data<-plot.data[!is.na(plot.data[,category.x]),]
  name<-unique(plot.data[,category.x])
  call<-c()
  for (i in name){
    num<-length(unique(plot.data[,x][plot.data[,category.x]==i]))
    call<-c(call,paste0(i,sep="\n(",num,sep=")"))
  }
  names(call)<-name
  
  
  if (length(c(category.x, category.y)) > 1){
    
    g <- g + facet_grid(paste(category.y, '~', category.x), space='free', scales='free') +
      theme(strip.text.x=element_text(size=strip.text.size, angle=facetangle),
            strip.text.y=element_text(size=strip.text.size, angle = 0))
    
  } else if(!is.null(category.x)){
    
    g <- g + facet_grid(paste('. ~', category.x),space='free', scales='free',labeller=as_labeller(call)) +
      theme(strip.text.x=element_text(size=strip.text.size, angle=facetangle))
    
  } else if(!is.null(category.y)){
    
    g <- g + facet_grid(paste(category.y, '~ .'), space='free', scales='free') +
      theme(strip.text.y=element_text(size=strip.text.size, angle = 0))
  }
  
  if(!is.null(extra.layers)) g <- g + extra.layers
  
  return(g)
}
