## Source correlation
source("helper.R")
library(reshape2)
library(ggplot2)
library(circlize)
library(RColorBrewer)

# genes.of.interest=c("cMET", "EGFR", "ER", "ERCC1", "Her2/Neu", "MGMT", "MLH1", "MSH2", "MSH6", "PD-1", "PD-1 IHC", "PD-L1", "PD-L1 IHC", "PGP", "PMS2", "PR", "PTEN", "RPM1", "SPACR Monoclonal", "SPARC Polyclonal", "TLE3", "TOP2A", "TOPO1", "TS", "TUBB3")
genes.of.interest=sort(as.character(unique(Actionability$Gene)))

my.molecular1=my.molecular[my.molecular$Biomarker%in%genes.of.interest,]
# my.molecular1=my.molecular

## Create a new dataset for building heatmap
molecular.new<-cbind(my.molecular1$PatientID,my.molecular1$Biomarker,my.molecular1$Technology,my.molecular1$Positive)
colnames(molecular.new)<-c("PatientID","Biomarker","Technology","Positive")
molecular.new<-data.frame(molecular.new)
molecular.new$PatientID<-sort(as.numeric(molecular.new$PatientID))
molecular.new$Technology<-as.character(molecular.new$Technology)
molecular.new$Biomarker<-as.character(molecular.new$Biomarker)
molecular.new$Positive<-as.numeric(as.character(molecular.new$Positive))


gene<-as.character(unique(molecular.new$Biomarker))
id<-unique(molecular.new$PatientID)
data.heatmap<-matrix(ncol=length(gene)+2)
colnames(data.heatmap)<-c("ID",gene,"Technology")
for (i in (1:length(id))){
  d.new<-matrix(ncol=length(gene)+2)
  colnames(d.new)<-c("ID",gene,"Technology")
  
  id.data<-molecular.new[molecular.new$PatientID==i,]
  
  if (length(unique(id.data$Technology))==1){
    
    for (j in (1:nrow(id.data))){
      g<-id.data$Biomarker[j]
      
      d.new[,g]<-id.data$Positive[j]
      
    }
    d.new[,'Technology']<-unique(id.data$Technology)
    d.new[,'ID']<-i
    data.heatmap<-rbind(data.heatmap,d.new)
    
    
  }
  else{
    
    t<-unique(id.data$Technology)
    for (k in t){
      t.data<-id.data[id.data$Technology==k,]
      for (l in (1:nrow(t.data))){
        gl<-t.data$Biomarker[l]
        d.new[,gl]<-t.data$Positive[l]
        
      }
      d.new[,'Technology']<-k
      d.new[,'ID']<-i
      data.heatmap<-rbind(data.heatmap,d.new)
    }
    
  }
  
}


data.heatmap<-data.frame(data.heatmap[2:nrow(data.heatmap),])
## Second type dataset
data.heatmap[,1:(ncol(data.heatmap)-1)]<-apply(data.heatmap[,1:(ncol(data.heatmap)-1)],2,function(x) as.numeric(as.character(x)))
# data.heatmap<-data.heatmap[,2:ncol(data.heatmap)]

a<-split(data.heatmap,data.heatmap$Technology)
correlation<-list()
for (i in 1:length(a)){
  c1<-list()
  for (j in 1:length(a)){
    common.id<-intersect(a[[i]]$ID,a[[j]]$ID)
    
    if (sum(common.id)==0){
      c1[[j]]<-NA
    }
    else{
    c1[[j]]<-cor(a[[i]][a[[i]]$ID%in%common.id,2:17],a[[j]][a[[j]]$ID%in%common.id,2:17],use="pairwise.complete.obs")

    }
    }
  names(c1)<-names(a)
  correlation[[i]]<-c1
  
}
names(correlation)<-names(a)


name<-sort(unlist(unname(sapply(names(correlation),function(x) x[which(sum(!is.na(correlation[[x]]))!=0)]))))
color.division<-unique(unlist(correlation))
color.pos<-color.division[color.division>0]
color.neg<-color.division[color.division<0]
pos<-quantile(color.pos,na.rm=TRUE)
neg<-quantile(color.neg,na.rm=TRUE)

plot.data.pre<-function(x,y){
  
  plot.data=correlation[[x]][[y]]
  name<-colnames(plot.data)
  plot.data[is.na(plot.data)]<-0

  
  if (x==y){
    diag(plot.data)<-0
    ##change all values <= 0.6 in absolute value to 0
    # plot.data[abs(plot.data) <= 0.5] <- 0

    if (sum(plot.data)==0){
      warning('There are no L1 and L2 biomarkers for this combination or cannot calculate correlation')
    }
    else{
      col.name<-name[apply(plot.data,2,function(x) !all(x==0))]
      plot.data <- plot.data[,apply(plot.data,2,function(x) !all(x==0))]
      if(is.null(dim(plot.data))){

        plot.data<-matrix(plot.data)
        dimnames(plot.data)<-list(name,col.name)
        row.name<-name[apply(plot.data,1,function(x) !all(x==0))]
        plot.data <- plot.data[apply(plot.data,1,function(x) !all(x==0)),]
        plot.data<-matrix(plot.data)
      }
      else{
        row.name1<-colnames(plot.data)
        row.name<-colnames(plot.data)
        col.name<-name[apply(plot.data,1,function(x) !all(x==0))]
        plot.data <- plot.data[apply(plot.data,1,function(x) !all(x==0)),]}
      if(is.null(dim(plot.data))){
        plot.data<-matrix(plot.data)
        dimnames(plot.data)<-list(row.name1,col.name)
      }
      else{
        if(is.null(rownames(plot.data))){
          dimnames(plot.data)<-list(row.name,col.name)
        }
        else{
          plot.data<-plot.data
        }
      }

}}
  else{

    if (sum(plot.data)==0){
      warning('There are no L1 and L2 biomarkers for this combination or cannot calculate correlation')
    }
    else{
      
      col.name<-name[apply(plot.data,2,function(x) !all(x==0))]
      plot.data <- plot.data[,apply(plot.data,2,function(x) !all(x==0))]
      if(is.null(dim(plot.data))){
        
        plot.data<-matrix(plot.data)
        dimnames(plot.data)<-list(name,col.name)
        row.name<-name[apply(plot.data,1,function(x) !all(x==0))]
        plot.data <- plot.data[apply(plot.data,1,function(x) !all(x==0)),]
        plot.data<-matrix(plot.data)
      }
      else{
        row.name1<-colnames(plot.data)
        row.name<-colnames(plot.data)
        col.name<-name[apply(plot.data,1,function(x) !all(x==0))]
        plot.data <- plot.data[apply(plot.data,1,function(x) !all(x==0)),]}
      if(is.null(dim(plot.data))){
        plot.data<-matrix(plot.data)
        dimnames(plot.data)<-list(row.name1,col.name)
        plot.data
      }
      else{
      if(is.null(rownames(plot.data))){
          dimnames(plot.data)<-list(row.name,col.name)
          plot.data
        }
        else{
          plot.data<-plot.data
          plot.data
        }    }}}}

circle.plot2<-function(x,y) {
  
plot.data=plot.data.pre(x,y)

if (!is.matrix(plot.data)){
  warning(plot.data)
}
else{

cols <- rev(brewer.pal(8,"RdYlBu"))


matCols <- matrix("#FFFFFF", nrow=nrow(plot.data), ncol=ncol(plot.data))
matCols[plot.data < 0] <- cols[4]
matCols[plot.data < -neg[2]] <- cols[3]
matCols[plot.data < neg[3]] <- cols[2]
matCols[plot.data < neg[4]] <- cols[1]

matCols[plot.data > 0] <- cols[5]
matCols[plot.data > pos[2]] <- cols[6]
matCols[plot.data > pos[3]] <- cols[7]
matCols[plot.data > pos[4]] <- cols[8]

if (x==y){
  

    i=1.3
    par(mar=rep(i,4), mgp=c(i*2.6,i, i*1.5), las=2,cex=0.5*i)
    set.seed(2)
    bgcolor<-rand_color(dim(plot.data)*2)
    circos.par("track.height"=0.1,cell.padding=c(0,0,0,0))
    chordDiagram(plot.data, annotationTrack = "grid",symmetric = TRUE,col=matCols)
    circos.trackPlotRegion(track.index = 1,ylim = c(0, 1),track.height=0.6,bg.col=bgcolor,panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      sector.name = get.cell.meta.data("sector.index")
      circos.text(mean(xlim), 1.3, sector.name, facing = "clockwise",
                  niceFacing = TRUE, adj = c(0.05,0),cex=1)
    }, bg.border = NA)
}
else{
  

    colnames(plot.data)<-paste(colnames(plot.data),strsplit(x,split="[ \\|,+-]+")[[1]][1],sep=" ")
    rownames(plot.data)<-paste(rownames(plot.data),strsplit(y,split="[ \\|,+-]+")[[1]][1],sep=" ")
    
    set.seed(2)
    bgcolor<-rand_color(dim(plot.data)*2)
    circos.par("track.height"=0.1,cell.padding=c(0,0,0,0))
    chordDiagram(plot.data, annotationTrack = "grid",symmetric = FALSE,col=matCols)
     circos.trackPlotRegion(track.index = 1,ylim = c(0, 1),track.height=0.6,bg.col=bgcolor,panel.fun = function(x, y) {
       xlim = get.cell.meta.data("xlim")
       sector.name = get.cell.meta.data("sector.index")
       circos.text(mean(xlim), 1.2, sector.name, facing = "clockwise",
                   niceFacing = TRUE, adj = c(0.05,0),cex=0.7)}, bg.border = NA)
}}}# here set bg.border to NA is important





