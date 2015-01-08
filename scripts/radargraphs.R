library(plotrix)
tissue <- read.csv("tissue_expression.csv",sep="\t",row.names=1)
tissue.labels<-rownames(tissue)


for (gene in names(tissue)){

  filename <- paste(gene,"_tissue_plot.jpg",sep="")
  jpeg(file=filename,width = 120, height = 120);
  polar.plot(log10(tissue[,gene]+1),,testpos,labels="",start=90,clockwise=TRUE,rp.type="poly",show.grid.labels=F,line.col="blue",lwd=5)
  dev.off()

  filename <- paste(gene,"_tissue_big_plot.jpg",sep="")
  jpeg(file=filename,width = 720, height = 640);
  tissue.labels[14]<-"\n\nMale accessory glands"
  polar.plot(log10(tissue[,gene]+1),,testpos,labels=tissue.labels,start=90,clockwise=TRUE,rp.type="poly",show.grid.labels=F,line.col="blue",lwd=2,main=gene)

  dev.off()
  
  filename <- paste(gene,"_tissue_scaled_big_plot.jpg",sep="")
  jpeg(file=filename,width = 720, height = 640);
  polar.plot(log10(tissue[,gene]+1),,testpos,labels=tissue.labels,start=90,clockwise=TRUE,rp.type="poly",show.grid.labels=F,line.col="blue",lwd=2,main=gene,radial.lim=pretty(c(0,4.5)))

  dev.off();  
}


