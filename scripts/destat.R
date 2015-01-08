cortable<-read.csv("pair_up_table_formatted.csv",sep="\t",header=T)
cahlist<-c(1:9,11:16)

genes <- vector()

kaksmean <- vector()
kakssd <- vector()

kamean <- vector()
kasd <- vector()

ksmean <- vector()
kssd <- vector()

zscmean <- vector()
zscsd <- vector()

kaksmax <- vector()

for (i in cahlist){

  filename <- paste("STAT_",i,"_ka_ks.txt",sep="");
  dataset1 <- read.table(filename,sep="\t",header=T);

  kaksmean <- c(kaksmean,mean(unlist(dataset1$Ka_Ks)));
  kamean <- c(kamean,mean(unlist(dataset1$D_n)));
  ksmean <- c(ksmean,mean(unlist(dataset1$D_s)));
  zscmean <- c(zscmean,mean(unlist(dataset1$Z.SCORE)));

  kakssd <- c(kakssd,sd(unlist(dataset1$Ka_Ks)));
  kasd <- c(kasd,sd(unlist(dataset1$D_n)));
  kssd <- c(kssd,sd(unlist(dataset1$D_s)));
  zscsd <- c(zscsd,sd(unlist(dataset1$Z.SCORE)));

  kaksmax <- c(kaksmax,max(unlist(dataset1$Ka_Ks)));

  
  genes <- c(genes,paste("CAH",i,sep=""));
}

df <- data.frame(kaks.mean=kaksmean,kaks.sd=kakssd,kaks.max=kaksmax,ka.mean=kamean,ka.sd=kasd,ks.mean=kamean,ks.sd=kssd,z.score.mean=zscmean,z.score.sd=zscsd,row.names=genes);


write.csv(df,file="kaks_descriptive_stats.csv");

