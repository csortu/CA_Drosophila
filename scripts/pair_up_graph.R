#cahlist<-c(1:9,11:15)
#cahlist<-c(13,15)
cahlist<-c(6,16)

for (i in cahlist){
  for (j in cahlist){
      if(i <j){
      	   filename<-paste("STAT_",i,"_ka_ks.txt",sep="");
      	   dataset1<-read.table(filename,sep="\t",header=T);
	   splist1<-levels(dataset1$SPECIES_NO.1);
      	   filename<-paste("STAT_",j,"_ka_ks.txt",sep="");
      	   dataset2<-read.table(filename,sep="\t",header=T);
	   splist2<-levels(dataset2$SPECIES_NO.1);
	   splist.common<-intersect(splist1,splist2);

           gene1 <- paste("CAH",i,sep="");
           gene2 <- paste("CAH",j,sep="");
           
           splist.combo <- combn(splist.common,2,simplify=T);

           kaks.set1<-vector();
           kaks.set2<-vector();
           
           ka.set1<-vector();
           ka.set2<-vector();
           
           ks.set1<-vector();
           ks.set2<-vector();
           
           z.score.set1<-vector();
           z.score.set2<-vector();
           
	   for(k in 1:(length(splist.combo)/2)){

	   	 sp1<-splist.combo[2*k -1]
	   	 sp2<-splist.combo[2*k]

                 kaks1<-0;
                 kaks2<-0;

                 ka1<-0;
                 ka2<-0;
                 
                 ks1<-0;
                 ks2<-0;
                 
                 z.score1<-0;
                 z.score2<-0;

		 if (sum(dataset1$SPECIES_NO.1==sp1 & dataset1$SPECIES_NO.2==sp2)>0){
		    kaks1<-dataset1[(dataset1$SPECIES_NO.1==sp1 & dataset1$SPECIES_NO.2==sp2),"Ka_Ks"][1];
		    ka1<-dataset1[(dataset1$SPECIES_NO.1==sp1 & dataset1$SPECIES_NO.2==sp2),"D_n"][1];
		    ks1<-dataset1[(dataset1$SPECIES_NO.1==sp1 & dataset1$SPECIES_NO.2==sp2),"D_s"][1];
		    z.score1<-dataset1[(dataset1$SPECIES_NO.1==sp1 & dataset1$SPECIES_NO.2==sp2),"Z.SCORE"][1];
                  }
		 if (sum(dataset1$SPECIES_NO.1==sp2 & dataset1$SPECIES_NO.2==sp1)>0){
		    kaks1<-dataset1[(dataset1$SPECIES_NO.1==sp2 & dataset1$SPECIES_NO.2==sp1),"Ka_Ks"][1];
		    ka1<-dataset1[(dataset1$SPECIES_NO.1==sp2 & dataset1$SPECIES_NO.2==sp1),"D_n"][1];
		    ks1<-dataset1[(dataset1$SPECIES_NO.1==sp2 & dataset1$SPECIES_NO.2==sp1),"D_s"][1];
		    z.score1<-dataset1[(dataset1$SPECIES_NO.1==sp2 & dataset1$SPECIES_NO.2==sp1),"Z.SCORE"][1];
		 }

		 if (sum(dataset2$SPECIES_NO.1==sp1 & dataset2$SPECIES_NO.2==sp2)>0){
		    kaks2<-dataset2[(dataset2$SPECIES_NO.1==sp1 & dataset2$SPECIES_NO.2==sp2),"Ka_Ks"][1];
		    ka2<-dataset2[(dataset2$SPECIES_NO.1==sp1 & dataset2$SPECIES_NO.2==sp2),"D_n"][1];
		    ks2<-dataset2[(dataset2$SPECIES_NO.1==sp1 & dataset2$SPECIES_NO.2==sp2),"D_s"][1];
		    z.score2<-dataset2[(dataset2$SPECIES_NO.1==sp1 & dataset2$SPECIES_NO.2==sp2),"Z.SCORE"][1];
		 }
		 if (sum(dataset2$SPECIES_NO.1==sp2 & dataset2$SPECIES_NO.2==sp1)>0){
		    kaks2<-dataset2[(dataset2$SPECIES_NO.1==sp2 & dataset2$SPECIES_NO.2==sp1),"Ka_Ks"][1];
		    ka2<-dataset2[(dataset2$SPECIES_NO.1==sp2 & dataset2$SPECIES_NO.2==sp1),"D_n"][1];
		    ks2<-dataset2[(dataset2$SPECIES_NO.1==sp2 & dataset2$SPECIES_NO.2==sp1),"D_s"][1];
		    z.score2<-dataset2[(dataset2$SPECIES_NO.1==sp2 & dataset2$SPECIES_NO.2==sp1),"Z.SCORE"][1];
		 }

                 if (kaks1 & kaks2){
#                   print(paste(sp1,sp2,"KaKs1:",kaks1,"KaKs2:",kaks2));
                   kaks.set1<-c(kaks.set1,kaks1);
                   kaks.set2<-c(kaks.set2,kaks2);

                   ks.set1<-c(ks.set1,ks1);
                   ks.set2<-c(ks.set2,ks2);

                   ka.set1<-c(ka.set1,ka1);
                   ka.set2<-c(ka.set2,ka2);

                   z.score.set1<-c(z.score.set1,z.score1);
                   z.score.set2<-c(z.score.set2,z.score2);

                 }                 

               }
# Plotting the things

                 jfname <- paste("correlations_gene_",i,"_",j,".jpg",sep="");
                 jpeg(filename = jfname, width = 960, height = 960, pointsize = 22, quality = 100,bg = "white");
                 par(mfrow=c(2,2)) 
                 
# Ka Ks graph
                 plot(kaks.set1~kaks.set2,xlab=gene1,ylab=gene2,main="Ka/Ks correlation");
                 slm <- lm(kaks.set1~kaks.set2-1);
                 abline(slm);
                 abline(a=0,b=1,lty=3);
                 
# Ka graph
                 plot(ka.set1~ka.set2,xlab=gene1,ylab=gene2,main="Ka correlation");
                 slm <- lm(ka.set1~ka.set2-1);
                 abline(slm);
                 abline(a=0,b=1,lty=3);
                 
# Ks graph
                 plot(ks.set1~ks.set2,xlab=gene1,ylab=gene2,main="Ks correlation");
                 slm <- lm(ks.set1~ks.set2-1);
                 abline(slm);
                 abline(a=0,b=1,lty=3);
                 
# z_score graph
                 plot(z.score.set1~z.score.set2,xlab=gene1,ylab=gene2,main="Z.Score correlation");
                 slm <- lm(z.score.set1~z.score.set2-1);
                 abline(slm);
                 abline(a=0,b=1,lty=3);


                 dev.off();
                 
	   


      }
  }  
}


