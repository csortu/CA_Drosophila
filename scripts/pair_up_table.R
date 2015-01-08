cahlist<-c(1:9,11:15)
#cahlist<-c(13,15)
#cahlist<-c(13,2)

print(paste("gene1", "gene2", "kaks cor", "kaks conf int 1", "kaks conf int 2", "kaks p-value", "ka cor", "ka conf int 1", "ka conf int 2", "ka p-value", "ks cor", "ks conf int 1", "ks conf int 2", "ks p-value", "z.score cor", "z.score conf int 1", "z.score conf int 2", "z.score p-value"));

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

          kaks1.norm<-0;
          kaks2.norm<-0;
          
          ka1.norm<-0;
          ka2.norm<-0;
                 
          ks1.norm<-0;
          ks2.norm<-0;
                 
          z.score1.norm<-0;
          z.score2.norm<-0;

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
# print statistics
# gene1, gene2, kaks cor, kaks conf int 1, kaks conf int 2, kaks p-value, ka cor, ka conf int 1, ka conf int 2, ka p-value, ks cor, ks conf int 1, ks conf int 2, ks p-value, z.score cor, z.score conf int 1, z.score conf int 2, z.score p-value

        kaks.cor <- cor.test(kaks.set1, kaks.set2)
        ka.cor <- cor.test(ka.set1, ka.set2)
        ks.cor <- cor.test(ks.set1, ks.set2)
        z.score.cor <- cor.test(z.score.set1, z.score.set2)
                 
        print(paste(gene1,gene2,kaks.cor$estimate,kaks.cor$conf.int[1],kaks.cor$conf.int[2],kaks.cor$p.value,ka.cor$estimate,ka.cor$conf.int[1],ka.cor$conf.int[2],ka.cor$p.value,ks.cor$estimate,ks.cor$conf.int[1],ks.cor$conf.int[2],ks.cor$p.value,z.score.cor$estimate,z.score.cor$conf.int[1],z.score.cor$conf.int[2],z.score.cor$p.value,sep=" "));
                 
        

      }
    }  
}


# read back the formatted data like this:
# cortable<-read.csv("pair_up_table_formatted.csv",sep="\t",header=T)
