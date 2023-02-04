######################################################
## Project: BACTERIAL HISTONES
## Compare a.a similarities within bdellovibrio species
## 
## 
## Aim : 
## compute all best reciprocal hits beteen 1 specie and the rest
## 
## 
## Date: March 22
## Author: Antoine Hocher
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.

Sys.setenv(PATH="/Users/ahocher/opt/miniconda3/bin:/Users/ahocher/opt/miniconda3/condabin:/opt/local/bin:/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin")
source("/Users/ahocher/Dropbox/Scripts/PLOT_LIBRARY/AH_plot_library.R")
setwd("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/BDELLOVIBRIO/RECIPROCALBLAST")
#Libraries : 
library(rtracklayer)
library(metablastr)
ProteomeDir="/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/BDELLOVIBRIO/PROTEINS/AnnotatedGenomes/"
# #For each proteome, compute the best reciprocal hits : 
setwd(ProteomeDir)

Genomes=read.table("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/BDELLOVIBRIO/PROTEINS/AnnotatedGenomes/Annotated_genomes_info.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

Prot_available=list.files(pattern = ".faa$")

for (Organism in Prot_available){

  Assembly=paste(strsplit(Organism,"_")[[1]][1:2],collapse = "_")
  Specie=Genomes[which(Genomes$Assembly==Assembly),]$Specie


     for(Organism2 in Prot_available){
       Assembly2=paste(strsplit(Organism2,"_")[[1]][1:2],collapse = "_")
       Specie2=Genomes[which(Genomes$Assembly==Assembly2),]$Specie
       

         if(Organism!=Organism2){
           if(file.exists(paste("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/BDELLOVIBRIO/RECIPROCALBLAST/",Assembly,"_vs_",Assembly2,"__",Specie,"_vs_",Specie2,"_rblast_eval_1e3.txt",sep=""))==F){


           Rhits=blast_best_reciprocal_hit(query=Organism,subject = Organism2,evalue = 0.001 ,search_type = "protein_to_protein",task = "blastp",is.subject.db = T,cores = 4)

          write.table(Rhits,paste("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/BDELLOVIBRIO/RECIPROCALBLAST/",Assembly,"_vs_",Assembly2,"__",Specie,"_vs_",Specie2,"_rblast_eval_1e3.txt",sep=""),row.names = F,sep="\t",quote=F)
           }
      }
  }
}






#2nd step : %identity in Bdellovibrio bacteriovorus HD100 vs others : 


setwd("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/BDELLOVIBRIO/RECIPROCALBLAST/")
FilesBRB=list.files(pattern="GCA_000196175.1_vs*")

BRB=read.table(FilesBRB[1],header=T,sep="\t",quote="",stringsAsFactors = F)
BRB=BRB[,1:3]
names(BRB)[2:3]=paste( names(BRB.2)[2:3],paste(strsplit(FilesBRB[1],"_")[[1]][4:5],collapse="_"),sep="_")

for(i in FilesBRB[2:length(FilesBRB)]){
  BRB.2=read.table(i,header=T,sep="\t",quote="",stringsAsFactors = F)
  BRB.2=BRB.2[,which(names(BRB.2)%in%c("query_id","subject_id","perc_identity"))]
  names(BRB.2)[2:3]=paste( names(BRB.2)[2:3],paste(strsplit(i,"_")[[1]][4:5],collapse="_"),sep="_")
  BRB=merge(BRB,BRB.2,by="query_id",all.x=T,all.y=T)
}


write.table(BRB,file="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/AADIVERGENCE/BdellovibrioHD100_vs_others_BRBH.txt",row.names=F,quote=F,sep="\t")


#Reload : 

BRB=read.delim("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/AADIVERGENCE/BdellovibrioHD100_vs_others_BRBH.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

#Compute average IDentity , and see where our gene of interest lies within the distribution of all identities : 

BRB$AvgID=apply(BRB[,grep(names(BRB),pattern="perc_identity_")],1,function(x) mean(x,na.rm=T))

BRB$HF=ifelse(BRB$query_id=="959_Bdellovibrio^bacteriovorus^HD100@CAE77736.1$GCA_000196175.1",1,0)
BRB$ZF_HF=ifelse(BRB$query_id=="959_Bdellovibrio^bacteriovorus^HD100@CAE80810.1$GCA_000196175.1",1,0)


Values=BRB[which(BRB$HF==1),grep(names(BRB),pattern="perc_identity_")]
AssemblyToRemove=names(Values[which(is.na(Values)==T)])
AssemblyToRemove=unlist(lapply(AssemblyToRemove,function(x) strsplit(x,"_")[[1]][4]))
AssemblyToRemove=paste("perc_identity_GCA_",AssemblyToRemove,sep="")
BRB.r=BRB[,-which(names(BRB)%in%AssemblyToRemove)]

BRB.r$AvgID=apply(BRB.r[,grep(names(BRB.r),pattern="perc_identity_")],1,function(x) mean(x,na.rm=T))


library(ggplot2);library(ggpubr)
Example=ggplot(data=BRB)+geom_violin(aes(x=1,y=AvgID))+geom_point(aes(x=1,y=BRB[which(BRB$HF==1),]$AvgID))+ylab("Average aa identity")+theme_pubclean()+theme(aspect.ratio=3)+ylim(c(0,100))

ggsave(Example,filename="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/AADIVERGENCE/BdellovibrioHD100_vs_others_BRBH.pdf")

Example=ggplot(data=BRB)+geom_violin(aes(x=1,y=AvgID))+geom_point(aes(x=1,y=BRB[which(BRB$ZF_HF==1),]$AvgID))+ylab("Average aa identity")+theme_pubclean()+theme(aspect.ratio=3)+ylim(c(0,100))

ggsave(Example,filename="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/AADIVERGENCE/BdellovibrioHD100_vs_others_BRBH_ZFHFBd3044.pdf")


#Only on the subset with values for histones (cleaner I suppose)
Example.r=ggplot(data=BRB.r)+geom_violin(aes(x=1,y=AvgID))+geom_point(aes(x=1,y=BRB.r[which(BRB.r$HF==1),]$AvgID))+ylab("Average aa identity")+theme_pubclean()+theme(aspect.ratio=3)+ylim(c(0,100))
ggsave(Example.r,filename="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/AADIVERGENCE/BdellovibrioHD100_vs_others_BRBH_Species_With_aBRBHagainstBd0055.pdf")

#adding Bd3044 : 
Example.r=ggplot(data=BRB.r)+geom_violin(aes(x=1,y=AvgID))+geom_point(aes(x=1,y=BRB.r[which(BRB.r$HF==1),]$AvgID))+geom_point(aes(x=1,y=BRB.r[which(BRB.r$ZF_HF==1),]$AvgID),col="blue")+ylab("Average aa identity")+theme_pubclean()+theme(aspect.ratio=3)+ylim(c(0,100))
ggsave(Example.r,filename="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/AADIVERGENCE/BdellovibrioHD100_vs_others_BRBH_Species_With_aBRBHagainstBd0055_with_bd3044_blue.pdf")






#Only on the subset of proteins with BRBH in all species: 

BRB.r$NbNA=apply(BRB.r[,grep(names(BRB.r),pattern="perc_identity_")],1,function(x) length(which(is.na(x)==T)))


#Graph part : 
Example.r2=ggplot(data=BRB.r2)+geom_violin(aes(x=1,y=AvgID))+geom_point(aes(x=1,y=BRB.r2[which(BRB.r2$HF==1),]$AvgID))+ylab("Average aa identity")+theme_pubclean()+theme(aspect.ratio=3)+ylim(c(0,100))

ggsave(Example.r2,filename="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/AADIVERGENCE/BdellovibrioHD100_vs_others_BRBH_Species_With_aBRBHagainstBd0055_and_BRBinallSpecies.pdf")


#Trying cumulative frequency plot : 

BRBorder=BRB.r[order(BRB.r$perc_identity_GCA_000317895.1),]
BRBorder=BRBorder[which(is.na(BRBorder$perc_identity_GCA_000317895.1)==F),]

BRBorder2=BRB.r[order(BRB.r$perc_identity_GCA_000786105.1),]
BRBorder2=BRBorder2[which(is.na(BRBorder2$perc_identity_GCA_000786105.1)==F),]

BRBorder3=BRB.r[order(BRB.r$perc_identity_GCA_002773975.1),]
BRBorder3=BRBorder3[which(is.na(BRBorder3$perc_identity_GCA_002773975.1)==F),]


ThreeSpecieExample=ggplot() + geom_step(data=BRBorder3,aes(x=perc_identity_GCA_002773975.1,y=..y..),stat="ecdf",col=col_vector[3],size=1)+ geom_step(data=BRBorder2,aes(x=perc_identity_GCA_000786105.1,y=..y..),stat="ecdf",col=col_vector[2],size=1)+ geom_step(data=BRBorder,aes(x=perc_identity_GCA_000317895.1,y=..y..),stat="ecdf",col=col_vector[1],size=1)+geom_point(data=BRBorder[which(BRBorder$HF==1),],aes(y=which(BRBorder$HF==1)/dim(BRBorder)[1],x=perc_identity_GCA_000317895.1),col=col_vector[1],size=2)+geom_point(data=BRBorder2[which(BRBorder2$HF==1),],aes(y=which(BRBorder2$HF==1)/dim(BRBorder2)[1],x=perc_identity_GCA_000786105.1),col=col_vector[2],size=2)+geom_point(data=BRBorder3[which(BRBorder3$HF==1),],aes(y=which(BRBorder3$HF==1)/dim(BRBorder3)[1],x=perc_identity_GCA_002773975.1),col=col_vector[3],size=2)+theme_pubr()+theme(aspect.ratio = 1/3)+xlab(" best reciprocal blast hits a.a identity (%)")+ylab("Fraction")+ scale_y_continuous(labels = function(x) paste0(x*100, "%"))
ggsave(ThreeSpecieExample,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/AADIVERGENCE/BdellovibrioHD100_vs_3specie_Focus_HF.pdf")


ThreeSpecieExample=ggplot() + geom_step(data=BRBorder3,aes(x=perc_identity_GCA_002773975.1,y=..y..),stat="ecdf",col=col_vector[3],size=1)+ geom_step(data=BRBorder2,aes(x=perc_identity_GCA_000786105.1,y=..y..),stat="ecdf",col=col_vector[2],size=1)+ geom_step(data=BRBorder,aes(x=perc_identity_GCA_000317895.1,y=..y..),stat="ecdf",col=col_vector[1],size=1)+geom_point(data=BRBorder[which(BRBorder$ZF_HF==1),],aes(y=which(BRBorder$ZF_HF==1)/dim(BRBorder)[1],x=perc_identity_GCA_000317895.1),col=col_vector[1],size=2)+geom_point(data=BRBorder2[which(BRBorder2$ZF_HF==1),],aes(y=which(BRBorder2$ZF_HF==1)/dim(BRBorder2)[1],x=perc_identity_GCA_000786105.1),col=col_vector[2],size=2)+geom_point(data=BRBorder3[which(BRBorder3$ZF_HF==1),],aes(y=which(BRBorder3$ZF_HF==1)/dim(BRBorder3)[1],x=perc_identity_GCA_002773975.1),col=col_vector[3],size=2)+theme_pubr()+theme(aspect.ratio = 1/3)+xlab(" best reciprocal blast hits a.a identity (%)")+ylab("Fraction")+ scale_y_continuous(labels = function(x) paste0(x*100, "%"))
ggsave(ThreeSpecieExample,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/AADIVERGENCE/BdellovibrioHD100_vs_3specie_Focus_ZF_HF.pdf")



#Specie name : 
"Bdellovibrio sp. CG10"
"Bdellovibrio sp. ArHS"
"Bdellovibrio bacteriovorus str. Tiberius"















