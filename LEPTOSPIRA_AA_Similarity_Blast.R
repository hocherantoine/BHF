######################################################
## Project: BACTERIAL HISTONES
## Compare a.a similarities within Leptospira species
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
setwd("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/LEPTOSPIRA/RECIPROCALBLAST/")

#Libraries : 
library(rtracklayer)
library(metablastr)

#List of genbank genomes for leptospira
ProteomeDir="/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/LEPTOSPIRA/genome_assemblies_genome_fasta/genome_assemblies_prot_fasta/ncbi-genomes-2022-06-17/"
# #For each proteome, compute the best reciprocal hits : 
setwd(ProteomeDir)

Genomes=read.delim("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/DB_BACT95_HCLUST0.5/PROTEOMES_FASTA/AnnotatedGenomes/Annotated_genomes_info.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

Prot_available=list.files(pattern = ".faa$")

for (Organism in Prot_available){
  
  Assembly=paste(strsplit(Organism,"_")[[1]][1:2],collapse = "_")
  Specie=Genomes[which(Genomes$Assembly==Assembly),]$Specie
  
  
  for(Organism2 in Prot_available){
    Assembly2=paste(strsplit(Organism2,"_")[[1]][1:2],collapse = "_")
    Specie2=Genomes[which(Genomes$Assembly==Assembly2),]$Specie
    
    
    if(Organism!=Organism2){
      if(file.exists(paste("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/LEPTOSPIRA/RECIPROCALBLAST/",Assembly,"_vs_",Assembly2,"__",Specie,"_vs_",Specie2,"_rblast_eval_1e3.txt",sep=""))==F){
        
        
        Rhits=blast_best_reciprocal_hit(query=Organism,subject = Organism2,evalue = 0.001 ,search_type = "protein_to_protein",task = "blastp",is.subject.db = F,cores = 4)
        
        write.table(Rhits,paste("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/LEPTOSPIRA/RECIPROCALBLAST/",Assembly,"_vs_",Assembly2,"__",Specie,"_vs_",Specie2,"_rblast_eval_1e3.txt",sep=""),row.names = F,sep="\t",quote=F)
      }
    }
  }
}






#2nd step : %identity in Bdellovibrio bacteriovorus HD100 vs others : 


setwd("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/LEPTOSPIRA/RECIPROCALBLAST/")
FilesBRB=list.files(pattern="GCA_000092565.1_vs*")

BRB=read.table(FilesBRB[1],header=T,sep="\t",quote="",stringsAsFactors = F)
BRB=BRB[,1:3]
names(BRB)[2:3]=paste( names(BRB)[2:3],paste(strsplit(FilesBRB[1],"_")[[1]][4:5],collapse="_"),sep="_")

for(i in FilesBRB[2:length(FilesBRB)]){
  BRB.2=read.table(i,header=T,sep="\t",quote="",stringsAsFactors = F)
  BRB.2=BRB.2[,which(names(BRB.2)%in%c("query_id","subject_id","perc_identity"))]
  names(BRB.2)[2:3]=paste( names(BRB.2)[2:3],paste(strsplit(i,"_")[[1]][4:5],collapse="_"),sep="_")
  BRB=merge(BRB,BRB.2,by="query_id",all.x=T,all.y=T)
}


write.table(BRB,file="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/AADIVERGENCE/LeptospiraInterrogansLai_GCA_000092565.1_vs_others_BRBH.txt",row.names=F,quote=F,sep="\t")


#Reload : 

BRB=read.delim("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/AADIVERGENCE/LeptospiraInterrogansLai_GCA_000092565.1_vs_others_BRBH.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

#Compute average IDentity , and see where our gene of interest lies within the distribution of all identities : 

BRB$AvgID=apply(BRB[,grep(names(BRB),pattern="perc_identity_")],1,function(x) mean(x,na.rm=T))

BRB$HF=ifelse(BRB$query_id=="AAN49657.2",1,0)


Values=BRB[which(BRB$HF==1),grep(names(BRB),pattern="perc_identity_")]

library(ggplot2);library(ggpubr)
Example=ggplot(data=BRB)+geom_violin(aes(x=1,y=AvgID))+geom_point(aes(x=1,y=BRB[which(BRB$HF==1),]$AvgID))+ylab("Average aa identity")+theme_pubclean()+theme(aspect.ratio=3)+ylim(c(0,100))

ggsave(Example,filename="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/AADIVERGENCE/LeptospiraInterrogans_vs_others_BRBH.pdf")


#Trying cumulative frequency plot : 

#Leptospira mayottensis
BRBorder=BRB[order(BRB$perc_identity_GCA_000306675.3),]
BRBorder=BRBorder[which(is.na(BRBorder$perc_identity_GCA_000306675.3)==F),]

#Leptospira hartskeerlii
BRBorder2=BRB[order(BRB$perc_identity_GCA_002811475.1),]
BRBorder2=BRBorder2[which(is.na(BRBorder2$perc_identity_GCA_002811475.1)==F),]

#Leptospira barantonii
BRBorder3=BRB[order(BRB$perc_identity_GCA_002811925.1),]
BRBorder3=BRBorder3[which(is.na(BRBorder3$perc_identity_GCA_002811925.1)==F),]

BRBExport=BRB[,which(names(BRB)%in%c("query_id","HF","subject_id_GCA_000306675.3","perc_identity_GCA_000306675.3","subject_id_GCA_002811475.1","perc_identity_GCA_002811475.1","subject_id_GCA_002811925.1","perc_identity_GCA_002811925.1"))]

ThreeSpecieExample=ggplot() + geom_step(data=BRBorder3,aes(x=perc_identity_GCA_002811925.1,y=..y..),stat="ecdf",col=col_vector[3],size=1)+ geom_step(data=BRBorder2,aes(x=perc_identity_GCA_002811475.1,y=..y..),stat="ecdf",col=col_vector[2],size=1)+ geom_step(data=BRBorder,aes(x=perc_identity_GCA_000306675.3,y=..y..),stat="ecdf",col=col_vector[1],size=1)+geom_point(data=BRBorder[which(BRBorder$HF==1),],aes(y=which(BRBorder$HF==1)/dim(BRBorder)[1],x=perc_identity_GCA_000306675.3),col=col_vector[1],size=2)+geom_point(data=BRBorder2[which(BRBorder2$HF==1),],aes(y=which(BRBorder2$HF==1)/dim(BRBorder2)[1],x=perc_identity_GCA_002811475.1),col=col_vector[2],size=2)+geom_point(data=BRBorder3[which(BRBorder3$HF==1),],aes(y=which(BRBorder3$HF==1)/dim(BRBorder3)[1],x=perc_identity_GCA_002811925.1),col=col_vector[3],size=2)+theme_pubr()+theme(aspect.ratio = 1/3)+xlab(" best reciprocal blast hits a.a identity (%)")+ylab("Fraction")+ scale_y_continuous(labels = function(x) paste0(x*100, "%"))
ggsave(ThreeSpecieExample,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/AADIVERGENCE/LeptospiraInterrogans_vs_3specie_Focus_HF.pdf")




#Exporting source data:
BRBExport=BRB[,which(names(BRB)%in%c("query_id","HF","subject_id_GCA_000306675.3","perc_identity_GCA_000306675.3","subject_id_GCA_002811475.1","perc_identity_GCA_002811475.1","subject_id_GCA_002811925.1","perc_identity_GCA_002811925.1"))]
write.table(BRBExport,"/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/AADIVERGENCE/Leptospira_Data_source_BRB_Figure.txt",row.names=F,sep="\t",quote=F)












