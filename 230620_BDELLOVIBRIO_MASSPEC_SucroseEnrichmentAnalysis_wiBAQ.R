#####################################################
## Project: BACTERIAL HISTONES
## MASS SPEC BDELLOVIBRIO ANALYSIS
## 
## 
## Aim : 
## Analyse Nucleoid enrichment experiments 
## 
## 
## Date: January 2021
## Amended 20 june 2023
## Author: Antoine Hocher
####################################################

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
source("/Users/ahocher/Dropbox/Scripts/PLOT_LIBRARY/AH_plot_library.R")
# 

library(DEP)
library("SummarizedExperiment") 
#Largely inspired from bellow for the dep analysis :
#https://www.csc.fi/documents/200270/382475/Proteomics_R_script/dfe2d119-1cc2-4228-8032-767e2bb24f5a



if(file.exists("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/BDELLOVIBRIO/NucEnrich/b023p051_AH_Bd1-3_Nuc-Top_proteinGroupsPFAMANNO_wiBAQ.txt")==F){
#Loading mass spectrometry data :
data=read.delim("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/Original_data_wiBAQ/b023p051_Bdello_with_iBAQ_proteinGroups.txt",header=T,sep="\t",stringsAsFactors = F)


#Removing e.coli proteins (not exactly contaminants as Bdellovibrio prey e.coli: 
if(length(grep("ECOBD",data$Fasta.headers))>0){
data=data[-grep("ECOBD",data$Fasta.headers),]}

#Because there has been some name formating issue
data$ID=unlist(lapply(data$Protein.IDs,function(x) paste(strsplit(x,"\\.1")[[1]][1],".1",collapse="",sep="")))


#To annotate PFAM dna binding domains : 
#I load a result file from a previous hmm search
library(rhmmer)
Bact_Res=as.data.frame(read_tblout("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/HMM_RESULTS/EBMC_Bact-vs-Pfam.txt"))


#Here formating and subselecting bdellovibrio data
Bact_Res$Taxid=unlist(lapply(Bact_Res$domain_name, function(x) strsplit(x,"_")[[1]][1]))
Bact_ResSub=Bact_Res[which(Bact_Res$Taxid=="264462"),]
rm(Bact_Res)

Bact_ResSub$ID=unlist(lapply(Bact_ResSub$domain_name, function(x) strsplit(x,"@")[[1]][2]))
Bact_ResSub$ID=gsub("\\|","",Bact_ResSub$ID)


#Hmm corresponding to PFAM DNA binding , or to NAPs :
PfamDNAbinding=read.table(file="/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/DBPome/PFAM_DNA_binding_list.txt",stringsAsFactors = F)$V1

PfamDNAbinding=c(PfamDNAbinding,"TFIIE_alpha")

NAPsHMM=c("Cren7","7kD_DNA_binding","CBFD_NFYB_HMF","Alba","Bac_DNA_binding","MC1","Histone")


#For each protein annoatting pfam hit
for( PFAM in unique(Bact_ResSub$query_name)){
  if(length(which(data$ID%in%Bact_ResSub[which(Bact_ResSub$query_name==PFAM),]$ID))>0){
    print(PFAM)
    data[,PFAM]=0
    data[which(data$ID%in%Bact_ResSub[which(Bact_ResSub$query_name==PFAM),]$ID),PFAM]=1}
}



#Annotating all proteins containing a DNA binding domain : 
data$PFAM_DNA_binding=0
data[which(data$ID%in%Bact_ResSub[which(Bact_ResSub$query_name%in%PfamDNAbinding),]$ID),]$PFAM_DNA_binding=1

#Annotating all NAPs
data$NAP=0
data[which(data$ID%in%Bact_ResSub[which(Bact_ResSub$query_name%in%NAPsHMM),]$ID),]$NAP=1


#Starting to merge technical duplicates

LFQ_columns <- grep("LFQ.", colnames(data)) # get LFQ column numbers

#Verify names : 
#colnames(data)[grep("LFQ.", colnames(data))]

data$LFQ.intensity.Bd1_Nuc=(data$LFQ.intensity.Bd1_Nuc_TR01+data$LFQ.intensity.Bd1_Nuc_TR02)/2
data$LFQ.intensity.Bd2_Nuc=(data$LFQ.intensity.Bd2_Nuc_TR01+data$LFQ.intensity.Bd2_Nuc_TR02)/2
data$LFQ.intensity.Bd3_Nuc=(data$LFQ.intensity.Bd3_Nuc_TR01+data$LFQ.intensity.Bd3_Nuc_TR02)/2

data$LFQ.intensity.Bd1_Top=(data$LFQ.intensity.Bd1_Top_TR01+data$LFQ.intensity.Bd1_Top_TR02)/2
data$LFQ.intensity.Bd2_Top=(data$LFQ.intensity.Bd2_Top_TR01+data$LFQ.intensity.Bd2_Top_TR02)/2
data$LFQ.intensity.Bd3_Top=(data$LFQ.intensity.Bd3_Top_TR01+data$LFQ.intensity.Bd3_Top_TR02)/2


data=data[,-LFQ_columns]

#Export
write.table(data,file="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/BDELLOVIBRIO/NucEnrich/b023p051_AH_Bd1-3_Nuc-Top_proteinGroupsPFAMANNO_wiBAQ.txt",sep="\t",quote=F,row.names=F)

}


#Re-load
data=read.delim("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/BDELLOVIBRIO/NucEnrich/b023p051_AH_Bd1-3_Nuc-Top_proteinGroupsPFAMANNO_wiBAQ.txt",header=T,sep="\t",quote="",stringsAsFactors = F)


LFQ_columns <- grep("LFQ.", colnames(data)) # get LFQ column numbers

data=data[,c(LFQ_columns,which(names(data)%in%c("PFAM_DNA_binding","ID","Protein.IDs","NAP","Reverse","Potential.contaminant")))]


#This is required for dep
data$Gene.names=data$ID




#Generate a SummarizedExperiment object by parsing condition information from the
experimental_design=as.data.frame(matrix(nrow=6,ncol=3))
names(experimental_design)=c("label","condition","replicate")
experimental_design$label=c("Bd1_Nuc","Bd1_Top","Bd2_Nuc","Bd2_Top","Bd3_Nuc","Bd3_Top")
experimental_design$condition=c("Nuc","Top","Nuc","Top","Nuc","Top")
experimental_design$replicate=c(1,1,2,2,3,3)


data_results <- LFQ(data, experimental_design,fun ="man",control = "Top" ,type = "control", alpha = 0.05, lfc = 0)

#Storing the results
AA=data_results$results

#Correcting p-value
AA$PadjBH=p.adjust(AA$Nuc_vs_Top_p.val, method="BH")


#Quick look at Bd0055
AA[which(AA$name=="CAE77736.1"),]

#Rename the ID column to avoid problem merging later on
names(AA)[which(names(AA)=="ID")]="FullEntry"

#Loading whole cell extract : 
WCEMspec=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/Original_data_wiBAQ/b023p050_Bdello_with_iBAQ_proteinGroups_SelectedColumns.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

#Some proteins have a double match to E.coli and B.bacteriovorus (example recA); because the condition we study barely has any E.coli left, we put forward the B.bacteriovorus annotation
WCEMspec$ID=WCEMspec$Protein.IDs
WCEMspec$ID=unlist(lapply(WCEMspec$ID,function(x) strsplit(x,";")[[1]][1]))

iBAQ_columns <- colnames(WCEMspec)[grep("iBAQ.", colnames(WCEMspec))] # get iBAQ column names



WCEMspec=WCEMspec[,which(names(WCEMspec)%in%c("ID",iBAQ_columns))]

#Merging with the iBAQ from nucleoid experiments : 

#Re-load the table
data=read.delim("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/BDELLOVIBRIO/NucEnrich/b023p051_AH_Bd1-3_Nuc-Top_proteinGroupsPFAMANNO_wiBAQ.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

iBAQ_columns <- colnames(data)[grep("iBAQ.", colnames(data))] # get iBAQ column names
data=data[,which(names(data)%in%c("ID",iBAQ_columns,"PFAM_DNA_binding"))]
#Averaging all replicates for Nuc and Top fractions
data$iBAQ.Bd1_Nuc_avg=(data$iBAQ.Bd1_Nuc_TR01+data$iBAQ.Bd1_Nuc_TR02+data$iBAQ.Bd2_Nuc_TR01+data$iBAQ.Bd2_Nuc_TR02+data$iBAQ.Bd3_Nuc_TR01+data$iBAQ.Bd3_Nuc_TR02)/6
data$iBAQ.Bd1_Top_avg=(data$iBAQ.Bd1_Top_TR01+data$iBAQ.Bd1_Top_TR02+data$iBAQ.Bd2_Top_TR01+data$iBAQ.Bd2_Top_TR02+data$iBAQ.Bd3_Top_TR01+data$iBAQ.Bd3_Top_TR02)/6


WCEMspec2=merge(WCEMspec,data,by="ID",all=T)
#Merging with the nucleoid enrichment experiments : 
WCEMspec3=merge(WCEMspec2,AA,by.x="ID",by.y="name",all=T)


#Just for the plot i put back the PFAM_DNA_binding annotation
AA$PFAM_DNA_binding=0
AA[which(AA$name%in%data[which(data$PFAM_DNA_binding==1),]$ID),]$PFAM_DNA_binding=1

#Seeing if we have a statistical enrichment for dna binding proteins

#AH 13.11.2022 : exporting boxplot figure : 
pdf("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/230620_BdellovibrioBActeriovorus_Nuc_vs_Top_Ratio_vs_PfamDNAbinding_Boxplot.pdf")
boxplot(AA$Nuc_vs_Top_ratio~AA$PFAM_DNA_binding)
dev.off()

t.test(AA$Nuc_vs_Top_ratio~AA$PFAM_DNA_binding)
wilcox.test(AA$Nuc_vs_Top_ratio~AA$PFAM_DNA_binding)

volcanoPlot=ggplot(data=AA)+geom_point(aes(x=Nuc_vs_Top_ratio,y=-log10(Nuc_vs_Top_p.val),color=as.factor(PFAM_DNA_binding)),size=0.5)+geom_point(data=AA[which(AA$name=="CAE77736.1"),],aes(x=Nuc_vs_Top_ratio,y=-log10(Nuc_vs_Top_p.val)),color="orange",size=2)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))
ggsave(plot = volcanoPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/230620_Bacterial_histones_phylogeny/Figures/MASS_SPEC/230620_Nuc_vs_Top_volcanoPlotUncorrectedpval.pdf")


volcanoPlot=ggplot(data=AA)+geom_point(aes(x=Nuc_vs_Top_ratio,y=-log10(PadjBH),color=as.factor(PFAM_DNA_binding)),size=0.5)+geom_point(data=AA[which(AA$name=="CAE77736.1"),],aes(x=Nuc_vs_Top_ratio,y=-log10(Nuc_vs_Top_p.val)),color="orange",size=2)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))
ggsave(plot = volcanoPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/230620_Bacterial_histones_phylogeny/Figures/MASS_SPEC/230620_Nuc_vs_Top_volcanoPlot.pdf")




EnrichmentPlot=ggplot(data=WCEMspec3)+geom_point(aes(x=Nuc_vs_Top_ratio,y=sqrt(LFQ.intensity.Bd1_Nuc),color=as.factor(PFAM_DNA_binding)),size=0.5)+geom_point(data=WCEMspec3[which(WCEMspec3$ID=="CAE77736.1"),],aes(x=Nuc_vs_Top_ratio,y=sqrt(LFQ.intensity.Bd1_Nuc)),color="orange",size=2)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))

ggsave(plot = EnrichmentPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/230620_Nuc_vs_Top_enrichment_orange_candidate_blueDNAbinding.pdf")






library(corrplot)
LFQ_columns <- grep("LFQ.", colnames(WCEMspec3)) # get LFQ column numbers

CorMat=WCEMspec3[,LFQ_columns]
CorMatRes=cor(CorMat,use="complete.obs")
pdf("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/230620_PearsonCorrealtion_MS_replicates.pdf",height=6,width=6)
corrplot(CorMatRes)
dev.off()

CorMatRes=cor(CorMat,method="spearman",use="complete.obs")

pdf("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/230620_SpearmanCorrealtion_MS_replicates.pdf",height=6,width=6)
corrplot(CorMatRes)
dev.off()


#Focusing on the histone fold protein : 

WCEMspec3$CBFD_NFYB_HMF=ifelse(WCEMspec3$ID=="CAE77736.1",1,0)
Example=ggplot(data=WCEMspec3)+geom_point(aes(x=LFQ.intensity.Bd1_Top,y=LFQ.intensity.Bd1_Nuc,color=CBFD_NFYB_HMF,size=as.factor(CBFD_NFYB_HMF)))+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+xlim(c(0,3e9))+ylim(c(0,3e9))+ggtitle("axis scale meant \nnot displaying one protein (intensity ~ 1e10")

ggsave(plot = Example,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/230620_Nuc_vs_Top_replicates1_linearscale.pdf")


Rep1=ggplot(data=WCEMspec3)+geom_point(aes(x=LFQ.intensity.Bd1_Top,y=LFQ.intensity.Bd1_Nuc,color=CBFD_NFYB_HMF,size=as.factor(CBFD_NFYB_HMF)))+scale_x_log10()+scale_y_log10()+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")
Rep2=ggplot(data=WCEMspec3)+geom_point(aes(x=LFQ.intensity.Bd2_Top,y=LFQ.intensity.Bd2_Nuc,color=CBFD_NFYB_HMF,size=as.factor(CBFD_NFYB_HMF)))+scale_x_log10()+scale_y_log10()+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")
Rep3=ggplot(data=WCEMspec3)+geom_point(aes(x=LFQ.intensity.Bd3_Top,y=LFQ.intensity.Bd3_Nuc,color=CBFD_NFYB_HMF,size=as.factor(CBFD_NFYB_HMF)))+scale_x_log10()+scale_y_log10()+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")


ggsave(plot = ggarrange(Rep1,Rep2,Rep3,ncol=3,nrow=1),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/230620_Nuc_vs_Top_replicates1_to3.pdf",height = 8,width=14)



#Localisation : 

LocalisationResult=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/EXTERNAL_DATASETS/PSORTDB_LOCALISATION_BDELLOVIBRIO/psortdb-results_Genbank.txt",header=T,sep="\t",stringsAsFactors = F,quote="")


LocalisationResult=LocalisationResult[,which(names(LocalisationResult)%in%c("ID","Final_Localization"))]
WCEMspec4=merge(WCEMspec3,LocalisationResult,by="ID",all.x=T)






#Exporting the table
WCEMspec4=WCEMspec4[order(-WCEMspec4$iBAQ.Bd1_Nuc_avg),]
WCEMspec4$name=WCEMspec4$Gene.names

write.table(WCEMspec4,"/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/BDELLOVIBRIO/NucEnrich/230620_b023p051_Bdellovibrio_Merged_WCE_and_SucroseProteomics_w_localisation_wiBAQ.txt",sep="\t",row.names=F)





