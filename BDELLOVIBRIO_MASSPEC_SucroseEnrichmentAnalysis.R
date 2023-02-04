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
## Author: Antoine Hocher
####################################################

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
source("/Users/ahocher/Dropbox/Scripts/PLOT_LIBRARY/AH_plot_library.R")
# 

library(DEP)
library("SummarizedExperiment") 
#Largely inspired from bellow for the dep analysis :
#https://www.csc.fi/documents/200270/382475/Proteomics_R_script/dfe2d119-1cc2-4228-8032-767e2bb24f5a



if(file.exists("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/BDELLOVIBRIO/NucEnrich/b023p051_AH_Bd1-3_Nuc-Top_proteinGroupsPFAMANNO.txt")==F){
#Loading mass spectrometry data :
data=read.delim("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/BDELLOVIBRIO/NucEnrich/b023p051_AH_Bd1-3_Nuc-Top_proteinGroups.txt",header=T,sep="\t",stringsAsFactors = F)


#Removing e.coli proteins (not exactly contaminants as Bdellovibrio prey e.coli: 
data=data[-grep("ECOBD",data$Fasta.headers),]

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

data$LFQ.intensity.Bd1_Nuc=(data$LFQ.intensity.Bd1_Nuc_01+data$LFQ.intensity.Bd1_Nuc_02)/2
data$LFQ.intensity.Bd2_Nuc=(data$LFQ.intensity.Bd2_Nuc_01+data$LFQ.intensity.Bd2_Nuc_02)/2
data$LFQ.intensity.Bd3_Nuc=(data$LFQ.intensity.Bd3_Nuc_01+data$LFQ.intensity.Bd3_Nuc_02)/2

data$LFQ.intensity.Bd1_Top=(data$LFQ.intensity.Bd1_Top_01+data$LFQ.intensity.Bd1_Top_02)/2
data$LFQ.intensity.Bd2_Top=(data$LFQ.intensity.Bd2_Top_01+data$LFQ.intensity.Bd2_Top_02)/2
data$LFQ.intensity.Bd3_Top=(data$LFQ.intensity.Bd3_Top_01+data$LFQ.intensity.Bd3_Top_02)/2


data=data[,-LFQ_columns]

#Export
write.table(data,file="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/BDELLOVIBRIO/NucEnrich/b023p051_AH_Bd1-3_Nuc-Top_proteinGroupsPFAMANNO.txt",sep="\t",quote=F,row.names=F)

}


#Re-load
data=read.delim("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/BDELLOVIBRIO/NucEnrich/b023p051_AH_Bd1-3_Nuc-Top_proteinGroupsPFAMANNO.txt",header=T,sep="\t",quote="",stringsAsFactors = F)


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
WCEMspec=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/BDELLOVIBRIO/b023p050_B_bacteriovorus_WCE_proteinGroups_PFAMANNO.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

#Some proteins have a double match to E.coli and B.bacteriovorus (example recA); because the condition we study barely has any E.coli left, we put forward the B.bacteriovorus annotation

WCEMspec$ID=unlist(lapply(WCEMspec$ID,function(x) strsplit(x,";")[[1]][1]))

LFQ_columns <- grep("LFQ.", colnames(WCEMspec)) # get LFQ column numbers


WCEMspec$LFQ.intensity.WCE1=(WCEMspec$LFQ.intensity.B.bacteriovorus_WCE1_1+WCEMspec$LFQ.intensity.B.bacteriovorus_WCE1_2)/2
WCEMspec$LFQ.intensity.WCE2=(WCEMspec$LFQ.intensity.B.bacteriovorus_WCE2_1+WCEMspec$LFQ.intensity.B.bacteriovorus_WCE2_2)/2

WCEMspec=WCEMspec[,-LFQ_columns]

WCEMspec=WCEMspec[,which(names(WCEMspec)%in%c("ID","LFQ.intensity.WCE1","LFQ.intensity.WCE2"))]

#Merging with the LFQ from nucleoid experiments : 
WCEMspec2=merge(WCEMspec,data,by="ID",all=T)
#Merging with the nucleoid enrichment experiments : 
WCEMspec3=merge(WCEMspec2,AA,by.x="ID",by.y="name",all=T)


#Just for the plot i put back the PFAM_DNA_binding annotation
AA$PFAM_DNA_binding=0
AA[which(AA$name%in%data[which(data$PFAM_DNA_binding==1),]$ID),]$PFAM_DNA_binding=1

#Seeing if we have a statistical enrichment for dna binding proteins

#AH 13.11.2022 : exporting boxplot figure : 
pdf("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/BdellovibrioBActeriovorus_Nuc_vs_Top_Ratio_vs_PfamDNAbinding_Boxplot.pdf")
boxplot(AA$Nuc_vs_Top_ratio~AA$PFAM_DNA_binding)
dev.off()

t.test(AA$Nuc_vs_Top_ratio~AA$PFAM_DNA_binding)
wilcox.test(AA$Nuc_vs_Top_ratio~AA$PFAM_DNA_binding)

volcanoPlot=ggplot(data=AA)+geom_point(aes(x=Nuc_vs_Top_ratio,y=-log10(Nuc_vs_Top_p.val),color=as.factor(PFAM_DNA_binding)),size=0.5)+geom_point(data=AA[which(AA$name=="CAE77736.1"),],aes(x=Nuc_vs_Top_ratio,y=-log10(Nuc_vs_Top_p.val)),color="orange",size=2)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))
ggsave(plot = volcanoPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Nuc_vs_Top_volcanoPlotUncorrectedpval.pdf")


volcanoPlot=ggplot(data=AA)+geom_point(aes(x=Nuc_vs_Top_ratio,y=-log10(PadjBH),color=as.factor(PFAM_DNA_binding)),size=0.5)+geom_point(data=AA[which(AA$name=="CAE77736.1"),],aes(x=Nuc_vs_Top_ratio,y=-log10(Nuc_vs_Top_p.val)),color="orange",size=2)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))
ggsave(plot = volcanoPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Nuc_vs_Top_volcanoPlot.pdf")

#RecA : 
ggplot(data=AA)+geom_point(aes(x=Nuc_vs_Top_ratio,y=-log10(PadjBH),color=as.factor(PFAM_DNA_binding)),size=0.5)+geom_point(data=AA[which(AA$name=="CAE78490.1"),],aes(x=Nuc_vs_Top_ratio,y=-log10(Nuc_vs_Top_p.val)),color="orange",size=2)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))



EnrichmentPlot=ggplot(data=WCEMspec3)+geom_point(aes(x=Nuc_vs_Top_ratio,y=sqrt(LFQ.intensity.Bd1_Nuc),color=as.factor(PFAM_DNA_binding)),size=0.5)+geom_point(data=WCEMspec3[which(WCEMspec3$ID=="CAE77736.1"),],aes(x=Nuc_vs_Top_ratio,y=sqrt(LFQ.intensity.Bd1_Nuc)),color="orange",size=2)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))

ggsave(plot = EnrichmentPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Nuc_vs_Top_enrichment_orange_candidate_blueDNAbinding.pdf")






library(corrplot)
LFQ_columns <- grep("LFQ.", colnames(WCEMspec3)) # get LFQ column numbers

CorMat=WCEMspec3[,LFQ_columns]
CorMatRes=cor(CorMat,use="complete.obs")
pdf("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/PearsonCorrealtion_MS_replicates.pdf",height=6,width=6)
corrplot(CorMatRes)
dev.off()

CorMatRes=cor(CorMat,method="spearman",use="complete.obs")

pdf("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/SpearmanCorrealtion_MS_replicates.pdf",height=6,width=6)
corrplot(CorMatRes)
dev.off()


#Focusing on the histone fold protein : 

WCEMspec3$CBFD_NFYB_HMF=ifelse(WCEMspec3$ID=="CAE77736.1",1,0)
Example=ggplot(data=WCEMspec3)+geom_point(aes(x=LFQ.intensity.Bd1_Top,y=LFQ.intensity.Bd1_Nuc,color=CBFD_NFYB_HMF,size=as.factor(CBFD_NFYB_HMF)))+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+xlim(c(0,3e9))+ylim(c(0,3e9))+ggtitle("axis scale meant \nnot displaying one protein (intensity ~ 1e10")

ggsave(plot = Example,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Nuc_vs_Top_replicates1_linearscale.pdf")


Rep1=ggplot(data=WCEMspec3)+geom_point(aes(x=LFQ.intensity.Bd1_Top,y=LFQ.intensity.Bd1_Nuc,color=CBFD_NFYB_HMF,size=as.factor(CBFD_NFYB_HMF)))+scale_x_log10()+scale_y_log10()+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")
Rep2=ggplot(data=WCEMspec3)+geom_point(aes(x=LFQ.intensity.Bd2_Top,y=LFQ.intensity.Bd2_Nuc,color=CBFD_NFYB_HMF,size=as.factor(CBFD_NFYB_HMF)))+scale_x_log10()+scale_y_log10()+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")
Rep3=ggplot(data=WCEMspec3)+geom_point(aes(x=LFQ.intensity.Bd3_Top,y=LFQ.intensity.Bd3_Nuc,color=CBFD_NFYB_HMF,size=as.factor(CBFD_NFYB_HMF)))+scale_x_log10()+scale_y_log10()+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")


ggsave(plot = ggarrange(Rep1,Rep2,Rep3,ncol=3,nrow=1),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Nuc_vs_Top_replicates1_to3.pdf",height = 8,width=14)



#Localisation : 

LocalisationResult=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/EXTERNAL_DATASETS/PSORTDB_LOCALISATION_BDELLOVIBRIO/psortdb-results_Genbank.txt",header=T,sep="\t",stringsAsFactors = F,quote="")


LocalisationResult=LocalisationResult[,which(names(LocalisationResult)%in%c("ID","Final_Localization"))]
WCEMspec4=merge(WCEMspec3,LocalisationResult,by="ID",all.x=T)


WCEMspec4$LFQ.intensity.Bd_Nuc_avg=(WCEMspec4$LFQ.intensity.Bd1_Nuc+WCEMspec4$LFQ.intensity.Bd2_Nuc+WCEMspec4$LFQ.intensity.Bd3_Nuc)/3




#Exporting the table
WCEMspec4=WCEMspec4[order(-WCEMspec4$LFQ.intensity.Bd_Nuc_avg),]
WCEMspec4$name=WCEMspec4$Gene.names

write.table(WCEMspec4,"/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/BDELLOVIBRIO/NucEnrich/b023p051_Bdellovibrio_Merged_WCE_and_SucroseProteomics_w_localisation.txt",sep="\t",row.names=F)




ggplot(data=WCEMspec4)+geom_point(aes(y=Nuc_vs_Top_ratio,x=sqrt(LFQ.intensity.Bd_Nuc_avg),color=as.factor(Final_Localization)),size=0.5)+geom_point(data=WCEMspec4[which(WCEMspec4$name=="CAE77736.1"),],aes(y=Nuc_vs_Top_ratio,x=sqrt(LFQ.intensity.Bd_Nuc_avg)),color="black",size=2)+theme_pubclean()+theme(aspect.ratio=1/4,legend.position = "")+scale_color_manual(values=col_vector)


LocaAndName=ggplot(data=WCEMspec4)+geom_point(aes(y=Nuc_vs_Top_ratio,x=log2(LFQ.intensity.Bd_Nuc_avg),color=as.factor(Final_Localization)),size=0.5)+geom_text_repel(data=WCEMspec4[which(WCEMspec4$name=="CAE77736.1"),],aes(y=Nuc_vs_Top_ratio,x=log2(LFQ.intensity.Bd_Nuc_avg),label=name),color="black",size=1,min.segment.length = 0.1)+theme_pubclean()+theme(aspect.ratio=1/1.61,legend.position = "")+scale_color_manual(values=col_vector)+xlab("Intensity, nucleoid fraction (log2)")+ylab("Enrichment \nNucleoid vs Top fraction (log2)")
ggsave(addSmallLegend(LocaAndName),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Nuc_vs_Top_enrichment_scatter_log2scale.pdf")


LocaAndName=ggplot(data=WCEMspec4)+geom_point(aes(y=Nuc_vs_Top_ratio,x=LFQ.intensity.Bd_Nuc_avg,color=as.factor(Final_Localization)),size=0.5)+geom_text_repel(data=WCEMspec4[which(WCEMspec4$name=="CAE77736.1"),],aes(y=Nuc_vs_Top_ratio,x=LFQ.intensity.Bd_Nuc_avg,label=name),color="black",size=1,min.segment.length = 0.1)+theme_pubclean()+theme(aspect.ratio=1/1.61,legend.position = "top")+scale_color_manual(values=col_vector)+xlab("Intensity, nucleoid fraction")+ylab("Enrichment \nNucleoid vs Top fraction (log2)")
ggsave(addSmallLegend(LocaAndName),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Nuc_vs_Top_enrichment_scatter_linearscale.pdf")


WCEMspec4=WCEMspec4[order(-WCEMspec4$LFQ.intensity.Bd_Nuc_avg),]
WCEMspec4$Description=gsub("\\[BdellovibriobacteriovorusHD100\\]","",WCEMspec4$ID)
DNABindingFocus=ggplot(data=WCEMspec4)+geom_point(aes(y=Nuc_vs_Top_ratio,x=log2(LFQ.intensity.Bd_Nuc_avg),color=as.factor(PFAM_DNA_binding)),size=0.5)+geom_text_repel(data=WCEMspec4[which(WCEMspec4$PFAM_DNA_binding==1 & WCEMspec4$Nuc_vs_Top_ratio>1 & WCEMspec4$PadjBH <0.01)[1:20],],aes(y=Nuc_vs_Top_ratio,x=log2(LFQ.intensity.Bd_Nuc_avg),label=Description),color="black",size=1,min.segment.length = 0.1)+theme_pubclean()+theme(aspect.ratio=1/1.61,legend.position = "top")+scale_color_manual(values=c("grey80","blue"))+xlab("Intensity, nucleoid fraction")+ylab("Enrichment \nNucleoid vs Top fraction (log2)")
ggsave(addSmallLegend(DNABindingFocus),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Nuc_vs_Top_enrichment_DNA_bindingFocus_scatter_log2scale.pdf")


DNABindingFocus=ggplot(data=WCEMspec4)+geom_point(aes(y=Nuc_vs_Top_ratio,x=LFQ.intensity.Bd_Nuc_avg,color=as.factor(PFAM_DNA_binding)),size=0.5)+geom_text_repel(data=WCEMspec4[which(WCEMspec4$PFAM_DNA_binding==1 & WCEMspec4$Nuc_vs_Top_ratio>1 & WCEMspec4$PadjBH <0.01)[1:20],],aes(y=Nuc_vs_Top_ratio,x=LFQ.intensity.Bd_Nuc_avg,label=Description),color="black",size=1,min.segment.length = 0.1)+theme_pubclean()+theme(aspect.ratio=1/1.61,legend.position = "top")+scale_color_manual(values=c("grey80","blue"))+xlab("Intensity, nucleoid fraction")+ylab("Enrichment \nNucleoid vs Top fraction")+xlim(c(0,3e9))
ggsave(addSmallLegend(DNABindingFocus),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Nuc_vs_Top_enrichment_DNA_bindingFocus_scatter.pdf")



#Showing Bd3044: 
DNABindingFocus=ggplot(data=WCEMspec4)+geom_point(aes(y=Nuc_vs_Top_ratio,x=LFQ.intensity.Bd_Nuc_avg,color=as.factor(PFAM_DNA_binding)),size=0.5)+geom_text_repel(data=WCEMspec4[which(WCEMspec4$PFAM_DNA_binding==1 & WCEMspec4$Nuc_vs_Top_ratio>1 & WCEMspec4$PadjBH <0.01)[1:20],],aes(y=Nuc_vs_Top_ratio,x=LFQ.intensity.Bd_Nuc_avg,label=Description),color="black",size=1,min.segment.length = 0.1)+geom_point(data=WCEMspec4[which(WCEMspec4$name=="CAE80810.1"),],aes(y=Nuc_vs_Top_ratio,x=LFQ.intensity.Bd_Nuc_avg),color="red",size=0.5)+theme_pubclean()+theme(aspect.ratio=1/1.61,legend.position = "top")+scale_color_manual(values=c("grey80","blue"))+xlab("Intensity, nucleoid fraction")+ylab("Enrichment \nNucleoid vs Top fraction")+xlim(c(0,3e9))
ggsave(addSmallLegend(DNABindingFocus),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Nuc_vs_Top_enrichment_DNA_bindingFocus_scatter_ShowingBd3044_inred.pdf")

#Showing the different annotated naps as well as Bd3044 : 
WCEMspec4[which(WCEMspec4$NAP==1),]

NAPFocus=ggplot(data=WCEMspec4)+geom_point(aes(y=Nuc_vs_Top_ratio,x=LFQ.intensity.Bd_Nuc_avg,color=as.factor(NAP)),size=0.5)+geom_text_repel(data=WCEMspec4[which(WCEMspec4$NAP==1 |WCEMspec4$name=="CAE80810.1"),],aes(y=Nuc_vs_Top_ratio,x=LFQ.intensity.Bd_Nuc_avg,label=Description),color="black",size=2,min.segment.length = 0.1)+theme_pubclean()+theme(aspect.ratio=1/1.61,legend.position = "top")+scale_color_manual(values=c("grey80","blue"))+xlab("Intensity, nucleoid fraction")+ylab("Enrichment \nNucleoid vs Top fraction")+xlim(c(0,3e9))+ggtitle( "CAE77736.1 = Bd0055,CAE80810.1=Bd3044, CAE79947.1= HU-2" )
ggsave(addSmallLegend(NAPFocus),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Nuc_vs_Top_enrichment_vs_NucleoidAbundancy_NAP_and_Bd3044.pdf")












ggplot(data=WCEMspec4)+geom_point(aes(y=Nuc_vs_Top_ratio,x=LFQ.intensity.Bd_Nuc_avg,color=as.factor(Final_Localization)),size=0.5)+geom_text_repel(data=WCEMspec4[which(WCEMspec4$name=="CAE77736.1"),],aes(y=Nuc_vs_Top_ratio,x=LFQ.intensity.Bd_Nuc_avg,label=name),color="black",size=1,min.segment.length = 0.1)+theme_pubclean()+theme(aspect.ratio=1/4,legend.position = "top")+scale_color_manual(values=col_vector)+xlim(0,3e9)





WCEMspec4$Final_Localization=factor(WCEMspec4$Final_Localization,levels = rev(c("Outer Membrane","Cytoplasmic Membrane","Extracellular","Unknown","Periplasmic","Cytoplasmic")))
LocalisationPlot=ggplot(data=WCEMspec4)+geom_violin(aes(y=Nuc_vs_Top_ratio,x=as.factor(Final_Localization),group=as.factor(Final_Localization),fill=as.factor(Final_Localization)),size=0,width=1)+geom_point(data=WCEMspec4[which(WCEMspec4$name=="CAE77736.1"),],aes(y=Nuc_vs_Top_ratio,x=as.factor(Final_Localization),group=as.factor(Final_Localization)),color="black",size=2) +theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_fill_manual(values=col_vector)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10))+xlab("")+coord_flip()
ggsave(LocalisationPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Nuc_vs_Top_enrichment_stratified_by_localisation.pdf")




WCEMspec4$Final_Localization=factor(WCEMspec4$Final_Localization,levels = rev(c("Outer Membrane","Cytoplasmic Membrane","Extracellular","Unknown","Periplasmic","Cytoplasmic")))
LocalisationPlot=ggplot(data=WCEMspec4)+geom_violin(aes(y=Top_vs_WCE_ratio,x=as.factor(Final_Localization),group=as.factor(Final_Localization),fill=as.factor(Final_Localization)),size=0,width=1)+geom_point(data=WCEMspec4[which(WCEMspec4$name=="CAE77736.1"),],aes(y=Top_vs_WCE_ratio,x=as.factor(Final_Localization),group=as.factor(Final_Localization)),color="black",size=2) +theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_fill_manual(values=col_vector)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10))+xlab("")+coord_flip()
ggsave(LocalisationPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Top_vs_WCE_enrichment_stratified_by_localisation.pdf")

WCEMspec4=WCEMspec4[order(-WCEMspec4$LFQ.intensity.Bd2_Top),]

ggplot(data=WCEMspec4[1:100,])+geom_bar(aes(x=as.factor(Final_Localization),group=as.factor(Final_Localization)),size=0,width=1)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_fill_manual(values=col_vector)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10))+xlab("")+coord_flip()


WCEMspec4=WCEMspec4[order(-WCEMspec4$LFQ.intensity.Bd2_Top),]

ggplot(data=WCEMspec4[1:100,])+geom_bar(aes(x=as.factor(Final_Localization),group=as.factor(Final_Localization)),size=0,width=1)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_fill_manual(values=col_vector)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10))+xlab("")+coord_flip()





AA[which(AA$name=="CAE79605.1"),]
AA[which(AA$name=="CAE77874.1"),]





#Checking percentile of Bd0055  :

AA=AA[order(-AA$),]

which(AA$name=="CAE77736.1")



data$Avg_LFQ_Nuc=(data$LFQ.intensity.Bd1_Nuc+data$LFQ.intensity.Bd2_Nuc+data$LFQ.intensity.Bd3_Nuc)/3
data=data[order(-data$Avg_LFQ_Nuc),]

which(data$ID=="CAE77736.1")
length(which(data$Avg_LFQ_Nuc>0))
19/1997

