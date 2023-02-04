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
library(rhmmer)
library(DEP)
library("SummarizedExperiment") 
library(metablastr)
source("/Users/ahocher/Dropbox/Scripts/ComputeProtStats_from_Fasta.R")

Sys.setenv(PATH="/Users/ahocher/opt/miniconda3/bin:/Users/ahocher/opt/miniconda3/condabin:/opt/local/bin:/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin")

#Largely inspired from bellow for the dep analysis :
#https://www.csc.fi/documents/200270/382475/Proteomics_R_script/dfe2d119-1cc2-4228-8032-767e2bb24f5a

#Loading mass spectrometry data :
data=read.delim("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/LEPTOSPIRA/INHOUSE/OriginalData/b023p057_proteinGroups_L_interrogans_Uniprot.txt",header=T,sep="\t",stringsAsFactors = F)
names(data)[1]="ID"

data$Avg_WCE_LFQ=(data$LFQ.intensity.INT_wce_1+data$LFQ.intensity.INT_wce_2)/2
data=data[order(-data$Avg_WCE_LFQ),]
head(data)
which(data$ID=="A0A2H1XGH2")
length(data$Avg_WCE_LFQ[which(data$Avg_WCE_LFQ!=0)])

#Obtaining the corresponding names : uniprot => Genbank
#
#
#Uniprot sequences  : 
Uniprot="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/LEPTOSPIRA/INHOUSE/OriginalData/UP000234460_uniprot_Leptospira_interrogans_serovar_Manilae_20210908.fasta"
UniprotFa=read.fasta(Uniprot)
Names=unlist(lapply(getName(UniprotFa),function(x) strsplit(x,"\\|")[[1]][2]))
write.fasta(UniprotFa,Names,file.out = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/LEPTOSPIRA/INHOUSE/OriginalData/UP000234460_uniprot_Leptospira_interrogans_serovar_Manilae_20210908_ShortNames.fasta")

#Reload the proper names
Uniprot="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/LEPTOSPIRA/INHOUSE/OriginalData/UP000234460_uniprot_Leptospira_interrogans_serovar_Manilae_20210908_ShortNames.fasta"

#Genbank sequences :

Genbank="/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/LEPTOSPIRA/Proteome_L_interrogans_Manilae_HP/GCA_001047655.1_ASM104765v1_protein.faa"

Corresp=blast_best_reciprocal_hit(Uniprot,Genbank,search_type = "protein_to_protein",task="blastp-fast")
Correspdt=as.data.frame(Corresp)
Correspdt=Correspdt[order(Correspdt$evalue),]
head(Correspdt)


##
##
###VERY IMPORTANT, DO NOT SKIP. 
###
#Verify that the vast majority has a 100% identity ( otherwise it's not a correspondence table..)
table(Correspdt$perc_identity)

Correspdt=Correspdt[,c(1:2)]
names(Correspdt)=c("UniprotID","ID")

#Merge to have uniprot and genbank id in the same table
data$UniprotID=data$ID

data2=merge(data,Correspdt,by="UniprotID",all.x=T)



#To annotate PFAM dna binding domains : 
#Because it got complicated I just recomputed pfam domains for each proteomes (it's faster than looking in the larger bacterial database)
#

#L.interrogans


setwd("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/LEPTOSPIRA/INHOUSE/OriginalData")


#system("hmmsearch --cut_ga --noali --tblout UP000234460_uniprot_Leptospira_interrogans_serovar_Manilae_20210908_ShortNames-vs-Pfam.txt /Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/PfamDb/Pfam-A.hmm UP000234460_uniprot_Leptospira_interrogans_serovar_Manilae_20210908_ShortNames.fasta > UP000234460_uniprot_Leptospira_interrogans_serovar_Manilae_20210908_ShortNames-vs-Pfam_Full.txt")



#Loading Hmmsearch results : 
Bact_Res=as.data.frame(read_tblout("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/LEPTOSPIRA/INHOUSE/OriginalData/UP000234460_uniprot_Leptospira_interrogans_serovar_Manilae_20210908_ShortNames-vs-Pfam.txt"))
Bact_Res$ID=Bact_Res$domain_name

#Hmm corresponding to PFAM DNA binding , or to NAPs :
PfamDNAbinding=read.table(file="/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/DBPome/PFAM_DNA_binding_list_no_AAA.txt",stringsAsFactors = F)$V1

PfamDNAbinding=c(PfamDNAbinding,"TFIIE_alpha")

NAPsHMM=c("Cren7","7kD_DNA_binding","CBFD_NFYB_HMF","Alba","Bac_DNA_binding","MC1","Histone")


#For each protein annotating pfam hit
for( PFAM in unique(Bact_Res$query_name)){
  if(length(which(data2$UniprotID%in%Bact_Res[which(Bact_Res$query_name==PFAM),]$ID))>0){
    print(PFAM)
    data2[,PFAM]=0
    data2[which(data2$UniprotID%in%Bact_Res[which(Bact_Res$query_name==PFAM),]$ID),PFAM]=1}
}



#Annotating all proteins containing a DNA binding domain : 
data2$PFAM_DNA_binding=0
data2[which(data2$UniprotID%in%Bact_Res[which(Bact_Res$query_name%in%PfamDNAbinding),]$ID),]$PFAM_DNA_binding=1

#Annotating all NAPs
data2$NAP=0
if(length(which(data2$UniprotID%in%Bact_Res[which(Bact_Res$query_name%in%NAPsHMM),]$ID))>0){
  
data2[which(data2$UniprotID%in%Bact_Res[which(Bact_Res$query_name%in%NAPsHMM),]$ID),]$NAP=1
}




#Export
write.table(data2,file="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/LEPTOSPIRA/INHOUSE/Leptospira_Interrogans_MassSpec_PFAMAnno.txt",sep="\t",quote=F,row.names=F)






#Re-load
data=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/LEPTOSPIRA/INHOUSE/Leptospira_Interrogans_MassSpec_PFAMAnno.txt",header=T,sep="\t",quote="",stringsAsFactors = F)

# get LFQ column numbers
LFQ_columns <- grep("LFQ.", colnames(data))

#Removed the pooled measurements
LFQ_columns=LFQ_columns[c(1:3,6:10)]

names(data)[LFQ_columns]=gsub("INT_","",names(data)[LFQ_columns])


#This is required for dep
data$Gene.names=data$UniprotID
data$Protein.IDs=data$UniprotID



ShortNames=gsub("LFQ.intensity.","",names(data)[LFQ_columns])

#Exp.design : 
experimental_design=as.data.frame(matrix(nrow=8,ncol=3))
names(experimental_design)=c("label","condition","replicate")
experimental_design$label=names(data)[LFQ_columns]
experimental_design$condition=unlist(lapply(ShortNames,function(x) strsplit(x,"_")[[1]][1]))
experimental_design$replicate=unlist(lapply(ShortNames,function(x) strsplit(x,"_")[[1]][2]))



data_results <- LFQ(data, experimental_design,fun ="man",control = "top" ,type = "all", alpha = 0.05, lfc = 0)

AA=data_results$results

#Computing another p-value correction, as the initial algorythm deems no gene sig, which is not correct and due to multiple testing pval correction :
AA$nuc_vs_top_p.val.FDR=p.adjust(AA$nuc_vs_top_p.val,method = "fdr")

#Annotating dna binding proteins
AA$PFAM_DNA_binding=0
AA[which(AA$ID%in%data[which(data$PFAM_DNA_binding==1),]$Protein.IDs),]$PFAM_DNA_binding=1

#Export : 
write.table(AA,file="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/LEPTOSPIRA/INHOUSE/Leptospira_Interrogans_MassSpec_Enrichment.txt",sep="\t",quote=F,row.names=F)





volcanoPlot=ggplot(data=AA)+geom_point(aes(x=nuc_vs_top_ratio,y=-log10(nuc_vs_top_p.val),color=as.factor(PFAM_DNA_binding)),size=0.5)+geom_point(data=AA[which(AA$name=="A0A2H1XGH2"),],aes(x=nuc_vs_top_ratio,y=-log10(nuc_vs_top_p.val)),color="orange",size=2)+geom_point(data=AA[which(AA$name=="A0A2H1XG06"),],aes(x=nuc_vs_top_ratio,y=-log10(nuc_vs_top_p.val)),color="red",size=2)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))

ggsave(plot = volcanoPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/LEPTOSPIRA/INTERROGANS/Nuc_vs_Top_volcanoPlotUncorrectedpval.pdf")


volcanoPlot=ggplot(data=AA)+geom_point(aes(x=nuc_vs_top_ratio,y=-log10(PadjBH),color=as.factor(PFAM_DNA_binding)),size=0.5)+geom_point(data=AA[which(AA$name=="A0A2H1XGH2"),],aes(x=nuc_vs_top_ratio,y=-log10(PadjBH)),color="orange",size=2)+geom_point(data=AA[which(AA$name=="A0A2H1XG06"),],aes(x=nuc_vs_top_ratio,y=-log10(PadjBH)),color="red",size=2)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))
ggsave(plot = volcanoPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/LEPTOSPIRA/INTERROGANS/Nuc_vs_Top_volcanoPlot.pdf")



WCEPlot=merge(AA,data,by.x="name",by.y="UniprotID",all=T)





EnrichmentPlot=ggplot(data=WCEPlot)+geom_point(aes(x=nuc_vs_top_ratio,y=sqrt((LFQ.intensity.wce_1+LFQ.intensity.wce_2)/2),color=as.factor(PFAM_DNA_binding.x)),size=0.5)+geom_point(data=WCEPlot[which(WCEPlot$name=="A0A2H1XGH2"),],aes(x=nuc_vs_top_ratio,y=sqrt((LFQ.intensity.wce_1+LFQ.intensity.wce_2)/2)),color="orange",size=2)+geom_point(data=WCEPlot[which(WCEPlot$name=="A0A2H1XG06"),],aes(x=nuc_vs_top_ratio,y=sqrt((LFQ.intensity.wce_1+LFQ.intensity.wce_2)/2)),color="red",size=2)+theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_color_manual(values=c("grey80","blue"))+ylab("WCE LFQ intensity")
ggsave(plot = EnrichmentPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/LEPTOSPIRA/INTERROGANS/Nuc_vs_Top_enrichment_vs_WCE_Intensity_orange_candidate_blueDNAbinding.pdf")



WCEPlot=WCEPlot[order(-WCEPlot$LFQ.intensity.nuc_1),]


#Plotting correlations between different conditions

library(corrplot)
LFQ_columns <- grep("LFQ.", colnames(data)) # get LFQ column numbers
LFQ_columns[c(1:3,6:10)]
CorMat=data[,LFQ_columns]
CorMatRes=cor(CorMat)
pdf("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/LEPTOSPIRA/INTERROGANS/Linterrogans_PearsonCorrealtion_MS_replicates.pdf",height=6,width=6)
corrplot(CorMatRes)
dev.off()

CorMatRes=cor(CorMat,method="spearman")

pdf("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/LEPTOSPIRA/INTERROGANS/Linterrogans_SpearmanCorrealtion_MS_replicates.pdf",height=6,width=6)
corrplot(CorMatRes)
dev.off()






#Localisation : 


LocalisationResult=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/EXTERNAL_DATASETS/PSORT_LEPTOSPIRA/INTERROGANS/psortdb-results_Linterro_AllGenbank.tsv",header=T,sep="\t",stringsAsFactors = F,quote="")

LocalisationResult$ID

WCEMspec3=merge(WCEPlot,LocalisationResult,by.x="GenbankID",by.y="ID",all.x=T)


WCEPlot=WCEMspec3

WCEPlot$LFQ.intensity.Nuc_avg=(WCEPlot$LFQ.intensity.nuc_1+WCEPlot$LFQ.intensity.nuc_2+WCEPlot$LFQ.intensity.nuc_3)/3

WCEPlot$LFQ.intensity.Top_avg=(WCEPlot$LFQ.intensity.top_1+WCEPlot$LFQ.intensity.top_2+WCEPlot$LFQ.intensity.top_3)/3



WCEPlot$LFQ.intensity.WCE_avg=(WCEPlot$LFQ.intensity.wce_1+WCEPlot$LFQ.intensity.wce_2)/2



LocaAndName=ggplot(data=WCEPlot)+geom_point(aes(y=nuc_vs_top_ratio,x=log2(LFQ.intensity.Nuc_avg),color=as.factor(Final_Localization)),size=0.5)+geom_text_repel(data=WCEPlot[which(WCEPlot$name=="A0A2H1XGH2"),],aes(y=nuc_vs_top_ratio,x=log2(LFQ.intensity.Nuc_avg),label=name),color="black",size=1,min.segment.length = 0.1)+theme_pubclean()+theme(aspect.ratio=1/1.61,legend.position = "")+scale_color_manual(values=col_vector)+xlab("Intensity, nucleoid fraction (log2)")+ylab("Enrichment \nNucleoid vs Top fraction (log2)")
ggsave(addSmallLegend(LocaAndName),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/LEPTOSPIRA/INTERROGANS/Nuc_vs_Top_enrichment_scatter_log2scale.pdf")


LocaAndName=ggplot(data=WCEPlot)+geom_point(aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg,color=as.factor(Final_Localization)),size=0.5)+geom_text_repel(data=WCEPlot[which(WCEPlot$name=="A0A2H1XGH2"),],aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg,label=name),color="black",size=1,min.segment.length = 0.1)+geom_text_repel(data=WCEPlot[which(WCEPlot$name=="A0A2H1XG06"),],aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg,label=name),color="black",size=1,min.segment.length = 0.1)+theme_pubclean()+theme(aspect.ratio=1/1.61,legend.position = "top")+scale_color_manual(values=col_vector)+xlab("Intensity, nucleoid fraction")+ylab("Enrichment \nNucleoid vs Top fraction (log2)")
ggsave(addSmallLegend(LocaAndName),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/LEPTOSPIRA/INTERROGANS/Nuc_vs_Top_enrichment_scatter_linearscale.pdf")


WCEPlot=WCEPlot[order(-WCEPlot$LFQ.intensity.Nuc_avg),]
DNABindingFocus=ggplot(data=WCEPlot)+geom_point(aes(y=nuc_vs_top_ratio,x=log2(LFQ.intensity.Nuc_avg),color=as.factor(PFAM_DNA_binding.x)),size=0.5)+geom_text_repel(data=WCEPlot[which(WCEPlot$PFAM_DNA_binding.x==1 & WCEPlot$nuc_vs_top_ratio>1 & WCEPlot$PadjBH <0.01)[1:20],],aes(y=nuc_vs_top_ratio,x=log2(LFQ.intensity.Nuc_avg),label=ID),color="black",size=1,min.segment.length = 0.1)+theme_pubclean()+theme(aspect.ratio=1/1.61,legend.position = "top")+scale_color_manual(values=c("grey80","blue"))+xlab("Intensity, nucleoid fraction")+ylab("Enrichment \nNucleoid vs Top fraction (log2)")
ggsave(addSmallLegend(DNABindingFocus),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/LEPTOSPIRA/INTERROGANS/Nuc_vs_Top_enrichment_DNA_bindingFocus_scatter_log2scale.pdf")







#Adding descriscription 


WCEPlot$Description=unlist(lapply(WCEPlot$Fasta.headers,function(x) strsplit(x,"\\|")[[1]][3]))
WCEPlot$Description=unlist(lapply(WCEPlot$Description,function(x) strsplit(x,"\\OS=Leptospira interrogans serovar Manilae")[[1]][1]))


DNABindingFocus=ggplot(data=WCEPlot)+geom_point(aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg,color=as.factor(PFAM_DNA_binding.x)),size=0.5)+geom_text_repel(data=WCEPlot[which(WCEPlot$PFAM_DNA_binding.x==1 & WCEPlot$nuc_vs_top_ratio>1 & WCEPlot$PadjBH <0.01)[1:20],],aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg,label=ID),color="black",size=1,min.segment.length = 0.1)+theme_pubclean()+theme(aspect.ratio=1/1.61,legend.position = "top")+scale_color_manual(values=c("grey80","blue"))+xlab("Intensity, nucleoid fraction")+ylab("Enrichment \nNucleoid vs Top fraction")+scale_x_sqrt()
ggsave(addSmallLegend(DNABindingFocus),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/LEPTOSPIRA/INTERROGANS/Nuc_vs_Top_enrichment_DNA_bindingFocus_scatter_.pdf")



DNABindingFocus=ggplot(data=WCEPlot)+geom_point(aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg,color=as.factor(PFAM_DNA_binding.x)),size=0.5)+geom_point(data=WCEPlot[which(WCEPlot$ID=="A0A2H1XGH2"),],aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg,color=as.factor(PFAM_DNA_binding.x)),size=0.5,color="orange")+geom_text_repel(data=WCEPlot[which(WCEPlot$PFAM_DNA_binding.x==1 & WCEPlot$nuc_vs_top_ratio>1 & WCEPlot$PadjBH <0.01)[1:10],],aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg,label=Description),color="black",size=1,min.segment.length = 0.1,max.overlaps=1000)+theme_pubclean()+theme(aspect.ratio=1/1.61,legend.position = "top")+scale_color_manual(values=c("grey80","blue"))+xlab("Intensity, nucleoid fraction")+ylab("Enrichment \nNucleoid vs Top fraction")+scale_x_sqrt()
ggsave(addSmallLegend(DNABindingFocus),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/LEPTOSPIRA/INTERROGANS/Nuc_vs_Top_enrichment_DNA_bindingFocus_scatter_with_description.pdf")


#This time plotting WCE abundancy ( makes more sense in fact)
DNABindingFocus=ggplot(data=WCEPlot)+geom_point(aes(y=nuc_vs_top_ratio,x=LFQ.intensity.WCE_avg,color=as.factor(PFAM_DNA_binding.x)),size=0.5)+geom_point(data=WCEPlot[which(WCEPlot$ID=="A0A2H1XGH2"),],aes(y=nuc_vs_top_ratio,x=LFQ.intensity.WCE_avg,color=as.factor(PFAM_DNA_binding.x)),size=0.5,color="orange")+geom_text_repel(data=WCEPlot[which(WCEPlot$PFAM_DNA_binding.x==1 & WCEPlot$nuc_vs_top_ratio>1 & WCEPlot$PadjBH <0.01)[1:10],],aes(y=nuc_vs_top_ratio,x=LFQ.intensity.WCE_avg,label=Description),color="black",size=1,min.segment.length = 0.1,max.overlaps=1000)+theme_pubclean()+theme(aspect.ratio=1/1.61,legend.position = "top")+scale_color_manual(values=c("grey80","blue"))+xlab("Intensity, Whole cell extract")+ylab("Enrichment \nNucleoid vs Top fraction")+scale_x_sqrt()
ggsave(addSmallLegend(DNABindingFocus),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/LEPTOSPIRA/INTERROGANS/Nuc_vs_Top_enrichment_DNA_bindingFocus_scatter_with_description_xlabWCEIntensity.pdf")




WCEPlot$Final_Localization=factor(WCEPlot$Final_Localization,levels = rev(c("Outer Membrane","Cytoplasmic Membrane","Extracellular","Unknown","Periplasmic","Cytoplasmic")))

LocalisationPlot=ggplot(data=WCEPlot)+geom_violin(aes(y=nuc_vs_top_ratio,x=as.factor(Final_Localization),group=as.factor(Final_Localization),fill=as.factor(Final_Localization)),size=0,width=1)+geom_point(data=WCEPlot[which(WCEPlot$name=="A0A2H1XGH2"),],aes(y=nuc_vs_top_ratio,x=as.factor(Final_Localization),group=as.factor(Final_Localization)),color="black",size=2)+geom_point(data=WCEPlot[which(WCEPlot$name=="A0A2H1XG06"),],aes(y=nuc_vs_top_ratio,x=as.factor(Final_Localization),group=as.factor(Final_Localization)),color="black",size=2) +theme_pubclean()+theme(aspect.ratio=1,legend.position = "")+scale_fill_manual(values=col_vector)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10))+xlab("")+coord_flip()
ggsave(LocalisationPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/LEPTOSPIRA/INTERROGANS/Nuc_vs_Top_enrichment_stratified_by_localisation.pdf")
















#Additional exploratory plot ( to see what is abundant and enriched in the chromatin fraction)

DNABindingFocus=ggplot(data=WCEPlot)+geom_point(aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg,color=as.factor(PFAM_DNA_binding.x)),size=0.5)+geom_point(data=WCEPlot[which(WCEPlot$ID=="A0A2H1XGH2"),],aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg,color=as.factor(PFAM_DNA_binding.x)),size=0.5,color="orange")+geom_text_repel(data=WCEPlot[which(WCEPlot$nuc_vs_top_ratio>1 & WCEPlot$PadjBH <0.01)[1:75],],aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg,label=Description),color="black",size=1,min.segment.length = 0.1,max.overlaps=1000)+theme_pubclean()+theme(aspect.ratio=1/1.61,legend.position = "top")+scale_color_manual(values=c("grey80","blue"))+xlab("Intensity, nucleoid fraction")+ylab("Enrichment \nNucleoid vs Top fraction")+scale_x_sqrt()
ggsave(addSmallLegend(DNABindingFocus),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/LEPTOSPIRA/INTERROGANS/Nuc_vs_Top_enrichment_scatter_with_description_Top75.pdf",height = 12,width=12)


#Focus on Scc : 
DNABindingFocus=ggplot(data=WCEPlot)+geom_point(aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg),size=0.5)+geom_point(data=WCEPlot[which(WCEPlot$ID=="A0A2H1XHS2"),],aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg),size=1,color="red")+geom_point(data=WCEPlot[which(WCEPlot$ID=="A0A2H1XGH2"),],aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg),size=1,color="orange")+geom_point(data=WCEPlot[which(WCEPlot$ID=="A0A2H1XG06"),],aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg),size=1,color="green")+geom_text_repel(data=WCEPlot[which(WCEPlot$ID%in%c("A0A2H1XHS2","A0A2H1XGH2","A0A2H1XG06")),],aes(y=nuc_vs_top_ratio,x=LFQ.intensity.Nuc_avg,label=Description),color="black",size=3,min.segment.length = 0.1,max.overlaps=1000)+theme_pubr()+theme(aspect.ratio=1/1.61,legend.position = "top")+scale_color_manual(values=c("grey80","blue"))+xlab("Intensity, nucleoid fraction")+ylab("Enrichment \nNucleoid vs Top fraction")+scale_x_sqrt()
ggsave(addSmallLegend(DNABindingFocus),filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/LEPTOSPIRA/INTERROGANS/Nuc_vs_Top_enrichment_scatter_with_HFOrange_H1Green_SccRed.pdf",height = 12,width=12)




##############################
#Exporting data : 
##############################

#Subselecting columns to export : 
LFQ_columns <- grep("LFQ.", colnames(WCEPlot)) # get LFQ column numbers
WCEPlot.e=WCEPlot[,c(which(names(WCEPlot)%in%c(names(AA),"Final_Localization","Protein.Name","Signal_Details","PFAM_DNA_binding.y","NAP")),LFQ_columns)]

#Export Table with all data : 

write.table(WCEPlot.e,"/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/LEPTOSPIRA/INHOUSE/b023p057_Leptospira_Interrogans_MassSpec_Enrichment_w_localisation.txt",sep="\t",row.names=F,quote=F)









#######################
#Ranked abundancy plot: 
#######################
WCEPlot.e=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/LEPTOSPIRA/INHOUSE/b023p057_Leptospira_Interrogans_MassSpec_Enrichment_w_localisation.txt",header=T,sep="\t",quote="",stringsAsFactors = F)


WCEPlot.e$LFQ.intensity.wce_avg=(WCEPlot.e$LFQ.intensity.wce_1+WCEPlot.e$LFQ.intensity.wce_2)/2

WCEPlot.e=WCEPlot.e[-which(is.na(WCEPlot.e$LFQ.intensity.wce_avg)==T),]
WCEPlot.e=WCEPlot.e[order(WCEPlot.e$LFQ.intensity.wce_avg),]

WCEPlot.e$NAPname=NA
WCEPlot.e[which(WCEPlot.e$name=="A0A2H1XHS2"),]$NAPname="scc"
WCEPlot.e[which(WCEPlot.e$name=="A0A2H1XKC1"),]$NAPname="DPS"
WCEPlot.e[which(WCEPlot.e$name=="A0A2H1XG06"),]$NAPname="H1p"
WCEPlot.e[which(WCEPlot.e$name=="A0A2H1XGH2"),]$NAPname="HFp"


MassSpecPlot=addSmallLegend(ggplot(data=WCEPlot.e)+geom_segment(data=WCEPlot.e[which(is.na(WCEPlot.e$NAPname)==F),],aes(x=which(is.na(WCEPlot.e$NAPname)==F),xend=which(is.na(WCEPlot.e$NAPname)==F),y=log2(LFQ.intensity.wce_avg),yend=-Inf,color=NAPname))+geom_point(aes(x=1:dim(WCEPlot.e)[1],y=log2(LFQ.intensity.wce_avg)),shape=20,size=1,stroke=0.8,color="#0F4159")+geom_point(data=WCEPlot.e[which(is.na(WCEPlot.e$NAPname)==F),],aes(x=which(is.na(WCEPlot.e$NAPname)==F),y=log2(LFQ.intensity.wce_avg),color=NAPname),shape=1,size=1.8,stroke=0.8)+theme_pubr()+theme(aspect.ratio = 1/6)+xlab("Rank")+ylab("Protein abundance,\n (log2(LFQ))")+ggtitle(expression(paste(italic("L. interrogans"), sep="")))+scale_color_manual(values = c("#7FC97F","red","#FDC086","#FFFF99","green","blue","brown","grey20")))


ggsave(plot=MassSpecPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Leptospira_interrogans_WCE_abundancy_inhousedata.pdf",width=6,height=4)



