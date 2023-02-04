#####################################################
## Project: bacterial histones
## Analysis of mass spectrometry data :
## Bdellovibrio bacteriovorus
## Linking to local id and quality check
## focus DNA binding proteins 
## Date: December 2020
## Author: Antoine Hocher
####################################################


Sys.setenv(PATH="/Users/ahocher/opt/miniconda3/bin:/Users/ahocher/opt/miniconda3/condabin:/opt/local/bin:/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin")
#LOADING LIRARIES : 
library(ggplot2)
library(seqinr)
library(metablastr)
library(ape)
library(seqinr)
library(ggpubr)
library(rhmmer)
source("/Users/ahocher/Dropbox/Scripts/PLOT_LIBRARY/AH_plot_library.R")
source("/Users/ahocher/Dropbox/Scripts/ComputeProtStats_from_Fasta.R")

#Import Mass Spectrometry data: 
MassSpec=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/BDELLOVIBRIO/b023p050_AH_Bdellovibrio_proteinGroups_20210115.txt",header=T,sep="\t",stringsAsFactors = F,quote="")
names(MassSpec)[1]="ID"
head(MassSpec)



#Annotate DNA binding and all PFAM hits
#Loading Pfam HMM search (note for github : would have to be re-run : hmmsearch of all PFAM-A vs bdellovibrio proteome, it's not key here it's just to annotate DNA binding proteins:
Bact_Res=as.data.frame(read_tblout("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/HMM_RESULTS/EBMC_Bact-vs-Pfam.txt"))

Bact_Res$Taxid=unlist(lapply(Bact_Res$domain_name, function(x) strsplit(x,"_")[[1]][1]))

Bact_ResSub=Bact_Res[which(Bact_Res$Taxid=="264462"),]
Bact_ResSub$ID=unlist(lapply(Bact_ResSub$domain_name, function(x) strsplit(x,"@")[[1]][2]))
Bact_ResSub$ID=gsub("\\|","",Bact_ResSub$ID)


#Hmm corresponding to PFAM DNA binding , or to NAPs :

#File obtained from: Malhotra, S. & Sowdhamini, R. Collation and analyses of DNA-binding protein domain families from sequence and structural databanks. Mol. Biosyst. 11, 1110–1118 (2015).
PfamDNAbinding=read.table(file="/Users/ahocher/Dropbox/Laboratory/ArchealChromatinProteom/HMMs/DBPome/PFAM_DNA_binding_list.txt",stringsAsFactors = F)$V1

PfamDNAbinding=c(PfamDNAbinding,"TFIIE_alpha")

NAPsHMM=c("Cren7","7kD_DNA_binding","CBFD_NFYB_HMF","Alba","Bac_DNA_binding","MC1","Histone")



for( PFAM in unique(Bact_ResSub$query_name)){
  if(length(which(MassSpec$ID%in%Bact_ResSub[which(Bact_ResSub$query_name==PFAM),]$ID))>0){
    print(PFAM)
    MassSpec[,PFAM]=0
    MassSpec[which(MassSpec$ID%in%Bact_ResSub[which(Bact_ResSub$query_name==PFAM),]$ID),PFAM]=1}
}

head(MassSpec)
MassSpec[which(MassSpec$CBFD_NFYB_HMF==1),"ID"]



#Annotating all proteins containing a DNA binding domain : 
MassSpec$PFAM_DNA_binding=0
MassSpec[which(MassSpec$ID%in%Bact_ResSub[which(Bact_ResSub$query_name%in%PfamDNAbinding),]$ID),]$PFAM_DNA_binding=1

MassSpec$NAP=0
MassSpec[which(MassSpec$ID%in%Bact_ResSub[which(Bact_ResSub$query_name%in%NAPsHMM),]$ID),]$NAP=1

# Export
write.table(MassSpec,file="/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/BDELLOVIBRIO/b023p050_B_bacteriovorus_WCE_proteinGroups_PFAMANNO.txt",col.names =T,sep="\t",row.names = F,quote=F)




MassSpec=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/PROTEOMES/BDELLOVIBRIO/b023p050_B_bacteriovorus_WCE_proteinGroups_PFAMANNO.txt",header=T,sep="\t",quote="",stringsAsFactors = F)


#Graph : NAP, ranked abundancy plot : 


MassSpec$NAPs=NA
MassSpec[which(MassSpec$ID=="CAE77736.1"),]$NAPs="CBFD_NFYB_HMF"
MassSpec[which(MassSpec$ID=="CAE79514.1"),]$NAPs="IHFa"
MassSpec[which(MassSpec$ID=="CAE78670.1"),]$NAPs="IHFb"
MassSpec[which(MassSpec$ID=="CAE78181.1"),]$NAPs="HUa"
MassSpec[which(MassSpec$ID=="CAE79947.1"),]$NAPs="HUb"
MassSpec[which(MassSpec$ID=="CAE81096.1"),]$NAPs="YbaB"
MassSpec[which(MassSpec$ID=="CAE80413.1"),]$NAPs="Dps"


#I could not find dps in our data

#Averaging technical and biological replicates : 
MassSpec$AvgLFQIntensity=(MassSpec$LFQ.intensity.B.bacteriovorus_WCE1_1+MassSpec$LFQ.intensity.B.bacteriovorus_WCE1_2+MassSpec$LFQ.intensity.B.bacteriovorus_WCE2_1+MassSpec$LFQ.intensity.B.bacteriovorus_WCE2_2)/4
MassSpec=MassSpec[order(MassSpec$AvgLFQIntensity),]


MassSpecPlot=addSmallLegend(ggplot(data=MassSpec)+geom_segment(data=MassSpec[which(is.na(MassSpec$NAPs)==F),],aes(x=which(is.na(MassSpec$NAPs)==F),xend=which(is.na(MassSpec$NAPs)==F),y=log2(AvgLFQIntensity),yend=-Inf,color=NAPs))+geom_point(aes(x=1:dim(MassSpec)[1],y=log2(AvgLFQIntensity)),shape=20,size=1,stroke=0.8,color="#0F4159")+geom_point(data=MassSpec[which(is.na(MassSpec$NAPs)==F),],aes(x=which(is.na(MassSpec$NAPs)==F),y=log2(AvgLFQIntensity),color=NAPs),shape=1,size=1.8,stroke=0.8)+theme_pubr()+theme(aspect.ratio = 1/6)+xlab("Rank")+ylab("Protein abundance,\n (log2(LFQ))")+ggtitle(expression(paste(italic("B.bacteriovorus"), sep="")))+scale_color_manual(values = c("#7FC97F","red","#FDC086","#FFFF99","green","blue","brown")))

ggsave(plot=MassSpecPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Bdellovibrio_Bacteriovorus_MassSpecInhouse.pdf",width=6,height=4)


#another version with Bd3044
MassSpec[which(MassSpec$ID=="CAE80810.1"),]$NAPs="Bd3044"


MassSpec$AvgLFQIntensity=(MassSpec$LFQ.intensity.B.bacteriovorus_WCE1_1+MassSpec$LFQ.intensity.B.bacteriovorus_WCE1_2+MassSpec$LFQ.intensity.B.bacteriovorus_WCE2_1+MassSpec$LFQ.intensity.B.bacteriovorus_WCE2_2)/4
MassSpec=MassSpec[order(MassSpec$AvgLFQIntensity),]

MassSpecPlot=addSmallLegend(ggplot(data=MassSpec)+geom_segment(data=MassSpec[which(is.na(MassSpec$NAPs)==F),],aes(x=which(is.na(MassSpec$NAPs)==F),xend=which(is.na(MassSpec$NAPs)==F),y=log2(AvgLFQIntensity),yend=-Inf,color=NAPs))+geom_point(aes(x=1:dim(MassSpec)[1],y=log2(AvgLFQIntensity)),shape=20,size=1,stroke=0.8,color="#0F4159")+geom_point(data=MassSpec[which(is.na(MassSpec$NAPs)==F),],aes(x=which(is.na(MassSpec$NAPs)==F),y=log2(AvgLFQIntensity),color=NAPs),shape=1,size=1.8,stroke=0.8)+theme_pubr()+theme(aspect.ratio = 1/6)+xlab("Rank")+ylab("Protein abundance,\n (log2(LFQ))")+ggtitle(expression(paste(italic("B.bacteriovorus"), sep="")))+scale_color_manual(values = c("#7FC97F","red","#FDC086","#FFFF99","green","blue","brown","grey20")))

ggsave(plot=MassSpecPlot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Bdellovibrio_Bacteriovorus_MassSpecInhouse_withBd3044.pdf",width=6,height=4)


MassSpec$Fasta.headers=gsub("\\[Bdellovibrio bacteriovorus HD100]","",MassSpec$Fasta.headers)
MassSpec$NormInt=100*MassSpec$AvgLFQIntensity/sum(MassSpec$AvgLFQIntensity)
MassSpec=MassSpec[order(-MassSpec$NormInt),]
MassSpec$Fasta.headers=factor(MassSpec$Fasta.headers,levels = MassSpec$Fasta.headers)
DNA_binding_plot=ggplot(data=MassSpec[which(MassSpec$PFAM_DNA_binding==1 | MassSpec$NAP==1)[1:30],])+geom_bar(aes(x=Fasta.headers,y=NormInt,fill=as.factor(NAP)),stat="identity", width=0.7)+theme_pubr()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+ylab("% of proteome")+ggtitle("Top30 of proteins containing a PFAM DNA binding domain")+scale_fill_manual(values = c("darkblue","green"))+theme(aspect.ratio = 1/3,legend.position="")+xlab("")

ggsave(plot=DNA_binding_plot,filename = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Figures/MASS_SPEC/Bdellovibrio_Bacteriovorus_MassSpecInhouse_DNA_binding_PFAM.pdf",width=6,height=12)


