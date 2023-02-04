


rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.

LocalisationResult=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/EXTERNAL_DATASETS/PSORT_LEPTOSPIRA/BIFLEXA/psortdb-results_biflexa.tsv",header=T,sep="\t",stringsAsFactors = F,quote="")
LocalisationResult2=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/EXTERNAL_DATASETS/PSORT_LEPTOSPIRA/BIFLEXA/psortdb-results_biflexa-2.tsv",header=T,sep="\t",stringsAsFactors = F,quote="")
LocalisationResult3=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/EXTERNAL_DATASETS/PSORT_LEPTOSPIRA/BIFLEXA/psortdb-results_biflexa-3.tsv",header=T,sep="\t",stringsAsFactors = F,quote="")

LocalisationResult=rbind(LocalisationResult,LocalisationResult2)
LocalisationResult=rbind(LocalisationResult,LocalisationResult3)


head(LocalisationResult)

LocalisationResult$SeqID=gsub("ref\\|","",LocalisationResult$SeqID)
library(rtracklayer)

#Load refseq first
INfo=readGFF("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/LEPTOSPIRA/GCF_000017685.1_ASM1768v1_genomic.gff")
INfo=as.data.frame(INfo)
head(INfo)
INfo=INfo[which(is.na(INfo$locus_tag)==F),]
for(i in INfo$locus_tag){
  INfo[which(INfo$locus_tag==i),]$old_locus_tag=  names(table(unlist(INfo[which(INfo$locus_tag==i),]$old_locus_tag)))[1]
}
head(INfo)

INfo$old_locus_tag

#Loading genbank gff second. 
INfo2=readGFF("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/LEPTOSPIRA/GCA_000017685.1_ASM1768v1_genomic.gff")
INfo2$locus_tag
INfo2=as.data.frame(INfo2[,which(names(INfo2)%in%c("locus_tag","protein_id"))])
INfo=as.data.frame(INfo[,which(names(INfo)%in%c("old_locus_tag","protein_id"))])

names(INfo)[1]="locus_tag"
names(INfo)[2]="ProteinIDRefseq"

INfo=INfo[which(is.na(INfo$ProteinIDRefseq)==F),]

INfo2=INfo2[which(is.na(INfo2$protein_id)==F),]

names(INfo2)[2]="ID"
INfo$locus_tag=as.character(INfo$locus_tag)
INfo2$locus_tag=as.character(INfo2$locus_tag)

Correspondence=merge.data.frame(INfo,INfo2,by="locus_tag")


Correspondence=Correspondence[which(is.na(Correspondence$locus_tag)==F & is.na(Correspondence$ID)==F),]


LocalisationResult$ProteinIDRefseq=LocalisationResult$SeqID

LocalisationResultAnno=merge(LocalisationResult,Correspondence,by="ProteinIDRefseq")

write.table(LocalisationResultAnno,"/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/EXTERNAL_DATASETS/PSORT_LEPTOSPIRA/BIFLEXA/psortdb-results_biflexa_AllGenbank.tsv",sep="\t",quote=F,row.names=F)




###PART II : Leptospira interrogans :
###

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.

LocalisationResult=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/EXTERNAL_DATASETS/PSORT_LEPTOSPIRA/INTERROGANS/psortdb-results_Linterro1.tsv",header=T,sep="\t",stringsAsFactors = F,quote="")
LocalisationResult2=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/EXTERNAL_DATASETS/PSORT_LEPTOSPIRA/INTERROGANS/psortdb-results_Linterro2.tsv",header=T,sep="\t",stringsAsFactors = F,quote="")
LocalisationResult3=read.table("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/EXTERNAL_DATASETS/PSORT_LEPTOSPIRA/INTERROGANS/psortdb-results_Linterro3.tsv",header=T,sep="\t",stringsAsFactors = F,quote="")

LocalisationResult=rbind(LocalisationResult,LocalisationResult2)
LocalisationResult=rbind(LocalisationResult,LocalisationResult3)


head(LocalisationResult)

LocalisationResult$SeqID=gsub("ref\\|","",LocalisationResult$SeqID)
library(rtracklayer)

#Load refseq first
INfo=readGFF("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/LEPTOSPIRA/genome_assemblies_genome_gff_ManilaeHP/GCF_001047655.1_ASM104765v1_genomic.gff")
INfo=as.data.frame(INfo)
head(INfo)
INfo=INfo[which(is.na(INfo$locus_tag)==F),]
for(i in INfo$locus_tag){
  INfo[which(INfo$locus_tag==i),]$old_locus_tag=  names(table(unlist(INfo[which(INfo$locus_tag==i),]$old_locus_tag)))[1]
}
head(INfo)

INfo$old_locus_tag

#Loading genbank gff second. 
INfo2=readGFF("/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/LEPTOSPIRA/genome_assemblies_genome_gff_ManilaeHP/GCA_001047655.1_ASM104765v1_genomic.gff")
INfo2$locus_tag
INfo2=as.data.frame(INfo2[,which(names(INfo2)%in%c("locus_tag","protein_id"))])
INfo=as.data.frame(INfo[,which(names(INfo)%in%c("old_locus_tag","protein_id"))])

names(INfo)[1]="locus_tag"
names(INfo)[2]="ProteinIDRefseq"

INfo=INfo[which(is.na(INfo$ProteinIDRefseq)==F),]

INfo2=INfo2[which(is.na(INfo2$protein_id)==F),]

names(INfo2)[2]="ID"
INfo$locus_tag=as.character(INfo$locus_tag)
INfo2$locus_tag=as.character(INfo2$locus_tag)

Correspondence=merge.data.frame(INfo,INfo2,by="locus_tag")


Correspondence=Correspondence[which(is.na(Correspondence$locus_tag)==F & is.na(Correspondence$ID)==F),]


LocalisationResult$ProteinIDRefseq=LocalisationResult$SeqID

LocalisationResultAnno=merge(LocalisationResult,Correspondence,by="ProteinIDRefseq")

write.table(LocalisationResultAnno,"/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/EXTERNAL_DATASETS/PSORT_LEPTOSPIRA/INTERROGANS/psortdb-results_Linterro_AllGenbank.tsv",sep="\t",quote=F,row.names=F)

