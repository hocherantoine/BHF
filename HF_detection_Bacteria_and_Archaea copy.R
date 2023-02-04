######################################################
## Project: BACTERIAL HISTONES
## Investigate HF in bacteria
## 
## 
## Aim : 
## Retrieve all HF protein from a proteome database
## 
## 
## 
## Author: Antoine Hocher
## Note : this version is for GitHub, but does contains some bits of code that have not made it in the manuscript at the time of uploading. The code is shown for reproducibility purposes, but to be run one would need to re-build the concat_bact95_hclust05.fasta database from the  list of ncbi genbank accession (too large to put on github)
#############################
#############################

#First part happens in the terminal, launched from R.
#
#clean up memory
rm(list = ls());gc()

#Gives the system path to R.
Sys.setenv(PATH="/Users/ahocher/opt/miniconda3/bin:/Users/ahocher/opt/miniconda3/condabin:/opt/local/bin:/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin")


#On bacterial genomes :

system("jackhmmer --tblout /Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/BactHistSeed_vs_concat_bact95_hclust05_Jackhmmtblout.txt /Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Seeds_sequences/BactHistSeed_for_Jack.faa /Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/DB_BACT95_HCLUST0.5/PROTEOMES_FASTA/AnnotatedGenomes/concat_bact95_hclust05.fasta > /Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/BactHistSeed_vs_concat_bact95_hclust05_Jackhmm.txt")

system("jackhmmer --tblout /Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/ArchtHistSeed_vs_concat_bact95_hclust05_Jackhmmtblout.txt /Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Seeds_sequences/ArchHistSeed_for_Jack.faa /Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/DB_BACT95_HCLUST0.5/PROTEOMES_FASTA/AnnotatedGenomes/concat_bact95_hclust05.fasta > /Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/ArchtHistSeed_vs_concat_bact95_hclust05_Jackhmm.txt")

system("jackhmmer --tblout /Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/BactDoubleHistSeed_vs_concat_bact95_hclust05_Jackhmmtblout.txt /Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Seeds_sequences/BactDoubleHistSeed_for_Jack.faa /Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/DB_BACT95_HCLUST0.5/PROTEOMES_FASTA/AnnotatedGenomes/concat_bact95_hclust05.fasta > /Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/BactDoubleHistSeed_vs_concat_bact95_hclust05_Jackhmm.txt")

system("jackhmmer --noali --tblout /Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/ArchDoubleHistSeed_vs_concat_bact95_hclust05_Jackhmmtblout.txt /Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Seeds_sequences/ArchDoubleHistSeed_for_Jack.faa /Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/DB_BACT95_HCLUST0.5/PROTEOMES_FASTA/AnnotatedGenomes/concat_bact95_hclust05.fasta > /Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/ArchDoubleHistSeed_vs_concat_bact95_hclust05_Jackhmm.txt")


##########################
#########PARTII
##########################


##
#Retrieving sequences to work on histone phylogeny :
#
library(rhmmer);library(seqinr)

##############
#Isolating the subset of sequences that make it : ie, filtering the hmmtblout files, and retrieve the a.a sequences of hits
###############

OutputSequences <- function(EvalT,HmmoIName,NumberOfDomains,PathToGenomes,GenomeInfoPath,PathToHmm){
  
  #Start
  Genomes=read.delim2(file = GenomeInfoPath,header=T,sep="\t",stringsAsFactors = F,quote = "")
  Hmm=as.data.frame(read_tblout(file = PathToHmm))
  Hmm$query_name=HmmoIName
  HmmoI=Hmm[which(Hmm$query_name==HmmoIName & Hmm$sequence_evalue<EvalT & Hmm$domain_number_reg%in%NumberOfDomains),]
  
  #Getting only the hits which have a defined number of domain 
  
  #clean up some memory space : 
  rm(Hmm)
  
  
  HmmoI$Taxid=unlist(lapply(HmmoI$domain_name,function(x) strsplit(x,split = "_")[[1]][1]))
  HmmoI$SpecieStrain=unlist(lapply(HmmoI$domain_name,function(x) paste(strsplit(strsplit(x,"\\@")[[1]][1],"_")[[1]][-1],collapse ="-")))
  HmmoI$Assembly=unlist(lapply(HmmoI$domain_name,function(x) strsplit(x,split = "\\$")[[1]][-1]))
  
    #EXPORT SEQUENCES NAMES : 
  SeqNames=paste(HmmoIName,"_Proteins.txt",sep="")
  write.table(HmmoI[!duplicated(HmmoI$domain_name),]$domain_name,SeqNames,sep="\t",quote=F,row.names = F,col.names = F)
  
  #EXTRACT AND DO SUBLIST OF SEQUENCE : 
  FastaNames=paste(HmmoIName,"_Proteins.fa",sep="")
  Command=paste("seqtk subseq ",dirname(PathToGenomes),"/",basename(PathToGenomes)," ",SeqNames," > ",FastaNames,sep="")
  system(Command)
  
  
  #This is simply to have a visual overview of potential outliers
  HmmSeq=read.fasta(FastaNames)
  pdf(paste(HmmoIName,"_Proteins_length.pdf",sep=""))
  hist(getLength(HmmSeq),breaks=length(HmmSeq),main=paste(HmmoIName,"Hits","Eval < ",EvalT),xlab="Protein Length")
  dev.off()
  
}


setwd("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Sequences/Final_Databases/")

EvalT=0.001
PathToGenomes="/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/DB_Archaea95_update012020/genome_assemblies_prot_fasta/ncbi-genomes-2020-05-14/AnnotatedGenomes/concatArchaeaEBMC2.faa"
GenomeInfoPath=paste(dirname(PathToGenomes),"/Annotated_genomes_info.txt",sep="")


OutputSequences(EvalT = EvalT,HmmoIName = "Bacterial_singlets_all_archaea",NumberOfDomains = c(1:10),PathToGenomes = PathToGenomes,GenomeInfoPath = GenomeInfoPath,PathToHmm = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/BactHistSeed_vs_concatArchaeaEBMC2_Jackhmmtblout.txt")


OutputSequences(EvalT = EvalT,HmmoIName ="Bacterial_doublets_all_archaea",NumberOfDomains = c(1:10),PathToGenomes = PathToGenomes,GenomeInfoPath = GenomeInfoPath,PathToHmm = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/BactDoubleHistSeed_vs__concatArchaeaEBMC2_Jackhmmtblout.txt")

OutputSequences(EvalT = EvalT,HmmoIName ="Archaeal_singlets_all_archaea",NumberOfDomains = c(1:10),PathToGenomes = PathToGenomes,GenomeInfoPath = GenomeInfoPath,PathToHmm = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/ArchtHistSeed_vs_concatArchaeaEBMC2_Jackhmmtblout.txt")



OutputSequences(EvalT = EvalT,HmmoIName ="Archaeal_doublets_all_archaea",NumberOfDomains = c(1:10),PathToGenomes = PathToGenomes,GenomeInfoPath = GenomeInfoPath,PathToHmm = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/ArchDoubleHistSeed_vs_concatArchaeaEBMC2_Jackhmmtblout.txt")



#Here we don't use jackhmm results but HMM models made obtained from PFAM which contain their own threshold (hence no further filtering on E-value - filtering at 1e-3 would likely remove some false positive but its minor)
#Because it's computed on cut-ga no EvalT used (ie set to 1)
#
OutputSequences(EvalT = 1,HmmoIName = "pfam_CBFD_NFYB_HMF_archaea",NumberOfDomains = c(1:10),PathToGenomes = PathToGenomes,GenomeInfoPath = GenomeInfoPath,PathToHmm = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/CBFD_NFYB_HMFPfam_vs_ArchaeaEBMC2_tblout.txt")

#DUF1931 :
OutputSequences(EvalT = 1,HmmoIName = "pfam_DUF1931_archaea",NumberOfDomains = c(1:10),PathToGenomes = PathToGenomes,GenomeInfoPath = GenomeInfoPath,PathToHmm = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/DUF1931Pfam_vs_ArchaeaEBMC2_tblout.txt")

#Histone :
OutputSequences(EvalT = 1,HmmoIName = "pfam_Histone_archaea",NumberOfDomains = c(1:10),PathToGenomes = PathToGenomes,GenomeInfoPath = GenomeInfoPath,PathToHmm = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/HistonePfam_vs_ArchaeaEBMC2_tblout.txt")







#For Bacteria HF proteins :
PathToGenomes="/Users/ahocher/Dropbox/Laboratory/Final_analysis_archaeal_chromatin/SEQ_DB/DB_BACT95_HCLUST0.5/PROTEOMES_FASTA/AnnotatedGenomes/concat_bact95_hclust05.fasta"
GenomeInfoPath=paste(dirname(PathToGenomes),"/Annotated_genomes_info.txt",sep="")

OutputSequences(EvalT = EvalT,HmmoIName = "Bacterial_singlets_bacteria",NumberOfDomains = c(1:10),PathToGenomes = PathToGenomes,GenomeInfoPath = GenomeInfoPath,PathToHmm = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/BactHistSeed_vs_concat_bact95_hclust05_Jackhmmtblout.txt")

OutputSequences(EvalT = EvalT,HmmoIName = "Bacterial_doublets_bacteria",NumberOfDomains = c(1:10),PathToGenomes = PathToGenomes,GenomeInfoPath = GenomeInfoPath,PathToHmm = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/BactDoubleHistSeed_vs_concat_bact95_hclust05_Jackhmmtblout.txt")

OutputSequences(EvalT = EvalT,HmmoIName = "Archaeal_singlets_bacteria",NumberOfDomains = c(1:10),PathToGenomes = PathToGenomes,GenomeInfoPath = GenomeInfoPath,PathToHmm = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/ArchtHistSeed_vs_concat_bact95_hclust05_Jackhmmtblout.txt")

OutputSequences(EvalT = EvalT,HmmoIName = "Archaeal_doublets_bacteria",NumberOfDomains = c(1:10),PathToGenomes = PathToGenomes,GenomeInfoPath = GenomeInfoPath,PathToHmm = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/ArchDoubleHistSeed_vs_concat_bact95_hclust05_Jackhmmtblout.txt")





#Because it's computed on cut-ga no EvalT used (ie set to 1)
OutputSequences(EvalT = 1,HmmoIName = "pfam_CBFD_NFYB_HMF_bacteria",NumberOfDomains = c(1:10),PathToGenomes = PathToGenomes,GenomeInfoPath = GenomeInfoPath,PathToHmm = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/CBFD_NFYB_HMFPfam_vs_BacteriaEBMCLarge_tblout.txt")

#DUF1931 :
OutputSequences(EvalT = 1,HmmoIName = "pfam_DUF1931_bacteria",NumberOfDomains = c(1:10),PathToGenomes = PathToGenomes,GenomeInfoPath = GenomeInfoPath,PathToHmm = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/DUF1931Pfam_vs_BacteriaEBMCLarge_tblout.txt")

#Histone :
OutputSequences(EvalT = 1,HmmoIName = "pfam_Histone_bacteria",NumberOfDomains = c(1:10),PathToGenomes = PathToGenomes,GenomeInfoPath = GenomeInfoPath,PathToHmm = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/JackHmmResults/HistonePfam_vs_BacteriaEBMCLarge_tblout.txt")







############Computing overlap plot overall size distribution :
############
setwd("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Sequences/Final_Databases/")

#All possible archaeal things :
CBFD_NFYB_HMF=read.fasta("pfam_CBFD_NFYB_HMF_archaea_Proteins.fa")
DUF1931=read.fasta("pfam_DUF1931_archaea_Proteins.fa")
Histone=read.fasta("pfam_Histone_archaea_Proteins.fa")
BactSingle=read.fasta("Bacterial_singlets_all_archaea_Proteins.fa")
BactDouble=read.fasta("Bacterial_doublets_all_archaea_Proteins.fa")
ArchSingle=read.fasta("Archaeal_singlets_all_archaea_Proteins.fa")
ArchDouble=read.fasta("Archaeal_doublets_all_archaea_Proteins.fa")

#Bacterial proteins
Bact_DUF1931=read.fasta("pfam_DUF1931_bacteria_Proteins.fa")
Bact_Histone=read.fasta("pfam_Histone_bacteria_Proteins.fa")
Bact_CBFD=read.fasta("pfam_CBFD_NFYB_HMF_bacteria_Proteins.fa")
Bact_BactSingle=read.fasta("Bacterial_singlets_bacteria_Proteins.fa")
Bact_BactDouble=read.fasta("Bacterial_doublets_bacteria_Proteins.fa")
Bact_ArchSingle=read.fasta("Archaeal_singlets_bacteria_Proteins.fa")
Bact_ArchDouble=read.fasta("Archaeal_doublets_bacteria_Proteins.fa")




AllHF=c(Bact_BactSingle,Bact_BactDouble,Bact_ArchSingle,Bact_ArchDouble,Bact_DUF1931,Bact_CBFD,Bact_Histone,ArchSingle,ArchDouble,BactSingle,BactDouble,CBFD_NFYB_HMF,DUF1931,Histone)

AllHF=AllHF[!duplicated(AllHF)]

#Plot length distribution 
postscript("All_HF_Archaea_and_Bacteria_Size_DistributionBeforeCutoff.eps")
hist(getLength(AllHF),breaks=1000)
dev.off()

#Order before export
AllHF=AllHF[order(getLength(AllHF))]
#Export
write.fasta(AllHF,names = getName(AllHF),"All_HF_archaea_and_bacteria.fa")

#Remove the * special character : 
system("tr -d '*' < All_HF_archaea_and_bacteria.fa > All_HF_archaea_and_bacterianoStar.fa")

#Export all less than 200bp : 
write.fasta(AllHF[which(getLength(AllHF)<200)],names = getName(AllHF[which(getLength(AllHF)<200)]),"All_HF_archaea_and_bacteria_less200.fa")

#Remove the * special character : 
system("tr -d '*' < All_HF_archaea_and_bacteria_less200.fa > All_HF_archaea_and_bacteria_less200_noStar.fa")



#Export separately Bacterial proteins 
AllHFBacteria=c(Bact_BactSingle,Bact_BactDouble,Bact_ArchSingle,Bact_ArchDouble,Bact_DUF1931,Bact_CBFD,Bact_Histone)
AllHFBacteria=AllHFBacteria[!duplicated(AllHFBacteria)]

AllHFBacteria=AllHFBacteria[order(getLength(AllHFBacteria))]
write.fasta(AllHFBacteria,names = getName(AllHFBacteria),"All_HF_bacteria.fa")
system("tr -d '*' < All_HF_bacteria.fa > All_HF_bacteria_noStar.fa")


write.fasta(AllHFBacteria[which(getLength(AllHFBacteria)<200)],names = getName(AllHFBacteria[which(getLength(AllHFBacteria)<200)]),"All_HF_less200_bacteria.fa")
system("tr -d '*' < All_HF_less200_bacteria.fa > All_HF_less200_bacteria_noStar.fa")






#Align : 
system("mafft-linsi --reorder /Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Sequences/Final_Databases/All_HF_less200_bacteria_noStar.fa > /Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Sequences/Final_Databases/All_HF_less200_bacteria_noStar_alignedlinsi.fa")

#Find best model : 
setwd("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Sequences/Final_Databases/")
system("/Users/ahocher/Dropbox/My\ Mac\ \(mrc-10681.local\)/Documents/Softwares/iqtree-2.0-rc1-MacOSX/bin/iqtree -s All_HF_less200_bacteria_noStar_alignedlinsi.fa -m MF")


#Format protein name for tree building : 
BacterialHF=read.fasta("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Sequences/Final_Databases/All_HF_less200_bacteria_noStar_alignedlinsi.fa")
BacterialHFNames=getName(BacterialHF)
BacterialHFNamesReduced=unlist(lapply(BacterialHFNames, function(x) strsplit(strsplit(x,"\\@")[[1]][2],"\\$")[[1]][1]))
write.fasta(BacterialHF,names = BacterialHFNamesReduced,"/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Sequences/Final_Databases/All_HF_less200_bacteria_noStar_alignedlinsiReducedName.fa")


#Compute tree using RaXmL : 
system("raxml-ng --all --msa All_HF_less200_bacteria_noStar_alignedlinsiReducedName.fa --seed 2 --threads 32 --force  --model LG+R6 --tree pars{10} --bs-trees 100")










#Rename Tree to enable iTol and genespy compatibility : 
library(ape)
BactTree <- read.tree("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Sequences/Final_Databases/TREES/All_HF_less200_bacteria_noStar_alignedlinsiReducedName.fa.raxml.support")


  BactTree$tip.label=paste("_",BactTree$tip.label,sep="")

write.tree(BactTree,"/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Sequences/Final_Databases/TREES/All_HF_less200_bacteria_noStar_alignedlinsiReducedName_RenamediTol.raxml.support")





#Export a file to derive gene context from genespy :
source("/Users/ahocher/Dropbox/Scripts/Export_Genspy_targets_from_fasta.R")

setwd("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/iTol_files/AllHF_Bact_Arch_FinalDb/")

ExportGeneSpyTargets("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Sequences/Final_Databases/All_HF_bacteria.fa")
ExportGeneSpyTargets("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Sequences/Final_Databases/All_HF_archaea_and_bacteria.fa")

source("/Users/ahocher/Dropbox/Scripts/ExportInfo_for_iTol.R")
WriteiTolReadyFiles(FastaFilePath = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Sequences/All_HF_archaea_and_bacteria.fa",TreeFilePath = "/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Sequences/All_HF_archaea_and_bacteria.fa_AlignedlinsiBMGEnoTrimTree.treefile")

#Exporting additional informations for iTol : 
source("/Users/ahocher/Dropbox/Scripts/table2itol-master/table2itol.R")

#Within the good folder : 
setwd("/Users/ahocher/Dropbox/Laboratory/Bacterial_Histones/Bacterial_histones_phylogeny/Sequences/Final_Databases/")

Histone_statsf$TID=paste("_",Histone_statsf$ProtName,sep="")

Histone_statsf$FullName=Histone_statsf$OriginalName
write.table(Histone_statsf[,which(names(Histone_statsf)%in%c("TID","order","phylum","class","FullName"))],"AllHF_Stats.txt",sep="\t",quote=F,row.names=F)

create_itol_files(infiles = "AllHF_Stats.txt",identifier = "TID",separator = "\t",conversion = 4,label = "FullName")

