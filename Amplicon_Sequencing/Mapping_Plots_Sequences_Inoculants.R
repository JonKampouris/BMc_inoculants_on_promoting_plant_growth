library(ggnetwork)
library(ggplot2) #ggplot2: Fancy graphs
library(ggpubr) #ggplot expansion for more fancy graphs
library(readxl) # upload files  
library(readr) # read the files
library(tidyr) # Data handling
library(RColorBrewer) # Colours for fancy graphs
library(tibble)# Data handling: rownames/columnanmes
library(stringr) # Manipulate strings
library(igraph)
library(ggraph)
library(dplyr)
library(reshape2)
library(ape)
library(msa)
library(phangorn)
library(ggtree)
library(seqinr)
source(paste0("~/functions_load.R"))


fungi1=  read.table(paste0("~/otu_RH.fungi.Ioannis.txt"), sep=" ")


fungi2= read.table(paste0("~/CT2020_RH.Int_CtrlvsInt_BM_summary_final_OTU.txt"), sep="\t") %>%filter(logFC>0)%>%filter(genus=="Trichoderma")%>%
  mutate(ID=str_replace(`Row.names`, "_.*", ""))




OTU_Table_count1=read.table(paste0("~/otu_RH.fungi.Ioannis.txt"), sep=" ")%>%
  rownames_to_column(var="Col") %>%full_join(select(fungi_RA, Col, OTU), by="Col")%>%
  select(-Col)%>%column_to_rownames(var="OTU")%>%t()%>%
  as.data.frame()


Fungi_responders= read_csv("~/Fungal_Responders_OTU_levels_MS3.csv")%>%
  mutate(OTU=str_replace(Taxa, ".*.;",""))%>%filter(coeff>0)

fungi_taxonomy=   read.table(paste0("~ITS2_CT_2020_TAX_reduced.txt"), sep="\t")%>%
  
  rownames_to_column(var="Col")%>%full_join(fungi_RA%>%select(Col, OTU))%>%filter(
    OTU!="NA"
  )%>%filter(OTU%in%c(Fungi_responders$OTU)&genus=="Trichoderma")


ITS_Trichoderma_OTUs_CT <- read_excel("~/ITS_Trichoderma_OTUs_CT.xlsx") %>%filter(OTU%in%c(fungi_taxonomy$Col))%>%select(Col=OTU, everything())%>%
  full_join(fungi_taxonomy[,c("Col", "OTU")])

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"OTU"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"Sequence"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

writeFasta(ITS_Trichoderma_OTUs_CT, "~/Trichoderma_responders.fa")
ReadFasta<-function(file) {
  # Read the file line by line
  fasta<-read_lines(file)
  # Identify header lines
  ind<-grep(">", fasta)
  # Identify the sequence lines
  s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
  # Process sequence lines
  seqs<-rep(NA, length(ind))
  for(i in 1:length(ind)) {
    seqs[i]<-paste(fasta[s$from[i]:s$to[i]], collapse = "")
    print(i)
  }
  # Create a data frame 
  DF<-data.frame(Sample=gsub(">", "", fasta[ind]), sequence=seqs)
  # Return the data frame as a result object from the function
  return(DF)
}


library(Biostrings)
library(vegan)
Trichoderma_Responders= readBStringSet("~/Trichoderma_responders.fa", format="fasta")
NCBI_Trichoderma_Sequences= readBStringSet("~/Trich.fasta", format = "fasta")
OMG16_Trichoderma_Sequences= readBStringSet("~/np_ITS_regions.fa", format = "fasta")
OMG16_combined= c(DNAStringSet(Trichoderma_Responders),DNAStringSet(NCBI_Trichoderma_Sequences),DNAStringSet(OMG16_Trichoderma_Sequences))


OMG16_alignment= msa(OMG16_combined, method = "ClustalW")
OMG16_alignment_seqinr <- msaConvert(OMG16_alignment, type = "seqinr::alignment")
OMG16_dist <- seqinr::dist.alignment(OMG16_alignment_seqinr,matrix = "identity")


OMG16_distNMDS= metaMDSiter(OMG16_dist,k=2)


info_tree_OMG16= data.frame(label=(names(OMG16_combined)), info=(names(OMG16_combined)), Category=(names(OMG16_combined)))%>%
  mutate(info=str_replace(info, 'OTU*.*',"OTU"))%>%  mutate(info=str_replace(info, 'NR*.*',"Trichoderma spp. Type Strains"))%>%
  mutate(info=str_replace(info, 'ITS.*',"OMG16"))
  
annotation=info_tree_OMG16%>%column_to_rownames(var="label")%>%select(`ITS sequences`=info)

anlist1=list( `ITS sequences` = c(OMG16= "black",OTU=  "darkred",`Trichoderma spp. Type Strains`= "skyblue4") )
htmpt= pheatmap::pheatmap( as.matrix(OMG16_dist), 
                    show_rownames=F,
                    show_colnames=F,cutree_cols=6, cutree_rows=6, annotation_col=annotation,  annotation_row=annotation,
                    legend = T,annotation_names_row=F, annotation_names_col=F,
                    fontsize = 30, annotation_colors=anlist1 )

bacilli_responders= read.csv(file = "~/Resdonders_ASVs_log_regression_MS3.csv")%>%filter(coeff>0)%>%
  filter(str_detect(Taxa, "Bacillus"))%>%mutate(ASV=str_replace(Taxa, ".*.*.ASV", "ASV"))
bac_seq <- read.csv("~/plastid_clean_low_read_clean_contaminant_free_ASV_Table_LTE_2020.csv", row.names = 1)%>%
  filter(ASV%in%c(bacilli_responders$ASV))
writeFasta( bac_seq%>%select(OTU=ASV, everything()) , "~/ASVs_Bacillus.fa")



ASV_Bacillus_Responders= readBStringSet("~/ASVs_Bacillus.fa", format="fasta")
NCBI_Bacillus_Sequences= readBStringSet("~/type_16SrRNA_regionBacillus_newfa.fasta", format = "fasta")
Bacillus_combined= c(DNAStringSet(ASV_Bacillus_Responders),DNAStringSet(NCBI_Bacillus_Sequences))

bacilli_alignment= msa(Bacillus_combined, method = "ClustalW")
bacilli_alignment_seqinr <- msaConvert(bacilli_alignment, type = "seqinr::alignment")
bacilli_dist <- seqinr::dist.alignment(bacilli_alignment_seqinr,matrix = "identity")
#bacilli_distNMDS= metaMDSiter(bacilli_dist, k=2)
info_tree_bacilli= data.frame(label=(names(Bacillus_combined)), info=(names(Bacillus_combined)), Category=(names(Bacillus_combined)))%>%
  mutate(info=str_replace(info, 'ASV*.*',"ASV"))%>% mutate(info=str_replace(info, '.*.s atrophaeus.*.',"Bacillus spp. Type Strains"))%>% 
  mutate(info=str_replace(info, 'NR*.*',"Bacillus spp. Type Strains"))%>%
  mutate(info=str_replace(info, 'Bac_.*',"B. atrophaeus Isolate"))%>% mutate(Category=str_replace(Category, 'Bac_atrophaeus*',"B. atrophaeus Isolate"))%>% 
  mutate(Category=str_replace(Category, '_Bacillus',""))%>%
  #mutate(Category=str_replace(Category, '.*Bacillus',"B."))%>%
  
  mutate(Category=str_replace(Category, 'strain.*.',""))%>%
  mutate(Category=str_replace(Category, '.riboso.*',""))%>%
  mutate(Category=str_replace(Category, 'Bac_*.',"Bacillus Atropheus"))%>%
  mutate(Category=str_replace(Category, 'NCPPB.*.',""))%>%
  mutate(Category=str_replace(Category, 'NBRC.*.',""))%>%
  mutate(Category=str_replace(Category, 'DSM.*',""))%>%
  mutate(Category=str_replace(Category, 'MLS10.*',""))%>%
  mutate(Category=str_replace(Category, 'ATCC.*',""))%>%
  mutate(Category=str_replace(Category, "]",""))%>%
  mutate(Category=str_replace(Category, "safensis.*","safensis"))%>%
  mutate(Category=str_replace(Category, "coahuilensis.*","coahuilensis"))%>%
  mutate(Category=str_replace(Category, "Cons.*","Consensus"))%>%
  mutate(info=str_replace(info, "Isolate","Abi03"))%>%column_to_rownames(var="label")%>%select(`16S rRNA gene`=info)


anlist1=list( `16S rRNA gene` = c(`B. atrophaeus Abi03`= "black",ASV=  "darkred",`Bacillus spp. Type Strains`= "skyblue4") )

htmpt= pheatmap::pheatmap( as.matrix(bacilli_dist), 
                           show_rownames=F,
                           show_colnames=F,cutree_cols=6, cutree_rows=6, annotation_col=info_tree_bacilli,  annotation_row=info_tree_bacilli,
                           legend = T,
                           fontsize = 30, annotation_names_row=F, annotation_names_col=F, treeheight_row =40, annotation_colors=anlist1)


p2=ggplotify::as.ggplot(htmpt)+ theme(legend.key.size = unit(13,"cm"))
