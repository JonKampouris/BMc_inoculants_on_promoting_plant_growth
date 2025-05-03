library(ggplot2) #ggplot2: Fancy graphs
library(ggpubr) #ggplot expansion for more fancy graphs
library(readxl)
library(readr)
library(dplyr)# Data handling
library(vegan)# Multivariate stat tool
library(tidyr) # Data handling
library(RColorBrewer) # Colours for fancy graphs
library(tibble)# Data handling: rownames/columnanmes
library(phyloseq) # tool for 16S amplicon analysis
library(reshape2)  # Data handling
library(dunn.test)# non-parametric post-hoc test
library(ape) # tool for phylogenetic trees
library(phangorn) #tool for phylogenetic trees
library(stringr)# string manipulation and text modification
library(plyr)
library(ARTool)
library(forcats)
library(tidyverse)
library(pairwiseAdonis)

#Script for alpha and beta-diversity for Rhizosphere and 
# Soil 16S rRNA amplicon: Dataset Bernburg 2020
# Assume that you have your DNA concentrations or the Copies per qpcr
ASVtable_working <- read.csv("~/plastid_clean_ASV_Table_LTE_2020.csv")
DNA_NG=readxl::read_excel("~/Di_Control_DNA_Conc_2019-2020.xlsx", 
sheet = "Dcont")
#Load the copies
Copies_NG=readxl::read_excel("~/Di_Control_DNA_Conc_2019-2020.xlsx", 
                          sheet = "Dcont2")%>%na.omit() %>%aggregate(Copies~Sample, mean)
Copies_DNA= full_join(DNA_NG, Copies_NG, by="Sample")
#In this dataset MW is  concentration and einwage the weight. 
set.seed(17011990)

#Function to normalize the ASV table 
normalise_100 <- function(x){100*(x/sum(x))}

#Load the Cleaned from plastid sequences table
ASVtable_working <- read.csv("~/plastid_clean_ASV_Table_LTE_2020.csv")
#Load the Taxonomy Table of the Rarefied ASV table. 
nrow(ASVtable_working%>%
       filter(., str_detect(Kingdom, "Archaea")))

ASVtax_table <- ASVtable_working[,2:9] %>%as.data.frame()
ASVtax_table <- ASVtax_table %>% mutate(Taxonomy=paste0(Kingdom,";",Phylum,";",Class,";",Order,";",Family,";",Genus))%>%
  mutate(Taxonomy =gsub("[[:space:]]", "_", Taxonomy))%>%
  mutate(Taxonomy=str_replace_all(Taxonomy, ".ncultured.*.", "Unclassified"))%>%
  mutate(Taxonomy=str_replace_all(Taxonomy, "uncultured*", "Unclassified"))%>%
  mutate(Taxonomy=str_replace(Taxonomy, "metagenome", "Unclassified"))

ASVtax_table <- select(ASVtax_table,ASV,Sequence ,Taxonomy)%>%
  separate(Taxonomy,
           into = c("Kingdom","Phylum","Class","Order","Family","Genus"), sep=";")
ASVtax_table[ASVtax_table=="NA"]="Unclassified"
ASVtax_table[ASVtax_table=="uncultured"]<-"Unclassified"
ASVtax_table[ASVtax_table=="metagenome"]<-"Unclassified"
rownames(ASVtax_table)=ASVtax_table$ASV
gathered_ASV_table=gather(ASVtable_working[c(1,3:114)], -ASV, key = "Sample", value = "Abundance")
#Assign metadata via sample names
metadata_assign=data.frame(Sample=unique(gathered_ASV_table$Sample))%>%
  mutate(Category=Sample)%>%
  mutate(Category=str_replace(Category, "FSRW", "Field-Soil/Root-Window;"))%>%
  mutate(Category=str_replace(Category, "RHRW", "Rhizosphere/Root-Window;"))%>%
  mutate(Category=str_replace(Category, "RW", "Root-Window;"))%>%
  mutate(Category=str_replace(Category, "FS", "Field-Soil;"))%>%
  mutate(Category=str_replace(Category, "RH", "Rhizosphere;"))%>%
  mutate(Category=str_replace(Category, ".T0", "T0-Soil;"))%>%
  mutate(., Category=str_replace(Category, ".raw", ""))%>%
  separate(., Category, into = c("Type", "Rest"), sep=";")%>%
  mutate(., Rest=str_replace(Rest, ".CT.Ext.C", "CT;Extensive;Ctrl;"))%>%
  mutate(., Rest=str_replace(Rest, ".CT.Ext.I", "CT;Extensive;BMc;"))%>%
  mutate(., Rest=str_replace(Rest, ".CTEx.I", "CT;Extensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTEx.C", "CT;Extensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTIn.I", "CT;Intensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTIn.C", "CT;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTIn.I", "CT;Intensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTIn.C", "CT;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTIn.C", "CT;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTIn.C", "CT;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".CT.Int.C", "CT;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".CT.Int.I", "CT;Intensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTEx.C", "CT;Extensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".CT.Int.I", "CT;Intensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, "CTExtT0-Soil.", "CT;Extensive;T0-Soil;")) %>%
  mutate(., Rest=str_replace(Rest, "CTIntT0-Soil.", "CT;Intensive;T0-Soil;")) %>%
  mutate(., Rest=str_replace(Rest, "MPExtT0-Soil.", "MP;Extensive;T0-Soil;")) %>%
  mutate(., Rest=str_replace(Rest, "MPIntT0-Soil.", "MP;Intensive;T0-Soil;")) %>%
  mutate(., Rest=str_replace(Rest, ".MP.Ext.C", "MP;Extensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".MP.Ext.I", "MP;Extensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPEx.I", "MP'Extensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPEx.C1", "MP;Extensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPIn.I", "MP;Intensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPIn.C", "MP;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPIn.I", "MP;Intensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPIn.C", "MP;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPIn.C", "MP;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPIn.C", "MP;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".MP.Int.C", "MP;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".MP.Int.I", "MP;Intensive;BMc;")) %>%
  separate(., Rest, into=(c("Practice1", "N-Fertilization_Intensity", "Treatment", "Block")), sep = ";")%>% 
  filter(., Practice1=="CT")


#Load the Taxonomy Table of the Rarefied ASV table. 
rownames(ASVtax_table)=ASVtax_table$ASV
gathered_ASV_table=gather(ASVtable_working[c(1,3:114)], -ASV, key = "Sample", value = "Abundance")
#Assign metadata via sample names 
#If you have any other metadata file please use it.
metadata_assign=data.frame(Sample=unique(gathered_ASV_table$Sample))%>%
  mutate(Category=Sample)%>%
  mutate(Category=str_replace(Category, "FSRW", "Root-Associated_Soil/Root-Window;"))%>%
  mutate(Category=str_replace(Category, "RHRW", "Rhizosphere/Root-Window;"))%>%
  mutate(Category=str_replace(Category, "RW", "Root-Window;"))%>%
  mutate(Category=str_replace(Category, "FS", "Root-Associated_Soil;"))%>%
  mutate(Category=str_replace(Category, "RH", "Rhizosphere;"))%>%
  mutate(Category=str_replace(Category, ".T0", "T0-Soil;"))%>%
  mutate(., Category=str_replace(Category, ".raw", ""))%>%
  separate(., Category, into = c("Type", "Rest"), sep=";")%>%
  mutate(., Rest=str_replace(Rest, ".CT.Ext.C", "CT;Extensive;Ctrl;"))%>%
  mutate(., Rest=str_replace(Rest, ".CT.Ext.I", "CT;Extensive;BMc;"))%>%
  mutate(., Rest=str_replace(Rest, ".CTEx.I", "CT;Extensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTEx.C", "CT;Extensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTIn.I", "CT;Intensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTIn.C", "CT;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTIn.I", "CT;Intensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTIn.C", "CT;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTIn.C", "CT;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTIn.C", "CT;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".CT.Int.C", "CT;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".CT.Int.I", "CT;Intensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".CTEx.C", "CT;Extensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".CT.Int.I", "CT;Intensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, "CTExtT0-Soil.", "CT;Extensive;T0-Soil;")) %>%
  mutate(., Rest=str_replace(Rest, "CTIntT0-Soil.", "CT;Intensive;T0-Soil;")) %>%
  mutate(., Rest=str_replace(Rest, "MPExtT0-Soil.", "MP;Extensive;T0-Soil;")) %>%
  mutate(., Rest=str_replace(Rest, "MPIntT0-Soil.", "MP;Intensive;T0-Soil;")) %>%
  mutate(., Rest=str_replace(Rest, ".MP.Ext.C", "MP;Extensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".MP.Ext.I", "MP;Extensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPEx.I", "MP'Extensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPEx.C1", "MP;Extensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPIn.I", "MP;Intensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPIn.C", "MP;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPIn.I", "MP;Intensive;BMc;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPIn.C", "MP;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPIn.C", "MP;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".MPIn.C", "MP;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".MP.Int.C", "MP;Intensive;Ctrl;")) %>%
  mutate(., Rest=str_replace(Rest, ".MP.Int.I", "MP;Intensive;BMc;")) %>%
  separate(., Rest, into=(c("Practice1", "N-Fertilization_Intensity", 
"Treatment", "Block")), sep = ";")%>% 
  filter(., Practice1=="CT")%>%
  mutate(Block=str_replace(Block,"..2020*.*",""))

FS_RS_metadata=filter(metadata_assign, 
Type%in%c( "Root-Associated_Soil","Rhizosphere"))
ASVtable=dplyr::select(ASVtable_working, ASV,
c(FS_RS_metadata$Sample))%>%na.omit()%>%column_to_rownames(var="ASV")
metadata_assign$Treatment=factor(metadata_assign$Treatment, c("Ctrl","BMc"))
metadata_assign$Block=as.numeric(metadata_assign$Block)
#Join the tables. You might need to change the code based on your metadata

#I trained my model on the with the Field Soil dataset, because I lack for the water content of rhizosphere

metadata2= full_join( metadata_assign, Copies_DNA%>%select(-Sample), 
by=c("Practice1", "N-Fertilization_Intensity", "Type", 
     "Treatment", "Block"))%>%na.omit()
ASVtable2= ASVtable_working[, c(3, 10:ncol(ASVtable_working))]
ASVtable2[,2:ncol(ASVtable2)]=apply(ASVtable2[,
2:ncol(ASVtable2)]+1,2,function (x) 100*(x/sum(x)))%>%na.omit()


correlations=function(x){
file2= cor.test(x, 100*(metadata2$MW)/metadata2$Weight, method = "spearman")
file3=data.frame( p=file2$p.value, rho=file2$estimate)
return(file3)
}
#Apply the spearman correlations in parrallel for all ASVs
#Use your sample name order from metadata (DNA/Gene copies concentration)  
#to ensure you have the correct order

spearman_correlations= apply(ASVtable2[,c("ASV",metadata2$Sample)]%>%column_to_rownames(var="ASV"),1, function(x) correlations(x))
#Create the data.frame from the list.
spearman_correlations_df <- as.data.frame(do.call(rbind, spearman_correlations))
#Find the contaminatns
contaminants=filter(spearman_correlations_df, p<0.05&rho<0)
#Check your contaminants
View(contaminants)
rownames(metadata2)<-metadata2$Sample
metadata2$W=metadata2$MW/metadata2$Weight

ASVtable3= ASVtable2%>%column_to_rownames(var="ASV")
#Get the Relative Abundance of contaminants
contaminants_RA= data.frame(contaminants= colSums(ASVtable3[rownames(contaminants),]))
distribution_of_contaminants= ggplot(contaminants_RA, aes(x=contaminants)) + geom_bar(stat = "density")+
  ggtitle(paste0("Mean: ", round(mean(contaminants_RA$contaminants),2), " Standard Deviation: ",
                 round(sd(contaminants_RA$contaminants),2))) + xlab("Contaminants RA (%)")
distribution_of_contaminants

contaminant_free=filter(ASVtable_working, !ASV%in%c(rownames(contaminants)))

ASV_counts=data.frame(Counts=rowSums(contaminant_free[,FS_RS_metadata$Sample]), ASV=contaminant_free$ASV)%>%
  filter(Counts>=10)
contaminant_free2=filter(ASVtable_working, ASV%in%c(ASV_counts$ASV))
contamination_info=data.frame( number_of_removed=nrow(contaminant_free2),
                               number_of_ASVs=nrow(ASVtable_working))
contamination_info2=as.data.frame( cbind(colSums(ASVtable_working[,FS_RS_metadata$Sample]),
      colSums(contaminant_free2[,FS_RS_metadata$Sample])))
colnames(contamination_info2)<-c("Merged", "Contaminant_free")

write.csv(contamination_info, file = "~/number_of_ASVs_contamination_passed_LTE_2020.csv")
write.csv(contamination_info2, file = "~/number_of_reads_contamination_passed_LTE_2020.csv")
write.csv(contaminant_free2, file ="~/plastid_clean_low_read_clean_contaminant_free_ASV_Table_LTE_2020.csv", row.names =T)

