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
library(phangorn) #tool for phylogenetic trees
library(stringr)# string manipulation and text modification
library(igraph) 
library(ggnetwork)
library(lme4) # mixed effect models 
library(lmerTest) # anova like  test for mixed effect models 
library(plyr)
library(ARTool)
library(forcats)
library(tidyverse)
library(ggplot2)
library(dplyr) #
library(plyr) #ggplot expansion for more fancy graphs
library(readxl)

#Diff abundance 2020
set.seed(17011990)
#Load the Cleaned from plastid sequences table


fungi_RA=   read.table(paste0("~/otu_RH.fungi.Ioannis.txt"), sep=" ")%>%t()

fungi_RA=apply(fungi_RA, 1, function(x) 100*(x)/sum(x))%>%as.data.frame()%>%
  rownames_to_column(var="Col")
fungi_RA$OTU=paste0("OTU", row.names(fungi_RA))

OTU_Table_count1=read.table(paste0("~/otu_RH.fungi.Ioannis.txt"), sep=" ")%>%
  rownames_to_column(var="Col") %>%full_join(select(fungi_RA, Col, OTU), by="Col")%>%
  select(-Col)%>%column_to_rownames(var="OTU")%>%t()%>%
  as.data.frame()


fungi_taxonomy=   read.table(paste0("~/ITS2_CT_2020_TAX_reduced.txt"), sep="\t")%>%

  rownames_to_column(var="Col")%>%full_join(fungi_RA%>%select(Col, OTU))%>%filter(
    OTU!="NA"
  )



fungi_taxonomy[fungi_taxonomy=="NA"]="Unclassified"
fungi_taxonomy[fungi_taxonomy=="uncultured"]<-"Unclassified"
fungi_taxonomy[fungi_taxonomy=="metagenome"]<-"Unclassified"

gathered_fungi_table=gather(fungi_RA, -OTU, key = "Sample", value = "Abundance")

#Assign metadata via sample names
metadata_assign=data.frame(Sample=unique(gathered_fungi_table$Sample))%>%
  mutate(Category=Sample)%>% 
  separate(., Category, into=(c("Type", "Practice1", "N-Fertilization_Intensity", "Treatment")), sep = "_")%>% 
  filter(., Practice1=="CT")

metadata_assign$Block=substr(metadata_assign$Treatment,nchar(metadata_assign$Treatment), nchar(metadata_assign$Treatment)+1)
metadata_assign$Treatment2=substr(metadata_assign$Treatment,1, nchar(metadata_assign$Treatment)-1)


metadata_assign$Treatment=factor(metadata_assign$Treatment2, c("Ctrl","BMc"))
RS_metadata=filter(metadata_assign, Type=="Rhizosphere")
fungi_taxonomy2=mutate(fungi_taxonomy%>%select(-Col), Genus=paste0(phylum,";",class,";",order,";",family,";",genus))%>%
  mutate(., Family=paste0(phylum,";",class,";",order,";",family))%>%
  mutate(., Order=paste0(phylum,";",class,";",order))%>%
  mutate(., Class=paste0(phylum,";",class))%>%mutate(., OTU2=paste0(genus,";",OTU))%>%
  mutate(Phylum=phylum)



input1=fungi_RA%>% select(-Col)%>% gather(-OTU, key="Sample", value = "Abundance")%>%
  full_join(fungi_taxonomy2, by="OTU")%>%na.omit()

OTU_Table_RA=dcast(input1, OTU2~Sample, sum, value.var = "Abundance")%>% column_to_rownames(var="OTU2")
OTU_Table=dcast(input1, OTU2~Sample, sum, value.var = "Abundance")%>% column_to_rownames(var="OTU2")%>% apply(.,2, function(x) 100*(x)/sum(x))%>%t()%>%na.omit()
Genera_Table=dcast(input1, Genus~Sample, sum, value.var = "Abundance")%>% column_to_rownames(var="Genus")%>% apply(.,2, function(x) 100*(x)/sum(x))%>%t()%>%na.omit()
Families_Table=dcast(input1, Family~Sample, sum, value.var = "Abundance")%>% column_to_rownames(var="Family")%>%apply(.,2,function(x) 100*(x)/sum(x))%>%t()%>%na.omit()
Order_Table=dcast(input1, Order~Sample, sum, value.var = "Abundance")%>% column_to_rownames(var="Order")%>%apply(.,2,function(x) 100*(x)/sum(x))%>%t()%>%na.omit()
Class_Table=dcast(input1, Class~Sample, sum, value.var = "Abundance")%>% column_to_rownames(var="Class")%>%apply(.,2,function(x) 100*(x)/sum(x))%>%t()%>%na.omit()
Phylum_Table=dcast(input1, Phylum~Sample, sum, value.var = "Abundance")%>% column_to_rownames(var="Phylum")%>%apply(.,2,function(x) 100*(x)/sum(x))%>%t()%>%na.omit()

metadata_assign$Treatment[is.na(metadata_assign$Treatment)]<-"BMc"
metadata_assign$Treatment2=as.numeric(as.factor(metadata_assign$Treatment))-1
Ext_RS_metadata=filter(metadata_assign, `N-Fertilization_Intensity`=="Ext")
Int_RS_metadata=filter(metadata_assign, `N-Fertilization_Intensity`=="Int")
multiple_glm_ext= function(x){
  
Table_glm=glm(Ext_RS_metadata$Treatment2~x, family=binomial(link="logit"))
file2=anova(Table_glm, test = "Chisq")
file3=data.frame(coeff=round(Table_glm$coefficients[2],2), p=file2$`Pr(>Chi)`[2])
return(file3)
}
multiple_glm_int= function(x){
  
  Table_glm=glm(Int_RS_metadata$Treatment2~x, family=binomial(link="logit"))
  file2=anova(Table_glm, test = "Chisq")
  file3=data.frame(coeff=round(Table_glm$coefficients[2],2), p=file2$`Pr(>Chi)`[2])
  return(file3)
}

Ext_multiple_glms_phyla= apply(Phylum_Table[Ext_RS_metadata$Sample,],2, function(x) multiple_glm_ext(x))
Ext_multiple_glms_phyla_df <- as.data.frame(do.call(rbind, Ext_multiple_glms_phyla))
Ext_multiple_glms_phyla_df[is.na(Ext_multiple_glms_phyla_df)]<-1
Ext_multiple_glms_phyla_df$padj=p.adjust(Ext_multiple_glms_phyla_df$p, method="BH")
Ext_multiple_glms_phyla_df=filter(Ext_multiple_glms_phyla_df)%>%mutate(Treatment="Ext")

means=cbind(Ext_RS_metadata%>%select(Treatment), Phylum_Table[Ext_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Ext_RS_metadata%>%select(Treatment), Phylum_Table[Ext_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
phyla_ext=as.data.frame(cbind(Ext_multiple_glms_phyla_df,
t(means[,rownames(Ext_multiple_glms_phyla_df)]),t(sdt[,rownames(Ext_multiple_glms_phyla_df)])))%>%
  rownames_to_column(var="Taxa")%>%mutate(Taxonomic_Level="Phylum", `N-Fertilization-Intensity`="Ext")
  
  


Ext_multiple_glms_class= apply(Class_Table[Ext_RS_metadata$Sample,],2, function(x) multiple_glm_ext(x))
Ext_multiple_glms_class_df <- as.data.frame(do.call(rbind, Ext_multiple_glms_class))
Ext_multiple_glms_class_df[is.na(Ext_multiple_glms_class_df)]<-1
Ext_multiple_glms_class_df$padj=p.adjust(Ext_multiple_glms_class_df$p, method="BH")
Ext_multiple_glms_class_df=filter(Ext_multiple_glms_class_df)%>%mutate(Treatment="Ext")

means=cbind(Ext_RS_metadata%>%select(Treatment), Class_Table[Ext_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Ext_RS_metadata%>%select(Treatment), Class_Table[Ext_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

Class_Ext=as.data.frame(cbind(Ext_multiple_glms_class_df,
                              t(means[,rownames(Ext_multiple_glms_class_df),drop = FALSE]),t(sdt[,rownames(Ext_multiple_glms_class_df),drop = FALSE])))%>%
  mutate(Taxonomic_Level="Class", `N-Fertilization-Extensity`="Ext")
Class_Ext=rownames_to_column(Class_Ext,var="Taxa")


Ext_multiple_glms_order= apply(Order_Table[Ext_RS_metadata$Sample,],2, function(x) multiple_glm_ext(x))
Ext_multiple_glms_order_df <- as.data.frame(do.call(rbind, Ext_multiple_glms_order))
Ext_multiple_glms_order_df[is.na(Ext_multiple_glms_order_df)]<-1
Ext_multiple_glms_order_df$padj=p.adjust(Ext_multiple_glms_order_df$p, method="BH")
Ext_multiple_glms_order_df=filter(Ext_multiple_glms_order_df)%>%mutate(Treatment="Ext")

means=cbind(Ext_RS_metadata%>%select(Treatment), Order_Table[Ext_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Ext_RS_metadata%>%select(Treatment), Order_Table[Ext_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

order_ext=as.data.frame(cbind(Ext_multiple_glms_order_df,
                              t(means[,rownames(Ext_multiple_glms_order_df)]),t(sdt[,rownames(Ext_multiple_glms_order_df)])))%>%
  mutate(Taxonomic_Level="Order", `N-Fertilization-Intensity`="Ext")
order_ext=rownames_to_column(order_ext,var="Taxa")



Ext_multiple_glms_family= apply(Families_Table[Ext_RS_metadata$Sample,],2, function(x) multiple_glm_ext(x))
Ext_multiple_glms_family_df <- as.data.frame(do.call(rbind, Ext_multiple_glms_family))
Ext_multiple_glms_family_df[is.na(Ext_multiple_glms_family_df)]<-1
Ext_multiple_glms_family_df$padj=p.adjust(Ext_multiple_glms_family_df$p, method="BH")
Ext_multiple_glms_family_df=filter(Ext_multiple_glms_family_df)%>%mutate(Treatment="Ext")


means=cbind(Ext_RS_metadata%>%select(Treatment), Families_Table[Ext_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Ext_RS_metadata%>%select(Treatment), Families_Table[Ext_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
Family_ext=as.data.frame(cbind(Ext_multiple_glms_family_df,
                              t(means[,rownames(Ext_multiple_glms_family_df)]),t(sdt[,rownames(Ext_multiple_glms_family_df)])))%>%
  mutate(Taxonomic_Level="Family", `N-Fertilization-Intensity`="Ext")
Family_ext=rownames_to_column(Family_ext,var="Taxa")


Ext_multiple_glms_genera= apply(Genera_Table[Ext_RS_metadata$Sample,],2, function(x) multiple_glm_ext(x))
Ext_multiple_glms_genera_df <- as.data.frame(do.call(rbind, Ext_multiple_glms_genera))
Ext_multiple_glms_genera_df[is.na(Ext_multiple_glms_genera_df)]<-1
Ext_multiple_glms_genera_df$padj=p.adjust(Ext_multiple_glms_genera_df$p, method="BH")
Ext_multiple_glms_genera_df=filter(Ext_multiple_glms_genera_df)%>%mutate(Treatment="Ext")

means=cbind(Ext_RS_metadata%>%select(Treatment), Genera_Table[Ext_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Ext_RS_metadata%>%select(Treatment), Genera_Table[Ext_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

Genera_ext=as.data.frame(cbind(Ext_multiple_glms_genera_df,
                               t(means[,rownames(Ext_multiple_glms_genera_df)]),t(sdt[,rownames(Ext_multiple_glms_genera_df)])))%>%
  mutate(Taxonomic_Level="Genera", `N-Fertilization-Intensity`="Ext")
Genera_ext=rownames_to_column(Genera_ext,var="Taxa")

counts_ext=OTU_Table_count1[Ext_RS_metadata$Sample,]%>% rownames_to_column(var="Sample")%>%
  gather(-Sample, key="OTU", value = "RA")%>% filter(RA>1)%>% group_by( OTU)%>% dplyr::summarise(count=n())%>%
  filter(count>=3)
counts_ext=filter(fungi_taxonomy2, OTU%in%c(counts_ext$OTU))

Ext_multiple_glms_OTU= apply(OTU_Table[Ext_RS_metadata$Sample,counts_ext$OTU2],2, function(x) multiple_glm_ext(x))
Ext_multiple_glms_OTU_df <- as.data.frame(do.call(rbind, Ext_multiple_glms_OTU))
Ext_multiple_glms_OTU_df[is.na(Ext_multiple_glms_OTU_df)]<-1
Ext_multiple_glms_OTU_df$padj=p.adjust(Ext_multiple_glms_OTU_df$p, method="BH")
Ext_multiple_glms_OTU_df=filter(Ext_multiple_glms_OTU_df)%>%mutate(Treatment="Ext")

means=cbind(Ext_RS_metadata%>%select(Treatment), OTU_Table[Ext_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Ext_RS_metadata%>%select(Treatment), OTU_Table[Ext_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

OTU_ext=as.data.frame(cbind(Ext_multiple_glms_OTU_df,
                               t(means[,rownames(Ext_multiple_glms_OTU_df)]),t(sdt[,rownames(Ext_multiple_glms_OTU_df)])))%>%
  mutate(Taxonomic_Level="OTU", `N-Fertilization-Intensity`="Ext")
OTU_ext=rownames_to_column(OTU_ext,var="Taxa")



Int_multiple_glms_phyla= apply(Phylum_Table[Int_RS_metadata$Sample,],2, function(x) multiple_glm_int(x))
Int_multiple_glms_phyla_df <- as.data.frame(do.call(rbind, Int_multiple_glms_phyla))
Int_multiple_glms_phyla_df[is.na(Int_multiple_glms_phyla_df)]<-1
Int_multiple_glms_phyla_df$padj=p.adjust(Int_multiple_glms_phyla_df$p, method="BH")
Int_multiple_glms_phyla_df=filter(Int_multiple_glms_phyla_df)%>%mutate(Treatment="Int")

means=cbind(Int_RS_metadata%>%select(Treatment), Phylum_Table[Int_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Int_RS_metadata%>%select(Treatment), Phylum_Table[Int_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

Phylum_Int=as.data.frame(cbind(Int_multiple_glms_phyla_df,
                               t(means[,rownames(Int_multiple_glms_phyla_df)]),t(sdt[,rownames(Int_multiple_glms_phyla_df)])))%>%
  mutate(Taxonomic_Level="Phylum", `N-Fertilization-Intensity`="Int")
Phylum_Int=rownames_to_column(Phylum_Int,var="Taxa")



Int_multiple_glms_class= apply(Class_Table[Int_RS_metadata$Sample,],2, function(x) multiple_glm_int(x))
Int_multiple_glms_class_df <- as.data.frame(do.call(rbind, Int_multiple_glms_class))
Int_multiple_glms_class_df[is.na(Int_multiple_glms_class_df)]<-1
Int_multiple_glms_class_df$padj=p.adjust(Int_multiple_glms_class_df$p, method="BH")
Int_multiple_glms_class_df=filter(Int_multiple_glms_class_df)%>%mutate(Treatment="Int")

means=cbind(Int_RS_metadata%>%select(Treatment), Class_Table[Int_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Int_RS_metadata%>%select(Treatment), Class_Table[Int_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

Class_Int=as.data.frame(cbind(Int_multiple_glms_class_df,
                               t(means[,rownames(Int_multiple_glms_class_df)]),t(sdt[,rownames(Int_multiple_glms_class_df)])))%>%
  mutate(Taxonomic_Level="Class", `N-Fertilization-Intensity`="Int")
Class_Int=rownames_to_column(Class_Int,var="Taxa")



Int_multiple_glms_order= apply(Order_Table[Int_RS_metadata$Sample,],2, function(x) multiple_glm_int(x))
Int_multiple_glms_order_df <- as.data.frame(do.call(rbind, Int_multiple_glms_order))
Int_multiple_glms_order_df[is.na(Int_multiple_glms_order_df)]<-1
Int_multiple_glms_order_df$padj=p.adjust(Int_multiple_glms_order_df$p, method="BH")
Int_multiple_glms_order_df=filter(Int_multiple_glms_order_df)%>%mutate(Treatment="Int")

means=cbind(Int_RS_metadata%>%select(Treatment), Order_Table[Int_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Int_RS_metadata%>%select(Treatment), Order_Table[Int_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

Order_Int=as.data.frame(cbind(Int_multiple_glms_order_df,
                               t(means[,rownames(Int_multiple_glms_order_df)]),t(sdt[,rownames(Int_multiple_glms_order_df)])))%>%
  mutate(Taxonomic_Level="Order", `N-Fertilization-Intensity`="Int")
Order_Int=rownames_to_column(Order_Int,var="Taxa")


Int_multiple_glms_family= apply(Families_Table[Int_RS_metadata$Sample,],2, function(x) multiple_glm_int(x))
Int_multiple_glms_family_df <- as.data.frame(do.call(rbind, Int_multiple_glms_family))
Int_multiple_glms_family_df[is.na(Int_multiple_glms_family_df)]<-1
Int_multiple_glms_family_df$padj=p.adjust(Int_multiple_glms_family_df$p, method="BH")
Int_multiple_glms_family_df=filter(Int_multiple_glms_family_df)%>%mutate(Treatment="Int")

means=cbind(Int_RS_metadata%>%select(Treatment), Families_Table[Int_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Int_RS_metadata%>%select(Treatment), Families_Table[Int_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

Family_Int=as.data.frame(cbind(Int_multiple_glms_family_df,
                               t(means[,rownames(Int_multiple_glms_family_df)]),t(sdt[,rownames(Int_multiple_glms_family_df)])))%>%
  mutate(Taxonomic_Level="Family", `N-Fertilization-Intensity`="Int")
Family_Int=rownames_to_column(Family_Int,var="Taxa")


Int_multiple_glms_genera= apply(Genera_Table[Int_RS_metadata$Sample,],2, function(x) multiple_glm_int(x))
Int_multiple_glms_genera_df <- as.data.frame(do.call(rbind, Int_multiple_glms_genera))
Int_multiple_glms_genera_df[is.na(Int_multiple_glms_genera_df)]<-1
Int_multiple_glms_genera_df$padj=p.adjust(Int_multiple_glms_genera_df$p, method="BH")
Int_multiple_glms_genera_df=filter(Int_multiple_glms_genera_df)%>%mutate(Treatment="Int")

means=cbind(Int_RS_metadata%>%select(Treatment), Genera_Table[Int_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Int_RS_metadata%>%select(Treatment), Genera_Table[Int_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

Genera_Int=as.data.frame(cbind(Int_multiple_glms_genera_df,
                               t(means[,rownames(Int_multiple_glms_genera_df)]),t(sdt[,rownames(Int_multiple_glms_genera_df)])))%>%
  mutate(Taxonomic_Level="Genera", `N-Fertilization-Intensity`="Int")
Genera_Int=rownames_to_column(Genera_Int,var="Taxa")

counts_int=OTU_Table_count1[Int_RS_metadata$Sample,]%>% as.data.frame()%>% rownames_to_column(var="Sample")%>%
gather(-Sample, key="OTU", value = "RA")%>% filter(RA>1)%>% group_by(OTU)%>%dplyr:: summarise(count=n())%>%
  filter(count>=3)
counts_int=filter(fungi_taxonomy2, OTU%in%c(counts_int$OTU))

Int_multiple_glms_OTU= apply(OTU_Table[Int_RS_metadata$Sample,(counts_int$OTU2)],2, function(x) multiple_glm_int(x))
Int_multiple_glms_OTU_df <- as.data.frame(do.call(rbind, Int_multiple_glms_OTU))
Int_multiple_glms_OTU_df[is.na(Int_multiple_glms_OTU_df)]<-1
Int_multiple_glms_OTU_df$padj=p.adjust(Int_multiple_glms_OTU_df$p, method="BH")
Int_multiple_glms_OTU_df=filter(Int_multiple_glms_OTU_df)%>%mutate(Treatment="Int")

means=cbind(Int_RS_metadata%>%select(Treatment), OTU_Table[Int_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Int_RS_metadata%>%select(Treatment), OTU_Table[Int_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

OTU_Int=as.data.frame(cbind(Int_multiple_glms_OTU_df,
                               t(means[,rownames(Int_multiple_glms_OTU_df)]),t(sdt[,rownames(Int_multiple_glms_OTU_df)])))%>%
  mutate(Taxonomic_Level="OTU", `N-Fertilization-Intensity`="Int")
OTU_Int=rownames_to_column(OTU_Int,var="Taxa")


colnames(Class_Ext)<-colnames(Class_Int)
complete_info_tax_levels= rbind(phyla_ext,Phylum_Int,Class_Ext, Class_Int,order_ext,Order_Int,Family_ext,Family_Int, Genera_Int,Genera_ext, OTU_Int,OTU_ext)

write_csv(complete_info_tax_levels%>%filter(.,padj<0.05), file = "~/Fungal_Responders_over_taxonomic_levels_MS3.csv")
write_csv(rbind(OTU_Int,OTU_ext)%>%filter(.,padj<0.05), file = "~/Fungal_Responders_OTU_levels_MS3.csv")

input_plots=complete_info_tax_levels%>%filter(.,padj<0.05)
input_plots$practice=input_plots$`N-Fertilization-Intensity`
input_plots$Inoculation=input_plots$coeff
input_plots$Inoculation[input_plots$Inoculation>0]<-"BMc"
input_plots$Inoculation[!input_plots$Inoculation=="BMc"]<-"Ctrl"

two_managements_phyla= ggplot(input_plots%>%filter(Taxonomic_Level=="Phylum"), aes(y=fct_reorder(Taxa,log2((Mean_BMc+1)/(Mean_Ctrl+1))),
x=log2((Mean_BMc+1)/(Mean_Ctrl+1)), shape=practice, colour=Inoculation)) +
  geom_point(size=12) +
  scale_colour_manual(name="Inoculation",values=c( "#ff9900","#146eb4"))+
  theme_pubclean() + 
  geom_vline(xintercept=0, linetype="dashed") + scale_shape_manual(name="", values = c(16,17))+
  
  theme( axis.text.y = element_text( size=25, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+xlab("log2FC(BMc/Ctrl)")+ylab("")+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+scale_x_reverse()+ 
  guides(fill = guide_legend(override.aes = list(shape = 21)))
two_managements_phyla

two_managements_classes= ggplot(input_plots%>%filter(Taxonomic_Level=="Class"), aes(y=fct_reorder(Taxa, log2((Mean_BMc+1)/(Mean_Ctrl+1))),
x=log2((Mean_BMc+1)/(Mean_Ctrl+1)), shape=practice, colour=Inoculation)) +
  geom_point(size=12) +
  scale_colour_manual(name="Inoculation",values=c( "#ff9900","#146eb4"))+
  theme_pubclean() + 
  geom_vline(xintercept=0, linetype="dashed") + scale_shape_manual(name="", values = c(16,17))+
  
  theme( axis.text.y = element_text( size=25, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+xlab("log2FC(BMc/Ctrl)")+ylab("")+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+scale_x_reverse()+ 
  guides(fill = guide_legend(override.aes = list(shape = 21)))
two_managements_classes



two_managements_order= ggplot(input_plots%>%filter(Taxonomic_Level=="Order"), aes(y=fct_reorder(Taxa,log2((Mean_BMc+1)/(Mean_Ctrl+1))),x=log2((Mean_BMc+1)/(Mean_Ctrl+1)), shape=practice, colour=Inoculation)) +
  geom_point(size=12) +
  scale_colour_manual(name="Inoculation",values=c( "#ff9900","#146eb4"))+
  theme_pubclean() + 
  geom_vline(xintercept=0, linetype="dashed") + scale_shape_manual(name="", values = c(16,17))+
  
  theme( axis.text.y = element_text( size=25, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+xlab("log2FC(BMc/Ctrl)")+ylab("")+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+scale_x_reverse()+ 
  guides(fill = guide_legend(override.aes = list(shape = 21)))
two_managements_order


two_managements_families= ggplot(input_plots%>%filter(Taxonomic_Level=="Family"), aes(y=fct_reorder(Taxa,log2((Mean_BMc+1)/(Mean_Ctrl+1))),x=log2((Mean_BMc+1)/(Mean_Ctrl+1)), shape=practice, colour=Inoculation)) +
  geom_point(size=12) +
  scale_colour_manual(name="Inoculation",values=c( "#ff9900","#146eb4"))+
  theme_pubclean() + 
  geom_vline(xintercept=0, linetype="dashed") + scale_shape_manual(name="", values = c(16,17))+
  
  theme( axis.text.y = element_text( size=15, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+xlab("log2FC(BMc/Ctrl)")+ylab("")+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+scale_x_reverse()+ 
  guides(fill = guide_legend(override.aes = list(shape = 21)))
two_managements_families

two_managements_genera= ggplot(input_plots%>%filter(Taxonomic_Level=="Genera"), aes(y=fct_reorder(Taxa,log2((Mean_BMc+1)/(Mean_Ctrl+1))),x=log2((Mean_BMc+1)/(Mean_Ctrl+1)), shape=practice, colour=Inoculation)) +
  geom_point(size=12) +
  scale_colour_manual(name="Inoculation",values=c( "#ff9900","#146eb4"))+
  theme_pubclean() + 
  geom_vline(xintercept=0, linetype="dashed") + scale_shape_manual(name="", values = c(16,17))+
  
  theme( axis.text.y = element_text( size=10, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+xlab("log2FC(BMc/Ctrl)")+ylab("")+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+scale_x_reverse()+ 
  guides(fill = guide_legend(override.aes = list(shape = 21)))
two_managements_genera

top15Genera= fungi_RA%>%
  as.data.frame()%>%
  gather(-OTU, -Col, key="Sample", value = "RA")%>%
  full_join(fungi_taxonomy2, by="OTU") %>%
  mutate(OTU=paste0(phylum,";",class,";",order,";",family, ";",genus)) %>% 
  full_join(metadata_assign, by="Sample")%>%
  dcast(OTU+Sample~., value.var = "RA", sum)%>%dplyr::select(., RA=".",everything())%>%
  dcast(OTU~., value.var = "RA", median)
top15Genera1=top15Genera[order(top15Genera$., decreasing= TRUE), ]
top15Genera2=top15Genera1[1:15,]

top15Genera_RA= fungi_RA%>%
  as.data.frame()%>%
  gather(-OTU, -Col, key="Sample", value = "RA")%>%
  full_join(fungi_taxonomy2, by="OTU") %>%
  mutate(OTU=paste0(Genus)) %>% filter(OTU%in%c((top15Genera2$OTU)))%>% 
  full_join(metadata_assign, by="Sample")%>%
  dcast(OTU+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%dplyr::select(RA=".", everything())%>%
  dcast(OTU+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  mean)%>%dplyr::select(mean=".", everything())%>%
  filter(OTU%in%top15Genera$OTU)
top15Genera_RA=fungi_RA%>%
  as.data.frame()%>%
  gather(-OTU, -Col, key="Sample", value = "RA")%>%
  full_join(fungi_taxonomy2, by="OTU") %>%
  mutate(OTU=paste0(Genus)) %>% filter(OTU%in%c((top15Genera2$OTU)))%>% 
  full_join(metadata_assign, by="Sample")%>%
  dcast(OTU+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%dplyr::select(RA=".", everything())%>%
  dcast(OTU+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  sd)%>%dplyr::select(sdt=".", everything())%>%full_join(top15Genera_RA)%>%
  filter(OTU%in%top15Genera$OTU)



top15Genera_RA_sign_stats=filter( complete_info_tax_levels, Taxa%in%c( top15Genera2$OTU))%>%select(OTU=Taxa, padj, coeff,  `N-Fertilization_Intensity`=`N-Fertilization-Intensity`)%>%
  full_join(top15Genera_RA)%>%mutate( `p-value_corrected_(BH)`=padj)
top15Genera_RA_sign_stats$padj[top15Genera_RA_sign_stats$padj<0.05]<-"p<0.05"
top15Genera_RA_sign_stats$padj[!top15Genera_RA_sign_stats$padj=="p<0.05"]<-"p>0.05"
top15Genera_RA_sign_stats$mean=round(top15Genera_RA_sign_stats$mean,3)

library(ggplot2)

top15Genera_RA_sign_stats_plot=ggplot(top15Genera_RA_sign_stats%>%mutate(`N-Fertilization_Intensity`=str_replace(`N-Fertilization_Intensity`, "ensive","")), aes(y=fct_reorder(OTU,mean), x=paste0( Treatment), fill=mean, shape=padj)) +facet_wrap(~paste(`N-Fertilization_Intensity`)) + geom_point(size=25)+
  scale_fill_distiller(palette = "Spectral", name="Relative_Abundance (%)") + scale_shape_manual(values = c(21, 24), name="p-value (BMc vs Ctrl)") +#scale_shape_manual(values = c(21,24), "") +
  theme_pubr(legend = "right") + 
  theme( axis.text.y = element_text( size=25, face = "bold", colour="black")) +
  theme( strip.text = element_text( size=45, face = "bold", colour="white")) +
  theme(strip.background =  element_rect( fill="darkblue")) +
  
  theme(axis.text.x = element_text(size = 35, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
  theme(legend.key.size = unit(2,"cm")) +xlab("")+ylab("")+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+geom_text(aes(label=round(mean,2)),size=8)

png("~/top15_Fungal_Genera_RA_plot_with_log_regression.png", width = 2500, height=1500)
print(top15Genera_RA_sign_stats_plot) 
dev.off()

write_excel_csv(top15Genera_RA_sign_stats, file = "~/top15_genera_RA_plot_with_log_regression.csv")


top15Classes=fungi_RA %>%
  as.data.frame()%>%gather(-OTU,-Col, key="Sample", value = "RA")%>%
  full_join(fungi_taxonomy2, by="OTU") %>%
  mutate(OTU=Class) %>% full_join(metadata_assign, by="Sample")%>%na.omit()%>%
  dcast(OTU+Sample~., value.var = "RA", sum)%>%dplyr::select(., RA=".",everything())%>% 
  dcast(OTU~., value.var = "RA", median)
top15Classes1=top15Classes[order(top15Classes$., decreasing= TRUE), ]
top15Classes2=top15Classes1[1:15,]

top15Classes_RA= fungi_RA%>%
  as.data.frame()%>%
  gather(-OTU,-Col, key="Sample", value = "RA")%>%
  full_join(fungi_taxonomy2, by="OTU") %>%
  mutate(OTU=Class) %>% filter(OTU%in%c((top15Classes2$OTU)))%>% 
  full_join(metadata_assign, by="Sample")%>%na.omit()%>%
  dcast(OTU+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%dplyr::select(RA=".", everything())%>%
  dcast(OTU+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  mean)%>%dplyr::select(mean=".", everything())



top15Classes_RA= fungi_RA%>%
  as.data.frame()%>% gather(-OTU,-Col, key="Sample", value = "RA")%>%
  full_join(fungi_taxonomy2, by="OTU") %>%
  mutate(OTU=paste0(Class)) %>% filter(OTU%in%c((top15Classes2$OTU)))%>% 
  full_join(metadata_assign, by="Sample")%>%na.omit()%>%
  dcast(OTU+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%dplyr::select(RA=".", everything())%>%
  dcast(OTU+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  sd)%>%dplyr::select(sdt=".", everything())%>% full_join(top15Classes_RA)


top15Classes_RA_sign_stats=filter( complete_info_tax_levels, Taxa%in%c( top15Classes2$OTU))%>%select(OTU=Taxa, padj, coeff,  `N-Fertilization_Intensity`=`N-Fertilization-Intensity`)%>%
  full_join(top15Classes_RA)%>%mutate( `p-value_corrected_(BH)`=padj)%>%filter(OTU!="NA")

top15Classes_RA_sign_stats$padj[top15Classes_RA_sign_stats$padj<0.05]<-"p<0.05"
top15Classes_RA_sign_stats$padj[!top15Classes_RA_sign_stats$padj=="p<0.05"]<-"p>0.05"


top15Classes_RA_sign_stats$mean=round(top15Classes_RA_sign_stats$mean,3)
top15Classes_RA_sign_stats_plot=ggplot(top15Classes_RA_sign_stats%>%mutate(`N-Fertilization_Intensity`=str_replace(`N-Fertilization_Intensity`, "ensive",".")), aes(y=fct_reorder(OTU,mean), x=paste0( Treatment), fill=mean, shape=padj)) +facet_wrap(~paste(`N-Fertilization_Intensity`)) + geom_point(size=25)+
  scale_fill_distiller(palette = "Spectral", name="Relative_Abundance (%)") + scale_shape_manual(values = c(21, 24), name="p-value (BMc vs Ctrl)") +#scale_shape_manual(values = c(21,24), "") +
  theme_pubr(legend = "right") + 
  theme( axis.text.y = element_text( size=25, face = "bold", colour="black")) +
  theme( strip.text = element_text( size=45, face = "bold", colour="white")) +
  theme(strip.background =  element_rect( fill="darkblue")) +
  
  theme(axis.text.x = element_text(size = 35, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
  theme(legend.key.size = unit(2,"cm")) +xlab("")+ylab("")+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
size = 2, linetype = "solid"))+geom_text(aes(label=round(mean,2)),size=8)
top15Classes_RA_sign_stats_plot

png("~/top15_fungal_Classes_RA_plot_with_log_regression.png", width = 1600, height=1500)
print(top15Classes_RA_sign_stats_plot) 
dev.off()
write_excel_csv(top15Classes_RA_sign_stats, file = "~/top15_fungal_classes_RA_plot_with_log_regression.csv")


top15Phyla= fungi_RA%>%
  as.data.frame()%>% 
  gather(-OTU,-Col, key="Sample", value = "RA")%>%
  full_join(fungi_taxonomy2, by="OTU") %>%
  mutate(OTU=paste0(Phylum)) %>% full_join(metadata_assign, by="Sample")%>%
  dcast(OTU+Sample~., value.var = "RA", sum)%>%dplyr::select(., RA=".",everything())%>%
  na.omit()%>%
  dcast(OTU~., value.var = "RA", median)
top15Phyla1=top15Phyla[order(top15Phyla$., decreasing= TRUE), ]
top15Phyla2=top15Phyla1[1:15,]

top15Phyla_RA= fungi_RA%>%
  as.data.frame()%>% 
  gather(-OTU,-Col, key="Sample", value = "RA")%>%
  full_join(fungi_taxonomy2, by="OTU") %>%
  mutate(OTU=paste0(Phylum)) %>% filter(OTU%in%c((top15Phyla2$OTU)))%>% 
  full_join(metadata_assign, by="Sample")%>%
  dcast(OTU+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%
  dplyr::select(RA=".", everything())%>%na.omit()%>%
  dcast(OTU+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  mean)%>%dplyr::select(mean=".", everything())

top15Phyla_RA= fungi_RA%>%
  as.data.frame()%>%
  gather(-OTU,-Col, key="Sample", value = "RA")%>%
  full_join(fungi_taxonomy2, by="OTU") %>%
  mutate(OTU=paste0(Phylum)) %>% filter(OTU%in%c((top15Phyla2$OTU)))%>% 
  full_join(metadata_assign, by="Sample")%>%
  dcast(OTU+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%dplyr::select(RA=".", everything())%>%na.omit() %>%
  dcast(OTU+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  sd)%>%dplyr::select(sdt=".", everything())%>%full_join(top15Phyla_RA)

top15Phyla_RA_sign_stats=filter( complete_info_tax_levels, Taxa%in%c( top15Phyla2$OTU))%>%select(OTU=Taxa, padj, coeff,  `N-Fertilization_Intensity`=`N-Fertilization-Intensity`)%>%
  full_join(top15Phyla_RA)%>%mutate( `p-value_corrected_(BH)`=padj)
top15Phyla_RA_fungi=top15Phyla_RA_sign_stats


top15Phyla_RA_sign_stats$padj[top15Phyla_RA_sign_stats$padj<0.05]<-"p<0.05"
top15Phyla_RA_sign_stats$padj[!top15Phyla_RA_sign_stats$padj=="p<0.05"]<-"p>0.05"


library(ggplot2)
top15Phyla_RA_sign_stats$mean=round(top15Phyla_RA_sign_stats$mean,3)
top15Phyla_RA_sign_stats_plot=ggplot(top15Phyla_RA_sign_stats%>%
mutate(`N-Fertilization_Intensity`=str_replace(`N-Fertilization_Intensity`, "ensive","."))%>%select(padj, everything()), aes(y=fct_reorder(OTU,mean),
x=paste0( Treatment), fill=mean, shape=padj)) +facet_wrap(~paste(`N-Fertilization_Intensity`)) + geom_point(size=25)+
  scale_fill_distiller(palette = "Spectral", name="Relative_Abundance (%)") + scale_shape_manual(values = c(21, 24), name="p-value (BMc vs Ctrl)") +#scale_shape_manual(values = c(21,24), "") +
  theme_pubr(legend = "right") + 
  theme( axis.text.y = element_text( size=25, face = "bold", colour="black")) +
  theme( strip.text = element_text( size=45, face = "bold", colour="white")) +
  theme(strip.background =  element_rect( fill="darkblue")) +
  
  theme(axis.text.x = element_text(size = 35, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
  theme(legend.key.size = unit(2,"cm")) +xlab("")+ylab("")+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+
  theme(legend.position="bottom", legend.direction = "horizontal",
        legend.box = "vertical",legend.justification = "center", legend.key.width =  unit(2.3,"cm")) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE)) +
  theme( axis.line = element_line(colour = "black", 
                                  size = 2, linetype = "solid"))+
  geom_text(aes(label=round(mean,2)),size=8)

png("~/top15_fungal_Phyla_RA_plot_with_log_regression.png", width = 1300, height=1500)
print(top15Phyla_RA_sign_stats_plot) 
dev.off()
write_excel_csv(top15Phyla_RA_sign_stats, file = "~/top15_fungal_phyla_RA_plot_with_log_regression.csv")

combined_ASVs=rbind(OTU_ext, OTU_Int)%>%filter(padj<0.05)%>%mutate(logFC=log((10^3*Mean_BMc+1)/(10^3*Mean_Ctrl+1)))%>%
  mutate(Inoculation=logFC)
combined_ASVs$Inoculation[combined_ASVs$Inoculation>0]<-"BMc"
combined_ASVs$Inoculation[!combined_ASVs$Inoculation=="BMc"]<-"Ctrl"

two_managements= ggplot(combined_ASVs,
aes(y=fct_reorder(Taxa,logFC),x=logFC, shape=`N-Fertilization-Intensity`, colour=Inoculation)) +
  geom_point(size=12) +
  scale_colour_manual(name="Inoculation",values=c( "#ff9900","#146eb4"))+
  theme_pubclean() + 
  geom_vline(xintercept=0, linetype="dashed") + scale_shape_manual(name="", values = c(16,17))+
  
  theme( axis.text.y = element_text( size=25, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+xlab("log2FC(BMc/Ctrl)")+ylab("")+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+scale_x_reverse()+ 
  guides(fill = guide_legend(override.aes = list(shape = 21)))



png("~/combined_fungal_log_plot_with_log.png", width = 2000, height=8000)
print(two_managements) 
dev.off()



top15OTU= fungi_RA%>%
  as.data.frame()%>%
  gather(-OTU,-Col, key="Sample", value = "RA")%>%
  full_join(fungi_taxonomy2, by="OTU") %>%
  mutate(OTU=paste0(phylum,";",class,";",order,";", family,";", genus,";", OTU)) %>% full_join(metadata_assign, by="Sample")%>%
  dcast(OTU2+Sample~., value.var = "RA", sum)%>%dplyr::select(., RA=".",everything())%>%na.omit()%>% 
  dcast(OTU2~., value.var = "RA", median)
top15OTU1=top15OTU[order(top15OTU$., decreasing= TRUE), ]
top15OTU2=top15OTU1[1:15,]

top15OTU_RA= fungi_RA %>%
  as.data.frame()%>% 
  gather(-OTU,-Col, key="Sample", value = "RA")%>%
  full_join(fungi_taxonomy2, by="OTU") %>%
  mutate(OTU=paste0(phylum,";",class,";",order,";", family,";", genus,";", OTU)) %>%
  filter(OTU2%in%c((top15OTU2$OTU2)))%>% 
  full_join(metadata_assign, by="Sample")%>%
  dcast(OTU+OTU2+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%dplyr::select(RA=".", everything())%>%
  na.omit()%>%
  dcast(OTU+OTU2+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  mean)%>%dplyr::select(mean=".", everything())

top15OTU_RA= fungi_RA%>%
  as.data.frame()%>%   gather(-OTU,-Col, key="Sample", value = "RA")%>%
  full_join(fungi_taxonomy2, by="OTU") %>%
  mutate(OTU=paste0(phylum,";",class,";",order,";", family,";", genus,";", OTU)) %>%
  filter(OTU2%in%c((top15OTU2$OTU2)))%>% 
  full_join(metadata_assign, by="Sample")%>%
  dcast(OTU+OTU2+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%
  dplyr::select(RA=".", everything())%>%
  na.omit()%>%
  dcast(OTU+OTU2+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  sd)%>%
  dplyr::select(sdt=".", everything())%>%full_join(top15OTU_RA)


top15OTU_RA_sign_stats=filter( rbind(OTU_ext,OTU_Int), 
Taxa%in%c( top15OTU2$OTU2))%>%select(OTU2=Taxa, padj, coeff,  `N-Fertilization_Intensity`=`N-Fertilization-Intensity`)%>%
  full_join(top15OTU_RA)%>%mutate( `p-value_corrected_(BH)`=padj)
top15OTU_RA_sign_stats$padj[top15OTU_RA_sign_stats$padj<0.05]<-"p<0.05"
top15OTU_RA_sign_stats$padj[!top15OTU_RA_sign_stats$padj=="p<0.05"]<-"p>0.05"

library(ggplot2)
top15OTU_RA_sign_stats$mean=round(top15OTU_RA_sign_stats$mean,3)
top15OTU_RA_sign_stats_plot=ggplot(top15OTU_RA_sign_stats%>%mutate(`N-Fertilization_Intensity`=str_replace(`N-Fertilization_Intensity`, "ensive",".")), 
aes(y=fct_reorder(OTU2,mean), x=paste0( Treatment), fill=mean, shape=padj)) +facet_wrap(~paste(`N-Fertilization_Intensity`)) +
  geom_point(size=26)+
  scale_fill_distiller(palette = "Spectral", name="Relative_Abundance (%)") + 
  scale_shape_manual(values = c(21, 24), name="p-value (BMc vs Ctrl)") +#scale_shape_manual(values = c(21,24), "") +
  theme_pubr(legend = "right") + 
  theme( axis.text.y = element_text( size=25, face = "bold", colour="black")) +
  theme( strip.text = element_text( size=45, face = "bold", colour="white")) +
  theme(strip.background =  element_rect( fill="darkblue")) +
  
  theme(axis.text.x = element_text(size = 35, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
  theme(legend.key.size = unit(2,"cm")) +xlab("")+ylab("")+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+
  geom_text(aes(label=round(mean,2)),size=8) +
  theme(legend.position="bottom", legend.direction = "horizontal",
        legend.box = "vertical",legend.justification = "center", legend.key.width =  unit(1.8,"cm")) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE)) +
  theme( axis.line = element_line(colour = "black", 
                                  size = 2, linetype = "solid"))

png("~/top15_fungal_OTUs_RA_plot_with_log_regression.png", width = 1500, height=1500)
print(top15OTU_RA_sign_stats_plot) 
dev.off()
write_excel_csv(top15OTU_RA_sign_stats, file = "~/top15_fungal_OTUs_RA_plot_with_log_regression.csv")

