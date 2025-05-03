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
    ASVtable_working <- read.csv("~/plastid_clean_low_read_clean_contaminant_free_ASV_Table_LTE_2020.csv", row.names = 1)
    #Function to calculate median and prevalence, so we avoid the 
    # confirmed DNA extraction kit contaminants such as Ralstonia pickettii
    
    #Load the Taxonomy Table of the Rarefied ASV table. 
    nrow(ASVtable_working%>%
           filter(., str_detect(Kingdom, "Archaea")))
    #No Archaea 
    ASVtax_table <- ASVtable_working[,2:9] %>%as.data.frame()
    ASVtax_table <- ASVtax_table %>% mutate(Taxonomy=paste0(Kingdom,";",Phylum,";",Class,";",Order,";",Family,";",Genus))%>%
      mutate(Taxonomy =gsub("[[:space:]]", "_", Taxonomy))%>%
      mutate(Taxonomy=str_replace_all(Taxonomy, ".ncultured.*.", "Unclassified"))%>%
      mutate(Taxonomy=str_replace_all(Taxonomy, "uncultured*", "Unclassified"))%>%
      mutate(Taxonomy=str_replace(Taxonomy, "metagenome", "Unclassified"))
    
ASVtax_table <- dplyr::select(ASVtax_table,ASV,Sequence ,Taxonomy)%>%separate(Taxonomy,
                                                                                  into = c("Kingdom","Phylum","Class","Order","Family","Genus"), sep=";")
ASVtax_table[ASVtax_table=="NA"]="Unclassified"
ASVtax_table[ASVtax_table=="uncultured"]<-"Unclassified"
ASVtax_table[ASVtax_table=="metagenome"]<-"Unclassified"
rownames(ASVtax_table)=ASVtax_table$ASV
gathered_ASV_table=gather(ASVtable_working[c(1,3:114)], -ASV, key = "Sample", value = "Abundance")

ASVtax_table[is.na(ASVtax_table)]<-"Unclassified"
for(i in 1:nrow(ASVtax_table)){
  gn =paste0(ASVtax_table[i,8]) 
  if (ASVtax_table[i,8]=="Unclassified"){
    gn =paste0("Unclassified_", ASVtax_table[i,7]) 
    if (ASVtax_table[i,7]=="Unclassified"){
      gn =paste0("Unclassified_", ASVtax_table[i,6]) 
      if (ASVtax_table[i,6]=="Unclassified"){
        gn =paste0("Unclassified_", ASVtax_table[i,5]) 
        if (ASVtax_table[i,5]=="Unclassified"){
          gn =paste0("Unclassified_", ASVtax_table[i,4]) 
          
          if (ASVtax_table[i,4]=="Unclassified"){
            gn =paste0("Unclassified_", ASVtax_table[i,3]) 
            
          }
        }
      }
    }
  }
  ASVtax_table[i,8]=paste0(gn)
  print(paste0("Round ",i))
}

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

metadata_assign$Treatment=factor(metadata_assign$Treatment, c("Ctrl","BMc"))
FS_RS_metada=filter(metadata_assign, Type=="Field-Soil"|Type=="Rhizosphere"&Sample!="FS.CT.Ext.C4..2020.trimmed")
ASVtable=dplyr::select(ASVtable_working, ASV, c(FS_RS_metada$Sample))%>%na.omit()%>%column_to_rownames(var="ASV")
FS_metada=filter(metadata_assign, Type=="Field-Soil"&Sample!="FS.CT.Ext.C4..2020.trimmed")
ASVtable=dplyr::select(ASVtable_working, ASV, c(FS_RS_metada$Sample))%>%na.omit()%>%column_to_rownames(var="ASV")
#Exclude ASVs with less than five reads
ASVs=ASVtable%>%dplyr::select(., -`FS.CT.Ext.C4..2020.trimmed`)
#Relative abundance transformation
ASV_100=apply(ASVs, 2, function (x)(x/sum(x))*100)
metadata_assign$Treatment2 = as.numeric(metadata_assign$Treatment)-1
#Differential Abundance##################
#First Analysis of Rhizosphere samples
RS_metadata=filter(metadata_assign, Type=="Rhizosphere")
ASVtax_table2=mutate(ASVtax_table, Genus=paste0(Phylum,";",Class,";",Order,";",Family,";",Genus))%>%
  mutate(., Family=paste0(Phylum,";",Class,";",Order,";",Family))%>%
  mutate(., Order=paste0(Phylum,";",Class,";",Order))%>%
  mutate(., Class=paste0(Phylum,";",Class))%>%mutate(., ASV2=paste0(Genus,";",ASV))



input1=ASVs%>%rownames_to_column(var="ASV")%>%
  gather(-ASV, key="Sample", value = "Abundance")%>%
  full_join(ASVtax_table2, by="ASV")

ASV_Table_count=dcast(input1, ASV2~Sample, sum, value.var = "Abundance")%>% column_to_rownames(var="ASV2")
ASV_Table=dcast(input1, ASV2~Sample, sum, value.var = "Abundance")%>% column_to_rownames(var="ASV2")%>% apply(.,2, function(x) 100*(x)/sum(x))%>%t()
Genera_Table=dcast(input1, Genus~Sample, sum, value.var = "Abundance")%>% column_to_rownames(var="Genus")%>% apply(.,2, function(x) 100*(x)/sum(x))%>%t()
Families_Table=dcast(input1, Family~Sample, sum, value.var = "Abundance")%>% column_to_rownames(var="Family")%>%apply(.,2,function(x) 100*(x)/sum(x))%>%t()
Order_Table=dcast(input1, Order~Sample, sum, value.var = "Abundance")%>% column_to_rownames(var="Order")%>%apply(.,2,function(x) 100*(x)/sum(x))%>%t()
Class_Table=dcast(input1, Class~Sample, sum, value.var = "Abundance")%>% column_to_rownames(var="Class")%>%apply(.,2,function(x) 100*(x)/sum(x))%>%t()
Phylum_Table=dcast(input1, Phylum~Sample, sum, value.var = "Abundance")%>% column_to_rownames(var="Phylum")%>%apply(.,2,function(x) 100*(x)/sum(x))%>%t()


Extensive_RS_metadata=filter(metadata_assign, Type=="Rhizosphere"& 
                               `N-Fertilization_Intensity`=="Extensive")
Intensive_RS_metadata=filter(metadata_assign, Type=="Rhizosphere"& 
                               `N-Fertilization_Intensity`=="Intensive")
metadata_assign$Treatment2 = as.numeric(metadata_assign$Treatment)-1

multiple_glm_ext= function(x){
  
Table_glm=glm(Extensive_RS_metadata$Treatment2~x, family=binomial(link="logit"))
file2=anova(Table_glm, test = "Chisq")
file3=data.frame(coeff=round(Table_glm$coefficients[2],2), p=file2$`Pr(>Chi)`[2])
return(file3)
}
multiple_glm_int= function(x){
  
  Table_glm=glm(Intensive_RS_metadata$Treatment2~x, family=binomial(link="logit"))
  file2=anova(Table_glm, test = "Chisq")
  file3=data.frame(coeff=round(Table_glm$coefficients[2],2), p=file2$`Pr(>Chi)`[2])
  return(file3)
}

extensive_multiple_glms_phyla= apply(Phylum_Table[Extensive_RS_metadata$Sample,],2, function(x) multiple_glm_ext(x))
extensive_multiple_glms_phyla_df <- as.data.frame(do.call(rbind, extensive_multiple_glms_phyla))
extensive_multiple_glms_phyla_df[is.na(extensive_multiple_glms_phyla_df)]<-1
extensive_multiple_glms_phyla_df$padj=p.adjust(extensive_multiple_glms_phyla_df$p, method="BH")
extensive_multiple_glms_phyla_df=filter(extensive_multiple_glms_phyla_df)%>%mutate(Treatment="Extensive")

means=cbind(Extensive_RS_metadata%>%select(Treatment), Phylum_Table[Extensive_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Extensive_RS_metadata%>%select(Treatment), Phylum_Table[Extensive_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
phyla_ext=as.data.frame(cbind(extensive_multiple_glms_phyla_df,
t(means[,rownames(extensive_multiple_glms_phyla_df)]),t(sdt[,rownames(extensive_multiple_glms_phyla_df)])))%>%
  rownames_to_column(var="Taxa")%>%mutate(Taxonomic_Level="Phylum", `N-Fertilization-Intensity`="Extensive")
  
  


extensive_multiple_glms_class= apply(Class_Table[Extensive_RS_metadata$Sample,],2, function(x) multiple_glm_ext(x))
extensive_multiple_glms_class_df <- as.data.frame(do.call(rbind, extensive_multiple_glms_class))
extensive_multiple_glms_class_df[is.na(extensive_multiple_glms_class_df)]<-1
extensive_multiple_glms_class_df$padj=p.adjust(extensive_multiple_glms_class_df$p, method="BH")
extensive_multiple_glms_class_df=filter(extensive_multiple_glms_class_df)%>%mutate(Treatment="Extensive")

means=cbind(Extensive_RS_metadata%>%select(Treatment), Class_Table[Extensive_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Extensive_RS_metadata%>%select(Treatment), Class_Table[Extensive_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

Class_Ext=as.data.frame(cbind(extensive_multiple_glms_class_df,
                              t(means[,rownames(extensive_multiple_glms_class_df),drop = FALSE]),t(sdt[,rownames(extensive_multiple_glms_class_df),drop = FALSE])))%>%
  mutate(Taxonomic_Level="Class", `N-Fertilization-Extensity`="Extensive")
Class_Ext=rownames_to_column(Class_Ext,var="Taxa")


extensive_multiple_glms_order= apply(Order_Table[Extensive_RS_metadata$Sample,],2, function(x) multiple_glm_ext(x))
extensive_multiple_glms_order_df <- as.data.frame(do.call(rbind, extensive_multiple_glms_order))
extensive_multiple_glms_order_df[is.na(extensive_multiple_glms_order_df)]<-1
extensive_multiple_glms_order_df$padj=p.adjust(extensive_multiple_glms_order_df$p, method="BH")
extensive_multiple_glms_order_df=filter(extensive_multiple_glms_order_df)%>%mutate(Treatment="Extensive")

means=cbind(Extensive_RS_metadata%>%select(Treatment), Order_Table[Extensive_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Extensive_RS_metadata%>%select(Treatment), Order_Table[Extensive_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
order_ext=as.data.frame(cbind(extensive_multiple_glms_order_df,
                              t(means[,rownames(extensive_multiple_glms_order_df)]),t(sdt[,rownames(extensive_multiple_glms_order_df)])))%>%
  mutate(Taxonomic_Level="Order", `N-Fertilization-Intensity`="Extensive")
order_ext=rownames_to_column(order_ext,var="Taxa")



extensive_multiple_glms_family= apply(Families_Table[Extensive_RS_metadata$Sample,],2, function(x) multiple_glm_ext(x))
extensive_multiple_glms_family_df <- as.data.frame(do.call(rbind, extensive_multiple_glms_family))
extensive_multiple_glms_family_df[is.na(extensive_multiple_glms_family_df)]<-1
extensive_multiple_glms_family_df$padj=p.adjust(extensive_multiple_glms_family_df$p, method="BH")
extensive_multiple_glms_family_df=filter(extensive_multiple_glms_family_df)%>%mutate(Treatment="Extensive")


means=cbind(Extensive_RS_metadata%>%select(Treatment), Families_Table[Extensive_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Extensive_RS_metadata%>%select(Treatment), Families_Table[Extensive_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
Family_ext=as.data.frame(cbind(extensive_multiple_glms_family_df,
                              t(means[,rownames(extensive_multiple_glms_family_df)]),t(sdt[,rownames(extensive_multiple_glms_family_df)])))%>%
  mutate(Taxonomic_Level="Family", `N-Fertilization-Intensity`="Extensive")
Family_ext=rownames_to_column(Family_ext,var="Taxa")


extensive_multiple_glms_genera= apply(Genera_Table[Extensive_RS_metadata$Sample,],2, function(x) multiple_glm_ext(x))
extensive_multiple_glms_genera_df <- as.data.frame(do.call(rbind, extensive_multiple_glms_genera))
extensive_multiple_glms_genera_df[is.na(extensive_multiple_glms_genera_df)]<-1
extensive_multiple_glms_genera_df$padj=p.adjust(extensive_multiple_glms_genera_df$p, method="BH")
extensive_multiple_glms_genera_df=filter(extensive_multiple_glms_genera_df)%>%mutate(Treatment="Extensive")

means=cbind(Extensive_RS_metadata%>%select(Treatment), Genera_Table[Extensive_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Extensive_RS_metadata%>%select(Treatment), Genera_Table[Extensive_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

Genera_ext=as.data.frame(cbind(extensive_multiple_glms_genera_df,
                               t(means[,rownames(extensive_multiple_glms_genera_df)]),t(sdt[,rownames(extensive_multiple_glms_genera_df)])))%>%
  mutate(Taxonomic_Level="Genera", `N-Fertilization-Intensity`="Extensive")
Genera_ext=rownames_to_column(Genera_ext,var="Taxa")

counts_ext=ASV_Table_count[,Extensive_RS_metadata$Sample]%>%t()%>% as.data.frame()%>% rownames_to_column(var="Sample")%>%
  gather(-Sample, key="ASV", value = "RA")%>% filter(RA>1)%>% group_by(ASV)%>% dplyr::summarise(count=n())%>%
  filter(count>=3)


extensive_multiple_glms_ASV= apply(ASV_Table[Extensive_RS_metadata$Sample,counts_ext$ASV],2, function(x) multiple_glm_ext(x))
extensive_multiple_glms_ASV_df <- as.data.frame(do.call(rbind, extensive_multiple_glms_ASV))
extensive_multiple_glms_ASV_df[is.na(extensive_multiple_glms_ASV_df)]<-1
extensive_multiple_glms_ASV_df$padj=p.adjust(extensive_multiple_glms_ASV_df$p, method="BH")
extensive_multiple_glms_ASV_df=filter(extensive_multiple_glms_ASV_df)%>%mutate(Treatment="Extensive")

means=cbind(Extensive_RS_metadata%>%select(Treatment), ASV_Table[Extensive_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Extensive_RS_metadata%>%select(Treatment), ASV_Table[Extensive_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

ASV_ext=as.data.frame(cbind(extensive_multiple_glms_ASV_df,
                               t(means[,rownames(extensive_multiple_glms_ASV_df)]),t(sdt[,rownames(extensive_multiple_glms_ASV_df)])))%>%
  mutate(Taxonomic_Level="ASV", `N-Fertilization-Intensity`="Extensive")
ASV_ext=rownames_to_column(ASV_ext,var="Taxa")



intensive_multiple_glms_phyla= apply(Phylum_Table[Intensive_RS_metadata$Sample,],2, function(x) multiple_glm_int(x))
intensive_multiple_glms_phyla_df <- as.data.frame(do.call(rbind, intensive_multiple_glms_phyla))
intensive_multiple_glms_phyla_df[is.na(intensive_multiple_glms_phyla_df)]<-1
intensive_multiple_glms_phyla_df$padj=p.adjust(intensive_multiple_glms_phyla_df$p, method="BH")
intensive_multiple_glms_phyla_df=filter(intensive_multiple_glms_phyla_df)%>%mutate(Treatment="Intensive")

means=cbind(Intensive_RS_metadata%>%select(Treatment), Phylum_Table[Intensive_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Intensive_RS_metadata%>%select(Treatment), Phylum_Table[Intensive_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

Phylum_Int=as.data.frame(cbind(intensive_multiple_glms_phyla_df,
                               t(means[,rownames(intensive_multiple_glms_phyla_df)]),t(sdt[,rownames(intensive_multiple_glms_phyla_df)])))%>%
  mutate(Taxonomic_Level="Phylum", `N-Fertilization-Intensity`="Intensive")
Phylum_Int=rownames_to_column(Phylum_Int,var="Taxa")



intensive_multiple_glms_class= apply(Class_Table[Intensive_RS_metadata$Sample,],2, function(x) multiple_glm_int(x))
intensive_multiple_glms_class_df <- as.data.frame(do.call(rbind, intensive_multiple_glms_class))
intensive_multiple_glms_class_df[is.na(intensive_multiple_glms_class_df)]<-1
intensive_multiple_glms_class_df$padj=p.adjust(intensive_multiple_glms_class_df$p, method="BH")
intensive_multiple_glms_class_df=filter(intensive_multiple_glms_class_df)%>%mutate(Treatment="Intensive")

means=cbind(Intensive_RS_metadata%>%select(Treatment), Class_Table[Intensive_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Intensive_RS_metadata%>%select(Treatment), Class_Table[Intensive_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

Class_Int=as.data.frame(cbind(intensive_multiple_glms_class_df,
                               t(means[,rownames(intensive_multiple_glms_class_df)]),t(sdt[,rownames(intensive_multiple_glms_class_df)])))%>%
  mutate(Taxonomic_Level="Class", `N-Fertilization-Intensity`="Intensive")
Class_Int=rownames_to_column(Class_Int,var="Taxa")



intensive_multiple_glms_order= apply(Order_Table[Intensive_RS_metadata$Sample,],2, function(x) multiple_glm_int(x))
intensive_multiple_glms_order_df <- as.data.frame(do.call(rbind, intensive_multiple_glms_order))
intensive_multiple_glms_order_df[is.na(intensive_multiple_glms_order_df)]<-1
intensive_multiple_glms_order_df$padj=p.adjust(intensive_multiple_glms_order_df$p, method="BH")
intensive_multiple_glms_order_df=filter(intensive_multiple_glms_order_df)%>%mutate(Treatment="Intensive")

means=cbind(Intensive_RS_metadata%>%select(Treatment), Order_Table[Intensive_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Intensive_RS_metadata%>%select(Treatment), Order_Table[Intensive_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

Order_Int=as.data.frame(cbind(intensive_multiple_glms_order_df,
                               t(means[,rownames(intensive_multiple_glms_order_df)]),t(sdt[,rownames(intensive_multiple_glms_order_df)])))%>%
  mutate(Taxonomic_Level="Order", `N-Fertilization-Intensity`="Intensive")
Order_Int=rownames_to_column(Order_Int,var="Taxa")


intensive_multiple_glms_family= apply(Families_Table[Intensive_RS_metadata$Sample,],2, function(x) multiple_glm_int(x))
intensive_multiple_glms_family_df <- as.data.frame(do.call(rbind, intensive_multiple_glms_family))
intensive_multiple_glms_family_df[is.na(intensive_multiple_glms_family_df)]<-1
intensive_multiple_glms_family_df$padj=p.adjust(intensive_multiple_glms_family_df$p, method="BH")
intensive_multiple_glms_family_df=filter(intensive_multiple_glms_family_df)%>%mutate(Treatment="Intensive")

means=cbind(Intensive_RS_metadata%>%select(Treatment), Families_Table[Intensive_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Intensive_RS_metadata%>%select(Treatment), Families_Table[Intensive_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

Family_Int=as.data.frame(cbind(intensive_multiple_glms_family_df,
                               t(means[,rownames(intensive_multiple_glms_family_df)]),t(sdt[,rownames(intensive_multiple_glms_family_df)])))%>%
  mutate(Taxonomic_Level="Family", `N-Fertilization-Intensity`="Intensive")
Family_Int=rownames_to_column(Family_Int,var="Taxa")


intensive_multiple_glms_genera= apply(Genera_Table[Intensive_RS_metadata$Sample,],2, function(x) multiple_glm_int(x))
intensive_multiple_glms_genera_df <- as.data.frame(do.call(rbind, intensive_multiple_glms_genera))
intensive_multiple_glms_genera_df[is.na(intensive_multiple_glms_genera_df)]<-1
intensive_multiple_glms_genera_df$padj=p.adjust(intensive_multiple_glms_genera_df$p, method="BH")
intensive_multiple_glms_genera_df=filter(intensive_multiple_glms_genera_df)%>%mutate(Treatment="Intensive")

means=cbind(Intensive_RS_metadata%>%select(Treatment), Genera_Table[Intensive_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Intensive_RS_metadata%>%select(Treatment), Genera_Table[Intensive_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

Genera_Int=as.data.frame(cbind(intensive_multiple_glms_genera_df,
                               t(means[,rownames(intensive_multiple_glms_genera_df)]),t(sdt[,rownames(intensive_multiple_glms_genera_df)])))%>%
  mutate(Taxonomic_Level="Genera", `N-Fertilization-Intensity`="Intensive")
Genera_Int=rownames_to_column(Genera_Int,var="Taxa")

counts_int=ASV_Table_count[Intensive_RS_metadata$Sample]%>%t()%>% as.data.frame()%>% rownames_to_column(var="Sample")%>%
gather(-Sample, key="ASV", value = "RA")%>% filter(RA>1)%>% group_by(ASV)%>% dplyr::summarise(count=n())%>%
  filter(count>=3)

intensive_multiple_glms_ASV= apply(ASV_Table[Intensive_RS_metadata$Sample,(counts_int$ASV)],2, function(x) multiple_glm_int(x))
intensive_multiple_glms_ASV_df <- as.data.frame(do.call(rbind, intensive_multiple_glms_ASV))
intensive_multiple_glms_ASV_df[is.na(intensive_multiple_glms_ASV_df)]<-1
intensive_multiple_glms_ASV_df$padj=p.adjust(intensive_multiple_glms_ASV_df$p, method="BH")
intensive_multiple_glms_ASV_df=filter(intensive_multiple_glms_ASV_df)%>%mutate(Treatment="Intensive")

means=cbind(Intensive_RS_metadata%>%select(Treatment), ASV_Table[Intensive_RS_metadata$Sample,])
means=aggregate(means, .~Treatment, mean)%>%mutate(Treatment=paste0("Mean_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")
sdt=cbind(Intensive_RS_metadata%>%select(Treatment), ASV_Table[Intensive_RS_metadata$Sample,])
sdt=aggregate(sdt, .~Treatment, sd)%>%mutate(Treatment=paste0("Sd_",Treatment))%>%as.data.frame()%>% column_to_rownames(var="Treatment")

ASV_Int=as.data.frame(cbind(intensive_multiple_glms_ASV_df,
                               t(means[,rownames(intensive_multiple_glms_ASV_df)]),t(sdt[,rownames(intensive_multiple_glms_ASV_df)])))%>%
  mutate(Taxonomic_Level="ASV", `N-Fertilization-Intensity`="Intensive")
ASV_Int=rownames_to_column(ASV_Int,var="Taxa")


colnames(Class_Ext)<-colnames(Class_Int)
complete_info_tax_levels= rbind(phyla_ext,Phylum_Int,Class_Ext, Class_Int,order_ext,Order_Int,Family_ext,Family_Int, Genera_Int,Genera_ext)

write_csv(complete_info_tax_levels%>%filter(.,padj<0.05), file = "MS3_results/Resdonders_over_taxonomic_levels_MS3.csv")
write_csv(rbind(ASV_Int, ASV_ext)%>%filter(.,padj<0.05), file = "MS3_results/Resdonders_ASVs_log_regression_MS3.csv")

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

top15Genera= ASV_100[, RS_metadata$Sample]%>%
  as.data.frame()%>% rownames_to_column(var="ASV")%>%
  gather(-ASV, key="Sample", value = "RA")%>%
  full_join(ASVtax_table, by="ASV") %>%
  mutate(ASV=paste0(Phylum,";",Class,";",Order,";",Family, ";",Genus)) %>% 
  full_join(RS_metadata, by="Sample")%>%
  dcast(ASV+Sample~., value.var = "RA", sum)%>%dplyr::select(., RA=".",everything())%>%
  dcast(ASV~., value.var = "RA", median)
top15Genera1=top15Genera[order(top15Genera$., decreasing= TRUE), ]
top15Genera2=top15Genera1[1:15,]

top15Genera_RA= ASV_100[, RS_metadata$Sample]%>%
  as.data.frame()%>% rownames_to_column(var="ASV")%>%
  gather(-ASV, key="Sample", value = "RA")%>%
  full_join(ASVtax_table, by="ASV") %>%
  mutate(ASV=paste0(Phylum,";",Class,";",Order,";",Family, ";",Genus)) %>% filter(ASV%in%c((top15Genera2$ASV)))%>% 
  full_join(RS_metadata, by="Sample")%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%dplyr::select(RA=".", everything())%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  mean)%>%dplyr::select(mean=".", everything())
top15Genera_RA= ASV_100[, RS_metadata$Sample]%>%
  as.data.frame()%>% rownames_to_column(var="ASV")%>%
  gather(-ASV, key="Sample", value = "RA")%>%
  full_join(ASVtax_table, by="ASV") %>%
  mutate(ASV=paste0(Phylum,";",Class,";",Order,";",Family, ";",Genus)) %>% filter(ASV%in%c((top15Genera2$ASV)))%>% 
  full_join(RS_metadata, by="Sample")%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%dplyr::select(RA=".", everything())%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  sd)%>%dplyr::select(sdt=".", everything())%>%full_join(top15Genera_RA)



top15Genera_RA_sign_stats=filter( complete_info_tax_levels, Taxa%in%c( top15Genera2$ASV))%>%select(ASV=Taxa, padj, coeff,  `N-Fertilization_Intensity`=`N-Fertilization-Intensity`)%>%
  full_join(top15Genera_RA)%>%mutate( `p-value_corrected_(BH)`=padj)
top15Genera_RA_sign_stats$padj[top15Genera_RA_sign_stats$padj<0.05]<-"p<0.05"
top15Genera_RA_sign_stats$padj[!top15Genera_RA_sign_stats$padj=="p<0.05"]<-"p>0.05"


library(ggplot2)

top15Genera_RA_sign_stats$mean=round(top15Genera_RA_sign_stats$mean,3)
top15Genera_RA_sign_stats_plot=ggplot(top15Genera_RA_sign_stats%>%mutate(`N-Fertilization_Intensity`=str_replace(`N-Fertilization_Intensity`, "ensive",".")), aes(y=fct_reorder(ASV,mean), x=paste0( Treatment), fill=mean, shape=padj)) +facet_wrap(~paste(`N-Fertilization_Intensity`)) + geom_point(size=25)+
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

png("MS3_results/top15_Genera_RA_plot_with_log_regression.png", width = 2100, height=1500)
print(top15Genera_RA_sign_stats_plot) 
dev.off()

write_excel_csv(top15Genera_RA_sign_stats, file = "MS3_results/top15_genera_RA_plot_with_log_regression.csv")


top15Classes= ASV_100[, RS_metadata$Sample]%>%
  as.data.frame()%>% rownames_to_column(var="ASV")%>%
  gather(-ASV, key="Sample", value = "RA")%>%
  full_join(ASVtax_table, by="ASV") %>%
  mutate(ASV=paste0(Phylum,";",Class)) %>% full_join(RS_metadata, by="Sample")%>%
  dcast(ASV+Sample~., value.var = "RA", sum)%>%dplyr::select(., RA=".",everything())%>% 
  dcast(ASV~., value.var = "RA", median)
top15Classes1=top15Classes[order(top15Classes$., decreasing= TRUE), ]
top15Classes2=top15Classes1[1:15,]

top15Classes_RA= ASV_100[, RS_metadata$Sample]%>%
  as.data.frame()%>% rownames_to_column(var="ASV")%>%
  gather(-ASV, key="Sample", value = "RA")%>%
  full_join(ASVtax_table, by="ASV") %>%
  mutate(ASV=paste0(Phylum,";",Class)) %>% filter(ASV%in%c((top15Classes2$ASV)))%>% 
  full_join(RS_metadata, by="Sample")%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%dplyr::select(RA=".", everything())%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  mean)%>%dplyr::select(mean=".", everything())



top15Classes_RA= ASV_100[, RS_metadata$Sample]%>%
  as.data.frame()%>% rownames_to_column(var="ASV")%>%
  gather(-ASV, key="Sample", value = "RA")%>%
  full_join(ASVtax_table, by="ASV") %>%
  mutate(ASV=paste0(Phylum,";",Class)) %>% filter(ASV%in%c((top15Classes2$ASV)))%>% 
  full_join(RS_metadata, by="Sample")%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%dplyr::select(RA=".", everything())%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  sd)%>%dplyr::select(sdt=".", everything())%>% full_join(top15Classes_RA)


top15Classes_RA_sign_stats=filter( complete_info_tax_levels, Taxa%in%c( top15Classes2$ASV))%>%select(ASV=Taxa, padj, coeff,  `N-Fertilization_Intensity`=`N-Fertilization-Intensity`)%>%
  full_join(top15Classes_RA)%>%mutate( `p-value_corrected_(BH)`=padj)

top15Classes_RA_sign_stats$padj[top15Classes_RA_sign_stats$padj<0.05]<-"p<0.05"
top15Classes_RA_sign_stats$padj[!top15Classes_RA_sign_stats$padj=="p<0.05"]<-"p>0.05"


top15Classes_RA_sign_stats$mean=round(top15Classes_RA_sign_stats$mean,3)
top15Classes_RA_sign_stats_plot=ggplot(top15Classes_RA_sign_stats%>%mutate(`N-Fertilization_Intensity`=str_replace(`N-Fertilization_Intensity`, "ensive",".")), aes(y=fct_reorder(ASV,mean), x=paste0( Treatment), fill=mean, shape=padj)) +facet_wrap(~paste(`N-Fertilization_Intensity`)) + geom_point(size=25)+
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

png("MS3_results/top15_Classes_RA_plot_with_log_regression.png", width = 1500, height=1500)
print(top15Classes_RA_sign_stats_plot) 
dev.off()
write_excel_csv(top15Classes_RA_sign_stats, file = "MS3_results/top15_classes_RA_plot_with_log_regression.csv")


top15Phyla= ASV_100[, RS_metadata$Sample]%>%
  as.data.frame()%>% rownames_to_column(var="ASV")%>%
  gather(-ASV, key="Sample", value = "RA")%>%
  full_join(ASVtax_table, by="ASV") %>%
  mutate(ASV=paste0(Phylum)) %>% full_join(RS_metadata, by="Sample")%>%
  dcast(ASV+Sample~., value.var = "RA", sum)%>%dplyr::select(., RA=".",everything())%>% 
  dcast(ASV~., value.var = "RA", median)
top15Phyla1=top15Phyla[order(top15Phyla$., decreasing= TRUE), ]
top15Phyla2=top15Phyla1[1:15,]

top15Phyla_RA= ASV_100[, RS_metadata$Sample]%>%
  as.data.frame()%>% rownames_to_column(var="ASV")%>%
  gather(-ASV, key="Sample", value = "RA")%>%
  full_join(ASVtax_table, by="ASV") %>%
  mutate(ASV=paste0(Phylum)) %>% filter(ASV%in%c((top15Phyla2$ASV)))%>% 
  full_join(RS_metadata, by="Sample")%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%dplyr::select(RA=".", everything())%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  mean)%>%dplyr::select(mean=".", everything())

top15Phyla_RA= ASV_100[, RS_metadata$Sample]%>%
  as.data.frame()%>% rownames_to_column(var="ASV")%>%
  gather(-ASV, key="Sample", value = "RA")%>%
  full_join(ASVtax_table, by="ASV") %>%
  mutate(ASV=paste0(Phylum)) %>% filter(ASV%in%c((top15Phyla2$ASV)))%>% 
  full_join(RS_metadata, by="Sample")%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%dplyr::select(RA=".", everything())%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  sd)%>%dplyr::select(sdt=".", everything())%>%full_join(top15Phyla_RA)


top15Phyla_RA_sign_stats=filter( complete_info_tax_levels, Taxa%in%c( top15Phyla2$ASV))%>%select(ASV=Taxa, padj, coeff,  `N-Fertilization_Intensity`=`N-Fertilization-Intensity`)%>%
  full_join(top15Phyla_RA)%>%mutate( `p-value_corrected_(BH)`=padj)

top15Phyla_RA_sign_stats$padj[top15Phyla_RA_sign_stats$padj<0.05]<-"p<0.05"
top15Phyla_RA_sign_stats$padj[!top15Phyla_RA_sign_stats$padj=="p<0.05"]<-"p>0.05"


library(ggplot2)
top15Phyla_RA_sign_stats$mean=round(top15Phyla_RA_sign_stats$mean,3)
top15Phyla_RA_sign_stats_plot=ggplot(top15Phyla_RA_sign_stats%>%mutate(`N-Fertilization_Intensity`=str_replace(`N-Fertilization_Intensity`, "ensive",".")), 
                                     aes(y=fct_reorder(ASV,mean), x=paste0( Treatment), fill=mean, shape=padj)) +facet_wrap(~paste(`N-Fertilization_Intensity`)) + geom_point(size=25)+
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

png("MS3_results/top15_Phyla_RA_plot_with_log_regression.png", width = 1200, height=1500)
print(top15Phyla_RA_sign_stats_plot) 
dev.off()
write_excel_csv(top15Phyla_RA_sign_stats, file = "MS3_results/top15_phyla_RA_plot_with_log_regression.csv")

combined_ASVs=rbind(ASV_ext, ASV_Int)%>%filter(padj<0.05)%>%mutate(logFC=log((10^3*Mean_BMc+1)/(10^3*Mean_Ctrl+1)))%>%
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



png("MS3_results/combined_edgeR_plot_with_log.png", width = 4500, height=3000)
print(two_managements) 
dev.off()



top15ASV= ASV_100[, RS_metadata$Sample]%>%
  as.data.frame()%>% rownames_to_column(var="ASV")%>%
  gather(-ASV, key="Sample", value = "RA")%>%
  full_join(ASVtax_table, by="ASV") %>%
  mutate(ASV=paste0(Phylum,";",Class,";",Order,";", Family,";", Genus,";", ASV)) %>% full_join(RS_metadata, by="Sample")%>%
  dcast(ASV+Sample~., value.var = "RA", sum)%>%dplyr::select(., RA=".",everything())%>% 
  dcast(ASV~., value.var = "RA", median)
top15ASV1=top15ASV[order(top15ASV$., decreasing= TRUE), ]
top15ASV2=top15ASV1[1:15,]

top15ASV_RA= ASV_100[, RS_metadata$Sample]%>%
  as.data.frame()%>% rownames_to_column(var="ASV")%>%
  gather(-ASV, key="Sample", value = "RA")%>%
  full_join(ASVtax_table, by="ASV") %>%
  mutate(ASV=paste0(Phylum,";",Class,";",Order,";", Family,";", Genus,";", ASV)) %>%
  filter(ASV%in%c((top15ASV2$ASV)))%>% 
  full_join(RS_metadata, by="Sample")%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%dplyr::select(RA=".", everything())%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  mean)%>%dplyr::select(mean=".", everything())

top15ASV_RA= ASV_100[, RS_metadata$Sample]%>%
  as.data.frame()%>% rownames_to_column(var="ASV")%>%
  gather(-ASV, key="Sample", value = "RA")%>%
  full_join(ASVtax_table, by="ASV") %>%
  mutate(ASV=paste0(Phylum,";",Class,";",Order,";", Family,";", Genus,";", ASV)) %>%
  filter(ASV%in%c((top15ASV2$ASV)))%>% 
  full_join(RS_metadata, by="Sample")%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment+Sample~., value.var = "RA",  sum)%>%dplyr::select(RA=".", everything())%>%
  dcast(ASV+`N-Fertilization_Intensity`+Treatment~., value.var = "RA",  sd)%>%dplyr::select(sdt=".", everything())%>%full_join(top15ASV_RA)


top15ASV_RA_sign_stats=filter( rbind(ASV_ext,ASV_Int), 
Taxa%in%c( top15ASV2$ASV))%>%select(ASV=Taxa, padj, coeff,  `N-Fertilization_Intensity`=`N-Fertilization-Intensity`)%>%
  full_join(top15ASV_RA)%>%mutate( `p-value_corrected_(BH)`=padj)
top15ASV_RA_sign_stats$padj[top15ASV_RA_sign_stats$padj<0.05]<-"p<0.05"
top15ASV_RA_sign_stats$padj[!top15ASV_RA_sign_stats$padj=="p<0.05"]<-"p>0.05"
top15ASV_RA_sign_stats$Taxonomy=top15ASV_RA_sign_stats$ASV
top15ASV_RA_sign_stats=  separate(top15ASV_RA_sign_stats,
           Taxonomy,into=c(
             "Phylum",
             "Class",
             "Order",
             "Family",
             "Genus",
             "ASV"),
           sep=";") %>%
  mutate(., Acronym=paste0( Genus,"_",ASV ))
for(i in 1:nrow(top15ASV_RA_sign_stats)){
  gn =paste0(top15ASV_RA_sign_stats[i,"Acronym"]) 
#  ASVtax_table[i,7]=paste0(gn)
  print(paste0("Round ",i))
  if (top15ASV_RA_sign_stats[i,"Genus"]=="Unclassified"){
    gn=paste0(top15ASV_RA_sign_stats[i,"Family"],"_",top15ASV_RA_sign_stats[i,"ASV"])
  }
  print(gn)
  top15ASV_RA_sign_stats[i,"Acronym"]=paste0(gn)
}


library(ggplot2)
top15ASV_RA_sign_stats=mutate(top15ASV_RA_sign_stats, Acronym=str_replace(Acronym, ".*.*Rhizobium", "Rhizobium"))
top15ASV_RA_sign_stats$mean=round(top15ASV_RA_sign_stats$mean,3)
top15ASV_RA_sign_stats_plot=ggplot(top15ASV_RA_sign_stats%>%mutate(`N-Fertilization_Intensity`=str_replace(`N-Fertilization_Intensity`, "ensive",".")), 
                                   aes(y=fct_reorder(Acronym,mean), x=paste0( Treatment), fill=mean, shape=padj)) +facet_wrap(~paste(`N-Fertilization_Intensity`)) +
  geom_point(size=25)+
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
  theme(legend.position="bottom", legend.direction = "horizontal",
        legend.box = "vertical",legend.justification = "center") +guides(shape=guide_legend(nrow=2,byrow=TRUE)) +
  theme( axis.line = element_line(colour = "black", 
size = 2, linetype = "solid"))+geom_text(aes(label=round(mean,2)),size=8)

png("MS3_results/top15_ASVs_RA_plot_with_log_regression.png", width = 1250, height=1500)
print(top15ASV_RA_sign_stats_plot) 
dev.off()
write_excel_csv(top15ASV_RA_sign_stats, file = "MS3_results/top15_ASVs_RA_plot_with_log_regression.csv")

