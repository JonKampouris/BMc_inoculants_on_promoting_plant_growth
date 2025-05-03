
library(ggplot2) #ggplot2: Fancy graphs
library(ggpubr) #ggplot expansion for more fancy graphs
library(readxl) # upload files  
library(readr) # read the files
library(tidyr) # Data handling
library(RColorBrewer) # Colours for fancy graphs
library(tibble)# Data handling: rownames/columnanmes
library(stringr) # Manipulate strings
library(igraph)
library(dplyr)
library(reshape2)
library(rstatix)

# Load the deltacts
deltacts= read_excel(paste0("~/LTE2020_geneExpression_deltaCT.xlsx"))
metadata_genes_expr <- read_excel("~/metadata_genes_expr.xlsx")%>%
  mutate(Gene=str_replace(Gene, "Zm", ""))


# Here I created fake groups, it can be removed if you have the Gene infomartion metadata
deltacts2=deltacts%>%gather(-TREATMENT, key="Gene", value = "CT")%>%
  mutate(Gene=str_replace(Gene, "Zm", ""))



deltacts2=separate(deltacts2, "TREATMENT", sep = "-", into = c("Tillage","Intensity", "BMc", "Block"))
deltacts2=as.data.frame(deltacts2)
library(rstatix)
stat.test <- deltacts2 %>%
  group_by(Intensity, Gene) %>%
  t_test(CT ~ BMc) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()%>%select(p.adj.signif, Gene, Intensity)
stat.test

#Since the data is already log2 transformed the FC can be compute with subtraction 
deltacts3=dcast(Gene+Intensity~BMc, data = deltacts2, mean, value.var = "CT")%>%mutate(FC=(I-C))
deltacts3=dcast(Gene+Intensity~BMc, data = deltacts2, mean, value.var = "CT")%>%mutate(FC=(I-C))

deltacts4=full_join(deltacts3,stat.test)
deltacts4$p.adj.signif[deltacts4$p.adj.signif=="ns"]<-""




#Enter the correct metadata here with one column the gene and the other collumn the metabolic traits
deltacts4=full_join(deltacts4, metadata_genes_expr)

heatmp_plot=
ggplot( deltacts4, aes(x=paste0(Intensity, "."), y=Gene, fill= (FC)))+
  scale_fill_distiller(name= expression(bold(Log2FC)),palette = "RdBu")+
  geom_tile()+geom_text(aes(label=p.adj.signif), size=10)+
  facet_grid(Annotation~.,scales = "free",space="free",switch = "x") + 
  theme(strip.placement = "outside",
        strip.text.y.right =  element_text(angle = 0,vjust=1,colour="black", size=30,face = "bold"),
        strip.background = element_rect(fill="white", colour="white", linewidth = 2),
        panel.spacing = unit(0.2,"cm"))+xlab("")+ylab("")+scale_y_discrete( position = "right")+
  
  theme( axis.text.y = element_text( size=25, face = "italic", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
  theme(panel.margin=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1) )+
  theme(legend.key.size = unit(3,"cm"))+geom_vline(xintercept = 1.5, size=1)
heatmp_plot
 
png(file="test_heatmap_genes.png", width = 1200, height = 1000)
print(heatmp_plot)
dev.off()
