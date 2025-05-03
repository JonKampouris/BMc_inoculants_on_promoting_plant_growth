library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)
library(ggpubr)
library(readr)



top15ASV_RA_sign_stats=
read.csv( file = "~/top15_ASVs_RA_plot_with_log_regression.csv")%>%mutate(Amplicon="16S rRNA gene")%>%mutate(Level="ASV")%>%
  mutate(ASV=paste0(Genus,";", ASV))

top15Phyla_RA_sign_stats=
read.csv( file = "~/top15_phyla_RA_plot_with_log_regression.csv")%>%mutate(Amplicon="16S rRNA gene")%>%mutate(Level="Phylum")

top15Fungal_ASVsRA_sign_stats=
  read.csv( file = "MS3_results/top15_fungal_OTUs_RA_plot_with_log_regression.csv")%>%mutate(Amplicon="ITS")%>%
  mutate(Level="OTU")%>%select(ASV=OTU2, everything())%>%mutate(ASV=str_replace(ASV, "OTU", "ASV"))
top15_fungal_Phyla_RA_sign_stats= read.csv( file = "MS3_results/top15_fungal_phyla_RA_plot_with_log_regression.csv")%>%mutate(Amplicon="ITS")%>%
  mutate(Level="Phylum")%>%select(ASV=OTU,  everything())

complete_data_ploting= rbind( top15ASV_RA_sign_stats%>%select(colnames(top15Phyla_RA_sign_stats)), 
       top15Fungal_ASVsRA_sign_stats%>%select(colnames(top15Phyla_RA_sign_stats)),
       top15Phyla_RA_sign_stats,top15_fungal_Phyla_RA_sign_stats%>%select(colnames(top15Phyla_RA_sign_stats))
)

complete_data_ploting$mean=round(complete_data_ploting$mean,2)
###############################################################
###################################################################
top15ASV_RA_sign_stats$coeff2=top15ASV_RA_sign_stats$coeff
top15ASV_RA_sign_stats$coeff2[top15ASV_RA_sign_stats$coeff2>1]="Increased in BMc (p<0.05)"
top15ASV_RA_sign_stats$coeff2[
top15ASV_RA_sign_stats$coeff2!="Increased in BMc (p<0.05)"]="Increased in Ctrl (p<0.05)"
for(i in 1:nrow(top15ASV_RA_sign_stats)){
 if( top15ASV_RA_sign_stats[i,"padj"]=="p>0.05"){
   top15ASV_RA_sign_stats[i,"coeff2"]="ns"
 }
}

top15ASV_RA_sign_stats$coeff2=factor(
  top15ASV_RA_sign_stats$coeff2,
  levels = c("Increased in BMc (p<0.05)","Increased in Ctrl (p<0.05)","ns")
)
  
ASV_plot=ggplot(top15ASV_RA_sign_stats%>%
                                       mutate(`N-Fertilization_Intensity`=str_replace(`N.Fertilization_Intensity`, "ensive",""))%>%
                        mutate(`N-Fertilization_Intensity`=str_replace(`N-Fertilization_Intensity`, "t","t."))%>%
                        mutate(ASV=str_replace(ASV, ".*.*.Rhi","Rhi"))%>%select(padj, everything()), aes(y=fct_reorder(ASV,mean),
                                                                                                         x=paste0( Treatment), fill=coeff2)) +
  facet_grid(~`N-Fertilization_Intensity`) + geom_tile(size=70, colour="black", linewidth=1, alpha=0.6)+
  scale_fill_manual(values = c("#ff9900",
                               "#146eb4",
                               
                               "#f7f7f7"),
                    name="")  + scale_alpha_manual(values = c(1, 0.4), name="") +#scale_shape_manual(values = c(21,24), "") +
  theme_pubr(legend = "right") + 
  theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
  theme( strip.text = element_text( size=65, face = "bold", colour="white")) +
  theme(strip.background =  element_rect( fill="darkblue")) +
  
  theme(axis.text.x = element_text(size = 55, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=75, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=75, face="bold", colour = "black"))+
  theme(legend.key.size = unit(2,"cm")) +xlab("")+ylab("")+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+
  theme(legend.position="bottom", legend.direction = "horizontal",
        legend.box = "vertical",legend.justification = "center", legend.key.width =  unit(2.3,"cm")) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE)) + guides(alpha="none")+
  theme( axis.line = element_line(colour = "black", 
                                  size = 2, linetype = "solid"))+
  geom_text(aes(label=paste0(round(top15ASV_RA_sign_stats$mean,2), "±", round(top15ASV_RA_sign_stats$sdt,2))),size=13)+
  ggtitle("B) Bacterial ASVs")+theme(title = element_text( size=75, face="bold", colour = "black"))

######################################################
######################################################
#paste0(round(top15Phyla_RA_sign_stats$mean,2), "±", round(top15Phyla_RA_sign_stats$sdt,2))

top15Phyla_RA_sign_stats$coeff2=top15Phyla_RA_sign_stats$coeff
top15Phyla_RA_sign_stats$coeff2[top15Phyla_RA_sign_stats$coeff2>1]="Increased in BMc (p<0.05)"
top15Phyla_RA_sign_stats$coeff2[
  top15Phyla_RA_sign_stats$coeff2!="Increased in BMc (p<0.05)"]="Increased in Ctrl (p<0.05)"
for(i in 1:nrow(top15Phyla_RA_sign_stats)){
  if( top15Phyla_RA_sign_stats[i,"padj"]=="p>0.05"){
    top15Phyla_RA_sign_stats[i,"coeff2"]="ns"
  }
}

top15Phyla_RA_sign_stats$coeff2=factor(
  top15Phyla_RA_sign_stats$coeff2,
  levels = c("Increased in BMc (p<0.05)","Increased in Ctrl (p<0.05)","ns")
)



Bac_phyla_plot=ggplot(top15Phyla_RA_sign_stats%>%
                  mutate(`N-Fertilization_Intensity`=str_replace(`N.Fertilization_Intensity`, "ensive",""))%>%
                  mutate(`N-Fertilization_Intensity`=str_replace(`N-Fertilization_Intensity`, "t","t."))%>%
                  mutate(ASV=str_replace(ASV, ".*.*.Rhi","Rhi"))%>%select(padj, everything()),
                  aes(y=fct_reorder(ASV,mean),
                  x=paste0( Treatment), fill=coeff2)) +
  
  # scale_fill_manual(name="",values=c( "#ff9900","#146eb4"))+
  # theme(strip.background = element_rect(colour="white", fill="white"))+
  
  
  facet_grid(~`N-Fertilization_Intensity`) + geom_tile(size=70, colour="black", linewidth=1, alpha=0.6)+
  scale_fill_manual(values = c("#ff9900",
                               "#146eb4",
                               
                               "#f7f7f7"),
                    name="")  + scale_alpha_manual(values = c(1, 0.4), name="") +#scale_shape_manual(values = c(21,24), "") +
  theme_pubr(legend = "right") + 
  theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
  theme( strip.text = element_text( size=65, face = "bold", colour="white")) +
  theme(strip.background =  element_rect( fill="darkblue")) +
  
  theme(axis.text.x = element_text(size = 55, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=75, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=75, face="bold", colour = "black"))+
  theme(legend.key.size = unit(2,"cm")) +xlab("")+ylab("")+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+
  theme(legend.position="bottom", legend.direction = "horizontal",
        legend.box = "vertical",legend.justification = "center", legend.key.width =  unit(2.3,"cm")) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE)) + guides(alpha="none")+
  theme( axis.line = element_line(colour = "black", 
                                  size = 2, linetype = "solid"))+
  geom_text(aes(label=paste0(round(top15Phyla_RA_sign_stats$mean,2), "±", round(top15Phyla_RA_sign_stats$sdt,2))),size=13)+
  ggtitle("A) Bacterial Phyla")+theme(title = element_text( size=75, face="bold", colour = "black"))

Bac_phyla_plot

#paste0(round(top15Phyla_RA_sign_stats$mean,2), "±", round(top15Phyla_RA_sign_stats$sdt,2))
#####################################################
#######################################################

top15Fungal_ASVsRA_sign_stats$coeff2=top15Fungal_ASVsRA_sign_stats$coeff
top15Fungal_ASVsRA_sign_stats$coeff2[top15Fungal_ASVsRA_sign_stats$coeff2>1]="Increased in BMc (p<0.05)"
top15Fungal_ASVsRA_sign_stats$coeff2[
  top15Fungal_ASVsRA_sign_stats$coeff2!="Increased in BMc (p<0.05)"]="Increased in Ctlr (p<0.05)"
for(i in 1:nrow(top15Fungal_ASVsRA_sign_stats)){
  if( top15Fungal_ASVsRA_sign_stats[i,"padj"]=="p>0.05"){
    top15Fungal_ASVsRA_sign_stats[i,"coeff2"]="ns"
  }
}

top15Fungal_ASVsRA_sign_stats$coeff2=factor(
  top15Fungal_ASVsRA_sign_stats$coeff2,
  levels = c("Increased in BMc (p<0.05)","Increased in Ctlr (p<0.05)","ns")
)



Fungal_ASVsplot=ggplot(top15Fungal_ASVsRA_sign_stats%>%
                  mutate(`N-Fertilization_Intensity`=str_replace(`N.Fertilization_Intensity`, "ensive",""))%>%
                  mutate(`N-Fertilization_Intensity`=str_replace(`N-Fertilization_Intensity`, "t","t."))%>%
                  mutate(ASV=str_replace(ASV, ".*.*.Rhi","Rhi"))%>%select(padj, everything()), aes(y=fct_reorder(ASV,mean),
                                                                                                   x=paste0( Treatment), fill=coeff2)) +
  facet_grid(~`N-Fertilization_Intensity`)+ geom_tile(size=70, colour="black", linewidth=1, alpha=0.6)+
  scale_fill_manual(values = c("#ff9900",
                               "#146eb4",
                               
                               "#f7f7f7"),
                    name="")  + scale_alpha_manual(values = c(1, 0.4), name="") +#scale_shape_manual(values = c(21,24), "") +
  theme_pubr(legend = "right") + 
  theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
  theme( strip.text = element_text( size=65, face = "bold", colour="white")) +
  theme(strip.background =  element_rect( fill="darkblue")) +
  
  theme(axis.text.x = element_text(size = 55, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=75, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=75, face="bold", colour = "black"))+
  theme(legend.key.size = unit(2,"cm")) +xlab("")+ylab("")+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+
  theme(legend.position="bottom", legend.direction = "horizontal",
        legend.box = "vertical",legend.justification = "center", legend.key.width =  unit(2.3,"cm")) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE)) + guides(alpha="none")+
  theme( axis.line = element_line(colour = "black", 
                                  size = 2, linetype = "solid"))+
  geom_text(aes(label=paste0(round(top15Fungal_ASVsRA_sign_stats$mean,2), "±", round(top15Fungal_ASVsRA_sign_stats$sdt,2))),size=13)+
  ggtitle("D) Fungal ASVs")+theme(title = element_text( size=75, face="bold", colour = "black"))


Fungal_ASVsplot

#####################################################
#######################################################

top15_fungal_Phyla_RA_sign_stats$coeff2=top15_fungal_Phyla_RA_sign_stats$coeff
top15_fungal_Phyla_RA_sign_stats$coeff2[top15_fungal_Phyla_RA_sign_stats$coeff2>1]="Increased in BMc (p<0.05)"
top15_fungal_Phyla_RA_sign_stats$coeff2[
  top15_fungal_Phyla_RA_sign_stats$coeff2!="Increased in BMc (p<0.05)"]="Increased in Ctlr (p<0.05)"
for(i in 1:nrow(top15_fungal_Phyla_RA_sign_stats)){
  if( top15_fungal_Phyla_RA_sign_stats[i,"padj"]=="p>0.05"){
    top15_fungal_Phyla_RA_sign_stats[i,"coeff2"]="ns"
  }
}

top15_fungal_Phyla_RA_sign_stats$coeff2=factor(
  top15_fungal_Phyla_RA_sign_stats$coeff2,
  levels = c("Increased in BMc (p<0.05)","Increased in Ctlr (p<0.05)","ns")
)



fungal_phyla_plot=ggplot(top15_fungal_Phyla_RA_sign_stats%>%
                  mutate(`N-Fertilization_Intensity`=str_replace(`N.Fertilization_Intensity`, "ensive",""))%>%
                  mutate(`N-Fertilization_Intensity`=str_replace(`N-Fertilization_Intensity`, "t","t."))%>%
                  mutate(ASV=str_replace(ASV, ".*.*.Rhi","Rhi"))%>%select(padj, everything()), aes(y=fct_reorder(ASV,mean),
                                                                                                   x=paste0( Treatment), fill=coeff2)) +
  facet_grid(~`N-Fertilization_Intensity`)+ geom_tile(size=70, colour="black", linewidth=1, alpha=0.6)+
  scale_fill_manual(values = c("#ff9900",
                               "#146eb4",
                               
                               "#f7f7f7"),
                    name="")  + scale_alpha_manual(values = c(1, 0.4), name="") +#scale_shape_manual(values = c(21,24), "") +
  theme_pubr(legend = "right") + 
  theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
  theme( strip.text = element_text( size=65, face = "bold", colour="white")) +
  theme(strip.background =  element_rect( fill="darkblue")) +
  
  theme(axis.text.x = element_text(size = 55, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=75, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=25, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=75, face="bold", colour = "black"))+
  theme(legend.key.size = unit(2,"cm")) +xlab("")+ylab("")+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+
  theme(legend.position="bottom", legend.direction = "horizontal",
        legend.box = "vertical",legend.justification = "center", legend.key.width =  unit(2.3,"cm")) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE)) + guides(alpha="none")+
  theme( axis.line = element_line(colour = "black", 
                                  size = 2, linetype = "solid"))+
  geom_text(aes(label=paste0(round(top15_fungal_Phyla_RA_sign_stats$mean,2), "±", round(top15_fungal_Phyla_RA_sign_stats$sdt,2))),size=13)+ 
  ggtitle("C) Fungal Phyla")+theme(title = element_text( size=75, face="bold", colour = "black"))
  
fungal_phyla_plot

library(ggpubr)
plots=ggarrange(
                Bac_phyla_plot ,ASV_plot,   
                fungal_phyla_plot,Fungal_ASVsplot,  common.legend = TRUE, legend = "bottom",align = "hv")


plots

png("~/top15_both_amplicons_RA_plot_with_log_regression.png",
    res=180, width = 10000, height=9000)
print(plots) 
dev.off()

