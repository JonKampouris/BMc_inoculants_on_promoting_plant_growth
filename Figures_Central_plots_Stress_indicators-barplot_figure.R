
# LINEAR MIXED MODEL ###
library(dplyr)
library(agricolae)
library(vegan)
library(car)

library(multcomp)
library(ggplot2)
library(ggpmisc)
library(broom)
library(lme4)
library(phyloseq)
library(MuMIn)
library(ggpubr)
# 
# install.packages(c("nlme",lib=.libPaths()[2]))
# install.packages("MuMIn")
 
#### START ####
df.StrI<-read.table("Stress_Ind_Central_and_RW_legend.txt",h=T)
df.StrI

str(df.StrI)
head(df.StrI)
tail(df.StrI)

# subsetting central dataset
df.StrI.central<-subset(df.StrI, PLOT=="central")
str(df.StrI.central)


df.StrI.central$TREATMENT = factor(df.StrI.central$TREATMENT, 
                                 levels=c("Ext-Ctrl","Ext-BMc", "Int-Ctrl","Int-BMc"))



df.StrI.central$Bms = factor(df.StrI.central$Bms, 
                           levels=c("Ctrl","BMc"))


# subsetting central shoot dataset
df.StrI.central.shoot<-subset(df.StrI.central, LOCATION =="shoot")

# # subsetting central root dataset
# df.StrI.central.root<-subset(df.StrI.central, LOCATION =="root")


##### ascorbate_peroxidase_activity Central plots - Shoots #####
df.APX<-df.StrI.central.shoot[,c(3:9,10)]
df.APX

APX.summary <- summarySE(df.APX, measurevar="APX", groupvars=c("FERTILIZATION","Bms"))
APX.summary

APX.summary2 <- APX.summary
APX.summary2$FERTILIZATION <- factor(APX.summary2$FERTILIZATION)

val.APX<-c("b","a","b", "a")

bp.APX.C.shoot<-ggplot(data = APX.summary2) +
  geom_bar(
    aes(
      y = APX,
      x = FERTILIZATION,
      fill = Bms,
    ),
    colour = "black",
    size = 0.2,
    position = "dodge",
    stat = "identity",
    width=.7
  ) +
  geom_errorbar(
    aes(
      ymin = APX - sd,
      ymax = APX + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  geom_text(
    aes(
      y = APX + sd,
      x = FERTILIZATION,
      label = val.APX,
      group = Bms
    ),
    size = rel(7),
    vjust = - 1.5,
    position = position_dodge(0.7)
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    breaks = pretty_breaks(),
    expand = expansion(mult = c(0,0.2))
  ) +
  labs(y = "APX [ng g-1 FW]")+
  scale_x_discrete(
    name = ""
  ) +
  scale_fill_manual(
    values=c("#146eb4","#ff9900")) +
  ggtitle("") +
  # general layout
  theme(
    axis.text.x = element_text(size=rel(2),color="#000000"),
    axis.text.y = element_text(size=rel(2),color="#000000"),
    axis.title.y.left = element_text(colour = 'black', size = rel(2)), # text y axis
    axis.line = element_line(colour = 'black', size = rel(0.4)),
    axis.ticks = element_line(colour = "black", size = rel(0.5)),
    panel.border = element_rect(colour = "black", fill = NA, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    legend.key=element_blank(),
    legend.position = "none", # the position has to be adapted in dependence of plot proportions. Check position only at exported plot
    legend.margin = margin(0, 0, 0, 0),
    #       legend.position = "right",
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.6, 'cm')) + ylab(bquote(APX~"["*unit~g^-1*" FW]"))

bp.APX.C.shoot



##### SOD Central plots - Shoots #####
df.SOD<-df.StrI.central.shoot[,c(3:9,15)]
df.SOD

SOD.summary <- summarySE(df.SOD, measurevar="SOD", groupvars=c("FERTILIZATION","Bms"))
SOD.summary

SOD.summary2 <- SOD.summary
SOD.summary2$FERTILIZATION <- factor(SOD.summary2$FERTILIZATION)

val.SOD<-c("b","a","b", "a")

bp.SOD.C.shoot<-ggplot(data = SOD.summary2) +
  geom_bar(
    aes(
      y = SOD,
      x = FERTILIZATION,
      fill = Bms,
    ),
    colour = "black",
    size = 0.2,
    position = "dodge",
    stat = "identity",
    width=.7
  ) +
  geom_errorbar(
    aes(
      ymin = SOD - sd,
      ymax = SOD + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  geom_text(
    aes(
      y = SOD + sd,
      x = FERTILIZATION,
      label = val.SOD,
      group = Bms
    ),
    size = rel(7),
    vjust = - 1.5,
    position = position_dodge(0.7)
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    breaks = pretty_breaks(),
    expand = expansion(mult = c(0,0.2))
  ) +
  labs(y = "SOD [ng g-1 FW]")+
  scale_x_discrete(
    name = ""
  ) +
  scale_fill_manual(
    values=c("#146eb4","#ff9900")) +
  ggtitle("") +
  # general layout
  theme(
    axis.text.x = element_text(size=rel(2),color="#000000"),
    axis.text.y = element_text(size=rel(2),color="#000000"),
    axis.title.y.left = element_text(colour = 'black', size = rel(2)), # text y axis
    axis.line = element_line(colour = 'black', size = rel(0.4)),
    axis.ticks = element_line(colour = "black", size = rel(0.5)),
    panel.border = element_rect(colour = "black", fill = NA, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    legend.key=element_blank(),
    legend.position = "none", # the position has to be adapted in dependence of plot proportions. Check position only at exported plot
    legend.margin = margin(0, 0, 0, 0),
    #       legend.position = "right",
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.6, 'cm')) + ylab(bquote(SOD~"["*unit~g^-1*" FW]"))

bp.SOD.C.shoot




##### H2O2 Central plots - Shoots #####
df.H2O2<-df.StrI.central.shoot[,c(3:9,12)]
df.H2O2

H2O2.summary <- summarySE(df.H2O2, measurevar="H2O2", groupvars=c("FERTILIZATION","Bms"))
H2O2.summary

H2O2.summary2 <- H2O2.summary
H2O2.summary2$FERTILIZATION <- factor(H2O2.summary2$FERTILIZATION)

val.H2O2<-c("a","b","a", "b")

bp.H2O2.C.shoot<-ggplot(data = H2O2.summary2) +
  geom_bar(
    aes(
      y = H2O2,
      x = FERTILIZATION,
      fill = Bms,
    ),
    colour = "black",
    size = 0.2,
    position = "dodge",
    stat = "identity",
    width=.7
  ) +
  geom_errorbar(
    aes(
      ymin = H2O2 - sd,
      ymax = H2O2 + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  geom_text(
    aes(
      y = H2O2 + sd,
      x = FERTILIZATION,
      label = val.H2O2,
      group = Bms
    ),
    size = rel(7),
    vjust = - 1.5,
    position = position_dodge(0.7)
  ) +
scale_y_continuous(
  limits = c(0, NA),
  breaks = pretty_breaks(),
  expand = expansion(mult = c(0,0.2))
) +
  labs(y = "H2O2 [ng g-1 FW]")+
  scale_x_discrete(
    name = ""
  ) +
  scale_fill_manual(
    values=c("#146eb4","#ff9900")) +
  ggtitle("") +
  # general layout
  theme(
    axis.text.x = element_text(size=rel(2),color="#000000"),
    axis.text.y = element_text(size=rel(2),color="#000000"),
    axis.title.y.left = element_text(colour = 'black', size = rel(2)), # text y axis
    axis.line = element_line(colour = 'black', size = rel(0.4)),
    axis.ticks = element_line(colour = "black", size = rel(0.5)),
    panel.border = element_rect(colour = "black", fill = NA, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    legend.key=element_blank(),
    legend.position = "none", # the position has to be adapted in dependence of plot proportions. Check position only at exported plot
    legend.margin = margin(0, 0, 0, 0),
    #       legend.position = "right",
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.6, 'cm')) + ylab(bquote(H2O2~"["*µmol~g^-1*" FW]"))

bp.H2O2.C.shoot




##### proline Central plots - Shoots #####
df.proline<-df.StrI.central.shoot[,c(3:9,14)]
df.proline

proline.summary <- summarySE(df.proline, measurevar="proline", groupvars=c("FERTILIZATION","Bms"))
proline.summary

proline.summary2 <- proline.summary
proline.summary2$FERTILIZATION <- factor(proline.summary2$FERTILIZATION)

# val.proline<-c("b","a","b", "a")

bp.proline.C.shoot<-ggplot(data = proline.summary2) +
  geom_bar(
    aes(
      y = proline,
      x = FERTILIZATION,
      fill = Bms,
    ),
    colour = "black",
    size = 0.2,
    position = "dodge",
    stat = "identity",
    width=.7
  ) +
  geom_errorbar(
    aes(
      ymin = proline - sd,
      ymax = proline + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  # geom_text(
  #   aes(
  #     y = proline + sd,
  #     x = FERTILIZATION,
  #     label = val.proline,
  #     group = Bms
  #   ),
  #   size = rel(7),
  #   vjust = - 1.5,
  #   position = position_dodge(0.7)
  # ) +
  scale_y_continuous(
    limits = c(0, NA),
    breaks = pretty_breaks(),
    expand = expansion(mult = c(0,0.2))
  ) +
  labs(y = "proline [ng g-1 FW]")+
  scale_x_discrete(
    name = ""
  ) +
  scale_fill_manual(
    values=c("#146eb4","#ff9900")) +
  ggtitle("") +
  # general layout
  theme(
    axis.text.x = element_text(size=rel(2),color="#000000"),
    axis.text.y = element_text(size=rel(2),color="#000000"),
    axis.title.y.left = element_text(colour = 'black', size = rel(2)), # text y axis
    axis.line = element_line(colour = 'black', size = rel(0.4)),
    axis.ticks = element_line(colour = "black", size = rel(0.5)),
    panel.border = element_rect(colour = "black", fill = NA, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    legend.key=element_blank(),
    legend.position = "none", # the position has to be adapted in dependence of plot proportions. Check position only at exported plot
    legend.margin = margin(0, 0, 0, 0),
    #       legend.position = "right",
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.6, 'cm')) + ylab(bquote(Proline~"["*mg^-1*" FW]"))

bp.proline.C.shoot





##### GB Central plots - Shoots #####
df.GB<-df.StrI.central.shoot[,c(3:9,11)]
df.GB

GB.summary <- summarySE(df.GB, measurevar="GB", groupvars=c("FERTILIZATION","Bms"))
GB.summary

GB.summary2 <- GB.summary
GB.summary2$FERTILIZATION <- factor(GB.summary2$FERTILIZATION)

 val.GB<-c("b","a","b", "a")

bp.GB.C.shoot<-ggplot(data = GB.summary2) +
  geom_bar(
    aes(
      y = GB,
      x = FERTILIZATION,
      fill = Bms,
    ),
    colour = "black",
    size = 0.2,
    position = "dodge",
    stat = "identity",
    width=.7
  ) +
  geom_errorbar(
    aes(
      ymin = GB - sd,
      ymax = GB + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  geom_text(
    aes(
      y = GB + sd,
      x = FERTILIZATION,
      label = val.GB,
      group = Bms
    ),
    size = rel(7),
    vjust = - 1.5,
    position = position_dodge(0.7)
  ) +
scale_y_continuous(
  limits = c(0, NA),
  breaks = pretty_breaks(),
  expand = expansion(mult = c(0,0.2))
) +
  labs(y = "GB [ng g-1 FW]")+
  scale_x_discrete(
    name = ""
  ) +
  scale_fill_manual(
    values=c("#146eb4","#ff9900")) +
  ggtitle("") +
  # general layout
  theme(
    axis.text.x = element_text(size=rel(2),color="#000000"),
    axis.text.y = element_text(size=rel(2),color="#000000"),
    axis.title.y.left = element_text(colour = 'black', size = rel(2)), # text y axis
    axis.line = element_line(colour = 'black', size = rel(0.4)),
    axis.ticks = element_line(colour = "black", size = rel(0.5)),
    panel.border = element_rect(colour = "black", fill = NA, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    legend.key=element_blank(),
    legend.position = "none", # the position has to be adapted in dependence of plot proportions. Check position only at exported plot
    legend.margin = margin(0, 0, 0, 0),
    #       legend.position = "right",
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.6, 'cm')) + ylab(bquote(GB~"["*mg^-1*" FW]"))

bp.GB.C.shoot






##### phenolics Central plots - Shoots #####
df.phenolics<-df.StrI.central.shoot[,c(3:9,13)]
df.phenolics

phenolics.summary <- summarySE(df.phenolics, measurevar="phenolics", groupvars=c("FERTILIZATION","Bms"))
phenolics.summary

phenolics.summary2 <- phenolics.summary
phenolics.summary2$FERTILIZATION <- factor(phenolics.summary2$FERTILIZATION)

# val.phenolics<-c("b","a","b", "a")

bp.phenolics.C.shoot<-ggplot(data = phenolics.summary2) +
  geom_bar(
    aes(
      y = phenolics,
      x = FERTILIZATION,
      fill = Bms,
    ),
    colour = "black",
    size = 0.2,
    position = "dodge",
    stat = "identity",
    width=.7
  ) +
  geom_errorbar(
    aes(
      ymin = phenolics - sd,
      ymax = phenolics + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  # geom_text(
  #   aes(
  #     y = phenolics + sd,
  #     x = FERTILIZATION,
  #     label = val.phenolics,
  #     group = Bms
  #   ),
  #   size = rel(7),
  #   vjust = - 1.5,
  #   position = position_dodge(0.7)
  # ) +
  scale_y_continuous(
    limits = c(0, NA),
    breaks = pretty_breaks(),
    expand = expansion(mult = c(0,0.2))
  ) +
  labs(y = "phenolics [ng g-1 FW]")+
  scale_x_discrete(
    name = ""
  ) +
  scale_fill_manual(
    values=c("#146eb4","#ff9900")) +
  ggtitle("") +
  # general layout
  theme(
    axis.text.x = element_text(size=rel(2),color="#000000"),
    axis.text.y = element_text(size=rel(2),color="#000000"),
    axis.title.y.left = element_text(colour = 'black', size = rel(2)), # text y axis
    axis.line = element_line(colour = 'black', size = rel(0.4)),
    axis.ticks = element_line(colour = "black", size = rel(0.5)),
    panel.border = element_rect(colour = "black", fill = NA, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    legend.key=element_blank(),
    legend.position = "none", # the position has to be adapted in dependence of plot proportions. Check position only at exported plot
    legend.margin = margin(0, 0, 0, 0),
    #       legend.position = "right",
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.6, 'cm')) + ylab(bquote(phenolics~"["*mg~gallic~acid~equivalent~g^-1*"]"))

bp.phenolics.C.shoot






##### total_antioxidants Central plots - Shoots #####
df.total_antioxidants<-df.StrI.central.shoot[,c(3:9,16)]
df.total_antioxidants

total_antioxidants.summary <- summarySE(df.total_antioxidants, measurevar="total_antioxidants", groupvars=c("FERTILIZATION","Bms"))
total_antioxidants.summary

total_antioxidants.summary2 <- total_antioxidants.summary
total_antioxidants.summary2$FERTILIZATION <- factor(total_antioxidants.summary2$FERTILIZATION)

 val.total_antioxidants<-c("b","a","b", "a")

bp.total_antioxidants.C.shoot<-ggplot(data = total_antioxidants.summary2) +
  geom_bar(
    aes(
      y = total_antioxidants,
      x = FERTILIZATION,
      fill = Bms,
    ),
    colour = "black",
    size = 0.2,
    position = "dodge",
    stat = "identity",
    width=.7
  ) +
  geom_errorbar(
    aes(
      ymin = total_antioxidants - sd,
      ymax = total_antioxidants + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  geom_text(
    aes(
      y = total_antioxidants + sd,
      x = FERTILIZATION,
      label = val.total_antioxidants,
      group = Bms
    ),
    size = rel(7),
    vjust = - 1.5,
    position = position_dodge(0.7)
  ) +
scale_y_continuous(
  limits = c(0, NA),
  breaks = pretty_breaks(),
  expand = expansion(mult = c(0,0.2))
) +
  labs(y = "Total antioxidants [%]")+
  scale_x_discrete(
    name = ""
  ) +
  scale_fill_manual(
    values=c("#146eb4","#ff9900")) +
  ggtitle("") +
  # general layout
  theme(
    axis.text.x = element_text(size=rel(2),color="#000000"),
    axis.text.y = element_text(size=rel(2),color="#000000"),
    axis.title.y.left = element_text(colour = 'black', size = rel(2)), # text y axis
    axis.line = element_line(colour = 'black', size = rel(0.4)),
    axis.ticks = element_line(colour = "black", size = rel(0.5)),
    panel.border = element_rect(colour = "black", fill = NA, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    legend.key=element_blank(),
    legend.position = "none", # the position has to be adapted in dependence of plot proportions. Check position only at exported plot
    legend.margin = margin(0, 0, 0, 0),
    #       legend.position = "right",
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.6, 'cm')) 

bp.total_antioxidants.C.shoot




###### ggarrange ######
ggarrange( bp.APX.C.shoot, bp.SOD.C.shoot, bp.H2O2.C.shoot, bp.H2O2.C.shoot, bp.proline.C.shoot, bp.GB.C.shoot, bp.total_antioxidants.C.shoot,bp.phenolics.C.shoot,
           labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
           ncol = 4, nrow = 2,   font.label = list(size = 18, color = "black", face = "bold"))



ggsave("barchart_stress_indicators_nolegend.pdf", plot = last_plot(),
       device = "pdf", units = c("cm"),width = 60, height = 32, dpi = 300) 





