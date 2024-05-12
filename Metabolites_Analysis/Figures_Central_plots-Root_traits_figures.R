
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
library(Rmisc)
library(scales)

# install.packages(c("nlme",lib=.libPaths()[2]))
# install.packages("MuMIn")

df.dc.RW.root<-read.table("Root_traits_LTE1_2020_RW_sorted_biomass.txt",h=T)
df.dc.RW.root

str(df.dc.RW.root)
head(df.dc.RW.root)
tail(df.dc.RW.root)


df.dc.RW.root$TREATMENT = factor(df.dc.RW.root$TREATMENT, 
                                   levels=c("Ext-Ctrl","Ext-BMc", "Int-Ctrl","Int-BMc"))


df.dc.RW.root$Bms = factor(df.dc.RW.root$Bms, levels=c("Ctrl","BMc"))


##### total_root_length RW plots - root #####
df.TRL<-df.dc.RW.root[,c(3:9,10)]
df.TRL

TRL.summary <- summarySE(df.TRL, measurevar="TRL", groupvars=c("FERTILIZATION","Bms"))
TRL.summary

TRL.summary2 <- TRL.summary
TRL.summary2$FERTILIZATION <- factor(TRL.summary2$FERTILIZATION)

# val.TRL<-c("b","a","b", "a")

bp.TRL.C.root<-ggplot(data = TRL.summary2) +
  geom_bar(
    aes(
      y = TRL,
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
      ymin = TRL - sd,
      ymax = TRL + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  # geom_text(
  #   aes(
  #     y = TRL + sd,
  #     x = FERTILIZATION,
  #     label = val.TRL,
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
  labs(y = "Total root length [cm]")+
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

bp.TRL.C.root
#   



##### avg_diameter RW plots - root #####

df.AVD<-df.dc.RW.root[,c(3:9,11)]
df.AVD

AVD.summary <- summarySE(df.AVD, measurevar="AVD", groupvars=c("FERTILIZATION","Bms"))
AVD.summary

AVD.summary2 <- AVD.summary
AVD.summary2$FERTILIZATION <- factor(AVD.summary2$FERTILIZATION)

val.AVD<-c("a","b","ab", "ab")

bp.AVD.C.root<-ggplot(data = AVD.summary2) +
  geom_bar(
    aes(
      y = AVD,
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
      ymin = AVD - sd,
      ymax = AVD + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
   geom_text(
     aes(
       y = AVD + sd,
       x = FERTILIZATION,
       label = val.AVD,
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
  labs(y = "Root diameter [mm]")+
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

bp.AVD.C.root
#   


##### AMF_colonization RW plots - root #####

df.AMF<-df.dc.RW.root[,c(3:9,12)]
df.AMF

AMF.summary <- summarySE(df.AMF, measurevar="AMF", groupvars=c("FERTILIZATION","Bms"))
AMF.summary

AMF.summary2 <- AMF.summary
AMF.summary2$FERTILIZATION <- factor(AMF.summary2$FERTILIZATION)

val.AMF<-c("a","a","b", "b")

bp.AMF.C.root<-ggplot(data = AMF.summary2) +
  geom_bar(
    aes(
      y = AMF,
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
      ymin = AMF - sd,
      ymax = AMF + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  geom_text(
    aes(
      y = AMF + sd,
      x = FERTILIZATION,
      label = val.AMF,
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
  labs(y = "AMF root colonization [%]")+
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

bp.AMF.C.root
#   




##### root_hair_size RW plots - root #####

df.RHL<-df.dc.RW.root[,c(3:9,14)]
df.RHL

RHL.summary <- summarySE(df.RHL, measurevar="RHL", groupvars=c("FERTILIZATION","Bms"))
RHL.summary

RHL.summary2 <- RHL.summary
RHL.summary2$FERTILIZATION <- factor(RHL.summary2$FERTILIZATION)

# val.RHL<-c("a","c","b", "c")

bp.RHL.C.root<-ggplot(data = RHL.summary2) +
  geom_bar(
    aes(
      y = RHL,
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
      ymin = RHL - sd,
      ymax = RHL + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  # geom_text(
  #   aes(
  #     y = RHL + sd,
  #     x = FERTILIZATION,
  #     label = val.RHL,
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
  labs(y = "Root hair length [mm]")+
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

bp.RHL.C.root



##### Plant biomass - shoots #####

df.SDM<-df.dc.RW.root[,c(3:9,15)]
df.SDM

SDM.summary <- summarySE(df.SDM, measurevar="SDM", groupvars=c("FERTILIZATION","Bms"))
SDM.summary

SDM.summary2 <- SDM.summary
SDM.summary2$FERTILIZATION <- factor(SDM.summary2$FERTILIZATION)

val.SDM<-c("c","b","c", "a")

bp.SDM.C.root<-ggplot(data = SDM.summary2) +
  geom_bar(
    aes(
      y = SDM,
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
      ymin = SDM - sd,
      ymax = SDM + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  geom_text(
    aes(
      y = SDM + sd,
      x = FERTILIZATION,
      label = val.SDM,
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
  labs(y = "Root hair length root colonization [%]")+
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
    legend.key.size = unit(0.6, 'cm')) + ylab(bquote(SDM~"["*g~plant^-1*"]"))

bp.SDM.C.root



##### ggarrange ####

ggarrange(bp.SDM.C.root, bp.TRL.C.root, bp.AVD.C.root,bp.RHL.C.root,bp.AMF.C.root,   
                              labels = c("A", "B", "C", "D", "E"),
                              ncol = 5, nrow = 1,  font.label = list(size = 18, color = "black", face = "bold"))



ggsave("barchart_root_traits.pdf", plot = last_plot(),
       device = "pdf", units = c("cm"),width = 60, height = 13, dpi = 300) 

