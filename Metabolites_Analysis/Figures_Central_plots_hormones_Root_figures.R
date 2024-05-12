
# LINEAR MIXED MODEL ###
library(dplyr)
library(agricolae)
library(vegan)
library(car)
library(Rmisc)
library(ggthemes)
library(multcomp)
library(ggplot2)
library(ggpmisc)
library(broom)
library(lme4)
library(phyloseq)
library(MuMIn)
library(ggpubr)
library(scales)

# install.packages(c("nlme",lib=.libPaths()[2]))
# install.packages("MuMIn")

###### START #####
df.dc<-read.table("Hormones_Central_and_RW_legend.txt",h=T)
df.dc

str(df.dc)
head(df.dc)
tail(df.dc)

# subsetting central dataset
df.dc.central<-subset(df.dc, PLOT=="central")
str(df.dc.central)

df.dc.central$TREATMENT = factor(df.dc.central$TREATMENT, 
                                 levels=c("Ext-Ctrl","Ext-BMc", "Int-Ctrl","Int-BMc"))

df.dc.central$Bms = factor(df.dc.central$Bms, 
                                 levels=c("Ctrl","BMc"))
# subsetting central shoot dataset

# subsetting central shoot dataset
df.dc.central.shoot<-subset(df.dc.central, LOCATION =="shoot")

# subsetting central root dataset
df.dc.central.root<-subset(df.dc.central, LOCATION =="root")

str(df.dc.central.root)


##### IAA Central plots - roots #####
df.IAA<-df.dc.central.root[,c(3:9,13)]
df.IAA

IAA.summary <- summarySE(df.IAA, measurevar="IAA", groupvars=c("FERTILIZATION","Bms"))
IAA.summary

IAA.summary2 <- IAA.summary
IAA.summary2$FERTILIZATION <- factor(IAA.summary2$FERTILIZATION)

val.IAA<-c("b","a","b", "a")

bp.IAA.C.root<-ggplot(data = IAA.summary2) +
  geom_bar(
    aes(
      y = IAA,
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
      ymin = IAA - sd,
      ymax = IAA + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  geom_text(
    aes(
      y = IAA + sd,
      x = FERTILIZATION,
      label = val.IAA,
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
  labs(y = "IAA [ng g-1 FW]")+
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
    legend.key.size = unit(0.6, 'cm')) + ylab(bquote(IAA~"["*ng~g^-1*" FW]"))

bp.IAA.C.root
#   



##### CK Central plots - roots #####
df.ck<-df.dc.central.root[,c(3:9,11)]
df.ck

ck.summary <- summarySE(df.ck, measurevar="CK", groupvars=c("FERTILIZATION","Bms"))
ck.summary

ck.summary2 <- ck.summary
ck.summary2$FERTILIZATION <- factor(ck.summary2$FERTILIZATION)

val.CK<-c("b","a","b", "a")

bp.CK.C.root<-ggplot(data = ck.summary2) +
  geom_bar(
    aes(
      y = CK,
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
      ymin = CK - sd,
      ymax = CK + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  geom_text(
    aes(
      y = CK + sd,
      x = FERTILIZATION,
      label = val.CK,
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
  labs(y = "CK [ng g-1 FW]")+
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
    axis.title.y.left = element_text(colour = 'black', size = rel(2)),
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
        legend.key.size = unit(0.6, 'cm')) + ylab(bquote(CK~"["*ng~g^-1*" FW]"))

bp.CK.C.root



##### GA Central plots - roots #####
df.GA<-df.dc.central.root[,c(3:9,12)]
df.GA

GA.summary <- summarySE(df.GA, measurevar="GA", groupvars=c("FERTILIZATION","Bms"))
GA.summary

GA.summary2 <- GA.summary
GA.summary2$FERTILIZATION <- factor(GA.summary2$FERTILIZATION)

# val.GA<-c("ab","a","b", "ab")

bp.GA.C.root<-ggplot(data = GA.summary2) +
  geom_bar(
    aes(
      y = GA,
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
      ymin = GA - sd,
      ymax = GA + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  # geom_text(
  #   aes(
  #     y = GA + sd,
  #     x = FERTILIZATION,
  #     label = val.GA,
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
  labs(y = "GA [ng g-1 FW]")+
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
    axis.title.y.left = element_text(colour = 'black', size = rel(2)),
    axis.line = element_line(colour = 'black', size = rel(0.4)),
    axis.ticks = element_line(colour = "black", size = rel(0.5)),
    panel.border = element_rect(colour = "black", fill = NA, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    legend.key=element_blank(),
legend.position = "none", # the position has to be adapted in dependence of plot proportions. CheGA position only at exported plot
    legend.margin = margin(0, 0, 0, 0),
    #       legend.position = "right",
    legend.text = element_text(size = rel(0.9)),
        legend.key.size = unit(0.6, 'cm')) +  ylab(bquote(GA~"["*ng~g^-1*" FW]"))

bp.GA.C.root

##### ABA Central plots - roots ####
df.aba<-df.dc.central.root[,c(3:9,10)]
df.aba

### barlplot with standard deviation!!!!!
# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
aba.summary <- summarySE(df.aba, measurevar="ABA", groupvars=c("FERTILIZATION","Bms"))
aba.summary

aba.summary2 <- aba.summary
aba.summary2$FERTILIZATION <- factor(aba.summary2$FERTILIZATION)

val.ABA<-c("a","b","a", "b")


bp.ABA.C.root<-ggplot(data = aba.summary2) +
  geom_bar(
    aes(
      y = ABA,
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
      ymin = ABA - sd,
      ymax = ABA + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  # geom_text(
  #   aes(
  #     y = ABA + sd,
  #     x = FERTILIZATION,
  #     label = val.ABA,
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
  labs(y = "ABA [ng g-1 FW]")+
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
    axis.title.y.left = element_text(colour = 'black', size = rel(2)),
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
    legend.key.size = unit(0.6, 'cm'))+ ylab(bquote(ABA~"["*ng~g^-1*" FW]"))

bp.ABA.C.root
#   




##### JA Central plots - roots #####
df.JA<-df.dc.central.root[,c(3:9,14)]
df.JA

JA.summary <- summarySE(df.JA, measurevar="JA", groupvars=c("FERTILIZATION","Bms"))
JA.summary

JA.summary2 <- JA.summary
JA.summary2$FERTILIZATION <- factor(JA.summary2$FERTILIZATION)

val.JA<-c("b","a","b", "a")

bp.JA.C.root<-ggplot(data = JA.summary2) +
  geom_bar(
    aes(
      y = JA,
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
      ymin = JA - sd,
      ymax = JA + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  geom_text(
    aes(
      y = JA + sd,
      x = FERTILIZATION,
      label = val.JA,
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
  labs(y = "JA [ng g-1 FW]")+
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
    axis.title.y.left = element_text(colour = 'black', size = rel(2)), 
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
        legend.key.size = unit(0.6, 'cm')) + ylab(bquote(JA~"["*ng~g^-1*" FW]"))
bp.JA.C.root




##### SA Central plots - roots #####
df.SA<-df.dc.central.root[,c(3:9,15)]
df.SA

SA.summary <- summarySE(df.SA, measurevar="SA", groupvars=c("FERTILIZATION","Bms"))
SA.summary

SA.summary2 <- SA.summary
SA.summary2$FERTILIZATION <- factor(SA.summary2$FERTILIZATION)

val<-c("c","a","d", "b")

bp.SA.C.root<-ggplot(data = SA.summary2) +
  geom_bar(
    aes(
      y = SA,
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
      ymin = SA - sd,
      ymax = SA + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
#   geom_text(
#     aes(
#       y = SA + sd,
#       x = FERTILIZATION,
#       label = val,
#       group = Bms
#     ),
#     size = rel(7),
#     vjust = - 1.5,
#     position = position_dodge(0.7))
# +
   scale_y_continuous(
    limits = c(0, NA),
    breaks = pretty_breaks(),
    expand = expansion(mult = c(0,0.2))
  ) +
  labs(y = "SA [ng g-1 FW]")+
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
    axis.title.y.left = element_text(colour = 'black', size = rel(2)),
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
        legend.key.size = unit(0.6, 'cm')) + ylab(bquote(SA~"["*ng~g^-1*" FW]"))



bp.SA.C.root


## ggarrange ####
figure.hormones<-ggarrange( bp.IAA.C.root, bp.CK.C.root, bp.GA.C.root, bp.ABA.C.root, bp.JA.C.root, bp.SA.C.root, 
                           labels = c("Indole-3-acetic acid", "    Cytokinin", "    Gibberellins", "Abscisic Acid", " Jasmonic acid", " Salicylic acid"),
                           ncol = 3, nrow = 2,   font.label = list(size = 18, color = "black", face = "bold"))
figure.hormones


ggarrange( bp.IAA.C.root, bp.CK.C.root, bp.GA.C.root, bp.ABA.C.root, bp.JA.C.root, bp.SA.C.root, 
           labels = c("Indole-3-acetic acid", "    Cytokinin", "    Gibberellins", "Abscisic Acid", " Jasmonic acid", " Salicylic acid"),
           ncol = 3, nrow = 2,   font.label = list(size = 18, color = "black", face = "bold"))

ggarrange( bp.IAA.C.root, bp.CK.C.root, bp.GA.C.root, bp.ABA.C.root, bp.JA.C.root, bp.SA.C.root, 
           labels = c("A", "B", "C", "D", "E", "F"),
           ncol = 3, nrow = 2,   font.label = list(size = 18, color = "black", face = "bold"))



ggsave("barchart_root_hormones_nolegend.pdf", plot = last_plot(),
       device = "pdf", units = c("cm"),width = 40, height = 32, dpi = 300) 


# figure.hormones<-ggarrange(bp.ABA.C.root, bp.CK.C.root, bp.GA.C.root, bp.IAA.C.root, bp.JA.C.root, bp.SA.C.root, 
#                            labels = c("Abscisic Acid", "    Citokinin", "    Gibberellins","Indole-3-acetic acid", " Jasmonic acid", " Salicylic acid"),
#                            ncol = 3, nrow = 2,   font.label = list(size = 18, color = "black", face = "bold"))
# figure.hormones
# 
# 
# figure.hormones<-ggarrange(bp.ABA.C.root, bp.CK.C.root, bp.GA.C.root, bp.IAA.C.root, bp.JA.C.root, bp.SA.C.root, 
#                            labels = c("Abscisic Acid", "    Citokinin", "    Gibberellins","Indole-3-acetic acid", "   Jasmonic acid", " Salicylic acid"),
#                            ncol = 3, nrow = 2,   font.label = list(size = 20, color = "black", face = "bold"))
# figure.hormones
# 
# 
# #####
# figure.ABA.GA<-ggarrange(bp.ABA.C.root, bp.GA.C.root, 
#           labels = c("Abscisic Acid", "Gibberellins"),
#           ncol = 2, nrow = 1,   font.label = list(size = 18, color = "black", face = "bold"))
# 
# figure.ABA.GA
# 
# figure.ABA.GA.JA<-ggarrange(bp.ABA.C.root, bp.GA.C.root, bp.JA.C.root,
#                          labels = c("Abscisic Acid", "Gibberellins", "Jasmonic acid"),
#                          ncol = 3, nrow = 1,   font.label = list(size = 18, color = "black", face = "bold"))
# 
# figure.ABA.GA.JA

############ plot to get the legend #### 

bp.IAA.C.root.legend<-ggplot(data = IAA.summary2) +
  geom_bar(
    aes(
      y = IAA,
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
      ymin = IAA - sd,
      ymax = IAA + sd,
      x = FERTILIZATION,
      group = Bms
    ),
    width = 0.2,
    size = 0.2,
    position = position_dodge(0.7))+
  geom_text(
    aes(
      y = IAA + sd,
      x = FERTILIZATION,
      label = val,
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
  labs(y = "IAA [ng g-1 FW]")+
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
    legend.position = c(0.10, 0.927), # the position has to be adapted in dependence of plot proportions. Check position only at exported plot
    legend.margin = margin(0, 0, 0, 0),
    #       legend.position = "right",
    legend.text = element_text(size = rel(1.4)),
    legend.key.size = unit(1.0, 'cm'))

bp.IAA.C.root.legend
#   
