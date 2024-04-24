
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

str(df.dc.central.shoot)


##### ABA Central plots - shoots ####
df.aba<-df.dc.central.shoot[,c(3:9,10)]
df.aba

### barlplot with standard deviation!!!!!
# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
aba.summary <- summarySE(df.aba, measurevar="ABA", groupvars=c("FERTILIZATION","Bms"))
aba.summary

aba.summary2 <- aba.summary
aba.summary2$FERTILIZATION <- factor(aba.summary2$FERTILIZATION)

val<-c("a","b","a", "b")


bp.ABA.C.shoot<-ggplot(data = aba.summary2) +
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
  geom_text(
    aes(
      y = ABA + sd,
      x = FERTILIZATION,
      label = val,
      group = Bms
    ),
    size = rel(3),
    vjust = - 1.5,
    position = position_dodge(0.7)
  ) +
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
    axis.text.x = element_text(size=rel(1.5),color="#000000"),
    axis.text.y = element_text(size=rel(1),color="#000000"),
    axis.title.y.left = element_text(colour = 'black', size = rel(1)),
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
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.9, 'cm'))

bp.ABA.C.shoot
#   


hist(df.dc.central.shoot$ABA)
shapiro.test(df.dc.central.shoot$ABA)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.central.shoot$ABA, rnorm(n=length(df.dc.central.shoot$ABA), mean=mean(df.dc.central.shoot$ABA), sd=sd(df.dc.central.shoot$ABA)))

# create ANOVA model
model.bp.ABA.C.shoot<-lm( ABA  ~ TREATMENT, data = df.dc.central.shoot)
anova(model.bp.ABA.C.shoot)

# staristical model
model.bp.ABA.C.shoot2<-aov( ABA  ~ FERTILIZATION * Bms, data = df.dc.central.shoot)
anova(model.bp.ABA.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.central.shoot$ABA, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.central.shoot$ABA)

hist(residuals(model.bp.ABA.C.shoot))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.ABA.C.shoot )
MSerror<-deviance(model.bp.ABA.C.shoot)/df
comparison <- LSD.test(model.bp.ABA.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.ABA.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison

a# log transformation
df.dc.central.shoot$ABA_log10<-log10(df.dc.central.shoot$ABA)
df.dc.central.shoot

bp.ABA_log10.C.shoot<-ggplot(df.dc.central.shoot,aes(TREATMENT, ABA_log10 ,  color=Bms)) + geom_boxplot()
bp.ABA_log10.C.shoot

hist(df.dc.central.shoot$ABA_log10)
shapiro.test(df.dc.central.shoot$ABA_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.central.shoot$ABA_log10, rnorm(n=length(df.dc.central.shoot$ABA_log10), mean=mean(df.dc.central.shoot$ABA_log10), sd=sd(df.dc.central.shoot$ABA_log10)))

# create ANOVA model
model.bp.ABA_log10.C.shoot<-lm( ABA_log10  ~ TREATMENT, data = df.dc.central.shoot)
anova(model.bp.ABA_log10.C.shoot)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.central.shoot$ABA_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.central.shoot$ABA_log10)

hist(residuals(model.bp.ABA_log10.C.shoot))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.ABA_log10.C.shoot )
MSerror<-deviance(model.bp.ABA_log10.C.shoot)/df
comparison <- LSD.test(model.bp.ABA_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.ABA_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison




##### CK Central plots - shoots #####
df.ck<-df.dc.central.shoot[,c(3:9,11)]
df.ck

ck.summary <- summarySE(df.ck, measurevar="CK", groupvars=c("FERTILIZATION","Bms"))
ck.summary

ck.summary2 <- ck.summary
ck.summary2$FERTILIZATION <- factor(ck.summary2$FERTILIZATION)

bp.CK.C.shoot<-ggplot(data = ck.summary2) +
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
      label = val,
      group = Bms
    ),
    size = rel(3),
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
    axis.text.x = element_text(size=rel(1.5),color="#000000"),
    axis.text.y = element_text(size=rel(1),color="#000000"),
    axis.title.y.left = element_text(colour = 'black', size = rel(1)),
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
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.9, 'cm'))

bp.CK.C.shoot




hist(df.dc.central.shoot$CK)
shapiro.test(df.dc.central.shoot$CK)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.central.shoot$CK, rnorm(n=length(df.dc.central.shoot$CK), mean=mean(df.dc.central.shoot$CK), sd=sd(df.dc.central.shoot$CK)))

# create ANOVA model
model.bp.CK.C.shoot<-lm( CK  ~ TREATMENT, data = df.dc.central.shoot)
anova(model.bp.CK.C.shoot)

# statistical model
model.bp.CK.C.shoot2<-lm( CK  ~ FERTILIZATION * Bms, data = df.dc.central.shoot)
anova(model.bp.CK.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.central.shoot$CK, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.central.shoot$CK)

hist(residuals(model.bp.CK.C.shoot))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.CK.C.shoot )
MSerror<-deviance(model.bp.CK.C.shoot)/df
comparison <- LSD.test(model.bp.CK.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.CK.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.dc.central.shoot$CK_log10<-log10(df.dc.central.shoot$CK)
df.dc.central.shoot

bp.CK_log10.C.shoot<-ggplot(df.dc.central.shoot,aes(TREATMENT, CK_log10 ,  color=Bms)) + geom_boxplot()
bp.CK_log10.C.shoot

hist(df.dc.central.shoot$CK_log10)
shapiro.test(df.dc.central.shoot$CK_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.central.shoot$CK_log10, rnorm(n=length(df.dc.central.shoot$CK_log10), mean=mean(df.dc.central.shoot$CK_log10), sd=sd(df.dc.central.shoot$CK_log10)))

# create ANOVA model
model.bp.CK_log10.C.shoot<-lm( CK_log10  ~ TREATMENT, data = df.dc.central.shoot)
anova(model.bp.CK_log10.C.shoot)

# statistical model
model.bp.CK_log10.C.shoot2<-lm( CK_log10  ~ FERTILIZATION * Bms, data = df.dc.central.shoot)
anova(model.bp.CK_log10.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.central.shoot$CK_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.central.shoot$CK_log10)

hist(residuals(model.bp.CK_log10.C.shoot))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.CK_log10.C.shoot )
MSerror<-deviance(model.bp.CK_log10.C.shoot)/df
comparison <- LSD.test(model.bp.CK_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.CK_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison



##### GA Central plots - shoots #####
df.GA<-df.dc.central.shoot[,c(3:9,12)]
df.GA

GA.summary <- summarySE(df.GA, measurevar="GA", groupvars=c("FERTILIZATION","Bms"))
GA.summary

GA.summary2 <- GA.summary
GA.summary2$FERTILIZATION <- factor(GA.summary2$FERTILIZATION)

val<-c("ab","a","b", "ab")

bp.GA.C.shoot<-ggplot(data = GA.summary2) +
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
  geom_text(
    aes(
      y = GA + sd,
      x = FERTILIZATION,
      label = val,
      group = Bms
    ),
    size = rel(3),
    vjust = - 1.5,
    position = position_dodge(0.7)
  ) +
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
    axis.text.x = element_text(size=rel(1.5),color="#000000"),
    axis.text.y = element_text(size=rel(1),color="#000000"),
    axis.title.y.left = element_text(colour = 'black', size = rel(1)),
    axis.line = element_line(colour = 'black', size = rel(0.4)),
    axis.ticks = element_line(colour = "black", size = rel(0.5)),
    panel.border = element_rect(colour = "black", fill = NA, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    legend.key=element_blank(),
    legend.position = c(0.10, 0.927), # the position has to be adapted in dependence of plot proportions. CheGA position only at exported plot
    legend.margin = margin(0, 0, 0, 0),
    #       legend.position = "right",
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.9, 'cm'))

bp.GA.C.shoot


hist(df.dc.central.shoot$GA)
shapiro.test(df.dc.central.shoot$GA)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.central.shoot$GA, rnorm(n=length(df.dc.central.shoot$GA), mean=mean(df.dc.central.shoot$GA), sd=sd(df.dc.central.shoot$GA)))

# create ANOVA model
model.bp.GA.C.shoot<-lm( GA  ~ TREATMENT, data = df.dc.central.shoot)
anova(model.bp.GA.C.shoot)

# staristical model
model.bp.GA.C.shoot2<-lm( GA  ~ FERTILIZATION * Bms, data = df.dc.central.shoot)
anova(model.bp.GA.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.central.shoot$GA, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.central.shoot$GA)

hist(residuals(model.bp.GA.C.shoot))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.GA.C.shoot )
MSerror<-deviance(model.bp.GA.C.shoot)/df
comparison <- LSD.test(model.bp.GA.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.GA.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.dc.central.shoot$GA_log10<-log10(df.dc.central.shoot$GA)
df.dc.central.shoot

bp.GA_log10.C.shoot<-ggplot(df.dc.central.shoot,aes(TREATMENT, GA_log10 ,  color=Bms)) + geom_boxplot()
bp.GA_log10.C.shoot

hist(df.dc.central.shoot$GA_log10)
shapiro.test(df.dc.central.shoot$GA_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.central.shoot$GA_log10, rnorm(n=length(df.dc.central.shoot$GA_log10), mean=mean(df.dc.central.shoot$GA_log10), sd=sd(df.dc.central.shoot$GA_log10)))

# create ANOVA model
model.bp.GA_log10.C.shoot<-lm( GA_log10  ~ TREATMENT, data = df.dc.central.shoot)
anova(model.bp.GA_log10.C.shoot)

model.bp.GA_log10.C.shoot2<-lm( GA_log10  ~ FERTILIZATION * Bms, data = df.dc.central.shoot)
anova(model.bp.GA_log10.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.central.shoot$GA_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.central.shoot$GA_log10)

hist(residuals(model.bp.GA_log10.C.shoot))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.GA_log10.C.shoot )
MSerror<-deviance(model.bp.GA_log10.C.shoot)/df
comparison <- LSD.test(model.bp.GA_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.GA_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison


##### IAA Central plots - shoots #####
df.IAA<-df.dc.central.shoot[,c(3:9,13)]
df.IAA

IAA.summary <- summarySE(df.IAA, measurevar="IAA", groupvars=c("FERTILIZATION","Bms"))
IAA.summary

IAA.summary2 <- IAA.summary
IAA.summary2$FERTILIZATION <- factor(IAA.summary2$FERTILIZATION)

val<-c("b","a","b", "a")

bp.IAA.C.shoot<-ggplot(data = IAA.summary2) +
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
    size = rel(3),
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
    axis.text.x = element_text(size=rel(1.5),color="#000000"),
    axis.text.y = element_text(size=rel(1),color="#000000"),
    axis.title.y.left = element_text(colour = 'black', size = rel(1)),
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
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.9, 'cm'))

bp.IAA.C.shoot
#   


hist(df.dc.central.shoot$IAA)
shapiro.test(df.dc.central.shoot$IAA)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.central.shoot$IAA, rnorm(n=length(df.dc.central.shoot$IAA), mean=mean(df.dc.central.shoot$IAA), sd=sd(df.dc.central.shoot$IAA)))

# create ANOVA model
model.bp.IAA.C.shoot<-lm( IAA  ~ TREATMENT, data = df.dc.central.shoot)
anova(model.bp.IAA.C.shoot)

# staristical model
model.bp.IAA.C.shoot2<-lm( IAA  ~ FERTILIZATION * Bms, data = df.dc.central.shoot)
anova(model.bp.IAA.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.central.shoot$IAA, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.central.shoot$IAA)

hist(residuals(model.bp.IAA.C.shoot))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.IAA.C.shoot )
MSerror<-deviance(model.bp.IAA.C.shoot)/df
comparison <- LSD.test(model.bp.IAA.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.IAA.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.dc.central.shoot$IAA_log10<-log10(df.dc.central.shoot$IAA)
df.dc.central.shoot

bp.IAA_log10.C.shoot<-ggplot(df.dc.central.shoot,aes(TREATMENT, IAA_log10 ,  color=Bms)) + geom_boxplot()
bp.IAA_log10.C.shoot

hist(df.dc.central.shoot$IAA_log10)
shapiro.test(df.dc.central.shoot$IAA_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.central.shoot$IAA_log10, rnorm(n=length(df.dc.central.shoot$IAA_log10), mean=mean(df.dc.central.shoot$IAA_log10), sd=sd(df.dc.central.shoot$IAA_log10)))

# create ANOVA model
model.bp.IAA_log10.C.shoot<-lm( IAA_log10  ~ TREATMENT, data = df.dc.central.shoot)
anova(model.bp.IAA_log10.C.shoot)

# statistical model
model.bp.IAA_log10.C.shoot2<-lm( IAA_log10  ~ FERTILIZATION * Bms, data = df.dc.central.shoot)
anova(model.bp.IAA_log10.C.shoot2)


#let's check if the data are now normal
# Density plot
ggdensity(df.dc.central.shoot$IAA_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.central.shoot$IAA_log10)

hist(residuals(model.bp.IAA_log10.C.shoot))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.IAA_log10.C.shoot )
MSerror<-deviance(model.bp.IAA_log10.C.shoot)/df
comparison <- LSD.test(model.bp.IAA_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.IAA_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison


##### JA Central plots - shoots #####
df.JA<-df.dc.central.shoot[,c(3:9,14)]
df.JA

JA.summary <- summarySE(df.JA, measurevar="JA", groupvars=c("FERTILIZATION","Bms"))
JA.summary

JA.summary2 <- JA.summary
JA.summary2$FERTILIZATION <- factor(JA.summary2$FERTILIZATION)

val<-c("a","bc","ab", "c")

bp.JA.C.shoot<-ggplot(data = JA.summary2) +
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
      label = val,
      group = Bms
    ),
    size = rel(3),
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
    axis.text.x = element_text(size=rel(1.5),color="#000000"),
    axis.text.y = element_text(size=rel(1),color="#000000"),
    axis.title.y.left = element_text(colour = 'black', size = rel(1)),
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
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.9, 'cm'))

bp.JA.C.shoot

#   
#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.central.shoot$JA, rnorm(n=length(df.dc.central.shoot$JA), mean=mean(df.dc.central.shoot$JA), sd=sd(df.dc.central.shoot$JA)))

# create ANOVA model
model.bp.JA.C.shoot<-lm( JA  ~ TREATMENT, data = df.dc.central.shoot)
anova(model.bp.JA.C.shoot)

# staristical model
model.bp.JA.C.shoot2<-lm( JA  ~ FERTILIZATION * Bms, data = df.dc.central.shoot)
anova(model.bp.JA.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.central.shoot$JA, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.central.shoot$JA)

hist(residuals(model.bp.JA.C.shoot))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.JA.C.shoot )
MSerror<-deviance(model.bp.JA.C.shoot)/df
comparison <- LSD.test(model.bp.JA.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.JA.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.dc.central.shoot$JA_log10<-log10(df.dc.central.shoot$JA)
df.dc.central.shoot

bp.JA_log10.C.shoot<-ggplot(df.dc.central.shoot,aes(TREATMENT, JA_log10 ,  color=Bms)) + geom_boxplot()
bp.JA_log10.C.shoot

hist(df.dc.central.shoot$JA_log10)
shapiro.test(df.dc.central.shoot$JA_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.central.shoot$JA_log10, rnorm(n=length(df.dc.central.shoot$JA_log10), mean=mean(df.dc.central.shoot$JA_log10), sd=sd(df.dc.central.shoot$JA_log10)))

# create ANOVA model
model.bp.JA_log10.C.shoot<-lm( JA_log10  ~ TREATMENT, data = df.dc.central.shoot)
anova(model.bp.JA_log10.C.shoot)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.central.shoot$JA_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.central.shoot$JA_log10)

hist(residuals(model.bp.JA_log10.C.shoot))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.JA_log10.C.shoot )
MSerror<-deviance(model.bp.JA_log10.C.shoot)/df
comparison <- LSD.test(model.bp.JA_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.JA_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison



##### SA Central plots - shoots #####
df.SA<-df.dc.central.shoot[,c(3:9,15)]
df.SA

SA.summary <- summarySE(df.SA, measurevar="SA", groupvars=c("FERTILIZATION","Bms"))
SA.summary

SA.summary2 <- SA.summary
SA.summary2$FERTILIZATION <- factor(SA.summary2$FERTILIZATION)

# val<-c("a","bc","ab", "c")

bp.SA.C.shoot<-ggplot(data = SA.summary2) +
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
#     size = rel(3),
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
    axis.text.x = element_text(size=rel(1.5),color="#000000"),
    axis.text.y = element_text(size=rel(1),color="#000000"),
    axis.title.y.left = element_text(colour = 'black', size = rel(1)),
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
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.9, 'cm'))

bp.SA.C.shoot



hist(df.dc.central.shoot$SA)
shapiro.test(df.dc.central.shoot$SA)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.central.shoot$SA, rnorm(n=length(df.dc.central.shoot$SA), mean=mean(df.dc.central.shoot$SA), sd=sd(df.dc.central.shoot$SA)))

# create ANOVA model
model.bp.SA.C.shoot<-lm( SA  ~ TREATMENT, data = df.dc.central.shoot)
anova(model.bp.SA.C.shoot)

# staristical model
model.bp.SA.C.shoot2<-lm( SA  ~ FERTILIZATION * Bms, data = df.dc.central.shoot)
anova(model.bp.SA.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.central.shoot$SA, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.central.shoot$SA)

hist(residuals(model.bp.SA.C.shoot))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.SA.C.shoot )
MSerror<-deviance(model.bp.SA.C.shoot)/df
comparison <- LSD.test(model.bp.SA.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.SA.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.dc.central.shoot$SA_log10<-log10(df.dc.central.shoot$SA)
df.dc.central.shoot

bp.SA_log10.C.shoot<-ggplot(df.dc.central.shoot,aes(TREATMENT, SA_log10 ,  color=Bms)) + geom_boxplot()
bp.SA_log10.C.shoot

hist(df.dc.central.shoot$SA_log10)
shapiro.test(df.dc.central.shoot$SA_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.central.shoot$SA_log10, rnorm(n=length(df.dc.central.shoot$SA_log10), mean=mean(df.dc.central.shoot$SA_log10), sd=sd(df.dc.central.shoot$SA_log10)))

# create ANOVA model
model.bp.SA_log10.C.shoot<-lm( SA_log10  ~ TREATMENT, data = df.dc.central.shoot)
anova(model.bp.SA_log10.C.shoot)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.central.shoot$SA_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.central.shoot$SA_log10)

hist(residuals(model.bp.SA_log10.C.shoot))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.SA_log10.C.shoot )
MSerror<-deviance(model.bp.SA_log10.C.shoot)/df
comparison <- LSD.test(model.bp.SA_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.SA_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison




## ggarrange ####
figure.hormones<-ggarrange( bp.IAA.C.shoot, bp.CK.C.shoot, bp.GA.C.shoot, bp.ABA.C.shoot, bp.JA.C.shoot, bp.SA.C.shoot, 
                           labels = c("Indole-3-acetic acid", "    Cytokinin", "    Gibberellins", "Abscisic Acid", " Jasmonic acid", " Salicylic acid"),
                           ncol = 3, nrow = 2,   font.label = list(size = 18, color = "black", face = "bold"))
figure.hormones

# figure.hormones<-ggarrange(bp.ABA.C.shoot, bp.CK.C.shoot, bp.GA.C.shoot, bp.IAA.C.shoot, bp.JA.C.shoot, bp.SA.C.shoot, 
#                            labels = c("Abscisic Acid", "    Citokinin", "    Gibberellins","Indole-3-acetic acid", " Jasmonic acid", " Salicylic acid"),
#                            ncol = 3, nrow = 2,   font.label = list(size = 18, color = "black", face = "bold"))
# figure.hormones
# 
# 
# figure.hormones<-ggarrange(bp.ABA.C.shoot, bp.CK.C.shoot, bp.GA.C.shoot, bp.IAA.C.shoot, bp.JA.C.shoot, bp.SA.C.shoot, 
#                            labels = c("Abscisic Acid", "    Citokinin", "    Gibberellins","Indole-3-acetic acid", "   Jasmonic acid", " Salicylic acid"),
#                            ncol = 3, nrow = 2,   font.label = list(size = 20, color = "black", face = "bold"))
# figure.hormones
# 
# 
# #####
# figure.ABA.GA<-ggarrange(bp.ABA.C.shoot, bp.GA.C.shoot, 
#           labels = c("Abscisic Acid", "Gibberellins"),
#           ncol = 2, nrow = 1,   font.label = list(size = 18, color = "black", face = "bold"))
# 
# figure.ABA.GA
# 
# figure.ABA.GA.JA<-ggarrange(bp.ABA.C.shoot, bp.GA.C.shoot, bp.JA.C.shoot,
#                          labels = c("Abscisic Acid", "Gibberellins", "Jasmonic acid"),
#                          ncol = 3, nrow = 1,   font.label = list(size = 18, color = "black", face = "bold"))
# 
# figure.ABA.GA.JA


