
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

# install.packages(c("nlme",lib=.libPaths()[2]))
# install.packages("MuMIn")


df.dc.RW.root<-read.table("Root_traits_LTE1_2020_RW_sorted_biomass.txt",h=T)
df.dc.RW.root

str(df.dc.RW.root)
head(df.dc.RW.root)
tail(df.dc.RW.root)


df.dc.RW.root$TREATMENT = factor(df.dc.RW.root$TREATMENT, 
                                   levels=c("Extensive_Ctrl","Extensive_BMc", "Intensive_Ctrl","Intensive_BMc"))

##### TRL RW plots - root #####
bp.TRL.C.root<-ggplot(df.dc.RW.root,aes(FERTILIZATION, TRL ,  fill=Bms)) + geom_boxplot()
bp.TRL.C.root

bp.TRL.C.root<-ggplot(df.dc.RW.root,aes(FERTILIZATION, TRL ,  fill=Bms)) + geom_boxplot()+
scale_fill_manual(values=c("#146eb4",  "#ff9900", "#146eb4" , "#ff9900"))+
  theme(legend.position = "none",
        plot.title = element_text(color = "black", size = 10, face = "bold"),
        axis.title.x = element_text(color = "white", size = 10),
        axis.title.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(
          color = "black",
          size = 10,
          angle = 0,
          hjust = .5,
          vjust = .5,
          face = "plain"
        ),
        axis.text.y = element_text(
          color = "black",
          size = 10,
          angle = 0,
          hjust = 1,
          vjust = 0,
          face = "plain"
        ),
        legend.title=element_text(size = 10), 
        legend.text=element_text(size = 10))+labs(y = "Total root lenght (cm)")

bp.TRL.C.root

hist(df.dc.RW.root$TRL)
shapiro.test(df.dc.RW.root$TRL)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.RW.root$TRL, rnorm(n=length(df.dc.RW.root$TRL), mean=mean(df.dc.RW.root$TRL), sd=sd(df.dc.RW.root$TRL)))

# create ANOVA model
model.bp.TRL.C.root<-lm( TRL  ~ TREATMENT, data = df.dc.RW.root)
anova(model.bp.TRL.C.root)

# statistical model
model.bp.TRL.C.root2<-lm( TRL  ~ FERTILIZATION * Bms, data = df.dc.RW.root)
anova(model.bp.TRL.C.root2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.RW.root$TRL, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.RW.root$TRL)

hist(residuals(model.bp.TRL.C.root))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.TRL.C.root )
MSerror<-deviance(model.bp.TRL.C.root)/df
comparison <- LSD.test(model.bp.TRL.C.root,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.TRL.C.root,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.dc.RW.root$TRL_log10<-log10(df.dc.RW.root$TRL)
df.dc.RW.root

bp.TRL_log10.C.root<-ggplot(df.dc.RW.root,aes(TREATMENT, TRL_log10 ,  color=Bms)) + geom_boxplot()
bp.TRL_log10.C.root

hist(df.dc.RW.root$TRL_log10)
shapiro.test(df.dc.RW.root$TRL_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.RW.root$TRL_log10, rnorm(n=length(df.dc.RW.root$TRL_log10), mean=mean(df.dc.RW.root$TRL_log10), sd=sd(df.dc.RW.root$TRL_log10)))

# create ANOVA model
model.bp.TRL_log10.C.root<-lm( TRL_log10  ~ TREATMENT, data = df.dc.RW.root)
anova(model.bp.TRL_log10.C.root)

# staristical model
model.bp.TRL_log10.C.root2<-lm( TRL_log10  ~ FERTILIZATION * Bms, data = df.dc.RW.root)
anova(model.bp.TRL_log10.C.root2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.RW.root$TRL_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.RW.root$TRL_log10)

hist(residuals(model.bp.TRL_log10.C.root))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.TRL_log10.C.root )
MSerror<-deviance(model.bp.TRL_log10.C.root)/df
comparison <- LSD.test(model.bp.TRL_log10.C.root,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.TRL_log10.C.root,"TREATMENT", df, MSerror, group=T)
comparison




##### AVD RW plots - root #####
bp.AVD.C.root<-ggplot(df.dc.RW.root,aes(TREATMENT, AVD ,  fill=Bms)) + geom_boxplot()
bp.AVD.C.root

bp.AVD.C.root<-ggplot(df.dc.RW.root,aes(TREATMENT, AVD ,  fill=Bms)) + geom_boxplot()+
scale_fill_manual(values=c("#146eb4",  "#ff9900", "#146eb4" , "#ff9900"))+
  theme(legend.position = "none",
        plot.title = element_text(color = "black", size = 10, face = "bold"),
        axis.title.x = element_text(color = "white", size = 10),
        axis.title.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(
          color = "black",
          size = 10,
          angle = 0,
          hjust = .5,
          vjust = .5,
          face = "plain"
        ),
        axis.text.y = element_text(
          color = "black",
          size = 10,
          angle = 0,
          hjust = 1,
          vjust = 0,
          face = "plain"
        ),
        legend.title=element_text(size = 10), 
        legend.text=element_text(size = 10))+labs(y = "Root diameter (mm)")

bp.AVD.C.root

hist(df.dc.RW.root$AVD)
shapiro.test(df.dc.RW.root$AVD)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.RW.root$AVD, rnorm(n=length(df.dc.RW.root$AVD), mean=mean(df.dc.RW.root$AVD), sd=sd(df.dc.RW.root$AVD)))

# create ANOVA model
model.bp.AVD.C.root<-lm( AVD  ~ TREATMENT, data = df.dc.RW.root)
anova(model.bp.AVD.C.root)

# statistical model
model.bp.AVD.C.root2<-lm( AVD  ~ FERTILIZATION * Bms, data = df.dc.RW.root)
anova(model.bp.AVD.C.root2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.RW.root$AVD, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.RW.root$AVD)

hist(residuals(model.bp.AVD.C.root))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.AVD.C.root )
MSerror<-deviance(model.bp.AVD.C.root)/df
comparison <- LSD.test(model.bp.AVD.C.root,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.AVD.C.root,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.dc.RW.root$AVD_log10<-log10(df.dc.RW.root$AVD)
df.dc.RW.root

bp.AVD_log10.C.root<-ggplot(df.dc.RW.root,aes(TREATMENT, AVD_log10 ,  color=Bms)) + geom_boxplot()
bp.AVD_log10.C.root

hist(df.dc.RW.root$AVD_log10)
shapiro.test(df.dc.RW.root$AVD_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.RW.root$AVD_log10, rnorm(n=length(df.dc.RW.root$AVD_log10), mean=mean(df.dc.RW.root$AVD_log10), sd=sd(df.dc.RW.root$AVD_log10)))

# create ANOVA model
model.bp.AVD_log10.C.root<-lm( AVD_log10  ~ TREATMENT, data = df.dc.RW.root)
anova(model.bp.AVD_log10.C.root)

# staristical model
model.bp.AVD_log10.C.root2<-lm( AVD_log10  ~ FERTILIZATION * Bms, data = df.dc.RW.root)
anova(model.bp.AVD_log10.C.root2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.RW.root$AVD_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.RW.root$AVD_log10)

hist(residuals(model.bp.AVD_log10.C.root))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.AVD_log10.C.root )
MSerror<-deviance(model.bp.AVD_log10.C.root)/df
comparison <- LSD.test(model.bp.AVD_log10.C.root,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.AVD_log10.C.root,"TREATMENT", df, MSerror, group=T)
comparison




##### AMF RW plots - root #####
bp.AMF.C.root<-ggplot(df.dc.RW.root,aes(TREATMENT, AMF ,  fill=Bms)) + geom_boxplot()
bp.AMF.C.root

bp.AMF.C.root<-ggplot(df.dc.RW.root,aes(TREATMENT, AMF ,  fill=Bms)) + geom_boxplot()+
scale_fill_manual(values=c("#146eb4",  "#ff9900", "#146eb4" , "#ff9900"))+
  theme(legend.position = "none",
        plot.title = element_text(color = "black", size = 10, face = "bold"),
        axis.title.x = element_text(color = "white", size = 10),
        axis.title.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(
          color = "black",
          size = 10,
          angle = 0,
          hjust = .5,
          vjust = .5,
          face = "plain"
        ),
        axis.text.y = element_text(
          color = "black",
          size = 10,
          angle = 0,
          hjust = 1,
          vjust = 0,
          face = "plain"
        ),
        legend.title=element_text(size = 10), 
        legend.text=element_text(size = 10))+labs(y = "AMF root colonization (%)")
bp.AMF.C.root

hist(df.dc.RW.root$AMF)
shapiro.test(df.dc.RW.root$AMF)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.RW.root$AMF, rnorm(n=length(df.dc.RW.root$AMF), mean=mean(df.dc.RW.root$AMF), sd=sd(df.dc.RW.root$AMF)))

# create ANOVA model
model.bp.AMF.C.root<-lm( AMF  ~ TREATMENT, data = df.dc.RW.root)
anova(model.bp.AMF.C.root)

# statistical model
model.bp.AMF.C.root2<-lm( AMF  ~ FERTILIZATION * Bms, data = df.dc.RW.root)
anova(model.bp.AMF.C.root2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.RW.root$AMF, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.RW.root$AMF)

hist(residuals(model.bp.AMF.C.root))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.AMF.C.root )
MSerror<-deviance(model.bp.AMF.C.root)/df
comparison <- LSD.test(model.bp.AMF.C.root,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.AMF.C.root,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.dc.RW.root$AMF_log10<-log10(df.dc.RW.root$AMF)
df.dc.RW.root

bp.AMF_log10.C.root<-ggplot(df.dc.RW.root,aes(TREATMENT, AMF_log10 ,  color=Bms)) + geom_boxplot()
bp.AMF_log10.C.root

hist(df.dc.RW.root$AMF_log10)
shapiro.test(df.dc.RW.root$AMF_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.RW.root$AMF_log10, rnorm(n=length(df.dc.RW.root$AMF_log10), mean=mean(df.dc.RW.root$AMF_log10), sd=sd(df.dc.RW.root$AMF_log10)))

# create ANOVA model
model.bp.AMF_log10.C.root<-lm( AMF_log10  ~ TREATMENT, data = df.dc.RW.root)
anova(model.bp.AMF_log10.C.root)

# staristical model
model.bp.AMF_log10.C.root2<-lm( AMF_log10  ~ FERTILIZATION * Bms, data = df.dc.RW.root)
anova(model.bp.AMF_log10.C.root2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.RW.root$AMF_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.RW.root$AMF_log10)

hist(residuals(model.bp.AMF_log10.C.root))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.AMF_log10.C.root )
MSerror<-deviance(model.bp.AMF_log10.C.root)/df
comparison <- LSD.test(model.bp.AMF_log10.C.root,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.AMF_log10.C.root,"TREATMENT", df, MSerror, group=T)
comparison




##### root_depth RW plots - root #####
bp.root_depth.C.root<-ggplot(df.dc.RW.root,aes(TREATMENT, root_depth ,  fill=Bms)) + geom_boxplot()
bp.root_depth.C.root

bp.root_depth.C.root<-ggplot(df.dc.RW.root,aes(TREATMENT, root_depth ,  fill=Bms)) + geom_boxplot()+
scale_fill_manual(values=c("#146eb4",  "#ff9900", "#146eb4" , "#ff9900"))+
  theme(legend.position = "none",
        plot.title = element_text(color = "black", size = 10, face = "bold"),
        axis.title.x = element_text(color = "white", size = 10),
        axis.title.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(
          color = "black",
          size = 10,
          angle = 0,
          hjust = .5,
          vjust = .5,
          face = "plain"
        ),
        axis.text.y = element_text(
          color = "black",
          size = 10,
          angle = 0,
          hjust = 1,
          vjust = 0,
          face = "plain"
        ),
        legend.title=element_text(size = 10), 
        legend.text=element_text(size = 10))+labs(y = "Root depth (cm)")

bp.root_depth.C.root

hist(df.dc.RW.root$root_depth)

shapiro.test(df.dc.RW.root$root_depth)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.RW.root$root_depth, rnorm(n=length(df.dc.RW.root$root_depth), mean=mean(df.dc.RW.root$root_depth), sd=sd(df.dc.RW.root$root_depth)))

# create ANOVA model
model.bp.root_depth.C.root<-lm( root_depth  ~ TREATMENT, data = df.dc.RW.root)
anova(model.bp.root_depth.C.root)

# statistical model
model.bp.root_depth.C.root2<-lm( root_depth  ~ FERTILIZATION * Bms, data = df.dc.RW.root)
anova(model.bp.root_depth.C.root2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.RW.root$root_depth, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.RW.root$root_depth)

hist(residuals(model.bp.root_depth.C.root))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.root_depth.C.root )
MSerror<-deviance(model.bp.root_depth.C.root)/df
comparison <- LSD.test(model.bp.root_depth.C.root,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.root_depth.C.root,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.dc.RW.root$root_depth_log10<-log10(df.dc.RW.root$root_depth)
df.dc.RW.root

bp.root_depth_log10.C.root<-ggplot(df.dc.RW.root,aes(TREATMENT, root_depth_log10 ,  color=Bms)) + geom_boxplot()
bp.root_depth_log10.C.root

hist(df.dc.RW.root$root_depth_log10)
shapiro.test(df.dc.RW.root$root_depth_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.RW.root$root_depth_log10, rnorm(n=length(df.dc.RW.root$root_depth_log10), mean=mean(df.dc.RW.root$root_depth_log10), sd=sd(df.dc.RW.root$root_depth_log10)))

# create ANOVA model
model.bp.root_depth_log10.C.root<-lm( root_depth_log10  ~ TREATMENT, data = df.dc.RW.root)
anova(model.bp.root_depth_log10.C.root)

# staristical model
model.bp.root_depth_log10.C.root2<-lm( root_depth_log10  ~ FERTILIZATION * Bms, data = df.dc.RW.root)
anova(model.bp.root_depth_log10.C.root2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.RW.root$root_depth_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.RW.root$root_depth_log10)

hist(residuals(model.bp.root_depth_log10.C.root))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.root_depth_log10.C.root )
MSerror<-deviance(model.bp.root_depth_log10.C.root)/df
comparison <- LSD.test(model.bp.root_depth_log10.C.root,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.root_depth_log10.C.root,"TREATMENT", df, MSerror, group=T)
comparison






##### RHL RW plots - root #####
bp.RHL.C.root<-ggplot(df.dc.RW.root,aes(TREATMENT, RHL,  fill=Bms)) + geom_boxplot()
bp.RHL.C.root

bp.RHL.C.root<-ggplot(df.dc.RW.root,aes(TREATMENT, RHL,  fill=Bms)) + geom_boxplot()+
  scale_fill_manual(values=c("#146eb4",  "#ff9900", "#146eb4" , "#ff9900"))+
  theme(legend.position = "none",
        plot.title = element_text(color = "black", size = 10, face = "bold"),
        axis.title.x = element_text(color = "white", size = 10),
        axis.title.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(
          color = "black",
          size = 10,
          angle = 0,
          hjust = .5,
          vjust = .5,
          face = "plain"
        ),
        axis.text.y = element_text(
          color = "black",
          size = 10,
          angle = 0,
          hjust = 1,
          vjust = 0,
          face = "plain"
        ),
        legend.title=element_text(size = 10), 
        legend.text=element_text(size = 10))+labs(y = "Root hairs size (mm)")

bp.RHL.C.root

hist(df.dc.RW.root$RHL)
shapiro.test(df.dc.RW.root$RHL)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.RW.root$RHL, rnorm(n=length(df.dc.RW.root$RHL), mean=mean(df.dc.RW.root$RHL), sd=sd(df.dc.RW.root$RHL)))

# create ANOVA model
model.bp.RHL.C.root<-lm( RHL  ~ TREATMENT, data = df.dc.RW.root)
anova(model.bp.RHL.C.root)

# statistical model
model.bp.RHL.C.root2<-lm( RHL  ~ FERTILIZATION * Bms, data = df.dc.RW.root)
anova(model.bp.RHL.C.root2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.RW.root$RHL, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.RW.root$RHL)

hist(residuals(model.bp.RHL.C.root))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.RHL.C.root )
MSerror<-deviance(model.bp.RHL.C.root)/df
comparison <- LSD.test(model.bp.RHL.C.root,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.RHL.C.root,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.dc.RW.root$RHL_log10<-log10(100*(df.dc.RW.root$RHL))
df.dc.RW.root

bp.RHL_log10.C.root<-ggplot(df.dc.RW.root,aes(TREATMENT, RHL_log10 ,  color=Bms)) + geom_boxplot()
bp.RHL_log10.C.root

hist(df.dc.RW.root$RHL_log10)
shapiro.test(df.dc.RW.root$RHL_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.dc.RW.root$RHL_log10, rnorm(n=length(df.dc.RW.root$RHL_log10), mean=mean(df.dc.RW.root$RHL_log10), sd=sd(df.dc.RW.root$RHL_log10)))

# create ANOVA model
model.bp.RHL_log10.C.root<-lm( RHL_log10  ~ TREATMENT, data = df.dc.RW.root)
anova(model.bp.RHL_log10.C.root)

# staristical model
model.bp.RHL_log10.C.root2<-lm( RHL_log10  ~ FERTILIZATION * Bms, data = df.dc.RW.root)
anova(model.bp.RHL_log10.C.root2)

#let's check if the data are now normal
# Density plot
ggdensity(df.dc.RW.root$RHL_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.dc.RW.root$RHL_log10)

hist(residuals(model.bp.RHL_log10.C.root))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.RHL_log10.C.root )
MSerror<-deviance(model.bp.RHL_log10.C.root)/df
comparison <- LSD.test(model.bp.RHL_log10.C.root,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.RHL_log10.C.root,"TREATMENT", df, MSerror, group=T)
comparison




##### ggarrange ####

figure.Root.traits<-ggarrange(bp.TRL.C.root, bp.AVD.C.root, bp.AMF.C.root,bp.root_depth.C.root,bp.RHL.C.root,   
                          labels = c("Total root lenght", "Root diameter", "AMF root colonization", "Root depth", "Root hair size"),
                          ncol = 5, nrow = 1,   font.label = list(size = 10, color = "black", face = "bold"))

figure.Root.traits


