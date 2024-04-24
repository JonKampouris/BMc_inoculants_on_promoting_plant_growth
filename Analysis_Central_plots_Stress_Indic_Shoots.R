
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
# subsetting central shoot dataset
df.StrI.central.shoot<-subset(df.StrI.central, LOCATION =="shoot")

# # subsetting central root dataset
# df.StrI.central.root<-subset(df.StrI.central, LOCATION =="root")


##### APX Central plots - Shoots #####
bp.APX.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, APX ,  fill=Bms)) + geom_boxplot()
bp.APX.C.shoot

bp.APX.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, APX ,  fill=Bms))+ 
geom_boxplot()+scale_fill_manual(values=c("#146eb4",  "#ff9900", "#146eb4" , "#ff9900"))+
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
        legend.text=element_text(size = 10))+labs(y = "APX (Unit g-1 FW)")
  
bp.APX.C.shoot

hist(df.StrI.central.shoot$APX)
shapiro.test(df.StrI.central.shoot$APX)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.StrI.central.shoot$APX, rnorm(n=length(df.StrI.central.shoot$APX), mean=mean(df.StrI.central.shoot$APX), sd=sd(df.StrI.central.shoot$APX)))

# create ANOVA model
model.bp.APX.C.shoot<-lm( APX  ~ TREATMENT, data = df.StrI.central.shoot)
anova(model.bp.APX.C.shoot)

# staristical model
model.bp.APX.C.shoot2<-lm( APX  ~ FERTILIZATION * Bms, data = df.StrI.central.shoot)
anova(model.bp.APX.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.StrI.central.shoot$APX, fill = "lightgray")
# QQ plot
ggqqplot(df.StrI.central.shoot$APX)

hist(residuals(model.bp.APX.C.shoot))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.APX.C.shoot )
MSerror<-deviance(model.bp.APX.C.shoot)/df
comparison <- LSD.test(model.bp.APX.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.APX.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.StrI.central.shoot$APX_log10<-log10(df.StrI.central.shoot$APX)
df.StrI.central.shoot

bp.APX_log10.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, APX_log10 ,  color=Bms)) + geom_boxplot()
bp.APX_log10.C.shoot

hist(df.StrI.central.shoot$APX_log10)
shapiro.test(df.StrI.central.shoot$APX_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.StrI.central.shoot$APX_log10, rnorm(n=length(df.StrI.central.shoot$APX_log10), mean=mean(df.StrI.central.shoot$APX_log10), sd=sd(df.StrI.central.shoot$APX_log10)))

# create ANOVA model
model.bp.APX_log10.C.shoot<-lm( APX_log10  ~ TREATMENT, data = df.StrI.central.shoot)
anova(model.bp.APX_log10.C.shoot)

#let's check if the data are now normal
# Density plot
ggdensity(df.StrI.central.shoot$APX_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.StrI.central.shoot$APX_log10)

hist(residuals(model.bp.APX_log10.C.shoot))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.APX_log10.C.shoot )
MSerror<-deviance(model.bp.APX_log10.C.shoot)/df
comparison <- LSD.test(model.bp.APX_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.APX_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison



##### GB Central plots - Shoots #####
bp.GB.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, GB ,  fill=Bms)) + geom_boxplot()
bp.GB.C.shoot

bp.GB.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, GB ,  fill=Bms)) + geom_boxplot()+
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
        legend.text=element_text(size = 10))+labs(y = "GB (mg g-1 FW)")
  
bp.GB.C.shoot

hist(df.StrI.central.shoot$GB)
shapiro.test(df.StrI.central.shoot$GB)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.StrI.central.shoot$GB, rnorm(n=length(df.StrI.central.shoot$GB), mean=mean(df.StrI.central.shoot$GB), sd=sd(df.StrI.central.shoot$GB)))

# create ANOVA model
model.bp.GB.C.shoot<-lm( GB  ~ TREATMENT, data = df.StrI.central.shoot)
anova(model.bp.GB.C.shoot)

# staristical model
model.bp.GB.C.shoot2<-lm( GB  ~ FERTILIZATION * Bms, data = df.StrI.central.shoot)
anova(model.bp.GB.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.StrI.central.shoot$GB, fill = "lightgray")
# QQ plot
ggqqplot(df.StrI.central.shoot$GB)

hist(residuals(model.bp.GB.C.shoot))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.GB.C.shoot )
MSerror<-deviance(model.bp.GB.C.shoot)/df
comparison <- LSD.test(model.bp.GB.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.GB.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.StrI.central.shoot$GB_log10<-log10(df.StrI.central.shoot$GB)
df.StrI.central.shoot

bp.GB_log10.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, GB_log10 ,  color=Bms)) + geom_boxplot()
bp.GB_log10.C.shoot

hist(df.StrI.central.shoot$GB_log10)
shapiro.test(df.StrI.central.shoot$GB_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.StrI.central.shoot$GB_log10, rnorm(n=length(df.StrI.central.shoot$GB_log10), mean=mean(df.StrI.central.shoot$GB_log10), sd=sd(df.StrI.central.shoot$GB_log10)))

# create ANOVA model
model.bp.GB_log10.C.shoot<-lm( GB_log10  ~ TREATMENT, data = df.StrI.central.shoot)
anova(model.bp.GB_log10.C.shoot)

#let's check if the data are now normal
# Density plot
ggdensity(df.StrI.central.shoot$GB_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.StrI.central.shoot$GB_log10)

hist(residuals(model.bp.GB_log10.C.shoot))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.GB_log10.C.shoot )
MSerror<-deviance(model.bp.GB_log10.C.shoot)/df
comparison <- LSD.test(model.bp.GB_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.GB_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison


##### H2O2 Central plots - Shoots #####
bp.H2O2.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, H2O2 ,  fill=Bms)) + geom_boxplot()
bp.H2O2.C.shoot

bp.H2O2.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, H2O2 ,  fill=Bms)) + geom_boxplot()+
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
        legend.text=element_text(size = 10))+labs(y = "H2O2 (µmol g-1 FW)")
bp.H2O2.C.shoot

hist(df.StrI.central.shoot$H2O2)
shapiro.test(df.StrI.central.shoot$H2O2)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.StrI.central.shoot$H2O2, rnorm(n=length(df.StrI.central.shoot$H2O2), mean=mean(df.StrI.central.shoot$H2O2), sd=sd(df.StrI.central.shoot$H2O2)))

# create ANOVA model
model.bp.H2O2.C.shoot<-lm( H2O2  ~ TREATMENT, data = df.StrI.central.shoot)
anova(model.bp.H2O2.C.shoot)

# staristical model
model.bp.H2O2.C.shoot2<-lm( H2O2  ~ FERTILIZATION * Bms, data = df.StrI.central.shoot)
anova(model.bp.H2O2.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.StrI.central.shoot$H2O2, fill = "lightgray")
# QQ plot
ggqqplot(df.StrI.central.shoot$H2O2)

hist(residuals(model.bp.H2O2.C.shoot))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.H2O2.C.shoot )
MSerror<-deviance(model.bp.H2O2.C.shoot)/df
comparison <- LSD.test(model.bp.H2O2.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.H2O2.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.StrI.central.shoot$H2O2_log10<-log10(df.StrI.central.shoot$H2O2)
df.StrI.central.shoot

bp.H2O2_log10.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, H2O2_log10 ,  color=Bms)) + geom_boxplot()
bp.H2O2_log10.C.shoot

hist(df.StrI.central.shoot$H2O2_log10)
shapiro.test(df.StrI.central.shoot$H2O2_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.StrI.central.shoot$H2O2_log10, rnorm(n=length(df.StrI.central.shoot$H2O2_log10), mean=mean(df.StrI.central.shoot$H2O2_log10), sd=sd(df.StrI.central.shoot$H2O2_log10)))

# create ANOVA model
model.bp.H2O2_log10.C.shoot<-lm( H2O2_log10  ~ TREATMENT, data = df.StrI.central.shoot)
anova(model.bp.H2O2_log10.C.shoot)

#let's check if the data are now normal
# Density plot
ggdensity(df.StrI.central.shoot$H2O2_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.StrI.central.shoot$H2O2_log10)

hist(residuals(model.bp.H2O2_log10.C.shoot))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.H2O2_log10.C.shoot )
MSerror<-deviance(model.bp.H2O2_log10.C.shoot)/df
comparison <- LSD.test(model.bp.H2O2_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.H2O2_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison




##### phenolics Central plots - Shoots #####
bp.phenolics.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, phenolics ,  fill=Bms)) + geom_boxplot()
bp.phenolics.C.shoot

bp.phenolics.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, phenolics ,  fill=Bms)) + geom_boxplot()+
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
        legend.text=element_text(size = 10))+labs(y = "Phenolics (mg gallic acid equivalent g-1 FW)")
bp.phenolics.C.shoot





hist(df.StrI.central.shoot$phenolics)
shapiro.test(df.StrI.central.shoot$phenolics)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.StrI.central.shoot$phenolics, rnorm(n=length(df.StrI.central.shoot$phenolics), mean=mean(df.StrI.central.shoot$phenolics), sd=sd(df.StrI.central.shoot$phenolics)))

# create ANOVA model
model.bp.phenolics.C.shoot<-lm( phenolics  ~ TREATMENT, data = df.StrI.central.shoot)
anova(model.bp.phenolics.C.shoot)

# staristical model
model.bp.phenolics.C.shoot2<-lm( phenolics  ~ FERTILIZATION * Bms, data = df.StrI.central.shoot)
anova(model.bp.phenolics.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.StrI.central.shoot$phenolics, fill = "lightgray")
# QQ plot
ggqqplot(df.StrI.central.shoot$phenolics)

hist(residuals(model.bp.phenolics.C.shoot))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.phenolics.C.shoot )
MSerror<-deviance(model.bp.phenolics.C.shoot)/df
comparison <- LSD.test(model.bp.phenolics.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.phenolics.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.StrI.central.shoot$phenolics_log10<-log10(df.StrI.central.shoot$phenolics)
df.StrI.central.shoot

bp.phenolics_log10.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, phenolics_log10 ,  color=Bms)) + geom_boxplot()
bp.phenolics_log10.C.shoot

hist(df.StrI.central.shoot$phenolics_log10)
shapiro.test(df.StrI.central.shoot$phenolics_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.StrI.central.shoot$phenolics_log10, rnorm(n=length(df.StrI.central.shoot$phenolics_log10), mean=mean(df.StrI.central.shoot$phenolics_log10), sd=sd(df.StrI.central.shoot$phenolics_log10)))

# create ANOVA model
model.bp.phenolics_log10.C.shoot<-lm( phenolics_log10  ~ TREATMENT, data = df.StrI.central.shoot)
anova(model.bp.phenolics_log10.C.shoot)

#let's check if the data are now normal
# Density plot
ggdensity(df.StrI.central.shoot$phenolics_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.StrI.central.shoot$phenolics_log10)

hist(residuals(model.bp.phenolics_log10.C.shoot))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.phenolics_log10.C.shoot )
MSerror<-deviance(model.bp.phenolics_log10.C.shoot)/df
comparison <- LSD.test(model.bp.phenolics_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.phenolics_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison



##### proline Central plots - Shoots #####
bp.proline.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, proline ,  fill=Bms)) + geom_boxplot()
bp.proline.C.shoot

bp.proline.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, proline ,  fill=Bms)) + geom_boxplot()+
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
        legend.text=element_text(size = 10))+labs(y = "Proline (mg g-1 FW)")
  
bp.proline.C.shoot

hist(df.StrI.central.shoot$proline)
shapiro.test(df.StrI.central.shoot$proline)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.StrI.central.shoot$proline, rnorm(n=length(df.StrI.central.shoot$proline), mean=mean(df.StrI.central.shoot$proline), sd=sd(df.StrI.central.shoot$proline)))

# create ANOVA model
model.bp.proline.C.shoot<-lm( proline  ~ TREATMENT, data = df.StrI.central.shoot)
anova(model.bp.proline.C.shoot)

# staristical model
model.bp.proline.C.shoot2<-lm( proline  ~ FERTILIZATION * Bms, data = df.StrI.central.shoot)
anova(model.bp.proline.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.StrI.central.shoot$proline, fill = "lightgray")
# QQ plot
ggqqplot(df.StrI.central.shoot$proline)

hist(residuals(model.bp.proline.C.shoot))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.proline.C.shoot )
MSerror<-deviance(model.bp.proline.C.shoot)/df
comparison <- LSD.test(model.bp.proline.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.proline.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.StrI.central.shoot$proline_log10<-log10(df.StrI.central.shoot$proline)
df.StrI.central.shoot

bp.proline_log10.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, proline_log10 ,  color=Bms)) + geom_boxplot()
bp.proline_log10.C.shoot

hist(df.StrI.central.shoot$proline_log10)
shapiro.test(df.StrI.central.shoot$proline_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.StrI.central.shoot$proline_log10, rnorm(n=length(df.StrI.central.shoot$proline_log10), mean=mean(df.StrI.central.shoot$proline_log10), sd=sd(df.StrI.central.shoot$proline_log10)))

# create ANOVA model
model.bp.proline_log10.C.shoot<-lm( proline_log10  ~ TREATMENT, data = df.StrI.central.shoot)
anova(model.bp.proline_log10.C.shoot)

#let's check if the data are now normal
# Density plot
ggdensity(df.StrI.central.shoot$proline_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.StrI.central.shoot$proline_log10)

hist(residuals(model.bp.proline_log10.C.shoot))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.proline_log10.C.shoot )
MSerror<-deviance(model.bp.proline_log10.C.shoot)/df
comparison <- LSD.test(model.bp.proline_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.proline_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison



##### SOD Central plots - Shoots #####
bp.SOD.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, SOD ,  fill=Bms)) + geom_boxplot()
bp.SOD.C.shoot

bp.SOD.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, SOD ,  fill=Bms)) + geom_boxplot()+
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
        legend.text=element_text(size = 10))+labs(y = "SOD (Unit g-1 FW)")
  
bp.SOD.C.shoot

hist(df.StrI.central.shoot$SOD)
shapiro.test(df.StrI.central.shoot$SOD)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.StrI.central.shoot$SOD, rnorm(n=length(df.StrI.central.shoot$SOD), mean=mean(df.StrI.central.shoot$SOD), sd=sd(df.StrI.central.shoot$SOD)))

# create ANOVA model
model.bp.SOD.C.shoot<-lm( SOD  ~ TREATMENT, data = df.StrI.central.shoot)
anova(model.bp.SOD.C.shoot)

# staristical model
model.bp.SOD.C.shoot2<-lm( SOD  ~ FERTILIZATION * Bms, data = df.StrI.central.shoot)
anova(model.bp.SOD.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.StrI.central.shoot$SOD, fill = "lightgray")
# QQ plot
ggqqplot(df.StrI.central.shoot$SOD)

hist(residuals(model.bp.SOD.C.shoot))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.SOD.C.shoot )
MSerror<-deviance(model.bp.SOD.C.shoot)/df
comparison <- LSD.test(model.bp.SOD.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.SOD.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.StrI.central.shoot$SOD_log10<-log10(df.StrI.central.shoot$SOD)
df.StrI.central.shoot

bp.SOD_log10.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, SOD_log10 ,  color=Bms)) + geom_boxplot()
bp.SOD_log10.C.shoot

hist(df.StrI.central.shoot$SOD_log10)
shapiro.test(df.StrI.central.shoot$SOD_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.StrI.central.shoot$SOD_log10, rnorm(n=length(df.StrI.central.shoot$SOD_log10), mean=mean(df.StrI.central.shoot$SOD_log10), sd=sd(df.StrI.central.shoot$SOD_log10)))

# create ANOVA model
model.bp.SOD_log10.C.shoot<-lm( SOD_log10  ~ TREATMENT, data = df.StrI.central.shoot)
anova(model.bp.SOD_log10.C.shoot)

#let's check if the data are now normal
# Density plot
ggdensity(df.StrI.central.shoot$SOD_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.StrI.central.shoot$SOD_log10)

hist(residuals(model.bp.SOD_log10.C.shoot))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.SOD_log10.C.shoot )
MSerror<-deviance(model.bp.SOD_log10.C.shoot)/df
comparison <- LSD.test(model.bp.SOD_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.SOD_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison



##### total_antioxidants Central plots - Shoots #####
bp.total_antioxidants.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, total_antioxidants ,  fill=Bms)) + geom_boxplot()
bp.total_antioxidants.C.shoot

bp.total_antioxidants.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, total_antioxidants ,  fill=Bms)) + geom_boxplot()+
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
        legend.text=element_text(size = 10))+labs(y = "Total antioxidants (%)")

bp.total_antioxidants.C.shoot

hist(df.StrI.central.shoot$total_antioxidants)
shapiro.test(df.StrI.central.shoot$total_antioxidants)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.StrI.central.shoot$total_antioxidants, rnorm(n=length(df.StrI.central.shoot$total_antioxidants), mean=mean(df.StrI.central.shoot$total_antioxidants), sd=sd(df.StrI.central.shoot$total_antioxidants)))

# create ANOVA model
model.bp.total_antioxidants.C.shoot<-lm( total_antioxidants  ~ TREATMENT, data = df.StrI.central.shoot)
anova(model.bp.total_antioxidants.C.shoot)

# staristical model
model.bp.total_antioxidants.C.shoot2<-lm( total_antioxidants  ~ FERTILIZATION * Bms, data = df.StrI.central.shoot)
anova(model.bp.total_antioxidants.C.shoot2)

#let's check if the data are now normal
# Density plot
ggdensity(df.StrI.central.shoot$total_antioxidants, fill = "lightgray")
# QQ plot
ggqqplot(df.StrI.central.shoot$total_antioxidants)

hist(residuals(model.bp.total_antioxidants.C.shoot))

# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.total_antioxidants.C.shoot )
MSerror<-deviance(model.bp.total_antioxidants.C.shoot)/df
comparison <- LSD.test(model.bp.total_antioxidants.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.total_antioxidants.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison

# log transformation
df.StrI.central.shoot$total_antioxidants_log10<-log10(df.StrI.central.shoot$total_antioxidants)
df.StrI.central.shoot

bp.total_antioxidants_log10.C.shoot<-ggplot(df.StrI.central.shoot,aes(TREATMENT, total_antioxidants_log10 ,  color=Bms)) + geom_boxplot()
bp.total_antioxidants_log10.C.shoot

hist(df.StrI.central.shoot$total_antioxidants_log10)
shapiro.test(df.StrI.central.shoot$total_antioxidants_log10)

#Kolmogorov-Smirnov Tests for normality
ks.test(df.StrI.central.shoot$total_antioxidants_log10, rnorm(n=length(df.StrI.central.shoot$total_antioxidants_log10), mean=mean(df.StrI.central.shoot$total_antioxidants_log10), sd=sd(df.StrI.central.shoot$total_antioxidants_log10)))

# create ANOVA model
model.bp.total_antioxidants_log10.C.shoot<-lm( total_antioxidants_log10  ~ TREATMENT, data = df.StrI.central.shoot)
anova(model.bp.total_antioxidants_log10.C.shoot)

#let's check if the data are now normal
# Density plot
ggdensity(df.StrI.central.shoot$total_antioxidants_log10, fill = "lightgray")
# QQ plot
ggqqplot(df.StrI.central.shoot$total_antioxidants_log10)

hist(residuals(model.bp.total_antioxidants_log10.C.shoot))


# create ANOVA post-hoc test model (with two tests, Fisher LSD and Tukey HSD)
df<-df.residual(model.bp.total_antioxidants_log10.C.shoot )
MSerror<-deviance(model.bp.total_antioxidants_log10.C.shoot)/df
comparison <- LSD.test(model.bp.total_antioxidants_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison
comparison <- HSD.test(model.bp.total_antioxidants_log10.C.shoot,"TREATMENT", df, MSerror, group=T)
comparison


##### ggarrange ####


#  BOXPLOTS arranged VERTICALLY

figure.hormones.vertical<-ggarrange(bp.APX.C.shoot, bp.GB.C.shoot, bp.H2O2.C.shoot, 
                                    bp.phenolics.C.shoot,bp.proline.C.shoot,bp.SOD.C.shoot,bp.total_antioxidants.C.shoot,labels = c("APX ", "GB", "H2O2", "Phenolics", "Proline", " SOD", "Total antioxidants"),
                                    ncol = 3, nrow = 3,   font.label = list(size = 10, color = "black", face = "bold"))
figure.hormones.vertical

figure.hormones.vertical.no.lettering<-ggarrange(bp.APX.C.shoot, bp.GB.C.shoot, bp.H2O2.C.shoot, 
                                                 bp.phenolics.C.shoot,bp.proline.C.shoot,bp.SOD.C.shoot,bp.total_antioxidants.C.shoot,
                                                 ncol = 3, nrow = 3)
figure.hormones.vertical.no.lettering


figure.hormones.vertical.lettering<-ggarrange(bp.APX.C.shoot, bp.GB.C.shoot, bp.H2O2.C.shoot, 
                                              bp.phenolics.C.shoot,bp.proline.C.shoot,bp.SOD.C.shoot,bp.total_antioxidants.C.shoot,
                                              labels = c("A", "B", "C","D", "E", "F"),
                                              ncol = 3, nrow = 3,   font.label = list(size = 16, color = "black", face = "bold"))
figure.hormones.vertical.lettering


#  BOXPLOTS arranged HORIZONTALY

figure.hormones.horizontal<-ggarrange(bp.APX.C.shoot, bp.GB.C.shoot, bp.H2O2.C.shoot, 
                                    bp.phenolics.C.shoot,bp.proline.C.shoot,bp.SOD.C.shoot,bp.total_antioxidants.C.shoot,labels = c("APX ", "GB", "H2O2", "Phenolics", "Proline", " SOD", "Total antioxidants"),
                                    ncol = 4, nrow = 2,   font.label = list(size = 10, color = "black", face = "bold"))
figure.hormones.horizontal

figure.hormones.horizontal.no.lettering<-ggarrange(bp.APX.C.shoot, bp.GB.C.shoot, bp.H2O2.C.shoot, 
                                                   bp.phenolics.C.shoot,bp.proline.C.shoot,bp.SOD.C.shoot,bp.total_antioxidants.C.shoot, 
                                                 ncol = 4, nrow = 2)
figure.hormones.horizontal.no.lettering


figure.hormones.horizontal.lettering<-ggarrange(bp.APX.C.shoot, bp.GB.C.shoot, bp.H2O2.C.shoot, 
                                                          bp.phenolics.C.shoot,bp.proline.C.shoot,bp.SOD.C.shoot,bp.total_antioxidants.C.shoot,
                                              labels = c("A", "B", "C","D", "E", "F", "G"),
                                              ncol = 4, nrow = 2,   font.label = list(size = 16, color = "black", face = "bold"))
figure.hormones.horizontal.lettering


# 
# 
# 
# figure.APX.SOD<-ggarrange(bp.APX.C.shoot, bp.SOD.C.shoot, 
#                          labels = c("Ascorbate peroxidase activity", "Superoxide dismutase"),
#                          ncol = 2, nrow = 1,   font.label = list(size = 10, color = "black", face = "bold"))
# 
# figure.APX.SOD
# 
# 
# 
# 
# figure.APX.SOD.H2O2.TOTAL_Antiox<-ggarrange(bp.APX.C.shoot, bp.SOD.C.shoot,bp.H2O2.C.shoot, bp.total_antioxidants.C.shoot,
#                             labels = c("Ascorbate peroxidase activity", "Superoxide dismutase", "      H202", "Total antioxidants"),
#                             ncol = 4, nrow = 1,   font.label = list(size = 10, color = "black", face = "bold"))
# 
# figure.APX.SOD.H2O2.TOTAL_Antiox
# 
