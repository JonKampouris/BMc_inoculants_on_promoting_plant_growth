# Package load ----------------------------------------------------------
pacman::p_load(tidyverse, # data import and handling
               conflicted, # handling function conflicts
               readODS, # import from ods files
               ComplexHeatmap, # create complex heatmap
               circlize, # color ramp
               fastcluster, # cluster algorithm
               cluster) # clustern

#### Resolve pakage conflicts ####
conflict_prefer("filter","dplyr")
conflict_prefer("selected", "dplyr")
conflict_prefer("select","dplyr")

BiocManager::install("ComplexHeatmap")

# Data import -------------------------------------------------------------

#### Set the WD using this  object ####
# change path if necessary 
setwd("./NAS/documents/work_jan/poster/microbe2022/Pricipitation/metabolites/")
# 

#### Check to make sure everything worked ####
getwd()


#### Import data ####
dat <- read_ods("LTE2020_Metabolites_RW_PCA_all.ods")
dat

remove.packages("cli")
install.packages("cli")

remove.packages("dplyr")
install.packages("dplyr")

library(tidyverse)

#### rename metabolites ###
dat <- dat %>% rename("6MeBOA" = "2,4-dihydroxy-7-methoxy-(2H)_1,4-benzoxazin-3 (4H)-one",
                      "4MeBOA" ="4-Methyl-2H-1,4-benzoxazin-3(4H)-one",
                      "6NH2BOA" ="6-Amino-2H-1,4-benzoxazin-3(4H)-one",
                      MBOA ="6-Methoxy-2-benzoxazolinone_(MBOA) ",
                      Asparagine = "Asparagin",
                      "2HBOA" = "Benzoic_acid",
                      "Caffeic acid" = "Caffeic_acid",
                      "Catechin hydrate" = "Catechin_Hydrate",
                      "cis-Aconitic acid" = "cis-aconitic_acid",
                      "Citric acid" ="citric_acid",
                      "Fumaric acid" = "fumaric_acid",
                      Glucose = Glucose,
                      Glutamine = Glutamine,
                      Glycine = Glycin,
                      "Malic acid" = "malic_acid",
                      "p-Coumaric acid" = "p-Coumaric_acid",
                      "Quercetin-Naringenin" = "Quercetin-Naringenin",
                      Serine = Serin,
                      "Succinic acid" = "succinic_acid",
                      "trans-Aconitic acid" = "trans-aconitic_acid",
                      "trans-Cinnamic acid" = "Trans-Cinnamic_acid",
                      Trehalose = Trehalose,
                      Tryptophan = Tryptophan)



# dat <- dat %>% rename(Asparagine = "Asparagin",
#                       "Caffeic acid" = "Caffeic_acid",
#                       "Catechin hydrate" = "Catechin_Hydrate",
#                       "cis-Aconitic acid" = "cis-aconitic_acid",
#                       "Citric acid" ="citric_acid",
#                       "Fumaric acid" = "fumaric_acid",
#                       Glucose = Glucose,
#                       Glutamine = Glutamine,
#                       Glycine = Glycin,
#                       "Malic acid" = "malic_acid",
#                       "p-Coumaric acid" = "p-Coumaric_acid",
#                       "Quercetin-Naringenin" = "Quercetin-Naringenin",
#                       Serine = Serin,
#                       "Succinic acid" = "succinic_acid",
#                       "trans-Aconitic acid" = "trans-aconitic_acid",
#                       "trans-Cinnamic acid" = "Trans-Cinnamic_acid",
#                       Trehalose = Trehalose,
#                       Tryptophan = Tryptophan)

# Data overview -----------------------------------------------------------

#### Data metrics ####
names(dat)
str(dat)

dat$TILLAGE <- factor(
  dat$TILLAGE,
  levels = c("CT", "MP"),
  labels = c("CT", "MP")
)

dat$FERTILIZATION <- factor(
  dat$FERTILIZATION,
  levels = c("Int", "Ext"),
  labels = c("Int", "Ext")
)

dat$BMS <- factor(
  dat$BMS,
  levels = c("control", "Bms"),
  labels = c("Control", "BMs")
)

dat$LOCATION <- factor(
  dat$LOCATION,
  levels = c("apical","basal","crown_roots"),
  labels = c("apical","basal","crown roots")
)




# Matrix ------------------------------------------------------------------
#### Put data in matrix for heatmap
## Scale Data 
# (Comment: I scaled because otherwise you only see differences in highly abundant metabolites.
# In the comment Guenther suggested to use real values, however than you have to adapt the color ramp.
# You can also try transforming and centering the data)
dat_mat <- as.matrix(scale(dat[,6:28], center = FALSE))
# dat_mat <- as.matrix(dat)

# transpose matrix
dat_mat <- t(dat_mat)

# select columns for annotation
row_Ano <- as.matrix(dat[,1:5])
#row_Ano <- t(row_Ano)


# Color -------------------------------------------------------------------
# Can be adapted
#f1 = col_fun = circlize::colorRamp2(c(-4, 0, 4), c("blue", "white", "red"), space = "RGB")
f2 = colorRamp2(seq(min(dat_mat), max(dat_mat), length = 3), c("blue", "#EEEEEE", "red"), space = "LAB")


# Row annotatiion ---------------------------------------------------------
col_BMs = columnAnnotation(
  Fertilization = row_Ano[,2], # Fertilization annotation
  Location = row_Ano[,4], # Annotation of Root location
  BMs=row_Ano[,3], # Annotation of BMs
  annotation_name_side="left", # Site to display annotation text
  show_annotation_name = TRUE, # Display annotation name
  gap = unit(1, 'mm'), # Space between annotation
  show_legend = c(TRUE), # Display Legend for annotation 
  col = list (Fertilization = c("Int" = "brown2", # colors for fertilization intensity
                                "Ext" = "forestgreen"),
              Location = c("crown roots" = "#7C3E66", # color for Location
                           "apical" = "#243A73",
                           "basal" = "#A5BECC"),
              BMs = c("Control" = "#146eb4", # color for BM
                      "BMs" = "#ff9900")))

# Draw Heatmap ------------------------------------------------------------
# alternative clustering algorithm, can be adjusted if needed
fh = function(x) fastcluster::hclust(dist(x))

# HEATMAP
ha <- Heatmap(dat_mat,
        top_annotation = col_BMs, # adds anotation to top (bottom also possible)
        col = f2, # adds color sheme for heatmap
        cluster_rows = function(m) as.dendrogram(diana(m)), #  cluster: divisive hierarchical clustering, DIANA = DIvisive ANAlysis
        clustering_distance_rows = function(x) as.dist((1-cor(t(x)))/2), # dendrogram distance
        column_order = sort(colnames(dat_mat)),
#        cluster_columns = fh, # cluster of columns
#        clustering_distance_columns = fh,
#        show_row_names = FALSE,
#        row_km = 5,
#        row_names_max_width = unit(15, "cm"),
        row_names_gp = gpar(fontsize= 18),
#        column_names_gp = gpar(fontsize= 18),
        show_heatmap_legend = TRUE
        )
      
ha  

draw(ha,heatmap_legend_side="top",annotation_legend_side ="top")
