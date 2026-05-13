library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(vegan)
library(ape)
library(broom)
library(patchwork)
library(effectsize)
library(rcompanion)
library(marginaleffects)
library(lme4)
####log linear regression
ILL30 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "ILL_MinF30", rowNames = T, skipEmptyRows = FALSE)
ILL35 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "ILL_MinF35", rowNames = T, skipEmptyRows = FALSE)
PB30 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "PB_MinF30", rowNames = T, skipEmptyRows = FALSE)
PB35 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "PB_MinF35", rowNames = T, skipEmptyRows = FALSE)

getSheetNames("Fungi OTU tables.xlsx")
ILL30 <- ILL30[-c(1:4),]
ILL35 <- ILL35[-c(1:4),]
PB30 <- PB30[-c(1:4),]
PB35 <- PB35[-c(1:4),]

PB30[] <- lapply(PB30, as.numeric)
PB35[] <-  lapply(PB35, as.numeric)
ILL30[] <-  lapply(ILL30, as.numeric)
ILL35[] <- lapply(ILL35, as.numeric)

PB30 <- PB30[rowSums(PB30)>0,]
PB35 <- PB35[rowSums(PB35)>0,]
ILL30 <- ILL30[rowSums(ILL30)>0,]
ILL35 <- ILL35[rowSums(ILL35)>0,]

#######
out.PB30  <- data.frame(seq = log(colSums(PB30 ), 2), richness = (colSums(PB30>0)))
out.PB35 <- data.frame(seq = log(colSums(PB35), 2), richness = (colSums(PB35>0)))
out.ILL30  <- data.frame(seq = log(colSums(ILL30), 2), richness = (colSums(ILL30>0)))
out.ILL35 <- data.frame(seq = log(colSums(ILL35), 2), richness = (colSums(ILL35>0)))

fit <- lm(richness ~ seq, data= out.PB30)
summary(fit)
performance::performance(fit)
out.PB30_re<-parameters::parameters(fit)
###
fit <- lm(richness ~ seq, data= out.PB35)
summary(fit)
performance::performance(fit)
out.PB35_re<-parameters::parameters(fit)
###

linear.function <-function(data1, data2, platform1, platform2){
  data1$Platform <- platform1
  data2$Platform <- platform2
  combined <- bind_rows(data1, data2)
  
  # Fit models separately and extract coefficients + adj.r.squared
  model_info <- combined %>%
    group_by(Platform) %>%
    do({
      m <- lm(richness ~ seq, data = .)
      tidy_m <- tidy(m)
      glance_m <- glance(m)
      data.frame(
        Intercept = tidy_m$estimate[tidy_m$term == "(Intercept)"],
        Slope     = tidy_m$estimate[tidy_m$term == "seq"],
        AdjR2     = glance_m$adj.r.squared
      )
    }) %>%
    ungroup()
  
  label_df <- combined %>%
    group_by(Platform) %>%
    dplyr::summarise(
      x = min(seq),         
      y = mean(richness)
    ) %>%
    left_join(model_info, by = "Platform") %>%
    mutate(
      label = sprintf("%s\nmodel: Radj² = %.3f\nrichness = %.2f + %.2f * log2(sequencing depth)",
                      Platform, AdjR2, Intercept, Slope)
    )
  
  my_colors <- setNames(
    c("#1B9E77", "#D95F02"),
    c(platform1, platform2)
  )
  
  pb_combined <- ggplot(combined, aes(x = seq, y = richness, color = Platform)) +
    geom_point(size = 3) +
    geom_abline(data = model_info, 
                aes(intercept = Intercept, slope = Slope, color = Platform),
                show.legend = FALSE, size = 0.9) +
    labs(y = "OTU richness",
         x = "log2(sequencing depth)") +
    geom_text(data = label_df,
              aes(x = x, y = y, label = label, color = Platform),
              inherit.aes = FALSE,
              hjust = 0, size = 3) +
    scale_color_manual(values = my_colors) +
    theme_classic()
  
  return(pb_combined)
}

pb_combined<-linear.function(data1=out.PB30, data2=out.PB35, platform1 = "PacBio 30", platform2 = "PacBio 35")

###
fit <- lm(richness ~ seq, data= out.ILL35)
summary(fit)
performance::performance(fit)
out.ILL35_re<-parameters::parameters(fit)
###
fit <- lm(richness ~ seq, data= out.ILL30)
summary(fit)
performance::performance(fit)
out.ILL30_re<-parameters::parameters(fit)
####
ILL_combined<-linear.function(out.ILL30, out.ILL35, "Illumina 30", "Illumina 35")

####artefact
PBC <- read.xlsx(xlsxFile = "Artefact accumul.xlsx", sheet = "PB_MinF_chim", rowNames = T, skipEmptyRows = FALSE)
PBL <- read.xlsx(xlsxFile = "Artefact accumul.xlsx", sheet = "PB_MinF_lowq", rowNames = T, skipEmptyRows = FALSE)
PBS <- read.xlsx(xlsxFile = "Artefact accumul.xlsx", sheet = "PB_MinF_switch", rowNames = T, skipEmptyRows = FALSE)
ILLC <- read.xlsx(xlsxFile = "Artefact accumul.xlsx", sheet = "ILL_MinF_chim", rowNames = T, skipEmptyRows = FALSE)
ILLL <- read.xlsx(xlsxFile = "Artefact accumul.xlsx", sheet = "ILL_MinF_lowq", rowNames = T, skipEmptyRows = FALSE)
ILLS <- read.xlsx(xlsxFile = "Artefact accumul.xlsx", sheet = "ILL_MinF_switch", rowNames = F, skipEmptyRows = FALSE)
a <- ILLS[duplicated(ILLS$OTU),]$OTU
ILLS <- ILLS[!ILLS$OTU %in% a, ]
row.names(ILLS)<- ILLS$OTU
ILLS <- ILLS[, -1]
getSheetNames("Artefact accumul.xlsx")

PBC30 <- PBC[,names(PBC) %in% names(PB30)]
PBC35 <- PBC[,names(PBC) %in% names(PB35)]
PBL30 <- PBL[,names(PBL) %in% names(PB30)]
PBL35 <- PBL[,names(PBL) %in% names(PB35)]
PBS30 <- PBS[,names(PBS) %in% names(PB30)]
PBS35 <- PBS[,names(PBS) %in% names(PB35)]

ILLC30 <- ILLC[,names(ILLC) %in% names(ILL30)]
ILLC35 <- ILLC[,names(ILLC) %in% names(ILL35)]
ILLL30 <- ILLL[,names(ILLL) %in% names(ILL30)]
ILLL35 <- ILLL[,names(ILLL) %in% names(ILL35)]
ILLS30 <- ILLS[,names(ILLS) %in% names(ILL30)]
ILLS35 <- ILLS[,names(ILLS) %in% names(ILL35)]

ILLL30 <- ILLL30[rowSums(ILLL30)>0,]
ILLL35 <- ILLL35[rowSums(ILLL35)>0,]
ILLC30 <- ILLC30[rowSums(ILLC30)>0,]
ILLC35 <- ILLC35[rowSums(ILLC35)>0,]
ILLS30 <- ILLS30[rowSums(ILLS30)>0,]
ILLS35 <- ILLS35[rowSums(ILLS35)>0,]

PBL30 <- PBL30[rowSums(PBL30)>0,]
PBL35 <- PBL35[rowSums(PBL35)>0,]
PBC30 <- PBC30[rowSums(PBC30)>0,]
PBC35 <- PBC35[rowSums(PBC35)>0,]
PBS30 <- PBS30[rowSums(PBS30)>0,]
PBS35 <- PBS35[rowSums(PBS35)>0,]
#######
out.PBL30  <- data.frame(seq = log(colSums(PBL30), 2), richness = (colSums(PBL30>0)))
out.PBL35 <- data.frame(seq = log(colSums(PBL35), 2), richness = (colSums(PBL35>0)))
out.PBC30  <- data.frame(seq = log(colSums(PBC30 ), 2), richness = (colSums(PBC30>0)))
out.PBC35 <- data.frame(seq = log(colSums(PBC35), 2), richness = (colSums(PBC35>0)))
out.PBS30  <- data.frame(seq = log(colSums(PBS30 ), 2), richness = (colSums(PBS30>0)))
out.PBS35 <- data.frame(seq = log(colSums(PBS35), 2), richness = (colSums(PBS35>0)))


out.ILLL30  <- data.frame(seq = log(colSums(ILLL30 ), 2), richness = (colSums(ILLL30>0)))
out.ILLL35 <- data.frame(seq = log(colSums(ILLL35), 2), richness = (colSums(ILLL35>0)))
out.ILLC30  <- data.frame(seq = log(colSums(ILLC30 ), 2), richness = (colSums(ILLC30>0)))
out.ILLC35 <- data.frame(seq = log(colSums(ILLC35), 2), richness = (colSums(ILLC35>0)))
out.ILLS30  <- data.frame(seq = log(colSums(ILLS30 ), 2), richness = (colSums(ILLS30>0)))
out.ILLS35 <- data.frame(seq = log(colSums(ILLS35), 2), richness = (colSums(ILLS35>0)))


fit <- lm(richness ~ seq, data= out.PBL30)
summary(fit)
performance::performance(fit)
out.PBL30_re <- parameters::parameters(fit)
###
fit <- lm(richness ~ seq, data= out.PBL35)
summary(fit)
performance::performance(fit)
out.PBL35_re <- parameters::parameters(fit)
####
PBL_combined <-linear.function(out.PBL30, out.PBL35, "PacBio lowS 30", "PacBio lowS 35")
#########################pbc
fit <- lm(richness ~ seq, data= out.PBC30)
summary(fit)
performance::performance(fit)
out.PBC30_re <- parameters::parameters(fit)
###
fit <- lm(richness ~ seq, data= out.PBC35)
summary(fit)
performance::performance(fit)
out.PBC35_re <- parameters::parameters(fit)
####
PBC_combined <-linear.function(out.PBC30, out.PBC35, "PacBio Chimera 30", "PacBio Chimera 35")

#########################pbs
fit <- lm(richness ~ seq, data= out.PBS30)
summary(fit)
performance::performance(fit)
out.PBS30_re <- parameters::parameters(fit)
###
out.PBS35 <- out.PBS35[is.finite(out.PBS35$seq), ]
fit <- lm(richness ~ seq, data= out.PBS35)
summary(fit)
performance::performance(fit)
out.PBS35_re <- parameters::parameters(fit)
############
PBS_combined  <-linear.function(out.PBS30, out.PBS35, "PacBio Switch 30", "PacBio Switch 35")
##illu
fit <- lm(richness ~ seq, data= out.ILLL30)
summary(fit)
performance::performance(fit)
out.ILLL30_re <- parameters::parameters(fit)
##
fit <- lm(richness ~ seq, data= out.ILLL35)
summary(fit)
performance::performance(fit)
out.ILLL35_re <- parameters::parameters(fit)
##
ILLL_combined  <-linear.function(out.ILLL30, out.ILLL35, "Illumina lowS 30", "Illumina lowS 35")

#########################
fit <- lm(richness ~ seq, data= out.ILLC30)
summary(fit)
performance::performance(fit)
out.ILLC30_re <- parameters::parameters(fit)

###
fit <- lm(richness ~ seq, data= out.ILLC35)
summary(fit)
performance::performance(fit)
out.ILLC35_re <- parameters::parameters(fit)
ILLC_combined  <-linear.function(out.ILLC30, out.ILLC35, "Illumina Chimera 30", "Illumina Chimera 35")

####ILLLS
fit <- lm(richness ~ seq, data= out.ILLS30)
summary(fit)
performance::performance(fit)
out.ILLS30_re <- parameters::parameters(fit)
###
fit <- lm(richness ~ seq, data= out.ILLS35)
summary(fit)
performance::performance(fit)
out.ILLS35_re <- parameters::parameters(fit)
ILLS_combined  <-linear.function(out.ILLS30, out.ILLS35, "Illumina Switch 30", "Illumina Switch 35")
####
(
  ILL_combined| pb_combined
) /
  (
    ILLL_combined| PBL_combined
  ) /
  (
    ILLC_combined| PBC_combined
  ) /
  (
    ILLS_combined| PBS_combined
  )
out.ILL30_re$combination <- "ILL30"
out.ILL35_re$combination<-  "ILL35"
out.ILLC30_re$combination<- "ILLChimera30"
out.ILLC35_re$combination<-  "ILLChimera35"
out.ILLL30_re$combination<- "ILLlowS30"
out.ILLL35_re$combination<-"ILLlowS35" 
out.ILLS30_re$combination<- "ILL Switch 30"
out.ILLS35_re$combination<- "ILL Switch 35"
out.PB30_re$combination<- "PacBio30"
out.PB35_re$combination<-   "PacBio35"
out.PBC30_re$combination<- "PacBio Chimera 30"
out.PBC35_re$combination<- "PacBio Chimera 35"
out.PBL30_re$combination<- "PacBio lowS 30"
out.PBL35_re$combination<- "PacBio lowS 35"
out.PBS30_re$combination<-  "PacBio switch 30"
out.PBS35_re$combination<- "PacBio switch 35"

all <- rbind(out.ILL30_re,out.ILL35_re,out.ILLC30_re,out.ILLC35_re, out.ILLL30_re,
             out.ILLL35_re,out.ILLS30_re,out.ILLS35_re,out.PB30_re,out.PB35_re,  
             out.PBC30_re,out.PBC35_re,out.PBL30_re,out.PBL35_re,out.PBS30_re, 
             out.PBS35_re)
write.csv(all, "log2.csv")
####betadiversity analysis
getSheetNames("Fungi OTU tables.xlsx")

PBMINF30 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "PB_MinF30", rowNames = T, skipEmptyRows = FALSE)
PBMINF35 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "PB_MinF35", rowNames = T, skipEmptyRows = FALSE)
PBSTD <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "PB_STDcyc", rowNames = T, skipEmptyRows = FALSE)
PBSTD30 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "PB_STD30", rowNames = T, skipEmptyRows = FALSE)
PBSTD35 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "PB_STD35", rowNames = T, skipEmptyRows = FALSE)
PBHYBR30 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "PB_Hybr30", rowNames = T, skipEmptyRows = FALSE)
PBHYBR35 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "PB_HYBR35", rowNames = T, skipEmptyRows = FALSE)
PBASV30 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "PacB_ASV30", rowNames = T, skipEmptyRows = FALSE)
PBASV35 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "PacB_ASV35", rowNames = T, skipEmptyRows = FALSE)

ILLMINF30 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "ILL_MinF30", rowNames = T, skipEmptyRows = FALSE)
ILLMINF35 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "ILL_MinF35", rowNames = T, skipEmptyRows = FALSE)
ILLSTD30 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "ILL_STD30", rowNames = T, skipEmptyRows = FALSE)
ILLSTD35 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "ILL_STD35", rowNames = T, skipEmptyRows = FALSE)
ILLHYBR30 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "ILL_hybr30", rowNames = T, skipEmptyRows = FALSE)
ILLHYBR35 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "ILL_hybr35", rowNames = T, skipEmptyRows = FALSE)
ILLASV30 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "ILL_ASV30", rowNames = T, skipEmptyRows = FALSE)
ILLASV35 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "ILL_ASV35", rowNames = T, skipEmptyRows = FALSE)

NCHYBR30 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "NC_hybr30", rowNames = T, skipEmptyRows = FALSE)
NCHYBR35 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "NC_hybr35", rowNames = T, skipEmptyRows = FALSE)
NCASV30 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "NC_ASV30", rowNames = T, skipEmptyRows = FALSE)
NCASV35 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "NC_ASV35", rowNames = T, skipEmptyRows = FALSE)

PRONHYBR30 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "Pron_hybr30", rowNames = T, skipEmptyRows = FALSE)
PRONHYBR35 <- read.xlsx(xlsxFile = "Fungi OTU tables.xlsx", sheet = "Pron_hybr35", rowNames = T, skipEmptyRows = FALSE)
PRONASV30 <- read.xlsx(xlsxFile = "Fungi OTU tables replacements.xlsx", sheet = "PRONAME_ASV30", rowNames = T, skipEmptyRows = FALSE)
PRONASV35 <- read.xlsx(xlsxFile = "Fungi OTU tables replacements.xlsx", sheet = "PRONAME_ASV35", rowNames = T, skipEmptyRows = FALSE)

#######
permanova <- function(data) {
  data <- as.data.frame(t(data))
  metadata <- data[,1:4]
  data <- data[,-c(1:4)]
  data[]<-lapply(data, as.numeric)
  data_hel <- decostand(data, method = "hellinger")
  distance<-vegdist(data_hel,method = "bray")
  permanova_result <- adonis2(
    distance ~ ecosys*design,
    data = metadata,
    permutations = 999,
    strata = metadata$plot
  )
  
  return(permanova_result)
}

ILLASV30_re <- permanova(ILLASV30)
ILLASV30_re$parameters <-row.names(ILLASV30_re)
ILLASV30_re$combinations <- "ILLASV30"

ILLASV35_re <- permanova(ILLASV35)
ILLASV35_re$parameters <-row.names(ILLASV35_re)
ILLASV35_re$combinations <- "ILLASV35"

ILLHYBR30_re <- permanova(ILLHYBR30)
ILLHYBR30_re$parameters <-row.names(ILLHYBR30_re)
ILLHYBR30_re$combinations <- "ILLHYBR30"

ILLHYBR35_re <- permanova(ILLHYBR35)
ILLHYBR35_re$parameters <-row.names(ILLHYBR35_re)
ILLHYBR35_re$combinations <- "ILLHYBR35"


ILLMINF30_re <- permanova(ILLMINF30)
ILLMINF30_re$parameters <-row.names(ILLMINF30_re)
ILLMINF30_re$combinations <- "ILLMINF30"


ILLMINF35_re <- permanova(ILLMINF35)
ILLMINF35_re$parameters <-row.names(ILLMINF35_re)
ILLMINF35_re$combinations <- "ILLMINF35"

ILLSTD30_re <- permanova(ILLSTD30)
ILLSTD30_re$parameters <-row.names(ILLSTD30_re)
ILLSTD30_re$combinations <- "ILLSTD30"

ILLSTD35_re <- permanova(ILLSTD35)
ILLSTD35_re$parameters <-row.names(ILLSTD35_re)
ILLSTD35_re$combinations <- "ILLSTD35"

NCASV30_re <- permanova(NCASV30)
NCASV30_re$parameters <-row.names(NCASV30_re)
NCASV30_re$combinations <- "NCASV30"

NCASV35_re <- permanova(NCASV35)
NCASV35_re$parameters <-row.names(NCASV35_re)
NCASV35_re$combinations <- "NCASV35"

NCHYBR30_re <- permanova(NCHYBR30)
NCHYBR30_re$parameters <-row.names(NCHYBR30_re)
NCHYBR30_re$combinations <- "NCHYBR30"

NCHYBR35_re <- permanova(NCHYBR35)
NCHYBR35_re$parameters <-row.names(NCHYBR35_re)
NCHYBR35_re$combinations <- "NCHYBR35"

PBASV30_re <- permanova(PBASV30)
PBASV30_re$parameters <-row.names(PBASV30_re)
PBASV30_re$combinations <- "PBASV30"

PBASV35_re <- permanova(PBASV35)
PBASV35_re$parameters <-row.names(PBASV35_re)
PBASV35_re$combinations <- "PBASV35"

PBHYBR30_re <- permanova(PBHYBR30)
PBHYBR30_re$parameters <-row.names(PBHYBR30_re)
PBHYBR30_re$combinations <- "PBHYBR30"

PBHYBR35_re <- permanova(PBHYBR35)
PBHYBR35_re$parameters <-row.names(PBHYBR35_re)
PBHYBR35_re$combinations <- "PBHYBR35"

PBMINF30_re <- permanova(PBMINF30)
PBMINF30_re$parameters <-row.names(PBMINF30_re)
PBMINF30_re$combinations <- "PBMINF30"

PBMINF35_re <- permanova(PBMINF35)
PBMINF35_re$parameters <-row.names(PBMINF35_re)
PBMINF35_re$combinations <- "PBMINF35"

PBSTD30_re <- permanova(PBSTD30)
PBSTD30_re$parameters <-row.names(PBSTD30_re)
PBSTD30_re$combinations <- "PBSTD30"

PBSTD35_re <- permanova(PBSTD35)
PBSTD35_re$parameters <-row.names(PBSTD35_re)
PBSTD35_re$combinations <- "PBSTD35"

PRONASV30_re <- permanova(PRONASV30)
PRONASV30_re$parameters <-row.names(PRONASV30_re)
PRONASV30_re$combinations <- "PRONASV30"

PRONASV35_re <- permanova(PRONASV35)
PRONASV35_re$parameters <-row.names(PRONASV35_re)
PRONASV35_re$combinations <- "PRONASV35"

PRONHYBR30_re <- permanova(PRONHYBR30)
PRONHYBR30_re$parameters <-row.names(PRONHYBR30_re)
PRONHYBR30_re$combinations <- "PRONHYBR30"

PRONHYBR35_re <- permanova(PRONHYBR35)
PRONHYBR35_re$parameters <-row.names(PRONHYBR35_re)
PRONHYBR35_re$combinations <- "PRONHYBR35"

all_result<-rbind(ILLASV30_re, ILLASV35_re, ILLHYBR30_re, ILLHYBR35_re,
                  ILLMINF30_re, ILLMINF35_re, ILLSTD30_re, ILLSTD35_re,
                  NCASV30_re, NCASV35_re, NCHYBR30_re, NCHYBR35_re,
                  PBASV30_re, PBASV35_re, PBHYBR30_re, PBHYBR35_re,
                  PBMINF30_re, PBMINF35_re, PBSTD30_re, PBSTD35_re,
                  PRONASV30_re, PRONASV35_re, PRONHYBR30_re,
                  PRONHYBR35_re)
####PBSTD
permanova <- function(data) {
  data <- as.data.frame(t(data))
  metadata <- data[,1:4]
  data <- data[,-c(1:4)]
  data[]<-lapply(data, as.numeric)
  data_hel <- decostand(data, method = "hellinger")
  distance<-vegdist(data_hel,method = "bray")
  permanova_result <- adonis2(
    distance ~ ecosys*design + cyc,
    data = metadata,
    permutations = 999,
    strata = metadata$plot
  )
  
  return(permanova_result)
}

PBSTD_re <- permanova(PBSTD)
PBSTD_re$parameters <-row.names(PBSTD_re)
PBSTD_re$combinations <- "PacB_STD_cyc"

write.csv(PBSTD_re, "PBSTD_re.csv")
write.csv(all_result, "all_result.csv")

#####################models##
platform<-read.csv("platform.information2.csv",header = TRUE,row.names = 1,sep = ",")
platform2 <- platform[!grepl("mix",platform$Plot),]
platform2 <- platform2[!grepl("NegC|PosC|MOCK",platform2$Plot),]
platform2 <- platform2[grepl("HYBR|ASV",platform2$Bioinf),]
##fungal OTU richness
platform2 <- platform2 %>%
  mutate(across(8:50, as.numeric))
mod <- lm(LOG_FungiIDSp ~ LOG_FungiSeq,data = platform2)
fungisp.r<-resid(mod)
mod <- lm(LOG_FungiGenera ~ LOG_FungiSeq,data = platform2)
fungige.r<-resid(mod)
mod <- lm(LOG_FungiOTU100.match ~ LOG_FungiSeq,data = platform2)
fungi100.r<-resid(mod)
mod <- lm(LOG_FungiOTU ~ LOG_FungiSeq,data = platform2)
fungiotu.r<-resid(mod)
mod <- lmer(
  fungiotu.r ~ 
    PCR_cycles*Design + Platform + Ecosystem+
    Design:Ecosystem +
    Platform:PCR_cycles+
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)

a<-anova(mod)
b<-performance::performance(mod)
c<-parameters::parameters(mod)
d<-eta_squared(mod, partial = TRUE)
d$group<-"OTU richnsss"
a$group<-"OTU richnsss"
b$group<-"OTU richnsss"
c$group<-"OTU richnsss"
##fungal species richness
mod <- lmer(
  fungisp.r ~ 
    PCR_cycles*Design + Platform + Ecosystem+
    Design:Ecosystem +
    Platform:PCR_cycles+
    Bioinf + (1 |Plot),
  data = platform2,
  REML = TRUE
)

a1<-anova(mod)
b1<-performance::performance(mod)
c1<-parameters::parameters(mod)
d1<-eta_squared(mod, partial = TRUE)
d1$group<-"sp richnsss"
a1$group<-"sp richnsss"
b1$group<-"sp richnsss"
c1$group<-"sp richnsss"

##fungal genera richness
mod <- lmer(
  fungige.r ~ 
    PCR_cycles*Design + Platform + Ecosystem+
    Design:Ecosystem +
    Platform:PCR_cycles+
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)

a2<-anova(mod)
b2<-performance::performance(mod)
c2<-parameters::parameters(mod)
d2<-eta_squared(mod, partial = TRUE)
d2$group<-"ge richnsss"
a2$group<-"ge richnsss"
b2$group<-"ge richnsss"
c2$group<-"ge richnsss"

##fungal 100% richness
mod <- lmer(
  fungi100.r ~ 
    PCR_cycles*Design + Platform + Ecosystem+
    Design:Ecosystem +
    Platform:PCR_cycles+
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)
a3<-anova(mod)
b3<-performance::performance(mod)
c3<-parameters::parameters(mod)
d3<-eta_squared(mod, partial = TRUE)
d3$group<-"100%richnsss"
a3$group<-"100%richnsss"
b3$group<-"100%richnsss"
c3$group<-"100%richnsss"

##fungal shannon
platform2$FungiShannon <- as.numeric(platform2$FungiShannon)
mod <- lmer(
  FungiShannon ~ 
    PCR_cycles*Design + Platform + Ecosystem+
    Design:Ecosystem +
    Platform:PCR_cycles+
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a4<-anova(mod)
b4<-performance::performance(mod)
c4<-parameters::parameters(mod)
d4<-eta_squared(mod, partial = TRUE)
d4$group<-"shannon"
a4$group<-"shannon"
b4$group<-"shannon"
c4$group<-"shannon"

##fungal evenness
platform2$FungiEvenness <- as.numeric(platform2$FungiEvenness)
mod <- lmer(
  FungiEvenness ~ 
    PCR_cycles*Design + Platform + Ecosystem+
    Design:Ecosystem +
    Platform:PCR_cycles+
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a5<-anova(mod)
b5<-performance::performance(mod)
c5<-parameters::parameters(mod)
d5<-eta_squared(mod, partial = TRUE)
d5$group<-"evenness"
a5$group<-"evenness"
b5$group<-"evenness"
c5$group<-"evenness"

####LR_ChimSeqs
mod <- lmer(
  LR_ChimSeqs ~ 
    PCR_cycles*Design + Platform + Ecosystem+
    Design:Ecosystem +
    Platform:PCR_cycles+
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a6<-anova(mod)
b6<-performance::performance(mod)
c6<-parameters::parameters(mod)
d6<-eta_squared(mod, partial = TRUE)
d6$group<-"LR_ChimSeqs"
a6$group<-"LR_ChimSeqs"
b6$group<-"LR_ChimSeqs"
c6$group<-"LR_ChimSeqs"

###LR_LowqSeqs
mod <- lmer(
  LR_LowqSeqs ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a7<-anova(mod)
b7<-performance::performance(mod)
c7<-parameters::parameters(mod)
d7<-eta_squared(mod, partial = TRUE)
d7$group<-"LR_LowqSeqs"
a7$group<-"LR_LowqSeqs"
b7$group<-"LR_LowqSeqs"
c7$group<-"LR_LowqSeqs"

###LR_OKSeqs
platform2$LR_OKSeqs <- as.numeric(platform2$LR_OKSeqs)
mod <- lmer(
  LR_OKSeqs ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2[!is.na(platform2$LR_OKSeqs),],
  REML = TRUE
)
a8<-anova(mod)
b8<-performance::performance(mod)
c8<-parameters::parameters(mod)
d8<-eta_squared(mod, partial = TRUE)
d8$group<-"LR_OKSeqs"
a8$group<-"LR_OKSeqs"
b8$group<-"LR_OKSeqs"
c8$group<-"LR_OKSeqs"

###TagswitchOTU% (bad model)
mod <- lmer(
  TagswitchOTU. ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a9<-anova(mod)
b9<-performance::performance(mod)
c9<-parameters::parameters(mod)
d9<-eta_squared(mod, partial = TRUE)
d9$group<-"TagswitchOTU%"
a9$group<-"TagswitchOTU%"
b9$group<-"TagswitchOTU%"
c9$group<-"TagswitchOTU%"

###TagswitchSeq% 
mod <- lmer(
  TagswitchSeq. ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a91<-anova(mod)
b91<-performance::performance(mod)
c91<-parameters::parameters(mod)
d91<-eta_squared(mod, partial = TRUE)
d91$group<-"TagswitchSeq%"
a91$group<-"TagswitchSeq%"
b91$group<-"TagswitchSeq%"
c91$group<-"TagswitchSeq%"

###PosCSwitchWeighted
mod <- lmer(
  PosCSwitchWeighted ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a92<-anova(mod)
b92<-performance::performance(mod)
c92<-parameters::parameters(mod)
d92<-eta_squared(mod, partial = TRUE)
d92$group<-"PosCSwitchWeighted"
a92$group<-"PosCSwitchWeighted"
b92$group<-"PosCSwitchWeighted"
c92$group<-"PosCSwitchWeighted"

###Residuals of LOG_EukOTU
mod <- lm(LOG_EukOTU ~ LOG_EukSeq,data = platform2)
EUKotu.r<-resid(mod)
mod <- lm(LOG_EukGenera ~ LOG_EukSeq,data = platform2)
EUKge.r<-resid(mod)

###
mod <- lmer(
  EUKotu.r ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)

a93<-anova(mod)
b93<-performance::performance(mod)
c93<-parameters::parameters(mod)
d93<-eta_squared(mod, partial = TRUE)
d93$group<-"EUKotu.r"
a93$group<-"EUKotu.r"
b93$group<-"EUKotu.r"
c93$group<-"EUKotu.r"

##EUK genus richness
mod <- lmer(
  EUKge.r ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)

a94<-anova(mod)
b94<-performance::performance(mod)
c94<-parameters::parameters(mod)
d94<-eta_squared(mod, partial = TRUE)
d94$group<-"EUKge.r"
a94$group<-"EUKge.r"
b94$group<-"EUKge.r"
c94$group<-"EUKge.r"
###
bll<-rbind(b,b1,b2,b3,b4,b5,b6,b7,b8,b9,b91,b92,b93,b94)
cll<-rbind(c,c1,c2,c3,c4,c5,c6,c7,c8,c9,c91,c92,c93,c94)
all<-rbind(a,a1,a2,a3,a4,a5,a6,a7,a8,a9,a91,a92,a93,a94)
dll<-rbind(d,d1,d2,d3,d4,d5,d6,d7,d8,d9,d91,d92,d93,d94)
write.csv(bll, "models.ASV.HYBR.Rsquare.csv")
write.csv(cll, "models.ASV.HYBR.parameters.csv")
write.csv(cbind(dll,all), "models.ASV.HYBR.significance.csv")
write.csv(dll, "models.ASV.HYBR.Rsquare.each.factor.csv")
#####full dataset
platform<-read.csv("platform.information2.csv",header = TRUE,row.names = 1,sep = ",")
platform2 <- platform[!grepl("mix",platform$Plot),]
platform2 <- platform2[!grepl("NegC|PosC|MOCK",platform2$Plot),]
##fungal OTU richness
platform2 <- platform2 %>%
  mutate(across(8:50, as.numeric))
mod <- lm(LOG_FungiIDSp ~ LOG_FungiSeq,data = platform2)
fungisp.r<-resid(mod)
mod <- lm(LOG_FungiGenera ~ LOG_FungiSeq,data = platform2)
fungige.r<-resid(mod)
mod <- lm(LOG_FungiOTU100.match ~ LOG_FungiSeq,data = platform2)
fungi100.r<-resid(mod)
mod <- lm(LOG_FungiOTU ~ LOG_FungiSeq,data = platform2)
fungiotu.r<-resid(mod)
mod <- lmer(
  fungiotu.r ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)

a<-anova(mod)
b<-performance::performance(mod)
c<-parameters::parameters(mod)
d<-eta_squared(mod, partial = TRUE)
d$group<-"OTU richnsss"
a$group<-"OTU richnsss"
b$group<-"OTU richnsss"
c$group<-"OTU richnsss"
##fungal species richness
mod <- lmer(
  fungisp.r ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)

a1<-anova(mod)
b1<-performance::performance(mod)
c1<-parameters::parameters(mod)
d1<-eta_squared(mod, partial = TRUE)
d1$group<-"sp richnsss"
a1$group<-"sp richnsss"
b1$group<-"sp richnsss"
c1$group<-"sp richnsss"

##fungal genera richness
mod <- lmer(
  fungige.r ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)

a2<-anova(mod)
b2<-performance::performance(mod)
c2<-parameters::parameters(mod)
d2<-eta_squared(mod, partial = TRUE)
d2$group<-"ge richnsss"
a2$group<-"ge richnsss"
b2$group<-"ge richnsss"
c2$group<-"ge richnsss"

##fungal 100% richness
mod <- lmer(
  fungi100.r ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)
a3<-anova(mod)
b3<-performance::performance(mod)
c3<-parameters::parameters(mod)
d3<-eta_squared(mod, partial = TRUE)
d3$group<-"100%richnsss"
a3$group<-"100%richnsss"
b3$group<-"100%richnsss"
c3$group<-"100%richnsss"

##fungal shannon
platform2$FungiShannon <- as.numeric(platform2$FungiShannon)
mod <- lmer(
  FungiShannon ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a4<-anova(mod)
b4<-performance::performance(mod)
c4<-parameters::parameters(mod)
d4<-eta_squared(mod, partial = TRUE)
d4$group<-"shannon"
a4$group<-"shannon"
b4$group<-"shannon"
c4$group<-"shannon"

##fungal evenness
platform2$FungiEvenness <- as.numeric(platform2$FungiEvenness)
mod <- lmer(
  FungiEvenness ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a5<-anova(mod)
b5<-performance::performance(mod)
c5<-parameters::parameters(mod)
d5<-eta_squared(mod, partial = TRUE)
d5$group<-"evenness"
a5$group<-"evenness"
b5$group<-"evenness"
c5$group<-"evenness"

####LR_ChimSeqs
mod <- lmer(
  LR_ChimSeqs ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a6<-anova(mod)
b6<-performance::performance(mod)
c6<-parameters::parameters(mod)
d6<-eta_squared(mod, partial = TRUE)
d6$group<-"LR_ChimSeqs"
a6$group<-"LR_ChimSeqs"
b6$group<-"LR_ChimSeqs"
c6$group<-"LR_ChimSeqs"

###LR_LowqSeqs
mod <- lmer(
  LR_LowqSeqs ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a7<-anova(mod)
b7<-performance::performance(mod)
c7<-parameters::parameters(mod)
d7<-eta_squared(mod, partial = TRUE)
d7$group<-"LR_LowqSeqs"
a7$group<-"LR_LowqSeqs"
b7$group<-"LR_LowqSeqs"
c7$group<-"LR_LowqSeqs"

###LR_OKSeqs
platform2$LR_OKSeqs <- as.numeric(platform2$LR_OKSeqs)
mod <- lmer(
  LR_OKSeqs ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2[!is.na(platform2$LR_OKSeqs),],
  REML = TRUE
)
a8<-anova(mod)
b8<-performance::performance(mod)
c8<-parameters::parameters(mod)
d8<-eta_squared(mod, partial = TRUE)
d8$group<-"LR_OKSeqs"
a8$group<-"LR_OKSeqs"
b8$group<-"LR_OKSeqs"
c8$group<-"LR_OKSeqs"

###TagswitchOTU% (bad model)
mod <- lmer(
  TagswitchOTU. ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a9<-anova(mod)
b9<-performance::performance(mod)
c9<-parameters::parameters(mod)
d9<-eta_squared(mod, partial = TRUE)
d9$group<-"TagswitchOTU%"
a9$group<-"TagswitchOTU%"
b9$group<-"TagswitchOTU%"
c9$group<-"TagswitchOTU%"

###TagswitchSeq% 
mod <- lmer(
  TagswitchSeq. ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a91<-anova(mod)
b91<-performance::performance(mod)
c91<-parameters::parameters(mod)
d91<-eta_squared(mod, partial = TRUE)
d91$group<-"TagswitchSeq%"
a91$group<-"TagswitchSeq%"
b91$group<-"TagswitchSeq%"
c91$group<-"TagswitchSeq%"

###PosCSwitchWeighted
mod <- lmer(
  PosCSwitchWeighted ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a92<-anova(mod)
b92<-performance::performance(mod)
c92<-parameters::parameters(mod)
d92<-eta_squared(mod, partial = TRUE)
d92$group<-"PosCSwitchWeighted"
a92$group<-"PosCSwitchWeighted"
b92$group<-"PosCSwitchWeighted"
c92$group<-"PosCSwitchWeighted"

###Residuals of LOG_EukOTU
platform2$LOG_EukSeq<-as.numeric(platform2$LOG_EukSeq)
mod <- lm(LOG_EukOTU ~ LOG_EukSeq,data = platform2)
EUKotu.r<-resid(mod)
mod <- lm(LOG_EukGenera ~ LOG_EukSeq,data = platform2)
EUKge.r<-resid(mod)

###
mod <- lmer(
  EUKotu.r ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)

a93<-anova(mod)
b93<-performance::performance(mod)
c93<-parameters::parameters(mod)
d93<-eta_squared(mod, partial = TRUE)
d93$group<-"EUKotu.r"
a93$group<-"EUKotu.r"
b93$group<-"EUKotu.r"
c93$group<-"EUKotu.r"

##EUK genus richness
mod <- lmer(
  EUKge.r ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)

a94<-anova(mod)
b94<-performance::performance(mod)
c94<-parameters::parameters(mod)
d94<-eta_squared(mod, partial = TRUE)
d94$group<-"EUKge.r"
a94$group<-"EUKge.r"
b94$group<-"EUKge.r"
c94$group<-"EUKge.r"
###
b8$ICC<-NA
bll<-rbind(b,b1,b2,b3,b4,b5,b6,b7,b8,b9,b91,b92,b93,b94)
cll<-rbind(c,c1,c2,c3,c4,c5,c6,c7,c8,c9,c91,c92,c93,c94)
all<-rbind(a,a1,a2,a3,a4,a5,a6,a7,a8,a9,a91,a92,a93,a94)
dll<-rbind(d,d1,d2,d3,d4,d5,d6,d7,d8,d9,d91,d92,d93,d94)
write.csv(bll, "models.full.Rsquare.csv")
write.csv(cll, "models.full.parameters.csv")
write.csv(cbind(dll,all), "models.full.significance.csv")
write.csv(dll, "models.full.Rsquare.each.factor.csv")
###PacBio and Illumina
platform<-read.csv("platform.information2.csv",header = TRUE,row.names = 1,sep = ",")
platform2 <- platform[!grepl("mix",platform$Plot),]
platform2 <- platform2[!grepl("NegC|PosC|MOCK",platform2$Plot),]
platform2 <- platform2[grepl("Illumina|PacBio",platform2$Platform),]
##fungal OTU richness
platform2 <- platform2 %>%
  mutate(across(8:50, as.numeric))
mod <- lm(LOG_FungiIDSp ~ LOG_FungiSeq,data = platform2)
fungisp.r<-resid(mod)
mod <- lm(LOG_FungiGenera ~ LOG_FungiSeq,data = platform2)
fungige.r<-resid(mod)
mod <- lm(LOG_FungiOTU100.match ~ LOG_FungiSeq,data = platform2)
fungi100.r<-resid(mod)
mod <- lm(LOG_FungiOTU ~ LOG_FungiSeq,data = platform2)
fungiotu.r<-resid(mod)
mod <- lmer(
  fungiotu.r ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)

a<-anova(mod)
b<-performance::performance(mod)
c<-parameters::parameters(mod)
d<-eta_squared(mod, partial = TRUE)
d$group<-"OTU richnsss"
a$group<-"OTU richnsss"
b$group<-"OTU richnsss"
c$group<-"OTU richnsss"
##fungal species richness
mod <- lmer(
  fungisp.r ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)

a1<-anova(mod)
b1<-performance::performance(mod)
c1<-parameters::parameters(mod)
d1<-eta_squared(mod, partial = TRUE)
d1$group<-"sp richnsss"
a1$group<-"sp richnsss"
b1$group<-"sp richnsss"
c1$group<-"sp richnsss"

##fungal genera richness
mod <- lmer(
  fungige.r ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)

a2<-anova(mod)
b2<-performance::performance(mod)
c2<-parameters::parameters(mod)
d2<-eta_squared(mod, partial = TRUE)
d2$group<-"ge richnsss"
a2$group<-"ge richnsss"
b2$group<-"ge richnsss"
c2$group<-"ge richnsss"

##fungal 100% richness
mod <- lmer(
  fungi100.r ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)
a3<-anova(mod)
b3<-performance::performance(mod)
c3<-parameters::parameters(mod)
d3<-eta_squared(mod, partial = TRUE)
d3$group<-"100%richnsss"
a3$group<-"100%richnsss"
b3$group<-"100%richnsss"
c3$group<-"100%richnsss"

##fungal shannon
platform2$FungiShannon <- as.numeric(platform2$FungiShannon)
mod <- lmer(
  FungiShannon ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a4<-anova(mod)
b4<-performance::performance(mod)
c4<-parameters::parameters(mod)
d4<-eta_squared(mod, partial = TRUE)
d4$group<-"shannon"
a4$group<-"shannon"
b4$group<-"shannon"
c4$group<-"shannon"

##fungal evenness
platform2$FungiEvenness <- as.numeric(platform2$FungiEvenness)
mod <- lmer(
  FungiEvenness ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a5<-anova(mod)
b5<-performance::performance(mod)
c5<-parameters::parameters(mod)
d5<-eta_squared(mod, partial = TRUE)
d5$group<-"evenness"
a5$group<-"evenness"
b5$group<-"evenness"
c5$group<-"evenness"

####LR_ChimSeqs
mod <- lmer(
  LR_ChimSeqs ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a6<-anova(mod)
b6<-performance::performance(mod)
c6<-parameters::parameters(mod)
d6<-eta_squared(mod, partial = TRUE)
d6$group<-"LR_ChimSeqs"
a6$group<-"LR_ChimSeqs"
b6$group<-"LR_ChimSeqs"
c6$group<-"LR_ChimSeqs"

###LR_LowqSeqs
mod <- lmer(
  LR_LowqSeqs ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a7<-anova(mod)
b7<-performance::performance(mod)
c7<-parameters::parameters(mod)
d7<-eta_squared(mod, partial = TRUE)
d7$group<-"LR_LowqSeqs"
a7$group<-"LR_LowqSeqs"
b7$group<-"LR_LowqSeqs"
c7$group<-"LR_LowqSeqs"

###LR_OKSeqs
platform2$LR_OKSeqs <- as.numeric(platform2$LR_OKSeqs)
mod <- lmer(
  LR_OKSeqs ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2[!is.na(platform2$LR_OKSeqs),],
  REML = TRUE
)
a8<-anova(mod)
b8<-performance::performance(mod)
c8<-parameters::parameters(mod)
d8<-eta_squared(mod, partial = TRUE)
d8$group<-"LR_OKSeqs"
a8$group<-"LR_OKSeqs"
b8$group<-"LR_OKSeqs"
c8$group<-"LR_OKSeqs"

###TagswitchOTU% (bad model)
mod <- lmer(
  TagswitchOTU. ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a9<-anova(mod)
b9<-performance::performance(mod)
c9<-parameters::parameters(mod)
d9<-eta_squared(mod, partial = TRUE)
d9$group<-"TagswitchOTU%"
a9$group<-"TagswitchOTU%"
b9$group<-"TagswitchOTU%"
c9$group<-"TagswitchOTU%"

###TagswitchSeq% 
mod <- lmer(
  TagswitchSeq. ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a91<-anova(mod)
b91<-performance::performance(mod)
c91<-parameters::parameters(mod)
d91<-eta_squared(mod, partial = TRUE)
d91$group<-"TagswitchSeq%"
a91$group<-"TagswitchSeq%"
b91$group<-"TagswitchSeq%"
c91$group<-"TagswitchSeq%"

###PosCSwitchWeighted
mod <- lmer(
  PosCSwitchWeighted ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + LOG_FungiSeq+ (1 | Plot),
  data = platform2,
  REML = TRUE
)
a92<-anova(mod)
b92<-performance::performance(mod)
c92<-parameters::parameters(mod)
d92<-eta_squared(mod, partial = TRUE)
d92$group<-"PosCSwitchWeighted"
a92$group<-"PosCSwitchWeighted"
b92$group<-"PosCSwitchWeighted"
c92$group<-"PosCSwitchWeighted"

###Residuals of LOG_EukOTU
platform2$LOG_EukSeq<-as.numeric(platform2$LOG_EukSeq)
mod <- lm(LOG_EukOTU ~ LOG_EukSeq,data = platform2)
EUKotu.r<-resid(mod)
mod <- lm(LOG_EukGenera ~ LOG_EukSeq,data = platform2)
EUKge.r<-resid(mod)

###
mod <- lmer(
  EUKotu.r ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)

a93<-anova(mod)
b93<-performance::performance(mod)
c93<-parameters::parameters(mod)
d93<-eta_squared(mod, partial = TRUE)
d93$group<-"EUKotu.r"
a93$group<-"EUKotu.r"
b93$group<-"EUKotu.r"
c93$group<-"EUKotu.r"

##EUK genus richness
mod <- lmer(
  EUKge.r ~ 
    PCR_cycles*Design + Platform + 
    Platform:PCR_cycles+Ecosystem+
    Design:Ecosystem +
    Bioinf + (1 | Plot),
  data = platform2,
  REML = TRUE
)

a94<-anova(mod)
b94<-performance::performance(mod)
c94<-parameters::parameters(mod)
d94<-eta_squared(mod, partial = TRUE)
d94$group<-"EUKge.r"
a94$group<-"EUKge.r"
b94$group<-"EUKge.r"
c94$group<-"EUKge.r"
###
b8$ICC<-NA
bll<-rbind(b,b1,b2,b3,b4,b5,b6,b7,b8,b9,b91,b92,b93,b94)
cll<-rbind(c,c1,c2,c3,c4,c5,c6,c7,c8,c9,c91,c92,c93,c94)
all<-rbind(a,a1,a2,a3,a4,a5,a6,a7,a8,a9,a91,a92,a93,a94)
dll<-rbind(d,d1,d2,d3,d4,d5,d6,d7,d8,d9,d91,d92,d93,d94)
write.csv(dll, "models.Illumina.PacBio.Rsquare.each.factor.csv")
write.csv(bll, "models.Illumina.PacBio.Rsquare.csv")
write.csv(cll, "models.Illumina.PacBio.parameters.csv")
write.csv(cbind(dll,all), "models.Illumina.PacBio.significance.csv")
################ecology##
platform<-read.csv("platform.information2.csv",header = TRUE,row.names = 1,sep = ",")
platform2 <- platform[!grepl("mix",platform$Plot),]
platform2 <- platform2[!grepl("NegC|PosC|MOCK",platform2$Plot),]
platform2 <- platform2[grepl("HYBR|ASV",platform2$Bioinf),]
##fungal OTU richness
platform2 <- platform2 %>%
  mutate(across(8:50, as.numeric))
mod <- lm(LOG_FungiIDSp ~ LOG_FungiSeq,data = platform2)
fungisp.r<-resid(mod)
mod <- lm(LOG_FungiGenera ~ LOG_FungiSeq,data = platform2)
fungige.r<-resid(mod)
mod <- lm(LOG_FungiOTU100.match ~ LOG_FungiSeq,data = platform2)
fungi100.r<-resid(mod)
mod <- lm(LOG_FungiOTU ~ LOG_FungiSeq,data = platform2)
fungiotu.r<-resid(mod)
mod <- lmer(
  fungiotu.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem
     +(1| Plot),
  data = platform2
)

a<-anova(mod)
b<-performance::performance(mod)
c<-parameters::parameters(mod)
d<-eta_squared(mod, partial = TRUE)
d$group<-"OTU richnsss"
a$group<-"OTU richnsss"
b$group<-"OTU richnsss"
c$group<-"OTU richnsss"

##fungal species richness
mod <- lmer(
  fungisp.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem + 
    Bioinf + (1|Plot),
  data = platform2,
  REML = TRUE
)
a1<-anova(mod)
b1<-performance::performance(mod)
c1<-parameters::parameters(mod)
d1<-eta_squared(mod, partial = TRUE)
d1$group<-"sp richnsss"
a1$group<-"sp richnsss"
b1$group<-"sp richnsss"
c1$group<-"sp richnsss"

##fungal genera richness
mod <- lmer(
  fungige.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem + (1|Plot),
  data = platform2
)

a2<-anova(mod)
b2<-performance::performance(mod)
c2<-parameters::parameters(mod)
d2<-eta_squared(mod, partial = TRUE)
d2$group<-"ge richnsss"
a2$group<-"ge richnsss"
b2$group<-"ge richnsss"
c2$group<-"ge richnsss"

##fungal 100% richness
mod <- lmer(
  fungi100.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem
    + (1|Plot),
  data = platform2,
  REML = TRUE
)
a3<-anova(mod)
b3<-performance::performance(mod)
c3<-parameters::parameters(mod)
d3<-eta_squared(mod, partial = TRUE)
d3$group<-"100%richnsss"
a3$group<-"100%richnsss"
b3$group<-"100%richnsss"
c3$group<-"100%richnsss"

##fungal shannon
mod <- lmer(
  FungiShannon ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    LOG_FungiSeq + (1|Plot),
  data = platform2
)
a4<-anova(mod)
b4<-performance::performance(mod)
c4<-parameters::parameters(mod)
d4<-eta_squared(mod, partial = TRUE)
d4$group<-"shannon"
a4$group<-"shannon"
b4$group<-"shannon"
c4$group<-"shannon"

##fungal evenness
platform2$FungiEvenness <- as.numeric(platform2$FungiEvenness)
mod <- lmer(
  FungiEvenness ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem +LOG_FungiSeq + (1|Plot),
  data = platform2
)
a5<-anova(mod)
b5<-performance::performance(mod)
c5<-parameters::parameters(mod)
d5<-eta_squared(mod, partial = TRUE)
d5$group<-"evenness"
a5$group<-"evenness"
b5$group<-"evenness"
c5$group<-"evenness"

###Residuals of LOG_EukOTU
mod <- lm(LOG_EukOTU ~ LOG_EukSeq,data = platform2)
EUKotu.r<-resid(mod)
mod <- lm(LOG_EukGenera ~ LOG_EukSeq,data = platform2)
EUKge.r<-resid(mod)
###
mod <- lmer(
  EUKotu.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem
    + (1|Plot),
  data = platform2
)

a93<-anova(mod)
b93<-performance::performance(mod)
c93<-parameters::parameters(mod)
d93<-eta_squared(mod, partial = TRUE)
d93$group<-"EUKotu.r"
a93$group<-"EUKotu.r"
b93$group<-"EUKotu.r"
c93$group<-"EUKotu.r"

##fungal species richness
mod <- lmer(
  EUKge.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem + 
    + (1|Plot),
  data = platform2
)
a94<-anova(mod)
b94<-performance::performance(mod)
c94<-parameters::parameters(mod)
d94<-eta_squared(mod, partial = TRUE)
d94$group<-"EUKge.r"
a94$group<-"EUKge.r"
b94$group<-"EUKge.r"
c94$group<-"EUKge.r"

bll<-rbind(b,b1,b2,b3,b4,b5,b93,b94)
cll<-rbind(c,c1,c2,c3,c4,c5,c93,c94)
all<-rbind(a,a1,a2,a3,a4,a5,a93,a94)
dll<-rbind(d,d1,d2,d3,d4,d5,d93,d94)
write.csv(dll, "ecology.ASV.HYBR.Rsquare.for.each.factor.csv")
write.csv(bll, "ecology.ASV.HYBR.Rsquare.csv")
write.csv(cll, "ecology.ASV.HYBR.parameters.csv")
write.csv(cbind(dll,all), "ecology.ASV.HYBR.significance.csv")
###full dataset
##
platform<-read.csv("platform.information2.csv",header = TRUE,row.names = 1,sep = ",")
platform2 <- platform[!grepl("mix",platform$Plot),]
platform2 <- platform2[!grepl("NegC|PosC|MOCK",platform2$Plot),]
##fungal OTU richness
platform2 <- platform2 %>%
  mutate(across(8:50, as.numeric))
mod <- lm(LOG_FungiIDSp ~ LOG_FungiSeq,data = platform2)
fungisp.r<-resid(mod)
mod <- lm(LOG_FungiGenera ~ LOG_FungiSeq,data = platform2)
fungige.r<-resid(mod)
mod <- lm(LOG_FungiOTU100.match ~ LOG_FungiSeq,data = platform2)
fungi100.r<-resid(mod)
mod <- lm(LOG_FungiOTU ~ LOG_FungiSeq,data = platform2)
fungiotu.r<-resid(mod)
mod <- lmer(
  fungiotu.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem + 
    + (1|Plot),
  data = platform2
)

a<-anova(mod)
b<-performance::performance(mod)
c<-parameters::parameters(mod)
d<-eta_squared(mod, partial = TRUE)
d$group<-"OTU richnsss"
a$group<-"OTU richnsss"
b$group<-"OTU richnsss"
c$group<-"OTU richnsss"

##fungal species richness
mod <- lmer(
  fungisp.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem + 
    (1|Plot),
  data = platform2,
  REML = TRUE
)
a1<-anova(mod)
b1<-performance::performance(mod)
c1<-parameters::parameters(mod)
d1<-eta_squared(mod, partial = TRUE)
d1$group<-"sp richnsss"
a1$group<-"sp richnsss"
b1$group<-"sp richnsss"
c1$group<-"sp richnsss"

##fungal genera richness
mod <- lmer(
  fungige.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem + (1|Plot),
  data = platform2
)

a2<-anova(mod)
b2<-performance::performance(mod)
c2<-parameters::parameters(mod)
d2<-eta_squared(mod, partial = TRUE)
d2$group<-"ge richnsss"
a2$group<-"ge richnsss"
b2$group<-"ge richnsss"
c2$group<-"ge richnsss"

##fungal 100% richness
mod <- lmer(
  fungi100.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem + 
    (1|Plot),
  data = platform2,
  REML = TRUE
)
a3<-anova(mod)
b3<-performance::performance(mod)
c3<-parameters::parameters(mod)
d3<-eta_squared(mod, partial = TRUE)
d3$group<-"100%richnsss"
a3$group<-"100%richnsss"
b3$group<-"100%richnsss"
c3$group<-"100%richnsss"
##fungal shannon
mod <- lmer(
  FungiShannon ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem + LOG_FungiSeq + (1|Plot),
  data = platform2
)
a4<-anova(mod)
b4<-performance::performance(mod)
c4<-parameters::parameters(mod)
d4<-eta_squared(mod, partial = TRUE)
d4$group<-"shannon"
a4$group<-"shannon"
b4$group<-"shannon"
c4$group<-"shannon"

##fungal evenness
mod <- lmer(
  FungiEvenness ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem + LOG_FungiSeq + (1|Plot),
  data = platform2
)
a5<-anova(mod)
b5<-performance::performance(mod)
c5<-parameters::parameters(mod)
d5<-eta_squared(mod, partial = TRUE)
d5$group<-"evenness"
a5$group<-"evenness"
b5$group<-"evenness"
c5$group<-"evenness"

###Residuals of LOG_EukOTU
mod <- lm(LOG_EukOTU ~ LOG_EukSeq,data = platform2)
EUKotu.r<-resid(mod)
mod <- lm(LOG_EukGenera ~ LOG_EukSeq,data = platform2)
EUKge.r<-resid(mod)
###
mod <- lmer(
  EUKotu.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem + (1|Plot),
  data = platform2
)

a93<-anova(mod)
b93<-performance::performance(mod)
c93<-parameters::parameters(mod)
d93<-eta_squared(mod, partial = TRUE)
d93$group<-"EUKotu.r"
a93$group<-"EUKotu.r"
b93$group<-"EUKotu.r"
c93$group<-"EUKotu.r"

##fungal species richness
mod <- lmer(
  EUKge.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem 
    + (1|Plot),
  data = platform2
)
a94<-anova(mod)
b94<-performance::performance(mod)
c94<-parameters::parameters(mod)
d94<-eta_squared(mod, partial = TRUE)
d94$group<-"EUKge.r"
a94$group<-"EUKge.r"
b94$group<-"EUKge.r"
c94$group<-"EUKge.r"

bll<-rbind(b,b1,b2,b3,b4,b5,b93,b94)
cll<-rbind(c,c1,c2,c3,c4,c5,c93,c94)
all<-rbind(a,a1,a2,a3,a4,a5,a93,a94)
dll<-rbind(d,d1,d2,d3,d4,d5,d93,d94)
write.csv(dll, "ecology.full.Rsquare.each.factor.csv")
write.csv(bll, "ecology.full.Rsquare.csv")
write.csv(cll, "ecology.full.parameters.csv")
write.csv(cbind(dll,all), "ecology.full.significance.csv")

###PacBio and Illumina
platform<-read.csv("platform.information2.csv",header = TRUE,row.names = 1,sep = ",")
platform2 <- platform[!grepl("mix",platform$Plot),]
platform2 <- platform2[!grepl("NegC|PosC|MOCK",platform2$Plot),]
platform2 <- platform2[grepl("Illumina|PacBio",platform2$Platform),]
platform2 <- platform2 %>%
  mutate(across(8:50, as.numeric))
mod <- lm(LOG_FungiIDSp ~ LOG_FungiSeq,data = platform2)
fungisp.r<-resid(mod)
mod <- lm(LOG_FungiGenera ~ LOG_FungiSeq,data = platform2)
fungige.r<-resid(mod)
mod <- lm(LOG_FungiOTU100.match ~ LOG_FungiSeq,data = platform2)
fungi100.r<-resid(mod)
mod <- lm(LOG_FungiOTU ~ LOG_FungiSeq,data = platform2)
fungiotu.r<-resid(mod)
mod <- lmer(
  fungiotu.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem + (1|Plot),
  data = platform2
)

a<-anova(mod)
b<-performance::performance(mod)
c<-parameters::parameters(mod)
d<-eta_squared(mod, partial = TRUE)
d$group<-"OTU richnsss"
a$group<-"OTU richnsss"
b$group<-"OTU richnsss"
c$group<-"OTU richnsss"

##fungal species richness
mod <- lmer(
  fungisp.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem +
    (1|Plot),
  data = platform2,
  REML = TRUE
)
a1<-anova(mod)
b1<-performance::performance(mod)
c1<-parameters::parameters(mod)
d1<-eta_squared(mod, partial = TRUE)
d1$group<-"sp richnsss"
a1$group<-"sp richnsss"
b1$group<-"sp richnsss"
c1$group<-"sp richnsss"

##fungal genera richness
mod <- lmer(
  fungige.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem + (1|Plot),
  data = platform2
)

a2<-anova(mod)
b2<-performance::performance(mod)
c2<-parameters::parameters(mod)
d2<-eta_squared(mod, partial = TRUE)
d2$group<-"ge richnsss"
a2$group<-"ge richnsss"
b2$group<-"ge richnsss"
c2$group<-"ge richnsss"

##fungal 100% richness
mod <- lmer(
  fungi100.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem +
    (1|Plot),
  data = platform2,
  REML = TRUE
)
a3<-anova(mod)
b3<-performance::performance(mod)
c3<-parameters::parameters(mod)
d3<-eta_squared(mod, partial = TRUE)
d3$group<-"100%richnsss"
a3$group<-"100%richnsss"
b3$group<-"100%richnsss"
c3$group<-"100%richnsss"

##fungal shannon
mod <- lmer(
  FungiShannon ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem +LOG_FungiSeq+ (1|Plot),
  data = platform2
)
a4<-anova(mod)
b4<-performance::performance(mod)
c4<-parameters::parameters(mod)
d4<-eta_squared(mod, partial = TRUE)
d4$group<-"shannon"
a4$group<-"shannon"
b4$group<-"shannon"
c4$group<-"shannon"
##fungal evenness
mod <- lmer(
  FungiEvenness ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem +LOG_FungiSeq+ (1|Plot),
  data = platform2
)
a5<-anova(mod)
b5<-performance::performance(mod)
c5<-parameters::parameters(mod)
d5<-eta_squared(mod, partial = TRUE)
d5$group<-"evenness"
a5$group<-"evenness"
b5$group<-"evenness"
c5$group<-"evenness"

###Residuals of LOG_EukOTU
mod <- lm(LOG_EukOTU ~ LOG_EukSeq,data = platform2)
EUKotu.r<-resid(mod)
mod <- lm(LOG_EukGenera ~ LOG_EukSeq,data = platform2)
EUKge.r<-resid(mod)
###
mod <- lmer(
  EUKotu.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem +(1|Plot),
  data = platform2
)

a93<-anova(mod)
b93<-performance::performance(mod)
c93<-parameters::parameters(mod)
d93<-eta_squared(mod, partial = TRUE)
d93$group<-"EUKotu.r"
a93$group<-"EUKotu.r"
b93$group<-"EUKotu.r"
c93$group<-"EUKotu.r"

##fungal species richness
mod <- lmer(
  EUKge.r ~ 
    PCR_cycles+Design + Platform + Bioinf+Ecosystem+
    Design:Ecosystem + (1|Plot),
  data = platform2
)
a94<-anova(mod)
b94<-performance::performance(mod)
c94<-parameters::parameters(mod)
d94<-eta_squared(mod, partial = TRUE)
d94$group<-"EUKge.r"
a94$group<-"EUKge.r"
b94$group<-"EUKge.r"
c94$group<-"EUKge.r"

bll<-rbind(b,b1,b2,b3,b4,b5,b93,b94)
cll<-rbind(c,c1,c2,c3,c4,c5,c93,c94)
all<-rbind(a,a1,a2,a3,a4,a5,a93,a94)
dll<-rbind(d,d1,d2,d3,d4,d5,d93,d94)
write.csv(dll, "ecology.Illumina.PacBio.Rsquare.each.factor.csv")
write.csv(bll, "ecology.Illumina.PacBio.Rsquare.csv")
write.csv(cll, "ecology.Illumina.PacBio.parameters.csv")
write.csv(cbind(dll,all), "ecology.Illumina.PacBio.significance.csv")
#########################ecology2##
platform<-read.csv("platform.information2.csv",header = TRUE,row.names = 1,sep = ",")
platform2 <- platform[!grepl("mix",platform$Soil_sample),]
platform2 <- platform2[!grepl("NegC|PosC|MOCK",platform2$Plot),]
##get the subsets
platform2 <- platform2 %>%
  mutate(across(8:50, as.numeric))
mod <- lm(LOG_FungiIDSp ~ LOG_FungiSeq,data = platform2)
platform2$fungisp.r<-resid(mod)
mod <- lm(LOG_FungiGenera ~ LOG_FungiSeq,data = platform2)
platform2$fungige.r<-resid(mod)
mod <- lm(LOG_FungiOTU100.match ~ LOG_FungiSeq,data = platform2)
platform2$fungi100.r<-resid(mod)
mod <- lm(LOG_FungiOTU ~ LOG_FungiSeq,data = platform2)
platform2$fungiotu.r<-resid(mod)
platform2$LOG_EukSeq<-as.numeric(platform2$LOG_EukSeq)
mod <- lm(LOG_EukOTU ~ LOG_EukSeq,data = platform2)
platform2$EUKotu.r<-resid(mod)
mod <- lm(LOG_EukGenera ~ LOG_EukSeq,data = platform2)
platform2$EUKge.r<-resid(mod)
platform2$group <- paste0(platform2$Platform, "_", platform2$Bioinf,"_",platform2$PCR_cycles)
data <-list()
for (i in unique(platform2$group)) {
  data[[i]] <- platform2[platform2$group %in% i,]
}
##fungal OTU richness
result <- data.frame()
result.r <- data.frame()
result.r2 <- data.frame()
result.r3 <- data.frame()
for (i in names(data)) {
  mod <- lmer(
    fungiotu.r ~ 
      Design*Ecosystem + (1|Plot),
    data = data[[i]]
  )
  a<-performance::performance(mod)
  a$group <- i
  b<-parameters::parameters(mod) 
  b$group <- i
  d<-eta_squared(mod, partial = TRUE)
  c<-anova(mod)
  result <- rbind(result,b)
  result.r <- rbind(result.r,a[,c(1:5,ncol(a))])
  result.r2 <- rbind(result.r2,d)
  result.r3 <- rbind(result.r3,c)
}
write.csv(result.r,"ecology2.rsquare.fungi.richness.csv")
write.csv(result.r3,"ecology2.significance.fungi.richness.csv")
write.csv(result.r2,"ecology2.rsquare.fungi.richness.each.factor.csv")
write.csv(result,"ecology2.parameters.fungi.richness.csv")
##fungal species richness
result <- data.frame()
result.r <- data.frame()
result.r2 <- data.frame()
result.r3 <- data.frame()
for (i in names(data)) {
  mod <- lmer(
    fungisp.r ~ 
      Design*Ecosystem + (1|Plot),
    data = data[[i]]
  )
  a<-performance::performance(mod)
  a$group <- i
  b<-parameters::parameters(mod) 
  b$group <- i
  d<-eta_squared(mod, partial = TRUE)
  c<-anova(mod)
  result <- rbind(result,b)
  result.r <- rbind(result.r,a[,c(1:5,ncol(a))])
  result.r2 <- rbind(result.r2,d)
  result.r3 <- rbind(result.r3,c)
}
write.csv(result.r,"ecology2.rsquare.fungi.sp.richness.csv")
write.csv(result.r3,"ecology2.significance.fungi.sp.richness.csv")
write.csv(result.r2,"ecology2.rsquare.fungi.sp.richness.each.factor.csv")
write.csv(result,"ecology2.parameters.fungi.sp.richness.csv")
##fungal genera richness
result <- data.frame()
result.r <- data.frame()
result.r2 <- data.frame()
result.r3 <- data.frame()
for (i in names(data)) {
  mod <- lmer(
    fungige.r ~ 
      Design*Ecosystem+ (1|Plot),
    data = data[[i]]
  )
  a<-performance::performance(mod)
  a$group <- i
  b<-parameters::parameters(mod) 
  b$group <- i
  d<-eta_squared(mod, partial = TRUE)
  c<-anova(mod)
  result <- rbind(result,b)
  result.r <- rbind(result.r,a[,c(1:5,ncol(a))])
  result.r2 <- rbind(result.r2,d)
  result.r3 <- rbind(result.r3,c)
}
write.csv(result.r,"ecology2.rsquare.fungi.ge.richness.csv")
write.csv(result.r3,"ecology2.significance.fungi.ge.richness.csv")
write.csv(result.r2,"ecology2.rsquare.fungi.ge.richness.each.factor.csv")
write.csv(result,"ecology2.parameters.fungi.ge.richness.csv")

##fungal 100% richness
result <- data.frame()
result.r <- data.frame()
result.r2 <- data.frame()
result.r3 <- data.frame()
for (i in names(data)) {
  mod <- lmer(
    LOG_FungiOTU100.match ~ 
      Design*Ecosystem+ (1|Plot),
    data = data[[i]]
  )
  a<-performance::performance(mod)
  a$group <- i
  b<-parameters::parameters(mod) 
  b$group <- i
  d<-eta_squared(mod, partial = TRUE)
  c<-anova(mod)
  result <- rbind(result,b)
  result.r <- rbind(result.r,a[,c(1:5,ncol(a))])
  result.r2 <- rbind(result.r2,d)
  result.r3 <- rbind(result.r3,c)
}
write.csv(result.r,"ecology2.rsquare.fungi.100.richness.csv")
write.csv(result.r3,"ecology2.significance.fungi.100.richness.csv")
write.csv(result.r2,"ecology2.rsquare.fungi.100.richness.each.factor.csv")
write.csv(result,"ecology2.parameters.fungi.100.richness.csv")
##fungal shannon
result <- data.frame()
result.r <- data.frame()
result.r2 <- data.frame()
result.r3 <- data.frame()
for (i in names(data)) {
  mod <- lmer(
    FungiShannon ~ 
      Design*Ecosystem + LOG_FungiSeq+ (1|Plot),
    data = data[[i]]
  )
  a<-performance::performance(mod)
  a$group <- i
  b<-parameters::parameters(mod) 
  b$group <- i
  d<-eta_squared(mod, partial = TRUE)
  c<-anova(mod)
  result <- rbind(result,b)
  result.r <- rbind(result.r,a[,c(1:5,ncol(a))])
  result.r2 <- rbind(result.r2,d)
  result.r3 <- rbind(result.r3,c)
}
write.csv(result.r,"ecology2.rsquare.fungi.shannon.csv")
write.csv(result.r3,"ecology2.significance.fungi.shannon.csv")
write.csv(result.r2,"ecology2.rsquare.fungi.shannon.each.factor.csv")
write.csv(result,"ecology2.parameters.fungi.shannon.csv")
##fungal evenness
result <- data.frame()
result.r <- data.frame()
result.r2 <- data.frame()
result.r3 <- data.frame()
for (i in names(data)) {
  mod <- lmer(
    FungiEvenness ~ 
      Design*Ecosystem + LOG_FungiSeq+ (1|Plot),
    data = data[[i]]
  )
  a<-performance::performance(mod)
  a$group <- i
  b<-parameters::parameters(mod) 
  b$group <- i
  d<-eta_squared(mod, partial = TRUE)
  c<-anova(mod)
  result <- rbind(result,b)
  result.r <- rbind(result.r,a[,c(1:5,ncol(a))])
  result.r2 <- rbind(result.r2,d)
  result.r3 <- rbind(result.r3,c)
}
write.csv(result.r,"ecology2.rsquare.fungi.evenness.csv")
write.csv(result.r3,"ecology2.significance.fungi.evenness.csv")
write.csv(result.r2,"ecology2.rsquare.fungi.evenness.each.factor.csv")
write.csv(result,"ecology2.parameters.fungi.evenness.csv")

###Residuals of LOG_EukOTU
result <- data.frame()
result.r <- data.frame()
result.r2 <- data.frame()
result.r3 <- data.frame()
for (i in names(data)) {
  mod <- lmer(
    EUKotu.r ~ 
      Design*Ecosystem+ (1|Plot),
    data = data[[i]]
  )
  a<-performance::performance(mod)
  a$group <- i
  b<-parameters::parameters(mod) 
  b$group <- i
  d<-eta_squared(mod, partial = TRUE)
  c<-anova(mod)
  result <- rbind(result,b)
  result.r <- rbind(result.r,a[,c(1:5,ncol(a))])
  result.r2 <- rbind(result.r2,d)
  result.r3 <- rbind(result.r3,c)
}
write.csv(result.r,"ecology2.rsquare.fungi.euk.richness.csv")
write.csv(result.r3,"ecology2.significance.fungi.euk.richness.csv")
write.csv(result.r2,"ecology2.rsquare.fungi.euk.richness.each.factor.csv")
write.csv(result,"ecology2.parameters.fungi.euk.richness.csv")

##EUK genus richness
result <- data.frame()
result.r <- data.frame()
result.r2 <- data.frame()
result.r3 <- data.frame()
for (i in names(data)) {
  mod <- lmer(
    EUKge.r ~ 
      Design*Ecosystem+ (1|Plot),
    data = data[[i]]
  )
  a<-performance::performance(mod)
  a$group <- i
  b<-parameters::parameters(mod) 
  b$group <- i
  d<-eta_squared(mod, partial = TRUE)
  c<-anova(mod)
  result <- rbind(result,b)
  result.r <- rbind(result.r,a[,c(1:5,ncol(a))])
  result.r2 <- rbind(result.r2,d)
  result.r3 <- rbind(result.r3,c)
}
write.csv(result.r,"ecology2.rsquare.fungi.euk.ge.richness.csv")
write.csv(result.r3,"ecology2.significance.fungi.euk.ge.richness.csv")
write.csv(result.r2,"ecology2.rsquare.fungi.euk.ge.richness.each.factor.csv")
write.csv(result,"ecology2.parameters.fungi.euk.ge.richness.csv")

