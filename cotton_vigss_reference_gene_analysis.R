## Script Title: Analysis of cotton reference genes for cotton-herbivore VIGS studies
## Author: Mason Clark
## Date Started: 2024-06-27
## Date Last Modified: 2025-06-19 (Modified for repository)
## Description: This script performs a comprehensive analysis of qPCR data to identify
##              suitable, stable reference genes in Upland cotton under cotton aphid herbivory stress 
##              and VIGS treatments over time.
##              This script includes data loading, preparation of data frames for the RQdeltaCT 
##              package used for all stability analyses, data filtering/grooming, 
##              expression analysis using linear mixed-effects models and data bootstrapping.

#### Libraries and Setup ####

# Load necessary R packages.
library(ggplot2)     # Data visualization and figure ploting.
library(ReadqPCR)    # For reading and managing qPCR data.
library(NormqPCR)    # For normalization of qPCR data.
library(RQdeltaCT)   # For delta-Ct method calculations.
library(ctrlGene)    # Another package for reference gene analysis.
library(dplyr)       # A grammar of data manipulation, providing a consistent set of verbs.
library(tidyverse)   # An opinionated collection of R packages designed for data science.
library(purrr)       # A toolkit for functional programming.
library(lme4)        # For fitting and analyzing linear mixed-effects models.
library(lmerTest)    # Provides p-values in summary tables of lmer objects.
library(emmeans)     # Estimated Marginal Means, for post-hoc comparisons.
library(car)         # Companion to Applied Regression, for statistical diagnostics.
library(effectsize)  # For calculating effect sizes.
library(RankAggreg)  # For ranking aggregation methods.

#### directory setting ####
dir <- c("/path/to/directory")
setwd(dir)

#### data files ####

# Raw Cts df
rawcts <- read.csv("final_cqs_with_exp_metadata.csv", header = TRUE)

# Renaming column for analysis
colnames(rawcts)[4] <- c("target")


# Renaming gene names per reviewer request 04/09/25
rawcts$target <- dplyr::recode(rawcts$target, "ACT7" = "GhACT7",
                               "UBQ7" = "GhUBQ7",
                               "UBQ14" = "GhUBQ14",
                               "PP2A1" = "GhPP2A1",
                               "TMN5" = "GhTMN5",
                               "TBL6" = "GhTBL6",
                               "HYD1" = "GhHYDRA1")

#Changing level order for figures
rawcts$target <- factor(rawcts$target, levels = c("GhACT7", "GhUBQ7", "GhUBQ14", "GhPP2A1", "GhTMN5", "GhTBL6", "GhHYDRA1"))

#Now import using CQdeltaCT function read_long
data.Ct <- read_Ct_long(path = "CqdeltaCT_bygroup_all-cqs_withTBL.txt",
                        sep="\t",
                        dec = ".",
                        skip = 0,
                        column.Sample = 1,
                        column.Gene = 3,
                        column.Ct = 4,
                        column.Group = 5)

#Recoding gene names for reviewer
data.Ct$Gene <- dplyr::recode(data.Ct$Gene, "ACT7" = "GhACT7",
                               "UBQ7" = "GhUBQ7",
                               "UBQ14" = "GhUBQ14",
                               "PP2A1" = "GhPP2A1",
                               "TMN5" = "GhTMN5",
                               "TBL6" = "GhTBL6",
                               "HYD1" = "GhHYDRA1")

#Create analysis-ready dataframe
data.Ct.ready <- make_Ct_ready(data = data.Ct, imput.by.mean.within.groups = TRUE)



#Get summary table
data.Ct.summary <- data.Ct %>%
  group_by(Gene) %>%
  summarise(mean = mean(Ct),
            sd = sd(Ct),
            min = min(Ct),
            max = max(Ct))


#### Filtering raw CTs for standard deviations < 0.5 ####


# Function to calculate SD and filter technical replicates
filter_replicates <- function(cq_values) {
  # If there are any NA values, the SD is NA
  if (any(is.na(cq_values))) {
    return(data.frame(Cq = cq_values, SD = rep(NA, length(cq_values))))
  }
  
  # Calculate SD for all combinations of 2 out of 3
  combs <- combn(cq_values, 2, simplify = FALSE)
  combs_sd <- sapply(combs, sd)
  
  # Check if any combination of 2 has SD < 0.5
  if (any(combs_sd < 0.5)) {
    # Return the combination with SD < 0.5 and the calculated SD
    best_comb <- combs[[which.min(combs_sd)]]
    best_comb_sd <- rep(min(combs_sd), 3)
    best_comb_sd[3] <- NA  # Add NA for the removed replicate
    return(data.frame(Cq = c(best_comb, NA), SD = best_comb_sd))
  } else {
    # Calculate the SD for all 3 replicates
    sd_value <- sd(cq_values)
    return(data.frame(Cq = cq_values, SD = rep(sd_value, length(cq_values))))
  }
}

# Apply the function to each sample and target group
filtered_data <- rawcts %>%
  group_by(Sample, Target) %>%
  nest() %>%
  mutate(Filtered = map(data, ~ filter_replicates(.x$Cq))) %>%
  unnest(cols = c(Filtered))

# Select and rename columns as required
filtered_data <- filtered_data %>%
  select(Sample, Target, Cq, SD) %>%
  arrange(Sample, Target)

#Save the final raw Cq filtered data
write.csv(filtered_data, "Filtered_Ct_values.csv", row.names = FALSE)

#Remove NAs
filter_data <- na.omit(filtered_data)

#Get mean Cq values and show s.d.
filtered_data_summarized <- filter_data %>%
      group_by(Sample, Target) %>%
      summarize(Cq = mean(Cq),
                  s.d. = mean(SD))

#Remove samples where s.d. > 0.5
filtered_cq_for_analysis <- filtered_data_summarized[filtered_data_summarized$s.d. < 0.5, ]

#Final filtered cq file for analysis
#Need to add group data manually for CQdeltaCT package
write.csv(filtered_cq_for_analysis, "filtered_cq_for_analysis.csv")


## See which samples are missing a gene sample
## All samples need equal number of detectors for analysis using normQPCR package
sample_gene_summary <- rawcts %>%
                       group_by(Sample) %>%
                       summarize(Count = n())

write.csv(sample_gene_summary, "sample_gene_instances_counts.csv")


### FUNCTIONS ####

#Export figures function
export_figure <- function(plot, filename, width, height){
ggplot2::ggsave(plot = plot, filename = filename, width = width, height = height, units = "cm", dpi = 300)
}

#### Ct and dCt calculations ####

#I need to use the df prepared by the package for the analysis:
#data.Ct.ready
#New df for analysis
data.Ct.ready.deltact <- data.Ct.ready[,c(1:3, 5:9)]
data.Ct.ready.deltact <- subset(data.Ct.ready.deltact, data.Ct.ready.deltact$Group=="complete")


#DeltaCt for GhGhACT7
data.Ct.ready.deltact.ACT <- delta_Ct(data = data.Ct.ready.deltact,
                                      normalise = TRUE, 
                                      ref = c("GhPP2A1", "GhTMN5", "GhUBQ7", "GhUBQ14", "GhTBL6"),
                                      transform = FALSE)

#Need to update column names of deltacCt to gene names
data.Ct.ready.deltact.ACT <- data.Ct.ready.deltact.ACT %>% rename(GhACT7 = dCt)

#DeltaCt for GhPP2A1
data.Ct.ready.deltact.GhPP2A1 <- delta_Ct(data = data.Ct.ready.deltact,
                                      normalise = TRUE, 
                                      ref = c("GhACT7", "GhTMN5", "GhUBQ7", "GhUBQ14", "GhTBL6"),
                                      transform = FALSE)

#Need to update column names of deltacCt to gene names
data.Ct.ready.deltact.GhPP2A1 <- data.Ct.ready.deltact.GhPP2A1 %>% rename(GhPP2A1 = dCt)

#DeltaCt for TMN5
data.Ct.ready.deltact.TMN <- delta_Ct(data = data.Ct.ready.deltact,
                                        normalise = TRUE, 
                                        ref = c("GhPP2A1", "GhACT7", "GhUBQ7", "GhUBQ14", "GhTBL6"),
                                        transform = FALSE)

#Need to update column names of deltacCt to gene names
data.Ct.ready.deltact.TMN <- data.Ct.ready.deltact.TMN %>% rename(GhTMN5 = dCt)

#DeltaCt for UBQ7
data.Ct.ready.deltact.UBQ7 <- delta_Ct(data = data.Ct.ready.deltact,
                                      normalise = TRUE, 
                                      ref = c("GhPP2A1", "GhTMN5", "GhACT7", "GhUBQ14", "GhTBL6"),
                                      transform = FALSE)

#Need to update column names of deltacCt to gene names
data.Ct.ready.deltact.UBQ7 <- data.Ct.ready.deltact.UBQ7 %>% rename(GhUBQ7 = dCt)

#DeltaCt for UBQ14
data.Ct.ready.deltact.UBQ14 <- delta_Ct(data = data.Ct.ready.deltact,
                                       normalise = TRUE, 
                                       ref = c("GhPP2A1", "GhTMN5", "GhACT7", "GhUBQ7", "GhTBL6"),
                                       transform = FALSE)

#Need to update column names of deltacCt to gene names
data.Ct.ready.deltact.UBQ14 <- data.Ct.ready.deltact.UBQ14 %>% rename(GhUBQ14 = dCt)

#DeltaCt for TBL6
data.Ct.ready.deltact.TBL <- delta_Ct(data = data.Ct.ready.deltact,
                                        normalise = TRUE, 
                                        ref = c("GhPP2A1", "GhTMN5", "GhACT7", "GhUBQ7", "GhUBQ14"),
                                        transform = FALSE)


#Need to update column names of deltacCt to gene names
data.Ct.ready.deltact.TBL <- data.Ct.ready.deltact.TBL %>% rename(GhTBL6 = dCt)

#Combine to final df
deltact.df <- data.Ct.ready.deltact.ACT %>%
            full_join(data.Ct.ready.deltact.GhPP2A1, by = c("Sample", "Group")) %>%
            full_join(data.Ct.ready.deltact.TMN, by = c("Sample", "Group")) %>%
            full_join(data.Ct.ready.deltact.UBQ7, by = c("Sample", "Group")) %>%
            full_join(data.Ct.ready.deltact.UBQ14, by = c("Sample", "Group")) %>%
            full_join(data.Ct.ready.deltact.TBL, by = c("Sample", "Group"))


#Final dCt data frame
deltact.df <- deltact.df %>%
              pivot_longer(cols = c(GhACT7, GhPP2A1, GhTMN5, GhUBQ7, GhUBQ14, GhTBL6),
                          names_to = "target", values_to = "dCT") 

#Now I need to join delta.df with rawct to add in the deltaCt value
raw.delta <- inner_join(rawcts, deltact.df, by=c("Sample","target"))

#Changing level order
raw.delta$target <- factor(raw.delta$target, levels = c("GhACT7", "GhUBQ7", "GhUBQ14", "GhPP2A1", "GhTMN5", "GhTBL6", "GhHYDRA1"))


# I first want to see if I can use HYD1 plants
# by making sure silencing hyd1 didn't have an effect on cts of reference genes
# Subset Cts by gene
tmn <- subset(raw.delta, raw.delta$target=="TMN5") # Do not need to recode names here
act <- subset(raw.delta, raw.delta$target=="GhACT7") # Do not need to recode names here
pp <- subset(raw.delta, raw.delta$target=="GhPP2A1") # Do not need to recode names here
U7 <- subset(raw.delta, raw.delta$target=="UBQ7") # Do not need to recode names here
U14 <- subset(raw.delta, raw.delta$target=="UBQ14") # Do not need to recode names here
tbl <- subset(raw.delta, raw.delta$target=="TBL6") # Do not need to recode names here


# Now determine whether the gene target (none, GFP, or Hydra1 had an effect on reference gene Cts)
lm1 <- lm(Cq ~ gene_target, data = tmn)
aov1 <- aov(lm1)
summary(aov1) #SILENCING GHHYDRA1 HAD NO EFFECT ON TMN5 delta-Ct

lm2 <- lm(Cq ~ gene_target, data = act)
aov2 <- aov(lm2)
summary(aov2) #SILENCING GHHYDRA1 HAD NO EFFECT ON GhACT7 delta-Ct

lm3 <- lm(Cq ~ gene_target, data = pp)
aov3 <- aov(lm3)
summary(aov3) #SILENCING GHHYDRA1 HAD NO EFFECT ON GhPP2A1 delta-Ct

lm4 <- lm(Cq ~ gene_target, data = U7)
aov4 <- aov(lm4)
summary(aov4) #SILENCING GHHYDRA1 HAD NO EFFECT ON UBQ7 delta-Ct

lm5 <- lm(Cq ~ gene_target, data = U14)
aov5 <- aov(lm5)
summary(aov5) #SILENCING GHHYDRA1 HAD NO EFFECT ON UBQ14 delta-Ct 

lm6 <- lm(Cq ~ gene_target, data = tbl)
aov6 <- aov(lm6)
summary(aov6) #SILENCING GHHYDRA1 HAD NO EFFECT ON TBL6 delta-Ct 

## RAW CTs PLOT
all_samples <- ggplot(subset(rawcts, target != "GhHYDRA1"), aes(x=target, y=Cq)) +
                geom_boxplot(lwd = 1.0, fill = "grey", alpha = 0.4) +
                xlab('Candidate reference gene') +
                ylab("Ct value") +
                theme_classic() +
                theme(axis.text=element_text(size=18,  color="black"),
                      axis.text.x=element_text(angle = 45, vjust=1, hjust=1, face="italic"),
                      legend.position="none",
                      axis.line=element_line(size=1.0),
                      axis.ticks=element_line(size=1.0),
                      axis.ticks.length=unit(.25, "cm"),
                      axis.title=element_text(size=18,face="bold"))


#Plotting Cts separated by the two time points
# GhACT7 is *** significance
# GhPP2A1 is ***
# TBL is ** significance
# Figure updated 
bytime <- ggplot(subset(raw.delta, target != "GhHYDRA1"), aes(x=target, y=dCT, fill=dpi)) +
            geom_boxplot(lwd = 1.5) +
            xlab('Candidate reference gene') +
            scale_fill_manual(name="Days-post infestation", labels=c("14 DPI", "21 DPI"), values=c("grey36", "grey76")) +
            scale_y_continuous(name = "dCt", limits=c(-10,10)) +
            theme_classic() +
            theme(axis.text=element_text(size=18,  color="black"),
                    axis.text.x=element_text(angle = 45, vjust=1, hjust=1, , face="italic"),
                    legend.position="none",
                    axis.line=element_line(size=1.0),
                    axis.ticks=element_line(size=1.0),
                    axis.ticks.length=unit(.25, "cm"),
                    axis.title=element_text(size=18,face="bold")) +
            geom_line(data=tibble(x=c(0.80,1.20), y=c(2.0,2.0)), aes(x=x,y=y), size = 2.0, inherit.aes = FALSE) + #GhACT7
            geom_text(data=tibble(x=1.0, y=2.25), aes(x=x,y=y, label= "***", fontface="bold"), size = 12, inherit.aes=FALSE) + #GhACT7
            geom_line(data=tibble(x=c(5.80,6.20), y=c(8.25,8.25)), aes(x=x,y=y), size = 2.0, inherit.aes = FALSE) + #TBL6
            geom_text(data=tibble(x=6.0, y=8.50), aes(x=x,y=y, label= "**", fontface="bold"), size = 12, inherit.aes=FALSE) + #TBL6
            geom_line(data=tibble(x=c(3.80,4.20), y=c(8.25,8.25)), aes(x=x,y=y), size = 2.0, inherit.aes = FALSE) + #GhPP2A1
            geom_text(data=tibble(x=4.0, y=8.50), aes(x=x,y=y, label= "***", fontface="bold"), size = 12, inherit.aes=FALSE) #GhPP2A1
  
  
#Cts plot separating groups VIGS or WT
infiltration_status <- ggplot(subset(raw.delta, target != "GhHYDRA1"), aes(x=target, y=dCT, fill=vigs_condition)) +
                          geom_boxplot(lwd = 1.5) +
                          xlab('Candidate reference gene') +
                          scale_fill_manual(name="Infiltration status", labels=c("Uninfiltrated", "VIGS"), values=c("grey36", "grey76")) +
                          ylab("dCt") +
                          theme_classic() +
                          theme(axis.text=element_text(size=18,  color="black"),
                                axis.text.x=element_text(angle = 45, vjust=1, hjust=1, face="italic"),
                                legend.position="none",
                                axis.line=element_line(size=1.0),
                                axis.ticks=element_line(size=1.0),
                                axis.ticks.length=unit(.25, "cm"),
                                axis.title=element_text(size=18,face="bold"))

#Cts plot separating groups by whether or not they were infested by aphids
#Just TMN **
infestation_status <- ggplot(subset(raw.delta, target != "GhHYDRA1"), aes(x=target, y=dCT, fill=factor(aphid_status, levels=c("uninfested","infested")))) +
                        geom_boxplot(lwd = 1.5) +
                        xlab('Candidate reference gene') +
                        scale_fill_manual(name="Infestation status", labels=c("Uninfested", "Infested"), values=c("grey36", "grey76")) +
                        scale_y_continuous(name = "dCt", limits=c(-10,10)) +
                        theme_classic() +
                        theme(axis.text=element_text(size=18,  color="black"),
                              axis.text.x=element_text(angle = 45, vjust=1, hjust=1, face="italic"),
                              legend.position="none",
                              axis.line=element_line(size=1.0),
                              axis.ticks=element_line(size=1.0),
                              axis.ticks.length=unit(.25, "cm"),
                              axis.title=element_text(size=18,face="bold")) +
                        geom_line(data=tibble(x=c(4.80,5.20), y=c(8.25,8.25)), aes(x=x,y=y), size = 2.0, inherit.aes = FALSE) +
                        geom_text(data=tibble(x=5.0, y=8.5), aes(x=x,y=y, label="**", fontface="bold"), size = 12, inherit.aes=FALSE)



#Analysis for significant differences in Ct values between treatments

#lmer on conditions 

#TMN5
tmn_model <- lmer(dCT ~ aphid_status + vigs_condition + dpi + (1|plant.ID), dat = tmn)
summary(tmn_model) #aphid status p=0.00056

# Checking model residuls for normality
plot(density(residuals(tmn_model))) #Looks fine
shapiro.test(residuals(tmn_model)) #Looks fine

#Checking homoscedasticity)
tmn_residuals <- residuals(tmn_model)
tmn_interaction <- interaction(tmn$aphid_status, tmn$vigs_condition, tmn$dpi)
leveneTest(tmn_residuals ~ tmn_interaction, center = "median") #p > 0.05

#GhACT7
GhACT7_model <- lmer(dCT ~ aphid_status + vigs_condition + dpi + (1|plant.ID), dat = act_subset)
summary(GhACT7_model) #dpi p < 0.0001

# Pull out residuals
act_residuals <- residuals(GhACT7_model)

# Checking homoscedasticity
act_interaction <- interaction(act$aphid_status, act$vigs_condition, act$dpi)
leveneTest(act_residuals ~ act_interaction, center = "median") #p > 0.05

#GhPP2A1
pp_model <- lmer(dCT ~ aphid_status + vigs_condition + dpi + (1 | plant.ID),
                 dat = pp)
summary(pp_model) #dpi p < 0.0001

#Checking homoscedasticity
pp_residuals <- residuals(pp_model) #residuals not normal
pp_interaction <- interaction(pp$aphid_status, pp$vigs_condition, pp$dpi)
leveneTest(pp_residuals ~ pp_interaction, center = "median") #p > 0.05

#UBQ7
u7_model <- lmer(dCT ~ aphid_status + vigs_condition + dpi + (1 | plant.ID),
                 dat = U7)
summary(u7_model) #none

#Checking homoscedasticity
u7_residuals <- residuals(u7_model) #bimodal
u7_interaction <- interaction(U7$aphid_status, U7$vigs_condition, U7$dpi)
leveneTest(u7_residuals ~ u7_interaction, center = "median") #p > 0.05

#UBQ14
u14_model <- lmer(dCT ~ aphid_status + vigs_condition + dpi + (1|plant.ID), dat = U14)
summary(u14_model) # none

#Checking homoscedasticity
u14_residuals <- residuals(u14_model) #Bimodal distribution of residuals
u14_interaction <- interaction(U14$aphid_status, U14$vigs_condition, U14$dpi)
leveneTest(u14_residuals ~ u14_interaction, center = "median") #p > 0.05

#TBL6
tbl_model <- lmer(dCT ~ aphid_status + vigs_condition + dpi + (1|plant.ID), dat = tbl)
summary(tbl_model) #dpi p=0.00369

#Checking homoscedasticity
tbl_residuals <- residuals(tbl_model) #Normal enough
tbl_interaction <- interaction(tbl$aphid_status, tbl$vigs_condition, tbl$dpi)
leveneTest(tbl_residuals ~ tbl_interaction, center = "median") #p > 0.05


#Need to perform p.adj for multiple testings corrections

#aphid status
aphid.ps <- c(0.00056, 0.731, 0.690, 0.379009, 0.0772, 0.61247)
aphid.p.adj <- p.adjust(aphid.ps, method = "BH")

#time
time.ps <- c(0.10427, 3.28e-08, 8.24e-05, 0.618747, 0.7288, 0.00369)
time.p.adj <- p.adjust(time.ps, method = "BH")

#Export Ct plots
export_figure(all_samples, "all_sample_cts.jpeg", 20, 16)
export_figure(bytime, "cts_by_time.jpeg", 20, 16)
export_figure(infiltration_status, "cts_by_infiltration_status.jpeg", 20, 16)
export_figure(infestation_status, "cts_by_infestation_status.jpeg", 20, 16)

##### geNorm M Stability Figures (data from CqdeltaCt) #####

#M stability data.frame
mstabdat <- read.csv("geNORM_ranking_results.csv", header = TRUE)
mstabdat$gene <- dplyr::recode(mstabdat$gene, "ACT7" = "GhACT7",
                               "UBQ7" = "GhUBQ7",
                               "UBQ14" = "GhUBQ14",
                               "ACT7" = "GhPP2A1",
                               "TMN5" = "GhTMN5",
                               "TBL6" = "GhTBL6",
                               "HYD1" = "GhHYDRA1",
                               "ACT7/PP2A1" = "GhACT7/GhPP2A1",
                               "PP2A1/ACT7" = "GhPP2A1/GhACT7")

## First plot sets have both conditions on same figure to compare differences in stability scores down the rank

# All conditions in one plot
# Did not use plot for manuscript
# mstab_plot1 <- ggplot(mstabdat, aes(x=number.remaining, y=m.score, group=condition)) +
#   geom_line(aes(linetype=condition), size = 2.0) +
#   scale_x_reverse() +
#   xlab("Remaining candidate reference genes") +
#   ylab("geNorm M stability score") +
#   geom_point(aes(shape=condition), size = 8) +
#   theme_classic() +
#   theme(axis.text=element_text(size=28, face = "bold", color="black"),
#         axis.line=element_line(size=1.5),
#         axis.ticks=element_line(size=1.5),
#         axis.ticks.length=unit(.25, "cm"),
#         legend.position = "top",
#         legend.text=element_text(size=rel(1.5)),
#         legend.title=element_text(size=rel(1.5)),
#         axis.title=element_text(size=34,face="bold"))

# #export plot
# export_figure(mstab_plot1, "mstability_scores_all_conditions_on_same_plot.jpeg", 26, 20)

#Data from just comparing infiltration status
mstab.infiltrationstatus.dat <- subset(mstabdat, mstabdat$contrast == "vigs_status")

#Plot stability comparing conditions
mstab_plot2 <- ggplot(mstab.infiltrationstatus.dat, aes(x=genes.left, y=genorm.score, linetype=condition, shape=condition)) +
  geom_line(size = 1.0) +
  scale_x_reverse() +
  xlab("Remaining candidate reference genes") +
  ylab("geNorm stability score (M)") +
  scale_y_continuous(limits = c(0,2.5), breaks = c(0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.50)) +
  scale_shape_discrete(name="Infiltration status", labels=c("Uninfiltrated", "VIGS")) +
  scale_linetype_discrete(name="Infiltration status", labels=c("Uninfiltrated", "VIGS")) +
  geom_point(size = 4) +
  theme_classic() +
  theme(axis.text=element_text(size=12, color="black"),
        axis.line=element_line(size=1.0),
        axis.ticks=element_line(size=1.0),
        axis.ticks.length=unit(.25, "cm"),
        legend.position = "top",
        legend.text=element_text(size=rel(1)),
        legend.title=element_text(size=rel(1)),
        axis.title=element_text(size=12,face="bold"))



#Export
export_figure(mstab_plot2, "mstability_scores_infiltratin_status.jpeg", 12, 14)

#### START HERE

#Repeat for infested plants
mstab.infestationstatus.dat <- subset(mstabdat, mstabdat$contrast=="infestation_status")

#Plot based on infestation status
mstab_plot3 <- ggplot(mstab.infestationstatus.dat, aes(x=genes.left, y=genorm.score, linetype=condition, shape=condition)) +
                      geom_line(size = 1.0) +
                      scale_x_reverse() +
                      xlab("Remaining candidate reference genes") +
                      ylab("geNorm stability score (M)") +
              scale_y_continuous(limits = c(0,2.5), breaks = c(0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.50)) +
                                  scale_shape_discrete(name="Infestation status", labels=c("Infested", "Uninfested")) +
                                  scale_linetype_discrete(name="Infestation status", labels=c("Infested", "Uninfested")) +
                                  geom_point(size = 4) +
                                  theme_classic() +
              theme(axis.text=element_text(size=12, color="black"),
                    axis.line=element_line(size=1.0),
                    axis.ticks=element_line(size=1.0),
                    axis.ticks.length=unit(.25, "cm"),
                    legend.position = "top",
                    legend.text=element_text(size=rel(1)),
                    legend.title=element_text(size=rel(1)),
                    axis.title=element_text(size=12,face="bold"))
                      


#Infestation status plot
export_figure(mstab_plot3, "mstability_scores_infestation_status.jpeg", 12, 14)


#Comparing Mstab between time points
#Repeat for infested plants
mstab.dpi.dat <- subset(mstabdat, mstabdat$contrast=="bytime")

#Plot based on infestation status
mstab_plot4 <- ggplot(mstab.dpi.dat, aes(x = genes.left,
                      y = genorm.score, linetype = condition, shape = condition)) +
                      geom_line(size = 1.0) +
                      scale_x_reverse() +
                xlab("Remaining candidate reference genes") +
                ylab("geNorm stability score (M)") +
                scale_y_continuous(limits = c(0,2.5), breaks = c(0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.50)) +
                            scale_shape_discrete(name="Days post-infestation", labels=c("14 DPI", "21 DPI")) +
                            scale_linetype_discrete(name="Days post-infestation", labels=c("14 DPI", "21 DPI")) +
                            geom_point(size = 4) +
                            theme_classic() +
              theme(axis.text=element_text(size=12, color="black"),
                    axis.line=element_line(size=1.0),
                    axis.ticks=element_line(size=1.0),
                    axis.ticks.length=unit(.25, "cm"),
                    legend.position = "top",
                    legend.text=element_text(size=rel(1)),
                    legend.title=element_text(size=rel(1)),
                    axis.title=element_text(size=12,face="bold"))

#Infestation status plot
export_figure(mstab_plot4, "mstability_scores_bytime.jpeg", 12, 14)




## Second plot sets will now show ranking of genes on condition-by basis


# DPI - 14
generanking_14dpi <- ggplot(subset(mstab.dpi.dat, condition %in% "14dpi"), 
       aes(x=factor(gene, levels=c("GhUBQ7", "GhUBQ14", "GhTBL6", "GhTMN5", "GhACT7/GhPP2A1")), y=genorm.score, group = condition)) +
      geom_line(size=1.0) +
      xlab("Gene ranking") +
      ylab("geNorm stability score (M)") +
      scale_y_continuous(limits = c(0,2), breaks = c(0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0)) +
      geom_point(size = 4) +
      theme_classic() +
      theme(axis.text=element_text(size=12, color="black"),
            axis.text.x=element_text(angle=45, vjust=1.0, hjust=1.0, face="italic"),
            axis.line=element_line(size=1.0),
            axis.ticks=element_line(size=1.0),
            axis.ticks.length=unit(.25, "cm"),
            legend.position = "top",
            legend.text=element_text(size=rel(2.5)),
            legend.title=element_text(size=rel(3)),
            axis.title=element_text(size=12,face="bold")) +
      geom_text(data=tibble(x=4.5, y=2.0), aes(x=x,y=y, label= "14 DPI", fontface="bold"), size = 8, inherit.aes=FALSE)


#Export 14dpi gene ranking figure
export_figure(generanking_14dpi, "mstability_scores_14dpi.jpeg", 12, 14)

# DPI - 21
generanking_21dpi <- ggplot(subset(mstab.dpi.dat, condition %in% "21dpi"), 
                            aes(x=factor(gene, levels=c("GhUBQ7", "GhUBQ14", "GhTBL6", "GhTMN5", "GhACT7/GhPP2A1")), y=genorm.score, group = condition)) +
                    geom_line(size=1.0) +
                    xlab("Gene ranking") +
                    ylab("geNorm stability score (M)") +
                    scale_y_continuous(limits = c(0,2.5), breaks = c(0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.50)) +
                    geom_point(size = 4) +
                    theme_classic() +
                    theme(axis.text=element_text(size=12, color="black"),
                          axis.text.x=element_text(angle=45, vjust=1.0, hjust=1.0, face="italic"),
                          axis.line=element_line(size=1.0),
                          axis.ticks=element_line(size=1.0),
                          axis.ticks.length=unit(.25, "cm"),
                          legend.position = "top",
                          legend.text=element_text(size=rel(2.5)),
                          legend.title=element_text(size=rel(3)),
                          axis.title=element_text(size=12,face="bold")) +
                    geom_text(data=tibble(x=4.5, y=2.5), aes(x=x,y=y, label= "21 DPI", fontface="bold"), size = 8, inherit.aes=FALSE)

#Export 21dpi gene ranking figure
export_figure(generanking_21dpi, "mstability_scores_21dpi.jpeg", 12, 14)


#By infiltration treatment

#Uninfiltrated
generanking_uninfiltrated <- ggplot(subset(mstabdat, condition %in% "uninfiltrated"), 
                              aes(x=factor(gene, levels=c("GhUBQ7", "GhUBQ14", "GhTBL6", "GhTMN5", "GhACT7/GhPP2A1")), y=genorm.score, group = condition)) +
                              geom_line(size=1.0) +
                              xlab("Gene ranking") +
                              ylab("geNorm stability score (M)") +
                              scale_y_continuous(limits = c(0,2.5), breaks = c(0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25, 2.5)) +
                              geom_point(size = 4) +
                              theme_classic() +
                              theme(axis.text=element_text(size=12, color="black"),
                                    axis.text.x=element_text(angle=45, vjust=1.0, hjust=1.0, face="italic"),
                                    axis.line=element_line(size=1.0),
                                    axis.ticks=element_line(size=1.0),
                                    axis.ticks.length=unit(.25, "cm"),
                                    legend.position = "top",
                                    legend.text=element_text(size=rel(2.5)),
                                    legend.title=element_text(size=rel(3)),
                                    axis.title=element_text(size=12,face="bold")) +
                              geom_text(data=tibble(x=4.4, y=2.5), aes(x=x,y=y, label= "Uninfiltrated", fontface="bold"), size = 8, inherit.aes=FALSE)

export_figure(generanking_uninfiltrated, "mstability_scores_uninfiltrated.jpeg", 12, 14)

#VIGS
generanking_VIGS <- ggplot(subset(mstabdat, condition %in% "vigs-infiltrated"), 
                                    aes(x=factor(gene, levels=c("GhUBQ7", "GhUBQ14", "GhTBL6", "GhTMN5", "GhACT7/GhPP2A1")), y=genorm.score, group = condition)) +
                    geom_line(size=1.0) +
                    xlab("Gene ranking") +
                    ylab("geNorm stability score (M)") +
                    scale_y_continuous(limits = c(0,2.5), breaks = c(0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25, 2.5)) +
                    geom_point(size = 4) +
                    theme_classic() +
                    theme(axis.text=element_text(size=12, color="black"),
                        axis.text.x=element_text(angle=45, vjust=1.0, hjust=1.0, face="italic"),
                        axis.line=element_line(size=1.0),
                        axis.ticks=element_line(size=1.0),
                        axis.ticks.length=unit(.25, "cm"),
                        legend.position = "top",
                        legend.text=element_text(size=rel(2.5)),
                        legend.title=element_text(size=rel(3)),
                        axis.title=element_text(size=12,face="bold")) +
                  geom_text(data=tibble(x=5.0, y=2.5), aes(x=x,y=y, label= "VIGS", fontface="bold"), size = 8, inherit.aes=FALSE)

export_figure(generanking_VIGS, "mstability_scores_VIGS.jpeg", 12, 14)

#Infestation status

#Uninfested
generanking_uninfested <- ggplot(subset(mstabdat, condition %in% "uninfested"), 
                           aes(x=factor(gene, levels=c("GhUBQ7", "GhUBQ14", "GhTBL6", "GhTMN5", "GhPP2A1/GhACT7")), y=genorm.score, group = condition)) +
                            geom_line(size=1.0) +
                            xlab("Gene ranking") +
                            ylab("geNorm stability score (M)") +
                            scale_y_continuous(limits = c(0,2.5), breaks = c(0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25, 2.5)) +
                            geom_point(size = 4) +
                            theme_classic() +
                          theme(axis.text=element_text(size=12, color="black"),
                                axis.text.x=element_text(angle=45, vjust=1.0, hjust=1.0, face="italic"),
                                axis.line=element_line(size=1.0),
                                axis.ticks=element_line(size=1.0),
                                axis.ticks.length=unit(.25, "cm"),
                                legend.position = "top",
                                legend.text=element_text(size=rel(2.5)),
                                legend.title=element_text(size=rel(3)),
                                axis.title=element_text(size=12,face="bold")) +
                            geom_text(data=tibble(x=4.5, y=2.5), aes(x=x,y=y, label= "Uninfested", fontface="bold"), size = 8, inherit.aes=FALSE)

export_figure(generanking_uninfested, "mstability_scores_uninfested.jpeg", 12, 14)

#Infested
generanking_infested <- ggplot(subset(mstabdat, condition %in% "infested"), 
                                 aes(x=factor(gene, levels=c("GhUBQ7", "GhUBQ14", "GhTBL6", "GhTMN5", "GhACT7/GhPP2A1")), y=genorm.score, group = condition)) +
                        geom_line(size=1.0) +
                        xlab("Gene ranking") +
                        ylab("geNorm stability score (M)") +
                        scale_y_continuous(limits = c(0,2.5), breaks = c(0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25, 2.5)) +
                        geom_point(size = 4) +
                        theme_classic() +
                        theme(axis.text=element_text(size=12, color="black"),
                              axis.text.x=element_text(angle=45, vjust=1.0, hjust=1.0, face = "italic"),
                              axis.line=element_line(size=1.0),
                              axis.ticks=element_line(size=1.0),
                              axis.ticks.length=unit(.25, "cm"),
                              legend.position = "top",
                              legend.text=element_text(size=rel(2.5)),
                              legend.title=element_text(size=rel(3)),
                              axis.title=element_text(size=12,face="bold")) +
                        geom_text(data=tibble(x=4.7, y=2.5), aes(x=x,y=y, label= "Infested", fontface="bold"), size = 8, inherit.aes=FALSE)


export_figure(generanking_infested, "mstability_scores_infested.jpeg", 12, 14)

#All conditions

#Infested
generanking_allconditions <- ggplot(subset(mstabdat, condition %in% "all"), 
                               aes(x=factor(gene, levels=c("UBQ7", "UBQ14", "TBL6", "TMN5", "GhACT7/GhPP2A1")), y=genorm.score, group = condition)) +
                            geom_line(size=2.0) +
                            xlab("Gene ranking") +
                            ylab("geNorm stability score (M)") +
                            scale_y_continuous(limits = c(0,2.5), breaks = c(0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25, 2.5)) +
                            geom_point(size = 8) +
                            theme_classic() +
                            theme(axis.text=element_text(size=28, face = "bold"),
                                  axis.text.x=element_text(angle=45, vjust=1.0, hjust=1.0, color="black"),
                                  axis.line=element_line(size=1.5),
                                  axis.ticks=element_line(size=1.5),
                                  axis.ticks.length=unit(.25, "cm"),
                                  legend.position = "top",
                                  legend.text=element_text(size=rel(2.5)),
                                  legend.title=element_text(size=rel(3)),
                                  axis.title=element_text(size=34,face="bold")) +
                          geom_text(data=tibble(x=4.2, y=2.4), aes(x=x,y=y, label= "All conditions", fontface="bold"), size = 12, inherit.aes=FALSE)

#Export that shit
export_figure(generanking_allconditions, "mstability_scores_all_conditions.jpeg", 22, 24)



#### NORMFINDER and geNORM using CQDELTACT PACKAGE ###

## I'll first analyze using geNorm of individual groups
## Then I'll run again comparing conditions

## GeNorm analysis

#All samples first
all.conditions_genorm_normfinder <- find_ref_gene(data = data.Ct.ready[,1:9],
              groups = c("complete"),
              candidates = c("GhACT7", "GhPP2A1","TMN5","UBQ7","UBQ14", "TBL6"),
              col = c("#66c2a5", "#fc8d62","#6A6599", "#D62728", "black", "darkblue"),
              angle = 90,
              dpi = 600,
              width = 15, 
              height = 15,
              save.to.tiff = FALSE,
              genorm.score = TRUE,
              norm.finder.score = TRUE)

#Saving plot as a new object
all.conditions <- all.conditions_genorm_normfinder[[1]]

#Final all.conditions Ct plot
all.conditions.plot <- all.conditions +
                        theme_classic() +
                        scale_x_discrete(labels=paste0("S", as.numeric(factor(data.Ct.ready$Sample)))) + 
                        xlab("Experiment samples") +
                        ylab("Ct value") + 
                        theme(axis.text.y=element_text(size=26, face = "bold", color="black"),
                              axis.text.x=element_text(size=12, face = "bold", color="black"),
                              axis.line=element_line(size=1),
                              axis.ticks=element_line(size=1),
                              axis.ticks.length=unit(.25, "cm"),
                              legend.position = "top",
                              legend.text=element_text(size=rel(2)),
                              axis.title=element_text(size=34,face="bold"))

#Save all.conditions export plot
export_figure(all.conditions.plot, "all_conditions_Cts_plot.jpeg",40,20)


#Now geNorm for all uninfiltrated plants
uninfiltrated_genorm_normfinder <- find_ref_gene(data = data.Ct.ready[,1:9],
                                                  groups = c("uninfiltrated"),
                                                  candidates = c("GhACT7", "GhPP2A1","TMN5","UBQ7","UBQ14", "TBL6"),
                                                  col = c("#66c2a5", "#fc8d62","#6A6599", "#D62728", "black", "darkblue"),
                                                  angle = 45,
                                                  dpi = 600,
                                                  width = 15, 
                                                  height = 15,
                                                  save.to.tiff = FALSE,
                                                  genorm.score = TRUE,
                                                  norm.finder.score = TRUE)




#Final uninfiltrated Ct plot
uninfiltrated_genorm_normfinder.ctplot <- uninfiltrated_genorm_normfinder[[1]] +
                                          theme_classic() +
                                          xlab("Experiment samples") +
                                          scale_x_discrete(labels = setNames(data.Ct.ready$sample.id, data.Ct.ready$Sample)) + 
                                          ylab("Ct value") +
                                          theme(axis.text.y=element_text(size=26, face = "bold", color="black"),
                                                axis.text.x=element_text(size=18, color="black", angle = 45),
                                                axis.line=element_line(size=1),
                                                axis.ticks=element_line(size=1),
                                                axis.ticks.length=unit(.25, "cm"),
                                                legend.position = "top",
                                                legend.text=element_text(size=rel(2)),
                                                axis.title=element_text(size=34,face="bold"))

#Export uninfiltrated plot
export_figure(uninfiltrated_genorm_normfinder.ctplot, "uninfiltrated_Cts_plot.jpeg", 40,20)

#Now VIGS-infiltrated
vigs.infiltrated_genorm_normfinder <- find_ref_gene(data = data.Ct.ready[,1:9],
                                                 groups = c("vigs_infiltrated"),
                                                 candidates = c("GhACT7", "GhPP2A1","TMN5","UBQ7","UBQ14", "TBL6"),
                                                 col = c("#66c2a5", "#fc8d62","#6A6599", "#D62728", "black", "darkblue"),
                                                 angle = 90,
                                                 dpi = 600,
                                                 width = 15, 
                                                 height = 15,
                                                 save.to.tiff = FALSE,
                                                 genorm.score = TRUE,
                                                 norm.finder.score = TRUE)

#VIGS-infiltrated Ct plot
vigs.infiltrated_genorm_normfinder.plot <- vigs.infiltrated_genorm_normfinder[[1]] +
                                            theme_classic() +
                                            xlab("Experiment samples") +
                                            ylab("Ct value") +
                                            scale_x_discrete(labels = setNames(data.Ct.ready$sample.id, data.Ct.ready$Sample)) +
                                            theme(axis.text.y=element_text(size=26, face = "bold", color="black"),
                                                  axis.text.x=element_text(size=18, color="black"),
                                                  axis.line=element_line(size=1),
                                                  axis.ticks=element_line(size=1),
                                                  axis.ticks.length=unit(.25, "cm"),
                                                  legend.position = "top",
                                                  legend.text=element_text(size=rel(2)),
                                                  axis.title=element_text(size=34,face="bold"))



#Export Vigs-infiltrated Ct plot
export_figure(vigs.infiltrated_genorm_normfinder.plot, "vigs-infiltrated_Cts_plot.jpeg", 50,20)


#Uninfested plants
uninfested_genorm_normfinder <- find_ref_gene(data = data.Ct.ready[,1:9],
                                                    groups = c("uninfested"),
                                                    candidates = c("GhACT7", "GhPP2A1","TMN5","UBQ7","UBQ14", "TBL6"),
                                                    col = c("#66c2a5", "#fc8d62","#6A6599", "#D62728", "black", "darkblue"),
                                                    angle = 45,
                                                    dpi = 600,
                                                    width = 15, 
                                                    height = 15,
                                                    save.to.tiff = FALSE,
                                                    genorm.score = TRUE,
                                                    norm.finder.score = TRUE)
#Uninfested plants plot
uninfested_genorm_normfinder.cqsplot <- uninfested_genorm_normfinder[[1]] +
                                        theme_classic() +
                                        xlab("Experiment samples") +
                                        ylab("Ct value") +
                                        scale_x_discrete(labels = setNames(data.Ct.ready$sample.id, data.Ct.ready$Sample)) +
                                        theme(axis.text.y=element_text(size=26, face = "bold", color="black"),
                                              axis.text.x=element_text(size=18, color="black"),
                                              axis.line=element_line(size=1),
                                              axis.ticks=element_line(size=1),
                                              axis.ticks.length=unit(.25, "cm"),
                                              legend.position = "top",
                                              legend.text=element_text(size=rel(2)),
                                              axis.title=element_text(size=34,face="bold"))


#Export uninfested plants Ct plot
export_figure(uninfested_genorm_normfinder.cqsplot, "uninfested_plants_Cts_plot.jpeg", 50,20)


#Uninfested plants
infested_genorm_normfinder <- find_ref_gene(data = data.Ct.ready[,1:9],
                                              groups = c("infested"),
                                              candidates = c("GhACT7", "GhPP2A1","TMN5","UBQ7","UBQ14", "TBL6"),
                                              col = c("#66c2a5", "#fc8d62","#6A6599", "#D62728", "black", "darkblue"),
                                              angle = 45,
                                              dpi = 600,
                                              width = 15, 
                                              height = 15,
                                              save.to.tiff = FALSE,
                                              genorm.score = TRUE,
                                              norm.finder.score = TRUE)




#Infested plants plot
infested_genorm_normfinder.cqsplot <- infested_genorm_normfinder[[1]] +
                                        theme_classic() +
                                        xlab("Experiment samples") +
                                        ylab("Ct value") +
                                        scale_x_discrete(labels = setNames(data.Ct.ready$sample.id, data.Ct.ready$Sample)) +
                                        theme(axis.text.y=element_text(size=26, face = "bold", color="black"),
                                              axis.text.x=element_text(size=18, color="black"),
                                              axis.line=element_line(size=1),
                                              axis.ticks=element_line(size=1),
                                              axis.ticks.length=unit(.25, "cm"),
                                              legend.position = "top",
                                              legend.text=element_text(size=rel(2)),
                                              axis.title=element_text(size=34,face="bold"))


#Export uninfested plants Ct plot
export_figure(infested_genorm_normfinder.cqsplot, "infested_plants_Cts_plot.jpeg", 50,20)


#14dpi plants
fourteendpi_genorm_normfinder <- find_ref_gene(data = data.Ct.ready[,1:9],
                                            groups = c("fourteen"),
                                            candidates = c("GhACT7", "GhPP2A1","TMN5","UBQ7","UBQ14", "TBL6"),
                                            col = c("#66c2a5", "#fc8d62","#6A6599", "#D62728", "black", "darkblue"),
                                            angle = 45,
                                            dpi = 600,
                                            width = 15, 
                                            height = 15,
                                            save.to.tiff = FALSE,
                                            genorm.score = TRUE,
                                            norm.finder.score = TRUE)
#21dpi plants
twentyonedpi_genorm_normfinder <- find_ref_gene(data = data.Ct.ready[,1:9],
                                               groups = c("twenty-one"),
                                               candidates = c("GhACT7", "GhPP2A1","TMN5","UBQ7","UBQ14", "TBL6"),
                                               col = c("#66c2a5", "#fc8d62","#6A6599", "#D62728", "black", "darkblue"),
                                               angle = 45,
                                               dpi = 600,
                                               width = 15, 
                                               height = 15,
                                               save.to.tiff = FALSE,
                                               genorm.score = TRUE,
                                               norm.finder.score = TRUE)


## NormFinder analysis
## This for comparing between groups, but CqdeltaCt has a separate Normfinder function
## To measure within-group gene stability

#VIGS condition
uninfiltated_v_infiltrated.normfinder <- find_ref_gene(data = data.Ct.ready[,1:9],
                                                       groups = c("uninfiltrated", "vigs_infiltrated"),
                                                       candidates = c("GhACT7", "GhPP2A1","GhTMN5","GhUBQ7","GhUBQ14", "GhTBL6"),
                                                       col = c("#66c2a5", "#fc8d62","#6A6599", "#D62728", "black", "darkblue"),
                                                       angle = 45,
                                                       dpi = 600,
                                                       width = 15, 
                                                       height = 15,
                                                       save.to.tiff = FALSE,
                                                       genorm.score = FALSE,
                                                       norm.finder.score = TRUE)

#Normfinder data
normfinderdf <- read.csv("normfinder_results.csv", header = TRUE)


# Recoding genenames for reviewers
normfinderdf$gene <- dplyr::recode(normfinderdf$gene, "ACT7" = "GhACT7",
                                                       "UBQ7" = "GhUBQ7",
                                                       "UBQ14" = "GhUBQ14",
                                                       "PP2A1" = "GhPP2A1",
                                                       "TMN5" = "GhTMN5",
                                                       "TBL6" = "GhTBL6",
                                                       "HYD1" = "GhHYDRA1")

#Normfinder figure for infiltration status
uninfiltated_v_infiltrated.normfinder.plot <- ggplot(subset(normfinderdf, conidition %in% "infiltration_status"), 
                                                     aes(x=factor(gene, levels=c("GhUBQ7", "GhUBQ14", "GhTBL6", "GhTMN5", "GhACT7", "GhPP2A1")), y=normfinder.score, group = conidition)) +
                                              geom_line(size=1.0) +
                                              xlab("Gene ranking") +
                                              ylab("NormFinder score") +
                                              geom_point(size = 4) +
                                              theme_classic() +
                                              theme(axis.text=element_text(size=12, color="black"),
                                                    axis.text.x=element_text(angle=45, vjust=1.0, hjust=1.0, face="italic"),
                                                    axis.line=element_line(size=1.0),
                                                    axis.ticks=element_line(size=1.0),
                                                    axis.ticks.length=unit(.25, "cm"),
                                                    legend.position = "top",
                                                    legend.text=element_text(size=rel(2.5)),
                                                    legend.title=element_text(size=rel(3)),
                                                    axis.title=element_text(size=12,face="bold")) +
                                              geom_text(data=tibble(x=4.5, y=0.5), aes(x=x,y=y, label= "Infiltration status", fontface="bold"), size = 6, inherit.aes=FALSE)


#Normfinder figure for infestation status
uninfested_v_infested.normfinder.plot <- ggplot(subset(normfinderdf, conidition %in% "infestation_status"), 
                                                     aes(x=factor(gene, levels=c("GhUBQ7", "GhUBQ14", "GhTBL6", "GhTMN5", "GhACT7", "GhPP2A1")), y=normfinder.score, group = conidition)) +
                                                geom_line(size=1.0) +
                                                xlab("Gene ranking") +
                                                ylab("NormFinder score") +
                                                geom_point(size = 4) +
                                                theme_classic() +
                                                theme(axis.text=element_text(size=12, color="black"),
                                                      axis.text.x=element_text(angle=45, vjust=1.0, hjust=1.0, face="italic"),
                                                      axis.line=element_line(size=1.0),
                                                      axis.ticks=element_line(size=1.0),
                                                      axis.ticks.length=unit(.25, "cm"),
                                                      legend.position = "top",
                                                      legend.text=element_text(size=rel(2.5)),
                                                      legend.title=element_text(size=rel(3)),
                                                      axis.title=element_text(size=12,face="bold")) +
                                                geom_text(data=tibble(x=4.5, y=0.5), aes(x=x,y=y, label= "Infestation status", fontface="bold"), size = 6, inherit.aes=FALSE)

#Normfinder figure for infestation status
time.normfinder.plot <- ggplot(subset(normfinderdf, conidition %in% "time"), aes(x=factor(gene, levels=c("GhTBL6", "GhUBQ7", "GhACT7", "GhUBQ14", "GhPP2A1", "GhTMN5")), y=normfinder.score, group = conidition)) +
                                                      geom_line(size=1.0) +
                                                      xlab("Gene ranking") +
                                                      ylab("NormFinder score") +
                                                      geom_point(size = 4) +
                                                      theme_classic() +
                                                      theme(axis.text=element_text(size=12, color="black"),
                                                            axis.text.x=element_text(angle=45, vjust=1.0, hjust=1.0, face="italic"),
                                                            axis.line=element_line(size=1.0),
                                                            axis.ticks=element_line(size=1.0),
                                                            axis.ticks.length=unit(.25, "cm"),
                                                            legend.position = "top",
                                                            legend.text=element_text(size=rel(2.5)),
                                                            legend.title=element_text(size=rel(3)),
                                                            axis.title=element_text(size=12,face="bold")) +
                                                      geom_text(data=tibble(x=4.5, y=0.7), aes(x=x,y=y, label="Time post-infestation", fontface="bold"), size = 6, inherit.aes=FALSE)


export_figure(uninfiltated_v_infiltrated.normfinder.plot, "infiltrationsstatus.normfinder.plot.jpeg", 12,10)
export_figure(uninfested_v_infested.normfinder.plot, "infestationsstatus.normfinder.plot.jpeg", 12,10)
export_figure(time.normfinder.plot, "time.normfinder.plot.jpeg", 12,10)

#Infestation condition
uninfested_v_infested.normfinder <- find_ref_gene(data = data.Ct.ready[,1:9],
                                                       groups = c("uninfested", "infested"),
                                                       candidates = c("GhACT7", "GhPP2A1","TMN5","UBQ7","UBQ14", "TBL6"),
                                                       col = c("#66c2a5", "#fc8d62","#6A6599", "#D62728", "black", "darkblue"),
                                                       angle = 45,
                                                       dpi = 600,
                                                       width = 15, 
                                                       height = 15,
                                                       save.to.tiff = FALSE,
                                                       genorm.score = FALSE,
                                                       norm.finder.score = TRUE)

#by time
#Infestation condition
fourteendpi_v_twentyonedpi.normfinder <- find_ref_gene(data = data.Ct.ready[,1:9],
                                                  groups = c("fourteen", "twenty-one"),
                                                  candidates = c("GhACT7", "GhPP2A1","TMN5","UBQ7","UBQ14", "TBL6"),
                                                  col = c("#66c2a5", "#fc8d62","#6A6599", "#D62728", "black", "darkblue"),
                                                  angle = 45,
                                                  dpi = 600,
                                                  width = 15, 
                                                  height = 15,
                                                  save.to.tiff = FALSE,
                                                  genorm.score = FALSE,
                                                  norm.finder.score = TRUE)








#### BESTKEEPER AND NORMFINDER PAIRWISE V ANALYSIS using CQDELTACT ####

#First Bestkeeper

#Preparing dataframe for bestkeeper analysis for all treatment groups
all.samples.bestkeeperdf <- subset(data.Ct.ready, data.Ct.ready$Group=="complete" )
all.samples.bestkeeperdf <- all.samples.bestkeeperdf[,c(2:3, 5:9)]
all.samples.bestkeeperdf <- as.data.frame(all.samples.bestkeeperdf)
all.samples.bestkeeperdf <- all.samples.bestkeeperdf %>% column_to_rownames(var = "Sample")
all.samples.bestkeeper.res <- bestKeeper(all.samples.bestkeeperdf)

#Preparing dataframe for bestkeeper analysis for uninfiltrated plants
uninfiltrated.bestkeeperdf <- subset(data.Ct.ready, data.Ct.ready$Group=="uninfiltrated")
uninfiltrated.bestkeeperdf <- uninfiltrated.bestkeeperdf[,c(2:3, 5:9)]
uninfiltrated.bestkeeperdf <- as.data.frame(uninfiltrated.bestkeeperdf)
uninfiltrated.bestkeeperdf <- uninfiltrated.bestkeeperdf %>% column_to_rownames(var = "Sample")
uninfiltrated.bestkeeperdf.res <- bestKeeper(uninfiltrated.bestkeeperdf)

#Preparing dataframe for bestkeeper analysis for infiltrated plants
vigsinfiltrated.bestkeeperdf <- subset(data.Ct.ready, data.Ct.ready$Group=="vigs_infiltrated")
vigsinfiltrated.bestkeeperdf <- vigsinfiltrated.bestkeeperdf[,c(2:3, 5:9)]
vigsinfiltrated.bestkeeperdf <- as.data.frame(vigsinfiltrated.bestkeeperdf)
vigsinfiltrated.bestkeeperdf <- vigsinfiltrated.bestkeeperdf %>% column_to_rownames(var = "Sample")
vigsinfiltrated.bestkeeperdf.res <- bestKeeper(vigsinfiltrated.bestkeeperdf)

#Preparing dataframe for bestkeeper analysis for uninfested plants
uninfested.bestkeeperdf <- subset(data.Ct.ready, data.Ct.ready$Group=="uninfested")
uninfested.bestkeeperdf <- uninfested.bestkeeperdf[,c(2:3, 5:9)]
uninfested.bestkeeperdf <- as.data.frame(uninfested.bestkeeperdf)
uninfested.bestkeeperdf <- uninfested.bestkeeperdf %>% column_to_rownames(var = "Sample")
uninfested.bestkeeperdf.res <- bestKeeper(uninfested.bestkeeperdf)

#Preparing dataframe for bestkeeper analysis for infested plants
infested.bestkeeperdf <- subset(data.Ct.ready, data.Ct.ready$Group=="infested")
infested.bestkeeperdf <- infested.bestkeeperdf[,c(2:3, 5:9)]
infested.bestkeeperdf <- as.data.frame(infested.bestkeeperdf)
infested.bestkeeperdf <- infested.bestkeeperdf %>% column_to_rownames(var = "Sample")
infested.bestkeeperdf.res <- bestKeeper(infested.bestkeeperdf)

#Preparing dataframe for bestkeeper analysis for infested plants
fourteendpi.bestkeeperdf <- subset(data.Ct.ready, data.Ct.ready$Group=="fourteen")
fourteendpi.bestkeeperdf <- fourteendpi.bestkeeperdf[,c(2:3, 5:9)]
fourteendpi.bestkeeperdf <- as.data.frame(fourteendpi.bestkeeperdf)
fourteendpi.bestkeeperdf <- fourteendpi.bestkeeperdf %>% column_to_rownames(var = "Sample")
fourteendpi.bestkeeperdf.res <- bestKeeper(fourteendpi.bestkeeperdf, ctVal = TRUE)

#Preparing dataframe for bestkeeper analysis for infested plants
twentyonedpi.bestkeeperdf <- subset(data.Ct.ready, data.Ct.ready$Group=="twenty-one")
twentyonedpi.bestkeeperdf <- twentyonedpi.bestkeeperdf[,c(2:3, 5:9)]
twentyonedpi.bestkeeperdf <- as.data.frame(twentyonedpi.bestkeeperdf)
twentyonedpi.bestkeeperdf <- twentyonedpi.bestkeeperdf %>% column_to_rownames(var = "Sample")
twentyonedpi.bestkeeperdf.res <- bestKeeper(twentyonedpi.bestkeeperdf)

#VariationV function uses same df as above

#VariationV of all samples
pv.all <- pairwiseV(all.samples.bestkeeperdf)

#VariationV of uninfiltrated samples
pv.uninfiltrated <- pairwiseV(uninfiltrated.bestkeeperdf)

#VariationV of vigs-infiltrated plants
pv.vigs <- pairwiseV(vigsinfiltrated.bestkeeperdf)

#VariationV of uninfested plants
pv.uninfesetd <- pairwiseV(uninfested.bestkeeperdf)

#VariationV of infested plants
pv.infested <- pairwiseV(infested.bestkeeperdf)


#### RankAgg ####

#Raw DF for Rank aggregation
rankagg.df <- read.csv("tot_gene_rankings.csv", header = TRUE)


## Subset for all samples
rankagg.df.all <- subset(rankagg.df, rankagg.df$condition=="all")

# BruteForce algorithm
bruteagg.df.all <- BruteAggreg(as.matrix(rankagg.df.all[3:8]), k=6, distance="Spearman")

# List for al raggr dfs
raggrlist <- list()

# Functon for pulling out rankings for ggplot
# This function was temporarily modified for 21 DPI only to reflect different order of data, data1, and BF to BF, data, data1
extract_rankings <- function(raggr_object){
  
  # get each set of rankings
  r1 <- as.character(raggr_object$lists[1,])
  r2 <- as.character(raggr_object$lists[2,])
  r3 <- as.character(raggr_object$lists[3,])
  
  #Combine into vector
  gene <- c(r1,r2,r3)
  
  # Create the rank column
  rank <- rep(1:6, 3)
  
  # Create the group column
  group <- c(rep("data", 6), rep("data1", 6), rep("data1", 6))
  
  # Combine into the final data frame
  final_df <- data.frame(gene = gene, rank = rank, group = group)
  
  mean_ranks <- final_df %>%
    filter(group %in% c("data", "data1")) %>%
    group_by(gene) %>%
    summarise(mean_rank = mean(rank, na.rm = TRUE)) %>%
    mutate(group = "mean") %>%
    rename(rank = mean_rank)
  
  final_df <- rbind(final_df, mean_ranks)
  
  raggrlist <<- append(raggrlist, list(final_df))
  
  # Display the resulting data frame
  return()
  }


# Subset for uninfiltrated samples
rankagg.df.uninfiltrated <- subset(rankagg.df, rankagg.df$condition=="uninfiltrated")

# BruteForce algorithm
bruteagg.df.uninfiltrated <- BruteAggreg(as.matrix(rankagg.df.uninfiltrated[3:8]), k=6, distance="Spearman")

# Subset for vigs-infiltrated samples

#Subset for all samples
rankagg.df.vigsinfiltrated <- subset(rankagg.df, rankagg.df$condition=="vigs-infiltrated")

#BruteForce algorithm
bruteagg.df.vigsinfiltrated <- BruteAggreg(as.matrix(rankagg.df.vigsinfiltrated[3:8]), k=6, distance="Spearman")

## Subset for uninfested samples

#Subset for all samples
rankagg.df.uninfested <- subset(rankagg.df, rankagg.df$condition=="uninfested")

#BruteForce algorithm
bruteagg.df.uninfested <- BruteAggreg(as.matrix(rankagg.df.uninfested[3:8]), k=6, distance="Spearman")

## Subset for infested samples

#Subset for all samples
rankagg.df.infested <- subset(rankagg.df, rankagg.df$condition=="infested")

#BruteForce algorithm
bruteagg.df.infested <- BruteAggreg(as.matrix(rankagg.df.infested[3:8]), k=6, distance="Spearman")

## Subset for 14-DPI samples

#Subset for all samples
rankagg.df.14dpi <- subset(rankagg.df, rankagg.df$condition=="14dpi")

#BruteForce algorithm
bruteagg.df.14dpi <- BruteAggreg(as.matrix(rankagg.df.14dpi[3:8]), k=6, distance="Kendall")

#Subset for all samples
rankagg.df.21dpi <- subset(rankagg.df, rankagg.df$condition=="21dpi")

#BruteForce algorithm
bruteagg.df.21dpi <- BruteAggreg(as.matrix(rankagg.df.21dpi[3:8]), k=6, distance="Spearman")


#
extract_rankings(bruteagg.df.all) #1
extract_rankings(bruteagg.df.uninfiltrated) #2
extract_rankings(bruteagg.df.vigsinfiltrated) #3
extract_rankings(bruteagg.df.uninfested) #4
extract_rankings(bruteagg.df.infested) #5
extract_rankings(bruteagg.df.14dpi) #6
extract_rankings(bruteagg.df.21dpi) #8 actually



#Plotting
brutaggr_allsamples <- ggplot(as.data.frame(raggrlist[1]), aes(x = factor(gene, levels = c("GhUBQ7", "GhUBQ14", "GhTBL6", "GhACT7", "GhTMN5", "GhPP2A1")), y = rank, group = group)) +
  geom_line(aes(color = group), linewidth = 1) +
  scale_color_manual(values = c("data" = "grey", "data1" = "grey", "BF" = "red", "mean"= "black")) +
  labs(y = "Rank", x = "", color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "italic"),
        legend.position = "top",
        panel.grid.minor = element_blank()) +
  geom_text(data = tibble(x = 6, y = 2), aes(x = x, y = y, label = "All samples", fontface = "bold"), size = 6, inherit.aes = FALSE)

export_figure(brutaggr_allsamples, "bruteagg.all.conditions.plot.jpeg", 24, 6)

#Plotting
brutaggr_uninfiltrated <- ggplot(as.data.frame(raggrlist[2]), aes(x = factor(gene, levels = c("GhUBQ7", "GhUBQ14", "GhACT7", "GhTBL6", "GhPP2A1", "GhTMN5")), y = rank, group = group)) +
  geom_line(aes(color = group), linewidth = 1) +
  scale_color_manual(values = c("data" = "grey", "data1" = "grey", "BF" = "red", "mean"= "black")) +
  labs(y = "Rank", x = "", color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "italic"),
        legend.position = "top",
        panel.grid.minor = element_blank()) +
  geom_text(data = tibble(x = 6, y = 2), aes(x = x, y = y, label = "Uninfiltrated", fontface = "bold"), size = 6, inherit.aes = FALSE)

export_figure(brutaggr_uninfiltrated, "bruteagg.uninfiltrated.plot.jpeg", 24, 6)

#Plotting
brutaggr_infiltrated <- ggplot(as.data.frame(raggrlist[3]), aes(x = factor(gene, levels = c("GhUBQ7", "GhUBQ14", "GhTBL6", "GhACT7","GhTMN5", "GhPP2A1")), y = rank, group = group)) +
  geom_line(aes(color = group), linewidth = 1) +
  scale_color_manual(values = c("data" = "grey", "data1" = "grey", "BF" = "red", "mean"= "black")) +
  labs(y = "Rank", x = "", color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "italic"),
        legend.position = "top",
        panel.grid.minor = element_blank()) +
  geom_text(data = tibble(x = 6, y = 2), aes(x = x, y = y, label = "VIGS", fontface = "bold"), size = 6, inherit.aes = FALSE)

export_figure(brutaggr_infiltrated, "bruteagg.infiltrated.plot.jpeg", 24, 6)

#Plotting
brutaggr_uninfested <- ggplot(as.data.frame(raggrlist[4]), aes(x = factor(gene, levels = c("GhUBQ7", "GhUBQ14", "GhACT7", "GhTBL6","GhPP2A1", "GhTMN5")), y = rank, group = group)) +
  geom_line(aes(color = group), linewidth = 1) +
  scale_color_manual(values = c("data" = "grey", "data1" = "grey", "BF" = "red", "mean"= "black")) +
  labs(y = "Rank", x = "", color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "italic"),
        legend.position = "top",
        panel.grid.minor = element_blank()) +
  geom_text(data = tibble(x = 6, y = 2), aes(x = x, y = y, label = "Uninfested", fontface = "bold"), size = 6, inherit.aes = FALSE)

export_figure(brutaggr_uninfested, "bruteagg.uninfested.plot.jpeg", 24, 6)

#Plotting
brutaggr_infested <- ggplot(as.data.frame(raggrlist[5]), aes(x = factor(gene, levels = c("GhUBQ7", "GhUBQ14", "GhTBL6", "GhACT7", "GhTMN5", "GhPP2A1")), y = rank, group = group)) +
  geom_line(aes(color = group), linewidth = 1) +
  scale_color_manual(values = c("data" = "grey", "data1" = "grey", "BF" = "red", "mean"= "black")) +
  labs(y = "Rank", x = "", color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "italic"),
        legend.position = "top",
        panel.grid.minor = element_blank()) +
  geom_text(data = tibble(x = 6, y = 2), aes(x = x, y = y, label = "Infested", fontface = "bold"), size = 6, inherit.aes = FALSE)

export_figure(brutaggr_infested, "bruteagg.infested.plot.jpeg", 24, 6)

#Plotting
brutaggr_14dpi <- ggplot(as.data.frame(raggrlist[6])[7:24,], aes(x = factor(gene, levels = c("GhUBQ7", "GhUBQ14", "GhTBL6", "GhTMN5", "GhACT7", "GhPP2A1")), y = rank, group = group)) +
  geom_line(aes(color = group), linewidth = 1) +
  scale_color_manual(values = c("data1" = "grey", "BF" = "red", "mean"= "black")) +
  labs(y = "Rank", x = "", color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "italic"),
        legend.position = "top",
        panel.grid.minor = element_blank()) +
  geom_text(data = tibble(x = 6, y = 2), aes(x = x, y = y, label = "14 DPI", fontface = "bold"), size = 6, inherit.aes = FALSE)

export_figure(brutaggr_14dpi, "bruteagg.df.14dpi.jpeg", 24, 6)

#Plotting
brutaggr_21dpi <- ggplot(as.data.frame(raggrlist[8]), aes(x = factor(gene, levels = c("GhUBQ7", "GhUBQ14", "GhTBL6", "GhTMN5", "GhACT7", "GhPP2A1")), y = rank, group = group)) +
  geom_line(aes(color = group), linewidth = 1) +
  scale_color_manual(values = c("data" = "grey", "data1" = "grey", "BF" = "red", "mean"= "black")) +
  labs(y = "Rank", x = "", color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "italic"),
        legend.position = "top",
        panel.grid.minor = element_blank()) +
  geom_text(data = tibble(x = 6, y = 2), aes(x = x, y = y, label = "21 DPI", fontface = "bold"), size = 6, inherit.aes = FALSE)

export_figure(brutaggr_21dpi, "bruteagg.df.21dpi.jpeg", 24, 6)

### BESTKEEPER FIGURES ####

#DF
bestkeeper.df <- read.csv("bestkeeper.csv", header = TRUE)
bestkeeper.df$gene <- dplyr::recode(bestkeeper.df$gene, "ACT7" = "GhACT7",
                                        "UBQ7" = "GhUBQ7",
                                        "UBQ14" = "GhUBQ14",
                                        "PP2A1" = "GhPP2A1",
                                        "TMN5" = "GhTMN5",
                                        "TBL6" = "GhTBL6")


#All conditions
bestkeeper.df.all <- subset(bestkeeper.df, bestkeeper.df$contrast=="all")
bestkeeper.df.all$gene <- factor(bestkeeper.df.all$gene, levels = c("GhUBQ7", "GhUBQ14", "GhACT7", "GhPP2A1", "GhTMN5", "GhTBL6"))
 
bestkeeper.df.all.plot <- ggplot(bestkeeper.df.all, aes(x=gene)) +
                          geom_bar(aes(y=sd.value), stat = "identity", color="darkgrey", alpha = 0.7) +
                          geom_line(aes(y=cv.value, group = condition), size = 1.0) +
                          geom_point(aes(y=cv.value), size = 2) +
                          xlab("Candidate reference gene") + 
                          scale_y_continuous(name = "CV value", sec.axis = sec_axis(transform=~.*1, name="SD value")) +
                          theme_classic() +
                          theme(axis.text=element_text(size=12, color="black"),
                                axis.text.x=element_text(angle=45,vjust = 1, hjust=1, face="italic"),
                                axis.line=element_line(size=1.0),
                                axis.ticks=element_line(size=1.0),
                                axis.ticks.length=unit(.25, "cm"),
                                legend.position = "top",
                                legend.text=element_text(size=rel(2.5)),
                                legend.title=element_text(size=rel(3)),
                                axis.title=element_text(size=12,face="bold")) +
                          geom_text(data=tibble(x=4.0, y=9.9), aes(x=x,y=y, label= "All conditions", fontface="bold"), size = 8, inherit.aes=FALSE)

#VIGS status
bestkeeper.df.vigs_status <- subset(bestkeeper.df, bestkeeper.df$contrast=="vigs_status")


#Uninfiltrated
bestkeeper.df.vigs_status.uninfiltrated.plot <- ggplot(bestkeeper.df.vigs_status[1:6,], aes(x=factor(bestkeeper.df.vigs_status[1:6,]$gene, levels = c("GhUBQ7", "GhUBQ14", "GhACT7", "GhPP2A1", "GhTMN5", "GhTBL6")))) + # nolint: line_length_linter.
  geom_bar(aes(y=sd.value), stat = "identity", color="darkgrey", alpha = 0.7) +
  geom_line(aes(y=cv.value, group = condition), size = 1.0) +
  geom_point(aes(y=cv.value), size = 2) +
  xlab("Candidate reference gene") + 
  scale_y_continuous(name = "CV value", sec.axis = sec_axis(transform=~.*1, name="SD value")) +
  theme_classic() +
  theme(axis.text=element_text(size=12, color="black"),
        axis.text.x=element_text(angle=45,vjust = 1, hjust=1, face = "italic"),
        axis.line=element_line(size=1.0),
        axis.ticks=element_line(size=1.0),
        axis.ticks.length=unit(.25, "cm"),
        legend.position = "top",
        legend.text=element_text(size=rel(2.5)),
        legend.title=element_text(size=rel(3)),
        axis.title=element_text(size=12,face="bold")) +
  geom_text(data=tibble(x=4.4, y=9.9), aes(x=x,y=y, label= "Uninfiltrated", fontface="bold"), size = 8, inherit.aes=FALSE)

#Infiltrated
bestkeeper.df.vigs_status.infiltrated.plot <- ggplot(bestkeeper.df.vigs_status[7:12,], aes(x=factor(bestkeeper.df.vigs_status[7:12,]$gene, levels = c("GhUBQ7", "GhUBQ14", "GhACT7", "GhPP2A1", "GhTMN5", "GhTBL6")))) +
  geom_bar(aes(y=sd.value), stat = "identity", color="darkgrey", alpha = 0.7) +
  geom_line(aes(y=cv.value, group = condition), size = 1.0) +
  geom_point(aes(y=cv.value), size = 2) +
  xlab("Candidate reference gene") + 
  scale_y_continuous(name = "CV value", sec.axis = sec_axis(transform=~.*1, name="SD value")) +
  theme_classic() +
  theme(axis.text=element_text(size=12, color="black"),
        axis.text.x=element_text(angle=45,vjust = 1, hjust=1, face="italic"),
        axis.line=element_line(size=1.0),
        axis.ticks=element_line(size=1.0),
        axis.ticks.length=unit(.25, "cm"),
        legend.position = "top",
        legend.text=element_text(size=rel(2.5)),
        legend.title=element_text(size=rel(3)),
        axis.title=element_text(size=12,face="bold")) +
  geom_text(data=tibble(x=5.5, y=9.9), aes(x=x,y=y, label= "VIGS", fontface="bold"), size = 8, inherit.aes=FALSE)


# Infestation status
bestkeeper.df.infestation_status <- subset(bestkeeper.df, bestkeeper.df$contrast=="infestation_status")

# Uninfested
bestkeeper.df.infestation_status.uninfested.plot <- ggplot(bestkeeper.df.infestation_status[1:6,], aes(x=factor(bestkeeper.df.infestation_status[1:6,]$gene, levels = c("GhUBQ7", "GhACT7", "GhUBQ14", "GhPP2A1", "GhTMN5", "GhTBL6")))) +
                                                    geom_bar(aes(y=sd.value), stat = "identity", color="darkgrey", alpha = 0.7) +
                                                    geom_line(aes(y=cv.value, group = condition), size = 1.0) +
                                                    geom_point(aes(y=cv.value), size = 2) +
                                                    xlab("Candidate reference gene") + 
                                                    scale_y_continuous(name = "CV value", sec.axis = sec_axis(transform=~.*1, name="SD value")) +
                                                    theme_classic() +
                                                    theme(axis.text=element_text(size=12, color="black"),
                                                            axis.text.x=element_text(angle=45, vjust = 1, hjust=1, face = "italic"),
                                                            axis.line=element_line(size=1.0),
                                                            axis.ticks=element_line(size=1.0),
                                                            axis.ticks.length=unit(.25, "cm"),
                                                            legend.position = "top",
                                                            legend.text=element_text(size=rel(2.5)),
                                                            legend.title=element_text(size=rel(3)),
                                                            axis.title=element_text(size=12,face="bold")) +
                                                    geom_text(data=tibble(x=4.5, y=9.9), aes(x=x,y=y, label= "Uninfested", fontface="bold"), size = 8, inherit.aes=FALSE)


#Infested
bestkeeper.df.infestation_status.infested.plot <- ggplot(bestkeeper.df.infestation_status[7:12,], aes(x=factor(bestkeeper.df.infestation_status[7:12,]$gene, levels = c("GhUBQ7", "GhUBQ14", "GhACT7", "GhPP2A1", "GhTMN5", "GhTBL6")))) +
  geom_bar(aes(y=sd.value), stat = "identity", color="darkgrey", alpha = 0.7) +
  geom_line(aes(y=cv.value, group = condition), size = 1.0) +
  geom_point(aes(y=cv.value), size = 2) +
  xlab("Candidate reference gene") + 
  scale_y_continuous(name = "CV value", sec.axis = sec_axis(transform=~.*1, name="SD value")) +
  theme_classic() +
  scale_x_discrete(labels= c("GhACT7" = expression(italic("GhACT7")),
                             "GhUBQ7" = expression(italic("GhUBQ7")),
                             "GhUBQ14" = expression(italic("GhUBQ14")),
                             "GhPP2A1" = expression(italic("GhPP2A1")),
                             "GhTMN5" = expression(italic("GhTMN5")),
                             "GhTBL6" = expression(italic("GhTBL6")))) +
  theme(axis.text=element_text(size=12, color="black"),
        axis.text.x=element_text(angle=45, vjust = 1, hjust=1),
        axis.line=element_line(size=1.0),
        axis.ticks=element_line(size=1.0),
        axis.ticks.length=unit(.25, "cm"),
        legend.position = "top",
        legend.text=element_text(size=rel(2.5)),
        legend.title=element_text(size=rel(3)),
        axis.title=element_text(size=12,face="bold")) +
  geom_text(data=tibble(x=5.0, y=9.9), aes(x=x,y=y, label= "Infested", fontface="bold"), size = 8, inherit.aes=FALSE)


# By time
bestkeeper.df.bytime <- subset(bestkeeper.df, bestkeeper.df$contrast=="bytime")

# 14 DPI
bestkeeper.df.bytime.fourteen.plot <- ggplot(bestkeeper.df.bytime[7:12,], aes(x=factor(bestkeeper.df.bytime[7:12,]$gene, levels=c("GhUBQ7", "GhUBQ14", "GhACT7", "GhPP2A1", "GhTMN5","GhTBL6")))) +
  geom_bar(aes(y=sd.value), stat = "identity", color="darkgrey", alpha = 0.7) +
  geom_line(aes(y=cv.value, group = condition), size = 1.0) +
  geom_point(aes(y=cv.value), size = 2) +
  xlab("Candidate reference gene") + 
  scale_y_continuous(name = "CV value", sec.axis = sec_axis(transform=~.*1, name="SD value")) +
  theme_classic() +
  theme(axis.text=element_text(size=12, color="black"),
        axis.text.x=element_text(angle=45,vjust = 1, hjust=1, face="italic"),
        axis.line=element_line(size=1.0),
        axis.ticks=element_line(size=1.0),
        axis.ticks.length=unit(.25, "cm"),
        legend.position = "top",
        legend.text=element_text(size=rel(2.5)),
        legend.title=element_text(size=rel(3)),
        axis.title=element_text(size=12,face="bold")) +
  geom_text(data=tibble(x=5.2, y=9.9), aes(x=x,y=y, label= "14 DPI", fontface="bold"), size = 8, inherit.aes=FALSE)

# 21 DPI
bestkeeper.df.bytime.twentyone.plot <- ggplot(bestkeeper.df.bytime[1:6,], aes(x=factor(bestkeeper.df.bytime[1:6,]$gene, levels=c("GhUBQ7", "GhUBQ14", "GhACT7", "GhPP2A1", "GhTMN5","GhTBL6")))) +
  geom_bar(aes(y=sd.value), stat = "identity", color="darkgrey", alpha = 0.7) +
  geom_line(aes(y=cv.value, group = condition), size = 1.0) +
  geom_point(aes(y=cv.value), size = 2) +
  xlab("Candidate reference gene") + 
  scale_y_continuous(name = "CV value", sec.axis = sec_axis(transform=~.*1, name="SD value")) +
  theme_classic() +
  theme(axis.text=element_text(size=12, color="black"),
        axis.text.x=element_text(angle=45,vjust = 1, hjust=1, face="italic"),
        axis.line=element_line(size=1.0),
        axis.ticks=element_line(size=1.0),
        axis.ticks.length=unit(.25, "cm"),
        legend.position = "top",
        legend.text=element_text(size=rel(2.5)),
        legend.title=element_text(size=rel(3)),
        axis.title=element_text(size=12,face="bold")) +
  geom_text(data=tibble(x=5.2, y=11.9), aes(x=x,y=y, label= "21 DPI", fontface="bold"), size = 8, inherit.aes=FALSE)


#Now saving plots
#All
export_figure(bestkeeper.df.all.plot, "bestkeeper.all.plot.jpeg", 12,1)
#Uninfiltrated
export_figure(bestkeeper.df.vigs_status.uninfiltrated.plot, "bestkeeper.uninfiltrated.plot.jpeg", 12,16)
#VIGSed
export_figure(bestkeeper.df.vigs_status.infiltrated.plot, "bestkeeper.VIGSed.plot.jpeg", 12,16)
#Uninfested
export_figure(bestkeeper.df.infestation_status.uninfested.plot, "bestkeeper.uninfested.plot.jpeg", 12,16)
#Infested
export_figure(bestkeeper.df.infestation_status.infested.plot, "bestkeeper.infested.plot.jpeg", 12,16)
#14 DPI
export_figure(bestkeeper.df.bytime.fourteen.plot, "bestkeeper.14dpi.plot.jpeg", 12,16)
#21 DPI
export_figure(bestkeeper.df.bytime.twentyone.plot, "bestkeeper.21dpi.plot.jpeg", 12,16)


### OTHER SUMMARY STATISTICS FOR MANUSCRIPT #####

#This includes HYD1 KD samples
#Want to exclude those samples for now for HYD1 Ct value average
summstats_bycondition <- data.Ct %>% 
                         group_by(Gene, Group) %>% 
                         summarise(mean = mean(Ct),
                                   sd=sd(Ct))

summstats_whole_without_HYD1KD <- subset(rawcts, rawcts$gene_target != "HYD1")  %>% 
  group_by(target) %>% 
  summarise(mean = mean(Cq),
            sd=sd(Cq),
            min=min(Cq),
            max=max(Cq))


### HYD1 RELATIVE EXPRESSION ANALYSIS - NORMALIZED TO GhACT7 AND GhPP2A1 ####

hyd1.expression <- read.csv("hyd1_normalization_GhACT7GhPP2A1.csv", header = TRUE)

#Subset by time: to be analyzed separately
hyd1.expression.14dpi <- subset(hyd1.expression,hyd1.expression$days.post.infestation=="fourteen")
hyd1.expression.21dpi <- subset(hyd1.expression,hyd1.expression$days.post.infestation=="twenty-one")

#Plot 14 dpi
hyd1.expression.14dpi.plot <- ggplot(hyd1.expression.14dpi, aes(x=infiltration.status, y = relative.expression)) +
  geom_boxplot(lwd=1.2) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Normalized GhHYDRA1 Expression") +
  theme(axis.text=element_text(size=26, face = "bold", color="black"),
        plot.margin = unit(c(1.0,0.5,0.5,0.5), "cm"),
        axis.line=element_line(size=1.5),
        axis.ticks=element_line(size=1.5),
        axis.ticks.length=unit(.25, "cm"),
        axis.title=element_text(size=28,face="bold"))

#Relevel infestation status
#Changing level order
hyd1.expression.14dpi$infestation.status <- factor(hyd1.expression.14dpi$infestation.status, levels = c("uninfested","infested"))

ylab <- expression(bold(paste("Relative ", bolditalic("GhHYDRA1 "), "expression")))

hyd1.expression.14dpi.barplot <- ggplot(hyd1.expression.14dpi, aes(x=infestation.status, y = relative.expression, fill = infiltration.status)) +
  geom_bar(stat="summary", position = position_dodge(), color="black", alpha=0.7, size=1.0) +
  geom_errorbar(stat="summary", position = position_dodge(0.9), width=0.4, size=1.0, color = "black") +
  scale_fill_manual(name="Infiltration status", labels=c("uninfiltrated", "VIGS"), values=c("grey36", "grey76")) +
  scale_y_continuous(limits = c(0,2.8), breaks = c(0,0.5,1.0,1.5,2.0,2.5)) +
  theme_classic() +
  xlab("Infestation status") +
  ylab(ylab) +
  theme(axis.text=element_text(size=12, color="black"),
        plot.margin = unit(c(1.0,0.5,0.5,0.5), "cm"),
        axis.line=element_line(size=1.0),
        axis.ticks=element_line(size=1.0),
        axis.ticks.length=unit(.25, "cm"),
        axis.title=element_text(size=12, face="bold")) +
  geom_line(data=tibble(x=c(0.75,1.25), y=c(0.9,0.9)), aes(x=x,y=y), size = 1.0, inherit.aes = FALSE) +
  geom_text(data=tibble(x=1.0, y=1.05), aes(x=x,y=y, label= "n.s.", fontface="bold"), size = 8, inherit.aes=FALSE) +
  geom_line(data=tibble(x=c(1.75,2.25), y=c(2.1,2.1)), aes(x=x,y=y), size = 1.0, inherit.aes = FALSE) +
  geom_text(data=tibble(x=2.0, y=2.25), aes(x=x,y=y, label= "n.s.", fontface="bold"), size = 8, inherit.aes=FALSE) +
  geom_line(data=tibble(x=c(1.0,2.0), y=c(2.4,2.4)), aes(x=x,y=y), size = 1.0, inherit.aes = FALSE) +
  geom_text(data=tibble(x=1.5, y=2.55), aes(x=x,y=y, label= "p = 0.0001", fontface="bold"), size = 8, inherit.aes=FALSE) +
  geom_text(data=tibble(x=0.70, y=2.7), aes(x=x,y=y, label= "14 DPI", fontface="bold"), size = 8, inherit.aes=FALSE)

# Now analysis of effect experimental conditions had relative expression
relative.expression.model <- lmer(relative.expression ~ infestation.status * days.post.infestation + infiltration.status + (1|plant.ID), dat=hyd1.expression)
anova(relative.expression.model)

# Reduced version of model
relative.expression.model_reduced <- lm(relative.expression ~ infestation.status * days.post.infestation + infiltration.status, dat=hyd1.expression)

anova(relative.expression.model, relative.expression.model_reduced)

# Checking normality of residuals
shapiro.test(residuals(relative.expression.model)) # p >0 .05
qqnorm(residuals(relative.expression.model))
qqline(residuals(relative.expression.model))

# Get interactions among conditions
relative.expression.model.interaction <- interaction(hyd1.expression$infestation.status, hyd1.expression$infiltration.status, hyd1.expression$days.post.infestation)

# Levene's test for homoscedasticity of residuals
leveneTest(residuals(relative.expression.model) ~ relative.expression.model.interaction, center = "median") # p > 0.05

# Shows significant effect of infestation status
# Now use eemeans to get estimated marginal means for the interaction
relative.expression.model.emm <- emmeans(relative.expression.model, ~ infestation.status | days.post.infestation)
pairs(relative.expression.model.emm, adjust = "bonferroni")
summary(pairs(relative.expression.model.emm, adjust = "bonferroni", infer = TRUE))

#Now run a cohen's d to standardize effect sizes using effectsize package
# Define the inputs

#14 DPI emmans pair results
estimate_14 <- 0.900
SE_14 <- 0.196
df_14 <- 30.7

#21 DPI emmans pair results
estimate_21 <- 0.451
SE_21 <- 0.233
df_21 <- 31.2

# Back-calculate pooled SD
pooled_sd_14 <- SE_14 * sqrt(df_14 + 1)
pooled_sd_21 <- SE_21 * sqrt(df_21 + 1)

# Calculate Cohen's d
cohens_d_14 <- estimate_14 / pooled_sd_14
cohens_d_21 <- estimate_21 / pooled_sd_21

# Print the results
cat("Cohen's d for Day 14:", cohens_d_14, "\n")
cat("Cohen's d for Day 21:", cohens_d_21, "\n")    

#Onto 21days DPI
#Relevel infestation status
#Changing level order
hyd1.expression.21dpi$infestation.status <- factor(hyd1.expression.21dpi$infestation.status, levels = c("uninfested","infested"))


hyd1.expression.21dpi.barplot <- ggplot(hyd1.expression.21dpi, aes(x=infestation.status, y = relative.expression, fill = infiltration.status)) +
  geom_bar(stat="summary", position = position_dodge(), color="black", alpha=0.7, size=1.0) +
  geom_errorbar(stat="summary", position = position_dodge(0.9), width=0.4, size=1.0, color = "black") +
  scale_fill_manual(name="Infiltration status", labels=c("uninfiltrated", "VIGS"), values=c("grey36", "grey76")) +
  theme_classic() +
  xlab("Infestation status") +
  ylab(ylab) +
  theme(axis.text=element_text(size=12, color="black"),
        plot.margin = unit(c(1.0,0.5,0.5,0.5), "cm"),
        axis.line=element_line(size=1.0),
        axis.ticks=element_line(size=1.0),
        axis.ticks.length=unit(.25, "cm"),
        axis.title=element_text(size=12,face="bold")) +
  geom_line(data=tibble(x=c(0.75,1.25), y=c(1.2,1.2)), aes(x=x,y=y), size = 1.0, inherit.aes = FALSE) +
  geom_text(data=tibble(x=1.0, y=1.3), aes(x=x,y=y, label= "n.s.", fontface="bold"), size = 8, inherit.aes=FALSE) +
  geom_line(data=tibble(x=c(1.75,2.25), y=c(1.6,1.6)), aes(x=x,y=y), size = 1.0, inherit.aes = FALSE) +
  geom_text(data=tibble(x=2.0, y=1.7), aes(x=x,y=y, label= "n.s.", fontface="bold"), size = 8, inherit.aes=FALSE) +
  geom_line(data=tibble(x=c(1.0,2.0), y=c(1.9,1.9)), aes(x=x,y=y), size = 1.0, inherit.aes = FALSE) +
  geom_text(data=tibble(x=1.5, y=2.0), aes(x=x,y=y, label= "p = 0.06", fontface="bold"), size = 8, inherit.aes=FALSE) +
  geom_text(data=tibble(x=0.70, y=2.1), aes(x=x,y=y, label= "21 DPI", fontface="bold"), size = 8, inherit.aes=FALSE)

export_figure(hyd1.expression.14dpi.barplot, "hyd1.expression.14dpi.barplot.jpeg", 18,20)
export_figure(hyd1.expression.21dpi.barplot, "hyd1.expression.21dpi.barplot.jpeg", 18,20)

#Sumstats of hyd1 expression
hyd1.expression.summstats <- hyd1.expression %>%
                              group_by(infestation.status, infiltration.status, days.post.infestation) %>%
                              summarise(mean = mean(relative.expression),
                                        sd = sd(relative.expression),
                                        n = n(),
                                        se = sd/sqrt(n))

#Sumstats of hyd1 expression
hyd1.expression.summstats1 <- hyd1.expression %>%
  group_by(infestation.status, days.post.infestation) %>%
  summarise(mean = mean(relative.expression),
            sd = sd(relative.expression),
            n = n(),
            se = sd/sqrt(n))

### HYD1 RELATIVE EXPRESSION ANALYSIS - NORMALIZED TO UBQ7 ####

#Df of HYD1 normalized to Ubq7
hyd1_normed_to_ubq7 <- read.csv("hyd1_normalized_to_ubq7.csv", header = TRUE)


#Sumstats of hyd1 expression
hyd1.expression.norm_Ubq7_sumstats <- hyd1_normed_to_ubq7 %>%
  group_by(infestation.status, days.post.infestation) %>%
  summarise(mean = mean(relative.expression),
            sd = sd(relative.expression),
            n = n(),
            se = sd/sqrt(n))

#Subset by time
hyd1_normed_to_ubq7.14dpi <- subset(hyd1_normed_to_ubq7, hyd1_normed_to_ubq7$days.post.infestation=="fourteen")
hyd1_normed_to_ubq7.21dpi <- subset(hyd1_normed_to_ubq7, hyd1_normed_to_ubq7$days.post.infestation=="twenty-one")

hyd1_normed_to_ubq7.14dpi$infestation.status <- factor(hyd1_normed_to_ubq7.14dpi$infestation.status, levels = c("uninfested","infested"))
hyd1_normed_to_ubq7.21dpi$infestation.status <- factor(hyd1_normed_to_ubq7.21dpi$infestation.status, levels = c("uninfested","infested"))


#Plot 14 DPI
hyd1_normed_to_ubq7.14dpi.plot <- ggplot(hyd1_normed_to_ubq7.14dpi, aes(x=infestation.status, y = relative.expression, fill = infiltration.status)) +
  geom_bar(stat="summary", position = position_dodge(), color="black", alpha=0.7, size=1.0) +
  geom_errorbar(stat="summary", position = position_dodge(0.9), width=0.4, size=1.0, color = "black") +
  scale_fill_manual(name="Infiltration status", labels=c("uninfiltrated", "VIGS"), values=c("grey36", "grey76")) +
  theme_classic() +
  xlab("Infestation status") +
  ylab(ylab) +
  theme(axis.text=element_text(size=12, color="black"),
        plot.margin = unit(c(1.0,0.5,0.5,0.5), "cm"),
        axis.line=element_line(size=1.0),
        axis.ticks=element_line(size=1.0),
        axis.ticks.length=unit(.25, "cm"),
        axis.title=element_text(size=12,face="bold")) +
  geom_line(data=tibble(x=c(1,2), y=c(12.5,12.5)), aes(x=x,y=y), size = 1.0, inherit.aes = FALSE) +
  geom_text(data=tibble(x=1.5, y=13.2), aes(x=x,y=y, label= "p = 0.06", fontface="bold"), size = 8, inherit.aes=FALSE) +
  geom_text(data=tibble(x=0.70, y=14), aes(x=x,y=y, label= "14 DPI", fontface="bold"), size = 8, inherit.aes=FALSE)

hyd1_normed_to_ubq7.21dpi.plot <- ggplot(hyd1_normed_to_ubq7.21dpi, aes(x=infestation.status, y = relative.expression, fill = infiltration.status)) +
  geom_bar(stat="summary", position = position_dodge(), color="black", alpha=0.7, size=1.0) +
  geom_errorbar(stat="summary", position = position_dodge(0.9), width=0.4, size=1.0, color = "black") +
  scale_fill_manual(name="Infiltration status", labels=c("uninfiltrated", "VIGS"), values=c("grey36", "grey76")) +
  theme_classic() +
  xlab("Infestation status") +
  ylab(ylab) +
  theme(axis.text=element_text(size=12,color="black"),
        plot.margin = unit(c(1.0,0.5,0.5,0.5), "cm"),
        axis.line=element_line(size=1.0),
        axis.ticks=element_line(size=1.0),
        axis.ticks.length=unit(.25, "cm"),
        axis.title=element_text(size=12,face="bold")) +
  geom_line(data=tibble(x=c(1,2), y=c(8)), aes(x=x,y=y), size = 1.0, inherit.aes = FALSE) +
  geom_text(data=tibble(x=1.5, y=8.5), aes(x=x,y=y, label= "p = 0.26", fontface="bold"), size = 8, inherit.aes=FALSE) +
  geom_text(data=tibble(x=0.80, y=10), aes(x=x,y=y, label= "21 DPI", fontface="bold"), size = 8, inherit.aes=FALSE)


# Now analysis of effect experimental conditions had relative expression
relative.expression.model1 <- lmer(relative.expression ~ infestation.status * days.post.infestation + infiltration.status + (1|plant.ID), dat=hyd1_normed_to_ubq7)
anova(relative.expression.model1)

relative.expression.model1_reduced <- lm(relative.expression ~ infestation.status * days.post.infestation + infiltration.status, dat=hyd1_normed_to_ubq7)

anova(relative.expression.model1, relative.expression.model1_reduced)


# Check normality of residuals
shapiro.test(residuals(relative.expression.model1)) # p < 0.05
qqnorm(residuals(relative.expression.model1))
qqline(residuals(relative.expression.model1))

# Get interactions among conditions
relative.expression.model1.interaction <- interaction(hyd1_normed_to_ubq7$infestation.status, hyd1_normed_to_ubq7$infiltration.status, hyd1_normed_to_ubq7$days.post.infestation)

# Levene's test for homoscedasticity of residuals
leveneTest(residuals(relative.expression.model1) ~ relative.expression.model1.interaction, center = "median") # p > 0.05


# Shows significant effect of infestation status
# Now use eemeans to get estimated marginal means for the interaction
relative.expression_ubq7.model.emm <- emmeans(relative.expression.model1, ~ infestation.status | days.post.infestation)
summary(pairs(relative.expression_ubq7.model.emm, adjust = "bonferroni", infer = TRUE))


#Results
export_figure(hyd1_normed_to_ubq7.14dpi.plot, "hyd1_normed_to_ubq7.14dpi.plot.jpeg", 18,20)
export_figure(hyd1_normed_to_ubq7.21dpi.plot, "hyd1_normed_to_ubq7.21dpi.plot.jpeg", 18,20)

#### HYDRA1 RESAMPLING #####

#14 DPI expression values normalized to GhACT7/GhPP2A1
GOI_best_ref <- hyd1.expression.14dpi$relative.expression

#14 DPI expression values normalized using Ubq7
GOI_worst_ref <- hyd1_normed_to_ubq7.14dpi$relative.expression


#Bootstrap function
# Bootstrap function to calculate standard deviation
bootstrap_sd <- function(data, n_bootstrap = 1000) {
  resampled_sd <- replicate(n_bootstrap, {
    sample_data <- sample(data, size = length(data), replace = TRUE)
    sd(sample_data)
  })
  return(resampled_sd)
}

#Run analysis
set.seed(123)  # Set seed for reproducibility
sd_best <- bootstrap_sd(GOI_best_ref)
sd_worst <- bootstrap_sd(GOI_worst_ref)

#combine into df
sd_data <- data.frame(
  sd_value = c(sd_best, sd_worst),
  group = rep(c("GhACT7/GhPP2A1", "Ubq7"), each = length(sd_best))
)

#Plot the data
density_plot <- ggplot(sd_data, aes(x = sd_value, fill = group)) +
                geom_density(alpha = 0.5, size = 1) +
                xlab("Standard deviation") +
                ylab("Density") +
                scale_y_continuous(limits = c(0,4.5), breaks = c(1,2,3,4)) +
                theme_classic() +
                scale_fill_discrete("Normalization gene(s)", labels = c(expression(italic("GhACT7")/italic("GhPP2A1")), expression(italic("GhUBQ7")))) +
                theme(axis.text=element_text(size=12,color="black"),
                      plot.margin = unit(c(1.0,0.5,0.5,0.5), "cm"),
                      axis.line=element_line(size=1.0),
                      axis.ticks=element_line(size=1.0),
                      axis.ticks.length=unit(.25, "cm"),
                      axis.title=element_text(size=12,face="bold"),
                      legend.position="top") +
                geom_line(data=tibble(x=c(0.75,4), y=c(3.6)), aes(x=x,y=y), size = 1.0, inherit.aes = FALSE) +
                geom_text(data=tibble(x=2.3, y=3.75), aes(x=x,y=y, fontface="bold"), label = expression(bold(bolditalic(p) ~ "< 0.0001")), size = 4, inherit.aes=FALSE)

#
t_test_result <- t.test(sd_best, sd_worst)
print(t_test_result)

export_figure(density_plot, "sd_plot_comparing_GhACT7GhPP2A1_ubq7.jpeg", 12,14)

