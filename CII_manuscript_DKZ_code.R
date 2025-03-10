# Stewart, Zachman et al. Complex II regulates purines in AML

# Extended Data Fig 1e
## files needed: Supplemental Table 1, sheet 1 == dataframe
dataframe %>% pivot_longer(cols = 2:3,
             names_to = 'type',
             values_to = 'Canonical Correlation') %>%
  ggplot(mapping = aes(x = complex, y = `Canonical Correlation`, fill = type)) + 
  geom_violin(alpha = 0.5) +
  theme_light() +
  ylim(0,1) +
  scale_x_discrete(labels = label_wrap(10)) +
  xlab("") + 
  theme(axis.title.x = element_text(family = "Arial", size = 16, color = 'black'),
        axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
        axis.title.y = element_text(family = "Arial", size = 16, color = 'black'),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        legend.text = element_text(family = "Arial", size = 16, color = 'black'),
        legend.title = element_blank()) 

# Extended Data Fig 1f
## files needed: Supplemental Table 1, sheet 1 == dataframe
oxphos_pathway_coessentiality %>%
  pivot_longer(cols = 2:3,
               names_to = 'type',
               values_to = 'Canonical Correlation') %>%
  ggplot(mapping = aes(x = `Canonical Correlation`, fill = type)) + 
  geom_density(alpha = 0.5) +
  theme_light() +
  xlab("Canonical Correlation") +
  ylab("Density Estimate") + 
  theme(axis.title.x = element_text(family = "Arial", size = 16, color = 'black'),
        axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
        axis.title.y = element_text(family = "Arial", size = 16, color = 'black'),
        axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
        legend.text = element_text(family = "Arial", size = 16, color = 'black'),
        legend.title = element_blank()) 

# Fig 1b
## files needed: Supplemental Table 1, sheet 1 == dataframe
etc_pathway_matrix <- dataframe %>%
  dplyr::select(!`Gene-Pathway`) %>%
  mutate(complex = gsub(' subunits', '', complex)) %>% 
  pivot_wider(names_from = complex, 
              values_from = `Complex-Pathway`) %>%
  na.omit() %>%
  column_to_rownames('pathway') %>% 
  as.matrix()

hm <- ComplexHeatmap::pheatmap(etc_pathway_matrix, 
                               cellwidth = 25, 
                               cluster_cols = FALSE,
                               cluster_rows = TRUE,
                               show_rownames = FALSE,
                               scale = 'column',
                               fontsize_col = 10,
                               treeheight_row = 30, 
                               angle_col = '0',
                               heatmap_legend_param = list(title = '', 
                                                           legend_direction = 'horizontal', 
                                                           legend_side = 'top'))

# Supplemental Table 3
## files needed: Supplemental Table 1, sheet 1 == dataframe
dataframe %>% 
  filter(complex %in% c('CI subunits', 'CII subunits', 'CIII subunits', 'CIV subunits', 'CV subunits')) %>%
  dplyr::select(!`Gene-Pathway`) %>%
  pivot_wider(names_from = complex, 
              values_from = `Complex-Pathway`) %>%
  arrange(desc(`CII subunits`)) %>%
  filter(`CII subunits` >= 0.4488337+0.1259972 &  # median +1SD for unique and strong CII correlations
           `CI subunits` <= 0.4488337+0.1259972 &
           `CIII subunits` <= 0.4488337+0.1259972 &
           `CIV subunits` <= 0.4488337+0.1259972 &
           `CV subunits` <= 0.4488337+0.1259972)

# Fig 5a and Extended Data Fig 7a
sdhb_laml <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/LAML/SDHB.csv') 
sdhb_laml$Kaplan.Meier.Class <- factor(sdhb_laml$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_acc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/ACC/SDHB.csv')
sdhb_acc$Kaplan.Meier.Class <- factor(sdhb_acc$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_blca <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/BLCA/SDHB.csv')
sdhb_blca$Kaplan.Meier.Class <- factor(sdhb_blca$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_brca <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/BRCA/SDHB.csv')
sdhb_brca$Kaplan.Meier.Class <- factor(sdhb_brca$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_cesc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/CESC/SDHB.csv')
sdhb_cesc$Kaplan.Meier.Class <- factor(sdhb_cesc$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_chol <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/CHOL/SDHB.csv')
sdhb_chol$Kaplan.Meier.Class <- factor(sdhb_chol$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_coad <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/COAD/SDHB.csv')
sdhb_coad$Kaplan.Meier.Class <- factor(sdhb_coad$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_dlbc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/DLBC/SDHB.csv')
sdhb_dlbc$Kaplan.Meier.Class <- factor(sdhb_dlbc$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_esca <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/ESCA/SDHB.csv')
sdhb_esca$Kaplan.Meier.Class <- factor(sdhb_esca$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_gbm <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/GBM/SDHB.csv')
sdhb_gbm$Kaplan.Meier.Class <- factor(sdhb_gbm$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_hnsc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/HNSC/SDHB.csv')
sdhb_hnsc$Kaplan.Meier.Class <- factor(sdhb_hnsc$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_kich <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/KICH/SDHB.csv')
sdhb_kich$Kaplan.Meier.Class <- factor(sdhb_kich$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_kirc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/KIRC/SDHB.csv')
sdhb_kirc$Kaplan.Meier.Class <- factor(sdhb_kirc$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_kirp <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/KIRP/SDHB.csv')
sdhb_kirp$Kaplan.Meier.Class <- factor(sdhb_kirp$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_lgg <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/LGG/SDHB.csv')
sdhb_lgg$Kaplan.Meier.Class <- factor(sdhb_lgg$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_lihc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/LIHC/SDHB.csv')
sdhb_lihc$Kaplan.Meier.Class <- factor(sdhb_lihc$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_luad <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/LUAD/SDHB.csv')
sdhb_luad$Kaplan.Meier.Class <- factor(sdhb_luad$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_lusc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/LUSC/SDHB.csv')
sdhb_lusc$Kaplan.Meier.Class <- factor(sdhb_lusc$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_meso <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/MESO/SDHB.csv')
sdhb_meso$Kaplan.Meier.Class <- factor(sdhb_meso$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_ov <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/OV/SDHB.csv')
sdhb_ov$Kaplan.Meier.Class <- factor(sdhb_ov$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_paad <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/PAAD/SDHB.csv')
sdhb_paad$Kaplan.Meier.Class <- factor(sdhb_paad$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_pcpg <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/PCPG/SDHB.csv')
sdhb_pcpg$Kaplan.Meier.Class <- factor(sdhb_pcpg$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_prad <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/PRAD/SDHB.csv')
sdhb_prad$Kaplan.Meier.Class <- factor(sdhb_prad$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_read <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/READ/SDHB.csv')
sdhb_read$Kaplan.Meier.Class <- factor(sdhb_read$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_sarc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/SARC/SDHB.csv')
sdhb_sarc$Kaplan.Meier.Class <- factor(sdhb_sarc$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_skcm <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/SKCM/SDHB.csv')
sdhb_skcm$Kaplan.Meier.Class <- factor(sdhb_skcm$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_stad <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/STAD/SDHB.csv')
sdhb_stad$Kaplan.Meier.Class <- factor(sdhb_stad$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_tgct <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/TGCT/SDHB.csv')
sdhb_tgct$Kaplan.Meier.Class <- factor(sdhb_tgct$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_thca <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/THCA/SDHB.csv')
sdhb_thca$Kaplan.Meier.Class <- factor(sdhb_thca$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_thym <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/THYM/SDHB.csv')
sdhb_thym$Kaplan.Meier.Class <- factor(sdhb_thym$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_ucec <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/UCEC/SDHB.csv')
sdhb_ucec$Kaplan.Meier.Class <- factor(sdhb_ucec$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_ucs <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/UCS/SDHB.csv')
sdhb_ucs$Kaplan.Meier.Class <- factor(sdhb_ucs$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdhb_uvm <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/UVM/SDHB.csv')
sdhb_uvm$Kaplan.Meier.Class <- factor(sdhb_uvm$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))

sdha_laml <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/LAML/SDHA.csv')
sdha_laml$Kaplan.Meier.Class <- factor(sdha_laml$Kaplan.Meier.Class, levels = c("low expression", 'high expression'))
sdha_laml %>% summarise(mean(Gene.Value))

sdha_acc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/ACC/SDHA.csv')

sdha_blca <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/BLCA/SDHA.csv')

sdha_brca <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/BRCA/SDHA.csv')

sdha_cesc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/CESC/SDHA.csv')

sdha_chol <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/CHOL/SDHA.csv')

sdha_coad <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/COAD/SDHA.csv')

sdha_dlbc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/DLBC/SDHA.csv')

sdha_esca <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/ESCA/SDHA.csv')

sdha_gbm <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/GBM/SDHA.csv')

sdha_hnsc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/HNSC/SDHA.csv')

sdha_kich <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/KICH/SDHA.csv')

sdha_kirc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/KIRC/SDHA.csv')

sdha_kirp <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/KIRP/SDHA.csv')

sdha_lgg <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/LGG/SDHA.csv')

sdha_lihc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/LIHC/SDHA.csv')

sdha_luad <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/LUAD/SDHA.csv')

sdha_lusc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/LUSC/SDHA.csv')

sdha_meso <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/MESO/SDHA.csv')

sdha_ov <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/OV/SDHA.csv')

sdha_paad <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/PAAD/SDHA.csv')

sdha_pcpg <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/PCPG/SDHA.csv')

sdha_prad <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/PRAD/SDHA.csv')

sdha_read <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/READ/SDHA.csv')

sdha_sarc <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/SARC/SDHA.csv')

sdha_skcm <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/SKCM/SDHA.csv')

sdha_stad <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/STAD/SDHA.csv')

sdha_tgct <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/TGCT/SDHA.csv')

sdha_thca <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/THCA/SDHA.csv')

sdha_thym <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/THYM/SDHA.csv')

sdha_ucec <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/UCEC/SDHA.csv')

sdha_ucs <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/UCS/SDHA.csv')

sdha_uvm <- read.csv('https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/csvs/rnaseq/UVM/SDHA.csv')

# Kaplan Meier analysis

tcga_surv <- sdhb_laml
model_surv <- survfit(Surv(Time, Censor) ~ Kaplan.Meier.Class, data = tcga_surv)
names(model_surv$strata) <- gsub("Kaplan.Meier.Class=", "", names(model_surv$strata))

plot <- ggsurvplot(
  model_surv,
  palette = ggsci::pal_npg(alpha = 1)(6),
  legend = "none",
  legend.title = element_blank(),
  surv.median.line = "none",
  data = tcga_surv, 
  pval = TRUE,
  pval.method = FALSE,
  pval.coord = c(0, 0.025),
  risk.table = TRUE,
  tables.theme = theme_cleantable(),
  xlim = c(0, max(tcga_surv$Time)+100),
  title = ''
) +
  ylab('Overall Survival') +
  xlab('Time (Days)')

plot

# Figure 5i
## Beat AML regression analysis
## dataframe needed from Beat AML project is SampleID, Inhibitor, AUC, Normalized SDHB Expression
data <- beataml_sdhb_clean %>% dplyr::select(SampleID, Inhibitor, Area_under_curve, Normalized_Exprs) 

lm_results <- data.frame(
  inhibitor = character(),
  p_value = numeric(),
  coefficient = numeric(),
  adj_p_value = numeric(),
  stringsAsFactors = FALSE
)

inhibitors <- unique(data$Inhibitor)

for (inhib in inhibitors) {
  inhib_data <- subset(data, Inhibitor == inhib)
  
  model <- lm(Area_under_curve ~ Normalized_Exprs, data = inhib_data)
  
  model_summary <- summary(model)
  p_value <- coef(model_summary)[2, 4]
  coefficient <- coef(model_summary)[2, 1]
  
  lm_results <- rbind(lm_results, data.frame(
    inhibitor = inhib,
    p_value = p_value,
    coefficient = coefficient,
    stringsAsFactors = FALSE
  ))
}

lm_results$adj_p_value <- p.adjust(lm_results$p_value, method = "fdr")

# Figure 5l
## Cox Proportional Hazard models used to create dataframe with HR and confidence intervals
## use SDHB expression data frames from above in 'data' entry
coxph(Surv(Time, Censor) ~ Kaplan.Meier.Class, data = sdhb_laml) %>%
  finalfit::fit2df()
