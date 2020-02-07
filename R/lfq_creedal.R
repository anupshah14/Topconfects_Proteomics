#########

# Script to apply topconfects on limma object for creedal et. al dataset

#########

library(DEP)
library(tidyverse)
# BiocManager::install("topconfects")
library(topconfects)
library(SummarizedExperiment)
library(limma)
source("functions.R")
set.seed(123454)

protein_creedal<-read.table("proteinGroups_creedal.txt", header = T, sep = "\t")

exp_design_creedal<-read.table("exp_design_creedal.txt", 
                            header = T, sep = "\t",
                            stringsAsFactors=FALSE)
data_creedal <- protein_creedal %>% 
  dplyr::filter(Reverse!="+", Potential.contaminant!="+",
                Only.identified.by.site!="+", Razor...unique.peptides>=2)

data_creedal<- DEP:::make_unique(data_creedal,"Gene.names","Protein.IDs",delim=";")

lfq_cols_creedal<-grep("LFQ.", colnames(data_creedal))

creedal_se<-DEP:::make_se(data_creedal,lfq_cols_creedal,exp_design_creedal)


creedal_filter<-DEP:::filter_missval(creedal_se, thr = 1)


creedal_imp<-DEP:::impute(creedal_filter,fun="man",shift=1.8,scale=0.3)

creedal_diff_limma<-test_limma(creedal_imp,type='all')

creedal_dep<-DEP:::add_rejections(creedal_diff_limma,alpha = 0.05,lfc = 1)

creedal_limma_results<-get_results_proteins(creedal_dep)


creedal_design<-model.matrix(~exp_design_creedal$condition)
creedal_design[,]

creedal_fit<-lmFit(assay(creedal_imp), design = creedal_design)

creedal_confects<-limma_confects(creedal_fit, coef=2, fdr = 0.05)
creedal_confects_trend<-limma_confects(creedal_fit, coef=2, fdr = 0.05, trend = TRUE)
confects_plot(creedal_confects,n=10)

se_genes_creedal<-creedal_confects$table %>% filter(confect!="NA") %>% 
  select(name) %>% as.vector(.)






creedal_fit_eb <- eBayes(creedal_fit)
creedal_eb_trend<-eBayes(creedal_fit, trend = TRUE)
top_creedal <- topTable(creedal_fit_eb, coef=2, n=Inf)
top_creedal_trend<-topTable(creedal_eb_trend, coef=2, n=Inf)

rank_rank_plot(creedal_confects$table$name, rownames(top_creedal), "limma_confects", "topTable")
rank_rank_plot(creedal_confects_trend$table$name, rownames(top_creedal_trend), 
               "limma_confects_trend", "topTable_trend")

head(assay(creedal_imp))


##### T test

creedal_t_test<-apply(assay(creedal_imp), 1, function (x) t.test(x[c(4,5,6)],
                                                           x[c(1,2,3)],
                                                           var.equal = FALSE))

effect_creedal<-unlist(lapply(creedal_t_test, function(x) x$estimate[1]- x$estimate[2]))
pvalue_creedal<-unlist(lapply(creedal_t_test, function(x) x$p.value))
fdr_creedal<-p.adjust(pvalue_creedal, method = "BH")
fdr_creedal<-sort(fdr_creedal, decreasing = FALSE)

creedal_std_error<-unlist(lapply(creedal_t_test, function(x) x$stderr))
creedal_df<-unlist(lapply(creedal_t_test, function(x) x$parameter["df"]))


t_test_confects<-normal_confects_names(effect_creedal, creedal_std_error, 
                                       creedal_df, fdr= 0.05, full = TRUE)

confects_plot(t_test_confects, n=10)

rank_rank_plot(t_test_confects$table$name, names(fdr_creedal), "t-test_confects", "t-test")

rank_rank_plot(creedal_confects$table$name, names(fdr_creedal), "topconfects", "t-test")


names(fdr_creedal)

t_test_confects$table$name


plotSA(creedal_fit)
plotSA(creedal_fit_eb)
plotSA(creedal_eb_trend)


creedal_confects$table %>% filter(confect!="NA") %>%
  ggplot(., aes(confect) ) +
  geom_histogram(bins=100)




plot_protein(creedal_dep, protein = head(creedal_confects$table$name,16), type = "boxplot") 


plot_protein(creedal_dep, protein = head(names(fdr_creedal),16), type = "boxplot")




creedal_confects$table %>% filter(confect!="NA") %>%
  filter(abs(effect) >= 1.00) %>%
  reshape2::melt(., id.vars=c("name"), measure.vars=c("effect", "confect")) %>%
  ggplot(., aes(value, fill= variable, color=variable)) +
  geom_density(alpha=0.8 ) + 
  labs(title="Creedal et. al (n=3)",
       x="log 2 Fold change",
       y= "Density")+ 
  scale_fill_manual(values=fill_color)+
  scale_color_manual(values=fill_color)+
  theme_DEP1()


rank_rank_plot(creedal_confects$table$name, names(fdr_creedal), "topconfects", "t-test",
               n=20)
