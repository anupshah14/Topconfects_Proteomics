#########

# Script to apply topconfects on limma object

#########

library(DEP)
library(tidyverse)
# BiocManager::install("topconfects")
library(topconfects)
library(SummarizedExperiment)
library(limma)
source("functions.R")
set.seed(123454)

protein_demo<-read.table("proteinGroups.txt", header = T, sep = "\t")

exp_design_demo<-read.table("exp_design.txt", 
                                     header = T, sep = "\t",
                                     stringsAsFactors=FALSE)
data_demo <- protein_demo %>% 
  dplyr::filter(Reverse!="+", Potential.contaminant!="+",
                Only.identified.by.site!="+", Razor...unique.peptides>=2)

data_demo<- DEP:::make_unique(data_demo,"Gene.names","Protein.IDs",delim=";")

lfq_cols_demo<-grep("LFQ.", colnames(data_demo))

demo_se<-DEP:::make_se(data_demo,lfq_cols_demo,exp_design_demo)


demo_filter<-DEP:::filter_missval(demo_se, thr = 5)


demo_imp<-DEP:::impute(demo_filter,fun="man",shift=1.8,scale=0.3)

demo_diff_limma<-test_limma(demo_imp,type='all')

demo_dep<-DEP:::add_rejections(demo_diff_limma,alpha = 0.05,lfc = 1)

demo_limma_results<-get_results_proteins(demo_dep)


demo_design<-model.matrix(~exp_design_demo$condition)
demo_design[,]


demo_fit<-lmFit(assay(demo_imp), design = demo_design)

demo_confects<-limma_confects(demo_fit, coef=2, fdr = 0.05)
demo_confects_trend<-limma_confects(demo_fit, coef=2, fdr = 0.05, trend = TRUE)
confects_plot(demo_confects,n=10)
plotSA(demo_fit)


demo_fit_eb <- eBayes(demo_fit)
demo_eb_trend<-eBayes(demo_fit, trend = TRUE)
top_demo <- topTable(demo_fit_eb, coef=2, n=Inf)
top_demo_trend<-topTable(demo_eb_trend, coef=2, n=Inf)

rank_rank_plot(demo_confects$table$name, rownames(top_demo), "limma_confects", "topTable")
rank_rank_plot(demo_confects_trend$table$name, rownames(top_demo_trend), 
               "limma_confects_trend", "topTable_trend")

head(assay(demo_imp))


##### T test

demo_t_test<-apply(assay(demo_imp), 1, function (x) t.test(x[c(2,4,6,8,10,12,14,16,18,20)],
                                                    x[c(1,3,5,7,9,11,13,15,17,19)],
                                                    var.equal = FALSE))

effect_demo<-unlist(lapply(demo_t_test, function(x) x$estimate[1]- x$estimate[2]))
pvalue_demo<-unlist(lapply(demo_t_test, function(x) x$p.value))
fdr_demo<-p.adjust(pvalue_demo, method = "BH")
fdr_demo<-sort(fdr_demo, decreasing = FALSE)

demo_std_error<-unlist(lapply(demo_t_test, function(x) x$stderr))
demo_df<-unlist(lapply(demo_t_test, function(x) x$parameter["df"]))


t_test_confects<-normal_confects_names(effect_demo, demo_std_error, demo_df, fdr= 0.05,
                                       full = TRUE)

confects_plot(t_test_confects, n=10)

rank_rank_plot(t_test_confects$table$name, names(fdr_demo), "t-test_confects", "t-test")

rank_rank_plot(demo_confects$table$name, names(fdr_demo), "limma_confects", 
               "t-test")
rank_rank_plot(demo_confects$table$name, t_test_confects$table$name, "limma_confects", 
               "t-test_confects")





rank_rank_plot(demo_confects$table$name, names(fdr_demo), "Topconfects", 
               "t-test")

names(fdr_demo)

t_test_confects$table$name


plotSA(demo_fit)
plotSA(demo_fit_eb)
plotSA(demo_eb_trend)


#### Plot confect values



demo_confects$table %>% filter(confect!="NA") %>%
  ggplot(., aes(confect) ) +
  geom_density()

demo_confects$table %>% filter(confect!="NA") %>%
  ggplot(., aes(effect) ) +
  geom_density()

demo_confects$table %>% 
  ggplot(., aes(effect) ) +
  geom_density()

##### Enrichment =====
#install.packages("gprofiler2")
library(gprofiler2)


## vector with proteins with confects value
de_demo<-demo_confects$table %>% filter(confect!="NA") %>%
  select(name)

demo_enriched<-gost(de_demo$name, organism = "hsapiens", ordered_query = TRUE,
     multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
     measure_underrepresentation = FALSE, evcodes = FALSE,
     user_threshold = 0.05, correction_method = "fdr",
     domain_scope = "annotated", custom_bg = NULL,
     numeric_ns = "", sources = "GO:CC")


### Confects top 20
demo_enriched$result %>% filter(significant==TRUE) %>% 
  head(., 20) %>%
  ggplot(., aes(term_name, -log10(p_value))) +
           geom_col() +
           coord_flip()


### t-test top 20
demo_fdr_enrichment<- gost(names(fdr_demo[fdr_demo<0.05]), organism = "hsapiens", ordered_query = TRUE,
                         multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                         measure_underrepresentation = FALSE, evcodes = FALSE,
                         user_threshold = 0.05, correction_method = "fdr",
                         domain_scope = "annotated", custom_bg = NULL,
                         numeric_ns = "", sources = "GO:CC")


demo_fdr_enrichment$result %>% filter(significant==TRUE) %>% 
  head(., 20) %>%
  ggplot(., aes(term_name, -log10(p_value))) +
  geom_col() +
  coord_flip()


plot_protein(demo_dep, protein = head(demo_confects$table$name,12), type = "boxplot") +
  scale_color_brewer(palette = "Paired")
  

plot_protein(demo_dep, protein = head(names(fdr_demo),12), type = "boxplot") + 
  scale_color_brewer(palette = "Paired")



demo_confects$table %>% filter(confect!="NA") %>% 
  filter(abs(effect) >= 1.00) %>% 
  ggplot(., aes(confect) ) +
  geom_density()


demo_confects$table %>% filter(confect!="NA") %>%
  filter(abs(effect) >= 1.00) %>%
  ggplot(., aes(effect) ) +
  geom_density()



demo_confects$table %>% filter(confect!="NA") %>%
  filter(abs(effect) >= 1.00) %>%
  reshape2::melt(., id.vars=c("name"), measure.vars=c("effect", "confect")) %>%
  ggplot(., aes(value, color= variable)) +
  geom_density( ) 



demo_confects$table %>% filter(confect!="NA") %>%
  filter(abs(effect) >= 1.00) %>%
  reshape2::melt(., id.vars=c("name"), measure.vars=c("effect", "confect")) %>%
  ggplot(., aes(value, fill= variable, color=variable)) +
  geom_density(alpha=0.5 ) + 
  labs(title=" Prostate cancer Dataset (n=10)",
       x="log 2 Fold change",
       y= "Density")+ 
  theme_DEP1()





#### Histogram combined
fill_color <- c("#40b8d0", "#b2d183")

demo_confects$table %>% filter(confect!="NA") %>%
  filter(abs(effect) >= 1.00) %>%
  reshape2::melt(., id.vars=c("name"), measure.vars=c("effect", "confect")) %>%
  ggplot(., aes(value, color= variable)) +
  geom_histogram( bins=100, fill="white" ) + 
  labs(title=" Prostate cancer Dataset (n=10)",
       x="log 2 Fold change",
       y= "Counts")+ 
  scale_color_manual(values=fill_color)+
  theme_DEP1()



rank_rank_plot(demo_confects$table$name, names(fdr_demo), "topconfects", 
               "t-test", n=20)



plot_protein(demo_dep, protein = head(demo_confects$table$name,9), type = "boxplot") +
  scale_color_brewer(palette = "Paired")

plot_protein(demo_dep, protein = head(names(fdr_demo),9), type = "boxplot") + 
  scale_color_brewer(palette = "Paired")


ggplot() +
  geom_point()+
  #xlim(-2, 5)+
  ylim(0,3.5)+
  geom_vline(xintercept = 0)+
  scale_x_continuous(limits=c(-2,4), breaks = seq(-2, 4, by = 1))+
  labs(x="Effect")+
  theme_bw()



