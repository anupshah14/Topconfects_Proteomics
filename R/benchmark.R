library(readxl)

benchmark_sheet<-readxl::read_excel("benchmark_tmt_dia.xls", sheet = "AbsDeltaFC" )

dia_quant<-readxl::read_excel("benchmark_tmt_dia.xls", sheet = "ProteinQuantities_BGS-DIA",
                              na="NA",
                              skip = 1)

tmt_quant<-readxl::read_excel("benchmark_tmt_dia.xls", sheet = "ProteinQuantities_BGS-TMT",
                              na="NA",
                              skip = 1)

exp_benchmark<-read.table("exp_design_benchmark.txt",
                          header = T, sep = "\t",
                          stringsAsFactors=FALSE)
exp_benchmark_sm<-read.table("exp_design_benchmark_s1s5.txt",
                          header = T, sep = "\t",
                          stringsAsFactors=FALSE)
benchmark_design<-model.matrix(~exp_benchmark$condition)
benchmark_design_sm<-model.matrix(~exp_benchmark_sm$condition)
benchmark_design[,]

dia_quant_full<-dia_quant[complete.cases(dia_quant[,2:11]),]
rownames(dia_quant_full)<-dia_quant_full$PROTEIN

benchmark_dia_fit<-lmFit(log2(dia_quant_full[,2:11]), design = benchmark_design)
benchmark_dia_sm_fit<-lmFit(log2_dia_sm, 
                            design = benchmark_design_sm)
benchmark_dia_confects<-limma_confects(benchmark_dia_fit, coef=5, fdr = 0.05 )
benchmark_dia_sm_confects<-limma_confects(benchmark_dia_sm_fit, coef=2, fdr = 0.05 )

head(benchmark_dia_confects$table)





### t-test

log2_dia_sm<-log2(dia_quant_full[,c(2,3,10,11)])
rownames(log2_dia_sm)<-dia_quant_full$PROTEIN

benchmark_dia_sm_t_test<-apply(log2_dia_sm, 1, 
                               function (x) t.test(x[c(1,2)],
                                               x[c(3,4)],
                                  var.equal = FALSE))

effect_benchmark_dia_sm<-unlist(lapply(benchmark_dia_sm_t_test, function(x) x$estimate[1]- x$estimate[2]))
pvalue_benchmark_dia_sm<-unlist(lapply(benchmark_dia_sm_t_test, function(x) x$p.value))
fdr_benchmark_dia_sm<-p.adjust(pvalue_benchmark_dia_sm, method = "BH")
fdr_benchmark_dia_sm<-sort(fdr_benchmark_dia_sm, decreasing = FALSE)
benchmark_dia_sm_std_error<-unlist(lapply(benchmark_dia_sm_t_test, function(x) x$stderr))
benchmark_dia_sm_df<-unlist(lapply(benchmark_dia_sm_t_test, function(x) x$parameter["df"]))


benchmark_dia_sm_t_test_confects<-normal_confects_names(effect_benchmark_dia_sm, 
                                                        benchmark_dia_sm_std_error, 
                                                        benchmark_dia_sm_df, fdr= 0.05, 
                                                    full = TRUE)


rank_rank_plot(benchmark_dia_sm_confects$table$name, names(fdr_benchmark_dia_sm),
               "topconfects", "t-test")
benchmark_dia_sm_confects$table$name %>% head()
names(fdr_benchmark_dia_sm) %>% head()


### get pus proteins


benchmark_dia_sm_confects$table %>% filter(confect!="NA") %>% select(name)

confects_plot(benchmark_dia_sm_confects, n=12)


### === TMT =====### 

tmt_quant_full<-tmt_quant[complete.cases(tmt_quant[,2:11]),]
rownames()
bm_tmt_exp<-read.table("tmt_exp_s1s5",
                             header = T, sep = "\t",
                             stringsAsFactors=FALSE)

bm_tmt_design_sm<-model.matrix(~bm_tmt_exp$condition)


tmt_sm<-tmt_quant_full[,c(8,11,4,10)]
tmt_sm<-log2(tmt_sm)
colnames(tmt_sm)
rownames(tmt_sm)<-tmt_quant_full$PROTEIN
benchmark_tmt_sm_fit<-lmFit(tmt_sm, 
                            design = bm_tmt_design_sm)

bm_tmt_sm_confects<-limma_confects(benchmark_tmt_sm_fit, coef = 2)

bm_tmt_sm_confects$table %>% filter(confect!="NA") %>% select(name)

confects_plot(bm_tmt_sm_confects, n=12)
