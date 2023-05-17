## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----message=FALSE, , results='hide', include=FALSE-------------------------------------------------------------------------------------------
pacman::p_load(dplyr,
               tibble,
               readr,
               stringr,
               data.table,
               ggplot2,
               sjPlot, 
               tidyr,
               ggpubr,
               ggmosaic)



## ---------------------------------------------------------------------------------------------------------------------------------------------
dat = fread('DGEresults/DGE_CPM_data.csv')

head(dat)

dat$GeneID %>% unique() %>% length()

anno = fread('DGEresults/annotation_combined.csv')


## ---------------------------------------------------------------------------------------------------------------------------------------------
expressed = dat %>% 
  as.data.frame() %>% 
  mutate(GeneID = as.numeric(GeneID)) %>% 
  select(GeneID) %>% 
  left_join(anno) %>%
  distinct(GeneID, Parent_haplotype) %>% 
  count(Parent_haplotype, name = 'Expressed_gene_count')

p = ggplot(expressed, aes(Parent_haplotype, Expressed_gene_count, fill = Parent_haplotype))+
  geom_col(color = 'black')+
  labs(y='Expressed genes')

p



## ---------------------------------------------------------------------------------------------------------------------------------------------
expressed


## ---------------------------------------------------------------------------------------------------------------------------------------------
impCG = fread('DGEresults/Limma_results_table_CG.csv') %>% 
  mutate(Gene_type = ifelse(logFC > 0, 'Susceptiblity', 'Resistance'),
         Parent_haplotype = factor(Parent_haplotype, levels = c('J.regia', 'J.microcarpa')))


## ---------------------------------------------------------------------------------------------------------------------------------------------
ASEsums = impCG %>% distinct(GeneID, Parent_haplotype) %>% count(Parent_haplotype, name = 'DEG_gene_count')

p1 = ggplot(ASEsums, aes(Parent_haplotype, DEG_gene_count, fill = Parent_haplotype))+
  geom_col(color = 'black')+
  theme(legend.position = 'none')+
  labs(y='DEGs')+
  ggtitle('Total CG DEGs')

p1



## ---------------------------------------------------------------------------------------------------------------------------------------------
sum(ASEsums$DEG_gene_count)


## ---------------------------------------------------------------------------------------------------------------------------------------------
impCG_pos = impCG %>%
  filter(logFC > 0) %>%
  distinct(GeneID, Parent_haplotype) %>% 
  count(Parent_haplotype, name = 'Susceptibility_genes')

sum(impCG_pos$n)

p2 = ggplot(impCG_pos, aes(Parent_haplotype, Susceptibility_genes, fill = Parent_haplotype))+
  geom_col(color = 'black')+
  theme(legend.position = 'none')+
  labs(y='DEGs')+
  ggtitle('CG DEGs (+)')

p2


## ---------------------------------------------------------------------------------------------------------------------------------------------
binom.test.v = Vectorize(binom.test)

impCG_pos_test = impCG_pos %>% 
  mutate(observed_sum = sum(Susceptibility_genes)) %>% 
  rename(observed = Susceptibility_genes) %>% 
  left_join(expressed, by = 'Parent_haplotype') %>% 
  left_join(ASEsums, by = 'Parent_haplotype') %>% 
  mutate(Analysis = 'CG_pos',
         Expressed_gene_count_sum = sum(Expressed_gene_count),
         DEG_gene_count_sum = sum(DEG_gene_count),
         Background_prob = DEG_gene_count/DEG_gene_count_sum,
         expected_observed = Background_prob*observed_sum,
         log_FE = log2(observed/expected_observed)) %>% 
  group_by(Parent_haplotype) %>% 
  mutate(Binom_p_val = binom.test.v(x = observed,
                                n = observed_sum,
                                p = Background_prob)[[3]],
         Binom_estimate = binom.test.v(x = observed,
                                n = observed_sum,
                                p = Background_prob)[[5]][[1]]) %>% 
  relocate(Analysis)



## ---------------------------------------------------------------------------------------------------------------------------------------------
impCG_neg = impCG %>%
  filter(logFC < 0) %>%
  distinct(GeneID, Parent_haplotype) %>% 
  count(Parent_haplotype, name = 'Resistance_genes')

p3 = ggplot(impCG_neg, aes(Parent_haplotype, Resistance_genes, fill = Parent_haplotype))+
  geom_col(color = 'black')+
  theme(legend.position = 'none')+
  labs(y='DEGs')+
  ggtitle('CG DEGs (-)')

p3


## ---------------------------------------------------------------------------------------------------------------------------------------------
impCG_neg_test = impCG_neg %>% 
  mutate(observed_sum = sum(Resistance_genes)) %>% 
  rename(observed = Resistance_genes) %>% 
  left_join(expressed, by = 'Parent_haplotype') %>% 
  left_join(ASEsums, by = 'Parent_haplotype') %>% 
  mutate(Analysis = 'CG_neg',
         Expressed_gene_count_sum = sum(Expressed_gene_count),
         DEG_gene_count_sum = sum(DEG_gene_count),
         Background_prob = DEG_gene_count/DEG_gene_count_sum,
         expected_observed = Background_prob*observed_sum,
         log_FE = log2(observed/expected_observed)) %>% 
  group_by(Parent_haplotype) %>% 
  mutate(Binom_p_val = binom.test.v(x = observed,
                                n = observed_sum,
                                p = Background_prob)[[3]],
         Binom_estimate = binom.test.v(x = observed,
                                n = observed_sum,
                                p = Background_prob)[[5]][[1]]) %>% 
  relocate(Analysis)


View(impCG_neg_test)


## ---------------------------------------------------------------------------------------------------------------------------------------------
t = impCG_pos %>% 
  left_join(impCG_neg) %>% 
  arrange(desc(Susceptibility_genes)) %>% 
  column_to_rownames(var = 'Parent_haplotype') %>% 
  t()

t


## ---------------------------------------------------------------------------------------------------------------------------------------------
t = impCG_pos %>% 
  left_join(impCG_neg) %>%
  arrange(desc(Parent_haplotype)) %>% 
  relocate(Resistance_genes, .before = Susceptibility_genes) %>% 
  column_to_rownames(var = 'Parent_haplotype') %>% 

  t()

t


## ---------------------------------------------------------------------------------------------------------------------------------------------
sum(t)


## ---------------------------------------------------------------------------------------------------------------------------------------------
mosaicplot(t,
           color = T)


## ---------------------------------------------------------------------------------------------------------------------------------------------
fisher_res = fisher.test(t)

fisher_res



## ---------------------------------------------------------------------------------------------------------------------------------------------
t = impCG_pos %>% 
  left_join(impCG_neg) %>%
  arrange(desc(Parent_haplotype)) %>% 
  relocate(Resistance_genes, .before = Susceptibility_genes) %>% 
  column_to_rownames(var = 'Parent_haplotype') %>% 
  t()

t


## ---------------------------------------------------------------------------------------------------------------------------------------------
sum(t)


## ---------------------------------------------------------------------------------------------------------------------------------------------
mosaicplot(t,
           color = T)


## ---------------------------------------------------------------------------------------------------------------------------------------------
fisher_res = fisher.test(t)

fisher_res



## ---------------------------------------------------------------------------------------------------------------------------------------------
uniq_impCG = impCG %>% 
  distinct(GeneID, Parent_haplotype, Gene_type) %>% 
  mutate(Pathogen = 'A. tumefaciens',
         Stats = paste('Fisher p-value: ',
                     round(fisher_res$p.val, 4),
                     '\n', 'Odds ratio: ',
                     round(fisher_res$estimate, 3), sep = ''))

A = ggplot(uniq_impCG)+
  geom_mosaic(aes(x = product(Gene_type, Parent_haplotype), fill = Parent_haplotype), color = 'black')+
  theme_grey(base_size = 14)+
    theme(legend.position = 'none',
          strip.text = element_text(face = 'italic'))+
  facet_wrap(~Pathogen+Stats)

A


## ---------------------------------------------------------------------------------------------------------------------------------------------
a = ggarrange(p2, p3, nrow = 2)

a


## ---------------------------------------------------------------------------------------------------------------------------------------------
impPHY = fread('DGEresults/Limma_results_table_PHY.csv') %>% 
 mutate(Gene_type = ifelse(logFC > 0, 'Susceptiblity', 'Resistance'),
         Parent_haplotype = factor(Parent_haplotype, levels = c('J.regia', 'J.microcarpa')))


## ---------------------------------------------------------------------------------------------------------------------------------------------
ASEsums = impPHY %>% distinct(GeneID, Parent_haplotype) %>% count(Parent_haplotype, name = 'DEG_gene_count')

p1 = ggplot(ASEsums, aes(Parent_haplotype, DEG_gene_count, fill = Parent_haplotype))+
  geom_col(color = 'black')+
  theme(legend.position = 'none')+
  labs(y='DEGs')+
  ggtitle('Total PHY DEGs')

p1



## ---------------------------------------------------------------------------------------------------------------------------------------------
sum(ASEsums$DEG_gene_count)


## ---------------------------------------------------------------------------------------------------------------------------------------------
impPHY_pos = impPHY %>%
  filter(logFC > 0) %>%
  distinct(GeneID, Parent_haplotype) %>% 
  count(Parent_haplotype, name = 'Susceptibility_genes')

sum(impPHY_pos$n)

p2 = ggplot(impPHY_pos, aes(Parent_haplotype, Susceptibility_genes, fill = Parent_haplotype))+
  geom_col(color = 'black')+
  theme(legend.position = 'none')+
  labs(y='DEGs')+
  ggtitle('PHY DEGs (+)')

p2


## ---------------------------------------------------------------------------------------------------------------------------------------------
binom.test.v = Vectorize(binom.test)

impPHY_pos_test = impPHY_pos %>% 
  mutate(observed_sum = sum(Susceptibility_genes)) %>% 
  rename(observed = Susceptibility_genes) %>% 
  left_join(expressed, by = 'Parent_haplotype') %>% 
  left_join(ASEsums, by = 'Parent_haplotype') %>% 
  mutate(Analysis = 'PHY_pos',
         Expressed_gene_count_sum = sum(Expressed_gene_count),
         DEG_gene_count_sum = sum(DEG_gene_count),
         Background_prob = DEG_gene_count/DEG_gene_count_sum,
         expected_observed = Background_prob*observed_sum,
         log_FE = log2(observed/expected_observed)) %>% 
  group_by(Parent_haplotype) %>% 
  mutate(Binom_p_val = binom.test.v(x = observed,
                                n = observed_sum,
                                p = Background_prob)[[3]],
         Binom_estimate = binom.test.v(x = observed,
                                n = observed_sum,
                                p = Background_prob)[[5]][[1]]) %>% 
  relocate(Analysis)



## ---------------------------------------------------------------------------------------------------------------------------------------------
impPHY_neg = impPHY %>%
  filter(logFC < 0) %>%
  distinct(GeneID, Parent_haplotype) %>% 
  count(Parent_haplotype, name = 'Resistance_genes')

p3 = ggplot(impPHY_neg, aes(Parent_haplotype, Resistance_genes, fill = Parent_haplotype))+
  geom_col(color = 'black')+
  theme(legend.position = 'none')+
  labs(y='DEGs')+
  ggtitle('PHY DEGs (-)')

p3


## ---------------------------------------------------------------------------------------------------------------------------------------------
impPHY_neg_test = impPHY_neg %>% 
  mutate(observed_sum = sum(Resistance_genes)) %>% 
  rename(observed = Resistance_genes) %>% 
  left_join(expressed, by = 'Parent_haplotype') %>% 
  left_join(ASEsums, by = 'Parent_haplotype') %>% 
  mutate(Analysis = 'PHY_neg',
         Expressed_gene_count_sum = sum(Expressed_gene_count),
         DEG_gene_count_sum = sum(DEG_gene_count),
         Background_prob = DEG_gene_count/DEG_gene_count_sum,
         expected_observed = Background_prob*observed_sum,
         log_FE = log2(observed/expected_observed)) %>% 
  group_by(Parent_haplotype) %>% 
  mutate(Binom_p_val = binom.test.v(x = observed,
                                n = observed_sum,
                                p = Background_prob)[[3]],
         Binom_estimate = binom.test.v(x = observed,
                                n = observed_sum,
                                p = Background_prob)[[5]][[1]]) %>% 
  relocate(Analysis)


View(impPHY_neg_test)


## ---------------------------------------------------------------------------------------------------------------------------------------------
t = impPHY_pos %>% 
  left_join(impPHY_neg) %>% 
  arrange(desc(Susceptibility_genes)) %>%
  column_to_rownames(var = 'Parent_haplotype') %>% 
  t()

t


## ---------------------------------------------------------------------------------------------------------------------------------------------
sum(t)


## ---------------------------------------------------------------------------------------------------------------------------------------------
mosaicplot(t,
           color = T)


## ---------------------------------------------------------------------------------------------------------------------------------------------
fisher_res = fisher.test(t)

fisher_res



## ---------------------------------------------------------------------------------------------------------------------------------------------
uniq_impPHY = impPHY %>% 
  distinct(GeneID, Parent_haplotype, Gene_type) %>% 
  mutate(Pathogen = 'Phytophthora .spp',
         Stats = paste('Fisher p-value: ',
                     round(fisher_res$p.val, 3),
                     '\n', 'Odds ratio: ',
                     round(fisher_res$estimate, 3), sep = ''))

B = ggplot(uniq_impPHY)+
  geom_mosaic(aes(x = product(Gene_type, Parent_haplotype), fill = Parent_haplotype), color = 'black')+
  labs(title = substitute(paste(italic('Phytophthora'))),
       tag = paste('.spp ',
                   'Fisher p-value: ',
                     round(fisher_res$p.val, 3),
                     ', ', 'Odds ratio: ',
                     round(fisher_res$estimate, 3), sep = ''))+
  theme(legend.position = 'none', plot.tag.position = c(0.56, 0.98))

B


## ---------------------------------------------------------------------------------------------------------------------------------------------
a = ggarrange(p2, p3, nrow = 2)

a


## ---------------------------------------------------------------------------------------------------------------------------------------------
impNEM = fread('DGEresults/Limma_results_table_NEM.csv') %>% 
 mutate(Gene_type = ifelse(logFC > 0, 'Susceptiblity', 'Resistance'),
         Parent_haplotype = factor(Parent_haplotype, levels = c('J.regia', 'J.microcarpa')))


## ---------------------------------------------------------------------------------------------------------------------------------------------
ASEsums = impNEM %>% distinct(GeneID, Parent_haplotype) %>% count(Parent_haplotype, name = 'DEG_gene_count')

p1 = ggplot(ASEsums, aes(Parent_haplotype, DEG_gene_count, fill = Parent_haplotype))+
  geom_col(color = 'black')+
  theme(legend.position = 'none')+
  labs(y='DEGs')+
  ggtitle('Total NEM DEGs')

p1



## ---------------------------------------------------------------------------------------------------------------------------------------------
sum(ASEsums$DEG_gene_count)


## ---------------------------------------------------------------------------------------------------------------------------------------------
impNEM_pos = impNEM %>%
  filter(logFC > 0) %>%
  distinct(GeneID, Parent_haplotype) %>% 
  count(Parent_haplotype, name = 'Susceptibility_genes')

sum(impNEM_pos$n)

p2 = ggplot(impNEM_pos, aes(Parent_haplotype, Susceptibility_genes, fill = Parent_haplotype))+
  geom_col(color = 'black')+
  theme(legend.position = 'none')+
  labs(y='DEGs')+
  ggtitle('NEM DEGs (+)')

p2


## ---------------------------------------------------------------------------------------------------------------------------------------------
binom.test.v = Vectorize(binom.test)

impNEM_pos_test = impNEM_pos %>% 
  mutate(observed_sum = sum(Susceptibility_genes)) %>% 
  rename(observed = Susceptibility_genes) %>% 
  left_join(expressed, by = 'Parent_haplotype') %>% 
  left_join(ASEsums, by = 'Parent_haplotype') %>% 
  mutate(Analysis = 'NEM_pos',
         Expressed_gene_count_sum = sum(Expressed_gene_count),
         DEG_gene_count_sum = sum(DEG_gene_count),
         Background_prob = DEG_gene_count/DEG_gene_count_sum,
         expected_observed = Background_prob*observed_sum,
         log_FE = log2(observed/expected_observed)) %>% 
  group_by(Parent_haplotype) %>% 
  mutate(Binom_p_val = binom.test.v(x = observed,
                                n = observed_sum,
                                p = Background_prob)[[3]],
         Binom_estimate = binom.test.v(x = observed,
                                n = observed_sum,
                                p = Background_prob)[[5]][[1]]) %>% 
  relocate(Analysis)



## ---------------------------------------------------------------------------------------------------------------------------------------------
impNEM_neg = impNEM %>%
  filter(logFC < 0) %>%
  distinct(GeneID, Parent_haplotype) %>% 
  count(Parent_haplotype, name = 'Resistance_genes')

p3 = ggplot(impNEM_neg, aes(Parent_haplotype, Resistance_genes, fill = Parent_haplotype))+
  geom_col(color = 'black')+
  theme(legend.position = 'none')+
  labs(y='DEGs')+
  ggtitle('NEM DEGs (-)')

p3


## ---------------------------------------------------------------------------------------------------------------------------------------------
impNEM_neg_test = impNEM_neg %>% 
  mutate(observed_sum = sum(Resistance_genes)) %>% 
  rename(observed = Resistance_genes) %>% 
  left_join(expressed, by = 'Parent_haplotype') %>% 
  left_join(ASEsums, by = 'Parent_haplotype') %>% 
  mutate(Analysis = 'NEM_neg',
         Expressed_gene_count_sum = sum(Expressed_gene_count),
         DEG_gene_count_sum = sum(DEG_gene_count),
         Background_prob = DEG_gene_count/DEG_gene_count_sum,
         expected_observed = Background_prob*observed_sum,
         log_FE = log2(observed/expected_observed)) %>% 
  group_by(Parent_haplotype) %>% 
  mutate(Binom_p_val = binom.test.v(x = observed,
                                n = observed_sum,
                                p = Background_prob)[[3]],
         Binom_estimate = binom.test.v(x = observed,
                                n = observed_sum,
                                p = Background_prob)[[5]][[1]]) %>% 
  relocate(Analysis)


View(impNEM_neg_test)


## ---------------------------------------------------------------------------------------------------------------------------------------------
t = impNEM_pos %>% 
  left_join(impNEM_neg) %>%
  arrange(desc(Susceptibility_genes)) %>%
  column_to_rownames(var = 'Parent_haplotype') %>% 
  t()

t


## ---------------------------------------------------------------------------------------------------------------------------------------------
sum(t)


## ---------------------------------------------------------------------------------------------------------------------------------------------
mosaicplot(t,
           color = T)


## ---------------------------------------------------------------------------------------------------------------------------------------------
fisher_res = fisher.test(t)

fisher_res



## ---------------------------------------------------------------------------------------------------------------------------------------------
uniq_impNEM = impNEM %>% 
  distinct(GeneID, Parent_haplotype, Gene_type) %>% 
  mutate(Pathogen = 'P. vulnus',
         Stats = paste('Fisher p-value: ',
                     round(fisher_res$p.val, 3),
                     '\n', 'Odds ratio: ',
                     round(fisher_res$estimate, 3), sep = ''))

C = ggplot(uniq_impNEM)+
  geom_mosaic(aes(x = product(Gene_type, Parent_haplotype), fill = Parent_haplotype), color = 'black')+
  labs(title = substitute(paste(italic('P. vulnus'))),
       tag = paste('Fisher p-value: ',
                     round(fisher_res$p.val, 3),
                     ', ', 'Odds ratio: ',
                     round(fisher_res$estimate, 3), sep = ''))+
  theme(legend.position = 'none', plot.tag.position = c(0.5,0.98))

C


## ---------------------------------------------------------------------------------------------------------------------------------------------
a = ggarrange(p2, p3, nrow = 2)

a


## ---------------------------------------------------------------------------------------------------------------------------------------------
all = list(uniq_impCG, uniq_impPHY, uniq_impNEM) %>% 
  bind_rows()

D = ggplot(all)+
  geom_mosaic(aes(x = product(Gene_type, Parent_haplotype),
                  fill = Parent_haplotype),
              offset = 0.03,
              color = 'black')+
  theme_grey(base_size = 11)+
    theme(legend.position = 'none',
          strip.text = element_text(face = 'italic'))+
  facet_wrap(~Pathogen+Stats)

D

sjPlot::save_plot('Final_Figs_and_Draft_Manuscript_SCRI_ROOT/CG_PHY_NEM_EBH_3.png', D, height = 12, width = 16)

