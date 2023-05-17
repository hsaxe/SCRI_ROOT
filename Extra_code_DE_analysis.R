# NEM Fit model for genes associated with NEM score
```{r}
mm_NEM = model.matrix(~RLN_2Y, data = metadata_S)

head(mm_NEM)
```

## Filtered mean-variance trend NEM
```{r}
y <- voom(d, mm_NEM, plot = T)
```

## Fitting linear models in limma with random effects NEM
```{r}
## Need to tell limma where the within class correlation is coming from
dupcor_NEM = duplicateCorrelation(y, mm_NEM, block = metadata_S$Hybrid)

## How correlated are the hybrid replicates on average?
consensus.corr.NEM = dupcor_NEM$consensus.correlation

consensus.corr.NEM

# lmFit fits a linear model using weighted least squares for each gene:
fit_NEM = lmFit(y, design = mm_NEM, block = metadata_S$Hybrid, correlation = consensus.corr.NEM) 

# Ebayes
BlockFit_NEM = eBayes(fit_NEM)
```

## Limma results NEM
```{r}
res_summaries_NEM = BlockFit_NEM %>% decideTests() %>% summary()

res_summaries_NEM
```

## Table of results NEM
```{r}
impNEM = topTable(BlockFit_NEM,
                  sort.by = "logFC",
                  p.value = 0.05,
                  adjust.method = "BH",
                  number = Inf) %>%
  rownames_to_column(var = 'GeneID') %>% 
  mutate(R = sqrt(t^2/(t^2 + 40)), AveExpr = 2^AveExpr)

dim(impNEM)
## adding Hybrid as a blocking variable reduced DEGs by several thousand

head(impNEM)
```

## PCA of NEM DEGs
```{r}
ids = impNEM %>% 
  select(GeneID)

PCA_NEM = cpmd %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'GeneID') %>% 
  right_join(ids) %>% 
  column_to_rownames(var = 'GeneID')

pca_plot_NEM = plot_pca(PCA_NEM, metadata,
                        join_by_name = 'Sample',
                        plotting_factors_in = 'col_names',
                        plotting_factors_name = Group, 
                        x = 'PC1',
                        y = 'PC2',
                        scale = T, 
                        center = T, 
                        color = 'RLN_2Y',
                        fill = 'RLN_2Y',
                        plot_type = '2D')

pca_plot_NEM$plot

```

## Merge annotation with results NEM
```{r}
impNEM = impNEM %>%
  left_join(annotation_combined, by = "GeneID")

head(impNEM)

fwrite(impNEM, 'DGEresults/Limma_results_table_NEM.csv')

impNEM = fread('DGEresults/Limma_results_table_NEM.csv')
```

## Read in NEM results
```{r}
impNEM = fread('DGEresults/Limma_results_table_NEM.csv')
```

## Volcano plot
```{r}
labs = impNEM %>% 
  distinct(logFC, adj.P.Val, name) %>% 
  slice_min(order_by = logFC, n = 5) %>% 
  rbind(impNEM %>% 
          distinct(logFC, adj.P.Val, name) %>% 
          slice_max(order_by = logFC, n = 5))


ggplot(impNEM %>% 
         distinct(logFC, adj.P.Val, name), aes(logFC, log10(adj.P.Val)*-1, color = logFC))+
  geom_point()+
  ggrepel::geom_label_repel(data = labs, aes(label = name), size = 3, color = 'black', box.padding = 0.4, label.padding = 0.1, max.overlaps = Inf)+
  geom_hline(yintercept = log10(0.05)*-1, linetype = 'dashed', color = 'red')+
  geom_text(aes(min(logFC), log10(0.05)*-1), label = 'FDR 0.05', vjust = -1)+
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black')+
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black')+
  # lims(y = c(-1, max(log10(plot_dat$fdr)*-1)))+
  scale_color_gradient2(low = 'blue', high = 'red')
```


## Positive NEM DEGs for GO
```{r}
NEM_GO_pos_Jm = impNEM %>%
  filter(logFC > 0) %>%
  distinct(GeneID, logFC, Parent_haplotype, `Jr-GeneID`) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.microcarpa") %>%
  select(`Jr-GeneID`, logFC) %>%
  mutate(`Jr-GeneID` = paste0("LOC", `Jr-GeneID`)) %>% 
  dplyr::rename(GeneID = `Jr-GeneID`)

length(unique(NEM_GO_pos_Jm$GeneID))

NEM_GO_pos_both = impNEM %>% 
  filter(logFC > 0) %>%
  distinct(GeneID, logFC, Parent_haplotype) %>%
  drop_na() %>%
  filter(Parent_haplotype == "J.regia") %>% 
  select(GeneID, logFC) %>%
  mutate(GeneID = paste0("LOC", GeneID)) %>% 
  rbind(NEM_GO_pos_Jm)

length(unique(NEM_GO_pos_both$GeneID))

head(NEM_GO_pos_both)

fwrite(NEM_GO_pos_both, "For_Panther/NEM_DEGs_pos_GO.csv", sep = '\t')

```

## Expression from alleles is causing duplication in GO results which means data loss as GO removes duplicates. However, these are not true duplicates.


## Negative NEM DEGs for GO
```{r}
NEM_GO_neg_Jm = impNEM %>%
  filter(logFC < 0) %>%
  distinct(GeneID, logFC, Parent_haplotype, `Jr-GeneID`) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.microcarpa") %>%
  select(`Jr-GeneID`, logFC) %>%
  mutate(`Jr-GeneID` = paste0("LOC", `Jr-GeneID`)) %>% 
  dplyr::rename(GeneID = `Jr-GeneID`) 


length(unique(NEM_GO_neg_Jm$GeneID))

NEM_GO_neg_both = impNEM %>% 
  filter(logFC < 0) %>%
  distinct(GeneID, logFC, Parent_haplotype) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.regia") %>% 
  select(GeneID, logFC) %>%
  mutate(GeneID = paste0("LOC", GeneID)) %>% 
  rbind(NEM_GO_neg_Jm) 

length(unique(NEM_GO_neg_both$GeneID))

head(NEM_GO_neg_both)

fwrite(NEM_GO_neg_both, "For_Panther/NEM_DEGs_neg_GO.csv", sep = '\t')
```

## Annotate all NEM DEGs with GO
```{r}
NEM_GO_ALL = NEM_GO_pos_both %>% 
  rbind(NEM_GO_neg_both) %>% 
  distinct()

fwrite(NEM_GO_ALL, "DGEresults/NEM_DEGs_all.csv", sep = '\t')
```


## Plot PCA to look in R
```{r}
pca2 = plot_pca(cpmd2, metadata_S2,
                join_by_name = 'Sample',
                plotting_factors_in = 'col_names',
                # Using 'Group' here becuase 'Hybrid' is already in metadata. Will cause error if I used 'Hybrid'
                plotting_factors_name = Group, 
                x = 'PC1',
                y = 'PC2',
                scale = T, 
                center = T, 
                color = 'RLN_2_3Y',
                fill = 'RLN_2_3Y',
                plot_type = '2D')
```

```{r}
pca2$plot
```

```{r}
c2 = plot_pca(cpmd2, metadata,
              # Having to use 'Hybrid' here to get it into the plotting data because the function drops everything but the join_by_name when summarizing for scatter
              join_by_name = 'Hybrid',
              plotting_factors_in = 'col_names',
              plotting_factors_name = Hybrid, 
              x = 'PC5',
              y = 'RLN_2_3Y',
              scale = T, 
              center = T, 
              color = 'Hybrid',
              plot_type = 'scatter',
              summarise_for_scatter = T)

c2$plot
```

# 2_3Y NEM Fit model for genes associated with NEM score
```{r}
mm_NEM_2_3Y = model.matrix(~RLN_2_3Y, 
                           data = metadata_S2)

head(mm_NEM_2_3Y)
```

## Filtered mean-variance trend NEM
```{r}
y_NEM_2_3Y <- voom(d2,
                   mm_NEM_2_3Y,
                   plot = T)
```

## Fitting linear models in limma with random effects NEM_2_3Y
```{r}
## Need to tell limma where the within class correlation is coming from
dupcor_NEM_2_3Y = duplicateCorrelation(y_NEM_2_3Y,
                                       mm_NEM_2_3Y, 
                                       block = metadata_S2$Hybrid)

## How correlated are the hybrid replicates on average?
consensus.corr.NEM_2_3Y = dupcor_NEM_2_3Y$consensus.correlation

consensus.corr.NEM_2_3Y

# lmFit fits a linear model using weighted least squares for each gene:
fit_NEM_2_3Y = lmFit(y_NEM_2_3Y,
                     design = mm_NEM_2_3Y,
                     block = metadata_S2$Hybrid,
                     correlation = consensus.corr.NEM_2_3Y) 

# Ebayes
BlockFit_NEM_2_3Y = eBayes(fit_NEM_2_3Y)
```

## Limma results NEM_2_3Y
```{r}
res_summaries_NEM_2_3Y = BlockFit_NEM_2_3Y %>% decideTests() %>% summary()

res_summaries_NEM_2_3Y
```

## Table of results NEM_2_3Y
```{r}
impNEM_2_3Y = topTable(BlockFit_NEM_2_3Y,
                       sort.by = "logFC",
                       p.value = 0.05,
                       adjust.method = "BH",
                       number = Inf) %>%
  rownames_to_column(var = 'GeneID') %>% 
  mutate(R = sqrt(t^2/(t^2 + 40)), AveExpr = 2^AveExpr)

dim(impNEM_2_3Y)
## adding Hybrid as a blocking variable reduced DEGs by several thousand

head(impNEM_2_3Y)
```

## PCA of NEM_2_3Y DEGs
```{r}
ids = impNEM_2_3Y %>% 
  select(GeneID)

PCA_NEM_2_3Y = cpmd2 %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'GeneID') %>% 
  right_join(ids) %>% 
  column_to_rownames(var = 'GeneID')

pca_plot_NEM_2_3Y = plot_pca(PCA_NEM_2_3Y, metadata,
                             join_by_name = 'Sample',
                             plotting_factors_in = 'col_names',
                             plotting_factors_name = Group, 
                             x = 'PC1',
                             y = 'PC2',
                             scale = T, 
                             center = T, 
                             color = 'RLN_3Y',
                             fill = 'RLN_3Y',
                             plot_type = '2D')

pca_plot_NEM_2_3Y$plot

```

## Merge annotation with results NEM_2_3Y
```{r}
impNEM_2_3Y = impNEM_2_3Y %>%
  left_join(annotation_combined, by = "GeneID")

head(impNEM_2_3Y)

fwrite(impNEM_2_3Y, 'DGEresults/Limma_results_table_NEM_2_3Y.csv')

impNEM_2_3Y = fread('DGEresults/Limma_results_table_NEM_2_3Y.csv')
```

## Read in NEM_2_3Y results
```{r}
impNEM_2_3Y = fread('DGEresults/Limma_results_table_NEM_2_3Y.csv')
```

## Volcano plot
```{r}
labs = impNEM_2_3Y %>% 
  distinct(logFC, adj.P.Val, name) %>% 
  slice_min(order_by = logFC, n = 5) %>% 
  rbind(impNEM_2_3Y %>% 
          distinct(logFC, adj.P.Val, name) %>% 
          slice_max(order_by = logFC, n = 5))


ggplot(impNEM_2_3Y %>% 
         distinct(logFC, adj.P.Val, name), aes(logFC, log10(adj.P.Val)*-1, color = logFC))+
  geom_point()+
  ggrepel::geom_label_repel(data = labs, aes(label = name), size = 3, color = 'black', box.padding = 0.4, label.padding = 0.1, max.overlaps = Inf)+
  geom_hline(yintercept = log10(0.05)*-1, linetype = 'dashed', color = 'red')+
  geom_text(aes(min(logFC), log10(0.05)*-1), label = 'FDR 0.05', vjust = -1)+
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black')+
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black')+
  # lims(y = c(-1, max(log10(plot_dat$fdr)*-1)))+
  scale_color_gradient2(low = 'blue', high = 'red')
```


## Positive NEM_2_3Y DEGs for GO
```{r}
NEM_2_3Y_GO_pos_Jm = impNEM_2_3Y %>%
  filter(logFC > 0) %>%
  distinct(GeneID, logFC, Parent_haplotype, `Jr-GeneID`) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.microcarpa") %>%
  select(`Jr-GeneID`, logFC) %>%
  mutate(`Jr-GeneID` = paste0("LOC", `Jr-GeneID`)) %>% 
  dplyr::rename(GeneID = `Jr-GeneID`)

length(unique(NEM_2_3Y_GO_pos_Jm$GeneID))

NEM_2_3Y_GO_pos_both = impNEM_2_3Y %>% 
  filter(logFC > 0) %>%
  distinct(GeneID, logFC, Parent_haplotype) %>%
  drop_na() %>%
  filter(Parent_haplotype == "J.regia") %>% 
  select(GeneID, logFC) %>%
  mutate(GeneID = paste0("LOC", GeneID)) %>% 
  rbind(NEM_2_3Y_GO_pos_Jm)

length(unique(NEM_2_3Y_GO_pos_both$GeneID))

head(NEM_2_3Y_GO_pos_both)

fwrite(NEM_2_3Y_GO_pos_both, "For_Panther/NEM_2_3Y_DEGs_pos_GO.csv", sep = '\t')

```

## Expression from alleles is causing duplication in GO results which means data loss as GO removes duplicates. However, these are not true duplicates.


## Negative NEM_2_3Y DEGs for GO
```{r}
NEM_2_3Y_GO_neg_Jm = impNEM_2_3Y %>%
  filter(logFC < 0) %>%
  distinct(GeneID, logFC, Parent_haplotype, `Jr-GeneID`) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.microcarpa") %>%
  select(`Jr-GeneID`, logFC) %>%
  mutate(`Jr-GeneID` = paste0("LOC", `Jr-GeneID`)) %>% 
  dplyr::rename(GeneID = `Jr-GeneID`) 


length(unique(NEM_2_3Y_GO_neg_Jm$GeneID))

NEM_2_3Y_GO_neg_both = impNEM_2_3Y %>% 
  filter(logFC < 0) %>%
  distinct(GeneID, logFC, Parent_haplotype) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.regia") %>% 
  select(GeneID, logFC) %>%
  mutate(GeneID = paste0("LOC", GeneID)) %>% 
  rbind(NEM_2_3Y_GO_neg_Jm) 

length(unique(NEM_2_3Y_GO_neg_both$GeneID))

head(NEM_2_3Y_GO_neg_both)

fwrite(NEM_2_3Y_GO_neg_both, "For_Panther/NEM_2_3Y_DEGs_neg_GO.csv", sep = '\t')
```

## Annotate all NEM_2_3Y DEGs with GO
```{r}
NEM_2_3Y_GO_ALL = NEM_2_3Y_GO_pos_both %>% 
  rbind(NEM_2_3Y_GO_neg_both) %>% 
  distinct()

fwrite(NEM_2_3Y_GO_ALL, "DGEresults/NEM_2_3Y_DEGs_all.csv", sep = '\t')
```