## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----echo=TRUE--------------------------------------------------------------------------------------------------------------------------------
pacman::p_load(base,
               edgeR,
               tidytable,
               data.table,
               dplyr,
               stats,
               ggpubr,
               tidyr,
               ggplot2,
               gridExtra,
               limma,
               OmicsAnalyst,
               sjPlot,
               tibble,
               tidytext,
               utils)


## ---------------------------------------------------------------------------------------------------------------------------------------------
dat = fread("C:/Users/hsaxe/OneDrive/Documents/ALAB/Transcriptome_data/Root/SCRI_ROOT_RNAseq_counts_combined_genomes.txt")

dat$GeneID = gsub("LOC", "", dat$GeneID)


head(dat)


## ---------------------------------------------------------------------------------------------------------------------------------------------
metadata = fread("C:/Users/hsaxe/OneDrive/Documents/ALAB/Transcriptome_data/Root/R/Phenotyping/SCRI/LongList3_2Y.csv", stringsAsFactors = T)

head(metadata)


## ---------------------------------------------------------------------------------------------------------------------------------------------
metadata = data.frame(Sample = colnames(dat)[colnames(dat) != 'GeneID']) %>%
  mutate(Hybrid = as.factor(gsub("\\-\\d$", "", Sample))) %>%
  left_join(metadata, by = c('Hybrid' = 'CAL:_Wip_ID'))
  
row.names(metadata) = metadata$Sample

metadata = metadata %>% 
  mutate(PHY_Avg = rowMeans(select(., matches('PCLR|PRLR'))))

head(metadata)

fwrite(metadata, 'SCRI_ROOT_RNAseq_metadata.csv')


## ---------------------------------------------------------------------------------------------------------------------------------------------
metaLong = metadata %>% 
  select(!Sample) %>% 
  select(Hybrid, CG_Avg., PHY_Avg, TwoY_RLN) %>% 
  distinct() %>% 
  pivot_longer(where(is.numeric), names_to = 'Trait')

p =ggplot(metaLong, aes(reorder_within(Hybrid, value, Trait), value, fill = Hybrid))+
  geom_col()+
  scale_x_reordered()+
  theme(axis.text.x = element_text(angle = 30))+
  facet_wrap(~Trait, scales = 'free', ncol = 4)
p

save_plot('DGEresults/metadata_plot.png', p, height = 10, width = 20)


## ---------------------------------------------------------------------------------------------------------------------------------------------
metadata_S = metadata %>%
  mutate_if(is.numeric, scale)

head(metadata_S)


## ---------------------------------------------------------------------------------------------------------------------------------------------
dat1 = dat %>%
  column_to_rownames(var =  "GeneID") %>%
  as.matrix()

## Do colnames in data match rownames in metadata? If they don't, use match(x,y) produces the order of y required to match the order of x

all(colnames(dat1) == rownames(metadata_S))
# Names match

# If they didn't match, use below code
## What order do rows of metadat need to be in to match colnames of dat1?
# match(colnames(dat1), rownames(metadata))

## Reset rownames
# metadata = metadata[match(colnames(dat1), rownames(metadata)),]

# all(colnames(dat1) == rownames(metadata))
# now they match


## ---------------------------------------------------------------------------------------------------------------------------------------------
annotation_Jm = fread("C:/Users/hsaxe/OneDrive/Documents/ALAB/Genome_info/Genomic_Annotation_2/Jm_x_Jr/Jm_x_Jr_Genomic_annotation.csv")

## Extract everything but class mRNA and other isoforms. This reduces duplication in the data
annotation_Jm = annotation_Jm %>%
  filter(feature != "mRNA", !grepl('\\sX[2-9]$|\\sX1[0-9]$', name)) %>%
  mutate(GeneID = as.character(GeneID)) %>% 
  mutate(Parent_haplotype = "J.microcarpa")

head(annotation_Jm)



## ---------------------------------------------------------------------------------------------------------------------------------------------
annotation_Jr = fread("C:/Users/hsaxe/OneDrive/Documents/ALAB/Genome_info/Genomic_Annotation_2/Jr/Jr_Genomic_annotation.csv")

## Extract everything but class mRNA and other isoforms. This reduces duplication in the data
annotation_Jr = annotation_Jr %>%
  filter(feature != "mRNA", !grepl('\\sX[2-9]$|\\sX1[0-9]$', name)) %>%
  mutate(GeneID = as.character(GeneID)) %>% 
  mutate(Parent_haplotype = "J.regia")

head(annotation_Jr)

fwrite(annotation_Jr %>% 
         select(GeneID) %>% 
         mutate(GeneID = paste0('LOC', GeneID)),
       'Annotate_annotation_with_GO.csv')



## ---------------------------------------------------------------------------------------------------------------------------------------------
annotation_combined = annotation_Jm %>%
  rbind(annotation_Jr, fill = T)

fwrite(annotation_combined, 'DGEresults/annotation_combined.csv')

BP_anno = annotation_Jm %>% 
  distinct(`Jr-GeneID`) %>% 
  rbind(distinct(annotation_Jr, GeneID), use.names = F) %>% 
  distinct(`Jr-GeneID`) %>% 
  mutate(`Jr-GeneID` = paste0('LOC', `Jr-GeneID`))

fwrite(BP_anno, 'GeneIDs_for_annotation_with_GO.csv')


## ---------------------------------------------------------------------------------------------------------------------------------------------
library(edgeR)

dds = DGEList(dat1)

dim(dds$counts)

## Calculate library normalization factors (does not do anything to data)
dds = calcNormFactors(dds)

## These are the size factors (normalization factors) for each sample
# dds$samples


## ---------------------------------------------------------------------------------------------------------------------------------------------
d = expression_filter(dds, DGEList = T, FilterFUN = max, FilterThreshold = 75)

## CPM normalized counts of all data
cpm = cpm(dds, prior.count = 2, log = F) 

cpm = cpm %>% 
  data.frame() %>% 
  rownames_to_column(var = 'GeneID') %>%  
  rename_with(~ gsub('X', '', gsub('\\.', '-', .x)))

fwrite(cpm, 'DGEresults/cpm_plotting_data.csv')

## CPM normalized counts of filtered data
cpmd = cpm(d)

fwrite(cpmd %>% data.frame() %>% rownames_to_column(var = 'GeneID'), 'DGEresults/DGE_CPM_data.csv')


## ---------------------------------------------------------------------------------------------------------------------------------------------
# d2 = expression_filter(dds, DGEList = T, FilterFUN = mean, FilterThreshold = 20)
# 
# fwrite(as.data.frame(d2$counts), 'DGEresults/CPM_20_filtered.csv', row.names = T)


## ---------------------------------------------------------------------------------------------------------------------------------------------
## MDS visualization
# plotMDS(d, col = as.numeric(metadata$Hybrid), cex=1)


## ---------------------------------------------------------------------------------------------------------------------------------------------
pca = plot_pca(cpmd, metadata,
               join_by_name = 'Sample',
               plotting_factors_in = 'col_names',
               # Using 'Group' here becuase 'Hybrid' is already in metadata. Will cause error if I used 'Hybrid'
               plotting_factors_name = Group, 
               x = 'PC1',
               y = 'PC2',
               scale = T, 
               center = T, 
               color = 'Hybrid',
               fill = 'Hybrid',
               plot_type = '2D')


## ---------------------------------------------------------------------------------------------------------------------------------------------
pca$plot


## ---------------------------------------------------------------------------------------------------------------------------------------------
# a = plot_pca(cpmd, metadata, sample_colname = 'Sample', sample_names_in = 'col_names', x = 'PC1', y = 'PC2', scale = T, center = T, color = 'CG_Avg.', fill = 'CG_Avg.', plot_type = '2D', group_by = 'Hybrid')
# 
# b = plot_pca(cpmd, metadata, sample_colname = 'Sample', sample_names_in = 'col_names', x = 'PC1', y = 'PC2', scale = T, center = T, color = 'Cinn.PRLR', fill = 'Cinn.PRLR', plot_type = '2D', group_by = 'Hybrid')
#           
# c = plot_pca(cpmd, metadata, sample_colname = 'Sample', sample_names_in = 'col_names', x = 'PC1', y = 'PC2', scale = T, center = T, color = 'TwoY_RLN', fill = 'TwoY_RLN', plot_type = '2D', group_by = 'Hybrid')
# 
# 
# save_plot('PCA_Hybrids.png', ggarrange(a, b, c, labels = c('A', 'B', 'C'), ncol = 3), base_height = 3, base_width = 12)



## ---------------------------------------------------------------------------------------------------------------------------------------------
# a = plot_pca(cpmd, metadata,
#                join_by_name = 'Sample',
#                plotting_factors_in = 'col_names',
#                plotting_factors_name = Hybrid, 
#                x = 'Hybrid',
#                y = 'PC2',
#                scale = T, 
#                center = T, 
#                color = 'Hybrid',
#                fill = 'CG_Avg.',
#                plot_type = 'boxplot')
# 
# b = plot_pca(cpmd, metadata,
#                join_by_name = 'Sample',
#                plotting_factors_in = 'col_names',
#                plotting_factors_name = Hybrid, 
#                x = 'Hybrid',
#                y = 'PC5',
#                scale = T, 
#                center = T, 
#                color = 'Hybrid',
#                fill = 'PHY_Avg',
#                plot_type = 'boxplot')
# 
# c = plot_pca(cpmd, metadata,
#                join_by_name = 'Sample',
#                plotting_factors_in = 'col_names',
#                plotting_factors_name = Hybrid, 
#                x = 'Hybrid',
#                y = 'PC1',
#                scale = T, 
#                center = T, 
#                color = 'Hybrid',
#                fill = 'TwoY_RLN',
#                plot_type = 'boxplot')



## ----fig.height=10, fig.width=7---------------------------------------------------------------------------------------------------------------
## Available plotting data: 
# [1] "Sample"      "Hybrid"      "CG_Avg."     "CG_Dec."     "Cinn.PCLR"   "Cinn.PRLR"  
#  [7] "Pini.PCLR"   "Pini.PRLR"   "PHY_Dec."    "TwoY_length" "TwoY_RLN"    "NEM_Dec."   # [13] "Unity"       "PHY_Avg" 

a = plot_pca(cpmd, metadata,
             # Having to use 'Hybrid' here to get it into the plotting data because the function drops everything but the join_by_name when summarizing for scatter
             join_by_name = 'Hybrid',
             plotting_factors_in = 'col_names',
             plotting_factors_name = Hybrid,
             x = 'PC2', y = 'CG_Avg.',
             scale = T,
             center = T,
             color = 'Hybrid',
             plot_type = 'scatter',
             summarise_for_scatter = T)

b = plot_pca(cpmd, metadata,
             # Having to use 'Hybrid' here to get it into the plotting data because the function drops everything but the join_by_name when summarizing for scatter
             join_by_name = 'Hybrid',
             plotting_factors_in = 'col_names',
             plotting_factors_name = Hybrid,
             x = 'PC5',
             y = 'PHY_Avg',
             scale = T,
             center = T,
             color = 'Hybrid',
             plot_type = 'scatter',
             summarise_for_scatter = T)

c = plot_pca(cpmd, metadata,
             # Having to use 'Hybrid' here to get it into the plotting data because the function drops everything but the join_by_name when summarizing for scatter
             join_by_name = 'Hybrid',
             plotting_factors_in = 'col_names',
             plotting_factors_name = Hybrid, 
             x = 'PC1',
             y = 'TwoY_RLN',
             scale = T, 
             center = T, 
             color = 'Hybrid',
             plot_type = 'scatter',
             summarise_for_scatter = T)

# d2 = plot_pca(cpmd, metadata,
#              join_by_name = 'Sample',
#              plotting_factors_in = 'col_names',
#              plotting_factors_name = Hybrid, 
#              x = 'PC5',
#              y = 'TwoY_length',
#              scale = T, 
#              center = T, 
#              color = 'Hybrid',
#              plot_type = 'scatter',
#              summarise_for_scatter = T)

arr = ggarrange(a$plot,
                b$plot,
                c$plot,
                # d2$plot + ggtitle('Trait: Length'),
                labels = c('A)', 'B)', 'C)'), nrow = 3, ncol = 1)

arr

sjPlot::save_plot('DGEresults/PCA_scatter_CG_PH_NEM.png', arr, height = 24, width = 14)


## ---------------------------------------------------------------------------------------------------------------------------------------------
# library(pheatmap)
# 
# cor = cor(cpmd)
# 
# pheatmap(cor, annotation = select(metadata,Hybrid))


## ---------------------------------------------------------------------------------------------------------------------------------------------
mm = model.matrix(~CG_Avg., data = metadata_S)

head(mm)


## ---------------------------------------------------------------------------------------------------------------------------------------------
y <- voom(d, mm, plot = T)


## ---------------------------------------------------------------------------------------------------------------------------------------------
# tmp <- voom(dds, mm, plot = T)


## ---------------------------------------------------------------------------------------------------------------------------------------------
## Need to tell limma where the within class correlation is coming from
dupcor_CG = duplicateCorrelation(y, mm, block = metadata$Hybrid)

## How correlated are the hybrid replicates on average?
consensus.corr.CG = dupcor_CG$consensus.correlation

consensus.corr.CG

# lmFit fits a linear model using weighted least squares for each gene:
fit = lmFit(y, design = mm, block = metadata$Hybrid, correlation = consensus.corr.CG) 


## ---------------------------------------------------------------------------------------------------------------------------------------------
BlockFit_CG = eBayes(fit)


## ---------------------------------------------------------------------------------------------------------------------------------------------
res_summaries_CG = BlockFit_CG %>% decideTests() %>% summary()

res_summaries_CG


## ---------------------------------------------------------------------------------------------------------------------------------------------
impCG = limma::topTable(BlockFit_CG, 
                        sort.by = "logFC",
                        p.value = 0.05, 
                        adjust.method = "BH",
                        number = Inf) %>%
  rownames_to_column(var = 'GeneID') %>% 
  mutate(R = sqrt(t^2/(t^2 + 40)), AveExpr = 2^AveExpr)

fwrite(impCG, 'DGEresults/unAnnotated_results_CG.csv')

rawCG = fread('DGEresults/unAnnotated_results_CG.csv')

dim(impCG)
## adding Hybrid as a blocking variable reduced DEGs by several thousand


## ---------------------------------------------------------------------------------------------------------------------------------------------
ids = impCG %>% 
  select(GeneID)

PCA_CG = cpmd %>% 
  data.frame() %>% 
  rownames_to_column(var = 'GeneID') %>% 
  right_join(ids) %>% 
  column_to_rownames(var = 'GeneID')

pca_plot_CG = plot_pca(PCA_CG, metadata,
               join_by_name = 'Sample',
               plotting_factors_in = 'col_names',
               plotting_factors_name = Group, 
               x = 'PC1',
               y = 'PC2',
               scale = T, 
               center = T, 
               color = 'CG_Avg.',
               fill = 'CG_Avg.',
               plot_type = '2D')

pca_plot_CG$plot



## ---------------------------------------------------------------------------------------------------------------------------------------------
impCG = impCG %>%

  left_join(annotation_combined, by = "GeneID")

head(impCG)

length(unique(impCG$GeneID))

fwrite(impCG, 'DGEresults/Limma_results_table_CG.csv')


## ---------------------------------------------------------------------------------------------------------------------------------------------
impCG = fread('DGEresults/Limma_results_table_CG.csv')


## ---------------------------------------------------------------------------------------------------------------------------------------------
labs = impCG %>% 
  distinct(logFC, adj.P.Val, name) %>% 
  slice_min(order_by = logFC, n = 5) %>% 
  rbind(impCG %>% 
  distinct(logFC, adj.P.Val, name) %>% 
    slice_max(order_by = logFC, n = 5))
      

ggplot(impCG %>% 
         distinct(logFC, adj.P.Val, name), aes(logFC, log10(adj.P.Val)*-1, color = logFC))+
  geom_point()+
   ggrepel::geom_label_repel(data = labs, aes(label = name), size = 3, color = 'black', box.padding = 0.4, label.padding = 0.1, max.overlaps = Inf)+
  geom_hline(yintercept = log10(0.05)*-1, linetype = 'dashed', color = 'red')+
  geom_text(aes(min(logFC), log10(0.05)*-1), label = 'FDR 0.05', vjust = -1)+
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black')+
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black')+
  # lims(y = c(-1, max(log10(plot_dat$fdr)*-1)))+
  scale_color_gradient2(low = 'blue', high = 'red')


## ---------------------------------------------------------------------------------------------------------------------------------------------
CG_GO_pos_Jm = impCG %>%
  filter(logFC > 0) %>%
  distinct(GeneID, logFC, Parent_haplotype, `Jr-GeneID`) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.microcarpa") %>%
  select(`Jr-GeneID`, logFC) %>%
  mutate(`Jr-GeneID` = paste0("LOC", `Jr-GeneID`)) %>% 
  dplyr::rename(GeneID = `Jr-GeneID`)

length(unique(CG_GO_pos_Jm$GeneID))

CG_GO_pos_both = impCG %>% 
  filter(logFC > 0) %>%
  distinct(GeneID, logFC, Parent_haplotype) %>%
  drop_na() %>%
  filter(Parent_haplotype == "J.regia") %>% 
  select(GeneID, logFC) %>%
  mutate(GeneID = paste0("LOC", GeneID)) %>% 
  rbind(CG_GO_pos_Jm)

length(unique(CG_GO_pos_both$GeneID))

head(CG_GO_pos_both)

fwrite(CG_GO_pos_both, "DGEresults/CG_DEGs_pos_GO.csv", sep = '\t')



## ---------------------------------------------------------------------------------------------------------------------------------------------
CG_GO_neg_Jm = impCG %>%
  filter(logFC < 0) %>%
  distinct(GeneID, logFC, Parent_haplotype, `Jr-GeneID`) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.microcarpa") %>%
  select(`Jr-GeneID`, logFC) %>%
  mutate(`Jr-GeneID` = paste0("LOC", `Jr-GeneID`)) %>% 
  dplyr::rename(GeneID = `Jr-GeneID`) 
 

length(unique(CG_GO_neg_Jm$GeneID))

CG_GO_neg_both = impCG %>% 
  filter(logFC < 0) %>%
  distinct(GeneID, logFC, Parent_haplotype) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.regia") %>% 
  select(GeneID, logFC) %>%
  mutate(GeneID = paste0("LOC", GeneID)) %>% 
  rbind(CG_GO_neg_Jm) 

length(unique(CG_GO_neg_both$GeneID))

head(CG_GO_neg_both)

fwrite(CG_GO_neg_both, "DGEresults/CG_DEGs_neg_GO.csv", sep = '\t')


## ---------------------------------------------------------------------------------------------------------------------------------------------
CG_GO_ALL = CG_GO_pos_both %>% 
  rbind(CG_GO_neg_both) %>% 
  distinct()

fwrite(CG_GO_ALL, "DGEresults/CG_DEGs_all.csv", sep = '\t')


## ---------------------------------------------------------------------------------------------------------------------------------------------
mm = model.matrix(~PHY_Avg, data = metadata_S)

head(mm)


## ---------------------------------------------------------------------------------------------------------------------------------------------
y <- voom(d, mm, plot = T)


## ---------------------------------------------------------------------------------------------------------------------------------------------
## Need to tell limma where the within class correlation is coming from
dupcor_PHY = duplicateCorrelation(y, mm, block = metadata$Hybrid)

## How correlated are the hybrid replicates on average?
consensus.corr.PHY = dupcor_PHY$consensus.correlation

consensus.corr.PHY

# lmFit fits a linear model using weighted least squares for each gene:
fit = lmFit(y, design = mm, block = metadata$Hybrid, correlation = consensus.corr.PHY) 

# Ebayes
BlockFit_PHY = eBayes(fit)


## ---------------------------------------------------------------------------------------------------------------------------------------------
res_summaries_PHY = BlockFit_PHY %>% decideTests() %>% summary()

res_summaries_PHY


## ---------------------------------------------------------------------------------------------------------------------------------------------
impPHY = topTable(BlockFit_PHY, sort.by = "logFC", p.value = 0.05, adjust.method = "BH", number = Inf) %>%
    rownames_to_column(var = 'GeneID') %>% 
  mutate(R = sqrt(t^2/(t^2 + 40)), AveExpr = 2^AveExpr)

dim(impPHY)
## adding Hybrid as a blocking variable reduced DEGs by several thousand

head(impPHY)


## ---------------------------------------------------------------------------------------------------------------------------------------------
ids = impPHY %>% 
  select(GeneID)

PCA_PHY = cpmd %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'GeneID') %>% 
  right_join(ids) %>% 
  column_to_rownames(var = 'GeneID')

pca_plot_PHY = plot_pca(PCA_PHY, metadata,
               join_by_name = 'Sample',
               plotting_factors_in = 'col_names',
               plotting_factors_name = Group, 
               x = 'PC1',
               y = 'PC2',
               scale = T, 
               center = T, 
               color = 'PHY_Avg',
               fill = 'PHY_Avg',
               plot_type = '2D')

pca_plot_PHY$plot



## ---------------------------------------------------------------------------------------------------------------------------------------------
impPHY = impPHY %>%

  left_join(annotation_combined, by = "GeneID")

head(impPHY)

fwrite(impPHY, 'DGEresults/Limma_results_table_PHY.csv')


## ---------------------------------------------------------------------------------------------------------------------------------------------
impPHY = fread('DGEresults/Limma_results_table_PHY.csv')


## ---------------------------------------------------------------------------------------------------------------------------------------------
labs = impPHY %>% 
  distinct(logFC, adj.P.Val, name) %>% 
  slice_min(order_by = logFC, n = 5) %>% 
  rbind(impPHY %>% 
  distinct(logFC, adj.P.Val, name) %>% 
    slice_max(order_by = logFC, n = 5))
      

ggplot(impPHY %>% 
         distinct(logFC, adj.P.Val, name), aes(logFC, log10(adj.P.Val)*-1, color = logFC))+
  geom_point()+
   ggrepel::geom_label_repel(data = labs, aes(label = name), size = 3, color = 'black', box.padding = 0.4, label.padding = 0.1, max.overlaps = Inf)+
  geom_hline(yintercept = log10(0.05)*-1, linetype = 'dashed', color = 'red')+
  geom_text(aes(min(logFC), log10(0.05)*-1), label = 'FDR 0.05', vjust = -1)+
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black')+
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black')+
  # lims(y = c(-1, max(log10(plot_dat$fdr)*-1)))+
  scale_color_gradient2(low = 'blue', high = 'red')


## ---------------------------------------------------------------------------------------------------------------------------------------------
PHY_GO_pos_Jm = impPHY %>%
  filter(logFC > 0) %>%
  distinct(GeneID, logFC, Parent_haplotype, `Jr-GeneID`) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.microcarpa") %>%
  select(`Jr-GeneID`, logFC) %>%
  mutate(`Jr-GeneID` = paste0("LOC", `Jr-GeneID`)) %>% 
  dplyr::rename(GeneID = `Jr-GeneID`)

length(unique(PHY_GO_pos_Jm$GeneID))

PHY_GO_pos_both = impPHY %>% 
  filter(logFC > 0) %>%
  distinct(GeneID, logFC, Parent_haplotype) %>%
  drop_na() %>%
  filter(Parent_haplotype == "J.regia") %>% 
  select(GeneID, logFC) %>%
  mutate(GeneID = paste0("LOC", GeneID)) %>% 
  rbind(PHY_GO_pos_Jm)

length(unique(PHY_GO_pos_both$GeneID))

head(PHY_GO_pos_both)

fwrite(PHY_GO_pos_both, "DGEresults/PHY_DEGs_pos_GO.csv", sep = '\t')



## ---------------------------------------------------------------------------------------------------------------------------------------------
PHY_GO_neg_Jm = impPHY %>%
  filter(logFC < 0) %>%
  distinct(GeneID, logFC, Parent_haplotype, `Jr-GeneID`) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.microcarpa") %>%
  select(`Jr-GeneID`, logFC) %>%
  mutate(`Jr-GeneID` = paste0("LOC", `Jr-GeneID`)) %>% 
  dplyr::rename(GeneID = `Jr-GeneID`) 
 

length(unique(PHY_GO_neg_Jm$GeneID))

PHY_GO_neg_both = impPHY %>% 
  filter(logFC < 0) %>%
  distinct(GeneID, logFC, Parent_haplotype) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.regia") %>% 
  select(GeneID, logFC) %>%
  mutate(GeneID = paste0("LOC", GeneID)) %>% 
  rbind(PHY_GO_neg_Jm) 

length(unique(PHY_GO_neg_both$GeneID))

head(PHY_GO_neg_both)

fwrite(PHY_GO_neg_both, "DGEresults/PHY_DEGs_neg_GO.csv", sep = '\t')


## ---------------------------------------------------------------------------------------------------------------------------------------------
PHY_GO_ALL = PHY_GO_pos_both %>% 
  rbind(PHY_GO_neg_both) %>% 
  distinct()

fwrite(PHY_GO_ALL, "DGEresults/PHY_DEGs_all.csv", sep = '\t')


## ---------------------------------------------------------------------------------------------------------------------------------------------
mm = model.matrix(~Cinn.PRLR, data = metadata_S)

head(mm)


## ---------------------------------------------------------------------------------------------------------------------------------------------
y <- voom(d, mm, plot = T)


## ---------------------------------------------------------------------------------------------------------------------------------------------
## Need to tell limma where the within class correlation is coming from
dupcor_P_cinn_root = duplicateCorrelation(y, mm, block = metadata$Hybrid)

## How correlated are the hybrid replicates on average?
consensus.corr.P_cinn_root = dupcor_P_cinn_root$consensus.correlation

consensus.corr.P_cinn_root

# lmFit fits a linear model using weighted least squares for each gene:
fit = lmFit(y, design = mm, block = metadata$Hybrid, correlation = consensus.corr.P_cinn_root) 

# Ebayes
BlockFit_P_cinn_root = eBayes(fit)


## ---------------------------------------------------------------------------------------------------------------------------------------------
res_summaries_P_cinn_root = BlockFit_P_cinn_root %>% decideTests() %>% summary()

res_summaries_P_cinn_root


## ---------------------------------------------------------------------------------------------------------------------------------------------
impP_cinn_root = topTable(BlockFit_P_cinn_root, 
                          sort.by = "logFC",
                          p.value = 0.05, 
                          adjust.method = "BH",
                          number = Inf) %>%
    rownames_to_column(var = 'GeneID') %>% 
  mutate(R = sqrt(t^2/(t^2 + 40)), AveExpr = 2^AveExpr)

dim(impP_cinn_root)
## adding Hybrid as a blocking variable reduced DEGs by several thousand

head(impP_cinn_root)


## ---------------------------------------------------------------------------------------------------------------------------------------------
ids = impP_cinn_root %>% 

  select(GeneID)

PCA_P_cinn_root = cpmd %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'GeneID') %>% 
  right_join(ids) %>% 
  column_to_rownames(var = 'GeneID')

pca_plot_P_cinn_root = plot_pca(PCA_P_cinn_root, metadata,
               join_by_name = 'Sample',
               plotting_factors_in = 'col_names',
               plotting_factors_name = Group, 
               x = 'PC1',
               y = 'PC2',
               scale = T, 
               center = T, 
               color = 'Cinn.PRLR',
               fill = 'Cinn.PRLR',
               plot_type = '2D')

pca_plot_P_cinn_root$plot



## ---------------------------------------------------------------------------------------------------------------------------------------------
impP_cinn_root = impP_cinn_root %>%

  left_join(annotation_combined, by = "GeneID")

head(impP_cinn_root)

fwrite(impP_cinn_root, 'DGEresults/Limma_results_table_P_cinn_root.csv')


## ---------------------------------------------------------------------------------------------------------------------------------------------
impP_cinn_root = fread('DGEresults/Limma_results_table_P_cinn_root.csv')


## ---------------------------------------------------------------------------------------------------------------------------------------------
labs = impP_cinn_root %>% 
  distinct(logFC, adj.P.Val, name) %>% 
  slice_min(order_by = logFC, n = 5) %>% 
  rbind(impP_cinn_root %>% 
  distinct(logFC, adj.P.Val, name) %>% 
    slice_max(order_by = logFC, n = 5))
      

ggplot(impP_cinn_root %>% 
         distinct(logFC, adj.P.Val, name), aes(logFC, log10(adj.P.Val)*-1, color = logFC))+
  geom_point()+
   ggrepel::geom_label_repel(data = labs, aes(label = name), size = 3, color = 'black', box.padding = 0.4, label.padding = 0.1, max.overlaps = Inf)+
  geom_hline(yintercept = log10(0.05)*-1, linetype = 'dashed', color = 'red')+
  geom_text(aes(min(logFC), log10(0.05)*-1), label = 'FDR 0.05', vjust = -1)+
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black')+
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black')+
  # lims(y = c(-1, max(log10(plot_dat$fdr)*-1)))+
  scale_color_gradient2(low = 'blue', high = 'red')


## ---------------------------------------------------------------------------------------------------------------------------------------------
P_cinn_root_GO_pos_Jm = impP_cinn_root %>%
  filter(logFC > 0) %>%
  distinct(GeneID, logFC, Parent_haplotype, `Jr-GeneID`) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.microcarpa") %>%
  select(`Jr-GeneID`, logFC) %>%
  mutate(`Jr-GeneID` = paste0("LOC", `Jr-GeneID`)) %>% 
  dplyr::rename(GeneID = `Jr-GeneID`)

length(unique(P_cinn_root_GO_pos_Jm$GeneID))

P_cinn_root_GO_pos_both = impP_cinn_root %>% 
  filter(logFC > 0) %>%
  distinct(GeneID, logFC, Parent_haplotype) %>%
  drop_na() %>%
  filter(Parent_haplotype == "J.regia") %>% 
  select(GeneID, logFC) %>%
  mutate(GeneID = paste0("LOC", GeneID)) %>% 
  rbind(P_cinn_root_GO_pos_Jm)

length(unique(P_cinn_root_GO_pos_both$GeneID))

head(P_cinn_root_GO_pos_both)

fwrite(P_cinn_root_GO_pos_both, "DGEresults/P_cinn_root_DEGs_pos_GO.csv", sep = '\t')



## ---------------------------------------------------------------------------------------------------------------------------------------------
P_cinn_root_GO_neg_Jm = impP_cinn_root %>%
  filter(logFC < 0) %>%
  distinct(GeneID, logFC, Parent_haplotype, `Jr-GeneID`) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.microcarpa") %>%
  select(`Jr-GeneID`, logFC) %>%
  mutate(`Jr-GeneID` = paste0("LOC", `Jr-GeneID`)) %>% 
  dplyr::rename(GeneID = `Jr-GeneID`) 
 

length(unique(P_cinn_root_GO_neg_Jm$GeneID))

P_cinn_root_GO_neg_both = impP_cinn_root %>% 
  filter(logFC < 0) %>%
  distinct(GeneID, logFC, Parent_haplotype) %>%
  drop_na() %>% 
  filter(Parent_haplotype == "J.regia") %>% 
  select(GeneID, logFC) %>%
  mutate(GeneID = paste0("LOC", GeneID)) %>% 
  rbind(P_cinn_root_GO_neg_Jm) 

length(unique(P_cinn_root_GO_neg_both$GeneID))

head(P_cinn_root_GO_neg_both)

fwrite(P_cinn_root_GO_neg_both, "DGEresults/P_cinn_root_DEGs_neg_GO.csv", sep = '\t')


## ---------------------------------------------------------------------------------------------------------------------------------------------
P_cinn_root_GO_ALL = P_cinn_root_GO_pos_both %>% 
  rbind(P_cinn_root_GO_neg_both) %>% 
  distinct()

fwrite(P_cinn_root_GO_ALL, "DGEresults/P_cinn_root_DEGs_all.csv", sep = '\t')


## ---------------------------------------------------------------------------------------------------------------------------------------------
mm = model.matrix(~TwoY_RLN, data = metadata_S)

head(mm)


## ---------------------------------------------------------------------------------------------------------------------------------------------
y <- voom(d, mm, plot = T)


## ---------------------------------------------------------------------------------------------------------------------------------------------
## Need to tell limma where the within class correlation is coming from
dupcor_NEM = duplicateCorrelation(y, mm, block = metadata$Hybrid)

## How correlated are the hybrid replicates on average?
consensus.corr.NEM = dupcor_NEM$consensus.correlation

consensus.corr.NEM

# lmFit fits a linear model using weighted least squares for each gene:
fit = lmFit(y, design = mm, block = metadata$Hybrid, correlation = consensus.corr.NEM) 

# Ebayes
BlockFit_NEM = eBayes(fit)


## ---------------------------------------------------------------------------------------------------------------------------------------------
res_summaries_NEM = BlockFit_NEM %>% decideTests() %>% summary()

res_summaries_NEM


## ---------------------------------------------------------------------------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------------------------------------------------------------------------
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
               color = 'TwoY_RLN',
               fill = 'TwoY_RLN',
               plot_type = '2D')

pca_plot_NEM$plot



## ---------------------------------------------------------------------------------------------------------------------------------------------
impNEM = impNEM %>%

  left_join(annotation_combined, by = "GeneID")

head(impNEM)

fwrite(impNEM, 'DGEresults/Limma_results_table_NEM.csv')

impNEM = fread('DGEresults/Limma_results_table_NEM.csv')


## ---------------------------------------------------------------------------------------------------------------------------------------------
impNEM = fread('DGEresults/Limma_results_table_NEM.csv')


## ---------------------------------------------------------------------------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------------------------------------------------------------------------
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

fwrite(NEM_GO_pos_both, "DGEresults/NEM_DEGs_pos_GO.csv", sep = '\t')



## ---------------------------------------------------------------------------------------------------------------------------------------------
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

fwrite(NEM_GO_neg_both, "DGEresults/NEM_DEGs_neg_GO.csv", sep = '\t')


## ---------------------------------------------------------------------------------------------------------------------------------------------
NEM_GO_ALL = NEM_GO_pos_both %>% 
  rbind(NEM_GO_neg_both) %>% 
  distinct()

fwrite(NEM_GO_ALL, "DGEresults/NEM_DEGs_all.csv", sep = '\t')


## ---------------------------------------------------------------------------------------------------------------------------------------------
save_plot('DGEresults/PCA_DEGs.png',
ggarrange(pca_plot_CG$plot, pca_plot_PHY$plot, pca_plot_NEM$plot, labels = c('CG)', 'PHY)', 'NEM)')),
height = 20, width = 24
)


## ---------------------------------------------------------------------------------------------------------------------------------------------
CG_pos = impCG %>%
  filter(logFC >= 0) %>% 
  distinct(GeneID) %>% 
  pull(GeneID)
  

PHY_pos = impPHY %>%
  filter(logFC >= 0) %>% 
  distinct(GeneID) %>% 
  pull(GeneID)

NEM_pos = impNEM %>%
  filter(logFC >= 0) %>% 
  distinct(GeneID)  %>% 
  pull(GeneID)

Length_pos = impLength %>%
  filter(logFC >= 0) %>% 
  distinct(GeneID)  %>% 
  pull(GeneID)

venn_pos = list(CG_pos = CG_pos, PHY_pos = PHY_pos, NEM_pos = NEM_pos)



## ---------------------------------------------------------------------------------------------------------------------------------------------
CG_neg = impCG %>%
  filter(logFC <= 0) %>% 
  distinct(GeneID) %>% 
  pull(GeneID)
  

PHY_neg = impPHY %>%
  filter(logFC <= 0) %>% 
  distinct(GeneID) %>% 
  pull(GeneID)

NEM_neg = impNEM %>%
  filter(logFC <= 0) %>% 
  distinct(GeneID)  %>% 
  pull(GeneID)

Length_pos = impLength %>%
  filter(logFC <= 0) %>% 
  distinct(GeneID)  %>% 
  pull(GeneID)

venn_neg = list(CG_neg = CG_neg, PHY_neg = PHY_neg, NEM_neg = NEM_neg)



## ---------------------------------------------------------------------------------------------------------------------------------------------
library(ggvenn)

a = ggvenn(venn_pos, text_size = 6)

b = ggvenn(venn_neg, text_size = 6)

p = tableGrob(summary)

c = grid.arrange(p)

ggarrange(a, b, labels = c('A)', 'B)'), label.y = 0.75)


save_plot('DGEresults/Venn_Pos_Neg.png',
          ggarrange(a, b, labels = c('A)', 'B)'), 
                    font.label = list(size = 25),
                    label.y = 0.75, 
                    heights = c(4,1)),
                    width = 45, 
                    height = 23)


## ---------------------------------------------------------------------------------------------------------------------------------------------
library(ggcorrplot)

# Get matrix for correlation
cordat =  metadata_S %>%
  select(CG_Avg., PHY_Avg, TwoY_RLN, TwoY_length) %>%
  scale() %>%
  as.matrix()

# Calculate p-values
corpmat1 = cor_pmat(cordat, method = 'pearson')

# Calculate correlation coefficients
cormat1 =cor(cordat, method = 'pearson')

# Plot. All relationships are significant
p = ggcorrplot(cormat1, type = 'lower', hc.order = T, lab = T, p.mat = corpmat1)

p

save_plot('DGEresults/Pheno_corrplot.png', p)

