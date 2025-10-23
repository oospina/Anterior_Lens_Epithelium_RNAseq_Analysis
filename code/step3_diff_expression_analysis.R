##
# Differential gene expression analysis
#

library('tidyverse')
library('DESeq2')
library('ComplexHeatmap')
library('VennDiagram')

# Read color palettes
col_pal = readRDS('./data/color_palettes.RDS')

# Read DESeq object
deseq_obj = readRDS('./data/eye_corneal_endothelium_updated_model_deseq_obj.RDS')

# Perform DE analysis
deseq_res = DESeq(deseq_obj)

# Extract metadata
sample_meta = as.data.frame(colData(deseq_obj))

# Get p-values for corneal endothelium vs naive endothelium
corn_vs_naiv_df = as.data.frame(results(deseq_res, contrast=c('tissue_type', 'corneal_endo', 'naive_endo')))

# Get p-values for lens epithelium vs naive endothelium
lens_vs_naiv_df = as.data.frame(results(deseq_res, contrast=c('tissue_type', 'lens_epith', 'naive_endo')))

# Get p-values for lens epiithelium vs corneal endothelium
lens_vs_corn_df = as.data.frame(results(deseq_res, contrast=c('tissue_type', 'lens_epith', 'corneal_endo')))

# Save results to Excel
openxlsx::write.xlsx(list(corneal_endo_vs_naive_endo=corn_vs_naiv_df %>% rownames_to_column('gene_symbol') %>% arrange(padj, desc(log2FoldChange)),
                          lens_epith_vs_naive_endo=lens_vs_naiv_df %>% rownames_to_column('gene_symbol') %>% arrange(padj, desc(log2FoldChange)),
                          lens_epith_vs_corneal_endo=lens_vs_corn_df %>% rownames_to_column('gene_symbol') %>% arrange(padj, desc(log2FoldChange))), 
                     file='./results/deg_all_comparisons_results.xlsx')

# Create table to do volcano plot
deg_df = corn_vs_naiv_df %>% add_column(compar='corneal_endo_vs_naive_endo') %>% rownames_to_column('gene_symbol') %>%
  bind_rows(., lens_vs_naiv_df %>% add_column(compar='lens_epith_vs_naive_endo') %>% rownames_to_column('gene_symbol')) %>%
  bind_rows(., lens_vs_corn_df %>% add_column(compar='lens_epith_vs_corneal_endo') %>% rownames_to_column('gene_symbol')) %>%
  mutate(signif=case_when(padj < 0.05 & log2FoldChange > 0 ~ 'Up',
                          padj < 0.05 & log2FoldChange < 0 ~ 'Down',
                          TRUE ~ 'Not_DE')) %>%
  filter(!is.na(padj)) %>%
  mutate(signif=factor(signif, levels=c('Up', 'Down', 'Not_DE'))) %>%
  mutate(genes_highlight=ifelse(gene_symbol %in% c('ATP1A1', 'TJP1', 'VIM', 'GJA1', 
                                                   'CRYAA', 'CRYAB', 'CRYBB2', 'FOXC1', 'LMX1B', 
                                                   'POU3F3', 'COL12A1', 'EMILIN1'), 
                                gene_symbol, NA)) # Genes to be highlighted for paper

# Make volcano plot  
volc_p = ggplot() +
  geom_point(data=deg_df, aes(x=log2FoldChange, y=-log10(pvalue), color=signif), size=0.5) +
  ggrepel::geom_text_repel(data=deg_df %>% filter(padj < 0.001 & abs(log2FoldChange) >= 5 & is.na(genes_highlight) & !str_detect(gene_symbol, '^LOC[0-9]+')),
                           aes(x=log2FoldChange, y=-log10(pvalue), color=signif, label=gene_symbol),
                           size=3, show.legend=F, force_pull=2, max.iter=100, max.overlaps=5) +
  ggrepel::geom_label_repel(data=deg_df %>% filter(!is.na(genes_highlight)), # HIGHLIHGTED GENES
                            aes(x=log2FoldChange, y=-log10(pvalue), color=signif, label=genes_highlight), fill=alpha(c("white"), 0.5), 
                            size=3, show.legend=F) +
  labs(title='DE genes corneal vs lens (DESeq2 w/ surrogate variables)',
       x='log2(fold change)', y='-log10(p-value)', color='DE (p<0.05)') +
  ylim(c(0, 175)) +
  scale_color_manual(values=c(Up='red', Down='blue', Not_DE='black')) +
  guides(color=guide_legend(override.aes=list(size=3))) +
  theme_bw() +
  theme(text=element_text(size=16)) +
  facet_wrap(~compar, scales='free')

# Save volcano plot
graphics.off()
pdf('./results/volcano_plot_all_comparisons.pdf', width=20, height=6)
print(volc_p)
dev.off()

rm(corn_vs_naiv_df, lens_vs_naiv_df, lens_vs_corn_df, volc_p) # Clean env

# Make heatmap
## Make annotation table
ann_df = sample_meta %>% 
  select(c('sex', 'tissue_type', 'study')) %>%
  mutate(sex=replace_na(sex, replace='unknown'))

# Get normalized counts
hm_mtx = assay(vst(deseq_obj, blind=FALSE))

# Plot LONG heatmap
## Define genes to plot (DEGs + highlighted genes)
genes_hm = deg_df %>% 
  filter(padj < 0.001 & abs(log2FoldChange) >= 10 & 
           !str_detect(gene_symbol, '^LOC[0-9]+|^SNOR[A-Z][0-9]+|^RN[UY][0-9]+|^TRN[A-Z]+|^SCARNA')) %>% # Remove gene classes 
  pull(gene_symbol) ## DE genes and strong logFC
genes_hm = unique(c(genes_hm, 
                    'CRYAA', 'CRYAB', 'CRYBB2', 'SLC4A11', 'ALCAM', 'ATP1A1', 'PITX2', 'PRDX6',
                    'CLDN11', 'AJAP1', 'TMEM204', 'TMEM178A', 'MKI67', 'CENPF', 
                    'PCNA', 'TP53', 'COL8A1', 'COL4A2', 'FOXC1', 'LMX1B', 'POU3F3',
                    'COL12A1', 'EMILIN1')) # Add genes of interest
hm_mtx_long = t(scale(t(hm_mtx[genes_hm, ]))) # Scale by row (gene)

## Match order of samples in metadata and expression data
ann_df = ann_df[match(colnames(hm_mtx_long), rownames(ann_df)), , drop=FALSE] # Match order of columns

## Create modified legend
lgd_list = list(
  Legend(labels=names(col_pal[['tissue_col_pal']]), 
         legend_gp=gpar(fill=as.vector(col_pal[['tissue_col_pal']]))),
  Legend(labels=names(col_pal[['sex_col_pal']]), 
         legend_gp=gpar(fill=as.vector(col_pal[['sex_col_pal']]))),
  Legend(labels=names(col_pal[['study_col_pal']]), 
         legend_gp=gpar(fill=as.vector(col_pal[['study_col_pal']]))))

## Combine legends
combined_legend = packLegend(lgd_list[[1]], lgd_list[[2]], lgd_list[[3]], direction="vertical")

## Create heatmap object
hm_ann = HeatmapAnnotation(df=ann_df, show_legend=FALSE,
                           col=list(
                             tissue_type=col_pal[['tissue_col_pal']],
                             sex=col_pal[['sex_col_pal']],
                             study=col_pal[['study_col_pal']]))

# Generate plot
hm_long = Heatmap(hm_mtx_long, 
                  name='Scaled gene\nexpression',
                  column_title='DEGs all comparisons (DESeq2 w/ surrogate variables)\nAdj. p-value<0.001; log2FC>=12 + genes of interest',
                  cluster_rows=TRUE, 
                  cluster_columns=TRUE, 
                  show_row_names=TRUE, 
                  show_column_names=TRUE, 
                  show_row_dend=FALSE,
                  bottom_annotation=hm_ann)

graphics.off()
pdf('./results/heatmap_deg_all_comparisons.pdf', height=14)
draw(hm_long, annotation_legend_list=NULL, heatmap_legend_list=combined_legend, 
     heatmap_legend_side="right", annotation_legend_side="right")
dev.off()

rm(hm_mtx, hm_long, hm_mtx_long, lgd_list, hm_ann, ann_df, combined_legend, genes_hm) # Clean env

# Upset plot
gene_collections = list(corn_vs_naiv_DEGs=deg_df %>% filter(compar == 'corneal_endo_vs_naive_endo' & padj < 0.05) %>% pull(gene_symbol),
                        lens_vs_naiv_DEGs=deg_df %>% filter(compar == 'lens_epith_vs_naive_endo' & padj < 0.05) %>% pull(gene_symbol),
                        lens_vs_corn_DEGs=deg_df %>% filter(compar == 'lens_epith_vs_corneal_endo' & padj < 0.05) %>% pull(gene_symbol),
                        corn_vs_naiv_NotDE=deg_df %>% filter(compar == 'corneal_endo_vs_naive_endo' & padj >= 0.05) %>% pull(gene_symbol),
                        lens_vs_naiv_NotDE=deg_df %>% filter(compar == 'lens_epith_vs_naive_endo' & padj >= 0.05) %>% pull(gene_symbol),
                        lens_vs_corn_NotDE=deg_df %>% filter(compar == 'lens_epith_vs_corneal_endo' & padj >= 0.05) %>% pull(gene_symbol))

comb_mat = make_comb_mat(gene_collections, mode='intersect')

up_p = UpSet(comb_mat[comb_degree(comb_mat) == 2], 
             set_order=c("corn_vs_naiv_DEGs",
                         "lens_vs_naiv_DEGs",
                         "lens_vs_corn_DEGs",
                         "corn_vs_naiv_NotDE",
                         "lens_vs_naiv_NotDE",
                         "lens_vs_corn_NotDE"),
             comb_order=order(comb_size(comb_mat[comb_degree(comb_mat) == 2])),
             bg_col=c("lightblue", "lightblue", "lightblue", "pink", "pink","pink"),
             top_annotation=upset_top_annotation(comb_mat[comb_degree(comb_mat) == 2], add_numbers=TRUE))

graphics.off()
pdf('./results/upset_deg_Notdeg_sets.pdf', height=5)
print(up_p)
dev.off()

rm(gene_collections, comb_mat, up_p) # Clean env

# Venn diagrams
gene_sets = deg_df %>% filter(padj < 0.05)
gene_sets = list(corn_UpDEG=gene_sets %>% filter(compar == 'corneal_endo_vs_naive_endo' & log2FoldChange > 0 |
                                                   compar == 'lens_epith_vs_corneal_endo' & log2FoldChange < 0) %>% pull(gene_symbol) %>% unique(),
                 corn_DownDEG=gene_sets %>% filter(compar == 'corneal_endo_vs_naive_endo' & log2FoldChange < 0 |
                                                     compar == 'lens_epith_vs_corneal_endo' & log2FoldChange > 0) %>% pull(gene_symbol) %>% unique(),
                 naiv_UpDEG=gene_sets %>% filter(compar == 'corneal_endo_vs_naive_endo' & log2FoldChange < 0 |
                                                   compar == 'lens_epith_vs_naive_endo' & log2FoldChange < 0) %>% pull(gene_symbol) %>% unique(),
                 naiv_DownDEG=gene_sets %>% filter(compar == 'corneal_endo_vs_naive_endo' & log2FoldChange > 0 |
                                                     compar == 'lens_epith_vs_naive_endo' & log2FoldChange > 0) %>% pull(gene_symbol) %>% unique(),
                 lens_UpDEG=gene_sets %>% filter(compar == 'lens_epith_vs_naive_endo' & log2FoldChange > 0 |
                                                   compar == 'lens_epith_vs_corneal_endo' & log2FoldChange > 0) %>% pull(gene_symbol) %>% unique(),
                 lens_DownDEG=gene_sets %>% filter(compar == 'lens_epith_vs_naive_endo' & log2FoldChange < 0 |
                                                     compar == 'lens_epith_vs_corneal_endo' & log2FoldChange < 0) %>% pull(gene_symbol) %>% unique())
## Up regulated genes
venn_up_p = ggvenn(gene_sets[grep('UpDEG', names(gene_sets))],
                   show_elements=FALSE,
                   stroke_size=0.3,
                   set_name_size=5,
                   fill_color=c("gold", "red", 'blue')) +
  ggtitle("Up-regulated genes in each tissue type\n(across all DE comparisons)")

## Down regulated genes
venn_down_p = ggvenn(gene_sets[grep('DownDEG', names(gene_sets))],
                   show_elements=FALSE,
                   stroke_size=0.3,
                   set_name_size=5,
                   fill_color=c("gold", "red", 'blue')) +
  ggtitle("Down-regulated genes in each tissue type\n(across all DE comparisons)")

graphics.off()
pdf('./results/venn_diagrams_dif_expr_genes.pdf', width=5, height=5)
print(venn_up_p)
print(venn_down_p)
dev.off()

