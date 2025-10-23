##
# Generation of figures for manuscript
#

library('tidyverse')
library('DESeq2')
library('ComplexHeatmap')
library('ggvenn')

# Read color palettes
col_pal = readRDS('./data/color_palettes.RDS')

# Read DESeq object
deseq_obj = readRDS('./data/eye_corneal_endothelium_updated_model_deseq_obj.RDS')


################################################################################
## FIGURE_03_A - Correlation plot 

# Extract gene expression data from DESeq object and adjust for surrogate variables
expr_mtx = assay(vst(deseq_obj, blind=FALSE))
expr_mtx = limma::removeBatchEffect(expr_mtx, covariates=as.data.frame(colData(deseq_obj))[, c('SV1', 'SV2')] )

# Change sample names
colnames(expr_mtx) = sapply(colnames(expr_mtx), function(x){
  if(grepl("corn_endo", x)){
    gsub('corn_endo', 'CEC', x)
  } else if(grepl("lens_epi", x)){
    gsub('lens_epi', 'ALCEC', x)
  } else if(grepl("naiv_endo", x)){
    gsub('naiv_endo', 'nCEC', x)
  }
}) %>% as.vector()

# Compute Spearman correlation
corr_mtx = cor(expr_mtx, method='spearman')

# Create annotation object for heatmap showing correlations
ann_df = as.data.frame(colData(deseq_obj)) %>% 
  rownames_to_column('sample_name') %>%
  dplyr::select('sample_name', 'sex', 'tissue_type') %>%
  mutate(clean_name=colnames(expr_mtx)) %>%
  column_to_rownames('clean_name')
  
hm_ann = HeatmapAnnotation(df=ann_df %>% dplyr::select(-c('sample_name')),
                           show_legend=FALSE,
                           col=list(tissue_type=col_pal[['tissue_col_pal']],
                                    sex=col_pal[['sex_col_pal']]))

# Create modified legend
lgd_list = list(
  Legend(labels=names(col_pal[['tissue_col_pal']]), 
         legend_gp=gpar(fill=as.vector(col_pal[['tissue_col_pal']]))),
  Legend(labels=names(col_pal[['sex_col_pal']]), 
         legend_gp=gpar(fill=as.vector(col_pal[['sex_col_pal']]))))

# Combine legends
comb_lgd = packLegend(lgd_list[[1]], lgd_list[[2]], direction="vertical")
                      
# Generate plot
hm_p = Heatmap(corr_mtx, 
               name="Spearman's\ncorrelation",
               column_title='Gene expression correlation across samples',
               cluster_rows=TRUE, 
               cluster_columns=TRUE, 
               show_row_names=TRUE, 
               show_column_names=TRUE, 
               col=circlize::colorRamp2(breaks=c(0.6, 0.7, 0.8, 0.9, 1), 
                                        colors=c("#084081", "#1B9E77", "#F5F5F5", "#F46D43", "#A50026")),
               bottom_annotation=hm_ann)

graphics.off()
pdf('./docs/figures_manuscript/FIGURE_03_A_SPEARMANCORR.pdf', width=6, height=5)
draw(hm_p, annotation_legend_list=NULL, heatmap_legend_list=comb_lgd, 
     heatmap_legend_side="right", annotation_legend_side="right")
dev.off()

rm(corr_mtx, hm_ann, comb_lgd, hm_p, lgd_list) # Clean env


################################################################################
## FIGURE_03_B - Volcano plots

# Read DE analysis results
fp = './results/deg_all_comparisons_results.xlsx'
compars = readxl::excel_sheets(fp)
deg_ls = lapply(compars, function(i){
  df_tmp = readxl::read_xlsx(fp, sheet=i) %>%
    add_column(compar=i)
  return(df_tmp)
})
names(deg_ls) = readxl::excel_sheets(fp)

# Create table to do volcano plot
deg_df = do.call('bind_rows', deg_ls) %>%
  mutate(signif=case_when(padj < 0.01 & log2FoldChange > 2 ~ 'Up',
                          padj < 0.01 & log2FoldChange < -2 ~ 'Down',
                          TRUE ~ 'Not_DE')) %>%
  filter(!is.na(padj)) %>%
  mutate(signif=factor(signif, levels=c('Up', 'Down', 'Not_DE'))) %>%
  mutate(genes_highlight=ifelse(gene_symbol %in% c('ATP1A1', 'TJP1', 'VIM', 'GJA1', 
                                                   'CRYAA', 'CRYAB', 'CRYBB2', 'FOXC1', 'LMX1B', 
                                                   'POU3F3', 'COL12A1', 'EMILIN1'), 
                                gene_symbol, NA)) # Genes Dr. Alonso wants highlighted

# Make volcano plot  
volc_p = ggplot() +
  geom_point(data=deg_df, aes(x=log2FoldChange, y=-log10(pvalue), color=signif), size=0.5) +
  ggrepel::geom_text_repel(data=deg_df %>% filter(padj < 0.01 & abs(log2FoldChange) > 2 & is.na(genes_highlight) & !str_detect(gene_symbol, '^LOC[0-9]+')),
                           aes(x=log2FoldChange, y=-log10(pvalue), color=signif, label=gene_symbol),
                           size=3, show.legend=F, force=0.2, force_pull=10, max.iter=100, max.overlaps=3) +
  ggrepel::geom_label_repel(data=deg_df %>% filter(!is.na(genes_highlight)), # HIGHLIHGTED GENES
                            aes(x=log2FoldChange, y=-log10(pvalue), color=signif, label=genes_highlight), fill=alpha(c("white"), 0.5), 
                            size=3, show.legend=F) +
  labs(title='Differentially expressed genes among tissue types\n(FDR 0.01; log2(fold change) > 2 or log2(fold change) < -2)',
       x='log2(fold change)', y='-log10(p-value)', color='DEGs') +
  ylim(c(0, 175)) +
  scale_color_manual(values=c(Up='red', Down='blue', Not_DE='black')) +
  guides(color=guide_legend(override.aes=list(size=3))) +
  theme_bw() +
  theme(text=element_text(size=16)) +
  facet_wrap(~compar, scales='free')

# Save volcano plot
graphics.off()
pdf('./docs/figures_manuscript/FIGURE_03_B_VOLCANO.pdf', width=20, height=6)
print(volc_p)
dev.off()

rm(fp, deg_ls, compars, volc_p) # Clean env


################################################################################
## FIGURE_03_C - UpSet plot

# Create gene collections to display
gene_collections = list(corn_vs_naiv_DEGs=deg_df %>% filter(compar == 'corneal_endo_vs_naive_endo' & padj < 0.01 & abs(log2FoldChange) > 2) %>% pull(gene_symbol),
                        lens_vs_naiv_DEGs=deg_df %>% filter(compar == 'lens_epith_vs_naive_endo' & padj < 0.01 & abs(log2FoldChange) > 2) %>%  pull(gene_symbol),
                        lens_vs_corn_DEGs=deg_df %>% filter(compar == 'lens_epith_vs_corneal_endo' & padj < 0.01 & abs(log2FoldChange) > 2) %>% pull(gene_symbol),
                        corn_vs_naiv_NotDE=deg_df %>% filter(compar == 'corneal_endo_vs_naive_endo' & padj >= 0.01) %>% pull(gene_symbol),
                        lens_vs_naiv_NotDE=deg_df %>% filter(compar == 'lens_epith_vs_naive_endo' & padj >= 0.01) %>% pull(gene_symbol),
                        lens_vs_corn_NotDE=deg_df %>% filter(compar == 'lens_epith_vs_corneal_endo' & padj >= 0.01) %>% pull(gene_symbol))

# Set interaction mode
comb_mat = make_comb_mat(gene_collections, mode='intersect')

# Select intereactions to show
comb_mask = comb_name(comb_mat) %in% c('111000', '110000', '101000', '011000', '000110', '000101', '000011', '000111')

# Make plot
up_p = UpSet(comb_mat[comb_mask], 
             set_order=c("corn_vs_naiv_DEGs",
                         "lens_vs_naiv_DEGs",
                         "lens_vs_corn_DEGs",
                         "corn_vs_naiv_NotDE",
                         "lens_vs_naiv_NotDE",
                         "lens_vs_corn_NotDE"),
             comb_order=order(comb_size(comb_mat[comb_mask])),
             bg_col=c("gray50", "gray50", "gray50", "pink", "pink","pink"),
             column_title='Differentially expressed genes (FDR < 0.01; log2FC > 2 or log2FC < -2)\nNot differentially expressed (FDR >= 0.01)', 
             top_annotation=upset_top_annotation(comb_mat[comb_mask], add_numbers=TRUE))

graphics.off()
pdf('./docs/figures_manuscript/FIGURE_03_C_UPSET.pdf', height=5)
print(up_p)
dev.off()

rm(comb_mat, up_p) # Clean env


################################################################################
## FIGURE_03_E - Heatmap DEGs

# Extract non-DE genes to interrogate GPT-4o
max_length = max(length(gene_collections[['corn_vs_naiv_NotDE']]), 
                 length(gene_collections[['lens_vs_naiv_NotDE']]), 
                 length(gene_collections[['lens_vs_corn_NotDE']]))
genes_gpt = list(cec_ncec=c(gene_collections[['corn_vs_naiv_NotDE']], rep(NA, max_length - length(gene_collections[['corn_vs_naiv_NotDE']]))),
                 alcec_ncec=c(gene_collections[['lens_vs_naiv_NotDE']], rep(NA, max_length - length(gene_collections[['lens_vs_naiv_NotDE']]))),
                 alcec_cec=c(gene_collections[['lens_vs_corn_NotDE']], rep(NA, max_length - length(gene_collections[['lens_vs_corn_NotDE']]))))
genes_gpt = as.data.frame(genes_gpt)
write.csv(genes_gpt, './results/genes_highlight_query_gpt.txt', row.names=FALSE, quote=FALSE)

# Genes highlighted by GPT
genes_gpt = c("COL4A1", "COL8A2", "KRT8", "KRT18", "VIM", 
              "ATP1A1", "ANXA1", "CDH2", "SLC4A4", "BMP4", 
              "CLDN10", "TJP1", "FN1", "SPARC", "GAPDH", 
              "HSPA1A", "PECAM1", "ITGA5", "ACTB", "EFEMP1")

# Other genes selected by Dr. Alonso
genes_high = c('CRYAA', 'CRYAB', 'CRYBB2', 'SLC4A11', 'ALCAM', 'ATP1A1', 'PITX2', 'PRDX6',
               'CLDN11', 'AJAP1', 'TMEM204', 'TMEM178A', 'MKI67', 'CENPF', 
               'PCNA', 'TP53', 'COL8A1', 'COL4A2', 'FOXC1', 'LMX1B', 'POU3F3',
               'COL12A1', 'EMILIN1', genes_gpt)

# Subset data to genes expressed over median expression value
genes_median = apply(expr_mtx, 1, median)
genes_median = genes_median[genes_median > median(genes_median)]
hm_mtx_scl = expr_mtx[names(genes_median), ]

# Transpose matrix and match order of columns with meta data
hm_mtx_scl = hm_mtx_scl[rownames(hm_mtx_scl) %in% unique(unlist(gene_collections)), ]
hm_mtx_scl = hm_mtx_scl[grep('^LOC[0-9]+', rownames(hm_mtx_scl), invert=TRUE), ]
hm_mtx_scl = t(scale(t(hm_mtx_scl))) # Scale by row (gene)
hm_mtx_scl = hm_mtx_scl[, match(rownames(ann_df), colnames(hm_mtx_scl)), drop=FALSE] 

# Generate heatmap annotation object
hm_ann = HeatmapAnnotation(df=ann_df %>% dplyr::select(-c('sample_name')),
                           show_legend=FALSE,
                           col=list(tissue_type=col_pal[['tissue_col_pal']],
                                    sex=col_pal[['sex_col_pal']]))

# Create modified legend
lgd_list = list(
  Legend(labels=names(col_pal[['tissue_col_pal']]), 
         legend_gp=gpar(fill=as.vector(col_pal[['tissue_col_pal']]))),
  Legend(labels=names(col_pal[['sex_col_pal']]), 
         legend_gp=gpar(fill=as.vector(col_pal[['sex_col_pal']]))))

# Combine legends
comb_lgd = packLegend(lgd_list[[1]], lgd_list[[2]], direction="vertical")

# Make annotation to highlight gene names
row_ann = rowAnnotation(foo=anno_mark(at=as.vector(na.omit(match(genes_high, rownames(hm_mtx_scl)))),
                                      labels=genes_high[genes_high %in% rownames(hm_mtx_scl)], 
                                      labels_gp=gpar(fontsize=10)))

# Generate plot
hm_p = Heatmap(hm_mtx_scl, 
               name='Scaled gene\nexpression',
               column_title='Genes expressed over median expression value',
               cluster_rows=TRUE, 
               cluster_columns=TRUE, 
               show_row_names=FALSE, 
               show_column_names=FALSE, 
               show_row_dend=FALSE,
               col=circlize::colorRamp2(c(-4.5, -1.5, 0, 1.5, 4.5), colors=c('darkblue', 'blue', 'white', 'red', 'darkred')),
               bottom_annotation=hm_ann,
               right_annotation=row_ann)

hm_p = draw(hm_p, annotation_legend_list=NULL, heatmap_legend_list=comb_lgd,
            heatmap_legend_side="right", annotation_legend_side="right", padding=unit(c(2, 2, 2, 20), "mm"))

graphics.off()
pdf('./docs/figures_manuscript/FIGURE_03_E_GENEEXPRHEATMAP.pdf')
print(hm_p)
dev.off()

rm(lgd_list, comb_lgd, hm_ann, row_ann, hm_mtx_scl, max_length, hm_p, genes_median, genes_gpt, genes_high) # Clean env


################################################################################
## FIGURE_03_D - Venn diagram

# Venn diagrams
gene_sets = deg_df %>% filter(padj < 0.01 & abs(log2FoldChange) > 2)
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
pdf('./docs/figures_manuscript/FIGURE_03_D_VENNDIAG.pdf', width=5, height=5)
print(venn_up_p)
print(venn_down_p)
dev.off()

rm(gene_sets, venn_down_p, venn_up_p) # Clean env


################################################################################
## FIGURE_03_F - Gene set enrichment

# Get pathways
fps = list.files('/projects/nagy_lab_projects/projects_sergio/eye_corneal_lens_comparison/data/', pattern='\\.gmt', full.names=TRUE)
pws_ls = lapply(c('go.bp', 'go.cc'), function(i){
  pws_raw = readLines(grep(i, fps, value=TRUE))
  pws = lapply(pws_raw, function(i){
    pw_tmp = unlist(strsplit(i, split='\\t'))
    pw_name_tmp = pw_tmp[1]
    pw_genes_tmp = pw_tmp[-c(1:2)]
    return(list(pw_name=pw_name_tmp,
                pw_genes=pw_genes_tmp))
  })
  rm(pws_raw)
  pws_names = c()
  for(i in 1:length(pws)){
    pws_names = append(pws_names, pws[[i]][['pw_name']])
    pws[[i]] = pws[[i]][['pw_genes']]
  }
  names(pws) = pws_names
  
  return(pws)
})
names(pws_ls) = c('go.bp', 'go.cc')

# Get size of pathway
pws_sz = lapply(names(pws_ls), function(i){
  ls_tmp = lapply(pws_ls[[i]], function(j){
    return(length(j))
  })
  return(ls_tmp)
})
names(pws_sz) = c('go.bp', 'go.cc')

# Read FGSEA results
fp = './results/fgsea_enrichment_scores_kegg_cc_bp_reactome.RDS'
fgsea_ls = readRDS(fp)

# Threshold
padj_thr = 0.05
nes_thr = 0.5

# Source pathways selected by GPT4o
source('./code/gpt4o_selected_gene_sets.R')

# Generate "bubble" plots of enrichment scores (CEC or ALCEC vs nCEC)
bp_pw_ls = lapply(c('go.bp', 'go.cc'), function(j){
  # Merge results for each tissue and pathway sizes
  df_tmp = as.data.frame(fgsea_ls[['corneal_endo_vs_naive_endo']][[3]][[j]]) %>%
    dplyr::select('pathway', padj_1='padj', NES_1='NES', size_1='size') %>%
    left_join(., as.data.frame(fgsea_ls[['lens_epith_vs_naive_endo']][[3]][[j]]) %>%
                dplyr::select('pathway', padj_2='padj', NES_2='NES', size_2='size'), by='pathway') %>%
    dplyr::filter(pathway %in% c(bp_merged, cc_merged)) %>%
    left_join(., as.data.frame(t(data.frame(pws_sz[[j]]))) %>%
                rownames_to_column('pathway') %>%
                dplyr::rename(pw_size=2), by='pathway') %>%
    mutate(size_1=(size_1/pw_size) * 100,
           size_2=(size_2/pw_size) * 100)
  
  # Subset to pathways to show
  df_tmp = df_tmp %>%
    arrange(desc(NES_1)) %>%
    select(-c('pw_size')) %>%
    mutate(pathway=str_replace(pathway, '^GOBP_|^GOCC_', '')) %>%
    mutate(pathway=case_when(padj_1 < padj_thr & padj_2 < padj_thr ~ paste0(pathway, '_1_2'),
                             padj_1 < padj_thr & padj_2 >= padj_thr  ~ paste0(pathway, '_1'),
                             padj_1 >= padj_thr & padj_2 < padj_thr  ~ paste0(pathway, '_2'),
                             is.na(padj_1) &  is.na(padj_2) ~ pathway,
                             TRUE ~ pathway)) %>%
    dplyr::select(-c('padj_1', 'padj_2')) %>%
    pivot_longer(cols=c('NES_1', 'NES_2', 'size_1', 'size_2', ), names_to=c(".value", "group"), names_sep = "_") %>%
    mutate(group=ifelse(group == 1, 'CEC', 'ALCEC')) %>%
    arrange(desc(NES))
  
  bp = ggplot(df_tmp) +
    geom_point(aes(x=NES, y=reorder(pathway, NES), color=NES, size=size)) +
    ylab('') + xlab("Normalized enrichment score (NES)") + 
    ggtitle("Positively enriched gene sets\n(logFC wighted by DGE p-values)") +
    scale_color_gradient(low="#2166AC", high="#B2182B")+
    theme(panel.background=element_rect(color='grey10', fill=NULL), 
          axis.text.y=element_text(size=8)) +
    facet_wrap(~group)
  
  
  return(bp)
})
names(bp_pw_ls) = c('go.bp', 'go.cc')

graphics.off()
pdf('./docs/figures_manuscript/FIGURE_03_F_GSEABUBBLE.pdf', width=10, height=5)
print(bp_pw_ls)
dev.off()

