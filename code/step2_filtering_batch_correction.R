##
# Assess batch effects in data set
#

library('tidyverse')
library('DESeq2')

# Read DESeq object
deseq_obj = readRDS('./data/eye_corneal_endothelium_deseq_obj.RDS')

# Keep genes expressed at least in half of the samples
genes_keep = rowSums(assay(deseq_obj) == 0) <= round(ncol(deseq_obj)*0.5, 0)
deseq_obj = deseq_obj[genes_keep, ]

rm(genes_keep) # Clean env

# Calculate library size factors
deseq_obj = estimateSizeFactors(deseq_obj)

# Estimate SVs
sva_res = sva::svaseq(dat=assay(vst(deseq_obj, blind=TRUE)),
                      mod=model.matrix(~tissue_type, data=as.data.frame(colData(deseq_obj))),
                      mod0=model.matrix(~1, data=as.data.frame(colData(deseq_obj))), n.sv=2)

# Extract SVs to add to DESeq object
sva_mtx = sva_res[['sv']]
colnames(sva_mtx) = paste0('SV', 1:ncol(sva_mtx))
colData(deseq_obj) = S4Vectors::DataFrame(cbind(colData(deseq_obj), sva_mtx))

rm(sva_res, sva_mtx) # Clean env

# Save DESeq object with SVs
saveRDS(deseq_obj, './data/eye_corneal_endothelium_w_SVs_deseq_obj.RDS')
deseq_obj = readRDS('./data/eye_corneal_endothelium_w_SVs_deseq_obj.RDS')

# Plot SVs
p1 = ggplot(as.data.frame(colData(deseq_obj))) +
  geom_point(aes(x=SV1, y=SV2, color=study, shape=tissue_type), size=3)

graphics.off()
pdf('./results/surrogate_vars_plot.pdf', height=4, width=5)
print(p1)
dev.off()

rm(p1) # Clean env

# Identify high variance genes
high_var_genes = apply(log2(assay(vst(deseq_obj, blind=FALSE))), 1, sd) %>%
  sort(decreasing=TRUE) %>%
  head(5000)

# Plot PCs before and after SV correction
pca_res = assay(vst(deseq_obj, blind=TRUE))[names(high_var_genes), ]
pca_res = prcomp(t(pca_res), center=TRUE, scale=TRUE)
expl_var = round(as.vector(summary(pca_res)[['importance']][2, 1:2])*100, 2)
## With correction
pca_res_corr = assay(vst(deseq_obj, blind=FALSE))[names(high_var_genes), ]
pca_res_corr = limma::removeBatchEffect(pca_res_corr, covariates=as.data.frame(colData(deseq_obj))[, c('SV1', 'SV2')] )
pca_res_corr = prcomp(t(pca_res_corr), center=TRUE, scale=TRUE)
expl_var_corr = round(as.vector(summary(pca_res_corr)[['importance']][2, 1:2])*100, 2)

pcp_1 = ggplot(as.data.frame(pca_res[['x']]) %>% select(1:2) %>% 
                 rownames_to_column('sample_name') %>%
                 left_join(., as.data.frame(colData(deseq_obj)) %>%
                             rownames_to_column('sample_name'), by='sample_name')) +
  geom_point(aes(x=PC1, y=PC2, color=study, shape=tissue_type), size=3) +
  labs(x=paste0('PC1 (', expl_var[1], ')'), y=paste0('PC2 (', expl_var[2], ')'))

pcp_2 = ggplot(as.data.frame(pca_res_corr[['x']]) %>% select(1:2) %>% 
                 rownames_to_column('sample_name') %>%
                 left_join(., as.data.frame(colData(deseq_obj)) %>%
                             rownames_to_column('sample_name'), by='sample_name')) +
  geom_point(aes(x=PC1, y=PC2, color=study, shape=tissue_type), size=3) +
  labs(x=paste0('PC1 (', expl_var_corr[1], ')'), y=paste0('PC2 (', expl_var_corr[2], ')'))

graphics.off()
pdf('./results/pca_pre_and_post_sva_correction.pdf', height=5, width=8)
print(ggpubr::ggarrange(pcp_1, pcp_2, ncol=2, common.legend=TRUE, legend='bottom'))
dev.off()

rm(pcp_1, pcp_2, expl_var, expl_var_corr, pca_res, pca_res_corr) # Clean env

# Save DESeq object with updated model
design(deseq_obj) = formula(~SV1+SV2+tissue_type)
saveRDS(deseq_obj, './data/eye_corneal_endothelium_updated_model_deseq_obj.RDS')

