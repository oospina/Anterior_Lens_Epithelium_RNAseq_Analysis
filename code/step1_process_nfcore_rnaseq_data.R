##
# Load STAR-Salmon-generated counts and preprocess for DGE analysis
#

library('tidyverse')
library('DESeq2')

# Read transcript to gene name table (from nfcore/rnaseq/salmon)
fp = './nfcore_rnaseq_run/nfcore_rnaseq_output/star_salmon/tx2gene.tsv'
tx2gene = readr::read_delim(fp, delim='\t')
# Collapse ambiguous transcript IDs to single gene symbol
tx2gene = tx2gene %>% mutate(gene_id=str_replace(gene_id, "_[0-9]{1,2}$", ''))
tx2gene = tx2gene[, c(1:2)]

# Read and process sample meta data
fp = './data/Cornea RNAseq samples.xlsx'
sample_meta = readxl::read_excel(fp, .name_repair='minimal') %>% 
  janitor::clean_names() %>%
  mutate(tissue_type=factor(case_when(str_detect(rn_aseq_samples, 'corneal') ~ 'corneal_endo',
                                      str_detect(rn_aseq_samples, 'lens') ~ 'lens_epith',
                                      TRUE ~ 'naive_endo'))) %>%
  mutate(study=factor(case_when(tissue_type == 'naive_endo' ~ 'tokuda_etal',
                                TRUE ~ 'nagy_lab')))
sample_meta = sample_meta[!is.na(sample_meta[[2]]), colSums(is.na(sample_meta)) != nrow(sample_meta)]
sample_meta = sample_meta[, -2] # Remove column with filename
sample_meta[['samplename']] = c("corn_endo_1", "corn_endo_2", "corn_endo_3", "corn_endo_4",
                                "lens_epi_1", "lens_epi_2", "lens_epi_3", "lens_epi_4",
                                "naiv_endo_1", "naiv_endo_2", "naiv_endo_3", "naiv_endo_4")
sample_meta = sample_meta %>% column_to_rownames('samplename')

rm(fp) # Clean env

# Pass salmon quant files to tximport
fp = list.files('./nfcore_rnaseq_run/nfcore_rnaseq_output/star_salmon/', recursive=TRUE, pattern='quant\\.sf', full.names=TRUE)
names(fp) = stringr::str_extract(fp, 'corn_endo_[1-4]|lens_epi_[1-4]|naiv_endo_[1-4]')
txi = tximport::tximport(fp, type='salmon', tx2gene=tx2gene)

rm(fp, tx2gene) # Clean env

# make sure order of samples is the same in metadata and counts
sample_meta = sample_meta[match(colnames(txi[['counts']]), rownames(sample_meta)), ]

# Create initial DESeq object (no batch correction)
deseq_obj = DESeqDataSetFromTximport(txi=txi,
                                     colData=sample_meta,
                                     design=~tissue_type)

saveRDS(deseq_obj, './data/eye_corneal_endothelium_deseq_obj.RDS')

