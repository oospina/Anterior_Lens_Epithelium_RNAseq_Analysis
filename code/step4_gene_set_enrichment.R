##
# Gene set enrichment analysis using FGSEA
#

library('tidyverse')
library('fgsea')

# Get hallmark pathways (gene sets downloaded as gmt files from https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)
fps = list.files('./data/', pattern='\\.gmt', full.names=TRUE)
pws_ls = lapply(c('kegg', 'go.bp', 'go.cc', 'reactome'), function(i){
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
names(pws_ls) = c('kegg', 'go.bp', 'go.cc', 'reactome')

rm(fps) # Clean env

# Get size of pathway
pws_sz = lapply(names(pws_ls), function(i){
  ls_tmp = lapply(pws_ls[[i]], function(j){
    return(length(j))
  })
  return(ls_tmp)
})
names(pws_sz) = c('kegg', 'go.bp', 'go.cc', 'reactome')

# Read DE results
compars = readxl::excel_sheets('./results/deg_all_comparisons_results.xlsx')
deg_ls = lapply(compars, function(i){
  df_tmp = readxl::read_excel('./results/deg_all_comparisons_results.xlsx', sheet=i)
})
names(deg_ls) = compars

# Put all DGE tests in single data frame and add rank statistic
## Reverse logFC for controls
## Also calculate weighted ranks
deg_ls = lapply(compars, function(i){
  df_tmp = deg_ls[[i]] %>%
    dplyr::select('gene_symbol', 'log2FoldChange', 'pvalue') %>%
    filter(!is.na(pvalue)) %>%
    add_column(compar=i) %>%
    separate(col=compar, into=c('group1', 'group2'), sep='_vs_', remove=TRUE) %>%
    mutate(pvalue=case_when(pvalue == 0 ~ 1e-310, TRUE ~ pvalue)) %>% # Add very small value to p-value in case of "0"
    mutate(rank_stat_grp1_w= -log10(pvalue) * -log2FoldChange) %>% # Revert logFC since "group1" is the baseline/ctrl
    mutate(rank_stat_grp2_w= -log10(pvalue) * log2FoldChange) %>%
    mutate(rank_stat_grp1= -1 * log2FoldChange) %>% # Revert logFC since "group1" is the baseline/ctrl
    mutate(rank_stat_grp2= log2FoldChange)

  return(df_tmp)
})
names(deg_ls) = compars

# Create named vectors with ranks (only logFC)
ranks_sets = lapply(compars, function(i){
  df_tmp = deg_ls[[i]]
  ls_tmp_1 = sort(setNames(df_tmp[['rank_stat_grp1']], nm=df_tmp[['gene_symbol']]), decreasing=TRUE)
  ls_tmp_2 = sort(setNames(df_tmp[['rank_stat_grp2']], nm=df_tmp[['gene_symbol']]), decreasing=TRUE)
  ls_tmp_3 = sort(setNames(df_tmp[['rank_stat_grp1_w']], nm=df_tmp[['gene_symbol']]), decreasing=TRUE)
  ls_tmp_4 = sort(setNames(df_tmp[['rank_stat_grp2_w']], nm=df_tmp[['gene_symbol']]), decreasing=TRUE)
  ls_tmp = list(ls_tmp_1, ls_tmp_2, ls_tmp_3, ls_tmp_4)
  names(ls_tmp) = c(unique(df_tmp[['group1']]),
                    unique(df_tmp[['group2']]),
                    paste0(unique(df_tmp[['group1']]), '_w'),
                    paste0(unique(df_tmp[['group2']]), '_w'))
  return(ls_tmp)
})
names(ranks_sets) = compars

# Run FGSEA
fgsea_ls = lapply(compars, function(i){
  ls_tmp1 = lapply(names(pws_ls), function(p){
    df_tmp = fgsea(pathways=pws_ls[[p]], stats=ranks_sets[[i]][[1]], minSize=3, scoreType='pos', nproc=4, eps=1e-100)
    return(df_tmp)
  })
  names(ls_tmp1) = names(pws_ls)

  ls_tmp2 = lapply(names(pws_ls), function(p){
    df_tmp = fgsea(pathways=pws_ls[[p]], stats=ranks_sets[[i]][[2]], minSize=3, scoreType='pos', nproc=4, eps=1e-100)
    return(df_tmp)
  })
  names(ls_tmp2) = names(pws_ls)

  ls_tmp3 = lapply(names(pws_ls), function(p){
    df_tmp = fgsea(pathways=pws_ls[[p]], stats=ranks_sets[[i]][[3]], minSize=3, scoreType='pos', nproc=4, eps=1e-100)
    return(df_tmp)
  })
  names(ls_tmp3) = names(pws_ls)

  ls_tmp4 = lapply(names(pws_ls), function(p){
    df_tmp = fgsea(pathways=pws_ls[[p]], stats=ranks_sets[[i]][[4]], minSize=3, scoreType='pos', nproc=4, eps=1e-100)
    return(df_tmp)
  })
  names(ls_tmp4) = names(pws_ls)

  ls_tmp = list(ls_tmp1, ls_tmp2, ls_tmp3, ls_tmp4)
  names(ls_tmp) = names(ranks_sets[[i]])

  return(ls_tmp)
})
names(fgsea_ls) = compars

rm(pws_ls, deg_ls, compars, ranks_sets) # Clean env

# Save scores
saveRDS(fgsea_ls, './results/fgsea_enrichment_scores_kegg_cc_bp_reactome.RDS')
fgsea_ls = readRDS('./results/fgsea_enrichment_scores_kegg_cc_bp_reactome.RDS')

# Generate "bubble" plots of enrichment scores
bp_list = lapply(names(fgsea_ls), function(i){
  bp_pw_ls = lapply(c('kegg', 'go.bp', 'go.cc', 'reactome'), function(j){
    # Too many pathways to visualize for GO.BP
    nes_thr = 0.9
    if(j == 'go.bp'){
      nes_thr = 0.95
    }
    # Plot title
    p_title = paste0("Positive enrichment using logFC wighted by DGE p-values\n(Only gene sets with enrichment score higher than quartile ", nes_thr, ')')
    
    df_tmp = as.data.frame(fgsea_ls[[i]][[3]][[j]]) %>%
      select('pathway', 'padj_1'='padj', 'NES_1'='NES', 'size_1'='size') %>%
      left_join(., as.data.frame(fgsea_ls[[i]][[4]][[j]]) %>%
                  select('pathway', 'padj_2'='padj', 'NES_2'='NES', 'size_2'='size'), by='pathway') %>%
      left_join(., as.data.frame(t(data.frame(pws_sz[[j]]))) %>%
                  rownames_to_column('pathway') %>%
                  dplyr::rename(pw_size=2), by='pathway') %>%
      mutate(size_1=(size_1/pw_size) * 100) %>%
      mutate(size_2=(size_2/pw_size) * 100) %>%
      select(-c('pw_size'))
 
    thr_score1 = quantile(df_tmp[['NES_1']], probs=nes_thr, na.rm=TRUE)
    thr_score2 = quantile(df_tmp[['NES_2']], probs=nes_thr, na.rm=TRUE)
    signif_pw = df_tmp %>% filter(padj_1 < 0.05 | padj_2 < 0.05) %>% pull(pathway) %>% unique()
    pw_high_nes =  df_tmp %>% filter(NES_1 < -thr_score1 | NES_1 > thr_score1 |
                                       NES_2 < -thr_score2 | NES_2 > thr_score2) %>% pull(pathway) %>% unique()
    
    df_tmp = df_tmp %>% filter(pathway %in% unique(c(signif_pw, pw_high_nes))) %>%
      mutate(mean_NES=rowMeans(select(., NES_1, NES_2))) %>%
      mutate(pathway=str_replace(pathway, '^REACTOME_|^GOBP_|^GOCC_|^KEGG_', '')) %>%
      mutate(pathway=case_when(is.na(padj_1) | is.na(padj_2) ~ paste0(pathway, '_-.-'),
                               padj_1 >= 0.05 & padj_2 >= 0.05 ~ paste0(pathway, '_-.-'),
                               padj_1 < 0.05 & padj_2 < 0.05 ~ paste0(pathway, '_*.*'),
                               padj_1 < 0.05 & padj_2 >= 0.05 ~ paste0(pathway, '_*.-'),
                               padj_1 >= 0.05 & padj_2 < 0.05 ~ paste0(pathway, '_-.*'))) %>%
      select(-c('padj_1', 'padj_2')) %>%
      pivot_longer(cols=c('NES_1', 'NES_2', 'size_1', 'size_2', ), names_to=c(".value", "group"), names_sep = "_") %>%
      mutate(group=ifelse(group == 1, gsub('_w$', '', names(fgsea_ls[[i]])[3]), gsub('_w$', '', names(fgsea_ls[[i]])[4]))) %>%
      arrange(desc(mean_NES))
    
    bp = ggplot(df_tmp) +
      geom_point(aes(x=NES, y=reorder(pathway, dplyr::desc(mean_NES)), color=NES, size=size)) +
      ylab('') + xlab("Normalized enrichment score (NES)") + ggtitle(p_title) +
      khroma::scale_color_YlOrBr() +
      theme(panel.background=element_rect(color='grey10', fill=NULL)) +
      facet_wrap(~group)
    
    return(bp)
  })
  names(bp_pw_ls) = c('kegg', 'go.bp', 'go.cc', 'reactome')
  
  return(bp_pw_ls)
})
names(bp_list) = names(fgsea_ls)

# Print plots to file
lapply(names(bp_list), function(i){
  graphics.off()
  pdf(paste0('./results/', i, '_fgsea_scores_kegg.pdf'), width=10)
  print(bp_list[[i]][['kegg']])
  dev.off()
  
  graphics.off()
  pdf(paste0('./results/', i, '_fgsea_scores_gobp.pdf'), width=20, height=160)
  print(bp_list[[i]][['go.bp']])
  dev.off()
  
  graphics.off()
  pdf(paste0('./results/', i, '_fgsea_scores_gocc.pdf'), width=18, height=60)
  print(bp_list[[i]][['go.cc']])
  dev.off()
  
  graphics.off()
  pdf(paste0('./results/', i, '_fgsea_scores_reactome.pdf'), width=18, height=60)
  print(bp_list[[i]][['reactome']])
  dev.off()
})

# Save results to spreadsheets
lapply(names(fgsea_ls), function(i){
  fgsea_ls_tmp = fgsea_ls[i]
  lapply(c('kegg', 'go.bp', 'go.cc', 'reactome'), function(j){
    
    # Non-weighted logFC FGSEA results
    fgsea_non_w = lapply(names(fgsea_ls_tmp), function(k){
      fgsea_ls_tmp2 = fgsea_ls_tmp[[k]][!grepl('_w$', names(fgsea_ls_tmp[[k]]))]
      names(fgsea_ls_tmp2) = gsub(paste0('^', i, '_'), '', names(fgsea_ls_tmp2))
      
      df_tmp1 = as.data.frame(fgsea_ls_tmp2[[1]][[j]]) %>%
        mutate(leadingEdge=map_chr(leadingEdge, ~ str_c(.x, collapse = "; "))) %>%
        mutate(leadingEdge=paste0('"', leadingEdge, '"'))
      df_tmp2 = as.data.frame(fgsea_ls_tmp2[[2]][[j]]) %>%
        mutate(leadingEdge=map_chr(leadingEdge, ~ str_c(.x, collapse = "; "))) %>%
        mutate(leadingEdge=paste0('"', leadingEdge, '"'))
      
      colnames(df_tmp1)[-1] = paste0(names(fgsea_ls_tmp2)[1], '_', colnames(df_tmp1)[-1])
      colnames(df_tmp2)[-1] = paste0(names(fgsea_ls_tmp2)[2], '_', colnames(df_tmp2)[-1])
      fgsea_ls_tmp2 = df_tmp1 %>% left_join(., df_tmp2, by='pathway')
      
      return(fgsea_ls_tmp2)
    })
    names(fgsea_non_w) = names(fgsea_ls_tmp)
    
    # Weighted logFC FGSEA results
    fgsea_w = lapply(names(fgsea_ls_tmp), function(k){
      fgsea_ls_tmp2 = fgsea_ls_tmp[[k]][grepl('_w$', names(fgsea_ls_tmp[[k]]))]
      names(fgsea_ls_tmp2) = gsub(paste0('^', i, '_'), '', names(fgsea_ls_tmp2))
      
      df_tmp1 = as.data.frame(fgsea_ls_tmp2[[1]][[j]]) %>%
        mutate(leadingEdge=map_chr(leadingEdge, ~ str_c(.x, collapse = "; "))) %>%
        mutate(leadingEdge=paste0('"', leadingEdge, '"'))
      df_tmp2 = as.data.frame(fgsea_ls_tmp2[[2]][[j]]) %>%
        mutate(leadingEdge=map_chr(leadingEdge, ~ str_c(.x, collapse = "; "))) %>%
        mutate(leadingEdge=paste0('"', leadingEdge, '"'))
      
      colnames(df_tmp1)[-1] = paste0(names(fgsea_ls_tmp2)[1], '_', colnames(df_tmp1)[-1])
      colnames(df_tmp2)[-1] = paste0(names(fgsea_ls_tmp2)[2], '_', colnames(df_tmp2)[-1])
      fgsea_ls_tmp2 = df_tmp1 %>% left_join(., df_tmp2, by='pathway')
      
      return(fgsea_ls_tmp2)
    })
    names(fgsea_w) = names(fgsea_ls_tmp)
    
    openxlsx::write.xlsx(fgsea_non_w, paste0('./results/', j, '_', i, '_gene_set_enrichment_logFC.xlsx'))
    openxlsx::write.xlsx(fgsea_w, paste0('./results/', j, '_', i, '_gene_set_enrichment_weighted_logFC.xlsx'))
  })
})

