volcano_plot <- function(genes = data.frame() ){
    
    genes$genes <- genes$gene
    genes$logFC <- genes$log2fc 
    genes$PValue <- genes$pval
    # https://www.biostars.org/p/309962/
    AdjustedCutoff <- 0.05
    LabellingCutoff <- 0.05
    FCCutoff <- 1
    neg_FCCutoff <- -1
    genes$Significance <- ""
   
    genes$Significance[(genes$PValue<AdjustedCutoff) & (genes$logFC>FCCutoff)] <- 'Elevated in COVID-19'
    genes$Significance[(genes$PValue<AdjustedCutoff) & (genes$logFC<neg_FCCutoff)] <- 'Reduced in COVID-19'
    table(genes$Significance)
    genes$Significance <- factor(genes$Significance, levels=c("", 
                                                              'Elevated in COVID-19', 
                                                              'Reduced in COVID-19'))


    ggplot(genes, aes(x = logFC, y = -log10(PValue))) +
      geom_point(aes(color = Significance)) +
      scale_color_manual(values = c("grey", "#7fbf7b", "#af8dc3")) +
      theme_bw(base_size = 25) + 
      theme(legend.title = element_blank(),
            text=element_text(size=20)) +
      xlim(-10, 10) + 
      xlab(expression('log'[2]~'(fold change)')) +
      ylab('-log'[10]~'(p-value)') + 
      geom_vline(xintercept=1, linetype='dotted') + 
      geom_vline(xintercept=-1, linetype='dotted') + 
      geom_hline(yintercept=1.3, linetype='dotted') + 
      guides(color=guide_legend(override.aes=list(size=5))) + ggtitle(genes$cell_type %>% unique)
}