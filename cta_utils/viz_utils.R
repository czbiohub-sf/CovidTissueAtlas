library(forcats)


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
      geom_point(aes(color = Significance ),size = 0.1) +
      scale_color_manual(values = c("grey", "#7fbf7b", "#af8dc3")) +
      theme_bw(base_size = 12) + 
      theme(legend.title = element_blank(),
            text=element_text(size=12)) +
      xlim(-10, 10) + 
      xlab(expression('log'[2]~'(fold change)')) +
      ylab('-log'[10]~'(p-value)') + 
      geom_vline(xintercept=1, linetype='dotted') + 
      geom_vline(xintercept=-1, linetype='dotted') + 
      geom_hline(yintercept=1.3, linetype='dotted') + 
      guides(color=guide_legend(override.aes=list(size=5))) + ggtitle(genes$cell_type %>% unique)
}

plot_DE_direction <- function(master_de = data.frame(), which_tissue = 'liver'){
    
    master_de %>% mutate(expr = ifelse(log2fc >0, 'up','down')) %>% dplyr::filter(tissue ==which_tissue) %>% 
        group_by(tissue, cell_type) %>% summarise(n_up = sum(expr=='up'), n_down = sum(expr=='down')) -> split_expression

    
    split_expression <- split_expression %>%  mutate(cell_type = fct_reorder(cell_type, n_up))
    
    p <- ggplot(data =split_expression) + 
                        geom_bar(aes(x = cell_type, y = n_up), stat ='identity', fill = "#7fbf7b") + 
                        geom_bar(aes(x = cell_type, y = -n_down), stat='identity',fill =  "#af8dc3" ) + 
                        geom_hline(yintercept = 0, color =c("black")) + 
                        coord_flip() + theme_minimal() + 
                        theme(text = element_text(size = 12)) + 
                        scale_fill_viridis(discrete = T, option = "E")  +
                        theme(legend.position = "bottom")  + ggtitle(which_tissue) + 
                        ylab('N DE genes (- down + up)') + xlab('Cell type')
    return(p)
}


# cell phone DB dotplot 

library(ggplot2)
dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename = 'plot.pdf',
                    width = 8,
                    height = 10,
                    means_path = './means.txt',
                    pvalues_path = './pvalues.txt',
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.pdf'
){

  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)

  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]

  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }

  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }

  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]

  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

  ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

#   if (output_extension == '.pdf') {
#       ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
#   }
#   else {
#       ggsave(filename, width = width, height = height, limitsize=F)
#   }
    return(plot.data)
}
