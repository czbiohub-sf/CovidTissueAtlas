
# Functions to parse the output from PathfindR for different organs 
library(stringr)
#Rank top 5 regarding of database 
plotTopPathwaysOrgan <- function(path_df = data.frame() , top_pathways = 10, min_pval = 0.001 , max_logp = 6 , 
                                max_FC = 10, max_name_length = 30 ){
    path_df %>% dplyr::filter(highest_p < min_pval ) %>% group_by(cell_type, database)  %>% 
        top_n(n = -top_pathways, wt = highest_p)  %>% 
        select(Term_Description, cell_type, Fold_Enrichment, highest_p, support, database) %>% 
        arrange(highest_p) -> ranked_paths 

    
    
    ranked_paths$logp <- -log10(ranked_paths$highest_p)
    ranked_paths$logp[ranked_paths$logp > max_logp] <- max_logp
    ranked_paths$Fold_Enrichment[ranked_paths$Fold_Enrichment > max_FC] <- max_FC

    ranked_paths$Term_full <- ranked_paths$Term_Description 
    ranked_paths$Term_Description <- str_trunc(ranked_paths$Term_Description , width = max_name_length)
    # Get the names, sorted first by lg, then by avg
    ranked_paths$Term_Description <- make.unique(ranked_paths$Term_Description )
    nameorder <- ranked_paths$Term_Description[order(ranked_paths$cell_type, ranked_paths$Fold_Enrichment)]

    # Turn name into a factor, with levels in the order of nameorder
    ranked_paths$Term_Description <- factor(ranked_paths$Term_Description, levels = nameorder)

    ranked_paths %>% ggplot(aes(x = Fold_Enrichment, y = Term_Description)) + 
        scale_colour_brewer(palette = "Set1") + 
        geom_point(aes(colour = cell_type, size = logp )) + theme_bw() + theme(text =element_text(size = 15)) + 
          theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed")
          ) + 
            facet_grid(cell_type ~ ., scales = "free_y", space = "free_y") + 
            theme(strip.background =element_rect(fill="white",colour = "white"))

}


# Filter out low-confidence results 
filterTopPathwaysOrgan <- function( path_df = data.frame(), top_pathways = 20, 
                              min_pval = 0.001, minFC = 5,
                              min_genes_used = 2, 
                              direction ='up'){
   
    if(direction =="up"){
              n_genes_used <- apply(!str_split(path_df$Up_regulated, ",", simplify = T) =="", 1, sum)
              path_df <- path_df %>% dplyr::rename(genes  = Up_regulated )
          }else if(direction =="down"){
              n_genes_used <- apply(!str_split(path_df$Down_regulated, ",", simplify = T) =="", 1, sum)
              path_df <- path_df %>% dplyr::rename(genes  = Down_regulated )
    }
    path_df$n_genes_used = n_genes_used
    path_df %>% dplyr::filter(highest_p < min_pval, Fold_Enrichment > minFC, n_genes_used >min_genes_used ) %>% group_by(cell_type, database)  %>% 
            select(Term_Description, cell_type, Fold_Enrichment, highest_p, support, database, n_genes_used, genes) %>% 
            arrange(highest_p) -> ranked_paths 

    ranked_paths$logp <- -log10(ranked_paths$highest_p)

    return(ranked_paths)
}

# Jan 6th 2021 
# match the genes in a given pathway set with the individual DE values
getGeneLists <- function(pathway_hits = data.frame(), which_class = 'endothelial'){
    tissues <- pathway_hits$tissue %>% unique 
    genes_lists = list() 
    for(t in tissues){
        gene_lists <- pathway_hits %>% dplyr::filter(tissue ==t) %>% pull(genes)
        tissue_genes <- str_split(paste(gene_lists, collapse =",", sep = "") , ",")[[1]] %>% trimws() %>% unique 

        genes_lists[[t]] <-  data.frame(tissue = t, gene = tissue_genes)

        }
        genes_df <- do.call(rbind, genes_lists)
        row.names(genes_df) <- 1:dim(genes_df)[1]
    
        FC_matrix <- genes_df %>% left_join(master_df %>% dplyr::filter(cell_class==which_class), by =c('gene', 'tissue')) %>% select(tissue, gene, log2fc) %>% 
                tidyr::spread(key = gene, value = log2fc, fill = 0)
        tissue_rows <- FC_matrix$tissue
        x <- FC_matrix %>% dplyr::select(-tissue) %>% as.matrix()
        row.names(x) <- tissue_rows
    
    return(x)
}

# spread the tidy dataframe of LF changes into a matrix form
makeFCmatrix <- function(
    cell_type = 'endothelial',
    term_grep_pattern = 'VEGF',
    which_direction = 'up', n_pathways = 20){

    path_df = loadPathways(this_sample = cell_type,direction = which_direction )
    ranked_paths = filterTopPathways(path_df, top_pathways = n_pathways, min_genes_used = 2)

    pathway_hits <- ranked_paths %>% dplyr::filter(grepl(Term_Description, pattern = term_grep_pattern, ignore.case =T)) %>% arrange(tissue) %>% 
                        dplyr::filter(n_genes_used > 2, highest_p < 1e-3) 

    fc_matrix <- getGeneLists(pathway_hits, which_class = cell_type )
    return(fc_matrix)
}

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}