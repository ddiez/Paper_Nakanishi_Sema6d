plot_cellector <- function(x, title=NULL) {
  (ggplot(x$data, aes(imagerow, imagecol, color=BLA)) +
     geom_point(size=.8) +
     scale_color_distiller(palette="Spectral") +
     geom_path(aes(imagerow, imagecol, group=cluster), data=x$cluster_data, inherit.aes=FALSE) +
     theme(aspect.ratio=1, axis.title.x=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), axis.line=element_blank()) +
     labs(y=title)) +
    ggplot(x$data, aes(imagerow, imagecol, color=CeA)) +
    geom_point(size=.8) +
    scale_color_distiller(palette="Spectral") +
    geom_path(aes(imagerow, imagecol, group=cluster), data=x$cluster_data, inherit.aes=FALSE) +
    theme(aspect.ratio=1) +
    NoAxes() +
    ggplot(x$data, aes(imagerow, imagecol, color=factor(x$cluster, labels=c("Other", "BLA", "CeA")))) +
    geom_point(size=.8) +
    scale_color_manual("Region", values=c("Other"="lightgrey", "BLA"="violetred", "CeA"="limegreen")) +
    geom_path(aes(imagerow, imagecol, group=cluster), data=x$cluster_data, inherit.aes=FALSE) +
    theme(aspect.ratio=1) +
    NoAxes()
}

plot_cellector2 <- function(x, title=NULL) {
  d <- x$data
  d_cluster <- x$cluster
  d_cluster_data <- x$cluster_data
  
  d$imagerow <- -d$imagerow
  #d$imagecol <- -d$imagecol
  d_cluster_data$imagerow <- -d_cluster_data$imagerow
  #d_cluster_data$imagecol <- -d_cluster_data$imagecol
  
  (ggplot(d, aes(imagerow, imagecol, color=BLA)) +
     geom_point(size=.8) +
     scale_color_distiller(palette="Spectral") +
     geom_path(aes(imagerow, imagecol, group=cluster), data=d_cluster_data, inherit.aes=FALSE) +
     theme(aspect.ratio=1, axis.title.x=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), axis.line=element_blank()) +
     labs(y=title)) +
    ggplot(d, aes(imagerow, imagecol, color=CeA)) +
    geom_point(size=.8) +
    scale_color_distiller(palette="Spectral") +
    geom_path(aes(imagerow, imagecol, group=cluster), data=d_cluster_data, inherit.aes=FALSE) +
    theme(aspect.ratio=1) +
    NoAxes() +
    ggplot(d, aes(imagerow, imagecol, color=factor(d_cluster, labels=c("Other", "BLA", "CeA")))) +
    geom_point(size=.8) +
    scale_color_manual("Region", values=c("Other"="lightgrey", "BLA"="violetred", "CeA"="limegreen")) +
    geom_path(aes(imagerow, imagecol, group=cluster), data=d_cluster_data, inherit.aes=FALSE) +
    theme(aspect.ratio=1) +
    NoAxes()
}


#' @rdname plot_spatial
#' @export
plot_spatial <- function(x, features=NULL, image=NULL, size=.5, rotate=FALSE) {
  UseMethod("plot_spatial")
}

#' @rdname plot_spatial
#' @export
plot_spatial.Seurat <- function(x, features=NULL, image=NULL,size=.5, rotate=FALSE) {
  images <- SeuratObject::Images(x)
  if (is.null(image))
    image <- images[1]
  else
    image <- intersect(image, images)
  if (length(image) == 0) stop("image ", image, " not found.")
  
  d <- SeuratObject::GetTissueCoordinates(x, image=image)
  if (!is.null(features)) {
    exprs <- SeuratObject::GetAssayData(x, slot="data")[features, , drop=FALSE]
    sel <- grep(paste0(image, "_"), colnames(exprs))
    d <- cbind(d, Matrix::t(exprs[, sel, drop=FALSE])) 
  }
  
  if (rotate) {
    d$x <- -d$imagerow
    d$y <- -d$imagecol
  }
  else {
    d$x <- d$imagecol
    d$y <- -d$imagerow
  }
  
  p <- ggplot(d, aes(x, y)) +
    geom_point(size=size, stroke=.1) +
    labs(title=image) +
    coord_fixed() #+
    #theme(axis.line=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank())
  
  if (is.null(features))
    p <- p + aes(color=I("black"))
  else
    p <- p + aes(color=.data[[features]]) +
      scale_color_gradientn(colors=Seurat:::SpatialColors(1000))
         
  #       , fill=.data[[feature]])) +
   # geom_point(shape=21, size=size, stroke=.1) +
  #  scale_fill_gradientn(colors=Seurat:::SpatialColors(1000)) + #, limits=c(NA, maxlim[feature])) +
   # labs(title=image, subtitle=feature) +
  p
}

calculate_cpdb_diff <- function(x, y, cells=NULL) {
  # P.values
  p_x <- as.matrix(x$m_pval)
  rownames(p_x) <- x$proteins
  
  p_y <- as.matrix(y$m_pval)
  rownames(p_y) <- y$proteins
  
  # Means
  m_x <- as.matrix(x$m_mean)
  rownames(m_x) <- x$proteins
  
  m_y <- as.matrix(y$m_mean)
  rownames(m_y) <- y$proteins
  
  # Cells and protein pairs.
  if (is.null(cells))
    cells <- intersect(colnames(p_x), colnames(p_y))
  proteins <- intersect(rownames(p_x), rownames(p_y))
  
  p_x <- p_x[proteins, cells]
  p_y <- p_y[proteins, cells]
  
  m_x <- m_x[proteins, cells]
  m_y <- m_y[proteins, cells]
  
  # Scale
  s_x <- t(scale(t(m_x)))
  s_y <- t(scale(t(m_y)))
  
  list(
    x=list(
      pval=p_x,
      mean=m_x,
      smean=s_x
    ),
    y=list(
      pval=p_y,
      mean=m_y,
      smean=s_y
    ),
    diff=m_x-m_y
  )
}

plot_cpdb_diff <- function(x, y, cutoff=0.05, diff=NULL, cells=NULL, ...) {
  res <- calculate_cpdb_diff(x, y, cells=cells)
  
  sel <- res$x$pval < cutoff | res$y$pval < cutoff
  sel_cells <- colSums(sel) > 0
  sel_proteins <- rowSums(sel) > 0
  
  m <- res$diff[sel_proteins, sel_cells]
  
  if (!is.null(diff)) {
    sel_cells <- colSums(abs(m) > diff) > 0
    sel_proteins <- rowSums(abs(m) > diff) > 0
    m <- m[sel_proteins, sel_cells]
  }
  
  Heatmap(m, name="diff", ...)
}

calculate_cpdb_cellxcell_score <- function(x, cutoff=0.05) {
  cells <- x$cells
  m_pval <- x$m_pval
  
  cells <- unique(unlist(strsplit(cells, "\\|")))
  scores <- matrix(0L, ncol=length(cells), nrow=length(cells), dimnames=list(cells, cells))
  
  d <- m_pval |> 
    pivot_longer(everything(), names_to="pair", values_to="pval")
  
  if (!is.null(cutoff)) 
    d <- d |> filter(pval < cutoff)
  
  d <- d |> separate("pair", c("cell_a", "cell_b"), sep="\\|")
  d <- d |> count(cell_a, cell_b) |> arrange(desc(n))
  
  for (k in seq_len(nrow(d))) {
    l <- d[k, ]
    scores[l[[1]], l[[2]]] <- scores[l[[1]], l[[2]]] + l[[3]]
    if (l[[1]] != l[[2]])
      scores[l[[2]], l[[1]]] <- scores[l[[2]], l[[1]]] + l[[3]] 
  }
  scores
}

#library(pheatmap)
heatmaps_plot = function(meta_file, pvalues_file, count_filename, log_filename, count_network_filename, interaction_count_filename, count_network_separator, interaction_count_separator, show_rownames = T, show_colnames = T,
                         scale="none", cluster_cols = T,border_color='white', cluster_rows = T, fontsize_row=11,
                         fontsize_col = 11, main = '',treeheight_row=0, family='Arial', treeheight_col = 0,
                         col1 = "dodgerblue4", col2 = 'peachpuff', col3 = 'deeppink4', meta_sep='\t', pvalues_sep='\t', pvalue=0.05, log=FALSE){
  #######   Network
  
  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)
  
  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]
  
  
  split_sep = '\\|'
  join_sep = '|'
  
  pairs1_all = unique(meta[,2])
  
  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
  
  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    
    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
    else
      count1 = length(unique(n1))
    
    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }
  
  all_count = all_count[-1,]
  #write.table(all_count, count_network_filename, sep=count_network_separator, quote=F, row.names = F)
  
  #######   count interactions
  
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))
    
  }
  
  if (any(count1)>0)
  {
    count_matrix = matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
    rownames(count_matrix)= unique(meta[,2])
    colnames(count_matrix)= unique(meta[,2])
    
    all_sum = rowSums(count_matrix)
    all_sum = cbind(names(all_sum), all_sum)
    #write.table(all_sum, file=interaction_count_filename, quote=F, sep=count_network_separator, row.names=F)
    
    col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
    
    if (log)
      count_matrix <- log(count_matrix + 1)
    
    pheatmap::pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
                       border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
                       main = main, treeheight_row = treeheight_row,color = col.heatmap, treeheight_col = treeheight_col)
    
    # pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
    #          border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
    #          main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = count_filename)
    # 
    # pheatmap(log(count_matrix+1), show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
    #          border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
    #          main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = log_filename)
  } else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
  }
}

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
  
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  
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
  
  # if (output_extension == '.pdf') {
  #     ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
  # }
  # else {
  #     ggsave(filename, width = width, height = height, limitsize=F)
  # }
}