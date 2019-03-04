# ShinyApp for heatmap
# rm(list=ls())
options(stringsAsFactors = F)

library(shiny)

# ==== Data ====
# orthologs
mah_genes = readRDS("data/mah_genes.RDS")
# sample annotation
lin2_pd = readRDS("data/lin2_pd.RDS")
# count matrices
lin2_probe_QNlog2trans = readRDS("data/lin2_probe_QNlog2trans.RDS")
lin2_gene_QNlog2trans = readRDS("data/lin2_gene_QNlog2trans.RDS")


# ==== Clustering ====
# DATA
## Gene Counts: Complete, Euclidean Distance
d_gene_counts=dist(t(lin2_gene_QNlog2trans))
hc_gene_counts = hclust(d_gene_counts)
## Probe Counts: Complete, Euclidean Distance
d_probe_counts=dist(t(lin2_probe_QNlog2trans))
hc_probe_counts = hclust(d_probe_counts)

library(rafalib)
library(grid)
library(pheatmap)
library(RColorBrewer)
# ==== Heatmaps Setup ====
# ## Edit body of pheatmap:::draw_colnames, customizing it to your liking
# draw_colnames_45 <- function (coln, gaps, ...) {
#   coord = pheatmap:::find_coordinates(length(coln), gaps)
#   x = coord$coord - 0.5 * coord$size
#   res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
#   return(res)}
# ## 'Overwrite' default draw_colnames with your own version
# assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))
# 
# ## edit pheatmap:::scale_vec_colours, to add gray to NA values
# # assign gray to NA entries in heatmap
# scale_vec_colours <- function (x, col = rainbow(10), breaks = NA)
# {
#   #return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
#   res = col[as.numeric(cut(x, breaks = breaks, include.lowest = T))]
#   res[is.na(res)] = "white"
#   return(res)
# }
# my_scale_vec_colours <- function(x, col = rainbow(10), breaks = NA){
#   res = col[as.numeric(cut(x, breaks = breaks, include.lowest = T))]
#   res[is.na(res)] = "white"
#   return(res)
# }
# #
# assignInNamespace(x="scale_vec_colours", value="my_scale_vec_colours",ns=asNamespace("pheatmap"))
# 

server <- function(input, output) {
  
  data_set = reactive({
    if(input$dat == "genes"){
      d_mat = as.matrix(d_gene_counts)
      dimnames(d_mat) = list(
        paste(substring(lin2_pd$species,1,1),lin2_pd$tissue,sep="."),
        paste(substring(lin2_pd$species,1,1),lin2_pd$tissue,sep=".")
      )
      hcout = hc_gene_counts
      MAIN = "Gene Counts"
    }else if(input$dat == "probes"){
      d_mat = as.matrix(d_probe_counts)
      dimnames(d_mat) = list(
        paste(substring(lin2_pd$species,1,1),lin2_pd$tissue,sep="."),
        paste(substring(lin2_pd$species,1,1),lin2_pd$tissue,sep=".")
      )
      hcout = hc_probe_counts
      MAIN = "Probe Counts"
    }
    return(
      list(
        d_mat = d_mat,
        hcout = hcout,
        MAIN = MAIN
      )
    )
  })
  # data set
  
  output$distPlot <- renderPlot({
    
    #pheatmap(data_set()$d_mat)
    COL = ifelse(input$color_pal == "pheatmap","RdYlBu",input$color_pal)
    
    if(input$mat_display == "lower_tri"){
      tmp_mat = data_set()$d_mat[data_set()$hcout$order,data_set()$hcout$order]
      tmp_mat[upper.tri(tmp_mat)] = NA
      tmp_mat=tmp_mat[
        match(colnames(data_set()$d_mat),colnames(tmp_mat)),
        match(colnames(data_set()$d_mat),colnames(tmp_mat))
        ]
      # plot heatmap
      pheatmap(
        mat=tmp_mat[,data_set()$hcout$order],
        cluster_rows=data_set()$hcout,
        cluster_cols=F,
        scale="none",
        color=colorRampPalette(rev(brewer.pal(n=7,name=COL)))(100),
        border_color="white",
        fontsize_number=9,
        fontsize=9,
        fontsize_row=7,
        fontsize_col=7,
        main=data_set()$MAIN
      )
    }else if(input$mat_display == "full"){
      # plot heatmap
      pheatmap(
        mat=data_set()$d_mat,
        cluster_rows=data_set()$hcout,
        cluster_cols=data_set()$hcout,
        scale="none",
        color=colorRampPalette(rev(brewer.pal(n=7,name=COL)))(100),
        border_color="white",
        fontsize_number=9,
        fontsize=9,
        fontsize_row=7,
        fontsize_col=7,
        main=data_set()$MAIN
      )
    }
  })
  
  
}

