#' @title TarNet
#' @description
#' A function used to plot a network that consists of microbes and metabolites in same pathway.This function offers a network of microbes and metabolites in same pathways.
#'
#' @param meta_mic2func the results of the complete nodescores getting from Nodescore function.
#' @param r correlation coefficient of nodes.Default value is 0.1.
#' @param IIScore the result using IIS() function.
#' @param p_adjust p value the have been adjusted.Default value is 0.05.
#'
#' @export TarNet
#' @import igraph
#' @import ggrepel
#' @import visNetwork
#' @importFrom htmlwidgets JS
#' @import ppcor
#' @examples
#' IIScore=IIS(Microbe_example,Meta_example,confounderData,groupInfo)
#' meta_mic2func=MMfunc(IIScore)
#' TarNet(IIScore,meta_mic2func)


TarNet<-function(IIScore,meta_mic2func,r = 0.1,p_adjust = 0.05){
  mic_meta_cor=IIScore[["mic_meta_cor"]]
  mic_meta_cor_p.adjust=IIScore[["mic_meta_cor_p.adjust"]]

  meta_in_func=meta_mic2func[[1]]
  mic_in_func=meta_mic2func[[2]]

  mic_meta_cor=mic_meta_cor[mic_in_func,meta_in_func]
  mic_meta_cor_p.adjust=mic_meta_cor_p.adjust[mic_in_func,meta_in_func]

  CorrDF <- function(cormat, pmat) {
    ut = matrix (TRUE, nrow = as.numeric(nrow (cormat)), ncol = as.numeric(ncol (pmat)))
    data.frame(
      from = rownames(cormat)[row(cormat)[ut]],
      to = colnames(cormat)[col(cormat)[ut]],
      r =(cormat)[ut],
      p.adjust = pmat[ut]
    )
  }

  mic_meta_cor_df <- CorrDF(mic_meta_cor,mic_meta_cor_p.adjust)
  if(!file.exists("./results/TarNet")){
    dir.create("./results/TarNet",recursive = T)
  }

  write.csv(mic_meta_cor_df,file='./results/TarNet/network_cor_results.csv',row.names = FALSE)

  mic_meta_cor_filtered <- mic_meta_cor_df[which(abs(mic_meta_cor_df$r) > r & mic_meta_cor_df$p.adjust < p_adjust),]

  if(length(mic_meta_cor_filtered$from)==0){
    print("the valus of r is too high,please see the file 'network_cor_results.csv' to set suitable values" )
  }

  centrality = "degree"
  nodeattrib <- data.frame(nodes = union(mic_meta_cor_filtered$from,mic_meta_cor_filtered$to))

  nodeattrib$group <- 0
  for (i in as.character(nodeattrib$nodes)){
    if (i %in% mic_in_func == TRUE){
      nodeattrib[nodeattrib$nodes == i,"group"] <- "Microbe"
    }else if(i %in% meta_in_func == TRUE){
      nodeattrib[nodeattrib$nodes == i,"group"] <- "Metabolite"
    }

  }

  # save file ---------------------------------------------------------------


  write.csv(nodeattrib,file='./results/TarNet/network_nodes_results.csv',row.names = FALSE)

  co_net <- graph_from_data_frame(mic_meta_cor_filtered,direct = F, vertices = nodeattrib)
  co_net=simplify(co_net, remove.multiple=T)
  E(co_net)[ r< 0 ]$color <- "green"
  E(co_net)[ r>0 ]$color <- "red"

  network=visIgraph(co_net)%>%
    visIgraphLayout(layout = "layout_with_fr", smooth = TRUE)%>%
    visNodes(color = list(hover = "green"))%>%visInteraction(hover = TRUE)%>%visLegend(useGroups = T,stepX = 70,stepY=70)%>%
    visGroups(groupname = 'Microbe',color="purple",shape="triangle")%>%
    visGroups(groupname = 'Metabolite',color="tomato",shape="dot")%>%
    visOptions(selectedBy = 'group')
  saveNetwork(network,file = "./results/TarNet/network.html")
  TarnetPairs <- mic_meta_cor_filtered
  return(TarnetPairs)
}
