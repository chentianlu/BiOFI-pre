#' @title IPS
#' @description
#' A function used to calculate the importance pair score (IPS) of a network based on the return of the TarNet function.
#'
#' @param TarnetPairs the results of microbial-metabolic pairs from TarNet function.
#' @param IIS the results from Nodescore function.
#'
#' @export IPS
#' @examples
#' IIScore=IIS(Microbe_example,Meta_example,confounderData,groupInfo)
#' meta_mic2func=MMfunc(IIScore)
#' TarnetPairs=TarNet(IIScore,meta_mic2func)
#' IPSres=IPS(TarnetPairs,IIScore)

IPS<-function(TarnetPairs,IIS){
  IISdf1 = IIS[["Node_Score"]][,1:2]
  colnames(IISdf1) <- c("to","to_IIS")
  IISdf2 = IIS[["Node_Score"]][,1:2]
  colnames(IISdf2) <- c("from","from_IIS")
  IPSres <- merge(TarnetPairs,IISdf1,all.x = T)
  IPSres <- merge(IPSres,IISdf2,all.x = T)
  IPSres$r_score <- abs(IPSres$r)
  IPSres$r_score <-  seq(nrow(IPSres),1,-1)
  IPSres$r_score <- (IPSres$r_score - 1)/(nrow(IPSres) - 1) * 2 + 1
  IPSres$ips <- IPSres$to_IIS + IPSres$from_IIS + IPSres$r_score
  IPSres <- IPSres[order(IPSres$ips,decreasing = T),]
  return(IPSres)
}
