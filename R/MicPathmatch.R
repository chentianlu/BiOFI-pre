#'@title MicPathmatch
#'@description
#' A function used to predict the functional pathway of microbes in the 16S or Metagenome data.
#' You can get two forms using this function to predict the possible next level strains and potental pathways of microbes you are interested in.
#' One is a form including predicted detail next level strains and pathways of your own microbes,
#' The other is a form including only your own microbes and corresponding pathways.
#' The format of microbe name is like this:
#' k__XXX; p__XXX; c__XXX; o__XXX; f__XXX; g__XXX; s__XXX;
#' k means kingdom;
#' p means phylum;
#' c means class;
#' o means order;
#' f means family;
#' g means genus;
#' s means species.
#'
#' @param Microbelist a list of microbe names in 16S or Metagenome data, and it does not contain the names of functional information in 16S or Metagenome data.
#'
#' @return Microberesult a list that contains two parts, one is a data frame that contains pathway and strains or sub-strains of microbes in the Microbelist
#'         and the other is a data frame that only contains Microbelist names that need to predict and relative pathways
#' @export MicPathmatch
#' @examples
#' Micrest=MicPathmatch(Microbe_path_example)

MicPathmatch<-function(Microbelist){
  microbename=Microbelist
  if(!length(microbename[grep("; ",microbename)])==length(microbename)){

    print("the list of microbes may contain inappropriate format,please check the format of microbes.")
    print("the format of microbes is 'k__XXX; p__XXX; c__XXX; o__XXX; f__XXX; g__XXX; s__XXX;'")
  }
  microbelist=c()
  B=data.frame()
  for(i in 1:length(microbename)){
    microbe=strsplit(microbename[i],'; ')
    microbenew=vector()
    for(i in 1:length(microbe[[1]])){
      if(nchar(microbe[[1]][i])==3){
        numb=i-1
        break}else{
          numb=length(microbe[[1]])
        }}
    for(i in 1:numb){
      a=microbe[[1]][i]
      microbenew=append(microbenew,a)
    }
    microbenew=paste(microbenew,collapse ='; ')
    microbelist=append(microbelist,microbenew)
  }

  for(i in 1:length(microbelist)){
    micname=strsplit(microbelist[i],'; ')
    num=length(micname[[1]])
    if(num==3){
      if(micname[[1]][2]=='p__unclassified Bacteria'){
        A=pathlist[pathlist$Class=='c__Abyssogena phaseoliformis symbiont [TAX:596095]',]

      }else{
        A=pathlist[pathlist$Class==micname[[1]][3],]
      }
      if(length(rownames(A)!=0)){
        A$name=microbelist[i]}
    }
    if(num==4){
      A=pathlist[pathlist$Order==micname[[1]][4],]
      if(length(rownames(A)!=0)){
        A$name=microbelist[i]}
    }
    if(num==5){
      A=pathlist[pathlist$Family==micname[[1]][5],]
      if(length(rownames(A)!=0)){
        A$name=microbelist[i]}
    }
    if(num==6){
      micname6=gsub('g__','',micname[[1]][6])
      A41=pathlist[pathlist$Family==micname[[1]][5],]

      A=A41[grepl(micname6,A41$Species),]
      if(length(rownames(A)!=0)){
        A$name=microbelist[i]}

    }
    if(num==7){
      micname7=gsub('s__','',micname[[1]][7])
      A51=pathlist[pathlist$Order==micname[[1]][4],]
      A52=A51[A51$Family==micname[[1]][5],]
      A53=pathlist[pathlist$Genus==micname[[1]][6],]
      A=pathlist[grepl(micname7,A53$Species),]
      if(length(rownames(A)!=0)){
        A$name=microbelist[i]}
    }

    if(is.null(A)){
      A=data.frame()
    }

    B=rbind(B,A)
  }
  microbetotle=B
  if(!file.exists("./results/Micpathmatch")){
    dir.create("./results/Micpathmatch",recursive = T)
  }
  B1=unique(B[,c(11,10)])
  microbepredict=B1
  write.csv(microbetotle,file='./results/Micpathmatch/Pathways of microbes(including prediction of strain and substrain).csv',row.names = FALSE)
  write.csv(microbepredict,file='./results/Micpathmatch/Pathways of microbes.csv',row.names=FALSE)
  Microberesult=list('Pathways of microbes(detail)'=microbetotle,'Pathways of microbes'=microbepredict)
  return(Microberesult)
}
