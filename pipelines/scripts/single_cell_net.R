csRenameOrth<-function(
  expQuery,
  expTrain,
  orthTable,
  speciesQuery='human',
  speciesTrain='mouse'){

  # what genes in the orth table are in the query table?
  rownames(orthTable)<-as.vector(orthTable[,speciesQuery])
  cgenes<-intersect(rownames(expQuery), rownames(orthTable))
  cat("query genes in ortholog table = ",length(cgenes),"\n")

  # of these, which are in training data?
  oTab<- orthTable[cgenes,]
  rownames(oTab) <- as.vector(oTab[,speciesTrain])
  ccGenes<- intersect(rownames(expTrain), rownames(oTab))
  cat("training genes in ortholog table and query data = ",length(ccGenes),"\n")
  # should put a check here for sufficient number of genes 
  oTab<-oTab[ccGenes,]
  qGenes<-as.vector(oTab[,speciesQuery])
  expQuery <- expQuery[qGenes,]
  rownames(expQuery) <- rownames(oTab) 
  expTrain <- expTrain[rownames(oTab),]
  list(expQuery=expQuery, expTrain=expTrain)
}


#' find best pairs
#'
#' find best pairs
#'
#' @param expDat expDat
#' @param cellLabels named vector of cell groups
#'
#' @return vector of pairs
#' 
#' @export
gnrBP<-function(
  expDat,
  cellLabels,
  topX=50){

  ans<-vector()
    myPatternG<-sc_sampR_to_pattern(as.character(cellLabels))
    for(i in seq(length(myPatternG))){
      cat(i,"\n")
      xres<-sc_testPattern(myPatternG[[i]], expDat=expDat)
      tmpAns<-findBestPairs(xres, topX)
      ans<-append(ans, tmpAns)
    }
    unique(ans)
}


#' @title
#' Find Classifier-Worthy Gene Candidates
#' @description
#' Find classifier-worthy genes for each cancer category to train the classifier
#'
#' @param expDat a matrix of normalized expression data from \code{\link{trans_prop}}
#' @param sampTab a dataframe of the sample table
#' @param dLevel a string indicating the column name in sample table that contains the cancer category
#' @param topX an integer indicating the number of top positive classification genes for each category to select for training. Will also select topX number of negative classification genes.
#' @param dThresh a number representing the detection threshold
#' @param alpha1 a number representing proportion of cells in which a gene must be considered detected (as defined in geneStats)
#' @param alpha2 a number representing lower proportion of cells for genes that must have higher expression level
#' @param mu a number represeting threshold for average expression level of genes passing the lower proportion criteria
#'
#' @return a list containing two lists: a list of classifier worthy genes named 'cgenes' and a list of cancer category named 'grps'
#' @export
findClassyGenes<-function
(expDat,
sampTab,
dLevel,
topX=25,
dThresh=0,
alpha1=0.05,
alpha2=.001,
mu=2)
{
  if((dLevel %in% colnames(sampTab)) == FALSE) {
      stop("Please enter the correct column name for sampTab that indicates the categories")
    }

  if((topX * 2 > nrow(expDat)) == TRUE) {
      stop(paste0("Please enter a topX value smaller than ", as.integer(nrow(expDat) / 2)))
    }

  gsTrain<-sc_statTab(expDat, dThresh=dThresh)
  ggenes<-sc_filterGenes(gsTrain, alpha1=alpha1, alpha2=alpha2, mu=mu)
  grps<-as.vector(sampTab[,dLevel])
  names(grps)<-rownames(sampTab)
  xdiff<-gnrAll(expDat[ggenes,], grps)
  cgenes<-lapply(xdiff, getClassGenes, topX=topX)
  cgenes2<-unique(unlist(cgenes))
  list(cgenes=cgenes2, grps=grps, cgenes_list = cgenes)
}

#' @title
#' Find Classy Genes
#' @description
#' Extract genes suitable for training classifier
#' @param diffRes a dataframe with pval, cval, holm, and rownames as the gene names
#' @param topX a number dicataing the number of genes to select for training classifier
#' @param bottom logic if true use the top x genes with - cvals
#'
#' @return a vector of genes that are good for training classifier for that category
getClassGenes<-function(diffRes, topX=25, bottom=TRUE) {
  #exclude NAs
  xi<-which(!is.na(diffRes$cval))
  diffRes<-diffRes[xi,] # exclude the NAs. Select the rows that does not have NAs

  diffRes<-diffRes[order(diffRes$cval, decreasing=TRUE),] #order based on classification value largest to lowest
  ans<-rownames(diffRes[1:topX,]) # get the top 20 genes

  if(bottom){
    ans<-append(ans, rownames( diffRes[nrow(diffRes) - ((topX-1):0),]))
  }

  # get the least differentially expressed genes as house holders
  # sameRes<-diffRes[order(abs(diffRes$cval), decreasing = FALSE), ]
  # ans<-append(ans, rownames(sameRes[1:topX, ]))

  #return
  ans
}

#' find genes higher in a cluster compared to all other cells
#'
#' ind genes higher in a cluster compared to all other cells
#'
#' @param expDat expDat
#' @param cellLabels named vector of cell groups
#'
#' @return list of diffExp data framnes
#' 
#' @export
gnrAll<-function(
  expDat,
  cellLabels){

  myPatternG<-sc_sampR_to_pattern(as.character(cellLabels))
  # sparse matrix?
  if(class(expDat)[1]!='matrix'){
    expTrans = Matrix::t(expDat)
  }
  else{
    expTrans = t(expDat)
  }
  specificSets<-lapply(myPatternG, sc_testPatternTrans, expDat=expTrans)
  cat("Done testing\n")

#  grpOrder<-myGrpSort(cellLabels)

#  specificSets[grpOrder]

  specificSets
}


makePairTab<-function(genes){
  pTab<-t(combn(genes, 2))
  colnames(pTab)<-c("genes1", "genes2")
  pTab<-cbind(pTab, pairName=paste(pTab[,1], "_",pTab[,2], sep=''))
  pTab
}


#' @title
#' Find the best gene pairs for training
#' @description
#' Find the gene pairs that most distinguish a cancer group from the rest
#'
#' @param expDat expDat
#' @param cell_labels named vector, value is grp, name is cell name
#' @param cgenes_list the list of labelled cgenes
#' @param topX number of genepairs for training
#' @param sliceSize the size of the slice for pair transform. Default at 5e3
#' @param quickPairs TRUE if wanting to select the gene pairs in a quick fashion
#'
#' @import parallel
#' @return vector of top gene-pair names
#'
#' @export
ptGetTop <-function(expDat, cell_labels, cgenes_list=NA, topX=50, sliceSize = 5e3, quickPairs = TRUE){
  if(!quickPairs){
    ans<-vector()
    genes<-rownames(expDat)
    
    ncores<-parallel::detectCores() # detect the number of cores in the system
    mcCores<-1
    if(ncores>1){
      mcCores<-ncores - 1
    }
    cat(ncores, "  --> ", mcCores,"\n")
    
    # make a data frame of pairs of genes that will be sliced later
    pairTab<-makePairTab(genes)
    
    if(topX > nrow(pairTab)) {
      stop(paste0("The data set has ", nrow(pairTab), " total combination of gene pairs. Please select a smaller topX."))
    }
    
    # setup tmp ans list of sc_testPattern
    cat("setup ans and make pattern\n")
    grps<-unique(cell_labels)
    myPatternG<-sc_sampR_to_pattern(as.character(cell_labels))
    statList<-list()
    for(grp in grps){
      statList[[grp]]<-data.frame()
    }
    
    # make the pairedDat, and run sc_testPattern
    cat("make pairDat on slice and test\n")
    nPairs = nrow(pairTab)
    cat("nPairs = ",nPairs,"\n")
    str = 1
    stp = min(c(sliceSize, nPairs)) # detect what is smaller the slice size or npairs
    
    while(str <= nPairs){
      if(stp>nPairs){
        stp <- nPairs
      }
      cat(str,"-", stp,"\n")
      tmpTab<-pairTab[str:stp,]
      tmpPdat<-ptSmall(expDat, tmpTab)
      
      if (Sys.info()[['sysname']] == "Windows") {
        tmpAns<-lapply(myPatternG, sc_testPattern, expDat=tmpPdat)
      }
      else {
        tmpAns<-parallel::mclapply(myPatternG, sc_testPattern, expDat=tmpPdat, mc.cores=mcCores) # this code cannot run on windows
      }
      
      for(gi in seq(length(myPatternG))){
        grp<-grps[[gi]]
        statList[[grp]]<-rbind( statList[[grp]],  tmpAns[[grp]])
      }
      
      
      str<-stp+1
      stp<-str + sliceSize - 1
    }
    
    cat("compile results\n")
    for(grp in grps){
      tmpAns<-findBestPairs(statList[[grp]], topX)
      ans<-append(ans, tmpAns)
    }
    return(unique(ans))
    
  }else{
    myPatternG<-sc_sampR_to_pattern(as.character(cell_labels))
    ans<-vector()
    
    for(cct in names(cgenes_list)){
      genes<-cgenes_list[[cct]]
      pairTab<-makePairTab(genes)
      
      nPairs<-nrow(pairTab)
      cat("nPairs = ", nPairs," for ", cct, "\n")
      
      tmpPdat<-ptSmall(expDat, pairTab)
      
      tmpAns<-findBestPairs( sc_testPattern(myPatternG[[cct]], expDat=tmpPdat), topX)
      ans<-append(ans, tmpAns)
    }
    
    return(unique(ans))
  }
}
ptSmall<-function
(expDat,
pTab){
  npairs = nrow(pTab)
  ans<-matrix(0, nrow=npairs, ncol=ncol(expDat))
  genes1<-as.vector(pTab[,"genes1"])
  genes2<-as.vector(pTab[,"genes2"])

    for(i in seq(nrow(pTab))){
      #cat(genes1[i], ": ", genes2[i],"\n")
      ans[i,]<-as.numeric(expDat[genes1[i],]>expDat[genes2[i],]) 
    }
  colnames(ans)<-colnames(expDat)
  rownames(ans)<-as.vector(pTab[,"pairName"])
  ans
}


#' makes complete gene-to-gene comparison
#'
#' @param expDat expDat
#'
#' @return list of list(genes=data frame, expDat=binary matrix)
#'
#' @export
pair_transform<-function # convert to a vector of length = length(vect)^2 - 1 /2
(expDat){
  ngenes<-nrow(expDat)
  genes<-rownames(expDat)
  ans<-matrix(0, nrow=ngenes*(ngenes-1)/2, ncol=ncol(expDat))
  pair_index<-1
  genes1<-vector()
  genes2<-vector()
  for(i in 1:ngenes){
    for(j in 1:ngenes){
      if(j>i){
        genes1<-append(genes1, genes[i])
        genes2<-append(genes2, genes[j])
        ans[pair_index,]<-as.numeric(expDat[i,]>expDat[j,]) 
        pair_index<-pair_index +1
      }
    }
  }
  colnames(ans)<-colnames(expDat)
  tList2 <- list(genes=data.frame(g1=genes1, g2=genes2), tDat=ans)

  pairNames<-paste(tList2[[1]][,1], "_",tList2[[1]][,2], sep='')

  pairDat<-tList2[[2]]
  rownames(pairDat)<-pairNames
  pairDat
}


#' makes complete gene-to-gene comparison
#'
#' @param expDat expDat
#' @param genePairs genePairs
#'
#' @return matrix indicating which gene of a pair is greater
#'
#' @export
query_transform<-function # convert to a vector of length = length(vect)^2 - 1 /2
(expDat,
genePairs #vector of strings indicating pairs to compare
){
  genes<-strsplit(genePairs, "_")
    ans<-matrix(0, nrow=length(genes), ncol=ncol(expDat))
  pair_index<-1
  genes1<-vector()
  genes2<-vector()
  for(i in seq(length(genes))){
    ans[i,]<-as.numeric(expDat[genes[[i]][1],]>expDat[genes[[i]][2],])
  }
  colnames(ans)<-colnames(expDat)
  rownames(ans)<-genePairs
  ans
}

#' finds the best pairs to use
#'
#' @param xdiff xdiff
#' @param n number of pairs
#' @param maxPer indicates the number of pairs that a gene is allowed to be in
#'
#' @return vector of good pairs
#'
#' @export
findBestPairs<-function # find best and diverse set of pairs
(xdiff,
n=50,
maxPer=3){

  
  xdiff<-xdiff[order(xdiff$cval, decreasing=TRUE),]
  genes<-unique(unlist(strsplit(rownames(xdiff), "_")))
  countList <- rep(0, length(genes))
  names(countList) <- genes	

  i<-0
  ans<-vector()
  xdiff_index<-1
  pair_names<-rownames(xdiff)
  while(i < n ){
    tmpAns<-pair_names[xdiff_index]
    tgp <- unlist(strsplit(tmpAns, "_"))
    if( (countList[ tgp[1] ] < maxPer) & (countList[ tgp[2] ] < maxPer )){
      ans<-append(ans, tmpAns)
      countList[ tgp[1] ] <- countList[ tgp[1] ]+ 1
      countList[ tgp[2] ] <- countList[ tgp[2] ]+ 1
      i<-i+1
    }
    xdiff_index <- xdiff_index + 1
  }
  ans
}


#' @export
addRandToSampTab<-function(classRes, sampTab, desc, id="cell_name"){
  cNames<-colnames(classRes)
  snames<-rownames(sampTab)

  rnames<-setdiff(cNames, snames)
  cat("number of random samples: ",length(rnames), "\n")

  stNew<-data.frame(rid=rnames, rdesc=rep("rand", length(rnames)))
  stTop<-sampTab[,c(id, desc)]
  colnames(stNew)<-c(id, desc)
  ans<-rbind(stTop, stNew)
  rownames(ans)<-colnames(classRes)
  ans
}

getTopGenePairs <- function(expr, labels, topX = 50, verbose = TRUE) {
  myPatternG <- sc_sampR_to_pattern(as.character(labels))
  genes <- rownames(expr)

  ans <- vector()

  for (cct in names(myPatternG)) {
    if (verbose) cat("Processing group:", cct, "\n")
    
    # 为每个 group 使用相同的全基因集组合
    pairTab <- makePairTab(genes)
    tmpPdat <- ptSmall(expr, pairTab)
    
    testRes <- sc_testPattern(myPatternG[[cct]], expDat = tmpPdat)
    topPairs <- findBestPairs(testRes, topX)
    ans <- append(ans, topPairs)
  }

  return(unique(ans))
}

run_singleCellNet<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run singleCellNet
  Wrapper script to run singleCellNet on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.csv(DataPath,row.names = 1)
  colnames(Data) <- gsub('_','.',colnames(Data), fixed = TRUE)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                              singleCellNet                                #
  #############################################################################
  library(singleCellNet)
  library(dplyr)
  True_Labels_singleCellNet <- list()
  Pred_Labels_singleCellNet <- list()
  Training_Time_singleCellNet <- list()
  Testing_Time_singleCellNet <- list()
  Data = t(as.matrix(Data))              # deals also with sparse matrix
  
  for(i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      DataTrain <- Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]]
      DataTest <- Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]]
    }
    else{
      DataTrain <- Data[,Train_Idx[[i]]]
      DataTest <- Data[,Test_Idx[[i]]]
    }
    
    start_time <- Sys.time()
    cgenes2<-findClassyGenes(DataTrain, data.frame(Annotation = Labels[Train_Idx[[i]]]), "Annotation")
    cgenesA<-cgenes2[['cgenes']]
    grps<-cgenes2[['grps']]
    DataTrain<-as.matrix(DataTrain[cgenesA,])
    xpairs<-getTopGenePairs(DataTrain, grps)
    pdTrain<-query_transform(DataTrain[cgenesA, ], xpairs)
    rf<-sc_makeClassifier(pdTrain[xpairs,], genes=xpairs, groups=grps)
    end_time <- Sys.time()
    Training_Time_singleCellNet[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    start_time <- Sys.time()
    DataTest<-query_transform(DataTest[cgenesA,], xpairs)
    classRes <-rf_classPredict(rf, DataTest)
    end_time <- Sys.time()
    Testing_Time_singleCellNet[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_singleCellNet[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_singleCellNet[i] <- list((rownames(classRes)[apply(classRes,2,which.max)])[1:length(Test_Idx[[i]])])
  }
  True_Labels_singleCellNet <- as.vector(unlist(True_Labels_singleCellNet))
  Pred_Labels_singleCellNet <- as.vector(unlist(Pred_Labels_singleCellNet))
  Training_Time_singleCellNet <- as.vector(unlist(Training_Time_singleCellNet))
  Testing_Time_singleCellNet <- as.vector(unlist(Testing_Time_singleCellNet))
  
  setwd(OutputDir)
  
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_singleCellNet,paste('singleCellNet_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_singleCellNet,paste('singleCellNet_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Training_Time_singleCellNet,paste('singleCellNet_',NumGenes,'_Training_Time.csv', sep = ''),row.names = FALSE)
    write.csv(Testing_Time_singleCellNet,paste('singleCellNet_',NumGenes,'_Testing_Time.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_singleCellNet,'singleCellNet_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_singleCellNet,'singleCellNet_Pred_Labels.csv',row.names = FALSE)
    write.csv(Training_Time_singleCellNet,'singleCellNet_Training_Time.csv',row.names = FALSE)
    write.csv(Testing_Time_singleCellNet,'singleCellNet_Testing_Time.csv',row.names = FALSE)
  }
}