setwd("/mnt/volume1/qyc/")
source("auto-annot/pipelines/scripts/singleR.R")
source("auto-annot/pipelines/scripts/single_cell_net.R")
source("auto-annot/pipelines/evaluate.R")

library(singleCellNet)
library(dplyr)
library(SingleR)
library(Seurat)

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

# Define paths
# train_data_path <- "/mnt/volume1/qyc/data/Intra-dataset/CellBench/10x_5cl/10x_5cl_data.csv"
# train_labels_path <- "/mnt/volume1/qyc/data/Intra-dataset/CellBench/10x_5cl/Labels.csv"
# test_data_path <- "/mnt/volume1/qyc/data/Intra-dataset/CellBench/CelSeq2_5cl/CelSeq2_5cl_data.csv"
# test_labels_path <- "/mnt/volume1/qyc/data/Intra-dataset/CellBench/CelSeq2_5cl/Labels.csv"
train_data_path <- "/mnt/volume1/qyc/data/Inter-dataset/PbmcBench/10Xv3/10Xv3_pbmc1.csv"
train_labels_path <- "/mnt/volume1/qyc/data/Inter-dataset/PbmcBench/10Xv3/10Xv3_pbmc1Labels.csv"
test_data_path <- "/mnt/volume1/qyc/data/Inter-dataset/PbmcBench/CEL-Seq/CL_pbmc1.csv"
test_labels_path <- "/mnt/volume1/qyc/data/Inter-dataset/PbmcBench/CEL-Seq/CL_pbmc1Labels.csv"

# Create output directories
results_dir <- "/mnt/volume1/qyc/auto-annot/result_across"
singler_dir <- file.path(results_dir, "SingleR")
scnet_dir <- file.path(results_dir, "singleCellNet")

dir.create(singler_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(scnet_dir, recursive = TRUE, showWarnings = FALSE)

# Function to run SingleR cross-dataset prediction
run_singleR_cross <- function(train_data_path, train_labels_path, test_data_path, test_labels_path, output_dir) {
  # Read data
  train_data <- read.csv(train_data_path, row.names = 1)
  train_labels <- as.matrix(read.csv(train_labels_path))
  test_data <- read.csv(test_data_path, row.names = 1)
  test_labels <- as.matrix(read.csv(test_labels_path))
  
  # Transpose data for SingleR
  train_data <- t(as.matrix(train_data))
  test_data <- t(as.matrix(test_data))
  
  # Run SingleR
  start_time <- Sys.time()
  singler <- SingleR(method = "single", 
                    test_data, 
                    train_data, 
                    train_labels[,1])
  end_time <- Sys.time()
  
  # Calculate time
  total_time <- as.numeric(difftime(end_time, start_time, units = 'secs'))
  
  # Get predictions
  true_labels <- test_labels[,1]
  pred_labels <- as.vector(singler$labels)
  
  # Save results
  setwd(output_dir)
  write.csv(true_labels, 'SingleR_True_Labels.csv', row.names = FALSE)
  write.csv(pred_labels, 'SingleR_Pred_Labels.csv', row.names = FALSE)
  write.csv(total_time, 'SingleR_Total_Time.csv', row.names = FALSE)
  
  # Evaluate results
  result <- evaluate(true_labels, pred_labels)
  print("SingleR Results:")
  print(result)
}

# Function to run singleCellNet cross-dataset prediction
run_singlecellnet_cross <- function(train_data_path, train_labels_path, test_data_path, test_labels_path, output_dir) {
  # Read data
  train_data <- read.csv(train_data_path, row.names = 1)
  colnames(train_data) <- gsub('_','.',colnames(train_data), fixed = TRUE)
  train_labels <- as.matrix(read.csv(train_labels_path))
  test_data <- read.csv(test_data_path, row.names = 1)
  colnames(test_data) <- gsub('_','.',colnames(test_data), fixed = TRUE)
  test_labels <- as.matrix(read.csv(test_labels_path))
  
  # Transpose data
  train_data <- t(as.matrix(train_data))
  test_data <- t(as.matrix(test_data))
  
  # Training phase
  start_time <- Sys.time()
  cgenes2 <- findClassyGenes(train_data, 
                           data.frame(Annotation = train_labels[,1]), 
                           "Annotation")
  cgenesA <- cgenes2[['cgenes']]
  grps <- cgenes2[['grps']]
  train_data <- as.matrix(train_data[cgenesA,])
  xpairs <- getTopGenePairs(train_data, grps)
  pdTrain <- query_transform(train_data[cgenesA, ], xpairs)
  rf <- sc_makeClassifier(pdTrain[xpairs,], genes=xpairs, groups=grps)
  end_time <- Sys.time()
  training_time <- as.numeric(difftime(end_time, start_time, units = 'secs'))
  
  # Testing phase
  start_time <- Sys.time()
  test_data <- query_transform(test_data[cgenesA,], xpairs)
  classRes <- rf_classPredict(rf, test_data)
  end_time <- Sys.time()
  testing_time <- as.numeric(difftime(end_time, start_time, units = 'secs'))
  
  # Get predictions
  true_labels <- test_labels[,1]
  pred_labels <- rownames(classRes)[apply(classRes,2,which.max)]
  
  # Save results
  setwd(output_dir)
  write.csv(true_labels, 'singleCellNet_True_Labels.csv', row.names = FALSE)
  write.csv(pred_labels, 'singleCellNet_Pred_Labels.csv', row.names = FALSE)
  write.csv(training_time, 'singleCellNet_Training_Time.csv', row.names = FALSE)
  write.csv(testing_time, 'singleCellNet_Testing_Time.csv', row.names = FALSE)
  
  # Evaluate results
  result <- evaluate(true_labels, pred_labels)
  print("singleCellNet Results:")
  print(result)
}

# Run predictions
print("Running SingleR cross-dataset prediction...")
run_singleR_cross(train_data_path, train_labels_path, test_data_path, test_labels_path, singler_dir)

# print("Running singleCellNet cross-dataset prediction...")
# run_singlecellnet_cross(train_data_path, train_labels_path, test_data_path, test_labels_path, scnet_dir) 