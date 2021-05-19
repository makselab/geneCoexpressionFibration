library(tidyr)
library(dplyr)
library(foreach)
library(corrplot)

# group id is used only for the purpose of making squares in conditions to use fiber conditions for regulators inside the fiber
# parameters
ecoli = T
useOverlapConditions = F
plotMixedCorrelation = F
inputFileId = 6

# folder
setwd("/home/ian/Dropbox (City College)/Research/PhD work/shared folders/EXPRESSION-2019/ECOLI_11/Github/")

getConditions <- function(fiberId, geneName) {
  # function to get conditions for the gene by fiber id and name
  # in case fiberId = -1, we use all possible conditions for the system
  if(fiberId == -1) {
    rawData <- systemConditions[, grep(geneName, colnames(systemConditions))]
    rawData <- data.frame(cbind(rownames(systemConditions), rawData), stringsAsFactors = F)
    rawData$rawData <- as.numeric(rawData$rawData)
    colnames(rawData) <- c("X", geneName)
  } else {
    if(ecoli) {
      rawData <- read.csv(paste("Blocks/Ecoli_WT/", fiberId, " _WT_conditions.csv", sep = ""), stringsAsFactors = F)
    } else {
      rawData <- read.csv(paste("Blocks/Bacilus_rel/", fiberId, " _conditions.csv", sep = ""), stringsAsFactors = F)
    }
    rawData <- rawData[, -c((ncol(rawData) - 2) : ncol(rawData))]
    columnId <- grep(geneName, colnames(rawData))
    if(length(columnId) == 0) {
      return(NULL)
    }
    rawData <- rawData[, c(1, columnId)]
    rawData <- arrange(rawData, X)
  }
  return(rawData)
}

getConditionTable <- function(genes) {
  # building table with fiber ids to choose conditions
  conditionTable <- expand.grid(x = genes$FiberId, y = genes$FiberId)
  conditionTable <- apply(conditionTable, 1, function(x) paste(x, collapse = "/"))
  
  conditionTable <- split(conditionTable, ceiling(1:length(conditionTable)/nrow(genes)))
  conditionTable <- foreach(i = 1:length(conditionTable), .combine = rbind) %do% {return(unlist(conditionTable[i]))}
  
  colnames(conditionTable) <- genes$Gene
  rownames(conditionTable) <- genes$Gene
  
  # now we want to make sure that regulators have fiber conditions inside the fiber
  systemFibers <- genes %>%
    group_by(GroupId) %>%
    summarise(FiberId = first(FiberId), Size = n()) %>%
    dplyr::select(2, 3) %>%
    ungroup()

  startId <- 1
  endId <- systemFibers$Size[1]

  i <- 1
  while(1) {
    conditionTable[startId:endId, startId:endId] <- paste(systemFibers$FiberId[i], systemFibers$FiberId[i], sep = "/")
    if(i == nrow(systemFibers)) {break}
    startId <- startId + systemFibers$Size[i]
    endId <- startId + systemFibers$Size[i + 1] - 1
    i <- i + 1
  }
  
  return(conditionTable)
}

getConditionsOfFullSystem <- function(rawExpression, genes) {
  # returns all conditions for all fibers in the system to use for regulators
  systemFiberIds <- unique(genes$FiberId)
  systemFiberIds <- systemFiberIds[systemFiberIds != -1]
  
  if(length(systemFiberIds) == 0) {
    systemConditions = colnames(rawExpression)
  } else {
    systemConditions <- foreach(i = systemFiberIds, .combine = c) %do% {
      fiberGene <- first(genes$Gene[genes$FiberId == i])
      fiber_conditions <- getConditions(i, fiberGene)
      return(fiber_conditions$X)
    }
  }
  systemConditions <- unique(systemConditions)
  systemConditions <- colnames(rawExpression) %in% systemConditions
  systemConditions <- rawExpression[, systemConditions]
  rownames(systemConditions) <- rawExpression[, 1]
  systemConditions <- as.data.frame(t(systemConditions))
}

getRawExpression <- function() {
  if(ecoli) {
    rawExpression <- read.table("ExpressionData/ECOMICS_WT_2019_2.txt", stringsAsFactors = F, sep = "\t", header = T)
  } else {
    rawExpression <- read.table("ExpressionData/Expression_bacillus.txt", stringsAsFactors = F, sep = "\t", header = T)
    rawExpression$ORF <- rawExpression$name
    colnames(rawExpression)[1] <- "GeneName"
  }
  return(rawExpression)
}

rawExpression <- getRawExpression()

if(ecoli) {
  genes <- read.table(paste("CrossCorrelationCodeParameters/Ecoli", inputFileId, ".txt", sep = ""), header = T, sep = "\t", stringsAsFactors = F)
} else {
  genes <- read.table(paste("CrossCorrelationCodeParameters/Bacillus", inputFileId, ".txt", sep = ""), header = T, sep = "\t", stringsAsFactors = F)
}

genes <- genes[genes$Gene %in% rawExpression[, 1], ]

conditionTable <- getConditionTable(genes)

# get all conditions for all fibers to use for regulators
if(!useOverlapConditions) {
  systemConditions <- getConditionsOfFullSystem(rawExpression, genes)
}

correlation <- foreach(i = 1:nrow(genes), .combine = rbind, .errorhandling = 'remove') %do% {
  foreach(j = 1:nrow(genes), .combine = cbind, .errorhandling = 'remove') %do% {
    # i = 1
    # j = 2
    xGeneName <- colnames(conditionTable)[j]
    xGeneFiberId <- unlist(strsplit(conditionTable[i, j], split = "/"))[1]
    
    yGeneName <- rownames(conditionTable)[i]
    yGeneFiberId <- unlist(strsplit(conditionTable[i, j], split = "/"))[2]
    
    conditionsX <- getConditions(fiberId = xGeneFiberId, geneName = xGeneName)
    conditionsY <- getConditions(fiberId = yGeneFiberId, geneName = yGeneName)
    if(is.null(conditionsX) | is.null(conditionsY)) {stop(paste("Genes", xGeneName, "and", yGeneName, "are skipped for the lack of conditions", sep = " "))}
    if(xGeneFiberId == yGeneFiberId) {
      correlationValue = cor(conditionsX[, 2], conditionsY[, 2], method = "pearson")
    } else {
      if(!useOverlapConditions) {
        overlapConditions <- unique(c(conditionsX$X, conditionsY$X))
        overlapConditions <- colnames(rawExpression) %in% overlapConditions
        
        conditionsX <- t(rawExpression[grep(paste("^", xGeneName, "$", sep = ""), rawExpression$GeneName), overlapConditions])
        conditionsY <- t(rawExpression[grep(paste("^", yGeneName, "$", sep = ""), rawExpression$GeneName), overlapConditions])
        
        correlationValue <- as.double(cor(conditionsX, conditionsY, method = "pearson"))
      } else {
        overlapConditions <- conditionsX$X[conditionsX$X %in% conditionsY$X]
        conditionsX <- conditionsX[conditionsX$X %in% overlapConditions, ]
        conditionsY <- conditionsY[conditionsY$X %in% overlapConditions, ]
        
        correlationValue = cor(conditionsX[, 2], conditionsY[, 2], method = "pearson")
      }
    }
    
    return(correlationValue)
  }
}

colnames(correlation) <- genes$Gene
rownames(correlation) <- genes$Gene

# col1 <- colorRampPalette(c("blue", "white", "red"))(100)
col1 <- colorRampPalette(c("blue", "white", "white", "red"))(100)
col1 <- colorRampPalette(c("blue", "white", "white", "white", "red"))(100)
p.mat <- as.matrix(1 - abs(correlation))

# svg(filename="composite.svg",
#     width = 5,
#     height = 5,
#     pointsize=12)

corrplot(as.matrix(correlation),
         # type = "lower",
         # type = "upper",
         method = "color",
         tl.col = "black", tl.srt = 45,
         tl.cex = 1.5,
         cl.cex = 1.5, cl.lim = c(-1, 1),
         cl.length = 6,
         col = col1,
         addgrid.col = "grey",
         # p.mat = as.matrix(p.mat),
         # sig.level = 1 - 0.6, pch.cex = 1
         )#, insig = "blank")

corrRect(table(genes$GroupId), col = "black")
# dev.off()