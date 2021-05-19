library(tidyr)
library(dplyr)
library(foreach)
library(doParallel)
library(corrplot)
library(EnvStats)
library(purrr)
library(stringr)
library(ggplot2)
library(MASS)
library(zoo)

getICVSample <- function(sampleExpression) {
  sampleExpression <- data.frame(t(sampleExpression))
  
  sampleExpression <- foreach(i = 1:ncol(sampleExpression), .combine = cbind) %do% {
    return(sample(sampleExpression[, i], nrow(sampleExpression), replace = T))
  }
  
  sampleICV <- apply(sampleExpression, 1, mean) / apply(sampleExpression, 1, sd)
  return(sampleICV)
}

getCorrelationMatrix <- function(nodeIds, fiberNodes) {
  expression <- rawExpression[nodeIds, -1]
  
  expression <- as.data.frame(t(expression), stringsAsFactors = F)
  colnames(expression) <- rawExpression[nodeIds, 1]
  expression$Experiment = row.names(expression)
  
  expression$Mean = apply(expression[, fiberNodes], 1, mean)
  expression$SD = apply(expression[, fiberNodes], 1, sd)
  expression$ICV = expression$Mean / expression$SD
  expression = expression[!is.na(expression$ICV), ]
  
  if(filteringMethod == 2) {
    sampleICV = getICVSample(sampleExpression = expression[, fiberNodes])
    expression$PValue = 1 - pnorm(expression$ICV, mean = mean(sampleICV), sd = sd(sampleICV))
  }
  
  expression$MishaelZscore = abs(expression$ICV - mean(expression$ICV)) / sd(expression$ICV)
  expression$MishaelPValue = pnorm(-expression$MishaelZscore)
  
  unfilteredExpression = expression
  pvalue = NA
  if(filteringMethod == 1) {
    expression <- expression %>%
      filter(ICV > mean(ICV))
    pvalue = geoMean(expression$MishaelPValue)
  }

  if(filteringMethod == 2) {
    expression <- expression %>%
      filter(ICV > 1) %>%
      filter(PValue < 0.05)
    pvalue = mean(expression$PValue)
  }
  
  expression <- dplyr::select(expression, 1:length(nodeIds))
  
  correlationMatrix <- cor(expression, method = "pearson")
  
  return(list("correlationMatrix" = correlationMatrix, "nexperiments" = nrow(expression), "pvalue" = pvalue))
}

plotCorrelationMatrix <- function(correlationMatrix, id, fiberNodes) {
  return()
  col1 <- colorRampPalette(c("blue", "white", "red"))(100)
  # col1 <- colorRampPalette(c("blue", "white", "white", "red"))(100)
  # col1 <- colorRampPalette(c("blue", "white", "white", "white", "red"))(100)
  p.mat <- as.matrix(1 - abs(correlationMatrix))
  
  svg(filename = paste("NewCorrelationPlots/correlation_", id, ".svg", sep = ""),
      width = 6,
      height = 5,
      pointsize=12)
  corrplot(as.matrix(correlationMatrix),
           method = "color",
           tl.col = "black", tl.srt = 45,
           tl.cex = 1.5,
           cl.cex = 1.5, cl.lim = c(-1, 1),
           cl.length = 6,
           col = col1,
           addgrid.col = "grey",
           p.mat = as.matrix(p.mat),
           sig.level = 1 - 0.6, pch.cex = 1)
  corrRect(fiberNodes, col = "black")
  dev.off()
}

getNodeIdsByName <- function(nodeNames) {
  nodeIds <- foreach(i = 1:length(nodeNames), .combine = c) %do% {
    grep(paste("^", nodeNames[i], "$", sep = ""), rawExpression[, 1], ignore.case = T)
  }
  return(nodeIds)
}

# to use dopar for bigger datasets we need to split blocks in the separate lines and put them together using foreach, don't remove errors
# do this later
getBlockSynchronizationData <- function(blocks, realData) {
  foreach(blockId = 1:nrow(blocks), .errorhandling = "remove") %do% {
    print(paste(blockId, nrow(blocks), sep = "/"))
    
    # make a list of nodes
    blockRegulators <- unlist(strsplit(blocks$Regulators[blockId], split = "; "))
    regulatorNodeIds <- getNodeIdsByName(blockRegulators)
    blockFibers <- unlist(strsplit(blocks$Fiber[blockId], split = "; "))
    fiberNodeIds <- getNodeIdsByName(blockFibers)
    blocks$RegulatorSize[blockId] <- length(regulatorNodeIds)
    blocks$FiberSize[blockId] <- length(fiberNodeIds)
    
    if(length(fiberNodeIds) < 2) {
      blocks$Synchronization[blockId] = -1
      blocks$numberOfConditions[blockId] = -1
      blocks$MeanCorrelation[blockId] = -1
      blocks$MeanRegulatorCorrelation[blockId] = -1
      return()
    }
    
    # get correlation matrix and get off diagonal terms to analyze
    fiberNodes = (length(regulatorNodeIds) + 1):(length(regulatorNodeIds) + length(fiberNodeIds))
    correlationData <- getCorrelationMatrix(nodeIds = c(regulatorNodeIds, fiberNodeIds), fiberNodes = fiberNodes)
    correlationMatrix <- correlationData$correlationMatrix

    blockCorrelation <- correlationMatrix[fiberNodes, fiberNodes]
    blockCorrelation <- blockCorrelation[row(blockCorrelation) != col(blockCorrelation)]
    
    # we are finding mean correlation between regulators and regulators and regulators and fibers to test the null-hypothesis of C(fiber) = 1 & C(regulator) = 0
    if(length(regulatorNodeIds) != 0) {
      regulatorNodeIds <- 1:length(regulatorNodeIds)
      regulatorCorrelation <- matrix(data = correlationMatrix[regulatorNodeIds, ],
                                     nrow = length(regulatorNodeIds), ncol = ncol(correlationMatrix))
      regulatorCorrelation <- regulatorCorrelation[row(regulatorCorrelation) != col(regulatorCorrelation)]
    } else {
      regulatorCorrelation = -1
    }
    
    # if any of block correlations is NA, correlation matrix didn't calculate well, so we assume that there is some problem with conditions and remove this element from count
    if(reduce(is.na(blockCorrelation), function(x, y) x | y)) {
      blocks$Synchronization[blockId] = -1
      blocks$numberOfConditions[blockId] = -1
      blocks$MeanCorrelation[blockId] = -1
      blocks$MeanRegulatorCorrelation[blockId] = -1
      next
    }
    
    if(prod(blockCorrelation > 0.6)) {
      blocks$Synchronization[blockId] = 1
    } else {
      blocks$Synchronization[blockId] = 0
    }
    
    blocks$numberOfConditions[blockId] = correlationData$nexperiments
    blocks$meanPvalue[blockId] = correlationData$pvalue
    blocks$MeanCorrelation[blockId] = mean(blockCorrelation)
    blocks$MeanRegulatorCorrelation[blockId] = mean(regulatorCorrelation)
    
    if(realData) {
      plotCorrelationMatrix(correlationMatrix, blockId, c(length(regulatorNodeIds), length(fiberNodeIds)))
    }
  }
  
  return(blocks)
}

generateRandomBlocks <- function(N, geneNames, uniform = T, maxRegulatorSize = 3, maxFiberSize = 10) {
  blocks = as.data.frame(matrix(0, nrow = N, ncol = 2))
  colnames(blocks) = c("Regulators", "Fiber")
  blocks$Regulators = ""
  
  # generating sizes of regulators and fibers following given ditribution
  if(uniform) {
    regulatorSizes = sample(x = 0:maxRegulatorSize, size = N, replace = T)
    fiberSizes = sample(x = 2:maxFiberSize, size = N, replace = T)
  } else {
    # all these numbers are coming from counting of the building block regulator and fiber sizes
    regulatorSizes = sample(x = 0:3, size = N, replace = T,
                            prob = c(0.198, 0.582, 0.176, 0.0440))
    fiberSizes     = sample(x = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14,
                                  16, 17, 18, 22, 23, 24, 25, 36, 47),
                            size = N, replace = T,
                            prob = c(0.27272727, 0.12500000, 0.11363636, 0.06818182, 0.09090909, 0.04545455,
                                     0.05681818, 0.02272727, 0.02272727, 0.01136364, 0.02272727, 0.03409091,
                                     0.01136364, 0.01136364, 0.01136364, 0.02272727, 0.01136364, 0.01136364,
                                     0.01136364, 0.01136364, 0.01136364))
  }
  
  blocks$Regulators <- foreach(i = 1:N, .combine = c) %do% {
    paste(sample(geneNames, regulatorSizes[i]), collapse = "; ")
  }
  
  blocks$Fiber <- foreach(i = 1:N, .combine = c) %do% {
    paste(sample(geneNames, fiberSizes[i]), collapse = "; ")
  }
  return(blocks)
}

readRealBlocks <- function() {
  if(ecoli == T) {
    blocks = read.table("Blocks/ecoli_blocks_cleaned.txt", sep = "\t", header = T, stringsAsFactors = F)
    # blocks = read.table("Blocks/ecoli_blocks_cleaned_carbon.txt", sep = "\t", header = T, stringsAsFactors = F)
    # blocks = read.table("Blocks/figureBlocks.txt", sep = "\t", header = T, stringsAsFactors = F)
  } else {
    blocks = read.table("Blocks/bacillus_blocks_cleaned.txt", sep = "\t", header = T, stringsAsFactors = F)
  }
  return(blocks)
}

readRawExpression <- function() {
  if(ecoli == T) {
    rawExpression <- read.table("ExpressionData/ECOMICS_WT_2019_2.txt", stringsAsFactors = F, sep = "\t", header = T, nrows = 4097)
  } else {
    rawExpression <- read.table("ExpressionData/Expression_bacillus.txt", stringsAsFactors = F, sep = "\t", header = T)
    rawExpression$ORF = rawExpression$name
    rawExpression <- rawExpression[, -ncol(rawExpression)]
    colnames(rawExpression)[1] <- "GeneName"
  }
  return(rawExpression)
}

#################################
####### Script parameters #######
#################################

filteringMethod = 3  # 1 (Mishael), 2 (Ian), other (no condition filtering)
NRandomBlocks = 10000
ecoli = T

#################################
### Actual script starts here ###
#################################

setwd("/home/ian/Dropbox (City College)/Research/PhD work/shared folders/EXPRESSION-2019/ECOLI_11/Github/")

rawExpression <- readRawExpression()

# # this is real data
blocks <- readRealBlocks()
tic()
blocks <- getBlockSynchronizationData(blocks, T)
toc()

# blocks <- filter(blocks, FiberSize <= 25)
# blocks2 <- filter(blocks2, FiberSize <= 25)
# blocks$nl = blocks2$nl
# blocks2 = read.table("ecoli_blocks.txt", sep = "\t", header = T, stringsAsFactors = F)
# blocks$BlockName <- blocks2$BlockName
# blocks <- read.table("ecoli_blocks_cleaned_no_operons_output.txt", sep = "\t", header = T, stringsAsFactors = F)

# this is random simulation for Ian's plot
# tic()
# randomBlocks <- generateRandomBlocks(N = NRandomBlocks, geneNames = rawExpression[, 1], uniform = T, maxFiberSize = max(blocks$FiberSize))
# randomBlocks <- getBlockSynchronizationData(randomBlocks, F)
# toc()

# write.table(randomBlocks, "PvalueRuns/bacillus_noFilter.txt", quote = F, row.names = F, sep = "\t")
if(ecoli) {
  if(filteringMethod == 1) {
    randomBlocks <- read.table("PvalueRuns/uniformMeanMishael.txt", sep = "\t", header = T, stringsAsFactors = F)
  }
  if(filteringMethod == 3) {
    randomBlocks <- read.table("PvalueRuns/uniformMeanNoFilter.txt", sep = "\t", header = T, stringsAsFactors = F)
  }
} else {
  if(filteringMethod == 1) {
    randomBlocks <- read.table("PvalueRuns/bacillus_ICV.txt", sep = "\t", header = T, stringsAsFactors = F)
  }
  if(filteringMethod == 3) {
    randomBlocks <- read.table("PvalueRuns/bacillus_noFilter.txt", sep = "\t", header = T, stringsAsFactors = F)
  }
}

randomBlocksMeans <- randomBlocks %>%
  filter(!is.na(MeanCorrelation)) %>%
  dplyr::select(c(-1, -2)) %>%
  group_by(FiberSize) %>%
  summarise_all(list(mean = mean, sd = sd))

blocks <- filter(blocks, Synchronization != -1 & FiberSize <= max(randomBlocks$FiberSize))
blocks$RandomCorrelation_mean <-
  unlist(
    foreach(size = blocks$FiberSize, combine = c) %do% {
      filter(randomBlocksMeans, FiberSize == size)$MeanCorrelation_mean
    })
blocks$RandomCorrelation_sd <-
  unlist(
    foreach(size = blocks$FiberSize, combine = c) %do% {
      filter(randomBlocksMeans, FiberSize == size)$MeanCorrelation_sd
    })

# blocks$Significant <- factor(as.integer(blocks$MeanCorrelation > blocks$RandomCorrelation_mean + 1.65 * blocks$RandomCorrelation_sd), levels = c(0, 0.5, 1))
blocks$Significant <- (as.integer(blocks$MeanCorrelation > blocks$RandomCorrelation_mean + 1.65 * blocks$RandomCorrelation_sd) * 2) + 16
blocks$Color = blocks$Significant

approximateMaxAndMin = F
if(approximateMaxAndMin) {
  randomBlocksMeans$MeanCorrelation_sd_max <- randomBlocksMeans$MeanCorrelation_mean + 2 * randomBlocksMeans$MeanCorrelation_sd
  randomBlocksMeans$MeanCorrelation_sd_min <- randomBlocksMeans$MeanCorrelation_mean - 2 * randomBlocksMeans$MeanCorrelation_sd
  
  # plot(log(randomBlocksMeans$FiberSize[1:15]), -log(randomBlocksMeans$MeanCorrelation_sd_max[1:15] - 0.1))
  # plot(log(randomBlocksMeans$FiberSize[1:15]), -log(-randomBlocksMeans$MeanCorrelation_sd_min[1:15]))
  
  model.max <- nls(MeanCorrelation_sd_max ~ beta + FiberSize^(1/alpha),
                   data = randomBlocksMeans[1:15,],
                   start = list(alpha = -1, beta = 0.1))
  model.min <- nls(MeanCorrelation_sd_min ~ - beta * FiberSize^(1/alpha),
                   data = randomBlocksMeans[1:15,],
                   start = list(alpha = -1, beta = 1))
  # summary(model.max)
  # summary(model.min)
  
  ggplot() +
    geom_point(data = randomBlocksMeans, mapping = aes(x = FiberSize, y = MeanCorrelation_sd_max), color = "blue") +
    geom_point(data = randomBlocksMeans, mapping = aes(x = FiberSize, y = MeanCorrelation_sd_min), color = "blue") +
    geom_point(mapping = aes(x = 1:50, y = (1:50) ^ (1/coef(model.max)[1]) + coef(model.max)[2]), color = "red") +
    geom_point(mapping = aes(x = 1:50, y = -coef(model.min)[2] * (1:50) ^ (1/coef(model.min)[1])), color = "red")
  
  ggplot(data = blocks) +
    geom_errorbar(aes(x = round(seq(from = 2, to = 50, length.out = nrow(blocks))),
                      ymin = -coef(model.min)[2] * round(seq(from = 2, to = 50, length.out = nrow(blocks))) ^ (1/coef(model.min)[1]),
                      ymax = round(seq(from = 2, to = 50, length.out = nrow(blocks))) ^ (1/coef(model.max)[1]) + coef(model.max)[2]),
                  colour = "green", width = .1) +
    geom_point(mapping = aes(x = FiberSize, y = MeanCorrelation)) +
    geom_line(data = blocks %>%
                 group_by(FiberSize) %>%
                 summarise(MeanCorrelation = mean(MeanCorrelation)),
               mapping = aes(x = FiberSize, y = MeanCorrelation), color = "blue") +
    geom_line(mapping = aes(x = FiberSize, y = RandomCorrelation_mean), color = "red")
} else {
  png(filename = "significance.png", width = 1800, height = 800)
  ggplot(data = blocks) +
    geom_hline(yintercept = 0, size = 3) +
    geom_errorbar(data = randomBlocksMeans,
                  aes(x = FiberSize,
                      ymin = MeanCorrelation_mean - 1.65 * MeanCorrelation_sd,
                      ymax = MeanCorrelation_mean + 1.65 * MeanCorrelation_sd),
                  colour = "green4", width = .5, size = 5) +
    geom_line(data = blocks %>%
                group_by(FiberSize) %>%
                summarise(MeanCorrelation = mean(MeanCorrelation)),
              mapping = aes(x = FiberSize,
                            y = rollapply(MeanCorrelation, c(1, 3, rep(5, length(MeanCorrelation) - 4), 3, 1),
                                          FUN = mean, fill = NA, align = "center")),
              color = "blue", size = 5) +
    geom_segment(aes(x = 2, y = 0.0685, xend = max(randomBlocks$FiberSize), yend = 0.0685), linetype = "dashed", color = "red", size = 5) +
    geom_line(mapping = aes(x = FiberSize, y = RandomCorrelation_mean), color = "red", size = 5) +
    geom_point(mapping = aes(x = FiberSize, y = MeanCorrelation), shape = blocks$Significant, size = 15, color = blocks$Color) +
    theme(text = element_text(size = 50), legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    xlab("Fiber Size") + ylab("Mean Correlation") +
    scale_x_continuous(breaks = 2:24, limits = c(1, 25), expand = c(0, 0)) +
    scale_y_continuous(breaks = 0:5 * 0.2, limits = c(-0.2, 1.2), expand = c(0, 0))
    # scale_x_continuous(breaks = 2:24, limits = c(1, 25), expand = c(0, 0)) +
    # scale_y_continuous(breaks = -2:3 * 0.2, limits = c(-0.4, 0.7), expand = c(0, 0))
  dev.off()
}

# size = 2
pvalueData <- foreach(size = 2:50, .combine = rbind) %do% {
  experimentalMean = mean(filter(blocks, FiberSize == size)$MeanCorrelation)
  randomMean = filter(randomBlocksMeans, FiberSize == size)$MeanCorrelation_mean
  randomSd = (filter(randomBlocksMeans, FiberSize == size)$MeanCorrelation_sd)
  numberOfBlocks = nrow(filter(blocks, FiberSize == size))
  pvalue = 1 - pnorm(experimentalMean,
            mean = randomMean,
            sd   = randomSd / numberOfBlocks ^ (1/2))
  return(c(size, experimentalMean, randomMean, randomSd, numberOfBlocks, pvalue))
}
pvalueData <- as.data.frame(pvalueData)
colnames(pvalueData) <- c("size", "experimentalMean", "randomMean", "randomSd", "numberOfBlocks", "pvalue")
rownames(pvalueData) <- NULL
pvalueData$pvalue <- round(pvalueData$pvalue, digits = 3)
pvalueData <- filter(pvalueData, !is.na(experimentalMean))
pvalueData$pvalue <- round(pvalueData$pvalue, digits = 3)

pvalueData %>%
  filter(pvalue <= 0.05) %>%
  summarise(SignificantBlocks = sum(numberOfBlocks))

blocks$Zscore <- (blocks$MeanCorrelation - blocks$RandomCorrelation_mean) / blocks$RandomCorrelation_sd
blocks$Pvalue <- round(1 -
                         pnorm(blocks$MeanCorrelation,
                               mean = blocks$RandomCorrelation_mean,
                               sd = blocks$RandomCorrelation_sd),
                       digits = 3)

#################################
#### Actual script ends here ####
#################################