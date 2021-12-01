### Initialisation of the system ----

library(tidyverse)
library(compiler)
library(rlang)
library(ggstatsplot)

rngFunc <- function(randDist)
{
    f <- function(n, rng)
    {
        runif(n)*(rng[2]-rng[1]) + rng[1]
    }
    if (randDist == "normal")
    {
        f <- function(n, rng)
        {
            m <- sum(rng)/2
            s <- (rng[2]-rng[1])/5
            rnorm(n, m, s)
        }
    }
    if (randDist == "exp")
    {
        f <- function(n, rng)
        {
            rexp(n)*(rng[2]-rng[1]) + rng[1]
        }
    }
    if (randDist == "poisson")
    {
        f <- function(n, rng)
        {
            p <- rpois(n, mean(x)*10%>%round)
            p*(rng[2]-rng[1]) + rng[1]
        }
    }
    return(f)
}

initFunc <- function(n = 300, nhalf = 100, nFactor = 100, tMax = 1000, 
                     pUnderRange = c(0,0.3), pOverRange = c(0,0.8),
                     pnormRange = c(0,0.5), randDist = "unif",
                     onThresh = 400, pBase = NULL)
{
    mn <- global_env()
    
    UnderGenes <- paste0("Down", 1:nhalf)
    OverGenes <- paste0("Up", 1:nhalf)
    NormalGenes <- paste0("Norm", 1:nhalf)
    
    # pUnderRange <- c(0, 0.3)
    # pOverRange <- c(0.0, 0.8)
    # pnormRange <- c(0,0.5)
    
    randFunc <- rngFunc(randDist)
        
    
    pUnder <- randFunc(nhalf, pUnderRange)
    pOver <- randFunc(nhalf, pOverRange)
    pNorm <- randFunc(nhalf, pnormRange)
    
    positions <- sample(1:n, n)
    
    positionsUnder <- positions[1:nhalf]
    positionsOver <- positions[(nhalf+1):(2*nhalf)]
    positionsNormal <- positions[(2*nhalf + 1):n]
    
    names(UnderGenes) <- positionsUnder
    names(OverGenes) <- positionsOver
    names(NormalGenes) <- positionsNormal
    
    genes <- c(UnderGenes, OverGenes, NormalGenes)
    or <- order(as.integer(names(genes)))
    genes <- genes[or]
    mn$pBase <- pBase
    if (is.null(pBase))
        mn$pBase <- c(pUnder, pOver, pNorm)[or]
    mn$factorCount <- rep(nFactor, n)
    names(mn$factorCount) <- names(genes)
    mn$recruitCount <- rep(0, n)
    names(mn$recruitCount) <- names(genes)
    mn$pNew <- pBase
    mn$expressionLevel <- rep(0, n)
    names(mn$expressionLevel) <- names(genes)
    
    mn$geneStatus <- rep(-1, n)
    names(mn$geneStatus) <- names(genes)
    
    mn$time <- paste0("T", 1:tMax)
    mn$factorDf <- data.frame(Gene = genes %>% str_remove_all("\\d"),T0 = factorCount)
    mn$recruitDf <- data.frame(Gene = genes %>% str_remove_all("\\d"),T0 = recruitCount)
    mn$expressionDf <- data.frame(Gene = genes %>% str_remove_all("\\d"),T0 = expressionLevel)
    mn$geneStatusDf <- data.frame(Gene = genes %>% str_remove_all("\\d"),T0 = geneStatus)
    mn$onThresh <- onThresh

}

reset <- function(nFactor = 100, onThresh = 100, tMax = 100, n = 300)
{#browser()
    mn <- global_env()
    nam <- names(mn$factorCount)
    mn$factorCount <- rep(nFactor, n)
    names(mn$factorCount) <- nam
    mn$recruitCount <- rep(0, n)
    names(mn$recruitCount) <- nam
    mn$pNew <- mn$pBase
    mn$expressionLevel <- rep(0, n)
    names(mn$expressionLevel) <- nam
    
    mn$geneStatus <- rep(-1, n)
    names(mn$geneStatus) <- nam
    genes <- mn$factorDf$Gene
    mn$time <- paste0("T", 1:tMax)
    mn$factorDf <- data.frame(Gene = genes %>% str_remove_all("\\d"),T0 = mn$factorCount)
    mn$recruitDf <- data.frame(Gene = genes %>% str_remove_all("\\d"),T0 = mn$recruitCount)
    mn$expressionDf <- data.frame(Gene = genes %>% str_remove_all("\\d"),T0 = mn$expressionLevel)
    mn$geneStatusDf <- data.frame(Gene = genes %>% str_remove_all("\\d"),T0 = mn$geneStatus)
    mn$onThresh <- onThresh
}
