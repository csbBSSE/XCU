### Author: Kishore Hari
### Description: Script containing functions to update the XCU system



### Updating the activation factors available. Input is a vector of factors the same
## length as the number of genes ----
factorUpdate <- function(factorCount)
{
    updOrder <- sample(1:length(factorCount), length(factorCount))
    fc <- factorCount
    f <- sapply(updOrder, function(x){
        neighbours <- c(x-1, x+1)
        if(x == 1) neighbours <- x+1
        if(x == length(fc)) neighbours <- x - 1
        k <- which(fc[x] < fc[neighbours])
        if (length(k) ==0)
            return()
        k <- neighbours[k]
        for (y in k)
        {#browser()
            p <- floor((fc[y] - fc[x])/4)
            fc[y] <<- fc[y] - p
            fc[x] <<- fc[x] + p
        }
    })
    #browser()
    fc
}
factorUpdate <- cmpfun(factorUpdate)


### Recruitment of factors onto the genes ----
recruit <- function(recruitCount, pBase, pNew,factorCount, thresh = 100)
{
    #browser()
    available <- which(factorCount >= 1)
    probs <- runif(length(available))
    recruitCount[available] <<- recruitCount[available] + (pNew[available] > probs)
    factorCount[available] <<- factorCount[available] - (pNew[available] > probs)
    pNew <<- pBase + 0.2*(recruitCount/(thresh + recruitCount))
    
}
recruit <- cmpfun(recruit)


### Production of RNA molecules for each gene ----
expression <- function(recruitCount, expressionLevel, geneStatus, thresh = 400,
                       deathProb = 0.35, koffRand = F)
{
    deathChoice <- runif(length(expressionLevel))<deathProb
    expressionLevel[geneStatus == 1] <<- expressionLevel[geneStatus == 1] + 1
    expressionLevel[deathChoice] <<- sapply(expressionLevel[deathChoice], 
                                            function(x){max(x-1,0)})
    p <- recruitCount/(thresh + recruitCount)
    probs <- runif(length(p))
    if (koffRand)
    {
        p[geneStatus == 1] <- 0.5
    }
    geneStatus <<- ifelse(p > probs, 1, -1)
    # change <- ifelse(p > probs, -1, 1)
    # geneStatus <<- geneStatus*change
}
expression <- cmpfun(expression)

### One step of the simulation ----
stepFunc <- function()
{
    #browser()
    mn <- global_env()
    recruit(mn$recruitCount, mn$pBase, mn$pNew,mn$factorCount)
    expression(mn$recruitCount, mn$expressionLevel, mn$geneStatus, 
               thresh = mn$onThresh, koffRand = mn$koffRand)
    factorCount <<- factorUpdate(mn$factorCount)
}
stepFunc <- cmpfun(stepFunc)
