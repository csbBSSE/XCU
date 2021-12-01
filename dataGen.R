### Generation of count matrix from simulations -----

source("init.R")
source("func.R")
source("plot_theme.R")

koffRand <- T
tMax <- 1000
onThresh <- 70
nFactor <- 1000
high <- c(0.3, 0.8)
med <- c(0.2, 0.5)
low <- c(0.0, 0.5)
nhalf <- 100

nFactorList <- c(100, 200, 500, 700, 1000)


initFunc(pUnderRange = med, pnormRange = med, pOverRange = med,
         n = 3*nhalf, nhalf = nhalf, nFactor = nFactor, onThresh = onThresh,
         tMax = tMax)

for (nFactor in nFactorList)
{
    
    expressioni <- expressionDf %>% select(Gene)
    bursti <- geneStatusDf 
    sapply(1:100, function(i){
        reset(nFactor, onThresh, tMax, n = 3*nhalf)
        dummy <- sapply(time, function(x){
            factorCount <<- stepFunc()
            factorDf[[x]] <<- factorCount
            recruitDf[[x]] <<- recruitCount
            expressionDf[[x]] <<- expressionLevel
            geneStatusDf[[x]] <<- geneStatus
        })
        expressioni[[paste0("E", i)]] <<- expressionDf %>% 
            select(Gene,T1000) %>% 
            # mutate(Gene = Gene %>% str_remove_all("\\d")) %>%
            select(T1000) %>% unlist
        bursti[[paste0("burst", i)]] <<- geneStatusDf %>% 
            # filter(str_detect(Gene, "U")) %>%
            unite("Status", -Gene, sep = "") %>%
            mutate(burstFreq = str_count(Status, "-11")/tMax) %>%
            select(Gene, burstFreq) %>% select(burstFreq) %>% unlist
    })
    
    expressioni <- expressioni %>% mutate(pBase = pBase) %>% 
        select(Gene, pBase, contains("E"))
    bursti <- bursti %>% mutate(pBase = pBase) %>% 
        select(Gene, pBase, contains("burst"))
    write.csv(expressioni, paste0("XaXa_", nFactor,".csv"), row.names = F)
    write.csv(bursti, paste0("XaXa_burst", nFactor, ".csv"), row.names = F)
}

initFunc(pUnderRange = low, pnormRange = med, pOverRange = high,
         n = 3*nhalf, nhalf = nhalf, nFactor = nFactor, onThresh = onThresh,
         tMax = tMax)

for (nFactor in nFactorList)
{
    
    expressioni <- expressionDf %>% select(Gene) 
    bursti <- geneStatusDf 
    sapply(1:100, function(i){
        reset(nFactor, onThresh, tMax, n = 3*nhalf)
        dummy <- sapply(time, function(x){
            factorCount <<- stepFunc()
            factorDf[[x]] <<- factorCount
            recruitDf[[x]] <<- recruitCount
            expressionDf[[x]] <<- expressionLevel
            geneStatusDf[[x]] <<- geneStatus
        })
        expressioni[[paste0("E", i)]] <<- expressionDf %>% 
            select(Gene,T1000) %>% 
            # mutate(Gene = Gene %>% str_remove_all("\\d")) %>%
            select(T1000) %>% unlist
        bursti[[paste0("burst", i)]] <<- geneStatusDf %>% 
            # filter(str_detect(Gene, "U")) %>%
            unite("Status", -Gene, sep = "") %>%
            mutate(burstFreq = str_count(Status, "-11")/tMax) %>%
            select(Gene, burstFreq) %>% select(burstFreq) %>% unlist
    })
    
    expressioni <- expressioni %>% mutate(pBase = pBase) %>% 
        select(Gene, pBase, contains("E"))
    bursti <- bursti %>% mutate(pBase = pBase) %>% 
        select(Gene, pBase, contains("burst"))
    write.csv(expressioni, paste0("XaXi_", nFactor,".csv"), row.names = F)
    write.csv(bursti, paste0("XaXi_burst", nFactor, ".csv"), row.names = F)
}

