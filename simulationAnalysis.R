source("plot_theme.R")

statTest <- function(x1, x2)
{
    df <- data.frame(Gene = geneList, y = x)
    
}

XaXaFiles <- list.files(".", "XaXa.*\\d.csv")
XaXiFiles <- list.files(".", "XaXi.*\\d.csv")

csvFiles <- c(XaXaFiles, XaXiFiles)

formatter <- sapply(csvFiles, function(x){
    df <- read.csv(x) %>% 
        arrange(Gene) %>% mutate(ID = 1:nrow(.)) %>%
        mutate(Gene = paste0(Gene, ID)) %>% select(-ID)
    write.csv(df, paste0("dir/",x),row.names = F)
})

geneIDs <- c(sample(1:200, 100), 201:300)
lowGenes <- geneIDs[1:100]
highGenes <- geneIDs[101:200]

## pair-wise statTest

factorCounts <- c(100, 200, 500, 700, 1000)
setwd("dir")
tTests <- sapply(factorCounts, function(x){
    dfXa <- read_csv(XaXaFiles[str_detect(XaXaFiles, paste0("_",x, ".csv"))], 
                     col_types = cols())
    dfXi <- read_csv(XaXiFiles[str_detect(XaXaFiles, paste0("_",x, ".csv"))], 
                     col_types = cols())
    cols <- colnames(dfXa)[str_detect(colnames(dfXa), "E")]
    d <- sapply(cols, function(y){
        df <- cbind(dfXi[[y]], dfXa[[y]])
        tTestUp <- t.test(df[highGenes, 1], df[highGenes, 2], alternative = "greater")
        tTestDo <- t.test(df[lowGenes, 1], df[lowGenes, 2], alternative = "two.sided")
        c(tTestUp$p.value, tTestUp$estimate, tTestDo$p.value, tTestDo$estimate)
    }) %>% t %>% data.frame %>% 
        set_names(c("upRegP", "upRegXi", "upRegXa","downRegP", "downRegXi", "downRegXa"))
    write.csv(d, paste0("ExpressionDiff_", x,".csv"), row.names = F)
    d <- sapply(cols, function(y){
        Xi <- dfXi[[y]]
        Xa <- dfXa[[y]]
        tTestDo <- t.test(Xa[highGenes], Xa[lowGenes], alternative = "greater")
        tTestUp <- t.test(Xi[highGenes], Xi[lowGenes], alternative = "greater")
        c(tTestUp$p.value, tTestUp$estimate, tTestDo$p.value, tTestDo$estimate)
    }) %>% t %>% data.frame %>% 
        set_names(c("XiP", "XiUp", "XiDown","XaP", "XaUp", "XaDown"))
    write.csv(d, paste0("ExpressionWithinDiff_", x,".csv"), row.names = F)
    
    
    
    dfXa <- read_csv(XaXaFiles[str_detect(XaXaFiles, paste0("burst",x, ".csv"))], 
                     col_types = cols())
    dfXi <- read_csv(XaXiFiles[str_detect(XaXaFiles, paste0("burst",x, ".csv"))], 
                     col_types = cols())
    cols <- paste0("burst", 1:100)
    d <- cols %>% 
        sapply(function(y){
            Xi <- dfXi[[y]]
            Xa <- dfXa[[y]]
            tTestDo <- t.test(Xa[highGenes], Xa[lowGenes], alternative = "greater")
            tTestUp <- t.test(Xi[highGenes], Xi[lowGenes], alternative = "greater")
            c(tTestUp$p.value, tTestUp$estimate, tTestDo$p.value, tTestDo$estimate)
        }) %>% t %>% data.frame %>% 
        set_names(c("XiP", "XiUp", "XiDown","XaP", "XaUp", "XaDown"))
    
    write.csv(d, paste0("burst_", x,".csv"), row.names = F)
})
expressionTest <- list.files(".", "ExpressionDiff")
burstTest <- list.files(".", "burst_")
expressionWithinTest <- list.files(".", "Within")
tTestPlots <- sapply(expressionTest, function(x){
    nFact <- str_extract(x, "\\d+")
    df <- read.csv(x)
    upFrac <- sum(df$upRegP > 0.05)/nrow(df)
    downFrac <- sum(df$downRegP > 0.05)/nrow(df)
    c(as.integer(nFact), upFrac, downFrac)
}) %>% t %>% data.frame %>% set_names(c("nFactor", "upRegulation", "downRegulation")) %>%
    gather(key = "Regulation", value = "Fraction", -nFactor)
ggplot(tTestPlots %>% 
           mutate(Regulation = ifelse(Regulation == "downRegulation", "Down", "Up")), 
       aes(x = nFactor, y = Fraction, color = Regulation)) + 
    geom_point(size = 2) + 
    geom_smooth(method = "lm")+
    labs(x = "# of factors", y = "Probability of\ninsignificant difference") +
    theme_Publication()
ggsave("ExpressionLevels.png", width = 5.5, height = 5.5)

tTestPlots <- sapply(expressionWithinTest, function(x){
    nFact <- str_extract(x, "\\d+")
    df <- read.csv(x)
    upFrac <- sum(df$XiP > 0.05)/nrow(df)
    downFrac <- sum(df$XaP > 0.05)/nrow(df)
    c(as.integer(nFact), upFrac, downFrac)
}) %>% t %>% data.frame %>% set_names(c("nFactor", "Xi", "Xa")) %>%
    gather(key = "Regulation", value = "Fraction", -nFactor)
ggplot(tTestPlots, aes(x = nFactor, y = Fraction, color = Regulation)) + 
    geom_point(size = 2) + 
    geom_smooth(method = "lm")+
    labs(x = "# of factors", y = "Probability of\ninsignificant difference") +
    theme_Publication()
ggsave("exprWithinDiff.png", width = 5, height = 5.5)

setwd("..")


### burstFreq diff ----

dfXa <- read.csv("XaXa_burst100.csv") %>% select(Gene, burst2) %>% 
    set_names(c("Gene", "burstFreq")) %>% slice(c(lowGenes, highGenes)) %>%
    mutate(Gene = str_remove_all(Gene, "\\d")) %>%
    mutate(Gene = ifelse(Gene == "Up", "Up", "Down"))
ggbetweenstats(dfXa, x = Gene, y = burstFreq ,type = "parametric",
              pairwise.comparisons = T, ylab = "Burst Frequency", title = "",
              point.path = F, pairwise.display = "all", results.subtitle = T,
              p.adjust.method = "none") +
    theme_Publication() + labs(y="Burst Frequency")
ggsave("XaXa_burstFreq.png", height = 6.5, width = 5.5)


dfXa <- read.csv("XaXi_burst100.csv") %>% select(Gene, burst2) %>% 
    set_names(c("Gene", "burstFreq")) %>% slice(c(lowGenes, highGenes)) %>%
    mutate(Gene = str_remove_all(Gene, "\\d")) %>%
    mutate(Gene = ifelse(Gene == "Up", "Up", "Down"))
ggbetweenstats(dfXa, x = Gene, y = burstFreq ,type = "parametric",
               pairwise.comparisons = T, ylab = "Burst Frequency", title = "",
               point.path = F, pairwise.display = "all", results.subtitle = T,
               p.adjust.method = "none") +
    theme_Publication() + labs(y="Burst Frequency")
ggsave("XaXi_burstFreq.png", height = 6.5, width = 5.5)


### expression Diff  -----
dfXa <- read.csv("XaXa_100.csv") %>% select(Gene, E2) %>% 
    set_names(c("Gene", "burstFreq")) %>% slice(c(lowGenes, highGenes)) %>%
    mutate(Gene = str_remove_all(Gene, "\\d")) %>%
    mutate(Gene = ifelse(Gene == "Up", "Up", "Down"))
ggbetweenstats(dfXa, x = Gene, y = burstFreq ,type = "parametric",
               pairwise.comparisons = T, ylab = "Burst Frequency", title = "",
               point.path = F, pairwise.display = "all", results.subtitle = T,
               p.adjust.method = "none") +
    theme_Publication() + labs(y="Expression")
ggsave("XaXa_expression.png", height = 6.5, width = 5.5)


dfXi <- read.csv("XaXi_100.csv") %>% select(Gene, E1) %>% 
    set_names(c("Gene", "burstFreq")) %>% slice(c(lowGenes, highGenes)) %>%
    mutate(Gene = str_remove_all(Gene, "\\d")) %>%
    mutate(Gene = ifelse(Gene == "Up", "Up", "Down"))
ggbetweenstats(dfXi, x = Gene, y = burstFreq ,type = "parametric",
               pairwise.comparisons = T, ylab = "Burst Frequency", title = "",
               point.path = F, pairwise.display = "all", results.subtitle = T,
               p.adjust.method = "none") +
    theme_Publication() + labs(y="Expression")
ggsave("XaXi_expression.png", height = 6.5, width = 5.5)

df <- dfXa %>% mutate(XaXi = dfXi$burstFreq) %>%
    set_names(c("Gene","XaXa", "XaXi")) %>%
    gather(key = "CellType", value = "Expression", -Gene)
ggbetweenstats(df %>% filter(Gene != "Up"), x = CellType, y = Expression ,type = "parametric",
               pairwise.comparisons = T, ylab = "Burst Frequency", title = "",
               point.path = F, pairwise.display = "all", results.subtitle = T,
               p.adjust.method = "none") +
    theme_Publication() + labs(y="Expression")
ggsave("XaXi_expression_Down.png", height = 6.5, width = 5.5)
