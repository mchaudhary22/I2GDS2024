<This is the README for the R code>

## data import

```r

#feature counts file location
hybridfeaturects = "~/featurecounts_hybrid.csv"

#sample information file
Metaseq_data_day1today10 <- "~/Metaseq_data_day1today10.csv"

#read in counts file
rawhybridcts <- read.csv(hybridfeaturects, header = TRUE)

#split file by arabidopsis and phelipanche (AT/PE) and convert rownames to gene ids
rawctsAt <- rawhybridcts %>% 
  filter(str_detect(rawhybridcts$Geneid, "^AT")) %>%
  column_to_rownames(var = "Geneid")

#phelipanche split
rawctsPe <- rawhybridcts%>% 
  filter(str_detect(rawhybridcts$Geneid, "^T99")) %>% 
  column_to_rownames(var = "Geneid") 

#remove rows with 0 counts
At <- rawctsAt[rowSums(rawctsAt)>=1,]
Pe <- rawctsPe[rowSums(rawctsPe)>=1,]

#convert to matrices
countDataAt <- as.matrix(At)
countDataPe <- as.matrix(Pe)

```

```{r}
###Load Sampledata

#load sample info file
coldata <- read.csv(Metaseq_data_day1today10, row.names=1)

view(coldata)

#factor by treatment types
coldata <- coldata[,c("Day","Treatment")]

#factor each treatment
coldata$Treatment <- factor(coldata$Treatment)
coldata$Day <- factor(coldata$Day)

```


```{r}
###Data Subsetting
##From meeting with David (9/5/24)

#subset the counts so it is only counts of phelipanche
Pesubset <- countDataPe %>% 
  as_tibble() %>% 
  dplyr::select(starts_with("D")& ends_with(c("P1","P2","P3","P4")))

#subset column data (sample info) so it only contains phelipanche
colpesubset <- coldata[coldata$Treatment == "Phelipanche",]

#repeat with arabidopsis
Atsubset <- countDataAt %>%
  as_tibble() %>%
  dplyr::select(starts_with("D") & ends_with(c("C1","C2","C3","C4","C5")))

colatsubset <- coldata[coldata$Treatment == "Control",]

```

```{r}
###Run DESeq2 on At reads
ddsAt<- DESeqDataSetFromMatrix(countData = Atsubset, colData = colatsubset, design = ~Day)
ddsAt<-DESeq(ddsAt)
res01<-results(ddsAt)
summary(res01)
reOrdered<-res01[order(res01$padj),]
rldAt<- rlog(ddsAt, blind=TRUE)
vsdAt<-varianceStabilizingTransformation(ddsAt, blind=TRUE)

meanSdPlot(assay(vsdAt))
```

```{r}
###Run DESeq2 on Pe reads
ddsPe<- DESeqDataSetFromMatrix(countData = Pesubset, colData = colpesubset, design = ~Day)
ddsPe<-DESeq(ddsPe)
res02<-results(ddsPe)
summary(res02)
reOrdered<-res02[order(res01$padj),]
rldAt<- rlog(ddsPe, blind=TRUE)
vsdAt<-varianceStabilizingTransformation(ddsPe, blind=TRUE)

meanSdPlot(assay(vsdPe))
```

```{r}
###Log Fold Change Shrinkage Visualization/Ranking
#From David's code

resultsNames(ddsPe)

#Using LFC shrinkage to compare genes by L2FC
lfcD2P <- lfcShrink(ddsPe, coef="Day_Day2_vs_Day1", type="apeglm")
lfcD3P <- lfcShrink(ddsPe, coef="Day_Day3_vs_Day1", type="apeglm")
lfcD5P <- lfcShrink(ddsPe, coef="Day_Day5_vs_Day1", type="apeglm")
lfcD7P <- lfcShrink(ddsPe, coef="Day_Day7_vs_Day1", type="apeglm")
lfcD10P <- lfcShrink(ddsPe, coef="Day_Day10_vs_Day1", type="apeglm")

#be sure to cite:  "Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences.    Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

##Convert LFC tables to tibbles

lfcD2P_tbl <- as_tibble(lfcD2P, rownames = "Gene_ID") %>% 
  mutate(group = "lfcD2P") 
lfcD3P_tbl <- as_tibble(lfcD3P, rownames = "Gene_ID") %>% 
  mutate(group = "lfcD3P") 
lfcD5P_tbl <- as_tibble(lfcD5P, rownames = "Gene_ID") %>% 
  mutate(group = "lfcD5P") 
lfcD7P_tbl <- as_tibble(lfcD7P, rownames = "Gene_ID") %>% 
  mutate(group = "lfcD7P") 
lfcD10P_tbl <- as_tibble(lfcD10P, rownames = "Gene_ID") %>% 
  mutate(group = "lfcD10P")

```



```{r}
###Joining Lfc/other values
##From David's code

all_lfc <- lfcD2P_tbl %>% 
  inner_join(lfcD3P_tbl, by = "Gene_ID") %>% 
  inner_join(lfcD5P_tbl, by = "Gene_ID") %>% 
  inner_join(lfcD7P_tbl, by = "Gene_ID") %>% 
  inner_join(lfcD10P_tbl, by = "Gene_ID") %>% 
  dplyr::select(-pvalue, -pvalue.x, -pvalue.y) %>% 
  mutate(
    lfc_D2P = log2FoldChange.x,
    lfc_D3P = log2FoldChange.y,
    lfc_D5P = log2FoldChange.x.x,
    lfc_D7P = log2FoldChange.y.y,
    lfc_D10P = log2FoldChange,
    padj_D2P = padj.x,
    padj_D3P = padj.y,
    padj_D5P = padj.x.x,
    padj_D7P = padj.y.y,
    padj_D10P = padj,
  ) %>% 
  view()

lfc_long <- all_lfc %>% 
  dplyr::select(1, 29:33) %>% 
  pivot_longer(
    cols = -Gene_ID,
    names_to = c("LF2C", "Group"),
    names_sep = "_",
    values_to = c("Expression")
)
pval_long <- all_lfc %>% 
  dplyr::select(1, 34:38) %>% 
  pivot_longer(
    cols = -Gene_ID,
    names_to = c("adjp", "Group"),
    names_sep = c("_"),
    values_to = c("Pvalue")
)

exp_tbl <- inner_join(lfc_long, pval_long) %>%
  dplyr::select(-LF2C, -adjp) %>% 
  view()

View(all_lfc)

```



```{r}
###Annotating where LF2c >= 2
##From David's code

exp_tbl <- exp_tbl %>% 
  mutate(
    ExpSig = case_when(Expression >= log(7.4) & Pvalue <= 0.05 ~ "Up-regulated",
                           Expression <= -log(7.4) & Pvalue <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
    ) %>%
  view()

exp_tbl$Group <- factor(exp_tbl$Group, levels = c("D2P", "D3P", "D5P", "D7P", "D10P")) #refactoring to fix ordering for graphs

```




```{r}
###MA plots with top gene labeling
##From David's Code

p3 <- ggplot(exp_tbl, aes(Expression, -log(Pvalue,10))) +
  geom_point(aes(color = ExpSig), size = 2/5) +
  facet_wrap(vars(Group)) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"Pvalue")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "darkred")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  theme_bw()
p3

top <- 10
top_genes <- bind_rows(
  exp_tbl %>% 
    filter(ExpSig == 'Up-regulated') %>% 
    arrange(Pvalue, desc(abs(Expression))) %>% 
    head(top),
  exp_tbl %>% 
    filter(ExpSig == 'Down-regulated') %>% 
    arrange(Pvalue, desc(abs(Expression))) %>% 
    head(top)
)


p3 <-  p3 +
  geom_label(data = top_genes,
                   mapping = aes(Expression, -log(Pvalue,10), label = Gene_ID),
                   size = 2)
p3

```

```
