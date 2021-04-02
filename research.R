
# Inspecting Counts
do1$counts[1:10, 1:5]
colnames(do1$counts)

# Retrieving Groups 

coldata$condition4[match(colnames(do1$counts), coldata$SampleName)]
df1 <- do1$counts
str(df1)

x <- rownames(df1)
df1 <- as.data.frame(df1)
head(df1)
df1$geneid <- rownames(df1)

df1long <- df1 %>% pivot_longer(1:14, names_to= "ids")

df1long
df1long$group <- coldata$condition4[match(df1long$ids, coldata$SampleName)]
df1long$id <- substr(df1long$ids, 1, 2)
head(df1long)
df1long$gene <- genes$SYMBOL[match(df1long$geneid, genes$ENSEMBL)]
head(df1long)  

range(df1long$value)
df1long[df1long$value>25000,]


# Retrieving Groups 

coldata$condition4[match(colnames(do2$counts), coldata$SampleName)]
df2 <- do2$counts
str(df2)

x <- rownames(df2)
df2 <- as.data.frame(df2)
head(df2)
df2$geneid <- rownames(df2)

df2long <- df2 %>% pivot_longer(1:14, names_to= "ids")

df2long
df2long$group <- coldata$condition4[match(df2long$ids, coldata$SampleName)]
df2long$id <- substr(df2long$ids, 1, 2)
head(df2long)
df2long$gene <- genes$SYMBOL[match(df2long$geneid, genes$ENSEMBL)]
head(df2long)  

alldf <- rbind(df1long, df2long)
alldf
alldf %>% filter(gene == "Peg3") %>% 
  ggplot(aes(x= group, y=value)) + geom_boxplot() + geom_point() 
alldf <- unique(alldf)
