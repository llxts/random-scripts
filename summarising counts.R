library(tidyverse)
library(stringr)

#get list of files in folder

filelist <- list.files(path = "./CountData", recursive = TRUE,
                       pattern = "\\.txt$", 
                       full.names = TRUE)

#map_df is from purr. The map_dfr() function will take a vector, apply a function to each element of it, and then return the results as a data frame bound row-by-row
#set_names set names in a vector from the filelist
#in map_df, you apply function read_table
#read_table reads the data per file in filelist, set colnames to false to prevent having the first line being a header, and id creates a new column with the identity of the sourcedata


df <- filelist %>%
  set_names(.) %>%
  map_df(read_table, col_names = FALSE, .id = "FileName")


df %>% separate(X1, c('CRE', 'BC'), sep = "_") -> df

colnames(df)[4] <- 'BC_count' #change colname to bc count

#formating for mpranalyze
df_mpra <- df
df_mpra$FileName <- gsub(pattern = "./CountData/", replacement = "", x = df_mpra$FileName, fixed = TRUE)
df_mpra$FileName <- gsub(pattern = "counts_bowtie_F268.txt", replacement = "", x = df_mpra$FileName, fixed = TRUE)
df_mpra$SamType <- lapply(X = df_mpra$FileName, FUN = function(t) if (grepl("DNA", t)) "DNA" else "RNA")
df_mpra$Organ <- lapply(X = df_mpra$FileName, FUN = function(t) if (grepl("H", t)) "Heart" else "Liver")
df_mpra$SampleNo <- str_extract(df_mpra$FileName, "\\d+")

#pivot table, grouping multiple variables
mpra_h <- df_mpra[-1] %>% group_by(Organ) %>%filter(grepl("Heart", Organ)) %>% pivot_wider(names_from = c(SamType, SampleNo, BC), values_from = BC_count, names_sep = ".", values_fill = 0)
mpra_h <- mpra_h[-2]
mpra_l <- df_mpra[-1] %>% group_by(Organ) %>%filter(grepl("Liver", Organ)) %>% pivot_wider(names_from = c(SamType, SampleNo, BC), values_from = BC_count, names_sep = ".", values_fill = 0)
mpra_l <- mpra_h[-2]


#group data by filename and cre to come up with the number of retrieved barcodes per cre per sample
df_BCperSam <- df %>% group_by(FileName, CRE)%>%summarize(retrievedBC = n(), .groups='drop', sumBC = sum(BC_count, na.rm = TRUE))

#now we are ready to do the calculations. as per MPRAflow, normalising the insert = (no. of bc +1)/(total number of bc retrieved + 1)/total counts * 10**6
#add = true to augment the grouping since group_by overrides the existing grouping variables

df_normCounts <- df_BCperSam %>% group_by(FileName)%>%mutate(totalCount = sum(sumBC)) %>% group_by(CRE, .add = TRUE) %>% summarise(normInsertCount = (sumBC + 1)/(retrievedBC+1)/(totalCount) * 10**6)
df_normCounts$SampleName <-  gsub(pattern = "./CountData/", replacement = "", x = df_normCounts$FileName, fixed = TRUE)
df_normCounts$SampleName <-  gsub(pattern = "counts_bowtie_F268.txt", replacement = "", x = df_normCounts$SampleName, fixed = TRUE)
df_normCounts$SamType <- lapply(X = df_normCounts$SampleName, FUN = function(t) if (grepl("DNA", t)) "DNA" else "RNA")
df_normCounts$Organ <- lapply(X = df_normCounts$SampleName, FUN = function(t) if (grepl("H", t)) "Heart" else "Liver")
df_normCounts$SampleNo <- str_extract(df_normCounts$SampleName, "\\d+")


#filter to separate heart and liver into diff files
df_Heart <- df_normCounts %>% group_by(Organ) %>%filter(grepl("Heart", Organ))
df_Liver <- df_normCounts %>% group_by(Organ) %>%filter(grepl("Liver", Organ))

#small working function to join both RNA and DNA, and remove the working files later
df_heart_DNA <- df_Heart %>% group_by(SamType) %>% filter(grepl("DNA", SamType)) %>% select(CRE, normInsertCount, SamType, SampleNo)
df_heart_RNA <- df_Heart %>% group_by(SamType) %>% filter(grepl("RNA", SamType)) %>% select(CRE, normInsertCount, SamType, SampleNo)
df_heart_comb <- left_join(df_heart_DNA, df_heart_RNA, by = c("CRE", "SampleNo"))
df_liver_DNA <- df_Liver %>% group_by(SamType) %>% filter(grepl("DNA", SamType)) %>% select(CRE, normInsertCount, SamType, SampleNo)
df_liver_RNA <- df_Liver %>% group_by(SamType) %>% filter(grepl("RNA", SamType)) %>% select(CRE, normInsertCount, SamType, SampleNo)
df_liver_comb <- left_join(df_liver_DNA, df_liver_RNA, by = c("CRE", "SampleNo"))

#rename Columns
colnames(df_heart_comb)[2] <-  "DNA normCount"
colnames(df_heart_comb)[5] <-  "RNA normCount"

colnames(df_liver_comb)[2] <-  "DNA normCount"
colnames(df_liver_comb)[5] <-  "RNA normCount"

#remove samtype
df_heart_comb<- df_heart_comb[,-3]
df_heart_comb<- df_heart_comb[,-5]

df_liver_comb<- df_liver_comb[,-3]
df_liver_comb<- df_liver_comb[,-5]

#calculate log2FC
df_heart_Final <- df_heart_comb  %>% mutate(log2 = log2(`RNA normCount`/`DNA normCount`))
df_liver_Final <- df_liver_comb  %>% mutate(log2 = log2(`RNA normCount`/`DNA normCount`))

#extract to excel for further analysis
write.csv(df_heart_Final, "df_heart_Final.csv")
write.csv(df_liver_Final, "df_liver_Final.csv")


#correlation plots
library(corrplot)

#format df for corrplot 
corrHeart <- df_heart_Final %>% pivot_longer(cols = c("DNA normCount", "RNA normCount", "log2"), names_to = "type", values_to= "normCount")
corrHeart$working <- paste(corrHeart$type, "_", corrHeart$SampleNo)
corrHeart <- corrHeart[-2:-3] %>% pivot_wider(names_from = working, values_from = normCount, values_fill = 0) %>% select(sort(names(.)))

heart_DNA_normcounts <- as.data.frame(corrHeart[2:12], row.names = corrHeart$CRE)
heart_DNA_normcounts <- log2(heart_DNA_normcounts[,-1])
heart_DNA_normcounts[heart_DNA_normcounts == "-Inf"]<- 0
corr_hd <- cor(heart_DNA_normcounts)

heart_RNA_normcounts <- read.delim("heart_RNA_normcounts.txt")
rownames(heart_RNA_normcounts) <- heart_RNA_normcounts$Row.Labels
heart_RNA_normcounts <- log2(heart_RNA_normcounts[,-1])
heart_RNA_normcounts[heart_RNA_normcounts == "-Inf"]<- 0
corr_hr <- cor(heart_RNA_normcounts)

liver_DNA_normcounts <- read.delim("liver_DNA_normcounts.txt")
rownames(liver_DNA_normcounts) <- liver_DNA_normcounts$Row.Labels
liver_DNA_normcounts <- log2(liver_DNA_normcounts[,-1])
liver_DNA_normcounts[liver_DNA_normcounts == "-Inf"]<- 0
corr_ld <- cor(liver_DNA_normcounts)

liver_RNA_normcounts <- read.delim("liver_RNA_normcounts.txt")
rownames(liver_RNA_normcounts) <- liver_RNA_normcounts$Row.Labels
liver_RNA_normcounts <- log2(liver_RNA_normcounts[,-1])
liver_RNA_normcounts[liver_RNA_normcounts == "-Inf"]<- 0
corr_lr <- cor(liver_RNA_normcounts)

heart_log2 <- read.delim("log2-heart.txt")
rownames(heart_log2) <- heart_log2$Row.Labels
heart_log2 <- heart_log2[,-1]
corr_logH <- cor(heart_log2)

liver_log2 <- read.delim("log2-liver.txt")
rownames(liver_log2) <- liver_log2$Row.Labels
liver_log2 <- liver_log2[,-1]
corr_logL <- cor(liver_log2)

corrplot.mixed(corr_logL, order = 'hclust')

# general function test 
#import list of BCs 
dftest <- read.table("bclist.txt", quote="\"", comment.char="")
for (i in filelist) {
  dftest <- read.table(i) %>% full_join(dftest, i, by = c('V1'))
}

#rename columns
colnamelist <- as.list(filelist)
colnamelist <- lapply(X = colnamelist, FUN = function(t) gsub(pattern = "./CountData/", replacement = "", x = t, fixed = TRUE))
colnamelist <- lapply(X = colnamelist, FUN = function(t) gsub(pattern = "counts_bowtie_F268.txt", replacement = "", x = t, fixed = TRUE))
colnames(dftest)[2:ncol(dftest)] <- colnamelist
colnames(dftest)[1] <- "CRE"


