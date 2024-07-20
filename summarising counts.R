library(tidyverse)
library(stringr)

#get list of files in folder

filelist <- list.files(path = "./Counts", recursive = TRUE,
                       pattern = "\\.txt$", 
                       full.names = TRUE)

#map_df is from purr. The map_dfr() function will take a vector, apply a function to each element of it, and then return the results as a data frame bound row-by-row
#set_names set names in a vector from the filelist
#in map_df, you apply function read_table
#read_table reads the data per file in filelist, set colnames to false to prevent having the first line being a header, and id creates a new column with the identity of the sourcedata


df <- filelist %>%
  set_names(.) %>%
  map_df(read_table, col_names = FALSE, .id = "FileName")


df %>% separate(X1, c('CRE', 'BC'), sep = "_") %>% drop_na() -> df


colnames(df)[4] <- 'BC_count' #change colname to bc count

#formating for mpranalyze
df_mpra <- df
df_mpra$FileName <- gsub(pattern = "./Counts/", replacement = "", x = df_mpra$FileName, fixed = TRUE)
df_mpra$FileName <- gsub(pattern = "counts.txt", replacement = "", x = df_mpra$FileName, fixed = TRUE)
df_mpra$SamType <- lapply(X = df_mpra$FileName, FUN = function(t) if (grepl("DNA", t)) "DNA" else "RNA")
df_mpra$Organ <- lapply(X = df_mpra$FileName, FUN = function(t) if (grepl("H", t)) "Heart" else "Liver")
df_mpra$SampleNo <- str_extract(df_mpra$FileName, "\\d+")

#pivot table, grouping multiple variables
mpra_dna <- df_mpra[-1] %>% group_by(SamType) %>%filter(grepl("DNA", SamType)) %>% pivot_wider(names_from = c(SamType, Organ, SampleNo, BC), values_from = BC_count, names_sep = ".", values_fill = 0)
mpra_dna <- mpra_dna[-2]
mpra_rna <- df_mpra[-1] %>% group_by(SamType) %>%filter(grepl("RNA", SamType)) %>% pivot_wider(names_from = c(SamType, Organ, SampleNo, BC), values_from = BC_count, names_sep = ".", values_fill = 0)
mpra_rna <- mpra_rna[-2]

write.csv(mpra_dna, "mpra-dna.csv", col.names = TRUE)
write.csv(mpra_rna, "mpra-rna.csv", col.names = TRUE)

#group data by filename and cre to come up with the number of retrieved barcodes per cre per sample
df %>% group_by(FileName, CRE)%>%summarize(retrievedBC = n(), .groups='drop', sumBC = sum(BC_count, na.rm = TRUE)) -> df_BCperSam

#now we are ready to do the calculations. as per MPRAflow, normalising the insert = (no. of bc +1)/(total number of bc retrieved + 1)/total counts * 10**6
#add = true to augment the grouping since group_by overrides the existing grouping variables

df_normCounts <- df_BCperSam %>% group_by(FileName)%>%mutate(totalCount = sum(sumBC)) %>% group_by(CRE, .add = TRUE) %>% summarise(normInsertCount = (sumBC + 1)/(retrievedBC+1)/(totalCount) * 10**6)
df_normCounts$SampleName <-  gsub(pattern = "./Counts/", replacement = "", x = df_normCounts$FileName, fixed = TRUE)
df_normCounts$SampleName <-  gsub(pattern = "counts.txt", replacement = "", x = df_normCounts$SampleName, fixed = TRUE)
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


#correlation plots
library(corrplot)

#format df for corrplot 
corrHeart <- df_heart_Final %>% pivot_longer(cols = c("DNA normCount", "RNA normCount", "log2"), names_to = "type", values_to= "normCount")
corrHeart$working <- paste(corrHeart$type, "_", corrHeart$SampleNo)
corrHeart <- corrHeart[-2:-3] %>% pivot_wider(names_from = working, values_from = normCount, values_fill = 0) %>% select(sort(names(.)))

corrLiver <- df_liver_Final %>% pivot_longer(cols = c("DNA normCount", "RNA normCount", "log2"), names_to = "type", values_to= "normCount")
corrLiver$working <- paste(corrLiver$type, "_", corrLiver$SampleNo)
corrLiver <- corrLiver[-2:-3] %>% pivot_wider(names_from = working, values_from = normCount, values_fill = 0) %>% select(sort(names(.)))

#Save files
write.csv(corrHeart, "df_heart_Final.csv")
write.csv(corrLiver, "df_liver_Final.csv")

heart_DNA_normcounts <- as.data.frame(corrHeart[2:12])
rownames(heart_DNA_normcounts) <- corrHeart$CRE
heart_DNA_normcounts <- log2(heart_DNA_normcounts)
heart_DNA_normcounts[heart_DNA_normcounts == "-Inf"]<- 0
corr_hd <- cor(heart_DNA_normcounts)

heart_log2_normcounts <- as.data.frame(corrHeart[13:23])
rownames(heart_log2_normcounts) <-  corrHeart$CRE
heart_log2_normcounts[is.na(heart_log2_normcounts)]<- 0
corr_hlog2 <- cor(heart_log2_normcounts)

heart_RNA_normcounts <- as.data.frame(corrHeart[24:34])
rownames(heart_RNA_normcounts) <- corrHeart$CRE
heart_RNA_normcounts <- log2(heart_RNA_normcounts)
heart_RNA_normcounts[heart_RNA_normcounts == "-Inf"]<- 0
heart_RNA_normcounts[is.na(heart_RNA_normcounts)]<- 0
corr_hr <- cor(heart_RNA_normcounts)

liver_DNA_normcounts <- as.data.frame(corrLiver[2:12])
rownames(liver_DNA_normcounts) <- corrLiver$CRE
liver_DNA_normcounts <- log2(liver_DNA_normcounts)
liver_DNA_normcounts[liver_DNA_normcounts == "-Inf"]<- 0
corr_ld <- cor(liver_DNA_normcounts)

liver_log2_normcounts <- as.data.frame(corrLiver[13:23])
rownames(liver_log2_normcounts) <- corrLiver$CRE
liver_log2_normcounts[is.na(liver_log2_normcounts)]<- 0
corr_llog2 <- cor(liver_log2_normcounts)

liver_RNA_normcounts <- as.data.frame(corrLiver[24:34])
rownames(liver_RNA_normcounts) <- corrLiver$CRE
liver_RNA_normcounts <- log2(liver_RNA_normcounts)
liver_RNA_normcounts[liver_RNA_normcounts == "-Inf"]<- 0
liver_RNA_normcounts[is.na(liver_RNA_normcounts)]<- 0
corr_lr <- cor(liver_RNA_normcounts)

corrplot.mixed(corr_hlog2, order = 'hclust') 




##### MPRAnalyze ######

library(MPRAnalyze)
ce.dnaCounts <- as.matrix(read.delim("dnacountsrm8.txt", row.names = 1, stringsAsFactors = TRUE))
ce.rnaCounts <- as.matrix(read.delim("rnacounts_indnarm8.txt", row.names = 1, stringsAsFactors = TRUE))
ce.dnaAnnot <- as.data.frame(read.delim("dna-annot-rm8.txt", row.names = 1, stringsAsFactors = TRUE))
ce.rnaAnnot <- as.data.frame(read.delim("rna-annot-rm8.txt", row.names = 1, stringsAsFactors = TRUE))
ce.dnaAnnot$replicate <- factor(ce.dnaAnnot$Rep)
ce.dnaAnnot$barcode <- factor(ce.dnaAnnot$BC)
# ce.dnaAnnot$seqbatch <- factor(ce.dnaAnnot$seqbatch)

ce.rnaAnnot$replicate <- factor(ce.rnaAnnot$Rep)
ce.rnaAnnot$barcode <- factor(ce.rnaAnnot$BC)
# ce.rnaAnnot$seqbatch <- factor(ce.rnaAnnot$seqbatch)


#reverse the levels
ce.dnaAnnot$Organ<-factor(ce.dnaAnnot$Organ, levels=rev(levels(ce.dnaAnnot$Organ)))
ce.rnaAnnot$Organ<-factor(ce.rnaAnnot$Organ, levels=rev(levels(ce.rnaAnnot$Organ)))



obj <- MpraObject(dnaCounts = ce.dnaCounts, rnaCounts = ce.rnaCounts, rnaAnnot = ce.rnaAnnot, dnaAnnot = ce.dnaAnnot)

## If the library factors are different for the DNA and RNA data, separate 
## estimation of these factors is needed. We can also change the estimation 
## method (Upper quartile by default)
obj <- estimateDepthFactors(obj, lib.factor = c("Rep"),
                            which.lib = "dna", 
                            depth.estimator = "uq")
obj <- estimateDepthFactors(obj, lib.factor = c("Organ", "Rep"),
                            which.lib = "rna", 
                            depth.estimator = "uq")

#Empirical testing, just to see which sequence has a regulatory function
obj1 <- analyzeQuantification(obj = obj,
                              dnaDesign = ~ Rep + Organ,
                              rnaDesign = ~ Organ)

##extract alpha values from the fitted model --> gives you transcription rate
alpha <- getAlpha(obj1, by.factor = "Organ")


##visualize the estimates
boxplot(alpha)

alpha2 <- alpha
alpha2[is.na(alpha2)] <- 0

## test
res.heart <- testEmpirical(obj = obj1,
                           statistic = alpha2$Heart)
res.liver <- testEmpirical(obj = obj1,
                           statistic = alpha2$Liver)
alpha2$Liver_pval <- res.liver$pval.mad
alpha2$Heart_pval <- res.heart$pval.mad

alpha2$Liver_pvaladj <- p.adjust(alpha2$Liver_pval,method = "fdr")
alpha2$Heart_pvaladj <- p.adjust(alpha2$Heart_pval,method = "fdr")

par(mfrow=c(2,2))

hist(res.heart$pval.mad, main="heart, all")
hist(res.liver$pval.mad, main="liver, all")



#comparing between Organ. the first factor level is the reference -- i.e. heart is the reference. so the more negative logFC, the more heart specific
obj2 <- analyzeComparative(obj = obj, 
                           dnaDesign = ~ BC + Rep + Organ, 
                           rnaDesign = ~ Rep + Organ , 
                           reducedDesign = ~ Rep)

# comparing using std errors
objse <- analyzeComparative(obj = obj, 
                            dnaDesign = ~ BC + Rep + Organ, 
                            rnaDesign = ~ Organ, 
                            fit.se=TRUE)

# this tests for the significance of the Organ (i.e. liver or heart) in the model, with Liver expression being used as the contrast and heart as ref i.e. Liver vs Heart
results.se <- testCoefficient(objse, "Organ", "Heart") # positive logFC indicates that the CRE is more expressed in the heart compared to the liver (logFC is UPREGULATION in Liver vs HEART)

#testing for differential activity
res <- testLrt(obj2)
summary (res)
res2<- res
res2[is.na(res2)] <- 0
res2$zscore <- (res2$logFC-mean(res2$logFC))/sd(res2$logFC)



## plot log Fold-Change
reslogfc = plot(density(res2$logFC))

## plot volcano (-log10 pval --> higher is better)
resvolcano = plot(res2$logFC, -log10(res2$pval))


##export result
write.csv(res2, "mpranalyzeresults_noseqbatch_indna_less10-rm8.csv", row.names=TRUE)
write.csv(alpha2, "mpranalyze_quant-alpha_indna_less10-rm8.csv", row.names=TRUE)


