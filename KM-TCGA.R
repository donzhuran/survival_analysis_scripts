
setwd("C:/Users/Korisnik/Desktop/Survival analysis")

##### Install dependencies
install.packages("https://github.com/compgenome365/TCGA-Assembler-2.git")
library(githubinstall)
#githubinstall("TCGA-Assembler-2")
library(assemblerr)

##### Download gene expression data from TCGA using TCGA-Assembler2
geneExp2 <- DownloadRNASeqData(cancerType = "BRCA", assayPlatform
                               = "gene.normalized_RNAseq",tissueType = "TP",
                               saveFolderName = ".")

##### Download clinical data


##### Merging gene expression and clinical data
setwd("C:/Users/Korisnik/Desktop/Survival analysis")
getwd()
TCGAGEdata <- read.csv(file="BRCA__gene.normalized_RNAseq__TP__20230211101034.txt",sep = "\t",row.names = 1, check.names = F)

# Import data
TCGAMergedataOriginal <- read.csv(file ="TCGA-2.csv", sep = "\t",check.names = F,row.names = 1) # check.names = F to not replace hyphen with a dot
dim(TCGAMergedataOriginal)

TCGAMergedata <- read.csv(file ="TCGA-3.csv", sep = "\t",check.names = F,row.names = 1)
dim(TCGAMergedata)
View(TCGAMergedata)
getwd()

table(TCGAMergedata$vital_status)
table(TCGAMergedata$`Overall survival`)
table(TCGAMergedata$race)

#surv_object <- Surv(time = TCGAMergedata$`Overall survival`, event = TCGAMergedata$vital_status)
#fit <- survfit(surv_object~TCGAMergedata$gender)



######Estimation of survival difference between white and black or African american

# exclude three race cathegories
dataRace=subset(TCGAMergedata, (TCGAMergedata$race!='not reported')&(TCGAMergedata$race!='american indian or alaska native')&(TCGAMergedata$race!='asian')) 
table(dataRace$race)

# Calculate survival analysis
surv_object <- Surv(time = dataRace$`Overall survival`, event = dataRace$vital_status)
fit <- survfit(surv_object~dataRace$race) 
ggsurvplot(fit, data=dataRace,pval=TRUE)
surv_diff1 <- survdiff(Surv(dataRace$`Overall survival`, dataRace$vital_status) ~ dataRace$race, data = dataRace)
surv_diff1 # p-value
ggsurvplot(fit,data= dataRace,
           pval=TRUE, conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

ggsurvplot(fit,data=dataRace,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

----------------------
#### Estimation of survival difference between early and late stages of breast cancer

table(TCGAMergedata$ajcc_pathologic_stage)

# merge stage data into "early" or "late stage"
dataStage=mutate(TCGAMergedata, ajcc_pathologic_stage=ifelse(((ajcc_pathologic_stage=="Stage I")|(ajcc_pathologic_stage=="Stage IA")|(ajcc_pathologic_stage=="Stage IB")|(ajcc_pathologic_stage=="Stage II")|(ajcc_pathologic_stage=="Stage IIA")|(ajcc_pathologic_stage=="Stage IIA")|(ajcc_pathologic_stage=="Stage IIB")),"Early-Stage","Late-Stage"))
table(dataStage$ajcc_pathologic_stage)

# Calculate survival analysis
surv_object1 <- Surv(time = dataStage$`Overall survival`, event = dataStage$vital_status)
fit <- survfit(surv_object1~dataStage$ajcc_pathologic_stage) 
ggsurvplot(fit, data=dataStage,pval=TRUE)
surv_diff2 <- survdiff(Surv(dataStage$`Overall survival`, dataStage$vital_status) ~ dataStage$ajcc_pathologic_stage, data = dataStage)
surv_diff2
ggsurvplot(fit,data= dataStage,
           pval=TRUE, conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

ggsurvplot(fit,data=dataStage,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


------------------------------
#### Effect of gene expression in survival #BRAF

# Create survival object
surv_object2 <- Surv(time = TCGAMergedata$`Overall survival`, event = TCGAMergedata$vital_status)

# check is there NA values in BRAF column
which(is.na(TCGAMergedata$`BRAF|673`)) 

# Define categories based on median values for BRAF expressions (cutoff is "median")
fit3 <- survfit(surv_object2~TCGAMergedata$`BRAF|673`>median(TCGAMergedata$`BRAF|673`))

# Visualize
ggsurvplot(fit3, data=TCGAMergedata,pval=TRUE)
#surv_diff3 <- survdiff(Surv(TCGAMergedata$`Overall survival`, TCGAMergedata$vital_status) ~ TCGAMergedata$`BRAF|673`, data = TCGAMergedata)
#surv_diff3


#### Effect of gene expression on survival #BRCA1
which(is.na(TCGAMergedata$`TP53|7157`))
fit5 <- survfit(surv_object2~TCGAMergedata$`BRCA1|672`>median(TCGAMergedata$`BRCA1|672`))
ggsurvplot(fit5, data=TCGAMergedata,pval=TRUE) 
