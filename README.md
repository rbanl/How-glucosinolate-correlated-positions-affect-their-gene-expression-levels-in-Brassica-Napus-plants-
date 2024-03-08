# How-glucosinolate-correlated-positions-affect-their-gene-expression-levels-in-Brassica-Napus-plants-

Full details of how my analysis was performed:

To test my hypothesis of an uneven distribution of FDR significantly correlated genes, I obtained all the genes that had a FDR value of less than 0.05 from the “final ” dataset I produced. I used the subset 4 times, 1 for each section I allocated to chromosome 1.This was done for the purpose of the analysis. Each subset a set number of positions. The positions represent the genes’ rank orders Each subset also contained code that specified which chromosome the data should be extracted from. Of the 27 significantly correlated genes found, 17 of them were located in the 1st  subset .6,3 and 1  more were found. My hypothesis is supported by both of these findings and most of the results appear to be consistent with research. Most active genes are found in euchromatin (Murakami, 2013) which is at the distal parts of chromosomes It is also not uncommon to find phenotypes like cancer being associated with particular regions on chromosomes. (Thomassen, Tan and Kruse, 2008). The 22 genes found in the first section would correspond to the beginning of the periphery of a chromosome.This skew in  significantly correlated gene distribution might be caused by that region of DNA having hypermutable sections , which could be full of minisatellites at a large number of gene positions. Microsatellites tend to be G-C rich  (Poulsen, Kahl and Weising, 1993) and therefore may be prone to a higher mutation rate. In the Poulsen study, it was shown that higher levels of some lead to intraspecies polymorphism. In B.Napus  , it is possible that higher levels a particular minisatellite increased the chances of proliferation of the 17 genes in a small  section.

In order to measure whether different parts of the chromosome were more had more gene expression , unpaired Wilcoxon tests were performed between all the territories in only chromosome 1 since there was only 1 significantly correlated gene on the 2nd chromosome. All the results were insignificant and the largest difference was found between territories 1 and 2 in chromosome 1, where the P-values were measured at  0.1171  and 0.5476. No test was performed between the 3rd and 4th territories as there was only one correlated gene in the 4th section. Average P-values  were also compared between the  genes in chromosomes 1 and 2 section

The non-significant p-value for the means between sections 3 and 4 may be  unreliable as because the data samples were very small . After finding most of the significant genes in one small part of the territories I would’ve expected large and significant differences in gene expression between the territories. There were 11 fewer genes in territory 2. In territory 2 there were a few genes with very high average expression levels. No very highly expressed genes were in territory 1.  In the 2nd territory the genes with the highest  expression levels had an average position of 1200. It may be considered that the genes in section 2 were on average much more influential than the genes in section 1 since there were no significant differences in expression even though there were less genes. It also could mean that the most influential genes are found further away from the periphery than expected. This also could be reasoned because even in territory 1 the 2 highest ranking genes are not too far away the average position of 1200 It also potentially means that section 1 had a high number of non-causally related genes to glucosinolate content and that section 2’s genes on average were more causally-related to glucosionolate content.One might also consider the second point because gene expression was so high.
Although the data doesn’t present the actual positions of the genes , it presents us with the rank order and the higher in rank a gene is the more likely it is to be away from the periphery of the chromosome 1.









PLEASE SEE BELOW FOR THE CODE:


Supplementary methods
###################################################################################################
#TThe distribution of  glucosinolate correlated genes  and differing expression levels in different
#chromosome territories in Brassica Napus  from the the OSR101 data set                                            #
###################################################################################################


####PREPARATION OF DATA FOR ANALYSIS


# Full data set of OSR101 reads per kio base per million aligned reads(RPKM) was loaded for analysis 
OSR101_RPKM <- read.delim("http://www-users.york.ac.uk/~ah1309/BigData/data/OSR101_RPKM.txt", header=TRUE)



#checked the class and dimensions of data
#of my data
class(OSR101_RPKM)
dim(OSR101_RPKM)



#represents the first 5 rows and first 6 columns of data
OSR101_RPKM[1:5,1:6]



#sets the gene names as the rows and  the next step involved checking data is alright
rownames(OSR101_RPKM) <- OSR101_RPKM$Gene
OSR101_RPKM[1:2,]



#removed the first column of data and created a new object called OSR101_RPKM2 that is checked in the next step
OSR101_RPKM2 <- OSR101_RPKM[,-1]
OSR101_RPKM2[1:5,1:5]


#checked if all data is in a numeric format
is.numeric(OSR101_RPKM2)


#converted data(currently in matrix form) to a numeric  form
OSR101_num <- as.matrix(sapply(OSR101_RPKM2, as.numeric))



#Checking data is numeric  and set row names of bngenenum2's to bngenenum 
is.numeric(OSR101_num)
row.names(OSR101_num) <- row.names(OSR101_RPKM2)



#checked the gene expressions levels
hist(OSR101_num[1,])
hist(OSR101_num)



#used to calculate mean expressin  values for each gene
rowmeans <- rowMeans(OSR101_num)


#subset expression is used to keep only the rows
keepmeans <- subset(OSR101_num, rowmeans >= 1)



#produces a summary of the data
summary(keepmeans)[,1:5]



#checking dimensions of keepmeans
dim(keepmeans)



#double checked that keepmeans object is properly produced.It is also saved here
write.table(keepmeans, "OSR101_num.txt", quote=F, sep="\t")



#saving the summary data
write.table(summary(keepmeans), "Summary.txt", quote=F, sep="\t")




#Read in saved file,then removed an unwanted column for the time being.In the last step I check the data
OSR101_fin<- read.table("OSR101_num.txt",header = T)
OSR101_final<- OSR101_fin[-1]
OSR101_final[1:5,1:5]



#Loaded the file containing  data for the glucosinolate content and lines.
gluc <- read.delim("http://www-users.york.ac.uk/~ah1309/BigData/data/Glucosinolates.txt", header=T)



#setting  the rownames from the Line column, then I load the data set
rownames(gluc) <- gluc$Line
gluc



#Transposing data after setting max.print so rows aren't omitted .In the last line 
# I check the data
options(max.print = 100000)
tOSR <- t(OSR101_final)
tOSR[1:5,1:6]


#set max print again , then merged data transposed data
#with glucosinolate data, then it's printed
options(max.print = 10000)
OSR_merge <- merge(gluc, tOSR, by="row.names")
OSR_merge[1:5,1:5]


#set the row names and removed columns
rownames(OSR_merge) <- OSR_merge$Line
OSR_merge <- OSR_merge[,-c(1:2)]
OSR_merge[1:5,1:5]


#use of a linear model to find any trend between trait and
#the RPKMs
lm0 <- lm(OSR_merge$Trait~OSR_merge[,2])



anova <- as.data.frame(anova(lm0)[1,])
#Intercept of regression line
intercept = as.data.frame(coefficients(lm0))[1,1]
#Gradient of regression line
gradient = as.data.frame(coefficients(lm0))[2,1]
#R2 (correlation coefficient)
R2 <- as.data.frame(summary(lm0)$r.squared)



# the different results are being joined to the same object
#using c function
result1 <- as.data.frame(c(anova, R2, intercept, gradient))
result1 



#counting number of columns in the merged data set and printing it 
numcol <- ncol(OSR_merge)
numcol



results <- as.data.frame(matrix(nrow = 0, ncol = 7996))


#ran for loop for 
for (i in 2:numcol){
  lm0 <- lm(OSR_merge$Trait~OSR_merge[,i])
  anova <- as.data.frame(anova(lm0)[1,])
  intercept = as.data.frame(coefficients(lm0))[1,1]
  gradient = as.data.frame(coefficients(lm0))[2,1]
  R2 <- as.data.frame(summary(lm0)$r.squared)
  result1 <- as.data.frame(c(anova, R2, intercept, gradient))
  
  colnames(results) <- colnames(result1)
  results <-rbind(results, result1)
}


results[1:5,]



# putting the row names back on and assigning column names to columns , then in th e last steps I'm checked
#results and then saved it 
rownames(results) <- colnames((OSR_merge[,2:numcol]))
colnames(results) <- c("Df", "Sum.Sq", "Mean.Sq", "F.value", "P.value", "R2", "Intercept", "Gradient")
results[1:10,]
write.table(results, "results.txt", quote=F, sep="\t")



#loading results that was saved in results.txt
results <- read.table("results.txt", header = TRUE)



#rows of results are  ordered by p-value significance and put in results_sorted. 
#The first 5 rows  of those are then printed
results_sorted <- results[order(results$P.value),]
results_sorted[1:5,]



#data containing the genes' positions on the chromosomes is loaded
OSR_chr <- read.delim("http://www-users.york.ac.uk/~ah1309/BigData/data/osr_dir.txt", row.names = 1)
OSR_chr <- read.delim("http://www-users.york.ac.uk/~ah1309/BigData/data/osr_dir.txt")





#OSR_chr is merged with results and row names are set again
results_OSR_chrom <- merge(OSR_chr, results, by="row.names")
results_OSR_chrom[1:5,]




#Created  2nd results_sorted which contains the p values ordered by their significance 
#in the merged data set from the previous step.In the 2nd step a column with the logged p-values is added
#and Data is checked in the last step.
results_sorted2 <- results_OSR_chrom[order(results_OSR_chrom$P.value),]
results_sorted2$logP <- -log10(results_sorted2$P.value)
results_sorted2[1:5,]



# HOW MANY GENE EXPRESSION MARKERS PASS THE MULTIPLE CORRECTIONS THRESHOLD?
###############################################################################
library(tidyverse)


#subset for results on chromosome 1 and 2 was created using Graph ==..
chr1 <- subset(results_sorted2, results_sorted2$Graph==1)
chr2<- subset(results_sorted2, results_sorted2$Graph==2)


#finding number of rows for chromosome 1 
nrow(chr1)


#making a bonferroni adjustment for the p-values to rule out type 1 errors
#, the bon adjustment is also logged to match the logP values printed 
bon <- -log10(0.05/5374)


#ggplot is loaded and then used to create  2 manhattan plots for chromosomes 1 and 2
#showing positions of genes according to their logged p-values,
#plots contain bonferroni threshold

library(ggplot2)


ggplot(chr1, aes(x= Position), y=logP))+ geom_point(size=0.6)+theme_classic() + geom_hline(yintercept = bon)

ggplot(chr2, aes(x=Position ), y=logP))+ geom_point(size=0.6)+theme_classic() + geom_hline(yintercept = bon)



#plot is created for gene expressions of section 1 values and position
plot(chrom1_sec1$Position, gene.data4_SEC1$Mean_gene_exprsn, ylab="Gene expression", xlab=paste("Position , Section 1 of chromosome1"
abline(lm1)


plot(chrom1_sec1$Position, gene.data4_SEC1$Mean_gene_exprsn, ylab="Gene expression", xlab=paste("Position , Section 1 of chromosome1", colnames(chrom1_sec1$Position)))

v=nrow(results_sorted2)
#Add new column with numbers from 1 to v
results_sorted2$Rank <- 1:v


#Correcting each p-value and converting them to FDR values then  addinga column of them to results_sorted2
results_sorted2$FDR  <-results_sorted2$P.value*(v/results_sorted2$F.value)



#changing the first column name to Gene in results_sorted2
colnames(results_sorted2)[1] <- c("Gene")
results_sorted2[1:5,]



#saving results_sorted2
write.table(results_sorted2, "Results_sorted 2 with FDR.txt", quote=F, sep="\t", row.names=F)



#WHAT IS THE  DISTRIBUTION OF THE  SIGNIFICANTLY ASSOCIATED  DISTRIBUTION OF GENES IN THE CHROMOSOMES LIKE?
#############################################################################################################


#storing saved results in final, then it's dimensions are double-checked
final<- read.delim("Results_sorted 2 with FDR.txt", header=T)
dim(final)

mean(chrom1_sec2$FDR) 
mean(finalchrom1_sec1$FDR)
mean(chrom1_sec3$FDR)
chrom1_sec4$FDR
getwd()
res11 <- wilcox.test(chsec ~ Psn, paired = TRUE)
C1 <- read.csv("C1.csv" , header = TRUE)
res11 <-wilcox.test(chsec, Psn, paired = TRUE)
wilcox.test(before, after, paired = TRUE)


#WHAT IS THE  DISTRIBUTION OF THE  SIGNIFICANTLY ASSOCIATED GENES IN THE CHROMOSOMES LIKE  ? 
####################################################################################################




#used subset function to extract genes from final data  which have an FDR of less than 0.05,then it's dimensions
#are checked
final2 <- subset(final, final$FDR<0.05)
dim(final2)




#created a subset of final2 named finalchrom1.It contains genes that are  only on  chromosome 1 
#and that have  significant FDR values.In the next step I found out how many genes there are.
finalchrom1 <- subset(final2, Graph == 1  )
nrow(finalchrom1)




#another subset is created to find out how many genes are within the range of position 1 to 2000 in 
# chromosome 1
chrom1_sec1 <- subset(finalchrom1,Position<1100)
nrow(chrom1_sec1)



#Some of the analysis performed on chrom1_sec1 was repeated for this and the next 2 sections 
#Section 2 of chromosome 1 
chrom1_sec2 <- subset(finalchrom1 ,Position>1100 & Position<3100)
nrow(chrom1_sec2)



#Section 3 of chromosome 1 
chrom1_sec3 <- subset(finalchrom1,Position>3099 &Position <6000)
nrow(chrom1_sec3)





#Section 4 of chromosome 1
chrom1_sec4 <- subset(finalchrom1,Position>5999)
nrow(chrom1_sec4)



#Chromosome 2 - Analysis was also repeated for chromosome 2 
finalchrom2 <- subset(final2, Graph == 2  )
nrow(finalchrom2)


#only 1 gene was present so no tests were ran 
finalchrom2$Position



#loaded gene sequences 
seqs <- read.delim("http://www-users.york.ac.uk/~ah1309/BigData/data/genes.txt")
options(stringsAsFactors = FALSE)


# checked integrity of gene sequences
seqs[1:5,]

#seet rownames of seqs data set
rownames(seqs) <- seqs$Gene



# extracted  section 1's genes and  inserted it into sec1genes as a vector so its genes could be printed
sec1genes <- as.vector(final$Gene[1:17])
sec1genes


#printed genes topseqs to display genes
topseqs <-  seqs[sec1genes,]
topseqs


# ARE THERE DIFFERENT EXPRESSION LEVELS IN DIFFERENT CHROMOSOMAL TERRITORIES?
####################################################################################################

#  extracted out gene expression values from main data set - OSR101_RPKM by row
gene.data2 <-  OSR101_RPKM[c(6271,6188,6315,2392,652,7496,352,6655,1341,5679,6819,5371,792,2288,3368,3613,5771),]


# calculated mean gene expression for each row  and added it  in the form of the 
#row_sums column at the end of the data frame
gene.data2$row_sum = apply(gene.data2[,-1], 1, mean)


#ran this code so no rows were emitted
options(max.print = 100000)


#checked data
gene.data2


#deleted rows that didnt contain mean gene expressions
gene.data3 <- gene.data2[ -c(1,2:102) ]

# changed row_sum column to  mean gene expression to  make results in  data frame clearer
gene.data4_SEC1 <- gene.data3 %>% 
  rename(
    Mean_gene_exprsn = row_sum
    
  )


#checked data
gene.data4_SEC1





#repeated process previously carried out for the 2nd chromosome territory

gene.data2_SEC2 <-  OSR101_RPKM[c(404,6818,4367,674,4891,3922),]


gene.data2_SEC2$row_sum = apply(gene.data2_SEC2[,-1], 1, mean)


gene.data2_SEC2


gene.data3_SEC2 <- gene.data2_SEC2[ -c(1,2:102)]



gene.data4_SEC2 <- gene.data3_SEC2 %>% 
  rename(
    Mean_gene_exprsn = row_sum
    
  )



gene.data4_SEC2



#loaded tidyverse package in order to perform wilcoxn test 
library(tidyverse)


#cheked normality of 1st and 2nd chromosome territories' gene expression so I could use the appropriate tests
shapiro.test(gene.data4_SEC1$Mean_gene_exprsn)
shapiro.test(gene.data4_SEC2$Mean_gene_exprsn)


#carried out wilcoxon test since gene.data4_SEC2 wasnt normally distributed 
wilcox.test(gene.data4_SEC1$Mean_gene_exprsn , gene.data4_SEC2$Mean_gene_exprsn, alternative = "two.sided")


#repeated process for chromosome territory 3 
gene.data2_SEC3 <-  OSR101_RPKM[c(1494,4627,2243),]



gene.data2_SEC3$row_sum = apply(gene.data2_SEC3[,-1], 1, sum)

options(max.print = 100000)

gene.data2_SEC3

gene.data3_SEC3 <- gene.data2_SEC3[ -c(1,2:102) ]



gene.data4_SEC3 <- gene.data3 %>% 
  rename(
   Mean_gene_exprsn = row_sum
    
  )

gene.data4_SEC3



#cheked normality of 2nd and 3rd chromosome territories' gene expression so I could use the appropriate tests
shapiro.test(gene.data4_SEC2$Mean_gene_exprsn)
shapiro.test(gene.data4_SEC3$Mean_gene_exprsn)


#carried out wilcoxon test since both are not normally distributed 
wilcox.test(gene.data4_SEC2$Mean_gene_exprsn , gene.data4_SEC3$Mean_gene_exprsn, alternative = "two.sided")

# wilcoxon test whcih compared the average FDR  of the genes in sections 1 and 2 
wilcox.test(chrom1_sec1$FDR , chrom1_sec2$FDR, alternative = "two.sided")



