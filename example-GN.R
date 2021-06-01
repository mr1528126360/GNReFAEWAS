# setwd("~/ESPCA-CODE")print
library(CpGassoc)
library(psych)
source('GN-ReFAEWAS.R')
gene_new <- read.csv("gene_new.txt", header=FALSE)
result.1_p2 <- read.csv("network.txt", header=FALSE)
input_phenotype <- read.csv("input_phenotype.txt", header=TRUE)
nn = 4
n = 400
myString <- "start!"
print(myString)
edges = list()
gene1 = data.matrix(gene_new)
gene = t(gene1)
re_av = array(0,dim=c(200))
re_av1 = array(0,dim=c(200))
for (i in 1:27955){
  edges[[i]] = c((result.1_p2[i,1]+1),(result.1_p2[i,2]+1))
}
out3 =  ESPCA(gene, k = 3, edges, k.group=n,we=0.3, t = 0.1)
myString <- "end!"
print(myString)
#write.csv (out3[["V"]], file ="ESPCA.csv")
df <- data.frame(matrix(unlist(out3[3]), nrow=3, byrow=T))
df1 = data.matrix(df)
df2 = t(df1)
S1 = abs(df2[,1])
S2 = order(S1,decreasing=TRUE)[1:n]
data1 <- gene_new[c(S2),1:ncol(gene_new)]
A = principal(data1)
PC = A[["residual"]]

S3 = abs(df2[,2])
S4 = order(S3,decreasing=TRUE)[1:n]
data2 <- gene_new[c(S4),1:ncol(gene_new)]
A2 = principal(data2)
PC2 = A2[["residual"]]

S5 = abs(df2[,3])
S6 = order(S5,decreasing=TRUE)[1:n]
data3 <- gene_new[c(S6),1:ncol(gene_new)]
A3 = principal(data3)
PC3 = A3[["residual"]]

print(PC[1:6,])
