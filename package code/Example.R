# setwd("~/ESPCA-CODE")print
library(CpGassoc)
library(psych)
library(calculatep.R)
gene_new <- read.csv("gene_new.txt", header=FALSE)
Network <- read.csv("network.txt", header=FALSE)
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
  edges[[i]] = c((Network[i,1]+1),(Network[i,2]+1))
}
out3 =  ESPCA(gene, k = 3, edges, k.group=n,we=0.3, t = 0.1)
myString <- "Sparse PCA calculate end!"
print(myString)
write.csv (out3[["V"]], file ="ESPCA.csv")
out <- data.frame(matrix(unlist(out3[3]), nrow=3, byrow=T))
myString <- "Uncorrected Result"
print(myString)
test1<-cpg.assoc(gene_new,input_phenotype$Class,large.data=FALSE)
plot(test1)
Cell_type_simulation = Calculate_p(gene_new,input_phenotype,out,n = 400,nn =4)
print(Cell_type_simulation)
myString <- "end!"
print(myString)


