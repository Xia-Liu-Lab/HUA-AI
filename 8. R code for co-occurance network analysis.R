# 20230101
# Ludi Liu
### --------------------------
### -------- sparCC ----------
### ---------------------------
# linux FastSpar
# correlation coefficient
fastspar --threshold 0.3 --otu_table count.tsv --correlation correlation.tsv --covariance covariance.tsv
# P value
mkdir bootstrap_countsfastspar_bootstrap --otu_table count.tsv --number 1000 --prefix bootstrap_counts/count
mkdir bootstrap_correlationparallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*
fastspar_pvalues --otu_table count.tsv --correlation correlation.tsv --prefix bootstrap_correlation/cor_count_ --permutations 1000 --outfile pvalues.tsv

### --------------------------
### -------- degree ----------
### ---------------------------
m <- read.csv()  #correlation coefficient matrix
#????һ??????ͼ
g <- make_undirected_graph(t(m[,1:2]))
plot(g)

# ------- 1.total degree---------
a <- data.frame(degree(g,mode="in"))
b <- data.frame(degree(g,mode="out"))
c <- data.frame(degree(g,mode="total"))
n <- data.frame(degree(g,normalized = T))

degree <- cbind(a,b,c,n)
colnames(degree) <- c("degree(in)","degree(out)","degree(total)","degree(normalized)")

# ------ 2.closeness degree------
cl <- data.frame(closeness(g,vids=V(g),mode = "total"))

colnames(cl) <- "closeness"
degree <- cbind(degree,cl)

write.csv(degree,"degree.csv" )