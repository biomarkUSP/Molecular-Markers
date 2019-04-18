
#'========================================================================================
#' Métodos utilizados para análise de matrizes de SNPs 
#'========================================================================================
#'
#' Autores:
#'        Ana Letycia B. Garcia
#'        Germano Martins F Costa Neto
#'        Nathalia Salgado Silva
#'        Rafael Massahiro Yassue
#'----------------------------------------------------------------------------------------

#' R packages

# pacotes necessarios
pkg <- c("snpReady","rrBLUP","BGGE","rrBLUP","BGLR","FactoMineR","adegenet","superheat")

# verifica instalacao
if(length(setdiff(pkg, rownames(installed.packages()))) > 0) {   
  install.packages(setdiff(pkg, rownames(installed.packages())))}   

# chama pacotes
sapply(pkg,library,character.only = TRUE) 

#' Data set
#'----------------------------------------------------------------------------------------

nomeoarquivo = "maize" #' usando arquivo maize

maize <- read.csv(paste("https://raw.githubusercontent.com/biomarkUSP/Datasets/master/",
                      nomeoarquivo,".csv",sep=""),header=T)

str(maize)


maize[1:10,1:10]
maize$IND = paste("ind",maize$IND,sep="")
maize.m <- data.frame(maize[,-2],row.names = 1)
head(maize.m)


superheat(t(maize.m))
superheat(cov(maize.m),row.dendrogram = T,col.dendrogram = T)

superheat(cov(t(maize.m)))

pca1 <- PCA(maize.m)
group <- HCPC(pca1)


# separando dados fenotipicos
y = data.frame(maize[,1:2],row.names = 1)
head(y)


## Segundo dataset
#'----------------------------------------------------------------------------------------

nomeoarquivo = "geno" #' usando arquivo geno

geno = read.csv(paste("https://raw.githubusercontent.com/biomarkUSP/Datasets/master/",
                      nomeoarquivo,".csv",sep=""),header=T)

str(geno)
dim(geno)

geno[1:10,1:10]
summary(geno)

SNP = data.frame(geno,row.names = 1)

# arquivo de mapa
nomeoarquivo = "map"
map = read.csv(paste("https://raw.githubusercontent.com/biomarkUSP/Datasets/master/",
                     nomeoarquivo,".csv",sep=""),header=T)
head(map)
tail(map)
rownames(map) <- map$marker
head(map)

str(map)
map$chr = as.factor(map$chr)

all(colnames(SNP) == rownames(map))

# usando snpReady
args(raw.data)


# Control de qualidade (QC)
SNP = as.matrix(SNP)

QC = raw.data(data = SNP, frame = "wide",
               hapmap = map, maf = 0.05, call.rate = 0.90, 
               base = FALSE, imput = TRUE, imput.type = "knni", plot = TRUE)

# Reporta o que foi feito
QC$report

# Matriz de marcadores "limpa"
M = QC$M.clean

QC$Hapmap


# Manualmente!

# imputacao
W  = matrix(0,dim(SNP)[1],dim(SNP)[2])


for (j in 1:ncol(SNP)) {
  # cat('j = ', j, '\n')
  W[, j] = ifelse(is.na(SNP[, j]), mean(SNP[, j], na.rm = TRUE), SNP[, j])
}

summary(W)

superheat(SNP)

# Control de qualidade

freq = colMeans(W)/2
maf = ifelse(freq > 0.5, 1 - freq, freq)

maf.index = which(maf < 0.05)
SNP2 <- W[, -maf.index]

dim(W)
dim(SNP)
dim(M)
dim(SNP2)

# Predicao Genomica

head(y)


# aplicando ao maize data set

SNP = as.matrix(maize.m)

QC = raw.data(data = SNP, frame = "wide",
              maf = 0.05, call.rate = 0.90, 
              base = FALSE, imput = TRUE, imput.type = "knni", plot = TRUE)

SNP=QC$M.clean
dim(maize.m)
dim(SNP)


# ridge regression
fit = mixed.solve(y = y$GY,Z = SNP)
summary(fit)

fit$Vu
fit$beta
fit$u
barplot(fit$u)

ETA.rr = list(SNP = list(X = SNP,model = "BRR"))
ETA.ba = list(SNP = list(X = SNP,model = "BayesA"))


fit2 = BGLR(y = y$GY,ETA = ETA.rr,verbose = F)

fit2$ETA$SNP

fit2$ETA$SNP$b

barplot(fit2$ETA$SNP$b)


fit3 = BGLR(y = y$GY,ETA = ETA.ba,verbose = F)

fit3$ETA$SNP$b
barplot(fit3$ETA$SNP$b)

cor(y$GY,fit3$yHat)


# Genomic relationship matrix

G = G.matrix(SNP, method = "VanRaden", format = "wide", plot = F)
G$Ga
G$Gd

superheat(G$Gd)


# Na unha!

SNP_s = scale(SNP, center = TRUE, scale = TRUE)

n = nrow(SNP_s) # numero de individuos!
m = ncol(SNP_s) # numero de marcas!


G = tcrossprod(SNP_s)/m
G = G + diag(n) * 0.001

superheat(G)

# GBLUP
fit4 = mixed.solve(y = y$GY,K = G)
fit4$Vu
fit4$u
barplot(fit4$u)

cor(y$GY,fit4$u)
