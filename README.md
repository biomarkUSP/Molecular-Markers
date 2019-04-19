# Marcadores moleculares (tradicionais e SNPs)


O que são marcadores moleculares? Quais os fundamentos das técnicas? Como e por que usá-los?


Após o desenvolvimento teórico sobre o que são marcadores genéticos, quais os protocolos laboratoriais para sua obtenção e processamento, iremos demonstrar algumas aplicações do uso de marcadores genéticos, com ênfase em melhoramento de plantas.

Os exemplos descritos nessa seção envolvem alguns métodos utilizados para a análise de matrizes de SNPs. Primeiramente falaremos sobre como analisar de forma simples uma estrutura de população com base em PCA. Em seguida, como fazer o controle de qualidade de dados moleculares. Por fim algumas aplicações simples com predição genômica para fins de orientar a seleção (ou descarte) de genótipos durante o processo de melhoramento genético de plantas.

Nesta nossa página irão encontrar esse e outros scripts, além dos datasets usados. Também disponibilizamos nessa seção os materiais que utilizamos para orientar nossa apresentação.

Essas aplicações tem como base modelos de genética biométrica, isto é, a área da genética que trabalha com modelagem matemática para descrever processos vinculados à herança de caracteres, a relação genótipo-fenótipo, a interação de genótipos por ambientes, entre outros.

No passado, os geneticistas e biometricistas buscavam compreender fenômenos genéticos com base apenas no fenótipo observado. Com o advento dos marcadores genéticos, e suas propriedades biométricas úteis, tais como alta herdabilidade (próxima a 1), foi possível utiliza-los para melhorar a modelagem da relação genótipo-fenótipo.

Além disso, os marcadores também podem usados para caracterização da estrutura de genética de populações, alavancando estudos de diversidade genética (e molecular), filogenia e evolução (em espécies cultivadas e populações naturais).


## Link para o video (youtube)

<https://www.youtube.com/watch?v=6NxcMr6hC40>

## Scripts (em R)

```R
# Métodos utilizados para análise de matrizes de SNPs

```R
# pacotes necessarios
pkg <- c("snpReady","rrBLUP","BGGE","rrBLUP","BGLR","FactoMineR","adegenet","superheat")

# verifica instalacao
if(length(setdiff(pkg, rownames(installed.packages()))) > 0) {   
  install.packages(setdiff(pkg, rownames(installed.packages())))}   

# chama pacotes
sapply(pkg,library,character.only = TRUE)
```

# Data sets

## Data set I: Maize (171 genótipos, 14 marcas)

```{r}

nomeoarquivo = "maize" #' usando arquivo maize

maize <- read.csv(paste("https://raw.githubusercontent.com/biomarkUSP/Datasets/master/",
                      nomeoarquivo,".csv",sep=""),header=T)

str(maize)


maize[1:10,1:10]
maize$IND = paste("ind",maize$IND,sep="")
maize.m = data.frame(maize[,-2],row.names = 1)
head(maize.m)
```

## Algumas análises exploratórias

```{r}
superheat(t(maize.m))
superheat(cov(maize.m),row.dendrogram = T,col.dendrogram = T)

superheat(cov(t(maize.m)))

pca1 = PCA(maize.m)
group = HCPC(pca1,graph = F)
```

- separando dados fenotipicos para futuras análises

```{r}
y = data.frame(maize[,1:2],row.names = 1)
head(y)

```


## Data set II: Maize (55 genótipos, 218 marcas)

```{r}

nomeoarquivo = "geno" #' usando arquivo geno

geno = read.csv(paste("https://raw.githubusercontent.com/biomarkUSP/Datasets/master/",
                      nomeoarquivo,".csv",sep=""),header=T)

str(geno)
dim(geno)

geno[1:10,1:10]
summary(geno)

SNP = data.frame(geno,row.names = 1)

```


### mapa genético

```{r}

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
```


## Controle de Qualidade

### Uso de snpReady

```{r}
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


```

### QC Manualmente!

```{r}

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
```


# Predicao Genômica

## Diferenças entre rrBLUP, BayesA e GBLUP
```{r}

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

```

_____________________________________________________________________________________________________


## Leituras sugeridas


### Conceitos básicos

Atul Grover & P. C. Sharma (2016) **Development and use of molecular markers: past and present** Critical Reviews in Biotechnology, 36:2, 290-302, <DOI:10.3109/07388551.2014.959891>

Oraguizie. et al. (2007) **Association Mapping in Plants**. Springer, 289p.

Turchetto-Zolet et al. (2017) **Marcadores Moleculares na Era Genômica: Metodologias e Aplicações**. E-book disponível em <https://www.sbg.org.br/sites/default/files/e_book_marcadores_moleculares_sbg_2017_final.pdf>

### Aplicações em melhoramento genético de plantas


Glowinski & Flint-Garcia (2017) **Germplasm Resources for Mapping Quantitative Traits in Maize** IN: The Maize Genome, Springer, Pages 143-159

Singh & Singh (2015) **Marker-Assisted Plant Breeding: Principles and Practices**, Springer,514p.
 
Song & Cregan (2017) **Classical and Molecular Mapping** IN: The Soybean Genome, Springer, Pages 41-56
 
Vieira,  M.L. C. et al.(2016) **Microsatellite markers: what they mean and why they are so useful**. Genetics and Molecular Biology, 39, 3, 312-328

### Em português

Ramalho M.A.P et al. 2012. **Genética na agropecuária**, 5ª Edição Revisada. Lavras: Editora UFLA. 566 p.


Para outras culturas, ver:

**Compendium of Plant Genomes**, Springer
disponível em <https://link.springer.com/search?facet-series=%2211805%22&facet-content-type=%22Book%22>
