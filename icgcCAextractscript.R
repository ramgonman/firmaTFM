
#### Se importa tabla "exp_seq.tsv" obtenida en julio 2023 del portal ICGC

arras <- read.csv("exp_seq.tsv", stringsAsFactors=FALSE, sep="\t")
donors <- !duplicated(arras$icgc_donor_id)
sum(donors)
donors0 <- arras$icgc_donor_id[donors]
rm(donors)

donorlist <- list()

for (i in donors0) {
donorlist[[i]] <- arras[arras$icgc_donor_id==i, c("gene_id","raw_read_count")]
}

all(donors0==names(donorlist))  ### TRUE

dimall <- lapply(donorlist, dim)

## Comprobamos si el orden de las gene_ID (ENSG) es el mismo que las muestras indicadas, se encuentran distintos números de gene_IDs (genes)
## en algunas muestras

checks0 <- vector()
for (i in donors0){
checks0[i] <- all(donorlist[[i]][,1]==donorlist[[2]][,1])
}

sum(checks0)

checks0b <- vector()
for (i in donors0){
checks0b[i] <- all(donorlist[[i]][,1]==donorlist[[1]][,1])
}

sum(checks0b)

### Obtenemos cuántos genes hay en cada muestra y los gene_ID de cada donante

veclen <- vector()

for (i in donors0){
veclen[i] <- dim(donorlist[[i]])[1]
}

table(veclen)

geneidlist <- list()
for (i in donors0){
geneidlist[[i]] <- donorlist[[i]][,1]
}

lapply(geneidlist, length)

### Se identifican los genes que están duplicados y se retienen todos los no duplicados
### La mayoría de muestras tienen 60490 genes únicos 

donorlistdup <- list()
for (i in donors0){
donorlistdup[[i]] <- donorlist[[i]][!duplicated(donorlist[[i]][,1]),]
}

lapply(donorlistdup, dim)
dimdup <- lapply(donorlistdup, dim)

vecdup <- vector()
for (i in donors0){
vecdup[i] <- dimdup[[i]][1]==60490
}

sum(vecdup)

vecdupNO <- vector()
for (i in donors0){
vecdupNO[i] <- dimdup[[i]][1]!=60490
}

### Se asignan rownames a cada muestra en la lista tras haber eliminado los duplicados.
### Al eliminarlos, escogemos en 20 donantes de la lista, que contienen recogidos más de
### un experimento de expresión, sólo uno de ellos. Se comprueba a posteriori tras haber 
### seleccionado los tumores primarios y ductales si se requiriera seleccionar uno diferente
### al escogido en caso de que estos donantes sean elegidos por ser tumores primarios
### y ductales (pero anticipamos que no fue así).

for (i in donors0){
rownames(donorlistdup[[i]]) <- donorlistdup[[i]][,1]
}

### Se identifican los genes comunes a todos las muestras, y se extraen los genes comunes 
### en una lista (donorlistdupcom) para poder disponer de una matriz con el mismo número 
### de genes por muestra.

commongeneid <- intersect(rownames(donorlistdup[[1]]), rownames(donorlistdup[[2]]))

donorlistdupcom <- list()

for (i in donors0){
donorlistdupcom[[i]] <- donorlistdup[[i]][commongeneid,]
}

all(names(donorlistdupcom)==names(donorlistdup))  ### TRUE
rm(i)

checks1 <- vector()
for (i in donors0){
checks1[i] <- all(rownames(donorlistdupcom[[i]])==commongeneid)
}
sum(checks1)

### Se pasa de la lista donorlistdupcom que contiene los genes en el mismo orden en todas las muestras
### y sin duplicados a una matriz de las mismas dimensiones

icgc_ca_raw <- matrix(NA,234,47989)
for (i in 1:234){
icgc_ca_raw[i,] <- donorlistdupcom[[i]][,2]
}

icgc_ca_rawt <- t(icgc_ca_raw)

rownames(icgc_ca_rawt) <- donorlistdupcom[[1]][,1]
colnames(icgc_ca_rawt) <- donors0

write.csv(icgc_ca_rawt, "icgc_ca_raw.csv")
save.image("icgcCAextract.RData")
savehistory()

