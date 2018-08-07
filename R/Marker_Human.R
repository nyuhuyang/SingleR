library(SingleR)
library(genefilter)
library(dplyr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
################################################
# Obtain ranked marker gene list from hpca data
################################################
data("hpca")
sort(unique(hpca$main_types))
types = FineTune_Human(as.character(hpca$types),
                       main.type = FALSE)
main_types = FineTune_Human(as.character(hpca$main_types),
                            main.type = TRUE)
Hpca = CreateSinglerReference(name = 'Hpca',
                             expr = hpca$data,
                             types = types, 
                             main_types = main_types)

hpca_main = FindAllMarkerGenes(object = Hpca, main.type = TRUE, Human =TRUE)
write.csv(hpca_main,file = "../SingleR/output/hpca_main.csv")
hpca_sub = FindAllMarkerGenes(object = Hpca, main.type = FALSE, Human =TRUE)
write.csv(hpca_sub,file = "../SingleR/output/hpca_sub.csv")

################################################
# Obtain ranked marker gene list from blueprint_encode data
################################################
# check blueprint_encode data==============================
data("blueprint_encode")
head(blueprint_encode$types)
dim(blueprint_encode$data)
head(blueprint_encode$data[,1:5])
table(is.na(blueprint_encode$data))
blueprint_encode$data[is.na(blueprint_encode$data)] = 0
head(colSums(blueprint_encode$data))
boxplot(blueprint_encode$data, main="blueprint_encode")#slow!
boxplot(blueprint_encode$data[,1:100])#slow!
# remove low quanlity blueprint_encode data
par(mfrow=c(2,1))
hist(colMeans(blueprint_encode$data),breaks=ncol(blueprint_encode$data))
quantile_75 <- apply(blueprint_encode$data,2,function(x) quantile(x,probs =0.75))
hist(quantile_75, breaks=ncol(blueprint_encode$data))
rm_samples <- names(quantile_75)[quantile_75<1]
rm_index <- which(colnames(blueprint_encode$data) %in% rm_samples)
blueprint_encode_rm <- blueprint_encode$data[,-rm_index]
quantile_75_new <- apply(blueprint_encode_rm,2,function(x) quantile(x,probs =0.75))
hist(quantile_75, breaks=ncol(blueprint_encode$data))
hist(quantile_75_new, breaks=ncol(blueprint_encode_rm),xlim = c(0,4.1))

par(mfrow=c(1,1))
boxplot(blueprint_encode_rm)#slow!
title(main="blueprint_encode_rm")

types = FineTune_Human(as.character(blueprint_encode$types[-rm_index]),
                       main.type = FALSE)
main_types = FineTune_Human(as.character(blueprint_encode$main_types[-rm_index]),
                            main.type = TRUE)

Blueprint_encode = CreateSinglerReference(name = 'Blueprint_encode',
                             expr = blueprint_encode_rm,
                             types = types, 
                             main_types = main_types)

blueprint_encode_main = FindAllMarkerGenes(object = Blueprint_encode, 
                                           main.type = TRUE, Human =TRUE)
write.csv(blueprint_encode_main,file = "../SingleR/output/blueprint_encode_main.csv")
blueprint_encode_sub = FindAllMarkerGenes(object = Blueprint_encode, 
                                          main.type = FALSE, Human =TRUE)
write.csv(blueprint_encode_sub,file = "../SingleR/output/blueprint_encode_sub.csv")

################################################
# merge ranked marker gene list hpca and blueprint_encode 
################################################
Iname1 = load(file='../SingleR/data/Hpca.RData')
Iname1
Iname2 = load(file='../SingleR/data/Blueprint_encode.RData')
Iname2
hpca_blue_encode <- merge(Hpca$data,Blueprint_encode$data,
                          by="row.names",all=FALSE)
rownames(hpca_blue_encode) = hpca_blue_encode$Row.names
hpca_blue_encode <- hpca_blue_encode[-which(colnames(hpca_blue_encode)=="Row.names")]
head(hpca_blue_encode[,1:5])
dim(hpca_blue_encode)
testMMM(hpca_blue_encode)
par(mfrow=c(1,1))
boxplot(hpca_blue_encode,main='hpca_blueprint_encode')#slow!

# Fine tune types
ref = CreateSinglerReference(name = 'hpca_blueprint_encode',
                             expr = as.matrix(hpca_blue_encode), # the expression matrix
                             types = c(Hpca$types,Blueprint_encode$types), 
                             main_types = c(Hpca$main_types,Blueprint_encode$main_types))

save(ref,file='../SingleR/data/ref_Human.RData') # it is best to name the object and the file with the same name.

Hpca_Blueprint_encode_main = FindAllMarkerGenes(object = ref, 
                                           main.type = TRUE, Human =TRUE)
write.csv(Hpca_Blueprint_encode_main,file = "../SingleR/output/Hpca_Blueprint_encode_main.csv")
Hpca_Blueprint_encode_sub = FindAllMarkerGenes(object = ref, 
                                          main.type = FALSE, Human =TRUE)
write.csv(Hpca_Blueprint_encode_sub,file = "../SingleR/output/Hpca_Blueprint_encode_sub.csv")