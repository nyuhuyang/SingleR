############################################
# combine immgen and mouse.rnaseq
############################################
library(SingleR)
library(genefilter)
source("../R/Seurat_functions.R")
#-------skip if ref data was saved-----
data("immgen")
data("mouse.rnaseq")
head(immgen$types)
head(immgen$main_types)
head(mouse.rnaseq$types)
dim(immgen$data)
dim(mouse.rnaseq$data)
head(immgen$data[,1:5])
head(mouse.rnaseq$data[,1:5])
anyNA(immgen$data)
anyNA(mouse.rnaseq$data)
testMMM(immgen$data)
testMMM(mouse.rnaseq$data)
par(mfrow=c(1,2))
boxplot(immgen$data, main="immgen") #slow!
boxplot(mouse.rnaseq$data, main="mouse.rnaseq")#slow!



# remove low quanlity mouse.rnaseq data
par(mfrow=c(2,1))
hist(colMeans(mouse.rnaseq$data),breaks=ncol(mouse.rnaseq$data))
quantile_75 <- apply(mouse.rnaseq$data,2,function(x) quantile(x,probs =0.75))
hist(quantile_75, breaks=ncol(mouse.rnaseq$data))
rm_samples <- names(quantile_75)[quantile_75<5]
rm_index <- which(colnames(mouse.rnaseq$data) %in% rm_samples)
mouse.rnaseq_rm <- mouse.rnaseq$data[,-rm_index]
par(mfrow=c(1,1))
boxplot(mouse.rnaseq_rm,main="mouse.rnaseq_rm")#slow!
testMMM(mouse.rnaseq_rm)
# merge
immgen_mouse.rnaseq <- merge(immgen$data,mouse.rnaseq_rm,
                            by="row.names",all=FALSE)
rownames(immgen_mouse.rnaseq) = immgen_mouse.rnaseq$Row.names
immgen_mouse.rnaseq <- immgen_mouse.rnaseq[-which(colnames(immgen_mouse.rnaseq)=="Row.names")]
head(immgen_mouse.rnaseq[,1:2])
dim(immgen_mouse.rnaseq)
testMMM(immgen_mouse.rnaseq)
boxplot(immgen_mouse.rnaseq, main="immgen_mouse.rnaseq")


name = 'immgen_mouse.rnaseq'
expr = as.matrix(immgen_mouse.rnaseq) # the expression matrix
types = as.character(c(immgen$types,mouse.rnaseq$types[-rm_index])) # a character list of the types. Samples from the same type should have the same name.
main_types = as.character(c(immgen$main_types,mouse.rnaseq$main_types[-rm_index])) # a character list of the main types. 

# Fine tune types
FineTune <- function(x){
        x = gsub("DC \\(DC\\)","Dendritic cells",x)
        return(x)
}
ref_immgen_mouse.rnaseq = CreateSinglerReference(name = "immgen_mouse.rnaseq",
                                      expr = as.matrix(immgen_mouse.rnaseq),
                                      types = FineTune(types), 
                                      main_types = FineTune(main_types))

save(ref_immgen_mouse.rnaseq,file='../SingleR/data/ref_Mouse.RData') # it is best to name the object and the file with the same name.
