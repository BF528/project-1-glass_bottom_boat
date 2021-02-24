# programmer project 1 bf528 glass bottom boat
library("affy")
library("affyPLM")
library("sva")
library("AnnotationDbi")
library("hgu133plus2.db")

data_missing_one <- ReadAffy(celfile.path = "/projectnb/bf528/users/glass_bottom_boat/project_1/samples/CEL_files")
#dim(data_missing_one)
#str(data_missing_one)

#just for the time being I know this is not what they're looking for
missing_file <- ReadAffy(celfile.path = "/projectnb/bf528/users/glass_bottom_boat/project_1/programmer/CEL2")
#missing_file <- ReadAffy(celfile.path = "/projectnb/bf528/users/glass_bottom_boat/project_1/samples/GSM971958_JS_04_U133_2.CEL.gz")

raw_data <- merge(data_missing_one, missing_file)

# normalise data
rma_data <- rma(raw_data)

write.csv(rma_data, "rma_data_norm.csv")

boxplot(raw_data, names = NULL, xlab = "Samples", ylab = "log-intensity", main = "Raw Data")
boxplot(rma_data, names = NULL, xlab = "Samples", ylab = "log-intensity", main = "Normalized Data")

###
##
# affyPLM
##
###

pset <- fitPLM(raw_data, normalize = TRUE, background = TRUE)

# RLE
rle <- RLE(pset, names = NULL, xlab = "Samples", ylab = "RLE values", main = "RLE boxplot")
rle_stats <- RLE(pset, type = "stats")
rle_plot <- RLE(pset, type = "density", names = NULL, xlab = "Samples", ylab = "RLE values", main = "RLE density plot")

# NUSE
nuse <- NUSE(pset, names = NULL, xlab = "Samples", ylab = "NUSE values", main = "NUSE boxplot")
nuse_stats <- NUSE(pset, type = "stats")
nuse_plot <- NUSE(pset, type = "density", names = NULL, xlab = "Samples", ylab = "NUSE values", main = "NUSE density plot")

###
##
# sva
##
###

# read clincal and batching annotation file
metadata <- read.csv(file = "/project/bf528/project_1/doc/proj_metadata.csv", header = TRUE)

#head(metadata)
#dim(metadata)
#names(metadata)

####
####

Pheno = pData(rma_data)
edata = exprs(rma_data)
mod = model.matrix(~as.factor(metadata$normalizationcombatmod), data = Pheno)
mod0 = model.matrix(~1, data = Pheno)

Batch = metadata$normalizationcombatbatch
modcombat = model.matrix(~1, data = Pheno)
combat_edata = ComBat(dat = edata, batch = Batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)

pValuesComBat = f.pvalue(combat_edata, mod, mod0)
qValuesComBat = p.adjust(pValuesComBat, method = "BH")

####
####

write.csv(combat_edata, "tmp.csv", row.names = TRUE, col.names = TRUE)
#this is analyst file

###
##
# PCA
##
###

transposed <- t(combat_edata)
num_transp2 <- scale(transposed, center = TRUE , scale = TRUE)
re_transposed <- t(num_transp2)

pc_analysis <- prcomp(re_transposed, center = FALSE, scale = FALSE)

#The rotation measure provides the principal component loading. 
#Each column of rotation matrix contains the principal component loading vector.

rot <- pc_analysis$rotation
rot.df <- data.frame(rot) # only way I could make the plot

#means
pca_mean <- pc_analysis$center
boxplot(pca_mean)
plot(pca_mean)
#standard dev
pca_sd <- pc_analysis$sdev
pca_sd
boxplot(pca_sd)


plot(rot.df$PC1, rot.df$PC2, xlab = "PC1", ylab = "PC2", main = "PCA")
boxplot(rot.df$PC1, rot.df$PC2)

# still need to look for outliers with $importance
imp <-rot.df$importance
plot(imp)

plot(summary(rot.df$PC1, rot.df$PC2[])$importance)
