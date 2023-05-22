setwd("/nfs/users/rg/dgarrido/PhD/projects/sqtlseeker/paper/simulations/real_data/gemma-nf")
pheno <- read.table("data/phenotypes.tsv", h = T, sep = "\t", check.names = F)
cov <- read.table("data/covariates.tsv", h = T)
cov <- cov[ ,1:4]
cov$gender <- as.numeric(as.factor(cov$gender))-1
comm.inds <- as.character(intersect(pheno$ID, cov$ID))
rownames(pheno) <- pheno$ID
rownames(cov) <- cov$ID
pheno <- pheno[comm.inds, ]
cov <- cov[comm.inds, ]
# identical(cov$ID, pheno$ID)
pheno$ID <- cov$ID <- NULL
fit <- lm(as.matrix(pheno) ~ . , data = cov)
resid <- data.frame(ID=rownames(resid),fit$residuals, check.names = F)
write.table(resid, "data/phenotypes.resid.tsv", sep="\t", quote=F, row.names = F)
