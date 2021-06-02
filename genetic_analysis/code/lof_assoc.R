
############  Run phenotype associations for pLoF of genes ############
my_formula<-"phenotype ~ lof_score + Age + Sex + Array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"

# data is a dataframe includnig Age, Sex, Array, top 10 PCs, processed cardiac imaging phenotype, and lof_score for each gene
# lof_score = 1, if the individual is the carrier of at least one pLoF allele
# lof_score = 0, if the individual is the carrier of none pLoF allele
lm(as.formula(my_formula), data = data)





