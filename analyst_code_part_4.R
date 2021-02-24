#path to sample data
sampledir = "/project/bf528/project_1/data/example_intensity_data.csv"
datadir = "/projectnb/bf528/users/glass_bottom_boat/project_1/programmer/exprs_data.csv"

#read in data
data = read.csv(datadir, header=TRUE, sep=",", row.names = 1)

#How many genes/probes are there?
totalgenes <- dim(data)[1]

###Noise Filtering and Dimensionality Reduction
#Part 1. Expressed in at least 20% of samples 
#For each gene(rows), at least 20% of the gene-expression values must be > log2(15))

threshold <- log2(15)
expressed20 <- data[rowSums(data > threshold) >= (0.2*ncol(data)), ]

#How many genes are left after filtering?
filteredgenes <- dim(expressed20)[1]

#Part 2. Have a variance significantly different from the median variance of all probe sets
#use threshold of p<0.01 

#find variance of each row, use to calculate test statistic
#degrees of freedom = N-1

dof <- ncol(expressed20)-1
expressed20$variance <- apply(expressed20,MARGIN = 1,var) 
expressed20$tstat <- (dof*expressed20$variance)/median(expressed20$variance)

#chi squared test
chsq<- qchisq(0.01, dof, lower.tail = F)
chi_filter<-subset(expressed20, tstat > chsq)

#write to file expression matrix for genes filtered out so far
write.csv(chi_filter, '/projectnb/bf528/users/glass_bottom_boat/project_1/analyst/expression_filter_4.2.csv')

#Part 3. Have a coefficient of variation > 0.186.
var_filter<-subset(chi_filter, variance > 0.186)

#Part 4. Write file containing the gene expression matrix for genes passing all three of the filters from 4.1, 4.2, and 4.3.
write.csv(var_filter, '/projectnb/bf528/users/glass_bottom_boat/project_1/analyst/expression_all_filters.csv')

#Final number of genes that pass all of the thresholds
final_genes <- nrow(var_filter)

#remove extra columns for making heatmap later
var_filter.2<-subset(var_filter, select = -c(variance, tstat))

