###Hierarchical Clustering and Subtype Discovery
#Perform hierarchical clustering on your fully filtered data matrix from Part 4.4
clust1<-hclust(dist(t(var_filter.2)))

#Part 2. Cut the dendrogram such that the samples are divided into 2 clusters
clusters<-cutree(clust1, 2)

#How many samples are in each cluster?
numclust1 <- sum(clusters==1)
numclust2 <- sum(clusters==2)

#Part 4: Identifying differentially expressed genes between the two clusters
#Use expression matrix from Part 4.4 and the cluster memberships from Part 5.2
ttest <- apply(as.matrix(var_filter.2),MARGIN=1,function(x) t.test(x=x[clusters==1],y=x[clusters==2]))

#Write out a dataframe containing the probeset ID, t-statistic, p-value, and adjusted p-value 
pvals <- sapply(ttest,function(x) x$p.value)
tstats <- sapply(ttest,function(x) x$statistic)
p_adjusted <- p.adjust(pvals ,method = "fdr")
t_p_df <- data.frame("Probeset_ID" = c(row.names(var_filter.2)),
                     tstats,pvals,p_adjusted)

#save data in csv file
write.csv(t_p_df,sep = ",",row.names = FALSE,file = "analysis5.4.csv")

#How many genes are differentially expressed at adjusted p<0.05 between the clusters for both lists?
diff_expressed <- t_p_df$Probeset_ID[t_p_df$p_adjusted<0.05]
print(c("Number of probes left with adjusted p value < 0.05:",length(diff_expressed)))


#For biologist: perform same analysis on expression matrix from 4.5 and provide as csv
ttest_bio <- apply(as.matrix(chi_filter),MARGIN=1,function(x) t.test(x=x[clusters==1],y=x[clusters==2]))
pvals_bio <- sapply(ttest_bio,function(x) x$p.value)
tstats_bio <- sapply(ttest_bio,function(x) x$statistic)
p_adjusted_bio <- p.adjust(pvals_bio,"fdr")
biologist_4.5 <- data.frame("Probeset_ID" = c(row.names(chi_filter)),
                            tstats_bio,pvals_bio,p_adjusted_bio)

diff_expressed_bio <- biologist_4.5$Probeset_ID[biologist_4.5$p_adjusted_bio<0.05]
numdiff <- length(diff_expressed_bio)
write.csv(biologist_4.5,row.names = F,file = "analysis5.6.csv",sep=",")

