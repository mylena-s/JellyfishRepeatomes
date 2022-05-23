
if(!require(psych)){install.packages("psych")}
if(!require(pwr)){install.packages("pwr")}
install.packages("pwr", repos="http://cran.r-project.org")

#power test of sample size
pwr.r.test(n= 11, power = 0.80, alternative = "two.sided")


#spearman rank-test with fdr correction on "order" abundance data
data1 <- read.csv("stats/Abundance_Order_data.csv")
x1 = data1[ , c(2:11)]
y1 = data1$Genome_size
spearman1 <- corr.test(x1, y1, method="spearman", adjust="fdr", normal = FALSE)
combined_pmatrix1<- lowerUpper(spearman1$p.adj, upper=spearman1$p, diff=FALSE)
write.csv(spearman1$r, file = "stats/Abundance_ORDER_SpearmanR_2.csv")
write.csv(spearman1$p.adj, file = "stats/Abundance_ORDER_Spearmanp_adjusted_2.csv")
write.csv(spearman1$p, file = "stats/Abundance_ORDER_Spearmanp_2.csv")
write.csv(combined_pmatrix1, file = "stats/Abundance_ORDER_Spearmanp_all2.csv")


#spearman rank-test with fdr correction on "type" abundance data
data2 <- read.csv("stats/Abundance_TYPE_data.csv")
x2 = data2$Genome_size
y2 = data2[,2:4]
spearman2 <- corr.test(x2, y2, method="spearman", adjust="fdr", normal = FALSE)
combined_pmatrix2 <- lowerUpper(spearman2$p.adj, upper=spearman2$p, diff=FALSE)

write.csv(spearman2$r, file = "stats/Abundance_TYPE_SpearmanR.csv")
write.csv(spearman2$p.adj, file = "stats/Abundance_TYPE_Spearmanp_adjusted.csv")
write.csv(spearman2$p, file = "stats/Abundance_TYPE_Spearmanp.csv")
write.csv(combined_pmatrix2, file = "stats/Abundance_TYPE_Spearmanp2_all.csv")

#spearman rank-test with fdr correction on TE "order" diversity data
data3 <- read.csv("stats/Diversity_data.csv")
x3 = data3$Genome_size
y3 = data3[,c(2:5)]
spearman3 <- corr.test(x3, y3, method="spearman", adjust="fdr", normal = FALSE)
combined_pmatrix3 <- lowerUpper(spearman3$p.adj, upper=spearman3$p, diff=FALSE)

write.csv(spearman3$r, file = "stats/Diversity_SpearmanR.csv")
write.csv(spearman3$p.adj, file = "stats/Diversity_Spearmanp_adjusted.csv")
write.csv(spearman3$p, file = "stats/Diversity_Spearmanp_diversity.csv")
write.csv(combined_pmatrix3, file = "stats/Diversity_Spearmanp_all.csv")

#mod = betareg("log(Genome_size)", data = data, link = "logit", link.phi= "log")
#summary(mod)

#spearman rank-test with fdr correction on "type" absolute data
data5 <- read.csv("stats/Absolute_TYPE_data.csv")
x5 = data5
y5 = data5
spearman5 <- corr.test(x5, y5, method="spearman", adjust="fdr", normal = FALSE)
combined_pmatrix5 <- lowerUpper(spearman5$p.adj, upper=spearman5$p, diff=FALSE)

write.csv(spearman5$r, file = "stats/absolute_TYPE_SpearmanR.csv")
write.csv(spearman5$p.adj, file = "stats/absolute_TYPE_Spearmanp_adjusted.csv")
write.csv(spearman5$p, file = "stats/absolute_TYPE_Spearmanp.csv")
write.csv(combined_pmatrix5, file = "stats/absolute_TYPE_Spearmanp2_all.csv")

#spearman rank-test with fdr correction on "order" on absolute abundance data
data4 <- read.csv("stats/absolute_Order_data.csv")
x4 = data4
y4 = data4
spearman4 <- corr.test(x4, y4, method="spearman", adjust="fdr", normal = FALSE)
combined_pmatrix4<- lowerUpper(spearman4$p.adj, upper=spearman4$p, diff=FALSE)
write.csv(spearman4$r, file = "stats/absolute_ORDER_SpearmanR.csv")
write.csv(spearman4$p.adj, file = "stats/absolute_ORDER_Spearmanp_adjusted.csv")
write.csv(spearman4$p, file = "stats/absolute_ORDER_Spearmanp.csv")
write.csv(combined_pmatrix4, file = "stats/absolute_ORDER_Spearmanp_all.csv")

