#############################################
# 1. Statistical hypothesis test in R
## Two sample t-test for one protein with one feature
#############################################

#load data from Section 2
load("Section2.RData")

# Let's start with one protein, named "sp|P44015|VAC2_YEAST"
oneproteindata <- iprg[iprg$Protein == "sp|P44015|VAC2_YEAST", ]

# Now, we'll perform a t-test whether protein `sp|P44015|VAC2_YEAST` has a change in abundance 
# between Condition 1 and Condition 2.

# If you want to see more details, 
?t.test

# First, get two conditions only, because t.test only works for two groups (conditions).
oneproteindata.condition12 <- oneproteindata[which(oneproteindata$Condition %in% 
                                                     c('Condition1', 'Condition2')), ]
unique(oneproteindata.condition12$Condition)
unique(oneproteindata$Condition)

# t test for different abundance (log2Int) between Groups (Condition)
result <- t.test(oneproteindata.condition12$Log2Intensity ~ oneproteindata.condition12$Condition,
                 var.equal=TRUE)
# show the summary of t-test including confidence level with 0.95
result

# We can redo the t-test and change the confidence level for the log2 fold change.

result.ci90 <- t.test(oneproteindata.condition12$Log2Intensity ~ oneproteindata.condition12$Condition, 
                      var.equal=TRUE, 
                      conf.level=0.9)
result.ci90

# Let's have a more detailed look at what information we can learn from the results our t-test. 

# name of output
names(result)

# mean for each group
result$estimate 

# log2 transformed fold change between groups : Disease-Healthy
result$estimate[1]-result$estimate[2]

# test statistic value, T value
result$statistic 

# standard error
(result$estimate[1]-result$estimate[2])/result$statistic

# degree of freedom
result$parameter 

# p value for two-sides testing
result$p.value 

# 95% confidence interval for log2 fold change
result$conf.int 

# p value calculation for one side
1-pt(result$statistic, result$parameter)

# p value for two sides, which is the same as pvalue from t test (result$p.value)
2*(1-pt(result$statistic, result$parameter))

# We can also manually compute our t-test statistic using the formulas we descibed above and 
# compare it with the `summaryresult` 

summaryresult

summaryresult12 <- summaryresult[1:2, ]

# test statistic, It is the same as 'result$statistic' above.
diff(summaryresult12$mean) # same as result$estimate[1]-result$estimate[2]
sqrt(sum(summaryresult12$sd^2/summaryresult12$length)) # same as stand error

diff(summaryresult12$mean)/sqrt(sum(summaryresult12$sd^2/summaryresult12$length))

#############Sample size calculation#######
install.packages("pwr")
library(pwr)

?pwr.t.test

# Significance level alpha
alpha <- 0.05

# Power = 1 - beta
power <- 0.95

# anticipated log2 fold change 
delta <- 1

# anticipated variability
sigma <- 1.5

# Effect size
d <- delta/sigma
pwr.t.test(d = 0, sig.level = alpha, power = power, type = 'two.sample')



# anticipated log2 fold change 
delta <- seq(0.1, 0.7, .1)
nd <- length(delta)

# anticipated variability
sigma <- seq(0.1,0.5,.1)
ns <- length(sigma)

# obtain sample sizes
samsize <- matrix(0, nrow=ns*nd, ncol = 3)
counter <- 0
for (i in 1:nd){
  for (j in 1:ns){
    result <- pwr.t.test(d = delta[i]/sigma[j],
                         sig.level = alpha, power = power,
                         type = "two.sample")
    counter <- counter + 1
    samsize[counter,1] <- delta[i]
    samsize[counter,2] <- sigma[j]
    samsize[counter,3] <- ceiling(result$n)
  }
}
colnames(samsize) <- c("fd","var","value")


library(ggplot2)
samsize <- as.data.frame(samsize)
samsize$var <- as.factor(samsize$var)
ggplot(data=samsize, aes(x=fd, y=value, group = var, colour = var)) +
  geom_line() +
  geom_point(size=2, shape=21, fill="white") +
  labs(title="Sig=0.05 Power=0.05", x="Anticipated log2 fold change", y='Sample Size (n)') +
  theme(plot.title = element_text(size=20, colour="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=13)) 

#############################################
# 2. Comparing proportion in R
## Two-Sample Proportion Test
#############################################
#Dataset is from nature paper: Proteogenomic characterization of human colon and rectal cancer (Zhang et al. 2014)
#Load in the TCGA colorectal cancer sample informtaion 
TCGA.CRC <- read.csv("TCGA_sample_information.csv")

#Invoke a spreadsheet-style data viewer on a matrix-like R object
View(TCGA.CRC)

#`colnames` is short for column names. 
#
colnames(TCGA.CRC)

# Select columns from TCGA dataset: 
# We are interested in the cancer type and history of colon polyps
TCGA.CRC.gc <- TCGA.CRC[, c('Cancer', 'history_of_colon_polyps')]
nrow(TCGA.CRC.gc)

###########Two-Sample Proportion Test################
#Generate 2-way contingency tables
ov <- table(TCGA.CRC.gc)
ov

#dotchart
dotchart(t(ov), xlab="Observed counts")

#Hypothesis: whether the proportion of patients who have history of colon polyps in the patients with colon cancer is different from that in the patients with rectal cancer
?prop.test
pt <- prop.test(ov)
pt

# name of output
names(pt)

# proportion in each group
pt$estimate 

# test statistic value
pt$statistic 

# degree of freedom
pt$parameter

# p value for two-sides testing
pt$p.value 

# double check p-value
1-pchisq(pt$statistic,pt$parameter)



#Fisher's Exact Test
?fisher.test

ft <- fisher.test(ov) 
ft

# odds ratio
ft$estimate 

#############################################
# 3. Saving your work 
#############################################

save.image(file='Section3.RData')
