##############################
##############################
## ABRF 2017 Satellite workshop - Hands-on 2 : R markdown and simple statistics in R
## Date: March 25, 2017
## Created by Meena Choi and Ting Huang
##############################
##############################


# Summary

# Calculate simple statistics and visualize them using ggplot2.
# Statistical hypothesis testing by t-test.
# Saving your work


#############################################
# 1. Basic statistical summaries in R
#############################################

# Load data from previous section
load(file = 'Section1.RData')

## 1.1 Calculate simple statistics

# Let's start with one protein as an example 
# and calculate the mean, standard deviation, standard error of the mean across all replicates per condition. 
# We then store all the computed statistics into a single summary data frame for easy access.

# We can use the **aggregate** function to compute summary statistics

# check what proteins are in dataset, show all protein names
unique(iprg$Protein)

# Let's start with one protein, named "sp|P44015|VAC2_YEAST"
oneproteindata <- iprg[iprg$Protein == "sp|P44015|VAC2_YEAST", ]

# If you want to see more details, 
?aggregate

### 1.1.1 Calculate mean per groups
# splits 'oneproteindata' into subsets by 'Condition', 
# then, compute 'FUN=mean' of 'log2Int'
sub.mean <- aggregate(Log2Intensity ~ Condition, data=oneproteindata, FUN=mean)
sub.mean

### 1.1.2 Calculate SD(standard deviation) per groups
# The same as mean calculation above. 'FUN' is changed to 'sd'.
sub.sd <- aggregate(Log2Intensity ~ Condition, data=oneproteindata, FUN=sd)
sub.sd

### 1.1.3 Count the number of observation per groups
# The same as mean calculation. 'FUN' is changed 'length'.
sub.len <- aggregate(Log2Intensity ~ Condition, data=oneproteindata, FUN=length)
sub.len

### 1.1.4 Calculate SE(standard error of mean) per groups
# SE = sqrt(s^2 / n)
sub.se <- sqrt(sub.sd$Log2Intensity^2/sub.len$Log2Intensity)
sub.se

# make the summary table including the results above (mean, sd, se and length).
summaryresult <- data.frame(Group=c("Condition1", "Condition2", "Condition3", "Condition4"),
                            mean=sub.mean$Log2Intensity,
                            sd=sub.sd$Log2Intensity, 
                            se=sub.se, 
                            length=sub.len$Log2Intensity)
summaryresult

## 1.2 Visualization with error bars for descriptive purpose
#‘error bars’ can have a variety of meanings or conclusions if what they represent is not precisely specified. 
# Below we provide some examples of which types of error bars are common. 
# We're using the summary of protein `sp|P44015|VAC2_YEAST` from the previous section and the `ggplot2` package as it provides a convenient way to make easily adaptable plots.

# Let's draw plots with mean and error bars

# means without any errorbar
ggplot(aes(x=Group, y=mean, colour=Group), data=summaryresult)+
  geom_point(size=3)

# Let's change a number of visual properties to make the plot more atttractive
# Let's change the labels of x-axis and y-axis and title: 
# add labs(title="Mean", x="Condition", y='Log2(Intensity)')
# Let's change background color for white : add theme_bw()
# Let's change size or color of labels of axes and title, text of x-axis : in theme
# Let's change the position of legend :'none' remove the legend
# Let's make the box for legend
# Let's remove the box for legend key.

ggplot(aes(x=Group, y=mean, colour=Group), data=summaryresult)+
  geom_point(size=3)+
  labs(title="Mean", x="Group", y='Log2(Intensity)')+
  theme_bw()+
  theme(plot.title = element_text(size=25, colour="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        legend.position = 'bottom',
        legend.background = element_rect(colour='black'),
        legend.key = element_rect(colour='white'))

# Very similar but now as a bar plot.
ggplot(aes(x=Group, y=mean, fill=Group), data=summaryresult)+
  geom_bar(position=position_dodge(), stat='identity')+
  scale_x_discrete('Group')+
  labs(title="Mean", x="Group", y='Log2(Intensity)')+
  theme_bw()+
  theme(plot.title = element_text(size=25, colour="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        legend.background = element_rect(colour='black'),
        legend.key = element_rect(colour='white'))

# For the sake of this tutorial, we'll continue adding error bars for different statistics with the point plots. 
# We'll leave it as an exercise to add error bars to the barplots. 
# Let's first add the standard deviation, then the standard error of the mean. Which one is smaller?

# mean with SD
ggplot(aes(x=Group, y=mean, colour=Group), data=summaryresult)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = mean + sd, ymin=mean - sd), width=0.1)+
  scale_x_discrete('Group')+
  labs(title="Mean with SD", x="Group", y='Log2(Intensity)')+
  theme_bw()+
  theme(plot.title = element_text(size=25, colour="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        legend.position = 'bottom',
        legend.background = element_rect(colour='black'),
        legend.key = element_rect(colour='white'))

# mean with SE
ggplot(aes(x=Group, y=mean, colour=Group), data=summaryresult)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = mean + se, ymin=mean - se), width=0.1)+
  labs(title="Mean with SE", x="Condition", y='Log2(Intensity)')+
  theme_bw()+
  theme(plot.title = element_text(size=25, colour="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        legend.position = 'bottom',
        legend.background = element_rect(colour='black'),
        legend.key = element_rect(colour='white'))

# The SE is narrow than the SD!

## 1.3 Calculate the confidence interval

# Now that we've covered the standard error of the mean and the standard deviation, 
# let's investigate how we can add custom confidence intervals (CI) for our measurement of the mean. 
# We'll add these CI's to the summary results we previously stored for protein `sp|P44015|VAC2_YEAST`

# Confidence interval : mean + or - (SE * alpha /2 { quantile of t distribution})$

# 95% confident interval
# Be careful for setting quantile for two-sided. need to divide by two for error.
# For example, 95% confidence interval, right tail is 2.5% and left tail is 2.5%.

summaryresult$ciw.lower.95 <- summaryresult$mean - qt(0.975,summaryresult$len)*summaryresult$se
summaryresult$ciw.upper.95 <- summaryresult$mean + qt(0.975,summaryresult$len)*summaryresult$se
summaryresult

# mean with 95% two-sided confidence interval
ggplot(aes(x=Group, y=mean, colour=Group), data=summaryresult)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = ciw.upper.95, ymin=ciw.lower.95), width=0.1)+
  labs(title="Mean with 95% confidence interval", x="Condition", y='Log2(Intensity)')+
  theme_bw()+
  theme(plot.title = element_text(size=25, colour="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        legend.position = 'bottom',
        legend.background = element_rect(colour='black'),
        legend.key = element_rect(colour='white'))

# Let's repeat that one more time for the 99% two-sided confidence interval. 

# mean with 99% two-sided confidence interval
summaryresult$ciw.lower.99 <- summaryresult$mean - qt(0.995,summaryresult$len)*summaryresult$se
summaryresult$ciw.upper.99 <- summaryresult$mean + qt(0.995,summaryresult$len)*summaryresult$se
summaryresult

ggplot(aes(x=Group, y=mean, colour=Group), data=summaryresult)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = ciw.upper.99, ymin=ciw.lower.99), width=0.1)+
  labs(title="Mean with 99% confidence interval", x="Condition", y='Log2(Intensity)')+
  theme_bw()+
  theme(plot.title = element_text(size=25, colour="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        legend.position = 'bottom',
        legend.background = element_rect(colour='black'),
        legend.key = element_rect(colour='white'))


#############################################
# 2. Saving your work 
#############################################

# Finally, we can save this whole session you worked so hard on! 
  
save.image(file='section2.RData')




