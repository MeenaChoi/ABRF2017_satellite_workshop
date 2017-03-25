##############################
##############################
## ABRF 2017 Satellite workshop - Hands-on 1 : Introduction to R and experimental design
## Date: March 25, 2017
## Created by Meena Choi and Ting Huang
##############################
##############################

# Summary

# Creating a new RStudio project
# Reading in data in R
# Data exploration, subsetting and replacement
# Visualizing data
# Select random sample and randomize MS run orders.
# Saving your work

#############################################
## 1. Create a new Rstudio project
#############################################

# From the menu, select **File > New Project...**, 
# then select **Existing Directory** and choose the directory where you downloaded this script and the example datasets for this tutorial. 
# All the output files we'll be creating in this tutorial will be saved in the 'working directory' that now has been set by Rstudio.
  
# Let's verify the working directory path with the get working directory command.
getwd()


#############################################
## 2. Reading in data
#############################################

iprg <- read.csv("iPRG_example_runsummary.csv")


#############################################
## 3. Data exploration
#############################################

#`class` shows the type of a variable, in this case a 'data.frame`.
class(iprg)

# `dim` shows the dimension of a data.frame, which are the number of rows and the number of columns
dim(iprg)

#`colnames` is short for column names. 
colnames(iprg)

#`head` shows the first 6 rows of data. Try `tail` to show the last 6 rows of data.
head(iprg)

# Let's explore the type of every column/variable and a summary for the value range for every column.
summary(iprg)

# Inspect the possible values for the `Conditions` and the `BioReplicate` (8th) column 
# using the named and numbered column selection syntax for data frames.
unique(iprg[, 'Condition'])
unique(iprg[, 4])

unique(iprg[, c('Condition', 'BioReplicate', 'Run')])

# Select subsets of rows from iprg dataset: 
# i.e we might be interested in working from here on only with the Condition1
# or all measurements on one particular MS run.

iprg.condition1 <- iprg[iprg$Condition == 'Condition1', ]
iprg.condition1.bio1 <- iprg[iprg$Condition == 'Condition1' & iprg$BioReplicate == '1', ]
nrow(iprg.condition1.bio1)

# subset of data for condition1 or condition2
iprg.condition1.2 <- iprg[iprg$Condition == 'Condition1' 
                          | iprg$Condition == 'Condition2', ]
nrow(iprg.condition1.2)
unique(iprg.condition1.2$Condition)

# subset of data for condition1 or condition2
iprg.condition1.2 <- iprg[which(iprg$Condition %in% c('Condition1', 'Condition2')), ]
nrow(iprg.condition1.2)
unique(iprg.condition1.2$Condition)


#############################################
## 4. Summarizing and Visualizing data
#############################################

### 4.1 Histogram
# Make a histogram of all the MS1 intensities, quantified by Skyline, for `iPRG_example`.

hist(iprg$Intensity)

# Our histogram looks quite skewed. How does this look on log-scale? 
# Do you recognize this distribution? The distribution for log2-transformed intensities looks very similar to the normal distribution. 

hist(iprg$Log2Intensity, 
     xlab="log2 transformed intensities", main="Histogram of iPRG data")

hist(log2(iprg$Intensity), 
     xlab="log2 transformed intensities", main="Histogram of iPRG data")

# We look at the summary for the log2-transformed values including the value for the mean.
summary(iprg$Log2Intensity)


### 4.2 Boxplot or box-and-whisker plot

# Boxplots are extremely useful becasue they allow us to quickly visualize the data distribution, 
# without making assumptions of the distribution type (non-parametric). 

# Let's make the boxplot with `ggplot2`, one of the most popular and powerful R packages 
# for making graphics. 

# Load ggplot2
library(ggplot2)

ggplot(aes_string(x='Run', y='Log2Intensity'), data=iprg)+
  geom_boxplot(aes_string(fill='Condition'))

# Let's rename all axis labels and title, and rotate the x-axis labels 90 degrees. 
# We can add those specifications using the `labs` and `theme` functions of the `ggplot2` package.

ggplot(aes_string(x='Run', y='Log2Intensity'), data=iprg)+
  geom_boxplot(aes_string(fill='Condition'))+
  labs(title='Log2 transformed intensity distribution per MS run', 
       y='Log2(Intensity)',
       x='MS run')+
  theme(axis.text.x=element_text(angle=90))

# And easily switch from a boxplot to a violin plot representation by changing the `geom` type. 

ggplot(aes_string(x='Run', y='Log2Intensity'), data=iprg)+
  geom_violin(aes_string(fill='Condition'))+
  labs(title='Log2 transformed intensity distribution per Subject', 
       y='Log2(Intensity)',
       x='MS run')+
  theme(axis.text.x=element_text(angle=90))


#############################################
# 5. Randomization
#############################################

## 5.1 Random selection of samples from a larger set

# This particular dataset contains a total of 10 subjects across conditions. 
# Suppose we label them from 1 to 10 
# and randomly would like to select 3 subjects we can do this using the `sample` function. 
# When we run `sample` another time, different subjects will be selected. Try this a couple times.

sample(10, 3)
sample(10, 3)

# Now suppose we would like to select the same randomly selected samples every time, 
# then we can use a random seed number.

set.seed(3728)
sample(10, 3)

set.seed(3728)
sample(10, 3)


## 5.2 Completely randomized order of MS runs 

# We can also create a random order using all elements of iPRG dataset. 
# Again, we can achieve this using `sample`, asking for exactly the amount of samples in the subset. 
# This time, each repetition gives us a different order of the complete set.

msrun <- unique(iprg$Run)
msrun

# randomize order among all 12 MS runs
sample(msrun, length(msrun))

# different order will be shown.
sample(msrun, length(msrun))


## 5.3 Randomized block design

## Allow to remove known sources of variability that you are not interested in.
## Group conditions into blocks such that the conditions in a block are as similar as possible.
## Randomly assign samples with a block.

# This particular dataset contains a total of 12 MS runs across 4 conditions, 
# 3 technical replicates per condition. Using the `block.random` function in the `psych` package, 
# we can achieve randomized block designs!

# use 'psych' package. Load the package first.
library(psych)

msrun <- unique(iprg[, c('Condition','Run')])
msrun

# 4 Conditions of 12 MS runs randomly ordered
block.random(n=12, c(Group=4))


#############################################
# 6. Saving your work 
#############################################

# You can save plots to a number of different file formats. 
# PDF is by far the most common format because it's lightweight, cross-platform 
# and scales up well but jpegs, pngs and a number of other file formats are also supported. 
# Let's redo the last barplot but save it to the file system this time. 

# Let's save the boxplot as pdf file.
pdf('boxplot_log2intensity_distribution_byMSrun.pdf', width=10, height=8)
ggplot(aes_string(x='Run', y='Log2Intensity'), data=iprg)+
  geom_boxplot(aes_string(fill='Condition'))+
  labs(title='Log2 transformed intensity distribution per MS run', 
       y='Log2(Intensity)',
       x='MS run')+
  theme(axis.text.x=element_text(angle=90))
dev.off() 

# Finally, we can save this whole session you worked so hard on! 
  
save.image(file='Day1.RData')

# let's give it a rest for today. Saving an .RData is the easiest way to pick up your work right where you left it!

rm(list=ls())
load(file = 'Day1.RData')



