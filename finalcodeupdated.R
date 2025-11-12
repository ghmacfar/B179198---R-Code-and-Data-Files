rm(list=ls())

#This loads the necessary packages needed to run this code. They need to be installed first!
library(ggplot2) 
library(data.table)
library(dplyr)
library(brunnermunzel)
library(ggpubr)


data <- read.csv("fulldata.csv") #loads the combined data into R under the name 'data'
newdata <- as.data.table(data)  #loads that data as a data table

zeroweekdata <- newdata[Biomarker %like% "0weeks"] #creates a new set of data only including rows from week 0

zeroweekdata$Sex <- as.factor(zeroweekdata$Sex) #makes sex into a categorical variable 

zeroweekmale <- zeroweekdata[Sex %like% "1"] #creates a new set of data only including males at week 0

zeroweekfemale <- zeroweekdata[Sex %like% "2"] #creates a new set of data only including females at week 0 

#These tests check if each biomarker is significantly different from a normal distribution
shapiro.test(zeroweekdata$IL.8) #Significantly different
shapiro.test(zeroweekdata$VEGF.A) #Not significantly different (but close)
shapiro.test(zeroweekdata$OPG) #Significantly different
shapiro.test(zeroweekdata$TGF.beta.1) #Significantly different
shapiro.test(zeroweekdata$IL.6) #Significantly different
shapiro.test(zeroweekdata$CXCL9) #Significantly different
shapiro.test(zeroweekdata$CXCL1) #Significantly different
shapiro.test(zeroweekdata$IL.18) #Not significantly different
shapiro.test(zeroweekdata$CSF.1) #Significantly different


#These density plots are used to compare the distribution of each biomarker between males and females
ggdensity(zeroweekmale$IL.8)  
ggdensity(zeroweekfemale$IL.8)

ggdensity(zeroweekmale$VEGF.A)
ggdensity(zeroweekfemale$VEGF.A)

ggdensity(zeroweekmale$OPG)
ggdensity(zeroweekfemale$OPG)

ggdensity(zeroweekmale$TGF.beta.1)
ggdensity(zeroweekfemale$TGF.beta.1)

ggdensity(zeroweekmale$IL.6)
ggdensity(zeroweekfemale$IL.6)

ggdensity(zeroweekmale$CXCL9)
ggdensity(zeroweekfemale$CXCL9)

ggdensity(zeroweekmale$CXCL1)
ggdensity(zeroweekfemale$CXCL1)

ggdensity(zeroweekmale$IL.18)
ggdensity(zeroweekfemale$IL.18)

ggdensity(zeroweekmale$CSF.1)
ggdensity(zeroweekfemale$CSF.1)

#These tests check if there is a significant difference in the distribution of biomarker levels between males and females
brunnermunzel.test(IL.8 ~ Sex, data = zeroweekdata)
brunnermunzel.test(VEGF.A ~ Sex, data = zeroweekdata) #Significant
brunnermunzel.test(OPG ~ Sex, data = zeroweekdata)
brunnermunzel.test(TGF.beta.1 ~ Sex, data = zeroweekdata) #Significant
brunnermunzel.test(IL.6 ~ Sex, data = zeroweekdata)
brunnermunzel.test(CXCL9 ~ Sex, data = zeroweekdata)
brunnermunzel.test(CXCL1 ~ Sex, data = zeroweekdata) #Significant
brunnermunzel.test(IL.18 ~ Sex, data = zeroweekdata)
brunnermunzel.test(CSF.1 ~ Sex, data = zeroweekdata) #Signifcant


#This code creates boxplots to visually represent the differences in biomarker levels between males and females. This
#is done for the four biomarkers (VEGF-A, TGF-beta-1, CXCL1, and CSF-1) found to be significant in 
#Brunner-Munzel testing.
ggplot(zeroweekdata, aes(x = factor(Sex, labels=c("Male","Female")), y = VEGF.A, fill = factor(Sex, labels=c("Male", "Female")))) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  labs(title = "VEGF-A levels by Sex at Week 0",
       x = "Sex", y = "VEGF-A level") +
  theme_bw(base_size=12) +
  scale_fill_manual(values = c("Male" = "lightblue", "Female" = "pink")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size=12, face = "bold")
  )


ggplot(zeroweekdata, aes(x = factor(Sex, labels=c("Male","Female")), y = TGF.beta.1, fill = factor(Sex, labels=c("Male", "Female")))) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  labs(title = "TGF-beta-1 levels by Sex at Week 0",
       x = "Sex", y = "TGF-beta-1 level") +
  theme_bw() +
  scale_fill_manual(values = c("Male" = "lightblue", "Female" = "pink")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size=12, face = "bold")
  )

ggplot(zeroweekdata, aes(x = factor(Sex, labels=c("Male","Female")), y = CXCL1, fill = factor(Sex, labels=c("Male", "Female")))) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  labs(title = "CXCL1 levels by Sex at Week 0",
       x = "Sex", y = "CXCL1 level") +
  theme_bw() +
  scale_fill_manual(values = c("Male" = "lightblue", "Female" = "pink")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size=12, face = "bold")
  )

ggplot(zeroweekdata, aes(x = factor(Sex, labels=c("Male","Female")), y = CSF.1, fill = factor(Sex, labels=c("Male", "Female")))) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  labs(title = "CSF-1 levels by Sex at Week 0",
       x = "Sex", y = "CSF-1 level") +
  theme_bw() +
  scale_fill_manual(values = c("Male" = "lightblue", "Female" = "pink")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size=12, face = "bold")
  )

#--------------------------------REGRESSION TIME--------------------------------------------------------------

regdata <- read.csv("inclusiondata.csv") #Loads the datasheet required to run the linear model, which contains
#biomarker levels and covariates for each patient at week 0

newreg <- as.data.table(regdata) #Stores that datasheet as a data table

twentyreg <- newreg[sample(.N, 23),] #Takes a random set of 23 rows out of the lm dataset. The number 23 was chosen
#because it is approximately 20% of 117, which is the number of rows in the dataset

eightyreg <- newreg[!twentyreg, on = names(newreg)] #Loads the 80% of rows that were not loaded into twentyreg into
#a new dataset

eightyreg$Sex <- as.factor(eightyreg$Sex) #Changes sex to a categorical variable in eightyreg 
eightyreg$Smoker <- as.factor(eightyreg$Smoker) #Changes smoker to a categorical variable in eightyreg

twentyreg$Sex <- as.factor(twentyreg$Sex) #Changes sex to a categorical variable in twentyreg
twentyreg$Smoker <- as.factor(twentyreg$Smoker) #Changes smoker to a categorical variable in twentyreg

#This runs the linear model with 12-month VAS as the response variable and biomarkers/covariates as the predictors
x <- lm(VAS12 ~ VAS0 + Smoker + Sex + Age + CSF.1 + IL.18 + CXCL1 + CXCL9 + IL.6 + TGF.beta.1 + OPG + VEGF.A + IL.8, data = eightyreg)

str(eightyreg) #This just checks if all the predictors are the right kind of variable

summary(x) #Provides a summary of the linear model

twentyreg$predicted_VAS12 <- predict(x, newdata=twentyreg) #Predicts 12-month VAS values using the linear model
#and stores the predicted values for each patient in a new column called "Preducted_VAS12"

#This creates a plot showing the correlation between predicted and actual 12-month VAS scores, 
#with a red line that represents a perfect correlation
ggplot(twentyreg, aes(x=predicted_VAS12, y=VAS12)) +
  geom_point(color="steelblue") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
theme_minimal()


#This just calculated the R-squared value associated with the correlation between predicted and actual 12-month 
#VAS scores, and stores it in a variable called "compare"
compare <- cor(twentyreg$VAS12, twentyreg$predicted_VAS12, use ="complete.obs")^2

compare #Displays the R-squared value between preducted and actual 12-month VAS



#This code is essentially doing the same as the above code associated with the regression, but runs it 1000
#times instead of only once and stores the R-squared value in a list called "rsquare" after each repetition. 
#This allows for an average R-squared between predicted and actual 12-month VAS to be calculated, to account for 
#the fact that the 80:20 split of the original dataset may have been skewed, which would subsequently skew the 
#linear model. 

iterations <- 1000
rsquare <- numeric(iterations)

for(i in 1:iterations){
  twentyreg <- newreg[sample(.N, 23),]
  eightyreg <- newreg[!twentyreg, on = names(newreg)]
  
  eightyreg$Sex <- as.factor(eightyreg$Sex)
  eightyreg$Smoker <- as.factor(eightyreg$Smoker)
  
  twentyreg$Sex <- as.factor(twentyreg$Sex)
  twentyreg$Smoker <- as.factor(twentyreg$Smoker)
  
  x <- lm(VAS12 ~ VAS0 + Smoker + Sex + Age + CSF.1 + IL.18 + CXCL1 + CXCL9 + IL.6 + TGF.beta.1 + OPG + VEGF.A + IL.8, data = eightyreg)
  
  twentyreg$predicted_VAS12 <- predict(x, newdata=twentyreg)
  
  rsquare[i] <- cor(twentyreg$VAS12, twentyreg$predicted_VAS12, use ="complete.obs")^2
  
}

mean(rsquare) #Finds the mean of R-squared values from each repititon
sd(rsquare) #Calculates the standard deviation of R-squared values from the 1000 repititions

