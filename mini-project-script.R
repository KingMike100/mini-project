#library("phyloseq"); packageVersion("phyloseq")
data <- read.csv("diabimmune_16s_t1d_metadata.csv");
summary(data)
View(data)
dim(data)
class(data)

table(data$Case_Control)
library(ggplot2)
#Bar charts represneting the count of cases in each category
table(data$Case_Control)
barplot(table(data$Case_Control))
barplot(table(data$Gender))
qplot(data$Case_Control, fill = data$Gender) + geom_bar()
qplot(data$Delivery_Route, fill = data$Gender) + geom_bar() + labs(title = "A bar graph showing the Delivery routes against Gender", x = "Delivery Route", y= "Number of people", fill = "Gender")

#association of independece between disease status and other variables
library(graphics)
table1 <- table(data$Case_Control, data$Gender)  #Disease status against the gender5
table2 <- table(data$Case_Control, data$Age_at_Collection)
table3 <- table(data$Case_Control, data$Delivery_Route)  #Disease status against delivery mode

View(table3)
mosaicplot(table1)
mosaicplot(table3)
mosaicplot(~ Case_Control + Delivery_Route, data = data, main = "Case Control Vs Delivery Route", shade = TRUE)
mosaicplot(~ Case_Control + Gender, data = data, main = "Case Control Vs Gender", shade = TRUE)
mosaicplot(~ Case_Control + Age_at_Collection, data = data, main = "Case Control Vs Age", shade = TRUE)

#Association/Dependency between Disease status & (delivery mode, gender, age)
chisq.test(table1)

#As the p-value 0.5796 is greater than the .05 significance level, we don not reject the 

chisq.test(table3)  #p value  3.949e-09 is less 

chisq.test(table2)






otuTable <- read.table("otu_table")  #Importing the OTU table
head(otuTable, n = 1)

taxaTable <- read.table("taxa_table") #Importing the Taxa table
head(taxaTable, n = 1)

dim(otuTable)

dim(taxaTable)

class(otuTable)

class(taxaTable)

#Converting the Taxa and OTU Table into a Matrix
mtaxaTable <- as.matrix(taxaTable)
motuTable <- as.matrix(otuTable)

class(mtaxaTable)
class(motuTable)

#Data Cleansing. This is done to have consistent data across all the matrices
#This will involve making sure that the OTU/taxa row names match. Currently they dont as taxa have a trailing ";"

#rownames(mtaxaTable)[rownames(mtaxaTable) == "4333897;"] = "4333897"

head(mtaxaTable, n=1)
tnames <- rownames(mtaxaTable)  #Extract rownames from the matrix
ntnames <- gsub(x = tnames, pattern = ";", replacement = "")  #Remove the ; from the extracted rownames
#ntnames
rownames(mtaxaTable) <- ntnames  #Set the new rownames
head(mtaxaTable, n=1)


library(phyloseq)
library(ggplot2)

#Tell phyloseq to load them into a phyloseq object
OTU = otu_table(motuTable, taxa_are_rows = TRUE)
TAX = tax_table(mtaxaTable)

#OTU
#TAX

#Generating the phyloseq object
physeq = phyloseq(OTU, TAX)
physeq

#Plotting the phyloseq
plot_bar(physeq, fill = "Family.")


plot_richness(physeq)  #Defauit plot produced by the plot_richness function

