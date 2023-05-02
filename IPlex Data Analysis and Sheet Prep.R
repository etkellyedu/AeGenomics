#IPlex Data Analysis and Sheet Prep ####
# Clean Up #################################################
#start fresh for reproducible code #
# Clear packages
p_unload(all)  # Easier: clears all add-ons
#detach(“package:datasets”, unload = TRUE)  # For base
# Clear console
#cat(“\014")) ctrl+L
#Clear env
rm(list=ls())
# Clear mind :)

library(pacman)



library(dplyr)
library(tidyr)
library("ggfortify")
library(ggplot2)
library(dplyr)
library(devtools)
# load the necessary packages
library(readr)
library(readxl)
library(stringr)
library(plotly)

# STEP 1: Turn genotype excel sheet into a csv ####
# notes on manual manipulations of sheets
# First plate did not have 1534
# one submission did not have 1534 primers


# Define the path to .xlsx files and the folder for the new csv sheets
pathin <- "Data/RawIplex/NewFiles/"
pathout <- "Data/GenotypeCSV/"

# Get a list of all .xlsx files in the directory
files <- list.files(pathin, pattern = "\\.xlsx$", full.names = TRUE)

# Loop through each file
for (file in files) {
  
  # Read the third sheet from the .xlsx file
  df <- read_excel(file, sheet = 3)
  
  # Extract the file name without the extension
  file_name <- substr(basename(file), 1, nchar(basename(file))-5)
  
  # Save the data frame as a .csv file in a new folder using the original file name
  write.csv(df, file.path(pathout, paste0(file_name, ".csv")))
}



#1A: Check if we can bind the csvs together ####
library(tidyverse)

# Set the directory where your CSV files are located
csv_dir <- "Data/GenotypeCSV/"

# Get a list of all CSV files in the directory
csv_files <- list.files(path = csv_dir, pattern = "*.csv", full.names = TRUE)

# Create an empty data frame to store the results
results_df <- data.frame(filename = character(), num_rows = numeric(), num_cols = numeric())

# Loop through each CSV file and add its filename, number of rows, and number of columns to the results data frame
for (csv_file in csv_files) {
  # Read in the CSV file
  csv_data <- read_csv(csv_file)
  
  # Get the number of rows and columns in the CSV file
  num_rows <- nrow(csv_data)
  num_cols <- ncol(csv_data)
  
  # Create a row for the results data frame with the filename, number of rows, and number of columns
  results_row <- data.frame(filename = csv_file, num_rows = num_rows, num_cols = num_cols)
  
  # Add the results row to the results data frame
  results_df <- bind_rows(results_df, results_row)
}

# Print the results data frame
print(results_df)



# STEP 2: Make genotype data file by combining all genotype CSVs ####
# specify the directory path
# old drive directory_path <- "Google Drive/Shared drives/Attardo Lab/Projects/Insecticide Resistance in Aedes aegypti/Data/iPlex Results /Iplex Genotype CSV/"
directory_path <- "Data/GenotypeCSV/"

#Note: our test has two plates. Create master sheets for each plate, then merge at the end

#Note::: We delete the 2018 run of Dinuba samples because we updated the SNPs to include 1534. 
#However, the rerun was missing the 1534 primers, so we had to rerun those samples (DSY) for a single SNP. That data sheet was fixed manually. 
# create an empty dataframe to store the combined data for W2
combined_data2 <- data.frame()

# loop through each file in the directory
for (file_name in list.files(directory_path)) {
  if (grepl("W2", file_name) && grepl(".csv", file_name)) {
    file_path <- paste0(directory_path, '/', file_name)
    # read the csv file
    df2 <- read.csv(file_path)
    combined_data2 <- rbind(combined_data2, df2)
  }
}

# create an empty dataframe to store the combined data for W1
combined_data1 <- data.frame()

# loop through each file in the directory
for (file_name in list.files(directory_path)) {
  if (grepl("W1", file_name) && grepl(".csv", file_name)) {
    file_path <- paste0(directory_path, '/', file_name)
    # read the csv file
    df1 <- read.csv(file_path)
    combined_data1 <- rbind(combined_data1, df1)
  }
}

#merge W1 and W2
merged_data <- merge(combined_data1,combined_data2, by = "SAMPLE_NAME")

# 2A: Cleanup ####

#drop well column
merged_data <- merged_data[, -c(2,3,26,27)]

#some Ts get turned to TRUE (logical), change all to character
merged_data[merged_data == "TRUE"] <- "T"

# change all columns to characters to be sure
merged_data <- merged_data %>% mutate_all(as.character)

#replace blanks with NA
merged_data <- merged_data %>% mutate_all(na_if,"")

# clear all variables except for merged_data
rm(list = ls(pattern = "[^merged_data]"))

#write csv
write.csv(merged_data, file = "IPlexMasterSheet_27Mar2023.csv")

# STEP 3: PCA Sheet Prep ####
combinedSNPs <- merged_data
  
combinedSNPs[combinedSNPs == "TRUE"] <- "T"
combinedSNPs["X1000854266"][combinedSNPs["X1000854266"] == "T"] <- "1"
combinedSNPs["X1055614166"][combinedSNPs["X1055614166"] == "T"] <- "1"
combinedSNPs["X1179585685"][combinedSNPs["X1179585685"] == "G"] <- "1"
combinedSNPs["X1224849446"][combinedSNPs["X1224849446"] == "G"] <- "1"
combinedSNPs["X2038960557"][combinedSNPs["X2038960557"] == "A"] <- "1"
combinedSNPs["X2107208294"][combinedSNPs["X2107208294"] == "T"] <- "1"
combinedSNPs["X2377059147"][combinedSNPs["X2377059147"] == "A"] <- "1"
combinedSNPs["X3052344181"][combinedSNPs["X3052344181"] == "T"] <- "1"
combinedSNPs["X3315939224"][combinedSNPs["X3315939224"] == "A"] <- "1"
combinedSNPs["vgsc_3315931548"][combinedSNPs["vgsc_3315931548"] == "C"] <- "1"
combinedSNPs["vgsc_3315931672"][combinedSNPs["vgsc_3315931672"] == "T"] <- "1"
combinedSNPs["vgsc_3315983763"][combinedSNPs["vgsc_3315983763"] == "C"] <- "1"
combinedSNPs["vgsc_3315999297"][combinedSNPs["vgsc_3315999297"] == "G"] <- "1"
combinedSNPs["vgsc_3316014588"][combinedSNPs["vgsc_3316014588"] == "A"] <- "1"
combinedSNPs["vgsc_3316080722"][combinedSNPs["vgsc_3316080722"] == "C"] <- "1"
combinedSNPs["X1000854391"][combinedSNPs["X1000854391"] == "G"] <- "1"
combinedSNPs["X1016293828"][combinedSNPs["X1016293828"] == "T"] <- "1"
combinedSNPs["X1027746708"][combinedSNPs["X1027746708"] == "C"] <- "1"
combinedSNPs["X1070561909"][combinedSNPs["X1070561909"] == "A"] <- "1"
combinedSNPs["X1101943648"][combinedSNPs["X1101943648"] == "G"] <- "1"
combinedSNPs["X2001844389"][combinedSNPs["X2001844389"] == "T"] <- "1"
combinedSNPs["X2186624519"][combinedSNPs["X2186624519"] == "G"] <- "1"
combinedSNPs["X2400630292"][combinedSNPs["X2400630292"] == "A"] <- "1"
combinedSNPs["X2426833801"][combinedSNPs["X2426833801"] == "A"] <- "1"
combinedSNPs["X2429563021"][combinedSNPs["X2429563021"] == "G"] <- "1"
combinedSNPs["X2472373785"][combinedSNPs["X2472373785"] == "C"] <- "1"
combinedSNPs["X3000389496"][combinedSNPs["X3000389496"] == "T"] <- "1"
combinedSNPs["X3019326496"][combinedSNPs["X3019326496"] == "T"] <- "1"
combinedSNPs["X3068430261"][combinedSNPs["X3068430261"] == "G"] <- "1"
combinedSNPs["X3233038191"][combinedSNPs["X3233038191"] == "G"] <- "1"
combinedSNPs["X3233420299"][combinedSNPs["X3233420299"] == "A"] <- "1"
combinedSNPs["X3235265896"][combinedSNPs["X3235265896"] == "G"] <- "1"
combinedSNPs["X3305442338"][combinedSNPs["X3305442338"] == "A"] <- "1"
combinedSNPs["X3385240055"][combinedSNPs["X3385240055"] == "A"] <- "1"
combinedSNPs["vgsc_3315983611"][combinedSNPs["vgsc_3315983611"] == "C"] <- "1"
combinedSNPs["X1000854266"][combinedSNPs["X1000854266"] == "C"] <- "3"
combinedSNPs["X1055614166"][combinedSNPs["X1055614166"] == "C"] <- "3"
combinedSNPs["X1179585685"][combinedSNPs["X1179585685"] == "A"] <- "3"
combinedSNPs["X1224849446"][combinedSNPs["X1224849446"] == "A"] <- "3"
combinedSNPs["X2038960557"][combinedSNPs["X2038960557"] == "G"] <- "3"
combinedSNPs["X2107208294"][combinedSNPs["X2107208294"] == "G"] <- "3"
combinedSNPs["X2377059147"][combinedSNPs["X2377059147"] == "C"] <- "3"
combinedSNPs["X3052344181"][combinedSNPs["X3052344181"] == "C"] <- "3"
combinedSNPs["X3315939224"][combinedSNPs["X3315939224"] == "C"] <- "3"
combinedSNPs["vgsc_3315931548"][combinedSNPs["vgsc_3315931548"] == "A"] <- "3"
combinedSNPs["vgsc_3315931672"][combinedSNPs["vgsc_3315931672"] == "C"] <- "3"
combinedSNPs["vgsc_3315983763"][combinedSNPs["vgsc_3315983763"] == "T"] <- "3"
combinedSNPs["vgsc_3315999297"][combinedSNPs["vgsc_3315999297"] == "T"] <- "3"
combinedSNPs["vgsc_3316014588"][combinedSNPs["vgsc_3316014588"] == "T"] <- "3"
combinedSNPs["vgsc_3316080722"][combinedSNPs["vgsc_3316080722"] == "A"] <- "3"
combinedSNPs["X1000854391"][combinedSNPs["X1000854391"] == "C"] <- "3"
combinedSNPs["X1016293828"][combinedSNPs["X1016293828"] == "C"] <- "3"
combinedSNPs["X1027746708"][combinedSNPs["X1027746708"] == "T"] <- "3"
combinedSNPs["X1070561909"][combinedSNPs["X1070561909"] == "G"] <- "3"
combinedSNPs["X1101943648"][combinedSNPs["X1101943648"] == "T"] <- "3"
combinedSNPs["X2001844389"][combinedSNPs["X2001844389"] == "G"] <- "3"
combinedSNPs["X2186624519"][combinedSNPs["X2186624519"] == "A"] <- "3"
combinedSNPs["X2400630292"][combinedSNPs["X2400630292"] == "T"] <- "3"
combinedSNPs["X2426833801"][combinedSNPs["X2426833801"] == "G"] <- "3"
combinedSNPs["X2429563021"][combinedSNPs["X2429563021"] == "A"] <- "3"
combinedSNPs["X2472373785"][combinedSNPs["X2472373785"] == "T"] <- "3"
combinedSNPs["X3000389496"][combinedSNPs["X3000389496"] == "A"] <- "3"
combinedSNPs["X3019326496"][combinedSNPs["X3019326496"] == "C"] <- "3"
combinedSNPs["X3068430261"][combinedSNPs["X3068430261"] == "A"] <- "3"
combinedSNPs["X3233038191"][combinedSNPs["X3233038191"] == "A"] <- "3"
combinedSNPs["X3233420299"][combinedSNPs["X3233420299"] == "T"] <- "3"
combinedSNPs["X3235265896"][combinedSNPs["X3235265896"] == "A"] <- "3"
combinedSNPs["X3305442338"][combinedSNPs["X3305442338"] == "G"] <- "3"
combinedSNPs["X3385240055"][combinedSNPs["X3385240055"] == "G"] <- "3"
combinedSNPs["vgsc_3315983611"][combinedSNPs["vgsc_3315983611"] == "T"] <- "3"
combinedSNPs["X1000854266"][combinedSNPs["X1000854266"] == "TC"] <- "2"
combinedSNPs["X1055614166"][combinedSNPs["X1055614166"] == "TC"] <- "2"
combinedSNPs["X1179585685"][combinedSNPs["X1179585685"] == "GA"] <- "2"
combinedSNPs["X1224849446"][combinedSNPs["X1224849446"] == "GA"] <- "2"
combinedSNPs["X2038960557"][combinedSNPs["X2038960557"] == "AG"] <- "2"
combinedSNPs["X2107208294"][combinedSNPs["X2107208294"] == "TG"] <- "2"
combinedSNPs["X2377059147"][combinedSNPs["X2377059147"] == "AC"] <- "2"
combinedSNPs["X3052344181"][combinedSNPs["X3052344181"] == "TC"] <- "2"
combinedSNPs["X3315939224"][combinedSNPs["X3315939224"] == "AC"] <- "2"
combinedSNPs["vgsc_3315931548"][combinedSNPs["vgsc_3315931548"] == "CA"] <- "2"
combinedSNPs["vgsc_3315931672"][combinedSNPs["vgsc_3315931672"] == "TC"] <- "2"
combinedSNPs["vgsc_3315983763"][combinedSNPs["vgsc_3315983763"] == "CT"] <- "2"
combinedSNPs["vgsc_3315999297"][combinedSNPs["vgsc_3315999297"] == "GT"] <- "2"
combinedSNPs["vgsc_3316014588"][combinedSNPs["vgsc_3316014588"] == "AT"] <- "2"
combinedSNPs["vgsc_3316080722"][combinedSNPs["vgsc_3316080722"] == "CA"] <- "2"
combinedSNPs["X1000854391"][combinedSNPs["X1000854391"] == "GC"] <- "2"
combinedSNPs["X1016293828"][combinedSNPs["X1016293828"] == "TC"] <- "2"
combinedSNPs["X1027746708"][combinedSNPs["X1027746708"] == "CT"] <- "2"
combinedSNPs["X1070561909"][combinedSNPs["X1070561909"] == "AG"] <- "2"
combinedSNPs["X1101943648"][combinedSNPs["X1101943648"] == "GT"] <- "2"
combinedSNPs["X2001844389"][combinedSNPs["X2001844389"] == "TG"] <- "2"
combinedSNPs["X2186624519"][combinedSNPs["X2186624519"] == "GA"] <- "2"
combinedSNPs["X2400630292"][combinedSNPs["X2400630292"] == "AT"] <- "2"
combinedSNPs["X2426833801"][combinedSNPs["X2426833801"] == "AG"] <- "2"
combinedSNPs["X2429563021"][combinedSNPs["X2429563021"] == "GA"] <- "2"
combinedSNPs["X2472373785"][combinedSNPs["X2472373785"] == "CT"] <- "2"
combinedSNPs["X3000389496"][combinedSNPs["X3000389496"] == "TA"] <- "2"
combinedSNPs["X3019326496"][combinedSNPs["X3019326496"] == "TC"] <- "2"
combinedSNPs["X3068430261"][combinedSNPs["X3068430261"] == "GA"] <- "2"
combinedSNPs["X3233038191"][combinedSNPs["X3233038191"] == "GA"] <- "2"
combinedSNPs["X3233420299"][combinedSNPs["X3233420299"] == "AT"] <- "2"
combinedSNPs["X3235265896"][combinedSNPs["X3235265896"] == "GA"] <- "2"
combinedSNPs["X3305442338"][combinedSNPs["X3305442338"] == "AG"] <- "2"
combinedSNPs["X3385240055"][combinedSNPs["X3385240055"] == "AG"] <- "2"
combinedSNPs["vgsc_3315983611"][combinedSNPs["vgsc_3315983611"] == "CT"] <- "2"
combinedSNPs["X1000854266"][combinedSNPs["X1000854266"] == "CT"] <- "2"
combinedSNPs["X1055614166"][combinedSNPs["X1055614166"] == "CT"] <- "2"
combinedSNPs["X1179585685"][combinedSNPs["X1179585685"] == "AG"] <- "2"
combinedSNPs["X1224849446"][combinedSNPs["X1224849446"] == "AG"] <- "2"
combinedSNPs["X2038960557"][combinedSNPs["X2038960557"] == "GA"] <- "2"
combinedSNPs["X2107208294"][combinedSNPs["X2107208294"] == "GT"] <- "2"
combinedSNPs["X2377059147"][combinedSNPs["X2377059147"] == "CA"] <- "2"
combinedSNPs["X3052344181"][combinedSNPs["X3052344181"] == "CT"] <- "2"
combinedSNPs["X3315939224"][combinedSNPs["X3315939224"] == "CA"] <- "2"
combinedSNPs["vgsc_3315931548"][combinedSNPs["vgsc_3315931548"] == "AC"] <- "2"
combinedSNPs["vgsc_3315931672"][combinedSNPs["vgsc_3315931672"] == "CT"] <- "2"
combinedSNPs["vgsc_3315983763"][combinedSNPs["vgsc_3315983763"] == "TC"] <- "2"
combinedSNPs["vgsc_3315999297"][combinedSNPs["vgsc_3315999297"] == "TG"] <- "2"
combinedSNPs["vgsc_3316014588"][combinedSNPs["vgsc_3316014588"] == "TA"] <- "2"
combinedSNPs["vgsc_3316080722"][combinedSNPs["vgsc_3316080722"] == "AC"] <- "2"
combinedSNPs["X1000854391"][combinedSNPs["X1000854391"] == "CG"] <- "2"
combinedSNPs["X1016293828"][combinedSNPs["X1016293828"] == "CT"] <- "2"
combinedSNPs["X1027746708"][combinedSNPs["X1027746708"] == "TC"] <- "2"
combinedSNPs["X1070561909"][combinedSNPs["X1070561909"] == "GA"] <- "2"
combinedSNPs["X1101943648"][combinedSNPs["X1101943648"] == "TG"] <- "2"
combinedSNPs["X2001844389"][combinedSNPs["X2001844389"] == "GT"] <- "2"
combinedSNPs["X2186624519"][combinedSNPs["X2186624519"] == "AG"] <- "2"
combinedSNPs["X2400630292"][combinedSNPs["X2400630292"] == "TA"] <- "2"
combinedSNPs["X2426833801"][combinedSNPs["X2426833801"] == "GA"] <- "2"
combinedSNPs["X2429563021"][combinedSNPs["X2429563021"] == "AG"] <- "2"
combinedSNPs["X2472373785"][combinedSNPs["X2472373785"] == "TC"] <- "2"
combinedSNPs["X3000389496"][combinedSNPs["X3000389496"] == "AT"] <- "2"
combinedSNPs["X3019326496"][combinedSNPs["X3019326496"] == "CT"] <- "2"
combinedSNPs["X3068430261"][combinedSNPs["X3068430261"] == "AG"] <- "2"
combinedSNPs["X3233038191"][combinedSNPs["X3233038191"] == "AG"] <- "2"
combinedSNPs["X3233420299"][combinedSNPs["X3233420299"] == "TA"] <- "2"
combinedSNPs["X3235265896"][combinedSNPs["X3235265896"] == "AG"] <- "2"
combinedSNPs["X3305442338"][combinedSNPs["X3305442338"] == "GA"] <- "2"
combinedSNPs["X3385240055"][combinedSNPs["X3385240055"] == "GA"] <- "2"
combinedSNPs["vgsc_3315983611"][combinedSNPs["vgsc_3315983611"] == "TC"] <- "2"


#manual adds

SNPs_NAana["X1103349772"][SNPs_NAana["X1103349772"] == "T"] <- "1"
SNPs_NAana["X1103349772"][SNPs_NAana["X1103349772"] == "TC"] <- "2"
SNPs_NAana["X1103349772"][SNPs_NAana["X1103349772"] == "C"] <- "3"

SNPs_NAana["X1188368236"][SNPs_NAana["X1188368236"] == "C"] <- "1"
SNPs_NAana["X1188368236"][SNPs_NAana["X1188368236"] == "CT"] <- "2"
SNPs_NAana["X1188368236"][SNPs_NAana["X1188368236"] == "T"] <- "3"


combinedSNPs["X1103349772"][combinedSNPs["X1103349772"] == "T"] <- "1"
combinedSNPs["X1103349772"][combinedSNPs["X1103349772"] == "TC"] <- "2"
combinedSNPs["X1103349772"][combinedSNPs["X1103349772"] == "C"] <- "3"

combinedSNPs["X1188368236"][combinedSNPs["X1188368236"] == "C"] <- "1"
combinedSNPs["X1188368236"][combinedSNPs["X1188368236"] == "CT"] <- "2"
combinedSNPs["X1188368236"][combinedSNPs["X1188368236"] == "T"] <- "3"

PCAsheet <- left_join(combinedSNPs, samplemeta)

# Filter out rows that contain "none" or "Florida" in a specific column
# Creating a list of values to filter by
filter_values <- c(0,1,2,3,4,5,6,7,8,10,11)
filter_values <- c(9,12,13,14,15,16,17,18,19,20,21,22,23,25,26,27,28,29,30,31)
# Using filter() function from dplyr package

PCA_filtered <- PCAsheet %>% filter(PopIDforsorting %in% filter_values)

#STEP 3A: Run PCA ####
#PCA comp
#How many columns? 47 
ncol(PCA_filtered)


#Make sure there are NO NAs in the annotation columns! Or else they will be dropped

#filter out features in STEP 5
SNPs<- PCA_filtered

SNPs_NA<-na.omit(SNPs)
SNPs_NAana <- left_join(SNPs_NA, samplemeta)
# check if data is numeric
is.numeric(SNPs_NAana[c(2:35)])
# convert column 2 to 35 to double
SNPs_NAana[, 2:35] <- lapply(SNPs_NAana[, 2:35], as.double)
#filter 

write.csv(SNPs_NAana, "Data/SNPsNAanaNewIntroLA")
results <- prcomp(SNPs_NAana[c(2:35)])

summary(results)

#PCAs
screeplot(results, type = "l", npcs = 15, main = "Screeplot of the first 10 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
cumpro <- cumsum(results$sdev^2 / sum(results$sdev^2))
plot(cumpro[0:15], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 6, col="blue", lty=5)
abline(h = 0.88759, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC6"),
       col=c("blue"), lty=5, cex=0.6)


plot(results$x[,1],results$x[,2], xlab="PC1 (25%)", ylab = "PC2 (8%)", main = "PC1 / PC2 - plot")


library("factoextra")
fviz_pca_ind(results, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = SNPs_NAana$SampleArea, 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = FALSE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Population") +
  ggtitle("2D PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))

#Plotly PCA
#try to make a plot in plotly 
explained_variance_ratio <- summary(results)[["importance"]]['Proportion of Variance',]
explained_variance_ratio <- 100 * explained_variance_ratio
components <- results[["x"]]
components <- data.frame(components)
components <- cbind(components, SNPs_NAana$BayGroup )
components$PC3 <- -components$PC3
components$PC2 <- -components$PC2
tot_explained_variance_ratio <- summary(results)[["importance"]]['Proportion of Variance',]
tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio)
tit = 'Total Explained Variance = 39.61'

fig3D <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, color = ~SNPs_NAana$BayGroup ) %>%
  add_markers(size = 12)

fig <- components %>%
  plot_ly()  %>%
  add_trace(
    type = 'splom',
    dimensions = list(
      list(label=paste('PC 1 (',toString(round(explained_variance_ratio[1],1)),'%)',sep = ''), values=~PC1),
      list(label=paste('PC 2 (',toString(round(explained_variance_ratio[2],1)),'%)',sep = ''), values=~PC2),
      list(label=paste('PC 3 (',toString(round(explained_variance_ratio[3],1)),'%)',sep = ''), values=~PC3),
      list(label=paste('PC 4 (',toString(round(explained_variance_ratio[4],1)),'%)',sep = ''), values=~PC4)),
    color = ~SNPs_NAana$BayGroup) %>%
  style(diagonal = list(visible = FALSE)) %>%
  layout(
    legend=list(title=list(text='color')),
    hovermode='closest',
    dragmode= 'select',
    plot_bgcolor='rgba(240,240,240, 0.95)',
    xaxis=list(domain=NULL, showline=F, zeroline=F, gridcolor='#ffff', ticklen=4),
    yaxis=list(domain=NULL, showline=F, zeroline=F, gridcolor='#ffff', ticklen=4),
    xaxis2=axis,
    xaxis3=axis,
    xaxis4=axis,
    yaxis2=axis,
    yaxis3=axis,
    yaxis4=axis )

fig
fig3D

#STEP 7 PED File prep####
#MAP/PED file prep

#merged_data
# Apply gsub() function to columns 2:4
merged_data[, 2:42] <- apply(merged_data[, 2:42], 2, function(x) {
  gsub("(?<!.)\\b(\\w)\\b(?!\\.)", "\\1\\1", x, perl = TRUE)
})

PEDfile<-  left_join(merged_data, samplemeta)

write.csv(PEDfile, "Data/Apr11pedfile.csv")
# STEP 4: Structure Sheet Prep ####
# Use `str_length` from `stringr` to determine the number of characters in each cell
# If the number of characters is 2, keep only the first character
# If the number of characters is 1, keep the character

#Drop the second letter
drop_second_letter <- function(i) {
  # Create a new data frame to store the modified data
  data_mod <- i
  
  # Loop through each column of the data frame, starting from the second column
  for (col_name in names(data_mod)[2:ncol(data_mod)]) {
    # Use the `str_length` function from the `stringr` library to get the length of each cell in the column
    cell_lengths <- str_length(data_mod[, col_name])
    
    # Use the `ifelse` function to only drop the second letter if the cell has two letters
    data_mod[, col_name] <- ifelse(cell_lengths == 2, str_sub(data_mod[, col_name], 1, 1), data_mod[, col_name])
  }
  
  # Return the modified data frame
  return(data_mod)
}

#Drop the first letter
drop_first_letter <- function(i) {
  # Create a new data frame to store the modified data
  data_mod <- i
  
  # Loop through each column of the data frame, starting from the second column
  for (col_name in names(data_mod)[2:ncol(data_mod)]) {
    # Use the `str_length` function from the `stringr` library to get the length of each cell in the column
    cell_lengths <- str_length(data_mod[, col_name])
    
    # Use the `ifelse` function to only drop the second letter if the cell has two letters
    data_mod[, col_name] <- ifelse(cell_lengths == 2, str_sub(data_mod[, col_name], 2, 2), data_mod[, col_name])
  }
  
  # Return the modified data frame
  return(data_mod)
}
#Use functions
SNP_1 <- drop_second_letter(merged_data)
SNP_2 <- drop_first_letter(merged_data)

#Bind together 
StructureSheet_Date <- bind_rows(SNP_1, SNP_2)

#merge with metadata
StructureSheet_WMeta <- left_join(StructureSheet_Date, samplemeta)

#write csv
write.csv(StructureSheet_WMeta, "Data/Structure_Apr12023.csv")

#NOTES from datalab meeting ####
#PCA- wide data, matrix 
#could go wide to long, convert, the go back wide 
#column names are a problem, have to understand how R would interpret the column ID

library(googlesheets4)
library(stringr)
SNPcode <- read_sheet("https://docs.google.com/spreadsheets/d/15kRDVvlwXJhp4x39YMiufcGrKmJ6u-bi/edit?usp=sharing&ouid=115204603573220813867&rtpof=true&sd=true")

SNPcode <- read_excel("/Users/Taylor1/Downloads/Copy of Aaeg-pop-SNPs_iPLEX_KDT_Positions.xlsx")


#for loop 
for (i in 1:nrow(SNPcode)){
  
 mysnp<- SNPcode$SNP_ID[i] 
 #replace any punctuation with nothing
 # any punctuation [[:punct:]]
 mysnp2 <- str_replace_all(mysnp, "[[:punct:]]", "")
 #puts pieces togetehr with no sepration
 mysnpname <- paste0("X", mysnp2)
 
 wildbase = SNPcode$`Ref(WT) allele`[i]
 mutantbase = SNPcode$`var (Mutant) allele`[i]
 
 #col_index<- colnames(SNP_1)==mysnpname ###dplyr option may be better 
 SNP_1 = SNP_1 %>% mutate(!!mysnpname:=case_when() )
 
 
}


#STEP 5: Call Rate Evaluation ####

# data is wide, make it long 
MergedAna<- left_join(merged_data, samplemeta)
long_alleleA <- pivot_longer(MergedAna, cols = -c(1,43:47), names_to = "SNP", values_to = "Genotype")

#PCA long allel
long_allele <- pivot_longer(PCA_filtered, cols = -c(1,43:47), names_to = "SNP", values_to = "Genotype")
long_allele <- long_allele %>% select(-c(3,4))

#calculate SNP Call Rate Per Sample
Sample_CallRate <- long_allele %>% 
  group_by(SAMPLE_NAME) %>%
  summarize(n_NA = sum(is.na(Genotype))) %>% 
  mutate(Call = (1-(n_NA/42))*100) 

#get sample N
TotalSample <- unique(long_allele$SAMPLE_NAME)  %>% 
  length()

#calculate Call Rate Per SNP
SNP_CallRate <- long_allele %>% 
  group_by(SNP) %>%
  summarize(n_NA = sum(is.na(Genotype))) %>% 
  mutate(Call = (1-(n_NA/836))*100) 

# lookup table questions
#check its a data frame
class(merged_data)

#cols have names, need this 
# base r way and dplyr way to do it

SNP_CallRate %>% 
  filter(Call > 90) %>%
  select(SNP)
#need a vector
#Select cols
SNindex <- SNP_CallRate$SNP[SNP_CallRate$Call>90]
filteredcol <- select(merged_data, SAMPLE_NAME, all_of(SNindex))

#update sheet 
PCA_filtered <- select(PCA_filtered, SAMPLE_NAME, all_of(SNindex))

#filter rows
SAMindex <- Sample_CallRate$SAMPLE_NAME[SNP_CallRate$Call>97]
#filtering rows here, not searching by name. Instead search by some logical test, like is it within samindex
filteredSAM <- filter(filteredcol, SAMPLE_NAME %in% SAMindex)

#new sheet 
PCA_filtered <- filter(PCA_filtered, SAMPLE_NAME %in% SAMindex)

# STEP 6: Allele Frequency Tables ####

#read in recoding sheet 
SNPcode <- read.csv("Data/VGSCSNPCode.csv") %>% 
  select(1,6:9,13) %>% 
  pivot_longer(cols = 2:5, 
               names_to = "Phenotype", 
               values_to = "Genotype")

#read in sample metadata
samplemeta <- read.csv("Data/IPlexMetaData.csv")

#Make VGSC SNPs list for search
VGSC_SNPs <- c('vgsc_3315931548','vgsc_3315931672', 'vgsc_3315983611',
         'vgsc_3315983763','vgsc_3315999297' ,'vgsc_3316014588',
         'vgsc_3316080722','X3315939224')

#filter by VGSC list
Allele_VGSC <- long_alleleA %>% 
  filter(SNP %in% VGSC_SNPs) %>%
  count(Genotype, SNP) 

# Replace NAs with No Call
Allele_VGSC$Genotype  <- replace_na(Allele_VGSC$Genotype, "No Call")

#Merge VGSC Key and filtered IPlex Data 
VGSC_table <- merge(Allele_VGSC, SNPcode)

#Merge metadata info
VGSC_table2 <- merge(VGSC_table, samplemeta)

#Pivot Wider
VGSCWide <- VGSC_table %>% select(3:5) %>%
  pivot_wider(names_from = Musca,
                           values_from= n)
  
#Pivot Wider
VGSCWide <- VGSC_table2 %>% select(3:5,6,9) %>%
  pivot_wider(names_from = Musca,
              values_from= n)


#Population grouping 
long_alleleMeta <- left_join(long_alleleA, samplemeta)
MetaSNP <- left_join(long_alleleMeta, SNPcode) %>% 
  filter(SNP %in% VGSC_SNPs)


Counts <- MetaSNP %>% 
  group_by(SampleArea) %>%
  count(Genotype, SNP, Phenotype, Musca)

widecounts <- Counts %>% pivot_wider(names_from = Musca,
                                                 values_from= n)


together <-left_join(merged_data, samplemeta)

write.csv(together)
