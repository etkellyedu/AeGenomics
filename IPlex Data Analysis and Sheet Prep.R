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

# STEP 1: Turn genotype excel sheet into a csv ####
# notes on manual manipulations of sheets
# First plate did not have 1534
# one submission did not have 1534 primers


# Define the path to .xlsx files and the folder for the new csv sheets
pathin <- "Data/RawIplex/"
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


# STEP 2: Make genotype data file by combining all genotype CSVs ####
# specify the directory path
# old drive directory_path <- "Google Drive/Shared drives/Attardo Lab/Projects/Insecticide Resistance in Aedes aegypti/Data/iPlex Results /Iplex Genotype CSV/"
directory_path <- "Data/GenotypeCSV/"

#Note: our test has two plates. Create master sheets for each plate, then merge at the end
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

# STEP 2A: Cleanup ####

#drop well column
merged_data <- merged_data[, -c(2,3,26,27)]
#some Ts get turned to TRUE (logical), change all to character
merged_data[merged_data == "TRUE"] <- "T"

#troubleshooting: IGNORE
mastersheet <- read_excel("/Users/Taylor1/Downloads/Master Data Sheet AEG+OC.xlsx")
colnames(mastersheet)[3] <- 'SAMPLE_NAME'
df1 <- select(mastersheet, 2,3)
df2 <- merged_data

# Join the two datasets on the matching column, "SAMPLE_NAME"
merged_df <- merge(df1, df2, by = "SAMPLE_NAME", all.x = TRUE)

# Find the rows with no match in dataset 1
no_match_df1 <- df1[!(df1$SAMPLE_NAME %in% merged_df$SAMPLE_NAME), ]

# Find the rows with no match in dataset 2
no_match_df2 <- df2[!(df2$SAMPLE_NAME %in% merged_df$SAMPLE_NAME), ]

# STEP 3: PCA Sheet Prep ####
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



# STEP 4: Structure Sheet Prep ####
# Use `str_length` from `stringr` to determine the number of characters in each cell
# If the number of characters is 2, keep only the first character
# If the number of characters is 1, keep the character

#Drop the second letter
drop_second_letter <- function(i) {
  # Create a new data frame to store the modified data
  data_mod <- i
  
  # Loop through each column of the data frame, starting from the third column
  for (col_name in names(data_mod)[3:ncol(data_mod)]) {
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
  
  # Loop through each column of the data frame, starting from the third column
  for (col_name in names(data_mod)[3:ncol(data_mod)]) {
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


#STEP 5: Allele Frequency Tables ####
# data is wide, make it long 
long_allele <- pivot_longer(merged_data, cols = -1, names_to = "SNP", values_to = "Genotype")



