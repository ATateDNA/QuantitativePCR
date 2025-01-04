#Anderson Tate
#Genomic Variation Laboratory
#12/17/2024
#Data Cleaning, Merging and Visualizing Quantification Amplification Results and Cq Values CSV files from the Bio-Rad CFX96 Maestro Software.
#----
  
library(tidyverse)
library(dplyr)
library(tidyr)
library(conflicted)
library(readxl)
library(writexl)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("summarize", "dplyr")

# Import data
data <- read_csv("/Users/andersontate/Documents/UCDavis/GVL/eDNA_USGS_VernalPools/LabNotebook/Spadefoot qPCR/USGS_SPHA_11/Anderson_2025-01-03_SPHA11_DM -  Quantification Amplification Results_FAM.csv") # Import the Quantification Amplification Results CSV file.

# Tidy data to long form from CFX96
datatib <- tibble(data) #Create a tibble

datatidy <- datatib %>% pivot_longer(cols = -Cycle,# Long form for tidy data. Be careful if it is a partial plate. You will still need to identify the wells. Don't specify wells if it is a full plate. 
                             names_to = "Well",
                             values_to = "RFU")

# Combining sample name data from plate map with qPCR output files by well
# Load sample names from the spreadsheet if using Excel otherwise move down to loading CSV files
# Replace 'samples.xlsx' with your actual file path and sheet name
#Excelsampledata <- read_excel("samples.xlsx", col_names = FALSE)

# Extract sample names into a vector when using Excel
#Excelsamplenames <- Excelsampledata[[1]]  # Assuming sample names are in the first column

# Load the sample names from the CSV file
platemap <- read_csv("/Users/andersontate/Documents/UCDavis/GVL/eDNA_USGS_VernalPools/LabNotebook/Spadefoot qPCR/USGS_SPHA_11/qPCR_PlateMap_SPHA11.csv")

# Load the Cq values BioRad Cq output CSV file
Cq <- read_csv("/Users/andersontate/Documents/UCDavis/GVL/eDNA_USGS_VernalPools/LabNotebook/Spadefoot qPCR/USGS_SPHA_11/Anderson_2025-01-03_SPHA11_DM -  Quantification Cq Results.csv")

#Subset Cq values to only wells and Cq values
Cq2 <- Cq[, c("Well", "Cq")]
#Remove leading zeroes in well names otherwise it won't join later to merged_data.
Cq2$Well <- gsub("^(\\w)(0+)(\\d+)$", "\\1\\3", Cq2$Well)

# Extract sample names into a vector
samplenames <- platemap[[2]]# Assuming sample names are in the second column. Taken from plates maps in BOX folder for each qPCR plate. Sample names are exported from Excel to CSV as above.
sitenames <- platemap[[3]]

# Create a 96-well plate layout
rows <- LETTERS[1:8]   # Rows A to H
cols <- 1:12           # Columns 1 to 12

# Well names for well groups
wellnames <- as.vector(outer(rows, cols, paste0))  # Generate "A1", "A2", ..., "H12"

# Split wells into groups of three organized by columns on the plate
wellgroups <- split(wellnames, ceiling(seq_along(wellnames) / 3))

plate <- expand.grid(Row = rows, Col = cols)  # Create all combinations of rows and columns
plate$Well <- paste0(plate$Row, plate$Col)    # Combine rows and columns into well names

# Group wells in sets of three by column
plate$Sample <- rep(samplenames, each = 3, length.out = nrow(plate))
plate$Site <- rep(sitenames, each = 3, length.out = nrow(plate))

# Sort by column order (column-major order), this may not be necessary
plate <- plate[order(plate$Col, plate$Row), ]

#Merge sample names and qPCR output
# Join datasets by the 'Well' column to add sample names
datamerged <- left_join(datatidy, plate, by = "Well")
# Join datasets by the 'Well' column to add Cq values
datamerged <- left_join(Cq2, datamerged, by = "Well")

# Perform a logical test if using a threshold Cq for detection status, Cq values < 40 become 1, > 40 & "NaN" become 0
datamerged$Cq_Less_Than_40 <- ifelse(datamerged$Cq == "NaN", 0, ifelse(datamerged$Cq < 40, 1, 0))

# Create a combined label for the legend using Well and Cq rounded to two decimals for visualization later
datamerged <- datamerged %>%
  mutate(
    Cq = round(Cq, 2),
    Well_Cq = paste0(Well, " (Cq: ", Cq, ")"))

#subset the data to only have three wells of the sample replicates. Change the wellgroups number for a certain set of three wells. Organized by column with three replicates in each. See wellgroups in the Environment for the corresponding numbers.
subset <- subset(datamerged, Well %in% wellgroups[[3]]) 

# Basic graph plotting amplification. Can use the subset here or the tidy format for more specific wells. Don't specify any cols if you want the entire dataset.
ggplot(data = subset, aes(x = Cycle, y = RFU, group = Well_Cq)) + 
  geom_line(aes(y = RFU, color = Well_Cq)) + 
  labs(title = subset$Sample, x = "Cycle Number", 
       y = "RFUs", color = "Replicate") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5))  # Center the title

# Summarize the data by sample
# Remove duplicate rows for each well. This removes the cycle records down to only the first cycle for each sample in each well.
dataunique <- datamerged %>%
  distinct(Well, .keep_all = TRUE)

# Calculate the sum of Cq_Less_Than_40 for each sample
samplesums <- dataunique %>%
  group_by(Sample, Site) %>%
  summarize(Sum_Cq_Less_Than_40 = sum(Cq_Less_Than_40))

# Merge the calculated sums back into data_unique
dataunique <- dataunique %>%
  left_join(samplesums, by = "Sample")

# Summarize the data by groups of 3 wells, this may not be necessary
#summarytable <- dataunique %>%
  #group_by(Group) %>%
  #summarize(
    #Samples = paste(unique(Sample), collapse = ", "),  # Combine sample names
    #Sum_Cq_Less_Than_40 = sum(Cq_Less_Than_40))
#Create a directory to write files into
#dir.create("OutputTables", recursive = TRUE)
# Write file
write.csv(samplesums, "/Users/andersontate/Documents/UCDavis/GVL/eDNA_USGS_VernalPools/Coding/RWorkDir/OutputTables/SampleCqSumsSPHA11.csv", row.names = FALSE)

    
