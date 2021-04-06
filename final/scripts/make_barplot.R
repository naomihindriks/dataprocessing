#!/usr/bin/env Rscript
library(ggplot2)

source("scripts/functions.R")

args = commandArgs(trailingOnly=TRUE)
# args[1], args[2], args[3], args[4] are required, the rest is optional

# args[1], args[2] and args[3] are expected to be the filenames of the feature-table (a tab seperated file),
#   taxonomy (another tab seperated file) and SraRunTable (a comma seperated file)

# args[4] is expected to be one of the taxonomic levels: "Kingdom", "Phylum", "Class", "Order", "Family", "Genus" or "Species"

# args[5] is the percentage at which a feature needs to be present in at least one of the samples at a specific time and treatment
#   default is 0.5

# args[6] is the name to which to save the plot that is made with this script
#   default is barplot.jpg

# args[7] and args[8] are the width and height (in pixels) of the plot that is made with this script
#   default is 600 for both.
#   if only args[7] the width and height will both be set to this specified number

# args[9] is a possible filter to set e.g. Kingdom=Bacteria the data will now be filtered so that
#   only features of the kingdom bacteria will be shown in the plot

taxonomic_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

if (length(args) < 4 | !(capitalize(args[4]) %in% taxonomic_levels)) {
  writeLines(paste("Please give at least 4 argument:\nthe first should be a path to the features file\nthe second should be a bath to the taxonomy file\nthe third should be a path to the SraRunTable file\nand the fourth should specify which taxonomic level to use, one of the following:", paste(taxonomic_levels, collapse = ", ")))
  quit(save = "no", status = 1)
} 
if(length(args) < 5) {
  args[5] = 0.5
}
if(length(args) < 6) {
  args[6] = "barplot.jpg"
}
if (length(args) < 7) {
  args[7] = 600
  args[8] = 600
}
if (length(args) < 8) {
  args[8] = args[7]
}

args[4] <- capitalize(args[4])

# load information from feature table (how many times each feature is found in each sample)
feature_table <- read.csv(
  args[1], 
  sep = "\t", 
  skip = 1, 
  header = TRUE
)

# load information from taxonomy (for every feature in feature table taxonomic assignment)
taxonomy <- read.csv(
  args[2], 
  sep = "\t", 
  header = TRUE
)

# load information metadata from samples, which individual does each sample belong to, which group and which time
samples <- read.csv(
  args[3],
  sep = ",", 
  header = TRUE
)

# Split the taxonomic assignment in taxonomy on ; to get assignment per taxonomic level
split_taxa <- strsplit(taxonomy$Taxon, ";")

# Create empty dataframe to fill with the taxonomic assignments for each taxonomic level 
taxon_dataframe <- data.frame()

# There is a taxonomic assignment for every feature in  split_taxa, loop through every feature
for (feature_taxa in split_taxa) {
  i<-0
  prev_name <- NULL
  # For every taxonomic level in feature_taxa remove the first letter and underscore. 
  # If there is no taxonomic assignment on this specific level set to NaN
  for (taxon_level in feature_taxa) {
    i <- i+1
    split_name <- str_trim(unlist(strsplit(taxon_level, "_")), side = "both")
    if(length(split_name) < 3) {
      feature_taxa[i] <- NaN
    } else {
      if(startsWith(split_name[3], "[") & endsWith(split_name[3], "]")) {
        split_name[3] <- substr(split_name[3], 2, nchar(split_name[3]) - 1)
      }
      if (split_name[1] == "s") {
        feature_taxa[i] <- paste(prev_name, split_name[3])
      } 
      else {
        feature_taxa[i] <- split_name[3]
      }
    }
    prev_name = feature_taxa[i]
  }
  # If the length of feature_taxa less than 7 that means the lower taxonomic levels 
  # aren't filled in and will be set to NaN
  if (length(feature_taxa) < 7) {
    feature_taxa <- c(feature_taxa, rep(NaN, 7 - length(feature_taxa)))
  }
  taxon_dataframe <- rbind(taxon_dataframe, feature_taxa)
}

# Set the column names of the taxon dataframe to the different taxonomic levels
colnames(taxon_dataframe) <- taxonomic_levels

# Add the information in the taxon_dataframe to the taxonomy dataframe
taxonomy <- cbind(taxonomy, taxon_dataframe)

# Cleanup
rm(taxon_dataframe, feature_taxa, i, split_name, split_taxa, taxon_level, prev_name)

# Combine feature_table dataframe with taxonomy dataframe
merged_data = merge.data.frame(
  feature_table, taxonomy, 
  by = 1
)

# Clean up
rm(feature_table, taxonomy)

if(!is.na(args[9])) {
  if(!(args[9] == "")) {
    x <- capitalize(unlist(strsplit(args[9], "=")))
    merged_data <- merged_data[merged_data[x[1]] == x[2],]
  }
}

sample_columns = colnames(merged_data)[startsWith(colnames(merged_data), "SRR")]

merged_data <- merged_data[c(sample_columns, args[4])]

merged_data <- aggregate(formula(paste(". ~ ", args[4])), merged_data, FUN = sum)

sample_name_columns <- (2:(length(sample_columns)))

# Clean up
rm(sample_columns)

data <- reshape(
  merged_data, 
  direction = "long",
  varying = list(names(merged_data)[sample_name_columns]),
  v.names = "Value",
  timevar = "Sample",
  times = names(merged_data)[sample_name_columns]
)

# Clean up
rm(sample_name_columns, merged_data)

# Change NA values in time to 0
samples["Time"][is.na(samples["Time"])] <- 0

treatments = unique(samples["Treatment"])[unique(samples["Treatment"]) != "None"]

for (treatment in treatments) {
  add_data <- samples[samples$Treatment == "None", ]
  add_data["Treatment"] <- treatment
  samples <- rbind(samples, add_data)
}

# Clean up
rm(treatment, treatments, add_data)

samples <- samples[samples$Treatment != "None",]

data <- merge(data, samples[c("Run", "horse", "Treatment", "Time")], by.x = "Sample", by.y = "Run")

# Clean up
rm(samples)

data <- data[c("Sample", args[4], "Value", "horse", "Treatment", "Time")]

# Create dataframe that has a row for each possible combination of time, horse and treatment
ntimes <- length(unique(data$Time))
nhorses <- length(unique(data$horse))
ntreatments <- length(unique(data$Treatment))

total <- ntimes * nhorses * ntreatments

x <- data.frame(Time = rep(unique(data$Time), total / ntimes))
y <- data.frame(horse = rep(unique(data$horse), each = ntimes, length.out = total))
z <- data.frame(Treatment = rep(unique(data$Treatment), each = ntimes * nhorses, length.out = total))
subset_dataframe <- cbind(x, y, z)

# Cleanup
rm(total, x, y, z, ntreatments, nhorses, ntimes)

# For each row in the created dataframe above run the get_taxonomies_to_keep function 
# and get filter to get all unique taxonomies.
# This will give a vector (called x) that holds all the operational taxonomic unit found 
# at or greater than the threshold given with args[5] at each time point for each culture condition and horse 
x <- apply(X = subset_dataframe, MARGIN = 1, FUN = get_taxonomies_to_keep)
x <- unique(unlist(x))

rm(subset_dataframe)

# Filter data to only contain the OTU's in x
data <- data[data[,args[4]] %in% x,]

horse_labels <- paste("horse", unique(data$horse))
names(horse_labels) <- 1:length(unique(data$horse))

plot <- ggplot(data, aes(fill=get(args[4]), y=Value, x=factor(Time))) +
  geom_bar(position="fill", stat="identity") +
  facet_grid(Treatment ~ horse, labeller = labeller(horse = horse_labels)) +
  labs(
    title = paste("Relative abundance at the ", args[4], " level"),
    y="Relative abundance", 
    x = "Time (hrs)",
    fill = args[4]
  )

if (length(x) <= 12) {
  plot <- plot + scale_fill_brewer(palette = "Paired")
}

# Save plot in jpeg file
jpeg(args[6], width = as.numeric(args[7]), height = as.numeric(args[8]))
plot
graphics.off()