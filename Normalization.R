#.....................................................
#
# Script for normalization of MZmine data
#
# Author: 				    Thomas F. Dyrlund
# Modification date: 	2016-01-08
#
# Changelog:    
#               1.1 - Seperated the normalization and calibration script
#               1.0 - Initial version (Ashfaq Ali and Tommi Suvitaival)
#
#.....................................................


# Data set information
csvSeperator 		    = ","
csvSampleData       = "Z:\\_Data\\LIP1\\0034_Epos\\Processed_Aug_2016\\EPOS_0034_areas_norm.csv"
csvStandards        = "Z:\\_Data\\LIP1\\0034_Epos\\Processed_Aug_2016\\standards_areas_2.csv"
csvLibrary          = "Z:\\_Data\\LIP1\\0036_SevaMeal\\Processed\\CSV Files\\MS Library - 2016-02-18.csv"
csvOutput           = "Z:\\_Data\\LIP1\\0034_Epos\\Processed_Aug_2016\\EPOS_0034_areas_normalized.csv"

# Specify the normalisers and their RT boundaries for the unknown features
opts <- list()
opts$unknown.normalization$boundaries <- c( 0, 6, 8.999999999, 21 )
opts$unknown.normalization$normalisers <- c( "05.LPC(17:0)", "09b.PC(16:0/d30/18:1)", "07.TG(16:0/16:0/16:0)-13C3" )

# Specify the text which is used to identify sample columns in the sample data
opts$identifier.observation.columns <- "Peak area"


#............................................................................
#             You don't need to modify anything below this line
#............................................................................

# Load the MZmine data sets
Data_filled_area = read.csv(csvSampleData, sep=csvSeperator, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
Data_standards_loaded = read.csv(csvStandards, sep=csvSeperator, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)

# Replace zeros with NAs
Data_filled_area[Data_filled_area==0] <- NA

# Remove empty columns (MZmine will add an empty column at the end of the CSV file)
Data_filled_area <- Data_filled_area[ , which( !apply( X=is.na( Data_filled_area ), MAR=2, FUN=all ) ) ]
Data_standards_loaded <- Data_standards_loaded[ , which( !apply( X=is.na( Data_standards_loaded ), MAR=2, FUN=all ) ) ]

# Import the MS library with the normalizer annotation information
#Library_data = read.csv(csvLibrary, sep=csvSeperator, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)

# Merge the data from the MS library with the MZmine data
# Keep only the "Compound Name" and "Normaliser" column from the MS Library
#Combined_Data <- merge(Library_data[,c("ID", "Name", "Normaliser")], Data_filled_area, by.x="ID", by.y ="ID", all.y = TRUE)
Combined_Data <-Data_filled_area
# Sanity check: Did the matching based on ID lead to matchin names?
print( all( Combined_Data$"Name.x"==Combined_Data$"Name.y", na.rm=TRUE ) )

# Remove the Name.y column and rename Name.x to Name
Combined_Data = Combined_Data[,!names(Combined_Data) %in% c("Name.y")]
names(Combined_Data)[names(Combined_Data)=="Name.x"] <- "Name"


#............................................................
# If needed, assign a normaliser to compounds based on RT
#............................................................

# Subset the features that have no assigned normaliser or which are unknown compounds
subset1 <- subset(Combined_Data, is.na( Combined_Data$Normaliser ) | Combined_Data$Normaliser=="" )

# Assign the RT-based normaliser
subset1$Normaliser <- cut( subset1$"row retention time", opts$unknown.normalization$boundaries, labels=opts$unknown.normalization$normalisers )

# Subset the features that do have an assigned normaliser
subset2 <- subset(Combined_Data, !is.na( Combined_Data$Normaliser ) & Combined_Data$Normaliser!="" )

# Re-combine the two subsetted data sets
Combined_Data <- rbind(subset1, subset2)

# Sort the features based on the normaliser
Combined_Data <- Combined_Data[ order( as.character( Combined_Data$Normaliser ) ), ]


#............................................................
# Match the internal standards data to the sample data
#............................................................

# Extract the names of the samples and sort them alphabetically
tmp <- colnames( Data_filled_area )[ grep( pattern=opts$identifier.observation.columns, x=colnames( Data_filled_area ) ) ]
names.samples <- sort(tmp)

# Sort the standard data based on the sample column name
standards <- Data_standards_loaded[ , c( colnames( Data_standards_loaded )[ 1:5 ], names.samples ) ]

# Set the row names equal to the standard compound names
rownames(standards) <- standards$Name


#............................................................
# Perform the compound normalization based on the normalizer
#............................................................

# Split the data set based on the normaliser
Data_x <- split(Combined_Data, Combined_Data$Normaliser)

# Sanity check: There should be no more than 6 standards used for normalization
print(length(Data_x) <= 6)

# Initialize the list for the normalized data set.
norm <- list()

# Loop through all the normalizers
for (i in 1:length(Data_x) ) {
    
    # Normalize the values against those of the standards
    norm[[ i ]] <- t( t( Data_x[[ i ]][ , names.samples ] )/ unlist( standards[ names( Data_x )[ i ], names.samples ] ) )
    
    # Recombine the sample columns
    norm[[ i ]] <- cbind( Data_x[[ i ]][ , -which( colnames( Data_x[[ i ]] ) %in% names.samples ) ], as.data.frame( norm[[ i ]] ) )
    
}

# Re-combine the normalized data sets
Normalized <- NULL
for ( i in 1:length( norm ) ) {
    Normalized <- rbind( Normalized, norm[[ i ]] )
}


#............................................................
# Save the data and perform a clean up
#............................................................

# Order the normalized data by row ID
Normalized = Normalized[order(Normalized$"row ID"), ]

# Write the output CSV file
write.csv(Normalized, file=csvOutput)

# Clean up
rm(csvSeperator, csvSampleData, csvStandards, csvLibrary, csvOutput, opts)
rm(Data_filled_area, Data_standards_loaded, Library_data, Combined_Data, subset1, subset2)
rm(tmp, standards, names.samples, Data_x, i,norm,Normalized)

