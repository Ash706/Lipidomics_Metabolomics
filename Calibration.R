#.....................................................
#
# Script for calibration of normalized MZmine data
#
# Author: 				    Thomas F. Dyrlund
# Modification date: 	2016-01-11
#
# Changelog:    
#               1.1 - Seperated the normalization and calibration script
#               1.0 - Initial version (Ashfaq Ali and Tommi Suvitaival)
#
#.....................................................


# Data set information
csvSeperator 		    = ","
csvSampleData       = "Z:\\_Data\\LIP1\\0034_Epos\\Processed_Aug_2016\\EPOS_0034_areas_normalized.csv"
csvCalibrationData  = "Z:\\_Data\\LIP1\\0034_Epos\\Processed_Aug_2016\\Calibration - Normalized2016-08-16.csv"
csvLibrary          = "Z:\\_Data\\LIP1\\0036_SevaMeal\\Processed\\CSV Files\\MS Library - 2016-02-18.csv"
csvOutput           = "Z:\\_Data\\LIP1\\0034_Epos\\Processed_Aug_2016\\Nomalized_calibrated_data.csv"

# Specify the calibrants and their RT boundaries for the unknown features
opts <- list()
opts$unknown.calibration$boundaries <- c( 0, 6,8.999999999, 21.0 )
opts$unknown.calibration$names.calibrants <- c( "LPC(16:0)", "PC(16:0e/18:1(9Z))", "TG(17:0/17:0/17:0)" )

# Specify the stock values of the calibrants
opts$stocks <- c( "TG(17:0/17:0/17:0)"=103.2/100, "LPC(16:0)"=99.1/100, "CE(18:2)"=103.6/100, "PC(16:0e/18:1(9Z))"=100/100 )

# Specify the text which is used to identify sample columns in the sample data
opts$identifier.observation.columns <- "Peak area"

opts$calibration.model.with.intercept <- FALSE

opts$calibrant.concentration.range.used <- c( 0, Inf ) # c( 101, Inf )

#............................................................................
#             You don't need to modify anything below this line
#............................................................................

# Load the MZmine data sets
Data_normalized = read.csv(csvSampleData, sep=csvSeperator, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
calibration_normalized = read.csv(csvCalibrationData, sep=csvSeperator, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)

# Replace zeros with NAs
Data_normalized[Data_normalized==0] <- NA

# Remove empty columns (MZmine will add an empty column at the end of the CSV file)
Data_normalized <- Data_normalized[ , which( !apply( X=is.na( Data_normalized ), MAR=2, FUN=all ) ) ]
calibration_normalized <- calibration_normalized[ , which( !apply( X=is.na( calibration_normalized ), MAR=2, FUN=all ) ) ]

# Import the MS library with the normalizer annotation information
#Library_data = read.csv(csvLibrary, sep=csvSeperator, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)

# Merge the data from the MS library with the MZmine data
# Keep only the "Compound Name" and "Normaliser" column from the MS Library
#Combined_Data <- merge(Library_data[,c("ID", "Name", "Calibrant")], Data_normalized, by.x="ID", by.y ="ID", all.y = TRUE)
Combined_Data <-Data_normalized
# Sanity check: Did the matching based on ID lead to matchin names?
print( all( Combined_Data$"Name.x"==Combined_Data$"Name.y", na.rm=TRUE ) )

# Remove the Name.y column and rename Name.x to Name
Combined_Data = Combined_Data[,!names(Combined_Data) %in% c("Name.y")]
names(Combined_Data)[names(Combined_Data)=="Name.x"] <- "Name"


#............................................................
# If needed, assign a calibrant to compounds based on RT
#............................................................

# Subset the features that have no assigned normaliser or which are unknown compounds
subset1 <- subset(Combined_Data, is.na( Combined_Data$Calibrant ) | Combined_Data$Calibrant=="" )

# Assign the RT-based normaliser
subset1$Calibrant <- cut( subset1$"row retention time", opts$unknown.calibration$boundaries, labels=opts$unknown.calibration$names.calibrants )

# Subset the features that do have an assigned normaliser
subset2 <- subset(Combined_Data, !is.na( Combined_Data$Calibrant ) & Combined_Data$Calibrant!="" )

# Re-combine the two subsetted data sets
Combined_Data <- rbind(subset1, subset2)

# Sort the features based on the calibrant
Combined_Data <- Combined_Data[ order( as.character( Combined_Data$Calibrant ) ), ]


#..................................................................
# Extract the concentration information for each calibrant sample
#..................................................................

# Set the row names equal to the calibrant names
rownames( calibration_normalized ) <- calibration_normalized$Name
#rownames( calibration_normalized ) <- calibration_normalized$"Compound Name"

# Extract only the sample columns
idx.sample.columns <- grep( pattern=opts$identifier.observation.columns, x=colnames( calibration_normalized ) ) # FIXME: Change to sample names.
calibration_normalized <- calibration_normalized[, idx.sample.columns ]

# Transpose the data frame
cal_Normalized_t <- as.data.frame( t( calibration_normalized ) )

# Extract the unique part of the name, where the calibrant concentration of the sample is defined
tmp <- unlist( stringr::str_extract( string=rownames( cal_Normalized_t ), pattern="CALIBCURVE.*ngml" ) )

# Remove the non-unique beginning and end of the text and convert the number into an integer
tmp <- unlist( stringr::str_replace( string=tmp, pattern="CALIBCURVE-", replacement="" ) )
tmp <- stringr::str_replace( string=tmp, pattern="ngml", replacement="" )
cal_Normalized_t$concentrations <- as.numeric( tmp )

# Optional
# Keep only calibration samples above a certain concentration
indecies = which( cal_Normalized_t$concentrations > 100 )
indecies = which( cal_Normalized_t$concentrations > opts$calibrant.concentration.range.used[ 1 ] & cal_Normalized_t$concentrations < opts$calibrant.concentration.range.used[ 2 ] )
cal_Normalized_t = cal_Normalized_t[ 6:11, ]

# Multiply concentrations with the stocks
concentrations.by.stocks <- cal_Normalized_t$concentrations %o% opts$stocks
rownames( concentrations.by.stocks ) <- rownames( cal_Normalized_t )

#..................................................................
# Perform the calibration of the data based on the calibrants
#..................................................................

# Initialize the result data frame
Calibrated_Data = Combined_Data
# idx.sample.columns.normalized.2 <- grep( pattern="Peak area", x=colnames( Combined_Data ) )
names.samples <- colnames( Combined_Data )[ grep( pattern="Peak area", x=colnames( Combined_Data ) ) ]

# Find the batch information from the date in the sample name. In this case, it is the third item when split by "_".
batch.of.sample.runs.by.date <- as.factor( unlist( lapply( X=stringr::str_split( string=names.samples, pattern="_" ), FUN=function( x ) { x[[ 3 ]] } ) ) )
batch.of.calibration.runs.by.date <- as.factor( unlist( lapply( X=stringr::str_split( string=rownames( cal_Normalized_t ), pattern="_" ), FUN=function( x ) { x[[ 3 ]] } ) ) )

batch.of.sample.runs.by.date <- as.factor(rep(1, length(names.samples)))
batch.of.calibration.runs.by.date <- as.factor(rep(1, nrow(cal_Normalized_t)))

# Replace all sample values with NA
Calibrated_Data[ , names.samples ] <- NA

# Initialize a list for the calibration models
calibration.models <- vector( mode="list", length=nlevels( batch.of.sample.runs.by.date ) )
names( calibration.models ) <- levels( batch.of.sample.runs.by.date )

tmp = NULL
calibrated.calibrants <- array( dim=dim( cal_Normalized_t ), dimnames=dimnames( cal_Normalized_t ) )

for ( j in 1:nlevels( batch.of.sample.runs.by.date ) ) { # Go through the batches.

  # Find the sample runs in batch 'j'.
  idx.sample.runs.batch.j <- which( batch.of.sample.runs.by.date==levels( batch.of.calibration.runs.by.date )[ j ] )
  # Find the calibration runs in batch 'j'.
  # Note: Using the batch levels from sample runs to make sure that the sample and calibration runs are from the same batch.
  idx.calibration.runs.batch.j <- which( batch.of.calibration.runs.by.date==levels( batch.of.sample.runs.by.date )[ j ] )
  
  calibration.models[[ j ]] <- vector( mode="list", length=ncol( concentrations.by.stocks ) )
  names( calibration.models[[ j ]] ) <- colnames( concentrations.by.stocks )
  
  for ( i in 1:ncol( concentrations.by.stocks ) ) { # Go through all calibrants (i.e. columns of the matrix)
    
    # Get the name of the calibrant
    name.calibrant.i <- colnames( concentrations.by.stocks )[ i ]
    
    # Find the peaks that are assigned to be calibrated with calibrant 'i'.
    idx.peaks.calibrant.i <- which( Combined_Data$Calibrant==name.calibrant.i )
    
    # Fit a inverse-weighted linear model from the observed values to the expected values in the calibration samples and regress the observations in the sample runs with the calibration model.
    if ( opts$calibration.model.with.intercept ) {
      # Fit
      calibration.models[[ j ]][[ i ]] <- lm( formula=concentrations.by.stocks[ idx.calibration.runs.batch.j, name.calibrant.i ] ~ cal_Normalized_t[ idx.calibration.runs.batch.j, name.calibrant.i ], weight=1/cal_Normalized_t[ idx.calibration.runs.batch.j, name.calibrant.i ] )
      # Regress
      Calibrated_Data[ idx.peaks.calibrant.i, names.samples[ idx.sample.runs.batch.j ] ] <- Combined_Data[ idx.peaks.calibrant.i, names.samples[ idx.sample.runs.batch.j ] ] * calibration.models[[ j ]][[ i ]]$coefficients[ 2 ] + calibration.models[[ j ]][[ i ]]$coefficients[ 1 ]
      # Sanity check: Calibrate the calibrant itself in the calibration samples.
      calibrated.calibrants[ idx.calibration.runs.batch.j, name.calibrant.i ] <- cal_Normalized_t[ idx.calibration.runs.batch.j, name.calibrant.i ] * calibration.models[[ j ]][[ i ]]$coefficients[ 2 ] + calibration.models[[ j ]][[ i ]]$coefficients[ 1 ]
    } else {
      # Fit
      calibration.models[[ j ]][[ i ]] <- lm( formula=concentrations.by.stocks[ idx.calibration.runs.batch.j, name.calibrant.i ] ~ cal_Normalized_t[ idx.calibration.runs.batch.j, name.calibrant.i ] - 1, weight=1/cal_Normalized_t[ idx.calibration.runs.batch.j, name.calibrant.i ] )
      # Regress
      Calibrated_Data[ idx.peaks.calibrant.i, names.samples[ idx.sample.runs.batch.j ] ] <- Combined_Data[ idx.peaks.calibrant.i, names.samples[ idx.sample.runs.batch.j ] ] * calibration.models[[ j ]][[ i ]]$coefficients
      # Sanity check: Calibrate the calibrant itself in the calibration samples.
      calibrated.calibrants[ idx.calibration.runs.batch.j, name.calibrant.i ] <- cal_Normalized_t[ idx.calibration.runs.batch.j, name.calibrant.i ] * calibration.models[[ j ]][[ i ]]$coefficients
    }
  }
}


# Sanity check: Calibrated calibration samples
matplot( x=calibrated.calibrants, type="b", xlab="Run index", ylab="Concentration (ng/ml)" )
legend( x="top", legend=colnames( calibrated.calibrants ), col=1:ncol( calibrated.calibrants ), lty=1:ncol( calibrated.calibrants ), pch=paste( 1:ncol( calibrated.calibrants ) ) )


batch.index.for.plot <- batch.of.calibration.runs.by.date
levels( batch.index.for.plot ) <- 1:nlevels( batch.index.for.plot )
batch.index.for.plot <- as.numeric( batch.index.for.plot )

for ( i in 1:length( calibration.models[[ 1 ]] ) ) { # Go through all calibrants.
  
  name.calibrant.i <- colnames( concentrations.by.stocks )[ i ]
  
  plotTitle = name.calibrant.i
  plot( x=cal_Normalized_t[ , name.calibrant.i ], y=concentrations.by.stocks[ , name.calibrant.i ], col=batch.index.for.plot, main=plotTitle, xlab="Normalized Ratios", ylab="Concentration (ng/ml)" )
  
  for ( j in 1:length( calibration.models ) ) { # Go through all batches.
    
    if ( opts$calibration.model.with.intercept ) {
      abline( coef=calibration.models[[ j ]][[ i ]]$coefficients, col=j, lty=j )
    } else {
      abline( a=0, b=calibration.models[[ j ]][[ i ]]$coefficients, col=j, lty=j )
    }
    
  }
  
  legend( x="topleft", legend=c( "Batch, r^2", paste( "r^2=", round( x=unlist( lapply( X=calibration.models, FUN=function( x ) { summary( x[[ i ]] )$r.squared } ) ), digits=3 ) ) ), col=c( NA, 1:length( calibration.models ) ), lty=c( NA, 1:length( calibration.models ) ))
  
}

# Sanity check: Plot the peak-wise average calibrated intensity versus retention time.

plot( x=Calibrated_Data$"row retention time", y=rowMeans( Calibrated_Data[ , names.samples ], na.rm=TRUE ), xlab="Retention time (min)", ylab="Concentration (ng/ml)" )

# Sanity check: Plot the peak-wise average calibrated intensity after log-transformation versus retention time.
plot( x=Calibrated_Data$"row retention time", y=rowMeans( log2( Calibrated_Data[ , names.samples ] ), na.rm=TRUE ), xlab="Retention time (min)", ylab="Log2(Concentration)" )

# Santity check: Verify that the number of NA values in the input and output data is equal
all(is.na(Combined_Data) == is.na(Calibrated_Data))


#............................................................
# Save the data and perform a clean up
#............................................................


# Order the calibrated data by row ID
Calibrated_Data = Calibrated_Data[order(Calibrated_Data$"row ID"), ]

# Write the output CSV file
write.csv( Calibrated_Data, file=csvOutput, row.names=FALSE, na="" )

# Clean up
rm(csvSeperator, csvSampleData, csvLibrary, csvOutput, opts)
rm(Combined_Data, subset1, subset2, cal_Normalized_t, Calibrated_Data)
rm(tmp, i, calibration_normalized, concentrations.by.stocks, Data_normalized)
rm(Library_data, calibration.models, csvCalibrationData)
rm(idx.peaks.calibrant.i, idx.sample.columns, names.samples, name.calibrant.i)