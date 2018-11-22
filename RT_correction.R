##
## Script for adjusting the retention times for a linear shift.
##
## Idea:
## -Run targeted identification of internal standard compounds from the internal standards runs.
## -Provide updated retention times found from the internal standards samples for 'RT.stds.new'.
## -Provide path to the reference library
## -Provide path to RT-corrected library files (full library, internal standards library, calibration library) to be written
##
## Notes
## -Creating the calibration library is unfinished in this script. It does not include all the respective standards.
## -TODO: Find standard compounds and updated retention times from the csv file of internal standards runs exported from MZmine.
##
## Version history:
## -08.12.2015: Added to the repository.
##
## Author:
## Tommi Suvitaival
## tsvv@steno.dk
## 11.11.2015
##

opts <- list()
opts$digits <- 2 # Number of digits in the result values.

opts$path.library <- file.path( "Z:", "_Data", "LIP1", "0001_Alspac", "MS_Library" ) ## CHANGE
opts$fname.library.reference <- "MS_Library-2015-11-05_from_TFRD--identified_only.txt" ## CHANGE
opts$fname.library.new <- "MS_Library-2015-11-05_from_TFRD--RT_shift_corrected.txt" ## CHANGE
opts$fname.library.new.standards.only <- "MS_Library-2015-11-05_from_TFRD--RT_shift_corrected--internal_standards_only.txt" ## CHANGE
opts$fname.library.new.calibration.only <- "MS_Library-2015-11-05_from_TFRD--RT_shift_corrected--calibration_only.txt" ## CHANGE
opts$fname.RT.correction.model <- "RT_correction-linear_model.RData"

# Specify the normalisers and their RT boundaries for the unknown peaks. This has to be specified here, since the calibrant CE(18:2) does not have an assigned normalizer but it is determined based on the retention time. -3.12.15
opts$unknown.normalization$boundaries <- c( 0, 6, 8.999999999, 21 )
opts$unknown.normalization$normalisers <- c( "05.LPC(17:0)", "09b.PC(16:0/d30/18:1)", "07.TG(16:0/16:0/16:0)-13C3" )

# Reference retention times from "MS Library - 2015-11-05".
# These are retention times from the TUMME data set.

# RT.reference <- c( 7.85, 7.16, 7.56, 7.73, 4.25, 6.31, 10.10, 5.67, 7.24, 7.24 )
# names( RT.reference ) <- c( "01.PE(17:0/17:0)", "02.SM(d18:1/17:0)", "03.Cer(d18:1/17:0)", "04.PC(17:0/17:0)", "05.LPC(17:0)", "06.PC(14:0/d13)", "07.TG(16:0/16:0/16:0)-13C3", "08.TG(8:0/8:0/8:0)-13C3", "09a.PC(16:0/d31/18:1)", "09b.PC(16:0/d30/18:1)" )

# Average retention times from the 0001_ALSPAC standard samples ("not OK" samples not included)

RT.stds.new <- c( 7.92, 7.23, 7.63, 7.80, 4.30, 6.38, 10.20, 5.71, 7.31, 7.31 ) ## CHANGE
names( RT.stds.new ) <- c( "01.PE(17:0/17:0)", "02.SM(d18:1/17:0)", "03.Cer(d18:1/17:0)", "04.PC(17:0/17:0)", "05.LPC(17:0)", "06.PC(14:0/d13)", "07.TG(16:0/16:0/16:0)-13C3", "08.TG(8:0/8:0/8:0)-13C3", "09a.PC(16:0/d31/18:1)", "09b.PC(16:0/d30/18:1)" )

# Load the reference library.

library.reference <- read.table( file = file.path( opts$path.library, opts$fname.library.reference ),
                                 header = TRUE,
                                 sep = "\t",
                                 dec = ".",
                                 stringsAsFactors = FALSE
)

# Read the column names separately so that they can be saved unchanged in the new file.

colnames.library <- read.table( file = file.path( opts$path.library, opts$fname.library.reference ),
                                header = FALSE,
                                sep = "\t",
                                dec = ".",
                                stringsAsFactors = FALSE,
                                nrows = 1
)
colnames.library[ which( is.na( colnames.library ) ) ] <- ""

# Find the corresponding standard in the reference library.

idx.stds.reference <- match( x=names( RT.stds.new ), table=library.reference$Name )

difference.stds.new.v.reference <- RT.stds.new - library.reference$RT[ idx.stds.reference ]


## Plot difference in the RT as a function of the reference RT.

opts.plot <- list()
opts.plot$col <- 1:length( difference.stds.new.v.reference )
opts.plot$pch <- 16

plot( x=library.reference$RT[ idx.stds.reference ], y=difference.stds.new.v.reference, col=opts.plot$col, pch=opts.plot$pch )
legend( x="topleft", legend=names( difference.stds.new.v.reference ), col=opts.plot$col, pch=opts.plot$pch )


## Initialize the new library.

library.new <- library.reference

# Move the reference RTs into a separate column.

library.new$RT.reference <- library.new$RT
library.new$RT[] <- NA

# Copy the new standard RT values to the new library.

library.new$RT[ idx.stds.reference ] <- RT.stds.new

## Fit a linear model for the difference in the RT.

# Fit the linear model from reference RT to new RT using the observed standard RTs.
lm.RT <- lm( formula = RT ~ RT.reference, data = library.new, subset = idx.stds.reference )

save( lm.RT, file=file.path( opts$path.library, opts$fname.RT.correction.model ) )

# Predict the new RT of the other items given the reference RT and round the predicted values with the specified accuracy.
library.new$RT[ -idx.stds.reference ] <- round( x=predict.lm( object=lm.RT, newdata=library.new )[ -idx.stds.reference ], digits=opts$digits )

library.new.txt <- rbind( c( colnames.library, "RT Reference" ), as.matrix( library.new ) )

# Remove the empty column from the array.
idx.non.missing.columns.library.new <- which( apply( X=!is.na( library.new.txt ), MAR=2, FUN=all ) )
library.new.txt <- library.new.txt[ , idx.non.missing.columns.library.new ]

# Copy standards to a separate standards-only file.
# Drop the ID column (first in the array).
library.new.txt.standards.only <- library.new.txt[ c( 1, idx.stds.reference+1 ), -1 ]
# library.new.txt.standards.only[ 1, which( library.new.txt.standards.only[ 1, ]=="Theoretical m/z" ) ] <- "m/z"
# library.new.txt.standards.only[ 1, which( library.new.txt.standards.only[ 1, ]=="RT" ) ] <- "Retention time (min)"
# library.new.txt.standards.only[ 1, which( library.new.txt.standards.only[ 1, ]=="Name" ) ] <- "Identity"


write.table( x=library.new.txt, file=file.path( opts$path.library, opts$fname.library.new ), col.names=FALSE, dec=".", na="", row.names=FALSE, sep=";" )

write.table( x=library.new.txt.standards.only, file=file.path( opts$path.library, opts$fname.library.new.standards.only ), col.names=FALSE, dec=".", na="", row.names=FALSE, sep=";" )


## Calibration set

# Specify the names of the compounds that will be used for calibration.
# names.calibration <- c( "PC(16:0e/18:1(9Z))", "LPC(16:0)", "TG(17:0/17:0/17:0)" )
names.calibration <- c( "PC(16:0e/18:1(9Z))", "LPC(16:0)", "TG(17:0/17:0/17:0)", "CE(18:2)" ) # Calibrant for CEs added 3.12.15

# Find the calibration compounds in the library.
idx.calibration.reference <- match( x=names.calibration, table=library.new$Name )

# Find the respective normalisers for the calibration compounds.
idx.calibration.reference <- c( idx.calibration.reference, match( x=library.new[ idx.calibration.reference, "Normaliser" ], table=library.new$Name ) )

# Check if there were some calibrants that did not have a normaliser specified (CE(18:2) does not have).
tmp <- which( is.na( idx.calibration.reference ) )

if ( length( tmp ) > 0 ) {
  for ( i in 1:length( tmp ) ) {
    # Find the retention time of the calibrant 'i' that does not have normaliser specified.
    RT.calibrant.i <- library.new[ idx.calibration.reference[ tmp-length( names.calibration ) ], "RT" ]
    # Find the RT-based normaliser for the calibrant.
    idx.normaliser.i <- which( RT.calibrant.i < opts$unknown.normalization$boundaries )[ 1 ] - 1
    # Find the entry of the RT-based normaliser in the library.
    idx.calibration.reference[ tmp ] <- match( x=opts$unknown.normalization$normalisers[ idx.normaliser.i ], table=library.new$Name )
    # FIXME: Add the RT-based normaliser to the normaliser columm.

  }
}

# Remove duplicates.
idx.calibration.reference <- unique( idx.calibration.reference )

# Write the calibration library.

library.new.txt.calibration.only <- library.new.txt[ c( 1, idx.calibration.reference+1 ), -1 ]

write.table( x=library.new.txt.calibration.only, file=file.path( opts$path.library, opts$fname.library.new.calibration.only ), col.names=FALSE, dec=".", na="", row.names=FALSE, sep=";" )
