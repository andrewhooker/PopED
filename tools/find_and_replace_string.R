## The following example demonstrates
## how a serach and replace string task 
## can be peformed with R across several files

## Create two text files with content 
filenames <- c( tempfile(), tempfile() )
for( f in  filenames ){
  cat("We wish you a Merry Christmas!\n\nBest regards\n", file=f)
}

## Replace Merry Christmas with Happy New Year
for( f in filenames ){
  
  x <- readLines(f)
  y <- gsub( "Merry Christmas", "Happy New Year", x )
  cat(y, file=f, sep="\n")
  
}

## Review output
for( f in filenames ){ 
  cat(readLines(f), sep="\n")
}