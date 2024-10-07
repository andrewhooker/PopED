#  Examples
#  find.largest.index("sfg","bpop")
#  find.largest.index("sfg","b")
#  find.largest.index("sfg","bocc",mat=T,mat.row=T)
#  find.largest.index("sfg","x")
#  find.largest.index("sfg","a")
# 
# txt <- capture.output(eval(parse(text="sfg")))
# txt <- grep(paste("^[^\\#]*bpop","\\[",sep=""),txt,value=T)
# txt
# ind <- gsub(paste("^[^\\#]*","bpop","\\[","(\\d+)\\].*",sep=""),"\\1",txt)
# max(as.numeric(ind))

# # Sample variables
# lab <- "bpop"
# txt <- "some text bpop[ 42 ] more text bpop[ 56 ] and bpop[ 78 ]"

find.largest.index <- function (func.str="sfg",lab="bpop",mat=F,mat.row=T) {
  
  if(is.function(func.str)){
    #txt <- capture.output(func.str)
    txt <- deparse(func.str)
  } else {
    #txt <- capture.output(eval(parse(text=func.str)))
    txt <- deparse(eval(parse(text=func.str)))
  }
  
  txt <- grep(paste("^[^\\#]*",lab,"\\[",sep=""),txt,value=T)
  if(length(txt)==0) return(0)
  if(length(txt)!=0 && !mat){
    
    #ind <- gsub(paste("^[^\\#]*",lab,"\\[\\s*(\\d+)\\s*\\].*",sep=""),"\\1",txt) # only last match per string
    
    # Construct the pattern
    pattern <- paste("^[^\\#]*", lab, "\\[\\s*(\\d+)\\s*\\].*", sep="")
    
    # Use gregexpr to find all matches
    matches <- gregexpr(paste(lab, "\\[\\s*(\\d+)\\s*\\]", sep=""), txt)
    
    # Extract the matches
    extracted <- regmatches(txt, matches)
    
    # Extract the numbers from the matches
    numbers <- sapply(extracted, function(x) sub(paste("^[^\\#]*", lab, "\\[\\s*(\\d+)\\s*\\].*", sep=""), "\\1", x))
    
    # turn a list of lists into a vector
    ind <- max(as.numeric(unlist(numbers)))
    if(is.null(ind)) ind <- 0
    return(ind)
  }  
  if(length(txt)!=0 && mat){
    #ind <- gsub(paste("^[^\\#]*",lab,"\\[\\s*(\\d+)\\s*,.*?\\].*$",sep=""),"\\1",txt) # only last match per string
    
    # Construct the pattern
    pattern <- paste(lab, "\\[\\s*(\\d+)?\\s*,\\s*(\\d+)?\\s*\\]", sep="")
    
    # Find all matches
    matches <- gregexpr(pattern, txt, perl=TRUE)
    
    # Extract the matches
    extracted <- regmatches(txt, matches)
    
    # Extract the numbers from each match
    numbers_rows <- sapply(extracted, function(x) sub(paste("^[^\\#]*", lab, "\\[\\s*(\\d+)?\\s*,\\s*(\\d+)?\\s*\\].*$", sep=""), "\\1", x))
    numbers_cols <- sapply(extracted, function(x) sub(paste("^[^\\#]*", lab, "\\[\\s*(\\d+)?\\s*,\\s*(\\d+)?\\s*\\].*$", sep=""), "\\2", x))
    
    # turn a list of lists into a vector
    ind_rows <- as.numeric(unlist(numbers_rows))
    ind_cols <- as.numeric(unlist(numbers_cols))
    
    if(is.null(ind_rows)) ind_rows <- 0
    if(is.null(ind_cols)) ind_cols <- 0
    
    if(mat.row) return(max(ind_rows))
    if(!mat.row) return(max(ind_cols))
    #return(c(rows=max(ind_rows),cols=max(ind_cols)))
  }  
  #if(length(txt)!=0 && mat && !mat.row)  ind <- gsub(paste("^[^\\#]*",lab,"\\[.*?,\\s*(\\d+)\\s*\\].*",sep=""),"\\1",txt)
  
  #max(as.numeric(ind))
  
}



