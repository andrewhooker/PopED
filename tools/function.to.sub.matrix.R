
function.to.sub.matrix <- function(file=NULL,file.str=NULL,mat){
    ##from="../warfarin/eval/"
    ##to="test.R"
    ##browser()
    ##match.call()[[1]]

    tmp <- NULL
    
    if(!is.null(file)) tmp <- readLines(file,warn=F)
    if(!is.null(file.str)) tmp <- file.str


    if(is.null(tmp)) stop("no file or file string detected")

    file.str <- tmp
    
    for(i in 1:length(file.str)){
        fun.str <- file.str[i]
        for(tmp.name in mat){
            if(length(grep(paste(tmp.name,"\\(",sep=""),fun.str)!=0)){
                drop.str=""
                if(length(grep(paste("=.*",tmp.name,"\\(\\s*,",sep=""),fun.str))!=0) drop.str=",drop=F"
                if(length(grep(paste("=.*",tmp.name,"\\(.*?,\\s*\\)",sep=""),fun.str))!=0) drop.str=",drop=F"
                if(length(grep(paste("=.*",tmp.name,"\\(.*?\\:.*\\)",sep=""),fun.str))!=0) drop.str=",drop=F"
                if(length(grep(paste("<-[^\\n]*",tmp.name,"\\(.*?\\:.*\\)",sep=""),fun.str))!=0) drop.str=",drop=F"
                if(tmp.name=="\\W+epsi") drop.str=""
                
                ## find the right parenthesis to substitute
                ##fun.str <- "y(:,(i-1)*NumEPS+1:i*NumEPS)=tmp(:,1:NumEPS);" ## for testing
                                        #fun.str <- "y(:,i-1*NumEPS+1:i*NumEPS)=tmp(:,1:NumEPS);" ## for testing
                ##tmp.name <- "\\s*y" ## for testing
                mat.str <- gsub(paste(".*",tmp.name,"\\((.*)",sep=""),"\\1",fun.str)
                mat.arg.str.split <- unlist(strsplit(mat.str,")"))
                mat.arg.str.split ## for testing
                mat.arg.str <- ""
                mat.arg.str <- mat.arg.str.split[1]
                mat.arg.str.split <- mat.arg.str.split[-1]
                test.ok <- FALSE
                while(!test.ok){
                    if(num.char("\\(",mat.arg.str)!=(num.char("\\)",mat.arg.str))){
                        if(length(mat.arg.str.split)==0){
                            test.ok <- TRUE
                            cat("return arguments for not right")
                            next
                        }
                        mat.arg.str <- paste(mat.arg.str,mat.arg.str.split[1],sep=")")
                        mat.arg.str.split <- mat.arg.str.split[-1]
                    } else {
                        test.ok <- TRUE
                    }
                }
                ##cat(mat.arg.str,"\n")
                n.right.paren <- num.char("\\)",mat.arg.str)
                ##n.right.paren
                for(k in 0:n.right.paren){
                    if(k==0) regex.tmp <- paste("[^\\)]*)\\)",sep="")
                    if(k>0) regex.tmp <- paste("[^\\)]*\\)",regex.tmp,sep="") 
                }
                ##regex.tmp
                regex.tmp.1 <- paste("(",tmp.name,")\\((",regex.tmp,sep="")
                ##regex.tmp.1
                regex.tmp.2 <- paste("\\1[\\2",drop.str,"]",sep="")
                
                fun.str <- gsub(regex.tmp.1,regex.tmp.2,fun.str)
                ##fun.str
                ##fun.str <- gsub(paste("(",tmp.name,")\\(([^\\)]*)\\)",sep=""),paste("\\1[\\2",drop.str,"]",sep=""),fun.str)
                ##fun.str <- gsub(paste(tmp.name,"(",mat.arg.str,")",sep=""),
                ##                paste(tmp.name,"[",mat.arg.str,drop.str,"]",sep=""),
                ##                fun.str,fixed=T)
                
                ##after.mat <- gsub(paste(".*",tmp.name,"(.*)",sep=""),"\\1",fun.str)
                ## count number of parenthesis
                ##left.paren.loc <- gregexpr("\\(",after.mat)
                ##n.left <- length(left.paren.loc[[1]])
                ##if(left.paren.loc[[1]]==-1) n.left <- 0
                ##right.paren.loc <- gregexpr("\\)",after.mat)
                ##n.right <- length(right.paren.loc[[1]])
                ##if(right.paren.loc[[1]]==-1) n.right <- 0
                ##n.right - n.left
            }
        }
        file.str[i] <- fun.str
    }

    if(!is.null(file)){
        out.file <- gsub(".R$",".old.R",file)
        i <- 0
        while(file.exists(out.file)){
            i <- i+1
            out.file <- gsub(".R$",paste(".old.",i,".R",sep=""),file)
        }
        file.copy(file,out.file)
        writeLines(file.str, con=file)
    }
    return(file.str)
}
