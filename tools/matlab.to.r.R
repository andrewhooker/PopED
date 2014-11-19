num.char <- function(char.str,full.str){
    tmp <- gregexpr(char.str,full.str)
    n.tmp <- length(tmp[[1]])
    if(tmp[[1]][1]==-1) n.tmp <- 0
    return(n.tmp)
}



matlab.to.r <- function(from, function.input = F, convert.input=F){
    ##from="../warfarin/eval/"
    ##to="test.R"
    ##browser()
    ##match.call()[[1]]
    
    file.str <- readLines(from,warn=F)

    more <- FALSE
    more.ret.args <- FALSE
    all.args.str.open <- FALSE

    line.skip <- 0
    ret.args <- ""
    
    for(i in 1:length(file.str)){
        fun.str <- file.str[i]
        #if(i==329) browser()

        ##########################
        ## code specific
        ##########################
        fun.str <- gsub("function plotLineSearch(numRows,optSum,optsw,detmf_vector,it_vector,xt_vector,x_vector,a_vector,m,gni,it,Engine)",
                        "function []=plotLineSearch(numRows,optSum,optsw,detmf_vector,it_vector,xt_vector,x_vector,a_vector,m,gni,it,Engine)",fun.str, fixed=T)
        fun.str <- gsub("init_vec = [b_ind' e0];","init_vec = [b_ind', e0];",fun.str, fixed=T)
        fun.str <- gsub("designsin = cell(1,0);%temporary solution","designsin = cell(1,0); % temporary solution",fun.str, fixed=T)
        fun.str <- gsub("val = deriv(j).hx;","val = deriv(j)$hx;",fun.str, fixed=T)
        fun.str <- gsub("[lgao kita inversiona aopt]","[lgao, kita, inversiona, aopt]",fun.str, fixed=T)
        fun.str <- gsub("[lgxto kitxt inversionxt xtopt]","[lgxto, kitxt, inversionxt, xtopt]",fun.str, fixed=T)
        fun.str <- gsub("[lgvaro kitvar inversionvar varopt]","[lgvaro, kitvar, inversionvar, varopt]",fun.str, fixed=T)
        fun.str <- gsub("[params params_cvs]","[params, params_cvs]",fun.str, fixed=T)
        fun.str <- gsub("[params params_cv]","[params, params_cv]",fun.str, fixed=T)
        fun.str <- gsub("[params param_cvs]","[params, param_cvs]",fun.str, fixed=T)
        fun.str <- gsub("[strPath strFunctionName]","[strPath, strFunctionName]",fun.str, fixed=T)
        
        fun.str <- gsub("function blockexp(fn,globalStructure)","function []=blockexp(fn,globalStructure)",fun.str, fixed=T)
        fun.str <- gsub("function blockfinal(fn,fmf,dmf,groupsize,ni,xt,x,a,bpop,d,docc,sigma,m,globalStructure)",
                        "function []=blockfinal(fn,fmf,dmf,groupsize,ni,xt,x,a,bpop,d,docc,sigma,m,globalStructure)",
                        fun.str, fixed=T)
        fun.str <- gsub("function blockopt(fn,globalStructure)","function []=blockopt(fn,globalStructure)",fun.str, fixed=T)
        fun.str <- gsub("function blockother(fn,globalStructure)","function []=blockother(fn,globalStructure)",fun.str, fixed=T)
        fun.str <- gsub("function blockoptwrt(fn,optsw)","function []=blockoptwrt(fn,optsw)",fun.str, fixed=T)
        fun.str <- gsub("function Dtrace(fn,it,ni,xtopt,xopt,aopt,gxt,ga,dmf,diff,ixt,ia,itvector,dmfvector,globalStructure)",
                        "function []=Dtrace(fn,it,ni,xtopt,xopt,aopt,gxt,ga,dmf,diff,ixt,ia,itvector,dmfvector,globalStructure)",
                        fun.str, fixed=T)
        fun.str <- gsub("function writet(f,ni,xt)","function []=writet(f,ni,xt)",fun.str, fixed=T)
        fun.str <- gsub("function write_iterationfile(","function []=write_iterationfile(",fun.str, fixed=T)
        fun.str <- gsub("fprintf(fn,'function [%s] = %s()\n\n',strReturn,strFunctionName);",
                        "fprintf(fn,'%s <- \\f\\u\\n\\c\\t\\i\\o\\n\\(\\)\\{\n\n',strFunctionName)",fun.str, fixed=T)
        
        fun.str <- gsub("fprintf(fn,'%s','end');","fprintf(fn,'%s %s', strReturn,'}')",fun.str, fixed=T)

        ##fun.str <- "[amount,globalStructure]=calculate_diff(@diff_eq1,[g(5), 0],xt,globalStructure,params);"
        fun.str <- gsub("(calculate_diff\\([^,]+\\,)\\s*\\[(.*)\\]","\\1cbind\\(\\2\\)",fun.str)
        ##fun.str
        
        ## change or remove lines then skip to end
        ## this is code specific
        if(line.skip>0){
            fun.str <- paste("##",fun.str)
            line.skip <- line.skip-1
            file.str[i] <- fun.str
            next
        }
        if(length(grep("bpop = params(1:nbpop); d=params((1+nbpop):(nbpop+nd)); covd=params((1+nbpop+nd):(nbpop+nd+ncovd));",
                       fun.str,fixed=T))!=0){ 
            fun.str <- gsub("bpop = params(1:nbpop); d=params((1+nbpop):(nbpop+nd)); covd=params((1+nbpop+nd):(nbpop+nd+ncovd));",
                            "bpop = params[1:nbpop,drop=F] \n d=params[(1+nbpop):(nbpop+nd),drop=F] \n covd=params[(1+nbpop+nd):(nbpop+nd+ncovd),drop=F]\n",
                            fun.str, fixed=T)
            file.str[i] <- fun.str
            next
        }
        if(length(grep(" docc=params((1+nbpop+nd+ncovd):(nbpop+nd+ncovd+ndocc)); covdocc=params((1+nbpop+nd+ncovd+ndocc):(nbpop+nd+ncovd+ndocc+ncovdocc));",
                       fun.str,fixed=T))!=0){ 
            fun.str <- gsub(" docc=params((1+nbpop+nd+ncovd):(nbpop+nd+ncovd+ndocc)); covdocc=params((1+nbpop+nd+ncovd+ndocc):(nbpop+nd+ncovd+ndocc+ncovdocc));",
                            " docc=params[(1+nbpop+nd+ncovd):(nbpop+nd+ncovd+ndocc),drop=F] \n covdocc=params[(1+nbpop+nd+ncovd+ndocc):(nbpop+nd+ncovd+ndocc+ncovdocc)] \n",
                            fun.str, fixed=T)
            file.str[i] <- fun.str
            next
        }
        if(length(grep(" sigma=params((1+nbpop+nd+ncovd+ndocc+ncovdocc):(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma)); covsigma=params((1+nbpop+nd+ncovd+ndocc+ncovdocc+nsigma):(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma+ncovsigma));",
                       fun.str,fixed=T))!=0){ 
            fun.str <- gsub(" sigma=params((1+nbpop+nd+ncovd+ndocc+ncovdocc):(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma)); covsigma=params((1+nbpop+nd+ncovd+ndocc+ncovdocc+nsigma):(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma+ncovsigma));",
                            " sigma=params[(1+nbpop+nd+ncovd+ndocc+ncovdocc):(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma),drop=F] \n covsigma=params[(1+nbpop+nd+ncovd+ndocc+ncovdocc+nsigma):(nbpop+nd+ncovd+ndocc+ncovdocc+nsigma+ncovsigma)] \n",
                            fun.str, fixed=T)
            file.str[i] <- fun.str
            next
        }

         if(length(grep("t= [bpop(:,1); d(:,1); ones(length(globalStructure.covd),1); docc(:,1); ones(length(globalStructure.covdocc),1); ones(length(sigma),1)]; %",
                       fun.str,fixed=T))!=0){ 
            fun.str <- gsub("t= [bpop(:,1); d(:,1); ones(length(globalStructure.covd),1); docc(:,1); ones(length(globalStructure.covdocc),1); ones(length(sigma),1)]; %",
                            "t= matrix(c(bpop[,1,drop=F], d[,1,drop=F], matrix(1,length(globalStructure$covd),1), docc[,1,drop=F], matrix(1,length(globalStructure$covdocc),1), matrix(1,length(sigma),1)),ncol=1,byrow=T) #",
                            fun.str, fixed=T)
            file.str[i] <- fun.str
            next
        }
        if(length(grep("warning('off','MATLAB:nearlySingularMatrix')",fun.str,fixed=T))!=0 ||
           length(grep("warning('on','MATLAB:nearlySingularMatrix')",fun.str,fixed=T))!=0){
            fun.str <- paste("##",fun.str)
            ##line.skip <- 2
            file.str[i] <- fun.str
            next
        }
        if(length(grep("if (~strcmp(strfgModelFilePath,''))",fun.str,fixed=T))!=0 ||
           length(grep("if (~strcmp(strUserDistFilePath,''))",fun.str,fixed=T))!=0 ||
           length(grep("if (~strcmp(strEDPenaltyFilePath,''))",fun.str,fixed=T))!=0 ||
           length(grep("if (~strcmp(strAutoCorrelationFilePath,''))",fun.str,fixed=T))!=0 ||
           length(grep("if (~strcmp(strffModelFilePath,''))",fun.str,fixed=T))!=0 ||
           length(grep("if (~strcmp(strErrorModelFilePath,''))",fun.str,fixed=T))!=0 ||
           length(grep("if (~strcmp(strErrorModelFilePath,''))",fun.str,fixed=T))!=0){

            fun.str <- paste("##",fun.str)
            line.skip <- 2
            file.str[i] <- fun.str
            next
        }
        if(length(grep("eval(sprintf('globalStructure.subffPointers.ff_pointer%d = eval(''@%s'');',i,strffModelFilename))",
                       fun.str,fixed=T))!=0){
            fun.str <- "globalStructure$subffPointers[paste('ff_pointer',i,sep='')] = strffModelFilename" 
            file.str[i] <- fun.str
            next
        }
        if(length(grep("%Set the model file string path",fun.str,fixed=T))!=0){
            fun.str <- paste("##",fun.str,"\n","globalStructure$model_file = popedInput$ff_file",sep="")
            line.skip <- 6
            file.str[i] <- fun.str
            next
        }
        if(length(grep("rand('state', globalStructure.dSeed)",fun.str,fixed=T))!=0){
            fun.str <- "set.seed(globalStructure$dSeed)"
            line.skip <- 1
            file.str[i] <- fun.str
            next
        }
        if(length(grep("This function fixes a bug in FreeMat 4.0",fun.str,fixed=T))!=0){
            fun.str <- paste("##",fun.str)
            line.skip <- 7
            file.str[i] <- fun.str
            next
        }
        if(length(grep("the grid is defined by ls_step_size",fun.str,fixed=T))!=0){
            fun.str <- paste("##", fun.str,"\n itvector <- c() \n dmfvector <- c()")
            file.str[i] <- fun.str
            next
        }

        ## ########################
        ## END code specific
        ## ########################


        ## add return argument to function end in middle of file
        if(length(grep("^function\\s+",fun.str))!=0 && ret.args!=""){
            ## add return argument to end of function
            ret.args.2 <- gsub("^list\\(","",ret.args)
            ret.args.2 <- gsub("\\)$","",ret.args.2)
            ret.args.2 <- gsub("([^,]*),","\\1=\\1,",ret.args.2)
            ret.args.2 <- gsub(",([^,]*)$",",\\1=\\1",ret.args.2)
            ret.args <- gsub("^list\\((.*)\\)",paste("list(",ret.args.2,")",sep=""),ret.args)
            
            
            only.spaces <- grep("^[[:space:]]*$",file.str[1:(i-1)])
            only.comments <- grep("^\\s*\\#.*",file.str[1:(i-1)])
            only.scrap <- c(only.spaces,only.comments)
            if(length(only.scrap)!=0){
                foo <- 1:length(file.str[1:(i-1)])
                foo <- foo[-only.scrap]
                last.line.with.code <- foo[length(foo)]
            } else {
                last.line.with.code <- length(file.str[i-1])
            }
            
            ##if(file.str[length(file.str)]=="}"){
            if(file.str[last.line.with.code]=="}"){
                ##file.str[length(file.str)]=paste("return(",ret.args,")",sep="")
                file.str[last.line.with.code]=paste("return(",ret.args,") \n}",sep="")
                ##file.str[length(file.str)+1]="}"
                ##file.str[last.line.with.code+1]="}"
            } else {
                                        #cat("return argument needs to be added")
                file.str[i-1]=paste(file.str[i-1],"\n return(",ret.args,")\n}",sep="")
            }
            ret.args <- ""
        }

        ## get return arguments
        ##tmp <- grep("function\\s*\\[(.*)\\]\\s*=",fun.str,value=T)
        tmp <- grep("function\\s+.*=",fun.str,value=T)
        ##if(length(tmp)!=0) ret.args <- gsub(".*\\[(.*)\\].*","\\1",tmp)
        if(length(tmp)!=0){
            ret.args <- gsub("function(.*)=.*","\\1",tmp)
            ret.args <- gsub("\\[","",ret.args)
            ret.args <- gsub("\\]","",ret.args)
            if(gregexpr(",",ret.args)[[1]][1]!=-1) ret.args <- paste("list(",ret.args,")",sep="")
        }

        tmp <- grep("function\\s+.*\\.\\.\\.",fun.str,value=T)
        ##if(length(tmp)!=0) ret.args <- gsub(".*\\[(.*)\\].*","\\1",tmp)
        if(more.ret.args){
            if(length(grep("\\]\\s*=.*",fun.str,value=T))!=0) more.ret.args <- F
            ret.args.more <- gsub("(.*)=.*","\\1",fun.str)
            ret.args.more <- gsub("\\]","",ret.args.more)
            ret.args.more <- gsub("(.*)\\.\\.\\.","\\1",ret.args.more)
            ret.args <- paste(ret.args,ret.args.more,sep="\n")
            if(!more.ret.args) ret.args <- paste(ret.args,")",sep="")
            if(!more.ret.args) fun.str <- gsub(".*=\\s*(.*)(\\(.*\\))","\\1 <- function\\2\\{",fun.str)
        }
        if(length(tmp)!=0){
            more.ret.args <- T
            ret.args <- gsub("function(.*)\\.\\.\\.","\\1",tmp)
            ret.args <- gsub("\\[","",ret.args)
            ret.args <- gsub("\\]","",ret.args)
            if(gregexpr(",",ret.args)[[1]][1]!=-1) ret.args <- paste("list(",ret.args,sep="")
            fun.str <- "\n"
        }

        
        ## replace function definition
        fun.str <- gsub("function\\s+(.*)\\s*=\\s*(.*)(\\(.*\\))","\\2 <- function\\3\\{",fun.str)
        ##fun.str <- gsub("function\\s+(.*)\\s*=\\s*(.*)(\\(.*\\))","\\2 <- function\\3\\{",fun.str)

        #subCell <- length(grep("\\w*\\{[^}]*}",fun.str))!=0

        
        ## replace some common things
        fun.str <- gsub("^\\%","\\#",fun.str)
        fun.str <- gsub("\\.\\*","*",fun.str)
        fun.str <- gsub("\\.\\^","^",fun.str)
        fun.str <- gsub("\\.\\/","/",fun.str)
        fun.str <- gsub("(\\s+?)\\%","\\1\\#",fun.str)
        ##fun.str <- gsub("(fprintf\\([^\\)]*?\\s+?)\\#","\\1\\%",fun.str)
        if(length(grep("fprintf",fun.str,fixed=T))!=0){
            fun.str <- gsub("(\\s+?)\\#","\\1\\%",fun.str)
        }
        fun.str <- gsub("(\\s)end","\\1}",fun.str)
        fun.str <- gsub("imag\\(","Im(",fun.str)
        fun.str <- gsub("^end","}",fun.str)
        fun.str <- gsub("~","!",fun.str)
        fun.str <- gsub("(zeros\\(.*\\))'(.*)","t\\(\\1)\\2",fun.str)
        #fun.str <- gsub("zeros\\(","matrix\\(0,",fun.str)
        fun.str <- gsub("cell\\((.*)\\)'(.*)","t(cell(\\1))\\2",fun.str)
        fun.str <- gsub("=\\{}","=cell(0,0)",fun.str)
        ##fun.str <- gsub("cell\\(","matrix\\(0,",fun.str) ## more needed here
        fun.str <- gsub("pargen\\((.*)\\)'","t(pargen(\\1))",fun.str)
        fun.str <- gsub("correlateSamples\\((.*)\\)'","t(correlateSamples(\\1))",fun.str)
        fun.str <- gsub("ones\\(","matrix\\(1,",fun.str)
        fun.str <- gsub("eye\\(","diag\\(1,",fun.str)
        ##fun.str <- gsub("isempty\\((.*?)\\)","any(dim(\\1)==0)",fun.str)
        fun.str <- gsub("strcat\\((.*)\\)","paste\\(\\1,sep=''\\)",fun.str)
        fun.str <- gsub("strcat\\((.*)\\)","paste\\(\\1,sep=''\\)",fun.str)
        fun.str <- gsub("(\\W*)true(\\W*)","\\1TRUE\\2",fun.str)
        fun.str <- gsub("fopen\\(","file\\(",fun.str)
        fun.str <- gsub("fclose\\(","close\\(",fun.str)

        fun.str <- gsub("(calculate_diff\\()\\@([^,]+)\\,","\\1\\'\\2\\'\\,",fun.str)
        
        fun.str <- gsub("(\\W)disp\\(","\\1print\\(",fun.str)
        
        ## replace 'for' command
        ##fun.str <- gsub("for\\s*(.*)=\\s*(.*):(.*)","for\\(\\1 in \\2:\\3\\)\\{",fun.str)
        fun.str <- gsub("for\\s*(.*)=\\s*(.*):([^\\#]+)(.*)","for\\(\\1 in \\2:\\3\\)\\{\\4",fun.str)
        
        ##if(i==231) browser()
        
        ## replace 'while' command
        fun.str <- gsub("while\\s*(.*)","while\\1{",fun.str)

        ## replace 'if' command
        fun.str <- gsub("^(\\s*)if\\s*([^\\#]+)(.*)","\\1if(\\2){\\3",fun.str)
        ##fun.str <- gsub("(\\s)if\\s*(.*)","\\1if\\(\\2\\)\\{",fun.str)
        ##fun.str <- gsub("^if\\s*(.*)","if\\(\\1\\)\\{",fun.str)
        fun.str <- gsub("(\\s)else","\\1\\} else \\{",fun.str)
        fun.str <- gsub("^else","\\} else \\{",fun.str)
        
        ## replace the '=' operator with '<-', not needed
        #fun.str <- gsub("==","@@",fun.str)
        #fun.str <- gsub("!=","@_@",fun.str)
        #fun.str <- gsub("=","<-",fun.str)
        #fun.str <- gsub("@@","==",fun.str)
        #fun.str <- gsub("@_@","!=",fun.str)
        
        ## replace structures with lists
        fun.str <- gsub("([[:alpha:]]+)\\.([[:alpha:]]+?)","\\1\\$\\2",fun.str)
        fun.str <- gsub("\\$m'",".R'",fun.str)
        fun.str <- gsub("\\$exe'","\\.exe'",fun.str)
        ##fun.str <- gsub("_",".",fun.str)

        ## replace reshape and diag
        ##fun.str <- "ir=reshape(ir,ns,1)"
        fun.str <- gsub("reshape\\(","reshape_matlab\\(",fun.str)
        fun.str <- gsub("diag\\(","diag_matlab\\(",fun.str)
        ##fun.str

        fun.str <- gsub("(\\W)trace\\(","\\1trace_matrix\\(",fun.str)
        
        ## build matricies in R
        fun.str <- gsub("\\[\\]","matrix(0,0,0)",fun.str)
        matrixStart <- (length(grep("\\[",fun.str))!=0 && length(grep("\\[.*\\]\\s*=",fun.str))==0
                        && length(grep("error\\('[^']*\\[",fun.str))==0)
        if(more){ ## continuation of matrix
            more <- length(grep("]",fun.str))==0 # multiline matrix
            if(more){  ## get matrix values 
                mat.vals <- fun.str
            } else {
                mat.vals <- gsub("(.*)\\].*","\\1",fun.str) 
            }
            n.rows.new <- length(gregexpr(";",mat.vals)[[1]]) ## get number of new rows
            n.rows <- n.rows+n.rows.new
            mat.vals <- gsub(";",",",mat.vals)
            mat.vals <- gsub("(\\d+)\\s","\\1,",mat.vals)
            if(more){
                commas <- length(grep(",",mat.vals))!=0
                if(!commas) mat.vals <- gsub("^([^#\n]*)","\\1,",mat.vals)
                fun.str <- mat.vals
            } else {
                tmp <- paste(mat.vals,"\\1",sep="")
                fun.str <- gsub(".*(\\].*)",tmp,fun.str)
            }
            tmp.1 <- "]"
            tmp <- paste("\\),nrow=",n.rows,",byrow=T\\)",sep="")
            if(transpose){
                tmp.1 <- "]'"
                tmp <- paste("\\),nrow=",n.rows,",byrow=T\\)\\)",sep="")
            }
            fun.str <- gsub(tmp.1,tmp,fun.str)
        }
        if(matrixStart){ # start of a matrix
            #if(i==62) browser()
            more <- length(grep("]",fun.str))==0 # multiline matrix
            transpose <- length(grep("]'",fun.str))!=0 # transpose matrix
            if(more){  ## get matrix values 
                mat.vals <- gsub(".*\\[(.*)","\\1",fun.str)
            } else {
                mat.vals <- gsub(".*\\[(.*)\\].*","\\1",fun.str) 
            }
            n.rows <- length(gregexpr(";",mat.vals)[[1]]) ## get number of rows in mat.vals
            mat.vals <- gsub(";",",",mat.vals) 
            mat.vals <- gsub("(\\d+)\\s","\\1,",mat.vals)
            commas <- length(grep(",",mat.vals))!=0
            spaces <- length(grep("\\w+\\s+\\w+",mat.vals))!=0
            if(!commas && spaces) mat.vals <- gsub("(\\w+)(\\s+)","\\1,\\2",mat.vals)
            if(more){
                commas <- length(grep(",",mat.vals))!=0
                if(!commas) mat.vals <- gsub("^([^#\n]*)","\\1,",mat.vals)
                tmp <- paste("matrix(c(",mat.vals,sep="")
                if(transpose) tmp <- paste("t(matrix(c(",mat.vals,sep="")
                fun.str <- gsub("\\[.*",tmp,fun.str)
            } else {
                tmp <- paste("matrix(c(",mat.vals,"\\1",sep="")
                if(transpose) tmp <- paste("t(matrix(c(",mat.vals,"\\1",sep="")
                fun.str <- gsub("\\[.*(\\].*)",tmp,fun.str)
            }
            tmp.1 <- "]"
            tmp <- paste("\\),nrow=",n.rows,",byrow=T\\)",sep="")
            if(transpose){
                tmp.1 <- "]'"
                tmp <- paste("\\),nrow=",n.rows,",byrow=T\\)\\)",sep="")
            }
            fun.str <- gsub(tmp.1,tmp,fun.str)
        }

        ## build return arguments for called functions  
        functionStart <- (length(grep("\\[.*\\].*=",fun.str))!=0)
        functionStart.comment <- (length(grep("\\#.*\\[.*\\].*=",fun.str))!=0)
        cont.str <- (length(grep("\\.\\.\\.",fun.str))!=0)
        if(functionStart && !functionStart.comment){
            ## get return arguments
            ret.str <- gsub(".*?\\[(.*?)].*","\\1",fun.str)
            #ret.str.splt <- strsplit(ret.str,",")[[1]]
            ret.str.split <- unlist(strsplit(ret.str,","))
            n.tmp <- length(ret.str.split)
            all.args.str <- ""
            for(j in 1:n.tmp){
                if(length(ret.str.split)==0) next
                arg.str <- ret.str.split[1]
                ret.str.split <- ret.str.split[-1]
                test.ok <- FALSE
                while(!test.ok){
                    if(num.char("\\(",arg.str)!=num.char("\\)",arg.str)){
                        if(length(ret.str.split)==0){
                            test.ok <- TRUE
                            cat("return arguments for not right")
                            next
                        }
                        arg.str <- paste(arg.str,ret.str.split[1],sep=",")
                        ret.str.split <- ret.str.split[-1]
                    } else {
                        test.ok <- TRUE
                        ##cat(arg.str,"\n")
                    }
                }
                arg.str <- paste(arg.str," <- returnArgs[[",j,"]]",sep="")
                all.args.str <- paste(all.args.str,arg.str,sep="\n")
            }
            
            ## add a list called returnArgs as the return argument
            tmp <- gsub("(.*?)\\[.*\\].*?=(.*)","\\2",fun.str)
            ##cat(tmp, " needs to be defined\n")
            fun.str <- gsub("(.*?)\\[.*\\].*?=(.*)","\\1 returnArgs <- \\2",fun.str)
            
            ## make retun arguments == to values in list
            if(!cont.str) {
                fun.str <- paste(fun.str,all.args.str,sep=" ")
            } else {
                all.args.str.open <- TRUE
            }
        }
        if(all.args.str.open && !cont.str){
            fun.str <- paste(fun.str,all.args.str,sep=" ")
            all.args.str.open <- FALSE
        }
        
        ## must be after matrix replacement
        fun.str <- gsub(";","",fun.str)
        fun.str <- gsub("mfilename\\(\\)","match.call\\(\\)\\[\\[1\\]\\]",fun.str)
        fun.str <- gsub("GetCalculationEngine\\(.*\\)","list\\(Type=1,Version=version\\$version\\.string\\)",fun.str)
        ##fun.str <- gsub("isfield\\(\\s*(.*)\\s*,\\s*'(.*)'\\s*\\)","!is.null\\(\\1\\$\\2\\)",fun.str)
        fun.str <- gsub("\\(\\s*:\\s*,","\\(,",fun.str)
        fun.str <- gsub("\\(\\s*:\\s*\\)","",fun.str)
        fun.str <- gsub("\\(([^,]*),:\\)","(\\1,)",fun.str)
        ##fun.str <- gsub("size\\((.*?),(.*?)\\)","dim(\\1)[\\2]",fun.str)
        fun.str <- gsub("strcmp\\((.*),(.*)\\)","\\1==\\2",fun.str)
        fun.str <- gsub("char\\(","as.character\\(",fun.str)
        fun.str <- gsub("(\\w+)\\{([^}]*)}","\\1[[\\2]]",fun.str)  ## replace sub-cell calls
        fun.str <- gsub("false","FALSE",fun.str)
        fun.str <- gsub("don\\'\\'t","do not",fun.str)
        fun.str <- gsub("error\\((.*)\\)","stop(sprintf(\\1))",fun.str)
        fun.str <- gsub("error (\\'.*\\')","stop(sprintf(\\1))",fun.str)
        fun.str <- gsub("\\.\\.\\.","",fun.str)

        ##-------------------------------------------
        ## code specific
        ##-------------------------------------------
        fun.str <- gsub("sum\\(100\\*clock\\)","as.integer(Sys.time())",fun.str)
        fun.str <- gsub("=ni'\\*","=t(ni)\\*",fun.str)
        fun.str <- gsub("=\\s+(x\\(.*\\))'","= t(\\1)",fun.str)
        fun.str <- gsub("=\\s+(a\\(.*\\))'","= t(\\1)",fun.str)
        fun.str <- gsub("(xt\\(i,1\\:ni\\(i\\)\\))'","t(\\1)",fun.str)
        fun.str <- gsub("(model_switch\\(i,1\\:ni\\(i\\)\\))'","t(\\1)",fun.str)
        fun.str <- gsub("(b_ind)'","t(\\1)",fun.str)
        fun.str <- gsub("(m1_tmp)'","t(\\1)",fun.str)
        fun.str <- gsub("\\+(f1)'","\\+t(\\1)",fun.str)
        fun.str <- gsub("\\*(l)'","\\*t(\\1)",fun.str)
        fun.str <- gsub("\\*(locc_tmp)'","\\*t(\\1)",fun.str)
        fun.str <- gsub("\\*(lh)'","\\*t(\\1)",fun.str)
        fun.str <- gsub("\\*(h)'","\\*t(\\1)",fun.str)
        fun.str <- gsub("(\\s+)covd'","\\1t(covd)",fun.str)
        fun.str <- gsub("(\\s+)covdocc'","\\1t(covdocc)",fun.str)
        fun.str <- gsub("(\\s+)covsigma'","\\1t(covsigma)",fun.str)
        fun.str <- gsub("\\*(l\\(,\\w+\\))'","\\*t(\\1)",fun.str)
        fun.str <- gsub("\\*(lh\\(,\\(i-1\\)\\*NumSigma\\+1\\:i\\*NumSigma\\))'","\\*t(\\1)",fun.str)
        fun.str <- gsub("\\*(lh\\(,\\(n-1\\)\\*NumSigma\\+1\\:n\\*NumSigma\\))'","\\*t(\\1)",fun.str)
        fun.str <- gsub("\\*(lh\\(,\\(m-1\\)\\*NumSigma\\+1\\:m\\*NumSigma\\))'","\\*t(\\1)",fun.str)
        fun.str <- gsub("\\*(locc\\[\\[\\w+\\]\\]\\(,\\w+\\))'","\\*t(\\1)",fun.str)
        fun.str <- gsub("\\*(h\\(,\\w+\\))'","\\*t(\\1)",fun.str)
        ##fun.str <- gsub("\\*()'","\\*t(\\1)",fun.str)
        fun.str <- gsub("\\*(tmp_lh)'","\\*t(\\1)",fun.str)
        fun.str <- gsub("\\*(tmp_lh_\\w)'","\\*t(\\1)",fun.str)
        fun.str <- gsub("b_ind\\%ind","bind\\#ind",fun.str)
        fun.str <- gsub(",true,",",TRUE,",fun.str)
        
        
        ## add source to file if using the fileparts command
        fp.found <- (length(grep("fileparts",fun.str))!=0)
        if(fp.found){
            extra.str <- gsub(".*fileparts(.*?)\\n.*","source\\1",fun.str)
            fun.str <- paste(extra.str,fun.str,sep="\n")
        }
        ## change anonoumous functions to first argument of do.call  
        if(length(grep("globalStructure.*pointer\\s*=\\s*eval",fun.str)!=0)){
            fun.str <- gsub("(.*)eval\\(sprintf.*?,(.*)\\)\\).*","\\1\\2",fun.str)
        }

        ## change references to sub-matricies
        matricies <- c("global_model_switch","globalStructure\\$gxt","mat","gmaxxt","gminxt",
                       "globalStructure\\$ga","gmaxa","gmina","globalStructure\\$gx","inters",
                       "design\\$G","design\\$Ga","design\\$Gx",
                       "globalStructure\\$gni",
                       "globalStructure\\$gbpop",
                       "globalStructure\\$gd","t\\(x","t\\(a","\\:ni","\\(ni","groupsize",
                       "model_switch","f1","b_global","m1_switch","notfixed_bpop","bpop_plus","bpop_minus",
                       "bpop","\\W+b","\\W+a","\\W+g","\\W+epsi","df_dbeta","m2_switch",
                       "gradff_switch","g_plus","g_minus"," dff_dg0","gradfg_switch","b_plus","b_minus","dfg_db0",
                       "e_plus","e_minus","dfeps_de0","b_ind_plus","b_ind_minus","temp","\\W+y","dv_dbeta","\\W+l",
                       "dv_db_new","\\W+lh","locc\\[\\[\\w+\\]\\]","notfixed_d","d_plus","d_minus","notfixed_covd",
                       "notfixed_docc","docc_plus","docc_minus","notfixed_covdoc","notfixed_sigma","sigma_plus","sigma_minus",
                       "sigma","$notfixed_covsigma","f2","mv","mm4","tmp_fim","m2_tmp","m3_tmp","xt_ind","xt_new","optsw",
                       "bpopdescr","ddescr","docc","sigma_d","itvector","dmfvector","\\s+d","\\=params","\\sparams","params_cv",
                       "\\(params",
                       "\\/params","var_derivative","param_vars","[^\\_\\w]params","m_perm","sn_perm",
                       "[^[:alnum:]\\_]xt","\\Wimft","[^[:alnum:]\\_]axt","\\Wxt_plus","\\Wir","\\Wgdmf","\\Wni_var",
                       "\\Wvaropt","\\Wavar","\\Wmaxvar","\\Wminvar","\\Wvaropto","\\Wnormgvar","\\Wnavar","\\Wgradvar",
                       "\\Wxtopt","\\Wminxt","\\Wmaxxt","\\Waopt","\\Wmina","\\Wmaxa","\\Wx","\\Wxtopto","\\Wbest\\_xt",
                       "\\Wa\\_plus","\\Waa","\\Wk_perm","\\Wline\\_opta","\\Wbest\\_a","\\Wdiscrete\\_val",
                       "\\Wline\\_optx","\\Wbest\\_x","\\Wbpop\\_gen")

        ##subMatricies <- paste(matricies,"\\(",sep="")
        for(tmp.name in matricies){
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
        if(length(grep("globalStructure$global_model_switch[1:globalStructure$m,1:globalStructure$maxni]=design$model_switch",
                       fun.str,fixed=T))!=0){
            fun.str <- paste("if(is.null(globalStructure$global_model_switch)){",
                             "globalStructure$global_model_switch=design$model_switch",
                             "} else {",
                             fun.str,
                             "}",
                             sep="\n")
        }
        if(length(grep("globalStructure$Ga=design$Ga",fun.str,fixed=T))!=0){
            fun.str <- paste("if(is.null(dim(design$Ga))){",
                             "design$Ga=matrix(data=design$Ga,nrow=1,byrow=TRUE)",
                             "}",
                             fun.str,
                             sep="\n")
        }
        if(length(grep("for(i in 1:size(design$Gx,1)){",fun.str,fixed=T))!=0){
            fun.str <- paste(fun.str,
                             "if(any(size(design$Gx)==0)) next",
                             sep="\n")
        }
        if(length(grep("for(i in 1:globalStructure$NumOcc){",fun.str,fixed=T))!=0){
            fun.str <- paste(fun.str,
                             "if(globalStructure$NumOcc==0) next",
                             sep="\n")
        }
        if(length(grep("for(i in 1:length(globalStructure$notfixed_covsigma)){",fun.str,fixed=T))!=0){
            fun.str <- paste(fun.str,
                             "if(any(size(globalStructure$notfixed_covsigma)==0)) next",
                             sep="\n")
        }
        
        if(length(grep("bpop_plus[i] = bpop_plus[i]+j*h",fun.str,fixed=T))!=0){
            fun.str <- "bpop_plus[i] = complex(real=bpop_plus[i],imaginary=h)"
        }
        if(length(grep("g_plus[i] = g_plus[i]+j*globalStructure$hlf",fun.str,fixed=T))!=0){
            fun.str <- "g_plus[i] = complex(real=g_plus[i],imaginary=globalStructure$hlf)"
        }
        if(length(grep("b_plus[i] = b_plus[i]+j*globalStructure$hlg",fun.str,fixed=T))!=0){
            fun.str <- "b_plus[i] = complex(real=b_plus[i],imaginary=globalStructure$hlg)"
        }
        if(length(grep("e_plus[i] = e_plus[i]+j*globalStructure$hle",fun.str,fixed=T))!=0){
            fun.str <- "e_plus[i] = complex(real=e_plus[i],imaginary=globalStructure$hle)"
        }

        
        ## fun.str <- gsub("","",fun.str,fixed=T)
        fun.str <- gsub("fprintf(fn,\\'\\%\\s \\<\\- function\\(\\)\n\n',strReturn,strFunctionName)\\{",
                        "fprintf(fn,'%s <- function()\n\n',strReturn,strFunctionName)",fun.str,fixed=T)
        fun.str <- gsub("matrix(c(),nrow=1,byrow=T)","matrix(0,0,0)",fun.str,fixed=T)
        fun.str <- gsub("hold on","##hold on",fun.str,fixed=T)
        fun.str <- gsub("hold off","##hold off",fun.str,fixed=T)
        fun.str <- gsub("^(\\s*)tic(\\s*)","\\1tic()\\2",fun.str)
        fun.str <- gsub("(\\s+)toc(\\s*)","\\1toc()\\2",fun.str)
        fun.str <- gsub("=toc(\\s*)","=toc()\\2",fun.str)
        fun.str <- gsub("if(tmpfile=='')==0{","if(!(tmpfile=='')){",fun.str,fixed=T)
        fun.str <- gsub("if(out_file=='')==0{","if(!(out_file=='')){",fun.str,fixed=T)
        fun.str <- gsub("else {%STORING DESIGNS FOR PARALLEL EXECUTION","else { #STORING DESIGNS FOR PARALLEL EXECUTION",fun.str,fixed=T)
        fun.str <- gsub("temp_detmf = designout[[it]].ofv","temp_detmf = designout[[it]]$ofv",fun.str,fixed=T)
        fun.str <- gsub("temp_mf = designout[[it]].FIM","temp_mf = designout[[it]]$FIM",fun.str,fixed=T)
        fun.str <- gsub("temp_detmf = designout[[itx]].ofv","temp_detmf = designout[[itx]]$ofv",fun.str,fixed=T)
        fun.str <- gsub("temp_mf = designout[[itx]].FIM","temp_mf = designout[[itx]]$FIM",fun.str,fixed=T)
        fun.str <- gsub("x = designsin[[itx]].x","x = designsin[[itx]]$x",fun.str,fixed=T)

        fun.str <- gsub("temp_detmf = designout[[ita]].ofv","temp_detmf = designout[[ita]]$ofv",fun.str,fixed=T)
        fun.str <- gsub("temp_mf = designout[[ita]].FIM","temp_mf = designout[[ita]]$FIM",fun.str,fixed=T)
        fun.str <- gsub("a = designsin[[ita]].a","a = designsin[[ita]]$a",fun.str,fixed=T)

        fun.str <- gsub("mf_tmp = designout[[it]].FIM","mf_tmp = designout[[it]]$FIM",fun.str,fixed=T)
        fun.str <- gsub("mf_plus = designout[[it]].FIM","mf_plus = designout[[it]]$FIM",fun.str,fixed=T)
        
        
        fun.str <- gsub("globalStructure$sigma = popedInput$design$sigma",
                        "globalStructure$sigma = if (is.null(dim(popedInput$design$sigma))) {matrix(popedInput$design$sigma,nrow=1,byrow=T)} else {popedInput$design$sigma}",fun.str,fixed=T)
        fun.str <- gsub("bfgs_init=matrix(c(),nrow=1,byrow=T)","bfgs_init=t(zeros(0,1))",fun.str,fixed=T)
        fun.str <- gsub("lgvar'","t(lgvar)",fun.str,fixed=T)
        fun.str <- gsub("a[i,,drop=F]'","t(a[i,,drop=F])",fun.str,fixed=T)
        fun.str <- gsub("a_plus[i,,drop=F]'","t(a_plus[i,,drop=F])",fun.str,fixed=T)
        fun.str <- gsub(" while(it<=globalStructure$sgit) && (abs(diff)>globalStructure$convergence_eps){",
                        " while((it<=globalStructure$sgit) && (abs(diff)>globalStructure$convergence_eps)){",fun.str,fixed=T)
        fun.str <- gsub("(1:numel(aopt))'","t(1:numel(aopt))",fun.str,fixed=T)

        fun.str <- gsub("xt_plus(i,1:ni[i,drop=F])'","t(xt_plus(i,1:ni[i,drop=F]))",fun.str,fixed=T)
        fun.str <- gsub("xt_plus[i,1:ni[i,drop=F],drop=F]'","t(xt_plus[i,1:ni[i,drop=F],drop=F])",fun.str,fixed=T)

        fun.str <- gsub("x0=x0'","x0=t(x0)",fun.str,fixed=T)
        fun.str <- gsub("l=l'","l=t(l)",fun.str,fixed=T)
        fun.str <- gsub("u=u'","u=t(u)",fun.str,fixed=T)

        fun.str <- gsub("bfgsb_min(f_handle, x_k,lb,ub,options)","bfgsb_min(f_name,f_options, x_k,lb,ub,options)",fun.str,fixed=T)

                
        fun.str <- gsub("(1:numel(xtopt))'","t(1:numel(xtopt))",fun.str,fixed=T)
        fun.str <- gsub("f_handle=@(x)calc_ofv_and_grad(x,optxt, opta, model_switch,aa,axt,globalStructure$groupsize,ni,xtopt,xopt,aopt,bpop,d,globalStructure$sigma,docc_full,globalStructure)",
                        "f_name <- 'calc_ofv_and_grad' \n f_options <- list(x,optxt, opta, model_switch,aa,axt,globalStructure$groupsize,ni,xtopt,xopt,aopt,bpop,d,globalStructure$sigma,docc_full,globalStructure)",
                        fun.str,fixed=T)
        fun.str <- gsub("can''t","can not",fun.str,fixed=T)
        fun.str <- gsub("designsout[[it]].","designsout[[it]]$",fun.str,fixed=T)
        fun.str <- gsub("designsin[[it]].","designsin[[it]]$",fun.str,fixed=T)
        fun.str <- gsub("clear itvector","itvector <- c()",fun.str,fixed=T)        
        fun.str <- gsub("clear dmfvector","dmfvector <- c()",fun.str,fixed=T)
        fun.str <- gsub("1:numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma",
                        "1:(numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)",fun.str,fixed=T)
        fun.str <- gsub("trace_matrix(reshape_matlab(m2_tmp[,m,drop=F],n,n)*invv*reshape_matlab(m2_tmp[,k,drop=F],n,n)*invv",
                        "trace_matrix(reshape_matlab(m2_tmp[,m,drop=F],n,n)%*%invv%*%reshape_matlab(m2_tmp[,k,drop=F],n,n)%*%invv",
                        fun.str,fixed=T)
        fun.str <- gsub("trace_matrix(reshape_matlab(m3_tmp[,m,drop=F],n,n)*invv*reshape_matlab(m3_tmp[,k,drop=F],n,n)*invv",
                        "trace_matrix(reshape_matlab(m3_tmp[,m,drop=F],n,n)%*%invv%*%reshape_matlab(m3_tmp[,k,drop=F],n,n)%*%invv",
                        fun.str,fixed=T)
        fun.str <- gsub("reshape_matlab(m3_tmp[,m,drop=F],n,n)*invv*reshape_matlab(m2_tmp[,k,drop=F],n,n)*invv",
                        "reshape_matlab(m3_tmp[,m,drop=F],n,n)%*%invv%*%reshape_matlab(m2_tmp[,k,drop=F],n,n)%*%invv",
                        fun.str,fixed=T)
        fun.str <- gsub("2*t(m1_tmp)*invv*m1_tmp","2*t(m1_tmp)%*%invv%*%m1_tmp",fun.str,fixed=T)
        fun.str <- gsub("t(f1)*f2*f1","t(f1)%*%f2%*%f1",fun.str,fixed=T)
        fun.str <- gsub("(ct1-1)*n+1:(ct1-1)*n+n,(ct2-1)*n+1:(ct2-1)*n+n",
                        "((ct1-1)*n+1):((ct1-1)*n+n),((ct2-1)*n+1):((ct2-1)*n+n)",fun.str,fixed=T)
        fun.str <- gsub("numnotfixed_bpop+1:numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma",
                        "(numnotfixed_bpop+1):(numnotfixed_bpop+numnotfixed_d+numnotfixed_covd+numnotfixed_docc+numnotfixed_covdocc+numnotfixed_sigma+numnotfixed_covsigma)",fun.str,fixed=T)
        fun.str <- gsub("n+1:n+n*n","(n+1):(n+n*n)",fun.str,fixed=T)
        fun.str <- gsub("h*sigma*t(h)","h%*%sigma%*%t(h)",fun.str,fixed=T)
        
        fun.str <- gsub("(t(lgvar)*lgvaro<0)","any(t(lgvar)%*%lgvaro<0)",fun.str,fixed=T)

        fun.str <- gsub("OFV of %d","OFV of %g",fun.str,fixed=T)


        fun.str <- gsub("if((size(globalStructure$prior_fim) == size(mft))){","if(all(size(globalStructure$prior_fim) == size(mft))){",fun.str,fixed=T)
        

        fun.str <- gsub("x[1:numel(xtopto[notfixed])]=matrix(0,0,0)","x=x[-c(1:numel(xtopto[notfixed]))]",fun.str,fixed=T)

        fun.str <- gsub("(size(globalStructure$prior_fim)==size(mft))","all(size(globalStructure$prior_fim)==size(mft))",fun.str,fixed=T)

        fun.str <- gsub("l*d*t(l)","l%*%d%*%t(l)",fun.str,fixed=T)
        fun.str <- gsub("imft[ct2,,drop=F]*ir[,ct2,drop=F]","imft[ct2,,drop=F]%*%ir[,ct2,drop=F]",fun.str,fixed=T)
        fun.str <- gsub("grad_ff_tmp*gradfg","grad_ff_tmp%*%gradfg",fun.str,fixed=T)
        fun.str <- gsub("zeros(size(b_ind))","zeros(size(b_ind)[1],size(b_ind)[2])",fun.str,fixed=T)
        fun.str <- gsub("zeros(size(bocc_ind))","zeros(size(bocc_ind)[1],size(bocc_ind)[2])",fun.str,fixed=T)
        fun.str <- gsub("globalStructure$gni[1:globalStructure$m,drop=F]","globalStructure$gni[1:globalStructure$m,,drop=F]",fun.str,fixed=T)
                
        fun.str <- gsub("locc_tmp*docc*t(locc_tmp)","locc_tmp%*%docc%*%t(locc_tmp)",fun.str,fixed=T)

        fun.str <- gsub("struct('factr',globalStructure$BFGSConvergenceCriteriaMinStep,'pgtol',globalStructure$BFGSProjectedGradientTol,'ftol',globalStructure$BFGSTolerancef,'gtol',globalStructure$BFGSToleranceg,'xtol',globalStructure$BFGSTolerancex)",
                        "list('factr'=globalStructure$BFGSConvergenceCriteriaMinStep,'pgtol'=globalStructure$BFGSProjectedGradientTol,'ftol'=globalStructure$BFGSTolerancef,'gtol'=globalStructure$BFGSToleranceg,'xtol'=globalStructure$BFGSTolerancex)",fun.str,fixed=T)
        
        fun.str <- gsub("lh*kron_tmp(d,sigma)*t(lh)","lh%*%kron_tmp(d,sigma)%*%t(lh)",fun.str,fixed=T)
        fun.str <- gsub("t(ni)*matrix(1,globalStructure$m,1)","t(ni)%*%matrix(1,globalStructure$m,1)",fun.str,fixed=T)


        fun.str <- gsub("in group nr #g","in group nr %g",fun.str,fixed=T)

        fun.str <- gsub("index(minxt==maxxt)=matrix(0,0,0)","index=index[minxt!=maxxt]",fun.str,fixed=T)
        fun.str <- gsub("index(mina==maxa)=matrix(0,0,0)","index=index[mina!=maxa]",fun.str,fixed=T)

        fun.str <- gsub("]].FIM","]]$FIM",fun.str,fixed=T)
        
        fun.str <- gsub("gradofv_xt[model_switch,axt,globalStructure$groupsize,ni,xtopt,xopt,aopt,bpop,d,globalStructure$sigma,docc_full,globalStructure]","gradofv_xt(model_switch,axt,globalStructure$groupsize,ni,xtopt,xopt,aopt,bpop,d,globalStructure$sigma,docc_full,globalStructure)",fun.str,fixed=T)

        

        fun.str <- gsub("getfulld(d(,2),globalStructure$covd)","getfulld(d[,2,drop=F],globalStructure$covd)",fun.str,fixed=T)
        

        if(basename(from)=="LinMatrixL_occ.m") fun.str <- gsub("y=0","y=zeros(length(xt_ind),0)",fun.str,fixed=T)
        if(basename(from)=="LinMatrixLH.m") fun.str <- gsub("c(t(b_ind) e0","c(t(b_ind), e0",fun.str,fixed=T)
        if(basename(from)=="LinMatrixLH.m") fun.str <- gsub("deriv(j).hx","deriv(j)$hx",fun.str,fixed=T)
        if(basename(from)=="LinMatrixLH.m") fun.str <- gsub("(i-1)*NumEPS+1:i*NumEPS","((i-1)*NumEPS+1):(i*NumEPS)",fun.str,fixed=T)
        ## if(length(grep("for(j in 1:size(design$Gx,2)){",fun.str,fixed=T))!=0){
        ##     fun.str <- paste(fun.str,
        ##                      "if(size(design$Gx,2)==0) next",
        ##                      sep="\n")
        ## }
                                 #   g/\/s//solve(,)/g
                                        #   g/fsolve('\(..*\)'/s//ms(~\1 /g
                                        #   g/param(\(..*\))/s//param[ \1 ] /g
                                        #   g/var(\(..*\))/s//var[ \1 ] /g
                                        #   g/mod1(\(..*\)/s//mod1[ \1 /g
        
        file.str[i] <- fun.str
    }

    ## add return argument to end of function
    ret.args.2 <- gsub("^list\\(","",ret.args)
    ret.args.2 <- gsub("\\)$","",ret.args.2)
    ret.args.2 <- gsub("([^,]*),","\\1=\\1,",ret.args.2)
    ret.args.2 <- gsub(",([^,]*)$",",\\1=\\1",ret.args.2)
    ret.args <- gsub("^list\\((.*)\\)",paste("list(",ret.args.2,")",sep=""),ret.args)


    only.spaces <- grep("^[[:space:]]*$",file.str)
    only.comments <- grep("^\\s*\\#.*",file.str)
    only.scrap <- c(only.spaces,only.comments)
    if(length(only.scrap)!=0){
        foo <- 1:length(file.str)
        foo <- foo[-only.scrap]
        last.line.with.code <- foo[length(foo)]
    } else {
        last.line.with.code <- length(file.str)
    }

    ##if(file.str[length(file.str)]=="}"){
    if(file.str[last.line.with.code]=="}"){
        ##file.str[length(file.str)]=paste("return(",ret.args,")",sep="")
        file.str[last.line.with.code]=paste("return(",ret.args,") \n}",sep="")
        ##file.str[length(file.str)+1]="}"
        ##file.str[last.line.with.code+1]="}"
    } else {
        #cat("return argument needs to be added")
        file.str[length(file.str)+1]=paste("return(",ret.args,")\n}",sep="")
    }

    ## create empty list for function input
    if(function.input){
        file.str.2 <- file.str
        file.str.2[length(file.str)+1] <- ""
        file.str.2[2] <- "popedInput <- list()"
        file.str.2[-2] <- file.str
        file.str <- file.str.2
    }
    if(convert.input){
        file.str.2 <- file.str
        file.str.2[length(file.str)+1] <- ""
        file.str.2[2] <- "globalStructure <- list()"
        file.str.2[-2] <- file.str
        file.str <- file.str.2
    }

    to <- gsub(".m$",".R",from)
    writeLines(file.str, con=to)
    ##file.str[1:10]
    ##file.str
    ##file.str[c(62,95:98)]
    return(file.str)
}


