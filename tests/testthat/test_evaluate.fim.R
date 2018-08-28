context("FIM calculation")

test_that("RSE from evaluate.fim", {
  
  source("examples_fcn_doc/examples_evaluate.fim.R")
  
  expected.reduced <- c(4.7,2.8,13.9,25.6,30.3,25.8,11.2)
  expected.full <- c(3.6, 2.6,4.8,26.3, 30.9, 26.6, 12.4)
  comp.red.1 <- get_rse(FIM.1,poped.db)
  comp.red.4 <- get_rse(FIM.4,poped.db,fim.calc.type=4)
  comp.full.0 <- get_rse(FIM.0,poped.db)
  comp.red.1.prior <- get_rse(FIM.1.prior, poped.db.prior)
  
  expect_true(all.equal(comp.red.1, expected.reduced,
                        check.attributes = FALSE, check.names = FALSE,
                        tolerance = 0.01))

  expect_true(all.equal(comp.full.0, expected.full,
                        check.attributes = FALSE, check.names = FALSE,
                        tolerance = 0.01))
  
  expect_true(all.equal(head(comp.red.4,-1), head(expected.reduced,-1),
                        check.attributes = FALSE, check.names = FALSE,
                        tolerance = 0.01))
  
  expect_true(all.equal(comp.red.1/sqrt(2), comp.red.1.prior))
  
})

test_that("det(FIM) using evaluate.fim, approx and derivative types", {
  
  source("examples_fcn_doc/warfarin_basic.R")
  
  ## check all FIM calculations
  df <- c()
  for(i in c(0,1,4,5,6,7)){
    for(j in c(0,1)){
      FIM <- evaluate.fim(poped.db,fim.calc.type=i,deriv.type=j) 
      tmp <- data.frame("fim.calc.type"= i, "deriv.type"=j, "det.FIM"=det(FIM))
      df <- rbind(df,tmp)
    }
  }
  #print(df,digits=3,row.names=F,print.gap=3)

  for(i in c(0,1,4,5,6,7)){
    expect_that(df[df$fim.calc.type==i & df$deriv.type==0,]$det.FIM, 
                equals(df[df$fim.calc.type==i & df$deriv.type==1,]$det.FIM,
                       tolerance=1e-5))
  }
  
  expect_that(df[df$fim.calc.type==0 & df$deriv.type==0,]$det.FIM, 
              equals(df[df$fim.calc.type==5 & df$deriv.type==0,]$det.FIM,
                     tolerance=1e-5))
  expect_that(df[df$fim.calc.type==0 & df$deriv.type==0,]$det.FIM, 
              equals(df[df$fim.calc.type==6 & df$deriv.type==0,]$det.FIM,
                     tolerance=1e-5))
  expect_that(df[df$fim.calc.type==1 & df$deriv.type==0,]$det.FIM, 
              equals(df[df$fim.calc.type==7 & df$deriv.type==0,]$det.FIM,
                     tolerance=1e-5))
  
})

test_that("ofv calculation", {
  
  source("examples_fcn_doc/warfarin_optimize.R")
  
  FIM <- evaluate.fim(poped.db) # new name for function needed
  
  expect_that(ofv_fim(FIM,poped.db,ofv_calc_type=1), equals(det(FIM)))
  expect_that(ofv_fim(FIM,poped.db,ofv_calc_type=2), equals(1/trace_matrix(inv(FIM))))
  expect_that(ofv_fim(FIM,poped.db,ofv_calc_type=4), equals(log(det(FIM)))) 
  expect_that(ofv_fim(FIM,poped.db,ofv_calc_type=7), equals(1/sum(get_rse(FIM,poped.db,use_percent=FALSE))))
})

  
test_that("internal FIM calculations", {
  
  source("examples_fcn_doc/warfarin_optimize.R")
  
  ind=1
  
  mf3_calc <- function(){
    mf3(model_switch=t(poped.db$design$model_switch[ind,,drop=FALSE]),
        xt=t(poped.db$design$xt[ind,,drop=FALSE]),
        x=zeros(0,1),
        a=t(poped.db$design$a[ind,,drop=FALSE]),
        bpop=poped.db$parameters$bpop[,2,drop=FALSE],
        d=poped.db$parameters$param.pt.val$d,
        sigma=poped.db$parameters$sigma,
        docc=poped.db$parameters$param.pt.val$docc,
        poped.db)
  }
  
  # in this simple case the full FIM is just the sum of the individual FIMs
  # and all the individual FIMs are the same
  
  poped.db$settings$iFIMCalculationType = 0
  output <- mf3_calc()
  expect_identical(det(output$ret*32),det(evaluate.fim(poped.db,fim.calc.type=0)))
  
  poped.db$settings$iFIMCalculationType = 1
  output <- mf3_calc()
  expect_identical(det(output$ret*32),det(evaluate.fim(poped.db,fim.calc.type=1)))
  
  poped.db$settings$iFIMCalculationType = 4
  output <- mf3_calc()
  expect_identical(det(output$ret*32),det(evaluate.fim(poped.db,fim.calc.type=4)))
  
  poped.db$settings$iFIMCalculationType = 5
  output <- mf3_calc()
  expect_identical(det(output$ret*32),det(evaluate.fim(poped.db,fim.calc.type=5)))
  
  poped.db$settings$iFIMCalculationType = 1
  output <- mf7(model_switch=t(poped.db$design$model_switch[ind,,drop=FALSE]),
                xt=t(poped.db$design$xt[ind,,drop=FALSE]),
                x=zeros(0,1),
                a=t(poped.db$design$a[ind,,drop=FALSE]),
                bpop=poped.db$parameters$bpop[,2,drop=FALSE],
                d=poped.db$parameters$param.pt.val$d,
                sigma=poped.db$parameters$sigma,
                docc=poped.db$parameters$param.pt.val$docc,
                poped.db)
  expect_identical(det(output$ret*32),det(evaluate.fim(poped.db,fim.calc.type=6)))
  
  poped.db$settings$iFIMCalculationType = 7
  output <- mf3_calc()
  expect_identical(det(output$ret*32),det(evaluate.fim(poped.db,fim.calc.type=7)))
  
})

test_that("FIM calculations are as expected", {
  source("examples_fcn_doc/warfarin_basic.R")
  expect_equal(det(evaluate.fim(poped.db,fim.calc.type=0)),1.220371e+24,tolerance=1e-5)
  expect_equal(det(evaluate.fim(poped.db,fim.calc.type=1)),5.996147e+22,tolerance=1e-5)
  expect_equal(det(evaluate.fim(poped.db,fim.calc.type=4)),2.398459e+21,tolerance=1e-5)
  expect_equal(det(evaluate.fim(poped.db,fim.calc.type=5)),1.220371e+24,tolerance=1e-5)
  expect_equal(det(evaluate.fim(poped.db,fim.calc.type=6)),1.220371e+24,tolerance=1e-5)
  expect_equal(det(evaluate.fim(poped.db,fim.calc.type=7)),5.996147e+22,tolerance=1e-5)
})

test_that("group FIM calculations work", {
  # define parameters
  sfg <- function(x,a,bpop,b,bocc){
    parameters=c(CL=bpop[1]*exp(b[1]),
                 V=bpop[2]*exp(b[2]),
                 KA=bpop[3]*exp(b[3]),
                 Favail=bpop[4],
                 DOSE=a[1])
    return(parameters)
  }
  
  ## -- Define initial design  and design space
  poped.db <- create.poped.database(ff_fun=ff.PK.1.comp.oral.sd.CL,
                                    fg_fun=sfg,
                                    fError_fun=feps.add.prop,
                                    bpop=c(CL=0.15, V=8, KA=1.0, Favail=1),
                                    notfixed_bpop=c(1,1,1,0),
                                    d=c(CL=0.07, V=0.02, KA=0.6),
                                    sigma=c(0.01,0.25),
                                    groupsize=32,
                                    xt=list(c(0.5,1,2,6),c(24,36,72,120)),
                                    minxt=0.01,
                                    maxxt=120,
                                    a=rbind(70,50),
                                    m=2,
                                    mina=0.01,
                                    maxa=100)
  
  group_fim <- extract_norm_group_fim(poped.db)
  
  expect_equal(det(group_fim[[1]]*poped.db$design$groupsize[1] + group_fim[[2]]*poped.db$design$groupsize[2]),
               det(evaluate.fim(poped.db)))
})

