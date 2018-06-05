context("Pre-Optimization")

test_that("get_par_and_space works", {
  
  
  sfg <- function(x,a,bpop,b,bocc){
    ## -- parameter definition function 
    parameters=c(CL=bpop[1]*exp(b[1]),
                 V=bpop[2]*exp(b[2]),
                 KA=bpop[3]*exp(b[3]),
                 Favail=bpop[4],
                 DOSE=a[1])
    return(parameters) 
  }
  
  ff <- function(model_switch,xt,parameters,poped.db){
    ##-- Model: One comp first order absorption
    with(as.list(parameters),{
      y=xt
      y=(DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*xt)-exp(-KA*xt))
      return(list(y=y,poped.db=poped.db))
    })
  }
  
  feps <- function(model_switch,xt,parameters,epsi,poped.db){
    ## -- Residual Error function
    ## -- Proportional + additive
    y <- ff(model_switch,xt,parameters,poped.db)[[1]] 
    y = y*(1+epsi[,1]) + epsi[,2]  
    return(list(y=y,poped.db=poped.db)) 
  }
  
  
  poped_db_1 <- create.poped.database( 
    ff_fun=ff,
    fg_fun = sfg,
    fError_fun = feps,
    bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
    notfixed_bpop=c(1,1,1,0),
    d=c(CL=0.07, V=0.02, KA=0.6), 
    sigma=c(0.01,0.25),
    xt=list(c(1,2,3,4,5),
            c(1,2,3,4)),
    groupsize=c(50,20),
    a=list(c(WT=70,DOSE=1000),
           c(DOSE=1000,WT=35))
    )
  
  df_1 <- par_and_space_tbl(poped_db_1)
  expect_equal(nrow(df_1),20)
  expect_true(all(df_1$fixed==TRUE))
  expect_true(all(df_1$par==df_1$lower))
  expect_true(all(df_1$par==df_1$upper))
  
  
  poped_db_2 <- create.poped.database(poped_db_1, maxxt=10,minxt=0)
  df_2 <- par_and_space_tbl(poped_db_2)
  #expect_true(all(df_2[df_2$type=="xt","fixed"]==FALSE))
  #expect_true(all(df_2[df_2$type!="xt","fixed"]==TRUE))
  
  library(dplyr)
  df_2 %>% filter(type=="xt") %>% select(fixed) %>% '=='(FALSE) %>% all() %>% expect_true()
  df_2 %>% filter(type!="xt") %>% select(fixed) %>% '=='(TRUE) %>% all() %>% expect_true()
  
  poped_db_3 <- create.poped.database(poped_db_2, 
                                      model_switch = list(c(1,1,1,2,2),
                                                          c(1,1,2,2)))
  df_3 <- par_and_space_tbl(poped_db_3)
  #expect_true(all(df_3[df_3$type=="xt","model"]==c(1,1,1,2,2,1,1,2,2)))
  df_3 %>% filter(type=="xt") %>% select(model) %>% '=='(c(1,1,1,2,2,1,1,2,2)) %>% all() %>% expect_true()
  
  poped_db_4 <- create.poped.database( 
    ff_fun=ff,
    fg_fun = sfg,
    fError_fun = feps,
    bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
    notfixed_bpop=c(1,1,1,0),
    d=c(CL=0.07, V=0.02, KA=0.6), 
    sigma=c(0.01,0.25),
    xt=list(c(1,2,3,4,5),
            c(1,2,3,4)),
    groupsize=c(50,20)
  )
  
  df_4 <- par_and_space_tbl(poped_db_4)
  expect_true(!any(df_4$type=="a"))
  
  
  poped_db_5 <- create.poped.database(poped_db_1, maxxt=10,minxt=0, discrete_xt = list(-1:13))
  df_5 <- par_and_space_tbl(poped_db_5)
  expect_true(all(df_5$allowed_values[[1]]==c(1:10)))
  
  poped_db_6 <- create.poped.database(poped_db_5, ourzero = 0, minxt=0)
  df_6 <- par_and_space_tbl(poped_db_6)
  expect_true(all(df_6$allowed_values[[1]]==c(0:10)))
  
  
  poped_db_7 <- create.poped.database(poped_db_1,a=list(c(70,1000),c(35,1000)))
  df_7 <- par_and_space_tbl(poped_db_7)
  df_7 %>% filter(type=="a") %>% select(name) %>% is.na() %>% all() %>% expect_true()

  # ds_2 <- create_design_space(design_1,maxni=10,maxxt=10,minxt=0)
  # 
  # ds_3 <- create_design_space(design_1,maxni=10,mingroupsize=20,maxxt=10,minxt=0)
  # 
  # ds_4 <- create_design_space(design_1,maxa=c(100,2000))
  # 
  # ds_5 <- create_design_space(design_1,mina=c(10,20))
  # 
  # design_2 <- create_design(xt=list(c(1,2,3,4,5),
  #                                   c(1,2,3,4)),
  #                           groupsize=c(50,20),
  #                           a=list(c(WT=70,DOSE=1000),
  #                                  c(WT=35,DOSE=1000)),
  #                           x=list(c(SEX=1,DOSE_discrete=100),
  #                                  c(SEX=2,DOSE_discrete=200)))
  # 
  # ds_6 <- create_design_space(design_2) 
  # 
  # ds_7 <- create_design_space(design_2,
  #                             x_space=list(SEX=c(1,2),
  #                                          DOSE_discrete=seq(100,400,by=20)))
  # 
  # ds_8 <- create_design_space(design_2,
  #                             x_space=list(SEX=c(1,2),
  #                                          DOSE_discrete=seq(100,400,by=20)),
  #                             grouped_xt=c(1,2,3,4,5))
  # 
  # ds_9 <- create_design_space(design_2,
  #                             x_space=list(SEX=c(1,2),
  #                                          DOSE_discrete=seq(100,400,by=20)),
  #                             use_grouped_xt=TRUE)
  # 
  # design_3 <- create_design(xt=list(c(1,2,3,4,5),
  #                                   c(1,2,3,4)),
  #                           groupsize=c(50,20),
  #                           a=list(c(WT=35,DOSE=1000)),
  #                           x=list(c(SEX=1,DOSE_discrete=100)))
  # 
  # ds_10 <- create_design_space(design_3,
  #                              x_space=list(SEX=c(1,2),DOSE_discrete=seq(100,400,by=20)),
  #                              use_grouped_a=TRUE)
  # 
  # ds_11 <- create_design_space(design_2,
  #                              x_space=list(SEX=c(1,2),DOSE_discrete=seq(100,400,by=20)),
  #                              grouped_a=list(c(1,2),c(3,2)))
  # 
  # ds_12 <- create_design_space(design_3,
  #                              x_space=list(SEX=c(1,2),DOSE_discrete=seq(100,400,by=20)),
  #                              use_grouped_x=TRUE)
  # 
  # ds_13 <- create_design_space(design_3,
  #                              x_space=list(SEX=c(1,2),DOSE_discrete=seq(100,400,by=20)),
  #                              grouped_x=list(c(1,2),c(3,2)))
  # 
  # seq_1 <- 1:10
  # ds_14 <- create_design_space(design_1,maxxt=10,minxt=0,
  #                              xt_space = list(seq_1,seq_1,seq_1,seq_1,seq_1))
  # ds_15 <- create_design_space(design_1,maxxt=10,minxt=0,xt_space = list(seq_1))
  # 
  # possible_values <- as.matrix(cbind(list(0:10),list(0:10),list(0:10),list(0:20),list(0:20)))
  # xt_space <- as.matrix(rbind(possible_values,possible_values))
  # 
  # ds_16 <- create_design_space(design_1,maxxt=10,minxt=0,xt_space = xt_space)
  # 
  # ds_17 <- create_design_space(design_1,a_space = list(1:100,seq(1000,100000,by=1000)))
  
})




