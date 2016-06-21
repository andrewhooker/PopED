library(PopED)

design_1 <- create_design(xt=list(c(1,2,3,4,5),
                                  c(1,2,3,4)),
                          groupsize=c(50,20),
                          a=list(c(WT=70,DOSE=1000),
                                 c(DOSE=1000,WT=35)))

ds_1 <- create_design_space(design_1)

ds_1_a <- create_design_space(design_1,our_zero = 1e-5)

ds_2 <- create_design_space(design_1,maxni=10,maxxt=10,minxt=0)

ds_3 <- create_design_space(design_1,maxni=10,mingroupsize=20,maxxt=10,minxt=0)

ds_4 <- create_design_space(design_1,maxa=c(100,2000))

ds_5 <- create_design_space(design_1,mina=c(10,20))

design_2 <- create_design(xt=list(c(1,2,3,4,5),
                                  c(1,2,3,4)),
                          groupsize=c(50,20),
                          a=list(c(WT=70,DOSE=1000),
                                 c(WT=35,DOSE=1000)),
                          x=list(c(SEX=1,DOSE_discrete=100),
                                 c(SEX=2,DOSE_discrete=200)))

ds_6 <- create_design_space(design_2) 

ds_7 <- create_design_space(design_2,
                            x_space=list(SEX=c(1,2),
                                         DOSE_discrete=seq(100,400,by=20)))

ds_8 <- create_design_space(design_2,
                            x_space=list(SEX=c(1,2),
                                         DOSE_discrete=seq(100,400,by=20)),
                            grouped_xt=c(1,2,3,4,5))

ds_9 <- create_design_space(design_2,
                            x_space=list(SEX=c(1,2),
                                         DOSE_discrete=seq(100,400,by=20)),
                            use_grouped_xt=TRUE)

design_3 <- create_design(xt=list(c(1,2,3,4,5),
                                  c(1,2,3,4)),
                          groupsize=c(50,20),
                          a=list(c(WT=35,DOSE=1000)),
                          x=list(c(SEX=1,DOSE_discrete=100)))

ds_10 <- create_design_space(design_3,
                             x_space=list(SEX=c(1,2),DOSE_discrete=seq(100,400,by=20)),
                             use_grouped_a=TRUE)

ds_11 <- create_design_space(design_2,
                             x_space=list(SEX=c(1,2),DOSE_discrete=seq(100,400,by=20)),
                             grouped_a=list(c(1,2),c(3,2)))

ds_12 <- create_design_space(design_3,
                             x_space=list(SEX=c(1,2),DOSE_discrete=seq(100,400,by=20)),
                             use_grouped_x=TRUE)

ds_13 <- create_design_space(design_3,
                             x_space=list(SEX=c(1,2),DOSE_discrete=seq(100,400,by=20)),
                             grouped_x=list(c(1,2),c(3,2)))

seq_1 <- 1:10
ds_14 <- create_design_space(design_1,maxxt=10,minxt=0,
                             xt_space = list(seq_1,seq_1,seq_1,seq_1,seq_1))
ds_15 <- create_design_space(design_1,maxxt=10,minxt=0,xt_space = list(seq_1))

possible_values <- as.matrix(cbind(list(0:10),list(0:10),list(0:10),list(0:20),list(0:20)))
xt_space <- as.matrix(rbind(possible_values,possible_values))

ds_16 <- create_design_space(design_1,maxxt=10,minxt=0,xt_space = xt_space)

ds_17 <- create_design_space(design_1,a_space = list(1:100,seq(1000,100000,by=1000)))

