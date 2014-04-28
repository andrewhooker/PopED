# warfarin optimize model

# Adding 10% Uncertainty to fixed effects log-normal (not Favail)
bpop_vals <- c(CL=0.15, V=8, KA=1.0, Favail=1)
bpop_vals_ed_ln <- cbind(ones(length(bpop_vals),1)*4, # log-normal distribution
                      bpop_vals,
                      ones(length(bpop_vals),1)*(bpop_vals*0.1)^2) # 10% of bpop value
bpop_vals_ed_ln["Favail",]  <- c(0,1,0)
bpop_vals_ed_ln


# with log-normal distributions
pars.ln <- pargen(par=bpop_vals_ed_ln,
               user_dist_pointer=NULL,
               sample_size=100,
               bLHS=1,
               sample_number=NULL,
               poped.db)
#looks ok
colMeans(pars.ln)
var(pars.ln)

mean.diff.ln <- (colMeans(pars.ln) - bpop_vals_ed_ln[,2])/bpop_vals_ed_ln[,2]*100
mean.diff.ln
var.diff.ln <- (diag(var(pars.ln)) - bpop_vals_ed_ln[,3])/bpop_vals_ed_ln[,3]*100
var.diff.ln



# Adding 10% Uncertainty to fixed effects normal-distribution (not Favail)
bpop_vals_ed_n <- cbind(ones(length(bpop_vals),1)*1, # log-normal distribution
                      bpop_vals,
                      ones(length(bpop_vals),1)*(bpop_vals*0.1)^2) # 10% of bpop value
bpop_vals_ed_n["Favail",]  <- c(0,1,0)
bpop_vals_ed_n

# with normal distributions
pars.n <- pargen(par=bpop_vals_ed_n,
               user_dist_pointer=NULL,
               sample_size=100,
               bLHS=1,
               sample_number=NULL,
               poped.db)

# Looks ok
colMeans(pars.n)
var(pars.n)

mean.diff.n <- (colMeans(pars.n) - bpop_vals_ed_n[,2])/bpop_vals_ed_n[,2]*100
mean.diff.n
var.diff.n <- (diag(var(pars.n)) - bpop_vals_ed_n[,3])/bpop_vals_ed_n[,3]*100
var.diff.n


# Adding 10% Uncertainty to fixed effects uniform-distribution (not Favail)
bpop_vals_ed_u <- cbind(ones(length(bpop_vals),1)*2, # uniform distribution
                        bpop_vals,
                        ones(length(bpop_vals),1)*(bpop_vals*0.1)) # 10% of bpop value
bpop_vals_ed_u["Favail",]  <- c(0,1,0)
bpop_vals_ed_u

# with normal distributions
pars.u <- pargen(par=bpop_vals_ed_u,
                 user_dist_pointer=NULL,
                 sample_size=100,
                 bLHS=1,
                 sample_number=NULL,
                 poped.db)

# Looks ok
mean.diff.u <- (colMeans(pars.u) - bpop_vals_ed_u[,2])/bpop_vals_ed_u[,2]*100
mean.diff.u
range.diff.u <- mean.diff.u
for(i in 1:4){
  range.diff.u[i] <- (diff(range(pars.u[,i])) - bpop_vals_ed_u[i,3])/bpop_vals_ed_u[i,3]*100
}
range.diff.u


