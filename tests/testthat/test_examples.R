context("Examples")

test_that("The Examples run", {
  
  if(skip_examples) skip("Examples with long run times")
  
  cur_dir <- getwd()
  run_dir <- file.path(cur_dir,'..','..','inst','examples')
  #setwd(run_dir)

  source(file.path(run_dir,"ex.1.a.PK.1.comp.oral.md.intro.R"),chdir=T)  
  source(file.path(run_dir,"ex.1.b.PK.1.comp.oral.md.re-parameterize.R"),chdir=T)  
  
  
  source(file.path(run_dir,"ex.1.c.PK.1.comp.oral.md.ODE.compiled.R"),chdir=T) 
  source(file.path(run_dir,"ex.2.a.warfarin.evaluate.R"),chdir=T)
  source(file.path(run_dir,"ex.2.b.warfarin.optimize.R"),chdir=T)

  source(file.path(run_dir,"ex.2.c.warfarin.ODE.compiled.R"),chdir=T)  
  expect_equal(design_ode[["ofv"]],design_ode_compiled[["ofv"]],tol=1e-2)

  source(file.path(run_dir,"ex.2.d.warfarin.ED.R"),chdir=T)
  source(file.path(run_dir,"ex.2.e.warfarin.Ds.R"),chdir=T)
  source(file.path(run_dir,"ex.3.a.PKPD.1.comp.oral.md.imax.D-opt.R"),chdir=T)
  source(file.path(run_dir,"ex.3.b.PKPD.1.comp.oral.md.imax.ED-opt.R"),chdir=T)
  source(file.path(run_dir,"ex.4.PKPD.1.comp.emax.R"),chdir=T)
  source(file.path(run_dir,"ex.5.PD.emax.hill.R"),chdir=T)
  source(file.path(run_dir,"ex.6.PK.1.comp.oral.sd.R"),chdir=T)
  source(file.path(run_dir,"ex.7.PK.1.comp.maturation.R"),chdir=T)
  source(file.path(run_dir,"ex.8.tmdd_qss_one_target_compiled.R"),chdir=T)
  source(file.path(run_dir,"ex.9.PK.2.comp.oral.md.ode.compiled.R"),chdir=T)

  source(file.path(run_dir,"ex.10.PKPD.HCV.compiled.R"),chdir=T)
  expect_equal(crit,crit_reference_reduced,tolerance=0.01,scale=crit_reference_reduced)
  for(i in 1:length(rse)) expect_equal(rse[i], rse_reference_reduced[i],tolerance=0.05,scale=rse_reference_reduced[i],check.names=F)
  expect_equal(crit_full,crit_reference_full,tolerance=0.01,scale=crit_reference_full)
  for(i in 1:length(rse_full)) expect_equal(rse_full[i],rse_reference_full[i],tolerance=0.05,scale=rse_reference_full[i],check.names=F)

  source(file.path(run_dir,"ex.10.PKPD.HCV.compiled.R"),chdir=T)
  source(file.path(run_dir,"ex.11.PK.prior.R"),chdir=T)
  source(file.path(run_dir,"ex.12.covariate.distributions.R"),chdir=T)
  source(file.path(run_dir,"ex.13.shrinkage.R"),chdir=T)
  source(file.path(run_dir,"ex.14.PK.IOV.R"),chdir=T)
  source(file.path(run_dir,"ex.15.full.covariance.matrix.R"),chdir=T)
  
  
  
  
  
  
  #setwd(cur_dir)
  
  
})