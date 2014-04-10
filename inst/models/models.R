bolus_1cpt_Vk<-function()
{
  form1<-paste("dose/V*(exp(-k*t))")
  form1<-parse(text=form1,n=-1)
  tf<-Inf
  return(list(c(form1),tf))
}