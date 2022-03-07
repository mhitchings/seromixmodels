y_name = 'logigg' # Outcome variable to fit to
b_vars = c() # Covariates on the log-odds of seropositivity
g_vars = c() # Covariates on the mean seroresponse among seropositives
re = 1 # Whether to include a household (or other level) random effect
model_nonresponse=1 #distribution of non-responders (0=normal,1=log-logistic)
model_response=0 #distribution of responders (0=normal,1=lognormal,2=gamma)
negprior=1 #whether to use pre-pandemic negatives to inform priors
prior_sd=0.05 # Limits of the uniform distribution around prior parameters
loq=1 # whether there is a limit of quantification
gd_prior = 2 # Prior mean of the difference between IgG among seronegatives and IgG among seropositives

require(rstan)
require(survival)
require(loo)
options(mc.cores = parallel::detectCores()) # checks number of cores without having later to specify the cores argument

baseFile = 'stanmod'

# Dataset: one row per individual, with an outcome (y_name), covariates (contained in b_vars and/or g_vars),
# household ID if appropriate

d_foranalysis = read.csv('/path/to/file/')

# If using pre-pandemic negatives, read in file. Also has outcome (y_name)
if (negprior==1) {
  neg=read.csv('/path/to/negative/controls/')
  
  # Fit a model to the pre-pandemic negatives, e.g. log-logistic
  neg$survd=1
  
  llogfit = survreg(Surv(logiggtrans,survd)~1,data=neg,dist='loglogistic')
  # Parameters describing distribution of negatives
  lambda_llog = exp(-coef(llogfit)[1])
  gamma_llog = exp(-log(llogfit$scale))
  
  
}

# Set the numerical value of the LLOQ on the scale of y
if (loq==1) {
  loq_n = log(21)
} else {
  loq_n = 0
}


d_foranalysis = d_foranalysis[!is.na(d_foranalysis[,y_name]),]

# Numerical house IDs from 1:num_hh
d_foranalysis = d_foranalysis %>% mutate(houseid_num = as.numeric(factor(houseid))) %>% arrange(houseid_num)
hhids = d_foranalysis$houseid_num
num_hh = length(unique(hhids))

# Matrix for covariates on mean seroresponse
if (length(g_vars)>0) {
  reg_vars_g=d_foranalysis[,g_vars]
} else {
  reg_vars_g=array(0,dim=c(nrow(d_foranalysis),0))
}

if (re==0) {
  data_prov_basic <-
    list(
      N=nrow(d_foranalysis), #the number of serosamples
      r_b=length(b_vars),# number of covariates on log-odds of seropositivity
      r_g=length(g_vars), # number of covariates on mean seroresponse
      y=d_foranalysis[,y_name], #antibody levels
      reg_vars_b=cbind(rep(1,nrow(d_foranalysis)),
                       d_foranalysis[,b_vars]
      ), #covariates on seroresponse probability
      reg_vars_g=reg_vars_g, #covariates on mean seroresponse
      model_nonresponse=model_nonresponse, #distribution of non-responders (0=normal,1=log-logistic)
      model_response=model_response, #distribution of responders (0=normal,1=lognormal,2=gamma)
      negprior=1, #whether to use pre-pandemic negatives to inform priors
      negprior_params = c(unname(lambda_llog),unname(log(gamma_llog-1))
                          ,prior_sd), #prior parameters (rate, shape, and variation of prior dist)
      loq=loq_n, #lower limit of quantification
      gammadiff_prior = gd_prior # prior mean of difference between seronegative seroresponse and seropositive
    )
  
  # Compile the stan model
  compiled_model_basic = stan_model(file="mixmod_withloq.stan")
  
  # Run the model
  system.time(
    mm_basic <-
      sampling(
        compiled_model_basic,
        data = data_prov_basic,
        chains = 4,
        iter = 5000,
        # warmup = ,
        thin = 1,
        control = list(adapt_delta = 0.80)
      )
  )
  mm_basic@stanmodel@dso <- new("cxxdso")
  saveRDS(mm_basic, file = paste0(baseFile,'.rds'))
}

if (re==1) {
  data_prov_re <-
    list(
      N=nrow(d_foranalysis), #the number of serosamples
      r_b=length(b_vars),# number of covariates on log-odds of seropositivity
      r_g=length(g_vars), # number of covariates on mean seroresponse
      y=d_foranalysis[,y_name], #antibody levels
      reg_vars_b=cbind(rep(1,nrow(d_foranalysis)),
                       d_foranalysis[,b_vars]
      ), #covariates on seroresponse probability
      reg_vars_g=reg_vars_g, #covariates on mean seroresponse
      model_nonresponse=model_nonresponse, #distribution of non-responders (0=normal,1=log-logistic)
      model_response=model_response, #distribution of responders (0=normal,1=lognormal,2=gamma)
      negprior=1, #whether to use pre-pandemic negatives to inform priors
      negprior_params = c(unname(lambda_llog),unname(log(gamma_llog-1))
                          ,prior_sd), #prior parameters (rate, shape, and variation of prior dist)
      num_hh=num_hh, # number of households
      hhids=hhids, # numerical house IDs from 1:num_hh
      loq=loq_n, #lower limit of quantification
      gammadiff_prior = gd_prior # prior mean of difference between seronegative seroresponse and seropositive
    )
  
  # Compile the stan model
  compiled_model_re = stan_model(file="mixmod_raneff_withloq.stan")
  
  # Run the model
  system.time(
    mm_re <-
      sampling(
        compiled_model_re,
        data = data_prov_re,
        chains = 4,
        iter = 5000,
        # warmup = ,
        thin = 1,
        control = list(adapt_delta = 0.80)
      )
  )
  mm_re@stanmodel@dso <- new("cxxdso")
  saveRDS(mm_re, file = paste0(baseFile,'.rds'))
}

# Do some post-processing
m = readRDS(paste0(baseFile,'.rds'))
print(m)
# Look at convergence of the model (Rhat values)
hist(summary(m)$summary[,"Rhat"],main=paste0(b_vars,collapse='.'))

# Calculate leave-one-out information criteria (LOO) using look package
looic=loo(m,pars='log_lik')

# Get matrix of parameters and summarize
pars = as.data.frame(as.matrix(m))
pars = pars[,c("sigma0","mean_nonresponse","sigma","gamma_diff",grep("beta_covs",colnames(pars),value=T),grep("gamma_covs",colnames(pars),value=T))]

if(length(b_vars)>0){covb_fornames = paste0('beta: ',b_vars)}else{covb_fornames = NULL}
if(length(g_vars)>0){covg_fornames = paste0('gamma: ',g_vars)}else{covg_fornames = NULL}
if (model_nonresponse!=1) {
  colnames(pars) = c('sd_nonresp','log(mean_nonresp)','sd_resp','mean_resp','p_resp',covb_fornames,covg_fornames)
} else {
  colnames(pars) = c('shape_llog','rate_llog','sd_resp','mean_resp','p_resp',covb_fornames,covg_fornames)
}

# Transform parameters
pars$sd_resp=exp(pars$sd_resp)
pars$p_resp=expit(pars$p_resp)
if (model_nonresponse!=1){
  pars$sd_nonresp=exp(pars$sd_nonresp)
  pars$`log(mean_resp)`=pars$`log(mean_nonresp)`+pars$`log(mean_resp)`
} else {
  pars$shape_llog=1+exp(pars$shape_llog)
  pars$`log(mean_resp)` = pi/pars$shape_llog / (pars$rate_llog*sin(pi/pars$shape_llog)) + pars$mean_resp
}

# Median and 95% credible interval of parameters
parfit = apply(pars,2,function(x) quantile(x,c(0.025,0.5,0.975)))
print(parfit)

# Make a table of the regression parameters
xb = round(exp(summary(m,pars=c('beta_covs'))$summary[-1,c("2.5%","50%","97.5%")]),2)
xg = round(summary(m,pars=c('gamma_covs'))$summary[,c("2.5%","50%","97.5%")],2)
x = rbind(xb,xg)

res_table_temp = data.frame('Var'=c(b_vars,paste0('agecat_',c(1,3:5)),paste0('timecat_',c(1,3:5))),
                            'Estimate'=paste0(x[,2],' (',x[,1],',',x[,3],')'))
res_table = rbind(res_table,res_table_temp)


# Simulate the density for the random effect model
# Average pdf for the sample is the mean of pdf across individuals
overall_density = function(data,y_name,b_vars,g_vars,parfit,loq_n,x) {
  
  y = data[,y_name]
  
  hhid=data[,'houseid_num']
  
  beta_covs=data[,grepl(paste0(b_vars,collapse='|'),colnames(data))]
  beta_covs = beta_covs[order(match(colnames(beta_covs),b_vars))]
  
  gamma_covs=data[,grepl(paste0(g_vars,collapse='|'),colnames(data))]
  gamma_covs = gamma_covs[order(match(colnames(gamma_covs),b_vars))]
  
  beta_hh = parfit[["beta0"]] + as.numeric(parfit[,paste0('hh',hhid)])
  beta_i = beta_hh + rowSums(data.frame(mapply(`*`,beta_covs,parfit[grepl('beta: ',names(parfit))])))
  
  gamma_i = parfit[['gamma0']]
  
  if (model_nonresponse!=1) {
    gamma_i = gamma_i + parfit[['mean_nonresp']]
  } else {
    gamma_i = gamma_i + pi/(1+exp(parfit[["shape_llog"]])) / (parfit[["rate_llog"]]*sin(pi/(1+exp(parfit[["shape_llog"]])))) 
  }
  
  gamma_i = gamma_i + rowSums(data.frame(mapply(`*`,gamma_covs,parfit[grepl('gamma: ',names(parfit))])))
  
  
  if (model_response==2) {
    gamma_shape = exp(parfit[["sigma"]])
    gamma_rate = exp(parfit[["sigma"]])/gamma_i
    y.seropos = expit(beta_i)*
      dgamma(x,shape=gamma_shape,rate=gamma_rate)
    p_resp_belowloq = expit(beta_i)*pgamma(loq_n,shape=gamma_shape,rate=gamma_rate)
  } else if (model_response==1) {
    logmean = log(gamma_i)
    logsd = exp(parfit[["sigma"]])
    y.seropos = sapply(1:length(logmean),function(i) expit(beta_i[i])*
                         dlnorm(x,meanlog=logmean[i],sdlog=logsd))
    p_resp_belowloq = expit(beta_i)*plnorm(loq_n,meanlog=logmean,sdlog=logsd)
  } else if (model_response==0) {
    sd_norm=exp(parfit[["sigma"]])
    y.seropos <- expit(beta_i)*dnorm(x,mean=gamma_i,sd=sd_norm)
    p_resp_belowloq = expit(beta_i)*pnorm(loq_n,mean=gamma_i,sd=sd_norm)
  }
  
  if (model_nonresponse==0) {
    y.seroneg = (1-expit(beta_i))*dnorm(x,mean=parfit[["mean_nonresp"]],sd=exp(parfit[["b0"]]))
    p_nonresp_belowloq=(1-expit(beta_i))*pnorm(loq_n,mean=parfit[["mean_nonresp"]],sd=exp(parfit[["b0"]]))
  } else if (model_nonresponse==1) {
    y.seroneg = sapply(expit(beta_i),function(p) (1-p)*parfit[["rate_llog"]]*(1+exp(parfit[["shape_llog"]]))*
                         (parfit[["rate_llog"]]*x)^(1+exp(parfit[["shape_llog"]])-1)/(1+(parfit[["rate_llog"]]*x)^(1+exp(parfit[["shape_llog"]])))^2)
    p_nonresp_belowloq = (1-expit(beta_i))*1/(1+(parfit[["rate_llog"]]*loq_n)^(-exp(parfit[["shape_llog"]])))
  }
  
  if (loq_n>0) {
    y.seropos[x==loq_n]=p_resp_belowloq
    y.seroneg[x==loq_n]=p_nonresp_belowloq
  }
  
  dens = rowMeans(y.seropos+y.seroneg)
  denspos=rowMeans(y.seropos)
  densneg=rowMeans(y.seroneg)
  
  return(list(dens,denspos,densneg))
}

# Values of antibody to estimate pdf
x = seq(max(min(d_foranalysis[,y_name][is.finite(d_foranalysis[,y_name])],na.rm=T),loq_n),max(d_foranalysis[,y_name][is.finite(d_foranalysis[,y_name])],na.rm=T),by=0.01)

nsim=1000
sim_dens = matrix(NA,nrow=length(x),ncol=nsim)
sim_denspos = matrix(NA,nrow=length(x),ncol=nsim)
sim_densneg = matrix(NA,nrow=length(x),ncol=nsim)

for (i in 1:nsim) {
  
  par_row = pars[sample.int(nrow(pars),1),]
  
  r = overall_density(d_foranalysis,y_name,b_vars,g_vars,par_row,loq_n,x)
  sim_dens[,i] = r[[1]]
  sim_denspos[,i]=r[[2]]
  sim_densneg[,i]=r[[3]]
  
}

dens_summ = apply(sim_dens,1,function(x) quantile(x,c(0.025,0.5,0.975)))
denspos_summ = apply(sim_denspos,1,function(x) quantile(x,c(0.025,0.5,0.975)))
densneg_summ = apply(sim_densneg,1,function(x) quantile(x,c(0.025,0.5,0.975)))

d_foranalysis[,y_name][d_foranalysis[,y_name]<loq_n]=loq_n

# Observed density
dx = density(na.omit(d_foranalysis[,y_name]))$x
dy = density(na.omit(d_foranalysis[,y_name]))$y

# Basic plot of observed vs. model-estimated probability density function
plot(dx[dx>loq_n],dy[dx>loq_n],type='l',
     xlim=c(loq_n,max(d_foranalysis[,y_name])),
     xlab='Antibody (units)',
     ylab='Density')

# Add probability of being below LLOQ
if (loq_n>0) {
  points(loq_n,mean(d_foranalysis[,y_name]==loq_n),pch=16)
}
lines(x[x>=loq_n],dens_summ[2,x>=loq_n],col='blue')
lines(x[x>=loq_n],dens_summ[1,x>=loq_n],col='blue',lty=2)
lines(x[x>=loq_n],dens_summ[3,x>=loq_n],col='blue',lty=2)
points(x[x<loq_n],dens_summ[2,x<loq_n],col='blue',pch=16)
lines(rep(x[x<loq_n],2),c(dens_summ[1,x<loq_n],dens_summ[3,x<loq_n]),col='blue')

lines(x[x>=loq_n],denspos_summ[2,x>=loq_n],col='green')
lines(x[x>=loq_n],denspos_summ[1,x>=loq_n],col='green',lty=2)
lines(x[x>=loq_n],denspos_summ[3,x>=loq_n],col='green',lty=2)
points(x[x<loq_n],denspos_summ[2,x<loq_n],col='green',pch=16)
lines(rep(x[x<loq_n],2),c(denspos_summ[1,x<loq_n],denspos_summ[3,x<loq_n]),col='green')

lines(x[x>=loq_n],densneg_summ[2,x>=loq_n],col='red')
lines(x[x>=loq_n],densneg_summ[1,x>=loq_n],col='red',lty=2)
lines(x[x>=loq_n],densneg_summ[3,x>=loq_n],col='red',lty=2)
points(x[x<loq_n],densneg_summ[2,x<loq_n],col='red',pch=16)
lines(rep(x[x<loq_n],2),c(densneg_summ[1,x<loq_n],densneg_summ[3,x<loq_n]),col='red')





