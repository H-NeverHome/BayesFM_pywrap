#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 09:31:56 2021

@author: Hauke Niehaus
"""


''' 
### Befa_Wrappy
A handy python wrapper for the befa-package 

# Features currently NOT implemented in Wrapper
-> simul.dedic.facmod (Generate synthetic data from a dedicated factor model)
-> Any plotting in R backend
-> Covariate processing

# befa-package
Original Author: Rémi Piatek (remi.piatek at gmail.com)
Documentation: https://cran.r-project.org/web/packages/BayesFM/index.html

# Installation
-> Backend is R -> requires seperate installation
-> Befa_Wrappy works (currently) only on Linux due to rpy2 (see https://rpy2.github.io/)


# Minimum Dependencies for Python
-> rpy2 (https://rpy2.github.io/)
-> numpy (https://numpy.org/)
-> pandas (https://pandas.pydata.org/)
-> arviz
-> xarray

# For Backend/R dependencies see (https://cran.r-project.org/web/packages/BayesFM/index.html)

'''

#https://support.rstudio.com/hc/en-us/articles/219949047-Installing-older-versions-of-packages
#https://cran.r-project.org/src/contrib/Archive/BayesFM/



##### import func

from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro

##### helper functions

def r_to_pd(smth_r):
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    import rpy2.robjects as ro
    
    with localconverter(ro.default_converter + pandas2ri.converter):
        smth_pd   = ro.conversion.rpy2py(smth_r)
    return smth_pd

def pd_to_r(smth_pd):
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    import rpy2.robjects as ro
    
    with localconverter(ro.default_converter + pandas2ri.converter):
        smth_r   = ro.conversion.py2rpy(smth_pd)
    return smth_r






#befa_run(data,mcmc_iter,factor_spec, beta_prior, cov_mat_prior, dir_prior)
class befa:
    def __init__(self, mcmc_iter,factor_spec, beta_prior, cov_mat_prior, dir_prior,seed):
        ## imports
        import rpy2.robjects as ro
        from rpy2.robjects import globalenv
        from rpy2.robjects.packages import importr
        import pandas as pd
        
        base = importr('base')
        utils = importr('utils')
        utils.chooseCRANmirror(ind=1)
        ro.r('if (require(BayesFM) == FALSE) {install.packages("BayesFM")}')
        befa_base = importr('BayesFM')
        befa_base = importr('psych')
        
        ## clear global env
        ro.r('rm(list = ls())')
        ## set seed
        self.seed_curr = int(seed)
        globalenv['seed_curr'] = self.seed_curr
        ro.r('set.seed(seed_curr)')
        ## mcmc params
        self.mcmc_iter = mcmc_iter
        self.mcmc_burnin = (mcmc_iter/10)
        self.mcmc_effective_nsamples = self.mcmc_iter-self.mcmc_burnin
        ## factor specs
        self.Nid = factor_spec[0] # min load
        self.Kmax = factor_spec[1] # max n fac
        
        ## beta prior on probability τ0 
        # aka that a manifest variable does not load on any factor
        self.beta_1 = beta_prior[0]
        self.beta_2 = beta_prior[1]
        
        ## inv wishh prior
        self.covmatprior = cov_mat_prior
        if isinstance(cov_mat_prior,bool) == True:
            None
        else:
            self.kappa0 = cov_mat_prior[0] #nu0/dof_inv_w
            self.xi0 = cov_mat_prior[1] #S0 /scale_inv_w
        # fitted results
        self.posterior = ...
              
    def befa_nfac_prior_sim(self, nvar, kappa_range):
        # befa_run(int , (kappa_min,kappa_max) )    
        from rpy2.robjects import globalenv
        import rpy2.robjects as ro
        import pandas as pd
        
        ro.r('set.seed(seed_curr)')
        ## mcmc run
        globalenv['sim_mcmc_nrep'] = self.mcmc_iter
        ## nr of vars in dataset
        globalenv['sim_nvar'] = nvar
        ## range of simulated kappas
        globalenv['min_k'] = kappa_range[0]
        globalenv['max_k'] = kappa_range[1]
        ro.r('kappa_sim <- seq(min_k,max_k,.5)')        

        ## factor specs
        globalenv['sim_Kmax'] = float(self.Kmax) # max n fac
        globalenv['sim_nid'] = float(self.Nid) # max n fac
        
        #run sim
        ro.r('''prior_nfac_sim_res <- BayesFM::simul.nfac.prior(nrep = sim_mcmc_nrep,
                                                            nvar = sim_nvar,                                         
                                                            Kmax = sim_Kmax,
                                                            Nid  = sim_nid,
                                                            kappa = c(kappa_sim))''')
        ro.r('print(summary(prior_nfac_sim_res))')
        ro.r('res_sim_nfac <- summary(prior_nfac_sim_res)')
        ## collect simulated kappa values
        raw_ind = [i for i in globalenv['kappa_sim']]  
        ## convert R matrices to DS
        res_raw = [pd.DataFrame(i, index =raw_ind) for i in r_to_pd(globalenv['res_sim_nfac'])]    
        res_nfac_DF = res_raw[0].copy()
        ## set factor amount as columns
        res_nfac_DF.columns = [int(i)+1 for i in res_nfac_DF.columns.tolist()]
        ## include acceptance rate
        res_nfac_DF['accept'] = res_raw[1]
        return res_nfac_DF
 
       
    def befa_covmat_prior_sim(self):
        #befa_run(R.DF,int    ,(nid, kmax), (kappa0, xi0), (nu0,S0)/bool, float)    
        ## imports
        import rpy2.robjects as ro
        from rpy2.robjects import globalenv
        from rpy2.robjects.packages import importr
        import pandas as pd
        import numpy as np
        from befa_rpy_class_import import r_to_pd,pd_to_r
        ro.r('set.seed(seed_curr)')
        ## mcmc run
        globalenv['sim_mcmc_nrep'] = self.mcmc_iter
        
        ## factor specs
        globalenv['sim_Kmax'] = float(self.Kmax) # max n fac
        
        ## priors on cov mat of latent factors
        if isinstance(self.covmatprior,bool) == True:
            print('running befa w uninform HW.prior prior')
            ro.r('''prior_sim_res <- BayesFM::simul.R.prior(nrep = sim_mcmc_nrep,
                                                          Kmax = sim_Kmax,                                        
                                                          HW.prior = TRUE)''')
            ro.r('sim_sum_maxcorr <- summary(prior_sim_res)$maxcor')
            ro.r('sim_sum_mineig <- summary(prior_sim_res)$mineig')            
            ro.r('sim_sum_maxcorr_stat <- sim_sum_maxcorr$stat')
            ro.r('sim_sum_mineig_stat <- sim_sum_mineig$stat')
            ro.r('plot(prior_sim_res)')
            ro.r('summary(prior_sim_res)')
            
            
            # collect prior specifications aka rownames
            ro.r('sim_rwn <- rownames(sim_sum_mineig_stat)')
            rw_names = np.array(globalenv['sim_rwn'])
            
            # stat variables aka column names
            ro.r('sim_clmn <- colnames(sim_sum_mineig_stat)')
            cl_mineig_stat = np.array(globalenv['sim_clmn'])
            
            maxcorr_res = pd.DataFrame(r_to_pd(globalenv['sim_sum_maxcorr_stat']),
                                        index = rw_names, 
                                        columns = cl_mineig_stat)
            
            mineig_res = pd.DataFrame(r_to_pd(globalenv['sim_sum_mineig_stat']),
                                        index = rw_names, 
                                        columns = cl_mineig_stat)
            
            return {'min_eig': mineig_res,
                    'max_corr': maxcorr_res}
            
    
    
        elif isinstance(self.covmatprior,tuple) == True:
            print('running befa w inf inv wish prior')
            #nu0/dof_inv_w
            
            globalenv['min_nu0'] = 5
            globalenv['max_nu0'] = 32
            ro.r('kappa_sim_curr <- seq(min_nu0,max_nu0,3)')   
            ro.r('nu00 <- c(kappa_sim_curr)')
            
            #ro.r('nu00 <- sim_Kmax + c(2, 5, 6, 7, 8, 9, 10, 17)')
            
            #S0 /scale_inv_w
            globalenv['sim_scale_inv_w']    = self.xi0
            ro.r('''prior_sim_res <- BayesFM::simul.R.prior(nrep = sim_mcmc_nrep,
                                                            Kmax = sim_Kmax,                                         
                                                            nu0 = nu00,
                                                            S0 = sim_scale_inv_w)''')
            ro.r('sim_sum_maxcorr <- summary(prior_sim_res)$maxcor')
            ro.r('sim_sum_mineig <- summary(prior_sim_res)$mineig')
            ro.r('sim_sum_maxcorr_stat <- sim_sum_maxcorr$stat')
            ro.r('sim_sum_mineig_stat <- sim_sum_mineig$stat')
            
            # collect prior specifications aka rownames
            ro.r('sim_rwn <- rownames(sim_sum_mineig_stat)')
            rw_names = list(globalenv['sim_rwn'])
            
            # stat variables aka column names
            ro.r('sim_clmn <- colnames(sim_sum_mineig_stat)')
            cl_mineig_stat = list(globalenv['sim_clmn'])
            
            maxcorr_res = pd.DataFrame(r_to_pd(globalenv['sim_sum_maxcorr_stat']),
                                        index = rw_names, 
                                        columns = cl_mineig_stat)
            
            mineig_res = pd.DataFrame(r_to_pd(globalenv['sim_sum_mineig_stat']),
                                        index = rw_names, 
                                        columns = cl_mineig_stat)
            ro.r('plot(prior_sim_res)')
            return {'min_eig': mineig_res,
                    'max_corr': maxcorr_res} 
        
 
    def make_inference_dat(self,raw_dat):
        import pandas as pd
        import numpy as np
        import xarray as xr
        import arviz as az
        nsamples = int(self.mcmc_effective_nsamples)
        results = self.posterior
        df = pd.DataFrame()
        df["chain"] = [0 for i in range(nsamples)]
        df["draw"] = np.arange(nsamples, dtype=int)
        df = df.set_index(["chain", "draw"])
        names_all = list(results.keys())

        for key in names_all: 
            dat_curr = df.copy()
            for i in results[key]:
                df[i] = list(results[key][i])

        xdata_curr = xr.Dataset.from_dataframe(df)
        dataset_curr = az.InferenceData(posterior=xdata_curr)
        return dataset_curr
        
    
    def fit_uninform(self):
        import rpy2.robjects as ro
        from rpy2.robjects import globalenv
        from rpy2.robjects.packages import importr
        import pandas as pd
        print('running befa w uninform HW.prior prior')
        ro.r('''befa_RES <- BayesFM::befa(model = data_curr,
                                              burnin = mcmc_burnin,
                                              iter = mcmc_niter,
                                              Nid = Nid,
                                              Kmax = Kmax,
                                              kappa0 = beta_1,
                                              xi0 = beta_2,                                         
                                              HW.prior = TRUE)''')
        ro.r('befa_RES <- BayesFM::post.column.switch(befa_RES)')
        ro.r('befa_RES <- BayesFM::post.sign.switch(befa_RES)')

        ### get summary
        #ro.r('''print(summary(befa_RES, what='all'))''')
        ro.r('''res_all_uninf <- summary(befa_RES,what = 'maxp')''')
        
        # factor loadings
        ro.r('''res_all_alpha_uninf <- res_all_uninf$alpha''')
        #ro.r('''print(summary(befa_RES,what = 'maxp')$alpha)''')
        fac_load = pd.DataFrame(r_to_pd(globalenv['res_all_alpha_uninf']))
        #print(fac_load)
        
        #  Idiosyncratic variances
        ro.r('''res_all_sigma_uninf<- res_all_uninf$sigma''')
        #ro.r('''print(res_all_sigma)''')
        ideo_var = pd.DataFrame(r_to_pd(globalenv['res_all_sigma_uninf']))
        
        #  Correlation matrix of the latent factors 
        ro.r('''res_all_R_uninf <- res_all_uninf$R''')
        #ro.r('''print(res_all_R)''')
        fac_corrmat = pd.DataFrame(r_to_pd(globalenv['res_all_R_uninf']))
        
        # factor solution
        ro.r('''res_all_nfac_uninf <- res_all_uninf$nfac''')
        #ro.r('''print(res_all_nfac)''')
        ro.r('rw_nfac_uninf <- rownames(res_all_nfac_uninf)')
        ro.r('print(rw_nfac_uninf)')
        
        n_factors = len(r_to_pd(globalenv['res_all_nfac_uninf']))
        if n_factors > 1:
            nfac_solution = pd.DataFrame(r_to_pd(globalenv['res_all_nfac_uninf']),
                                index = list(r_to_pd(globalenv['rw_nfac_uninf'])))
            nfac_solution['nfac'] = list(nfac_solution.index)
            nfac_solution.columns = ['prob','nfac']
        else:
            fac_load['dedic'].max()
            nfac_solution = pd.DataFrame(index = [fac_load['dedic'].max()])
            nfac_solution['prob'] = r_to_pd(globalenv['res_all_nfac_uninf'])
            nfac_solution['nfac'] = [fac_load['dedic'].max()]
            nfac_solution.columns = ['prob','nfac']
        
        
        res_total_uninf = {'fac_load':fac_load,
                      'ideo_var':ideo_var,
                      'fac_corrmat':fac_corrmat,
                      'nfac_solution':nfac_solution}
        
        # get all keys in results
        ro.r('''names_befa_RES_uninf <- names(befa_RES)''')
        names_all = [ i for i in globalenv['names_befa_RES_uninf']]
        dict_res = {}
        for var in names_all:
            # string literals// f-function
            ro.r(f"curr_dat_uninf <- befa_RES${var}")
            ro.r(f"curr_col_uninf <- colnames(befa_RES${var})")
            if var not in ['nfac','MHacc']:
                colmns = [i for i in globalenv['curr_col_uninf']]
                dict_res[var] = pd.DataFrame(r_to_pd(globalenv['curr_dat_uninf']),
                                              columns = colmns).astype(float)
            
            else:
                colmns = [var]
                dict_res[var] = pd.DataFrame(r_to_pd(globalenv['curr_dat_uninf']),
                                              columns = colmns).astype(float)

        
        self.posterior = dict_res
        inf_dat_az = self.make_inference_dat(dict_res)
        
        #self.posterior_inf_data = inf_dat_az  
        return {'results_summary':res_total_uninf,
                'results_raw_mcmc_DF': dict_res,
                'results_inference_dat':inf_dat_az,
                'variables': names_all}         
        
        
        
        self.posterior = res_total_uninf
    
    
    def fit_inform(self):
        import rpy2.robjects as ro
        from rpy2.robjects import globalenv
        from rpy2.robjects.packages import importr
        import pandas as pd                
        print('running befa w inf inv wish prior')
        #nu0/dof_inv_w
        globalenv['dof_inv_w']      = self.kappa0 
        #S0 /scale_inv_w
        globalenv['scale_inv_w']    = self.xi0
        ro.r('''befa_RES_inf <- BayesFM::befa(model = data_curr,
                                              burnin = mcmc_burnin,
                                              iter = mcmc_niter,
                                              Nid = Nid,
                                              Kmax = Kmax,
                                              kappa0 = beta_1,
                                              xi0 = beta_2,                                         
                                              nu0 = dof_inv_w,
                                              S0 = scale_inv_w)''')
        
        ro.r('befa_RES_inf <- BayesFM::post.column.switch(befa_RES_inf)')
        ro.r('befa_RES_inf <- BayesFM::post.sign.switch(befa_RES_inf)')
        #ro.r('''plot(summary(befa_RES_inf, what='hppm'))''')
        ### get summary
        #ro.r('''print(summary(befa_RES_inf, what='all'))''')
        
        ro.r('''res_all_inf <- summary(befa_RES_inf, what = 'all')''')
        
        # factor loadings
        ro.r('''res_all_alpha_inf <- res_all_inf$alpha''')
        #ro.r('''print(res_all_alpha)''')
        fac_load = pd.DataFrame(r_to_pd(globalenv['res_all_alpha_inf']))
        
        #  Idiosyncratic variances
        ro.r('''res_all_sigma_inf <- res_all_inf$sigma''')
        #ro.r('''print(res_all_sigma)''')
        ideo_var = pd.DataFrame(r_to_pd(globalenv['res_all_sigma_inf']))
        
        #  Correlation matrix of the latent factors 
        ro.r('''res_all_R_inf <- res_all_inf$R''')
        #ro.r('''print(res_all_R)''')
        fac_corrmat = pd.DataFrame(r_to_pd(globalenv['res_all_R_inf']))
        
        # factor solution
        ro.r('''res_all_nfac_inf <- res_all_inf$nfac''')
        #ro.r('''print(res_all_nfac)''')
        ro.r('rw_nfac <- rownames(res_all_nfac_inf)')
        ro.r('print(rw_nfac)')
        
        n_factors = len(r_to_pd(globalenv['res_all_nfac_inf']))
        if n_factors > 1:
            nfac_solution = pd.DataFrame(r_to_pd(globalenv['res_all_nfac_inf']),
                                index = list(r_to_pd(globalenv['rw_nfac'])))
            nfac_solution['nfac'] = list(nfac_solution.index)
            nfac_solution.columns = ['prob','nfac']
        else:
            fac_load['dedic'].max()
            nfac_solution = pd.DataFrame(index = [fac_load['dedic'].max()])
            nfac_solution['prob'] = r_to_pd(globalenv['res_all_nfac_inf'])
            nfac_solution['nfac'] = [fac_load['dedic'].max()]
            nfac_solution.columns = ['prob','nfac']
        
        
        res_total_inf = {'fac_load':fac_load,
                      'ideo_var':ideo_var,
                      'fac_corrmat':fac_corrmat,
                      'nfac_solution':nfac_solution}   

        ro.r('''names_befa_RES_inf <- names(befa_RES_inf)''')
        names_all = [ i for i in globalenv['names_befa_RES_inf']]
        dict_res = {}
        for var in names_all:
            # string literals// f-function
            ro.r(f"curr_dat_inf <- befa_RES_inf${var}")
            ro.r(f"curr_col_inf <- colnames(befa_RES_inf${var})")
            if var not in ['nfac','MHacc']:
                colmns = [i for i in globalenv['curr_col_inf']]
                dict_res[var] = pd.DataFrame(r_to_pd(globalenv['curr_dat_inf']),
                                              columns = colmns).astype(float)
            
            else:
                colmns = [var]
                dict_res[var] = pd.DataFrame(r_to_pd(globalenv['curr_dat_inf']),
                                              columns = colmns).astype(float)
        
        
        self.posterior = dict_res
        inf_dat_az = self.make_inference_dat(dict_res)    
        self.posterior_inf_data = inf_dat_az                     
        return {'results_summary':res_total_inf,
                'results_raw_mcmc_DF': dict_res,
                'results_inference_dat':inf_dat_az,
                'variables': names_all}    



    def befa_fit_new(self, data): 
       #befa_run(R.DF,int    ,(nid, kmax), (kappa0, xi0), (nu0,S0)/bool, float)    
       import rpy2.robjects as ro
       from rpy2.robjects import globalenv
       from rpy2.robjects.packages import importr
       import pandas as pd
       ro.r('set.seed(seed_curr)')
       ## mcmc run
       globalenv['mcmc_burnin'] = self.mcmc_burnin
       globalenv['mcmc_niter'] = self.mcmc_iter-(self.mcmc_iter/10)
       
       ## factor specs
       globalenv['Nid'] = float(self.Nid) # min load
       globalenv['Kmax'] = float(self.Kmax) # max n fac
       
       ## beta prior on probability τ0 
       # aka that a manifest variable does not load on any factor
       globalenv['beta_1'] = self.beta_1 
       globalenv['beta_2'] = self.beta_2
       
       # (converted) data
       globalenv['data_curr'] = data
       
       ## priors on cov mat of latent factors
       if isinstance(self.covmatprior,bool) == True:
           res_1 = self.fit_uninform()
           return res_1
       elif isinstance(self.covmatprior,tuple) == True:           
           res_2 = self.fit_inform()
           return res_2
           



############## OLD


    
    # def befa_fit(self, data):
    #     #befa_run(R.DF,int    ,(nid, kmax), (kappa0, xi0), (nu0,S0)/bool, float)    
    #     import rpy2.robjects as ro
    #     from rpy2.robjects import globalenv
    #     from rpy2.robjects.packages import importr
    #     import pandas as pd
    #     ro.r('set.seed(seed_curr)')
    #     ## mcmc run
    #     globalenv['mcmc_burnin'] = self.mcmc_burnin
    #     globalenv['mcmc_niter'] = self.mcmc_iter-(self.mcmc_iter/10)
        
    #     ## factor specs
    #     globalenv['Nid'] = float(self.Nid) # min load
    #     globalenv['Kmax'] = float(self.Kmax) # max n fac
        
    #     ## beta prior on probability τ0 
    #     # aka that a manifest variable does not load on any factor
    #     globalenv['beta_1'] = self.beta_1 
    #     globalenv['beta_2'] = self.beta_2
        
    #     # (converted) data
    #     globalenv['data_curr'] = data
        
    #     ## priors on cov mat of latent factors
    #     if isinstance(self.covmatprior,bool) == True:
    #         print('running befa w uninform HW.prior prior')
    #         ro.r('''befa_RES <- BayesFM::befa(model = data_curr,
    #                                               burnin = mcmc_burnin,
    #                                               iter = mcmc_niter,
    #                                               Nid = Nid,
    #                                               Kmax = Kmax,
    #                                               kappa0 = beta_1,
    #                                               xi0 = beta_2,                                         
    #                                               HW.prior = TRUE)''')
    #         ro.r('befa_RES <- BayesFM::post.column.switch(befa_RES)')
    #         ro.r('befa_RES <- BayesFM::post.sign.switch(befa_RES)')

    #         ### get summary
    #         ro.r('''print(summary(befa_RES, what='all'))''')
    #         ro.r('''res_all_uninf <- summary(befa_RES,what = 'all')''')
            
    #         # factor loadings
    #         ro.r('''res_all_alpha_uninf <- res_all_uninf$alpha''')
    #         #ro.r('''print(res_all_alpha)''')
    #         fac_load = pd.DataFrame(r_to_pd(globalenv['res_all_alpha_uninf']))
            
    #         #  Idiosyncratic variances
    #         ro.r('''res_all_sigma_uninf<- res_all_uninf$sigma''')
    #         #ro.r('''print(res_all_sigma)''')
    #         ideo_var = pd.DataFrame(r_to_pd(globalenv['res_all_sigma_uninf']))
            
    #         #  Correlation matrix of the latent factors 
    #         ro.r('''res_all_R_uninf <- res_all_uninf$R''')
    #         #ro.r('''print(res_all_R)''')
    #         fac_corrmat = pd.DataFrame(r_to_pd(globalenv['res_all_R_uninf']))
            
    #         # factor solution
    #         ro.r('''res_all_nfac_uninf <- res_all_uninf$nfac''')
    #         #ro.r('''print(res_all_nfac)''')
    #         ro.r('rw_nfac_uninf <- rownames(res_all_nfac_uninf)')
    #         ro.r('print(rw_nfac_uninf)')
            
    #         n_factors = len(r_to_pd(globalenv['res_all_nfac_uninf']))
    #         if n_factors > 1:
    #             nfac_solution = pd.DataFrame(r_to_pd(globalenv['res_all_nfac_uninf']),
    #                                 index = list(r_to_pd(globalenv['rw_nfac_uninf'])))
    #             nfac_solution['nfac'] = list(nfac_solution.index)
    #             nfac_solution.columns = ['prob','nfac']
    #         else:
    #             fac_load['dedic'].max()
    #             nfac_solution = pd.DataFrame(index = [fac_load['dedic'].max()])
    #             nfac_solution['prob'] = r_to_pd(globalenv['res_all_nfac_uninf'])
    #             nfac_solution['nfac'] = [fac_load['dedic'].max()]
    #             nfac_solution.columns = ['prob','nfac']
            
            
    #         res_total_uninf = {'fac_load':fac_load,
    #                       'ideo_var':ideo_var,
    #                       'fac_corrmat':fac_corrmat,
    #                       'nfac_solution':nfac_solution} 



    #     elif isinstance(self.covmatprior,tuple) == True:
    #         print('running befa w inf inv wish prior')
    #         #nu0/dof_inv_w
    #         globalenv['dof_inv_w']      = self.kappa0 
    #         #S0 /scale_inv_w
    #         globalenv['scale_inv_w']    = self.xi0
    #         ro.r('''befa_RES_inf <- BayesFM::befa(model = data_curr,
    #                                               burnin = mcmc_burnin,
    #                                               iter = mcmc_niter,
    #                                               Nid = Nid,
    #                                               Kmax = Kmax,
    #                                               kappa0 = beta_1,
    #                                               xi0 = beta_2,                                         
    #                                               nu0 = dof_inv_w,
    #                                               S0 = scale_inv_w)''')
    
    #         ro.r('befa_RES_inf <- BayesFM::post.column.switch(befa_RES_inf)')
    #         ro.r('befa_RES_inf <- BayesFM::post.sign.switch(befa_RES_inf)')
    #         #ro.r('''plot(summary(befa_RES_inf, what='hppm'))''')
    #         ### get summary
    #         ro.r('''print(summary(befa_RES_inf, what='all'))''')
            
    #         ro.r('''res_all_inf <- summary(befa_RES_inf, what = 'all')''')
            
    #         # factor loadings
    #         ro.r('''res_all_alpha_inf <- res_all_inf$alpha''')
    #         #ro.r('''print(res_all_alpha)''')
    #         fac_load = pd.DataFrame(r_to_pd(globalenv['res_all_alpha_inf']))
            
    #         #  Idiosyncratic variances
    #         ro.r('''res_all_sigma_inf <- res_all_inf$sigma''')
    #         #ro.r('''print(res_all_sigma)''')
    #         ideo_var = pd.DataFrame(r_to_pd(globalenv['res_all_sigma_inf']))
            
    #         #  Correlation matrix of the latent factors 
    #         ro.r('''res_all_R_inf <- res_all_inf$R''')
    #         #ro.r('''print(res_all_R)''')
    #         fac_corrmat = pd.DataFrame(r_to_pd(globalenv['res_all_R_inf']))
            
    #         # factor solution
    #         ro.r('''res_all_nfac_inf <- res_all_inf$nfac''')
    #         #ro.r('''print(res_all_nfac)''')
    #         ro.r('rw_nfac <- rownames(res_all_nfac_inf)')
    #         ro.r('print(rw_nfac)')
            
    #         n_factors = len(r_to_pd(globalenv['res_all_nfac_inf']))
    #         if n_factors > 1:
    #             nfac_solution = pd.DataFrame(r_to_pd(globalenv['res_all_nfac_inf']),
    #                                 index = list(r_to_pd(globalenv['rw_nfac'])))
    #             nfac_solution['nfac'] = list(nfac_solution.index)
    #             nfac_solution.columns = ['prob','nfac']
    #         else:
    #             fac_load['dedic'].max()
    #             nfac_solution = pd.DataFrame(index = [fac_load['dedic'].max()])
    #             nfac_solution['prob'] = r_to_pd(globalenv['res_all_nfac_inf'])
    #             nfac_solution['nfac'] = [fac_load['dedic'].max()]
    #             nfac_solution.columns = ['prob','nfac']
            
            
    #         res_total_inf = {'fac_load':fac_load,
    #                       'ideo_var':ideo_var,
    #                       'fac_corrmat':fac_corrmat,
    #                       'nfac_solution':nfac_solution} 
            
            
            

    #     if isinstance(self.covmatprior,bool) == True:
    #         # get all keys in results
    #         ro.r('''names_befa_RES_uninf <- names(befa_RES)''')
    #         names_all = [ i for i in globalenv['names_befa_RES_uninf']]
    #         dict_res = {}
    #         for var in names_all:
    #             # string literals// f-function
    #             ro.r(f"curr_dat_uninf <- befa_RES${var}")
    #             ro.r(f"curr_col_uninf <- colnames(befa_RES${var})")
    #             if var not in ['nfac','MHacc']:
    #                 colmns = [i for i in globalenv['curr_col_uninf']]
    #                 dict_res[var] = pd.DataFrame(r_to_pd(globalenv['curr_dat_uninf']),
    #                                               columns = colmns).astype(float)
                
    #             else:
    #                 colmns = [var]
    #                 dict_res[var] = pd.DataFrame(r_to_pd(globalenv['curr_dat_uninf']),
    #                                               columns = colmns).astype(float)

    #         inf_dat_az = self.make_inference_dat(dict_res)
    #         self.posterior = dict_res
    #         self.posterior_inf_data = inf_dat_az  
    #         return {'results_summary':res_total_uninf,
    #                 'results_raw_mcmc_DF': dict_res,
    #                 'results_inference_dat':inf_dat_az,
    #                 'variables': list(dict_res.keys())}       

            
    #     elif isinstance(self.covmatprior,tuple) == True:
    #         # get all keys in results
    #         ro.r('''names_befa_RES_inf <- names(befa_RES_inf)''')
    #         names_all = [ i for i in globalenv['names_befa_RES_inf']]
    #         dict_res = {}
    #         for var in names_all:
    #             # string literals// f-function
    #             ro.r(f"curr_dat_inf <- befa_RES_inf${var}")
    #             ro.r(f"curr_col_inf <- colnames(befa_RES_inf${var})")
    #             if var not in ['nfac','MHacc']:
    #                 colmns = [i for i in globalenv['curr_col_inf']]
    #                 dict_res[var] = pd.DataFrame(r_to_pd(globalenv['curr_dat_inf']),
    #                                               columns = colmns).astype(float)
                
    #             else:
    #                 colmns = [var]
    #                 dict_res[var] = pd.DataFrame(r_to_pd(globalenv['curr_dat_inf']),
    #                                               columns = colmns).astype(float)

    #         inf_dat_az = self.make_inference_dat(dict_res)    
    #         self.posterior = dict_res
    #         self.posterior_inf_data = inf_dat_az                     
    #         return {'results_summary':res_total_inf,
    #                 'results_raw_mcmc_DF': dict_res,
    #                 'results_inference_dat':inf_dat_az,
    #                 'variables': list(dict_res.keys())}  
 

    # def plot_posterior(self,what):
    #     import arviz as az
    #     import pandas as pd
    #     for var in what:
    #         if var in ['nfac','dedic']:
    #             az.plot_posterior(self.posterior_inf_data,
    #                               var_names = [i for i in self.posterior[var]],
    #                               kind='hist',
    #                               point_estimate=None,
    #                               hdi_prob='hide')            
    #         # else:
    #         #     az.plot_posterior(self.posterior,
    #         #                       var_names = [i for i in df if'alpha' in i],
    #         #                       point_estimate='mean')
                
            
        
        
         
            
         
            
         
            
         
            
         
            
         
            
         
            
         
            
         
            
         
            
         