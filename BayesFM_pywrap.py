

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






#befa_run(data,mcmc_iter,factor_spec, beta_prior, cov_mat_prior, dir_prior, savepath, rnd_seed)
class befa:
    def __init__(self, mcmc_iter,factor_spec, beta_prior, cov_mat_prior, dir_prior,savepath,seed):
        ## imports
        import rpy2.robjects as ro
        from rpy2.robjects import globalenv
        from rpy2.robjects.packages import importr
        import pandas as pd
        # set path to output folder != working directory
        self.savepath = savepath
        
        
        base = importr('base')
        utils = importr('utils')
        utils.chooseCRANmirror(ind=1)
        ro.r('if (require(BayesFM) == FALSE) {install.packages("BayesFM")}')
        befa_base = importr('BayesFM')
        befa_base = importr('psych')
        
        ## clear global env
        ro.r('rm(list = ls())')
        
        ### safe parameters
        params = {}
        
        ## set seed
        self.seed = int(seed)

        
        globalenv['seed_curr'] = self.seed

        ro.r('set.seed(seed_curr)')
        
        params['seed'] = self.seed

        ## mcmc params
        self.mcmc_iter = mcmc_iter
        self.mcmc_burnin = (mcmc_iter/10)
        self.mcmc_effective_nsamples = self.mcmc_iter-self.mcmc_burnin
        
        params['mcmc_iter'] = mcmc_iter
        params['mcmc_burnin'] = (mcmc_iter/10)
        
        ## factor specs
        self.Nid = factor_spec[0] # min load
        self.Kmax = factor_spec[1] # max n fac
        
        params['Nid'] = factor_spec[0]
        params['Kmax'] = factor_spec[1]
        
        ## beta prior on probability τ0 
        # aka that a manifest variable does not load on any factor
        self.beta_1 = beta_prior[0]
        self.beta_2 = beta_prior[1]
        
        params['kappa0']= beta_prior[0]
        params['xi0']   = beta_prior[1]
        
        ## Concentration parameters of the Dirichlet prior distribution on the indicators.
        params['kappa'] = dir_prior
        self.kappa = dir_prior
        
        ## inv wishh prior on covariance matrix
        self.covmatprior = cov_mat_prior
        if isinstance(cov_mat_prior,bool) == True:
            params['nu0_S0'] = 'HW_prior' 
        else:
            self.nu0 = cov_mat_prior[0] #nu0/dof_inv_w
            self.S0 = cov_mat_prior[1] #S0 /scale_inv_w
            params['nu0'] = cov_mat_prior[0]
            params['S0'] = cov_mat_prior[1]
        
        ## save parameters
        self.params = params
        
        
        # fitted results
        self.posterior = ...
        self.summary_fac_load = ...


 
             
    def befa_nfac_prior_sim(self, nvar, kappa_range):
        # befa_run(int , (kappa_min,kappa_max) )    
        from rpy2.robjects import globalenv
        import rpy2.robjects as ro
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        ro.r('set.seed(seed_curr)')
        ## mcmc run
        globalenv['sim_mcmc_nrep'] = 25000
        ## nr of vars in dataset
        globalenv['sim_nvar'] = nvar
        
        # input kappas direct
        kappa_series = pd.Series(kappa_range)
        globalenv['kappa_sim'] = pd_to_r(kappa_series)
      
        ## factor specs
        globalenv['sim_Kmax'] = float(self.Kmax) # max n fac
        globalenv['sim_nid'] = float(self.Nid) # max n fac
        
        #run sim
        ro.r('''prior_nfac_sim_res <- BayesFM::simul.nfac.prior(nrep = sim_mcmc_nrep,
                                                            nvar = sim_nvar,                                         
                                                            Kmax = sim_Kmax,
                                                            Nid  = sim_nid,
                                                            kappa = c(kappa_sim))''')
        #ro.r('print(summary(prior_nfac_sim_res))')
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
        self.nfac_sim = res_nfac_DF
        
        res_nfac_DF.to_html(self.savepath + '/nfac_sim.html')
        
        plt.clf()
        plt.figure(figsize=(30,15))
        hm = sns.heatmap(res_nfac_DF.round(2),
                      annot = True,
                      square = True,
                      center=0.,
                      cmap="coolwarm",
                      linewidths = .9)
        hm.set_yticklabels(hm.get_yticklabels(),
                          rotation=0)

        path_hm = self.savepath + "/nfac_sim.png"
        plt.savefig(str(path_hm),
                    dpi=500,
                    format = 'png')
        
        return res_nfac_DF
 


       
    def befa_covmat_prior_sim(self):
        #befa_run(R.DF,int    ,(nid, kmax), (kappa0, xi0), (nu0,S0)/bool, float)    
        ## imports
        import rpy2.robjects as ro
        from rpy2.robjects import globalenv
        from rpy2.robjects.packages import importr
        import pandas as pd
        import numpy as np
        
        #from befa_rpy_class_import import r_to_pd,pd_to_r
        
        ro.r('set.seed(seed_curr)')
        ## mcmc run/ only 50k -> >100k gives weird error
        globalenv['sim_mcmc_nrep'] = 25000
        
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
            res_sim_covmat = {'min_eig': mineig_res,
                    'max_corr': maxcorr_res}
            self.covmat_sim = res_sim_covmat
            
            maxcorr_res.to_html(self.savepath + '/maxcorr_sim.html')
            mineig_res.to_html(self.savepath + '/mineig_sim.html')
            return res_sim_covmat
            
    
    
        elif isinstance(self.covmatprior,tuple) == True:
            print('running befa w inf inv wish prior')
            #nu0/dof_inv_w
            
            globalenv['min_nu0'] = 10
            globalenv['max_nu0'] = 32
            ro.r('kappa_sim_curr <- seq(min_nu0,max_nu0,3)')   
            ro.r('nu00 <- c(kappa_sim_curr)')
            

            #S0 /scale_inv_w
            globalenv['sim_scale_inv_w']    = self.S0
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
            res_sim_covmat = {'min_eig': mineig_res,
                              'max_corr': maxcorr_res}
            self.covmat_sim = res_sim_covmat
            maxcorr_res.to_html(self.savepath + '/maxcorr_sim.html')
            mineig_res.to_html(self.savepath + '/mineig_sim.html')
            return res_sim_covmat



    def inf_data(self,raw_dat):
        import pandas as pd
        import numpy as np
        import xarray as xr
        import arviz as az
        nsamples = int(self.mcmc_effective_nsamples)
        results = self.posterior
        df_mcmc = pd.DataFrame()
        df_mcmc["chain"] = [0 for i in range(nsamples)]
        df_mcmc["draw"] = np.arange(nsamples, dtype=int)
        df_mcmc = df_mcmc.set_index(["chain", "draw"])
        names_all = list(results.keys())

        # for each sampled variable
        res_inf_dat = {}

        for key in names_all: 
            dat_curr = df_mcmc.copy()
            curr_res = results[key]
            curr_vars = list(results[key].columns)
            for i in curr_vars:
                dat_curr[i] = list(curr_res[i].copy())
                
            xdata_curr = xr.Dataset.from_dataframe(dat_curr)
            dataset_curr = az.InferenceData(posterior=xdata_curr)
            res_inf_dat[key] = dataset_curr.copy()
        
        return res_inf_dat
    
 
       
    
    def fit_uninform(self):
        # Import
        import rpy2.robjects as ro
        from rpy2.robjects import globalenv
        from rpy2.robjects.packages import importr
        import pandas as pd
        
        # set seed 
        globalenv['seed_curr'] = self.seed 
        ro.r('set.seed(seed_curr)')  
        
        # fit befa model
        print('running befa w uninform HW.prior prior')
        ro.r('''befa_RES <- BayesFM::befa(model = data_curr,
                                              burnin = mcmc_burnin,
                                              iter = mcmc_niter,
                                              Nid = Nid,
                                              Kmax = Kmax,
                                              kappa = kappa,
                                              kappa0 = beta_1,
                                              xi0 = beta_2,                                         
                                              HW.prior = TRUE)''')
        
        ## Column and sign switching issues
        ro.r('befa_RES <- BayesFM::post.column.switch(befa_RES)')
        ro.r('befa_RES <- BayesFM::post.sign.switch(befa_RES)')

        ### get summary
        ro.r('''res_all_uninf <- summary(befa_RES,what = 'maxp')''')
        
        # factor loadings
        ro.r('''res_all_alpha_uninf <- res_all_uninf$alpha''')
        fac_load = pd.DataFrame(r_to_pd(globalenv['res_all_alpha_uninf']))

        
        #  Idiosyncratic variances
        ro.r('''res_all_sigma_uninf<- res_all_uninf$sigma''')
        ideo_var = pd.DataFrame(r_to_pd(globalenv['res_all_sigma_uninf']))
        
        #  Correlation matrix of the latent factors 
        ro.r('''res_all_R_uninf <- res_all_uninf$R''')
        fac_corrmat = pd.DataFrame(r_to_pd(globalenv['res_all_R_uninf']))
        
        # factor solution
        ro.r('''res_all_nfac_uninf <- res_all_uninf$nfac''')
        ro.r('rw_nfac_uninf <- rownames(res_all_nfac_uninf)')
        #ro.r('print(rw_nfac_uninf)')
        
        ## get and organize factor solution & probs
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
        self.summary_fac_load =    res_total_uninf
        
        # get all keys in results
        ro.r('''names_befa_RES_uninf <- names(befa_RES)''')
        names_all = [ i for i in globalenv['names_befa_RES_uninf']]
        dict_res = {}
        
        for var in names_all:
            # string literals// f-function
            ro.r(f"curr_dat_uninf <- befa_RES${var}")
            ro.r(f"curr_col_uninf <- colnames(befa_RES${var})")
            
            # seperate discrete and continous
            if var not in ['nfac','MHacc']:
                colmns = [i for i in globalenv['curr_col_uninf']]
                dict_res[var] = pd.DataFrame(r_to_pd(globalenv['curr_dat_uninf']),
                                              columns = colmns).astype(float)
            
            else:
                colmns = [var]
                dict_res[var] = pd.DataFrame(r_to_pd(globalenv['curr_dat_uninf']),
                                              columns = colmns).astype(float)

        # save posterior
        self.posterior = dict_res
        
        #convert to arviz & save
        inf_dat_az = self.inf_data(dict_res)
        self.posterior_inf_data = inf_dat_az  
        
        return {'results_summary':res_total_uninf,
                'results_raw_mcmc_DF': dict_res,
                'results_inference_dat':inf_dat_az,
                'variables': names_all,
                'seed' : ... }         
        
        


    
    def fit_inform(self):
        import rpy2.robjects as ro
        from rpy2.robjects import globalenv
        from rpy2.robjects.packages import importr
        import pandas as pd   
        from befa_rpy_class_import import r_to_pd,pd_to_r
                   
        print('running befa w inf inv wish prior')
        globalenv['seed_curr'] = self.seed 
        ro.r('set.seed(seed_curr)')  
        #nu0/dof_inv_w
        globalenv['dof_inv_w']      = self.nu0 
        
        #S0 /scale_inv_w
        globalenv['scale_inv_w']    = self.S0
        
        ## fit model
        ro.r('''befa_RES_inf <- BayesFM::befa(model = data_curr,
                                              burnin = mcmc_burnin,
                                              iter = mcmc_niter,
                                              Nid = Nid,
                                              Kmax = Kmax,
                                              kappa = kappa,
                                              kappa0 = beta_1,                                             
                                              xi0 = beta_2,                                         
                                              nu0 = dof_inv_w,
                                              S0 = scale_inv_w)''')
        
        ro.r('befa_RES_inf <- BayesFM::post.column.switch(befa_RES_inf)')
        ro.r('befa_RES_inf <- BayesFM::post.sign.switch(befa_RES_inf)')
        
        ro.r('''res_all_inf <- summary(befa_RES_inf, what = 'maxp')''')
        
        # factor loadings
        ro.r('''res_all_alpha_inf <- res_all_inf$alpha''')
        fac_load = pd.DataFrame(r_to_pd(globalenv['res_all_alpha_inf']))
        
        #  Idiosyncratic variances
        ro.r('''res_all_sigma_inf <- res_all_inf$sigma''')
        ideo_var = pd.DataFrame(r_to_pd(globalenv['res_all_sigma_inf']))
        
        #  Correlation matrix of the latent factors 
        ro.r('''res_all_R_inf <- res_all_inf$R''')
        fac_corrmat = pd.DataFrame(r_to_pd(globalenv['res_all_R_inf']))
        
        # factor solution
        ro.r('''res_all_nfac_inf <- res_all_inf$nfac''')
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
        self.summary_fac_load =    res_total_inf

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
        inf_dat_az = self.inf_data(dict_res)    
        self.posterior_inf_data = inf_dat_az
                  
        return {'results_summary':res_total_inf,
                'results_raw_mcmc_DF': dict_res,
                'results_inference_dat':inf_dat_az,
                'variables': names_all,
                'seed' : int(r_to_pd(globalenv['seed_curr'])) }    




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
        
        ## Kappa prior
        globalenv['kappa'] = self.kappa       
        # (converted) data
        globalenv['data_curr'] = data
       
        ## priors on cov mat of latent factors
        if isinstance(self.covmatprior,bool) == True:
            res_1 = self.fit_uninform()
            self.summary_fac_load = res_1['results_summary']['fac_load']
            self.pkldat = res_1
            safe_html = self.savepath
            for summ in res_1['results_summary']:
                res_1['results_summary'][summ].to_html(safe_html + '/' + summ + '.html' )
            return res_1
        elif isinstance(self.covmatprior,tuple) == True:           
            res_2 = self.fit_inform()
            self.summary_fac_load = res_2['results_summary']['fac_load']
            self.pkldat = res_2
            safe_html = self.savepath
            for summ in res_2['results_summary']:
                res_2['results_summary'][summ].to_html(safe_html + '/' + summ + '.html' )
            return res_2




    def plot_posterior(self, what):
        import arviz as az
        import matplotlib.pyplot as plt
        from pathlib import Path
        import os
        plt.ioff()
        # 'what' should contain a subset (list) of model variables like 'alpha' etc.
        # 'path' contains the (not relative) save location/folder path
        curr_inference_dat = self.posterior_inf_data
        res_path = Path(self.savepath+'/posterior_plots')
        print('Producing posterior plots')

        if res_path.exists() == True:
            None
        else:
            os.makedirs(res_path)
        
        for var in what:
            if var in ['dedic','nfac','MHacc']:
                print('plotting ' + f"{var}")
                kind = 'hist'
                hdi = 'hide'
                plot = az.plot_posterior(curr_inference_dat[var],
                                         kind=kind,
                                         figsize=(20,15),
                                         point_estimate = None)
                curr_path1 = res_path / f"{var}.png"
                plt.savefig(curr_path1,
                            dpi=500,
                            format = 'png')                
                
            if var=='alpha':
                print('plotting ' + f"{var}")
                #### maek extra directory for alpha
                alpha_path = Path(str(res_path) +'/'+'alpha')
                if alpha_path.exists() == True:
                    None
                else:
                    os.makedirs(alpha_path)
              
                curr_infdat = curr_inference_dat[var]                
                alpha_idx = list(curr_infdat.to_dict()['posterior'].keys())
                for alpha in alpha_idx:
                    plot_alpha = az.plot_posterior(curr_infdat,
                                                   var_names = alpha,
                                                   hdi_prob=.95,
                                                   rope=[-.1,.1],
                                                   ref_val = .0,
                                                   figsize=(20,15),
                                                   point_estimate='mean')               
                    file_alpha = str(alpha_path) +'/' +alpha+'.png'
                    plt.savefig(file_alpha,
                                dpi=500,
                                format = 'png')
                    plt.clf()  
                    
                    
            elif var == 'sigma':
                print('plotting ' + f"{var}")
                #### maek extra directory for alpha
                sigma_path = Path(str(res_path) +'/'+'sigma')
                if sigma_path.exists() == True:
                    None
                else:
                    os.makedirs(sigma_path)
              
                curr_infdat = curr_inference_dat[var]                
                sigma_idx = list(curr_infdat.to_dict()['posterior'].keys())
                for sigma in sigma_idx:
                    plot_alpha = az.plot_posterior(curr_infdat,
                                                   var_names = sigma,
                                                   hdi_prob=.95,
                                                   rope=[-.1,.1],
                                                   ref_val = .0,
                                                   figsize=(20,15),
                                                   point_estimate='mean')               
                    file_sigma = str(sigma_path) +'/' + sigma +'.png'
                    plt.savefig(file_sigma,
                                dpi=500,
                                format = 'png')
                    plt.clf()                 
                
                
            if var=='R':
                #### maek extra directory for alpha
                print('plotting ' + f"{var}")
                r_path = Path(str(res_path) +'/'+'R')
                if r_path.exists() == True:
                    None
                else:
                    os.makedirs(r_path)
              
                curr_infdat = curr_inference_dat[var]                
                r_idx = list(curr_infdat.to_dict()['posterior'].keys())
                for r in r_idx:
                    plot_r = az.plot_posterior(curr_infdat,
                                               var_names = r,
                                               hdi_prob=.95,
                                               rope=[-.1,.1],
                                               ref_val = .0,
                                               figsize=(20,15),
                                               point_estimate='mean')               
                    file_r = str(r_path) +'/' +r+'.png'
                    plt.savefig(file_r,
                                dpi=500,
                                format = 'png')
                    plt.clf()     
                    
        return curr_inference_dat
            
            


    def plot_traces(self,what):
        import arviz as az
        import matplotlib.pyplot as plt
        from pathlib import Path
        import os
        import matplotlib as mpl
        plt.ioff()
        
        # 'what' should contain a subset (list) of model variables like 'alpha' etc.
        # 'path' contains the (not relative) save location/folder path
        print('Producing traceplots')
        curr_inference_dat = self.posterior_inf_data
        res_path = Path(self.savepath+'/trace_plots')

        if res_path.exists() == True:
            None
        else:
            os.makedirs(res_path)
        
        for var in what:   

            if var=='alpha':
                print('plotting ' + f"{var}")
                #### maek extra directory for alpha
                alpha_path = Path(str(res_path) +'/'+'alpha')
                if alpha_path.exists() == True:
                    None
                else:
                    os.makedirs(alpha_path)
              
                curr_infdat = curr_inference_dat[var]                
                alpha_idx = list(curr_infdat.to_dict()['posterior'].keys())
                # for alpha in alpha_idx:
                    
                #     plot_alpha = az.plot_trace(curr_infdat,
                #                                var_names = alpha,
                #                                 figsize=(15,40),
                #                                 rug=True,
                #                                 combined = False,
                #                                 legend=True)               
                #     file_alpha = str(alpha_path) +'/' +alpha+'.png'
                #     plt.savefig(file_alpha,
                #                 dpi=500,
                #                 format = 'png')
                #     plt.clf() 
                plot_alpha = az.plot_trace(curr_infdat,
                                            figsize=(15,40),
                                            rug=False,
                                            legend=True)               
                file_alpha = str(alpha_path) +'/plot_alpha.png'
                plt.savefig(file_alpha,
                            dpi=250,
                            format = 'png')
                plt.clf()                     
                    
                    
            elif var == 'sigma':
                print('plotting ' + f"{var}")
                #### maek extra directory for alpha
                sigma_path = Path(str(res_path) +'/'+'sigma')
                if sigma_path.exists() == True:
                    None
                else:
                    os.makedirs(sigma_path)
              
                curr_infdat = curr_inference_dat[var]                
                sigma_idx = list(curr_infdat.to_dict()['posterior'].keys())
                # for sigma in sigma_idx:
                #     plot_sigma = az.plot_trace(curr_infdat,
                #                                var_names = sigma,
                #                                figsize=(20,15),
                #                                rug=True,
                #                                combined = False,
                #                                legend=True)              
                #     file_sigma = str(sigma_path) +'/' + sigma +'.png'
                #     plt.savefig(file_sigma,
                #                 dpi=500,
                #                 format = 'png')
                #     plt.clf()                 
                plt.clf()
                plot_sigma = az.plot_trace(curr_infdat,
                                           figsize=(15,40),
                                           rug=False,
                                           legend=True)              
                file_sigma = str(sigma_path) +'/plot_sigma.png'
                plt.savefig(file_sigma,
                            dpi=250,
                            format = 'png')
               
                
            elif var=='R':
                #### maek extra directory for alpha
                print('plotting ' + f"{var}")
                r_path = Path(str(res_path) +'/'+'R')
                if r_path.exists() == True:
                    None
                else:
                    os.makedirs(r_path)
              
                curr_infdat = curr_inference_dat[var]                
                r_idx = list(curr_infdat.to_dict()['posterior'].keys())
                # for r in r_idx:
                #     plot_r = az.plot_trace(curr_infdat,
                #                            var_names = r,
                #                            figsize=(20,15),
                #                            rug=True,
                #                            combined = False,
                #                            legend=True)              
                #     file_r = str(r_path) +'/' +r+'.png'
                #     plt.savefig(file_r,
                #                 dpi=500,
                #                 format = 'png')
                #     plt.clf()
                
                plt.clf()
                plot_r = az.plot_trace(curr_infdat,
                                       figsize=(15,40),
                                       rug=False,
                                       legend=True)              
                file_r = str(r_path) +'/plot_r.png'
                plt.savefig(file_r,
                            dpi=250,
                            format = 'png')

                
            else:
                print('thats all folks')
                    
        return curr_inference_dat

        


    def plot_factorload_summary(self):
        import seaborn as sns
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        from pathlib import Path
        import os      
        print('Producing factorload summary')
        res_path = Path(self.savepath+'/plots')

        if res_path.exists() == True:
            None
        else:
            os.makedirs(res_path)
            
        #sns.set(rc={'figure.figsize':(15,20)})
        what = ['dedic','mean']
        summ_res = self.summary_fac_load[what].copy()

        clms = list(summ_res.copy().T.columns)
        indx = np.unique([int(i) for i in summ_res['dedic']])
        indx_len = len(indx)
        plot_alpha = pd.DataFrame(index= indx,
                                  columns = list(summ_res.index))

        for clm in clms:
            curr_insert = pd.Series(index = indx, dtype=int)
            curr_A = summ_res.copy().T[clm].tolist()
            curr_insert.loc[int(curr_A[0])] = curr_A[1]
            plot_alpha[clm] = list(curr_insert)
        plt.clf()
        plt.figure(figsize=(30,15))
        hm = sns.heatmap(plot_alpha,
                      annot = True,
                      square = True,
                      center=0.,
                      cmap="coolwarm",
                      linewidths = .9,
                      xticklabels = [i for i in plot_alpha])
        hm.set_yticklabels(hm.get_yticklabels(),
                          rotation=0)

        path_hm = res_path / "factorload.png"
        plt.savefig(str(path_hm),
                    dpi=500,
                    format = 'png')

        return summ_res




    def pickle_results(self):
        
        import pickle
        safepath = self.savepath
        data_name = None
        if data_name == None:
            if isinstance(self.covmatprior,bool) == False:
                res_1 = self.pkldat
                res_1['covmat_sim']     = self.covmat_sim
                res_1['nfac_sim']   = self.nfac_sim 
                pickle.dump(res_1, open(safepath+'/'+"befa_res_uninf.p", "wb" ) )
        
            elif isinstance(self.covmatprior,tuple) == True:           
                res_2 = self.pkldat
                res_2['covmat_sim'] = self.covmat_sim
                res_2['nfac_sim']   = self.nfac_sim 
                pickle.dump(res_2, open( safepath+'/'+"befa_res_inf.p", "wb" ) )

        else:
            res_3 = self.pkldat
            res_3['covmat_sim'] = self.covmat_sim
            res_3['nfac_sim']   = self.nfac_sim 
            pickle.dump(res_3, open( safepath+'/'+ data_name +".p", "wb" ) )
            




