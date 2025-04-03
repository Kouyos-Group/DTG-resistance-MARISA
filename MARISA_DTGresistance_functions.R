
# Main model
# Forcing function, smooth, or linear
smooth_f = function(x, a, b, c, d){
  if(a == b){
    if(x < a){
      return(c) 
    }
    else{
      return(d) 
    }
  }
  else{
    if(x < a){
      return(c) 
    }
    else if(x > b){
      return(d) 
    }
    else{
      return(c+(d-c)*(3*((x-a)/(b-a))^2-2*((x-a)/(b-a))^3))
    }
  }
}

lin_f = function(x, a, b, c, d){
  if(a==b){
    if(x<a){
      return(c)
    }
    else{
      return(d)
    }
  }
  else{
    if(x<a){
      return(c)
    }
    else if(x>b){
      return(d)
    }
    else{
      return(c+(d-c)*(x-a)/(b-a))
    }
  }
}

#position in vector according to position in the four dimension of array

# i for Carestage & regimen
# j for CD4 level
# k for sex
# l for resistance genotype


pos = function(i, j, k, l){
  return(2 + i + (j-1)*25 + (k-1)*25*4 + (l-1)*25*4*2)
}

#position for vectors
posa = function(v1, v2, v3, v4){
  r = rep(0,length(v1)*length(v2)*length(v3)*length(v4))
  for(i in 1:length(v1)){
    for(j in 1:length(v2)){
      for(k in 1:length(v3)){
        for(l in 1:length(v4)){
          r[i +
              (j-1)*length(v1) +
              (k-1)*length(v1)*length(v2) +
              (l-1)*length(v1)*length(v2)*length(v3)] = 2 + v1[i]+(v2[j]-1)*25+(v3[k]-1)*25*4+(v4[l]-1)*25*4*2
        }
      }
    }
  }
  return(r)
}



#initial distribution of individual across CD4 classes
repart_f = function(k){
  x1=1/(1+k+k^2+k^3)
  x2=k*x1
  x3=k*x2
  x4=k*x3
  return(c(x1,x2,x3,x4))
}

#initial vector, in 2005
x_start2005 = function(nr_of_genotypes){
  #Previously estimated parameters
  dim_genotype = nr_of_genotypes
  rate_ratio = 0.5
  k1 = 1
  k2 = 2
  k3 = 2
  
  #Starting point
  undiag_start = c(1308090, 1875356)/1000  # in men, in women (72% of total)
  diag_start = c(355630, 795094)/1000      # in men, in women (26% of total)
  treat_start = c(32928, 55927)/1000       # in men, in women (2% of total)
  
  start_value = array(0, dim = c(25, 4, 2, dim_genotype))
  
  # For carestage I, all cd4, both gender, both NNRTI resistance
  start_value[1, 1:4, 1:2, 1] = array(rep(repart_f(k1), 2),             dim = c(4, 2, 1)) *  # same number in all CD4
    array(rep(undiag_start, each = 4),      dim = c(4, 2, 1)) *  # Nr (per 1000) of undiagnosed (I) PLHIV
    array(rep(c(0.99), each = 4*2),         dim = c(4, 2, 1))    # 1% NNRTI Resistant, 99% susceptible
  
  start_value[1, 1:4, 1:2, 3] = array(rep(repart_f(k1), 2),              dim = c(4, 2, 1)) *  # same number in all CD4
    array(rep(rep(undiag_start, each = 4)), dim = c(4, 2, 1)) *  # Nr (per 1000) of undiagnosed (I) PLHIV
    array(rep(c(0.01), each = 4*2),      dim = c(4, 2, 1))    # 1% NNRTI Resistant, 99% susceptible
  
  
  start_value[3, 1:4, 1:2, 1] = array(rep(repart_f(k2), 2),         dim=c(4, 2, 1)) *      # x2 in each CD4 stage
    array(rep(diag_start, each = 4),    dim=c(4, 2, 1)) *      # Nr (per 1000) of diagnosed (D) PLHIV
    array(rep(c(0.99), each = 4*2),     dim=c(4, 2, 1))        # 1% NNRTI Resistant, 99% susceptible
  
  start_value[3, 1:4, 1:2, 3] = array(rep(repart_f(k2), 2),         dim=c(4, 2, 1)) *      # x2 in each CD4 stage
    array(rep(diag_start, each = 4),    dim=c(4, 2, 1)) *      # Nr (per 1000) of diagnosed (D) PLHIV
    array(rep(c(0.01), each = 4*2),     dim=c(4, 2, 1))        # 1% NNRTI Resistant, 99% susceptible
  
  
  
  start_value[4, 1:4, 1:2, 1] = array(rep(repart_f(k3), 2),          dim=c(4, 2, 1)) *     # x2 in each CD4 stage
    array(rep(treat_start, each = 4),    dim=c(4, 2, 1)) *     # Nr (per 1000) of treated (T1) PLHIV
    array(rep(c(0.99), each = 4*2),      dim=c(4, 2, 1))       # 1% NNRTI Resistant, 99% susceptible
  start_value[4, 1:4, 1:2, 3] = array(rep(repart_f(k3), 2),          dim=c(4, 2, 1)) *     # x2 in each CD4 stage
    array(rep(treat_start, each = 4),    dim=c(4, 2, 1)) *     # Nr (per 1000) of treated (T1) PLHIV
    array(rep(c(0.01), each = 4*2),      dim=c(4, 2, 1))       # 1% NNRTI Resistant, 99% susceptible
  
  start_value = c(c(0, 0), as.vector(start_value))
  
  return(start_value)
}


#initial vector, in 2020
x_start2020 = function(x, p_DTG, nr_of_genotypes){
  dim_genotype = nr_of_genotypes
  #p_DTG = proportion of women eligible for DTG
  x <- array(x[-(1:2)], dim = c(25, 4, 2, dim_genotype))
  
  #Diagnostic
  diag_m <- array(x[2, 1:4, 1, 1:dim_genotype]) + array(x[3, 1:4, 1, 1:dim_genotype])
  diag_w <- array(x[2, 1:4, 2, 1:dim_genotype]) + array(x[3, 1:4, 2, 1:dim_genotype])
  
  #NNRTI according to cd4, sex, and resistance
  nnrti_m <- array(x[4:6, 1:4, 1, 1:dim_genotype]) + array(x[7:9, 1:4, 1, 1:dim_genotype])
  nnrti_w <- array(x[4:6, 1:4, 2, 1:dim_genotype]) + array(x[7:9, 1:4, 2, 1:dim_genotype])
  
  #diagnostic
  x[2, 1:4, 1, 1:dim_genotype] = diag_m
  x[3, 1:4, 1, 1:dim_genotype] = 0
  x[2, 1:4, 2, 1:dim_genotype] = diag_w * p_DTG
  x[3, 1:4, 2, 1:dim_genotype] = diag_w * (1-p_DTG)
  
  #men: all eligible for DTG
  x[4:6, 1:4, 1, 1:dim_genotype] = 0
  x[7:9, 1:4, 1, 1:dim_genotype] = nnrti_m
  x[c(10:12, 20, 21), 1:4, 1, 1:dim_genotype] = 0 #0, as no DTG in 2019
  
  x[4:6, 1:4, 2, 1:dim_genotype] = nnrti_w * (1-p_DTG)
  x[7:9, 1:4, 2, 1:dim_genotype] = nnrti_w * p_DTG
  
  x[c(10:12, 20, 21), 1:4, 2, 1:dim_genotype] = 0 #0, as no DTG in 2019
  
  x = c(c(0, 0), as.numeric(x))
  return(x)
}


# Model used in "Predicting emergent Dolutegravir resistance in South Africa: A modelling study"
MARISA_DTGRESISTANCE = function(t, # time
                                y,         # system state { susceptible , infected , recovered }
                                theta,     # parameters { transmission rate , recovery rate }
                                matrices,  # matrices { acquisition, reversion and effects of genotypes }
                                x_r,       # real valued fixed data
                                x_i){
  
  ## INPUT PARAMETERS ####
  dtg_1st =              theta[1];
  dtg_switch =           theta[2];
  dtg_eff =              theta[4];
  p_msm =                theta[5];
  risk_r =               theta[6];
  hiv_msm =              theta[7];
  var_switch_DTG_S =     theta[8];   
  var_switch_DTG_F =     theta[9];  
  DTG_start_year =       theta[10];  # Year of DTG introduction
  switch_decrease =      theta[11];  # decrease in switching rate compared to the one from IeDEA-SA
  DTG_delay_F_to_S =     theta[12];  
  DTG_delay_DTG_to_PI =  theta[13];
  PIswitch_change =      theta[14];  # T/F: whether the fitted DTG to PI switching rates should be overwritten (DTG_to_PI_counterfactual: vector of length 4 (per CD4 level) with switch rates)
  optimize_NRTI =        theta[16];  # whether NRTI are optimized (i.e., no resistance) at 1: not optimized, 2: in fail, 3: in suppressed, 4: in fail and suppressed.
  DTGresiswitchPI =      theta[17];
  avg_months_Frecent =   theta[18];
  rescue_rate =          theta[19];  # rate at which failing PI goes to rescue therapy (where they are suppressed). Set to 0 if not using.
  sens_Npop =            theta[20];
  sens_inf =             theta[21];
  sens_diag =            theta[22];
  sens_artini =          theta[23];
  sens_mort =            theta[24];
  failratePImodif =      theta[25];
  avg_months_Finterm = 12
  LC_rate = 109/1000/12
  LCtoCare = 1/22.8 
  
  ## INPUT MATRICES ####
  A =                         matrices[[1]];
  R =                         matrices[[2]];
  E =                         matrices[[3]];
  P_NNRTI =                   matrices[[4]];  # matrix describing genotypes assigned to NNRTI
  P_NRTI =                    matrices[[5]];  # matrix describing genotypes assigned to NRTI
  P_DTG =                     matrices[[6]];  # matrix describing genotypes assigned to DTG
  genotypes =                 matrices[[7]];  # vector describing possible genotypes
  alpha_g =                   matrices[[8]];
  E_DTG_low =                 matrices[[9]];
  E_DTG_mid =                 matrices[[10]];
  E_DTG_high =                matrices[[11]];
  M_DTG_atleastlow =          matrices[[12]];
  M_DTG_atleastmid =          matrices[[13]];
  M_DTG_atleasthigh =         matrices[[14]];
  DTG_to_PI_counterfactual =  matrices[[15]];
  Tr =                        matrices[[16]]; # Mutation during transmission: Tr[i][j] = P(genotype changes from i to j during transmission)
  # j-th row of Tr is the distribution of genotypes resulting from infection with genotype j
  M_nrtisus =                 matrices[[17]]; # Matrix used in counterfactual szenario where NRTI are optimized at switch
  
  
  # Check that all columns sum to one
  stopifnot(all(rowSums(Tr)==1))
  stopifnot(all(rowSums(M_nrtisus)==1))
  # Extract number of phenotypes
  dim_genotype = nrow(A)
  
  # E2: effect of treatment for S <-> F
  # Comment: in the first 3 months, impact is higher 
  E2 <- E
  E2[E2[1:dim_genotype, 1]>1, 1:2] <- 3.24  # switch from the 1.92 (NNRTI < 3m) to the 3.24 (NNRTI > 3m) 
  
  ## DETERMINE POSITIONS ####
  #  1  -  2  -  3  -  4  -  5  -  6  -   7   -  8
  # 138 - 140 - 148 - 155 - 230 - 263 - NNRTI - NRTI
  NNRTI_res_pos          = which(substring(genotypes, first = 7, last = 7) == "1")
  NRTI_res_pos           = which(substring(genotypes, first = 8, last = 8) == "1")
  DTG_res_pos            = which(substring(genotypes, first = 1, last = 6) != "000000")
  DTG_res_pos_low        = which(substring(genotypes, first = 1, last = 6) %in% E_DTG_low)
  DTG_res_pos_mid        = which(substring(genotypes, first = 1, last = 6) %in% E_DTG_mid)
  DTG_res_pos_high       = which(substring(genotypes, first = 1, last = 6) %in% E_DTG_high)
  DTG_res_pos_midandhigh = which(substring(genotypes, first = 1, last = 6) %in% c(E_DTG_mid, E_DTG_high))
  #pos for dummy 221011
  pos_dtg_nrti           = DTG_res_pos[which(DTG_res_pos %in% NRTI_res_pos)]
  pos_nnrti_nrti         = NNRTI_res_pos[which(NNRTI_res_pos %in% NRTI_res_pos)]
  pos_dtg_nnrti_nrti     = DTG_res_pos[which(DTG_res_pos %in% NNRTI_res_pos & DTG_res_pos %in% NRTI_res_pos)]
  
  all_pos = c(1:25);                         # position of all HIV stages
  inf_pos = c(2,3, 4,6, 7,9, 10,12, 13,15, 16:24);  # position of infectious HIV stages knowing their HIV-positive status
  
  
  ## CALCULATE VARIABLES ####
  # reshape y_matrix to array format
  ypos <- array(y[3:(25 * 4 * 2 * dim_genotype + 2)], c(25, 4, 2, dim_genotype))
  
  corr_factor_incidence = (p_msm*risk_r * hiv_msm + 2.0 * (1 - p_msm)) / (0.05 * 0.008 / 0.003 + 2.0 * 0.95)
  inf1 = sens_inf * c( 0.003 * risk_r * p_msm * hiv_msm / corr_factor_incidence,  # sens_inf is for in- or decreasing transmission in revisions
                       0.003 * (1 - p_msm) / corr_factor_incidence,
                       0.003 * (1 - p_msm) / corr_factor_incidence,
                       0.0);
  inf2 = 3.318999;
  inf3 = 0.5440886;
  
  #Diagnosis
  diag1 = 273.6363;                                                         # number of months before diagnosis (without opportunistic disease - oi) in 2005: 22 years
  diag2 = smooth_f(2005.0 + t/12.0, 2005.0, 2015.0, 1.0, 7.710578);         # fold increase between 2005 and 2015: 3 years
  diag_women = 1.25;                                                        # increase of diagnosis rate (without oi) for women
  oi_inc = c(0.05/12.0, 0.12/12.0, 0.27/12.0, 0.9/12.0);                    # oi incidence by cd4
  oi_test = smooth_f(2005 + t/12.0, 2005.0, 2015.0, 0.2, 0.8);              # proportion of oi test
  preg_inc = c(1.0 * 23.0 / (12.0 * 1000.0),
               0.96 * 23.0 / (12.0 * 1000.0), 
               0.87 * 23.0/(12.0*1000.0), 
               0.74 * 23.0/(12.0 * 1000.0));     # incidence of pregnancy by cd4
  preg_test = smooth_f(2005 + t / 12.0, 2005.0, 2010.0, 0.5, 0.98);         # proportion of pregnancy test
  
  #Treatment
  p_dtg = c(1.0, theta[3]);                                 # proportion of women opting for DTG
  p_tdf = 1.0;                                              # proportion of people opting/being prescribed tdf (over azt)
  rate_treat = 0.001879695 * sens_artini;                   # free treatment parameter rates fixed in project 1 |  sens_artini is for revisions
  t_cd4_elig = c(0.4, 0.5, 0.7, 1.0);                       # relative initiation rates by cd4 in 2020 (increase to become identical from 2022, Treat-All policy)
  t_cd4_treat_all = c(1.0/0.4, 1.0/0.5, 1.0/0.7, 1.0);
  t_cd4_elig_year1 = c(2015.5, 2013.5, 2009.5, 2001.0);
  t_cd4_elig_year2 = c(2016.5, 2015.5, 2012.5, 2004.0);
  t_1st = c();
  for(i in 1:4){
    t_1st[i] = rate_treat * smooth_f(2005.0 + t/12.0, 2001.0, 2012.0, 1.0, 200.0/12.0) *       # Increase in treatment rate between 2001 and 2012
      lin_f(2005.0 + t/12.0, t_cd4_elig_year1[i],  t_cd4_elig_year2[i], 0.0, t_cd4_elig[i]) *  # Change in eligibIlity criteria
      lin_f(2005.0 + t/12.0, 2017.0, 2022.0, 1.0, t_cd4_treat_all[i]);                         # Increase in treatment rate due to Treat-All policy
  }
  switchingPI = c(0, 0, 0, 0)  # used as replacement for t_switch in counterfactual scenario
  
  t_switch = c(4.35/531.0 /switch_decrease,
               4.35/427.0 /switch_decrease,
               4.35/294.0 /switch_decrease,
               4.35/189.0 /switch_decrease);   # Switching rate to PI, failing individuals (either DTG-ineligible on NNRTI or DTG-eligible on DTG)
  t_switch_elig = c();                         # Switch rate to PI, DTG-eligible
  t_1st_NNRTI_inel = c();
  t_1st_NNRTI_elig = c();
  t_1st_DTG = c();
  t_switch_DTG_S = c();
  t_switch_DTG_F = c();
  for(i in 1:4){
    t_1st_DTG [i]= t_1st[i] * smooth_f(2005.0 + t/12.0, DTG_start_year, DTG_start_year, 0.0, 1.0) * dtg_1st; #Treatment initiation rate, DTG
    t_1st_NNRTI_inel[i] = t_1st[i];
    t_1st_NNRTI_elig[i] = t_1st[i] - t_1st_DTG[i];
    t_switch_elig[i] = t_switch[i] * smooth_f(2005.0 + t/12.0, DTG_start_year, DTG_start_year, 1.0, 1.0 - dtg_1st); #switch always 0 after 2019, except if dtg_1st=0 (no DTG introduction) in which case it remains at 1
    t_switch_DTG_S[i] = 1.0/12.0 * dtg_switch * var_switch_DTG_S * smooth_f(2005.0 + t/12.0, DTG_start_year, DTG_start_year, 0.0, 1.0);
    if(dtg_1st == 0.0) t_switch_DTG_F[i] = 0.0;
    if(dtg_1st == 1.0 & dtg_switch == 0.0) t_switch_DTG_F[i] = t_switch[i];
    if(dtg_switch == 1.0) t_switch_DTG_F[i] = 1.0/12.0 * var_switch_DTG_F;
    t_switch_DTG_F[i] = t_switch_DTG_F[i] * smooth_f(2005.0 + t/12.0, DTG_start_year, DTG_start_year, 0.0, 1.0);
  }
  
  # Mortality
  rate_death = 0.1632000 * sens_mort;    # Mortality rate per month per 1000 people (reference, supp people with cd4 > 500) | sens_mort is for revisions
  t_prov = c()
  mort_approx = c()
  for(l in 1:37){
    t_prov = t-l+1.0
    mort_approx[l] = rate_treat * smooth_f(2005.0 + t_prov/12.0, 2001.0, 2012.0, 1.0, 200.0/12.0) *  # Increase in treatment rate between 2001 and 2012
      lin_f(2005.0 + t_prov/12.0, t_cd4_elig_year1[4],  t_cd4_elig_year2[4], 0.0, t_cd4_elig[4]) *   # Change in eligibility criteria
      lin_f(2005.0 + t_prov/12.0, 2017.0, 2022.0, 1.0, t_cd4_treat_all[4])
  }
  p1 = 1.0 - 0.27*exp(-0.05*(mean(mort_approx) - 0.005219697)/rate_treat)  #p1 =0.865181; # proportion of people with cd4>50 among people with cd4<200
  mu = matrix(c(1.57 * rate_death/1000.0, 2.0 * rate_death/1000.0, 4.57 * rate_death/1000.0, (p1 * 40.9 + (1 - p1) * 134.4) * rate_death/1000.0, # relative mortality for untreated, by cd4
                (0.9 * 1 + 0.1 * 3.92) * rate_death/1000.0, (0.9 * 1.26 + 0.1 * 3.92) * rate_death/1000.0,                                     # for people starting treatment, by cd4
                (0.9 * 1.94 + 0.1 * 4.28) * rate_death/1000.0, (0.9 * (p1*8.3+(1-p1)*41.7) + 0.1 * (p1*11.8+(1-p1)*59.7)) * rate_death/1000.0,
                1.0 * rate_death/1000.0, 1.26 * rate_death/1000.0, 1.94 * rate_death/1000.0, (p1*8.3+(1-p1)*41.7) * rate_death/1000.0,         # for suppressed people, by cd4
                3.92 * rate_death/1000.0, 3.92 * rate_death/1000.0, 4.28 * rate_death/1000.0, (p1*11.8+(1-p1)*59.7) * rate_death/1000.0),      # for people failing treatment, by cd4
              nrow = 4, byrow = T);
  
  # CD4 progression, for untreated, treated with NNRTI, treated with PI
  cd4_untreated = c(4.35 / 260.0, 4.35 / 156.0, 4.35 / 182.0); #untreated
  cd4_j = c()
  cd4  = matrix(c(4.35 / 206.0, 4.35 / 133.0, 4.35 / 263.0,    # NNRTI, start treatment, right
                  4.35 / 69.0, 4.35 / 70.0, 4.35 / 79.0,       # start treatment, left
                  4.35 / 73.0, 4.35 / 61.0, 4.35 / 41.0,       # suppressed
                  4.35 / 77.0, 4.35 / 66.0, 4.35 / 96.0,       # failing
                  4.35 / 140.0, 4.35 / 98.0, 4.35 / 145.0,     # PI,start treatment, right
                  4.35 / 68.0, 4.35 / 81.0, 4.35 / 177.0,      # start treatment, left
                  4.35 / 72.0, 4.35 / 59.0, 4.35 / 32.0,       # suppressed
                  4.35 / 61.0, 4.35 / 65.0, 4.35 / 69.0),      # failing
                nrow = 2, byrow = T);
  
  # Treatment suppression and failure rates, for NNRTI, DTG and PI, by CD4
  treat_j = c()
  treat = matrix(c(2.0 * 4.35 / 30.0,   2.0 * 4.35 / 30.0,   2.0 * 4.35 / 31.0,   2.0 * 4.35 / 34.0,   # NNRTI & DTG, T to S
                   2.0 * 4.35 / 203.0,  2.0 * 4.35 / 198.0,  2.0 * 4.35 / 164.0,  2.0 * 4.35 / 112.0,  # T to F
                   4.35 / 767.0,        4.35 / 582.0,        4.35 / 270.0,        4.35 / 96.0,         # S to F
                   4.35 / 28.0,         4.35 / 56.0,         4.35 / 62.0,         4.35 / 79.0,         # F to S   
                   
                   2.0 * 4.35 / 33.0,   2.0 * 4.35 / 33.0,   2.0 * 4.35 / 35.0,   2.0 * 4.35 / 43.0,   # PI
                   2.0 * 4.35 / 124.0,  2.0 * 4.35 / 122.0,  2.0 * 4.35 / 103.0,  2.0 * 4.35 / 66.0,   # T to F
                   4.35 / 267.0,        4.35 / 178.0,        4.35 / 174.0,        4.35 / 83.0,         # S to F
                   4.35 / 10.0,         4.35 / 20.0,         4.35 / 24.0,         4.35 / 51.0),        # F to S 
                 nrow = 2, byrow = T);

  # CD4 progression and ART efficacy by ART:
  #   1: CD4 progression/efficacy equivalent to NNRTI
  #   2: equivalent to PI
  art_eq  = c(1, 1, 1, 2);  # NNRTI(inel), NNRTI(el), DTG, PI
  art_j = c()
  
  
  # N values by gender, used to defined susceptible by gender y[1], y[2]
  N_value = sens_Npop * c(39845943.06, 40029729.5, 40544491.9, 41221762.41, 42376196.21, 43678397.16, 44530813.78, 44878159.38, 46527810.82, 47609855.86, 48395856.04, 49250847.33, 50105838.63, 50960829.93)/1000
  N_v <- smooth.spline(1:14, N_value, df = 5, all.knots = FALSE)
  N_value_f = function(x){return(predict(N_v, x/12)$y * c(0.5, 0.5))}
  N <- N_value_f(t + 1)
  
  
  # susceptible by gender
  y[1] = N[1] - sum(ypos[all_pos, 1:4, 1, 1:dim_genotype]);
  y[2] = N[2] - sum(ypos[all_pos, 1:4, 2, 1:dim_genotype]);
  
  
  ## DIFFERENTIAL EQUATIONS ####
  dy_dt_pos <- array(0, c(25, 4, 2, dim_genotype))

  dummy_rates_F1NNRTI = rep(0, 4)
  dummy_rates_F1DTG = rep(0, 12)
  dummy_switch_to_PI = rep(0, 6)
  dummy_NRswitching = rep(0, 5)
  dummy_GRTnumber = rep(0, 4)
  dummy_epidemic = rep(0, 3) # infected male, infected female, mortality
  
  # Newly infected
  dy_dt_pos[1, 1, 1, 1:dim_genotype] = dy_dt_pos[1, 1, 1, 1:dim_genotype] + y[1] / (y[1] + sum(ypos[all_pos, 1:4, 1, 1:dim_genotype])) * t(Tr) %*% (
    inf1[1] * inf2 * (colSums(ypos[1, 1:4, 1, 1:dim_genotype]) + inf3 * colSums(ypos[inf_pos ,1:4, 1, 1:dim_genotype], dims=2)) +
      inf1[2] * inf2 * (colSums(ypos[1, 1:4, 2, 1:dim_genotype]) + inf3 * colSums(ypos[inf_pos, 1:4, 2, 1:dim_genotype], dims=2))
  );
  dy_dt_pos[1, 1, 2, 1:dim_genotype] = dy_dt_pos[1, 1, 2, 1:dim_genotype] + y[2] / (y[2] + sum(ypos[all_pos, 1:4, 2, 1:dim_genotype])) * t(Tr) %*% (
    inf1[3] * inf2 * (colSums(ypos[1, 1:4, 1, 1:dim_genotype]) + inf3 * colSums(ypos[inf_pos, 1:4, 1, 1:dim_genotype], dims=2)) + 
      inf1[4] * inf2 * (colSums(ypos[1, 1:4, 2, 1:dim_genotype]) + inf3 * colSums(ypos[inf_pos, 1:4, 2, 1:dim_genotype], dims=2))
  )
  
  
  ### dummy epidemic newly infected male & female
  dummy_epidemic[1] = dummy_epidemic[1] +
    y[1] / (y[1] + sum(ypos[all_pos, 1:4, 1, 1:dim_genotype])) * sum(t(Tr) %*% (
      inf1[1] * inf2 * (colSums(ypos[1, 1:4, 1, 1:dim_genotype]) + inf3 * colSums(ypos[inf_pos ,1:4, 1, 1:dim_genotype], dims=2)) +
        inf1[2] * inf2 * (colSums(ypos[1, 1:4, 2, 1:dim_genotype]) + inf3 * colSums(ypos[inf_pos, 1:4, 2, 1:dim_genotype], dims=2)))) 
  
  dummy_epidemic[2] = dummy_epidemic[2] +
    y[2] / (y[2] + sum(ypos[all_pos, 1:4, 2, 1:dim_genotype])) * sum(t(Tr) %*% (
      inf1[3] * inf2 * (colSums(ypos[1, 1:4, 1, 1:dim_genotype]) + inf3 * colSums(ypos[inf_pos, 1:4, 1, 1:dim_genotype], dims =2)) +
        inf1[4] * inf2 * (colSums(ypos[1, 1:4, 2, 1:dim_genotype]) + inf3 * colSums(ypos[inf_pos, 1:4, 2, 1:dim_genotype], dims=2))));
  
  
  # Diagnosis
  for(j in 1:4){
    #Diagnosis
    #men
    dy_dt_pos[1, j, 1, ] = dy_dt_pos[1, j, 1, ] - sens_diag * (diag2 / diag1 + oi_inc[j] * oi_test) * ypos[1, j, 1, ];
    #women
    dy_dt_pos[1, j, 2, ] = dy_dt_pos[1, j, 2, ] - sens_diag * (diag2 / diag1 * diag_women + oi_inc[j] * oi_test + preg_inc[j] * preg_test) * ypos[1, j, 2, ];
    
    #DTG-eligible, men
    dy_dt_pos[2, j, 1, ] = dy_dt_pos[2, j, 1, ] + p_dtg[1] * sens_diag * (diag2 / diag1 + oi_inc[j] * oi_test) * ypos[1, j, 1, ];
    #DTG-ineligible, men
    dy_dt_pos[3, j, 1, ] = dy_dt_pos[3, j, 1, ] + (1 - p_dtg[1]) * sens_diag * (diag2 / diag1 + oi_inc[j] * oi_test) * ypos[1, j, 1, ];
    #DTG-eligible, women
    dy_dt_pos[2, j, 2, ] = dy_dt_pos[2, j, 2, ] + p_dtg[2] * sens_diag * (diag2 / diag1 * diag_women + oi_inc[j] * oi_test + preg_inc[j] * preg_test) * ypos[1, j, 2, ];
    #DTG-ineligible, women
    dy_dt_pos[3, j, 2, ] = dy_dt_pos[3, j, 2, ] + (1 - p_dtg[2]) * sens_diag * (diag2 / diag1 * diag_women + oi_inc[j] * oi_test + preg_inc[j] * preg_test) * ypos[1, j, 2, ];
    
    
    
    #Treatment initiation
    
    #NNRTI: 1) NNRTI+TDF for DTG-inel
    #       2) NNRTI+TDF for DTG-elig
    dy_dt_pos[4, j, , ] = dy_dt_pos[4, j, , ] + t_1st_NNRTI_inel[j] * ypos[3, j, , ];
    dy_dt_pos[7, j, , ] = dy_dt_pos[7, j, , ] + t_1st_NNRTI_elig[j] * ypos[2, j, , ];
    
    #DTG: DTG+TDF
    dy_dt_pos[10, j, , ] = dy_dt_pos[10, j, , ] + t_1st_DTG[j] * ypos[2, j, , ];
    dummy_NRswitching[5] = dummy_NRswitching[5] + sum(t_1st_DTG[j] * ypos[2, j, , ])
    
    dy_dt_pos[3, j, , ] = dy_dt_pos[3, j, , ] - t_1st_NNRTI_inel[j] * ypos[3, j, , ];
    dy_dt_pos[2, j, , ] = dy_dt_pos[2, j, , ] - (t_1st_NNRTI_elig[j] + t_1st_DTG[j]) * ypos[2, j, , ];
    
    
    #Treatment switch to PI after failure
    #Switch from NNRTI to PI, DTG-ineligible
    dy_dt_pos[13, j, , ] = dy_dt_pos[13, j, , ] + t_switch[j] * ypos[6, j, , ];
    dy_dt_pos[6, j, , ]  = dy_dt_pos[6, j, , ] - t_switch[j] * ypos[6, j, , ];
    
    #Switch from NNRTI to PI, DTG-eligible
    dy_dt_pos[13, j, , ] = dy_dt_pos[13, j, , ] + t_switch_elig[j] * ypos[9, j, , ];
    dy_dt_pos[9, j, , ]  = dy_dt_pos[9, j, , ] - t_switch_elig[j] * ypos[9, j, , ];
    dummy_rates_F1NNRTI[1] = t_switch_elig[1]
    
    #Switch from DTG to PI
    # PIswitch_change T/F whether the switch rates should be shifted (input from theta)
    # DTG_to_PI_counterfactual: vector with rates (input from matrices
    DTG_to_PI_counterfactual_start = 2024
    
    if(t >= ((DTG_to_PI_counterfactual_start-2020)*12) + 180){switchingPI = DTG_to_PI_counterfactual} else{switchingPI = t_switch}
    dy_dt_pos[13, j, , ] = dy_dt_pos[13, j, , ] + 
      DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[j], t_switch[j]) * ypos[12, j, , ] + 
      DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[j], t_switch[j]) * ypos[20, j, , ] +
      DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[j], t_switch[j]) * ypos[21, j, , ];
    dy_dt_pos[12, j, , ] = dy_dt_pos[12, j, , ] - DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[j], t_switch[j]) * ypos[12, j, , ];
    dy_dt_pos[20, j, , ] = dy_dt_pos[20, j, , ] - DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[j], t_switch[j]) * ypos[20, j, , ];
    dy_dt_pos[21, j, , ] = dy_dt_pos[21, j, , ] - DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[j], t_switch[j]) * ypos[21, j, , ];
    
    # 1: default, 
    # 2: in f>6m (20, 21) within 6months
    # 3: in f>6m (20, 21) within 12months
    # 4: in f>1.5y (21) within 6months
    # 5: in f>1.5y (21) within 12months
    if(DTGresiswitchPI > 1 & t >= ((DTG_to_PI_counterfactual_start-2020)*12) + 180){
      GRTswitch = c(NA,
                    1/6, 
                    1/12, 
                    1/6, 
                    1/12)
      switchrate = GRTswitch[DTGresiswitchPI]
      GRTstages = list(NA, 
                       c(20, 21),
                       c(20, 21),
                       21, 
                       21)
      
      for(s in GRTstages[[DTGresiswitchPI]]){
        dy_dt_pos[13, , , DTG_res_pos_midandhigh] = dy_dt_pos[13, , , DTG_res_pos_midandhigh] + switchrate * ypos[s, , , DTG_res_pos_midandhigh];
        dy_dt_pos[ s, , , DTG_res_pos_midandhigh] = dy_dt_pos[ s, , , DTG_res_pos_midandhigh] - switchrate * ypos[s, , , DTG_res_pos_midandhigh];
      }
    }
    
    
    dummy_rates_F1DTG[1] = DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[1], t_switch[1])
    dummy_rates_F1DTG[2] = DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[2], t_switch[2])
    dummy_rates_F1DTG[3] = DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[3], t_switch[3])
    dummy_rates_F1DTG[4] = DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[4], t_switch[4])
    
    # c("All", "DTGmid+", "CD4_1", "CD4_2", "CD4_3", "CD4_4")
    dummy_switch_to_PI[1] = dummy_switch_to_PI[1] + 
      sum(DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[j], t_switch[j]) * ypos[12, j, , ]) + 
      sum(t_switch[j] * ypos[6, j, , ]) +
      sum(t_switch_elig[j] * ypos[9, j, , ])
    dummy_switch_to_PI[2] = dummy_switch_to_PI[2] +
      sum(DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[j], t_switch[j]) * ypos[12, j, , DTG_res_pos_midandhigh]) + 
      sum(t_switch[j] * ypos[6, j, , DTG_res_pos_midandhigh]) +
      sum(t_switch_elig[j] * ypos[9, j, , DTG_res_pos_midandhigh])
    dummy_switch_to_PI[3] = ifelse(j==1, sum(DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[j], t_switch[j]) * ypos[12, j, , ]) + 
                                     sum(t_switch[j] * ypos[6, j, , ]) +
                                     sum(t_switch_elig[j] * ypos[9, j, , ]), 
                                   dummy_switch_to_PI[3])
    dummy_switch_to_PI[4] = ifelse(j==2, sum(DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[j], t_switch[j]) * ypos[12, j, , ]) + 
                                     sum(t_switch[j] * ypos[6, j, , ]) +
                                     sum(t_switch_elig[j] * ypos[9, j, , ]), 
                                   dummy_switch_to_PI[4])
    dummy_switch_to_PI[5] = ifelse(j==3, sum(DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[j], t_switch[j]) * ypos[12, j, , ]) + 
                                     sum(t_switch[j] * ypos[6, j, , ]) +
                                     sum(t_switch_elig[j] * ypos[9, j, , ]), 
                                   dummy_switch_to_PI[5])
    dummy_switch_to_PI[6] = ifelse(j==4, sum(DTG_delay_DTG_to_PI * ifelse(PIswitch_change == 1, switchingPI[j], t_switch[j]) * ypos[12, j, , ]) + 
                                     sum(t_switch[j] * ypos[6, j, , ]) +
                                     sum(t_switch_elig[j] * ypos[9, j, , ]), 
                                   dummy_switch_to_PI[6])
    
    
    
    #Treatment switch to DTG,  #HERE FOR NRTI OPTIMIZATION AT SWITCH FOR COUNTERFACTUAL
    if(optimize_NRTI %in% c(3, 4)){
      for(g in 1:dim_genotype){
        from = genotypes[g]
        if(substring(from, first = 8, last = 8) == "0"){
          dy_dt_pos[11, j, , g] = dy_dt_pos[11, j, , g] + t_switch_DTG_S[j] * ypos[8, j, , g]
          dy_dt_pos[ 8, j, , g] = dy_dt_pos[ 8, j, , g] - t_switch_DTG_S[j] * ypos[8, j, , g]
          dummy_NRswitching[3] = dummy_NRswitching[3] + sum(t_switch_DTG_S[j] * ypos[8, j, , g])
          dummy_NRswitching[4] = dummy_NRswitching[4] + sum(t_switch_DTG_S[j] * ypos[8, j, , g]) * ifelse(g %in% NRTI_res_pos, 1, 0)
        } else{
          dy_dt_pos[11, j, , g-1] = dy_dt_pos[11, j, , g-1] + t_switch_DTG_S[j] * ypos[8, j, , g]
          dy_dt_pos[ 8, j, , g]   = dy_dt_pos[ 8, j, , g]   - t_switch_DTG_S[j] * ypos[8, j, , g]
          dummy_NRswitching[3] = dummy_NRswitching[3] + sum(t_switch_DTG_S[j] * ypos[8, j, , g])
          dummy_NRswitching[4] = dummy_NRswitching[4] + sum(t_switch_DTG_S[j] * ypos[8, j, , g]) * ifelse(g %in% NRTI_res_pos, 1, 0)
        }
      }
    } else{
      dy_dt_pos[11, j, , ] = dy_dt_pos[11, j, , ] + t_switch_DTG_S[j] * ypos[8, j, , ]
      dy_dt_pos[ 8, j, , ] = dy_dt_pos[ 8, j, , ] - t_switch_DTG_S[j] * ypos[8, j, , ]
      dummy_NRswitching[3] = dummy_NRswitching[3] + sum(t_switch_DTG_S[j] * ypos[8, j, , ])
      dummy_NRswitching[4] = dummy_NRswitching[4] + sum(t_switch_DTG_S[j] * ypos[8, j, , NRTI_res_pos])
    }
    
    if(optimize_NRTI %in% c(2, 4)){
      for(g in 1:dim_genotype){
        from = genotypes[g]
        if(substring(from, first = 8, last = 8) == "0"){
          dy_dt_pos[10, j, , g] = dy_dt_pos[10, j, , g] + t_switch_DTG_F[j] * ypos[9, j, , g]
          dy_dt_pos[ 9, j, , g] = dy_dt_pos[ 9, j, , g] - t_switch_DTG_F[j] * ypos[9, j, , g]
          dummy_NRswitching[1] = dummy_NRswitching[1] + sum(t_switch_DTG_F[j] * ypos[9, j, , g])
          dummy_NRswitching[2] = dummy_NRswitching[2] + sum(t_switch_DTG_F[j] * ypos[9, j, , g]) * ifelse(g %in% NRTI_res_pos, 1, 0)
        } else{
          dy_dt_pos[10, j, , g-1] = dy_dt_pos[10, j, , g-1] + t_switch_DTG_F[j] * ypos[9, j, , g]
          dy_dt_pos[ 9, j, , g]   = dy_dt_pos[ 9, j, , g]   - t_switch_DTG_F[j] * ypos[9, j, , g]
          dummy_NRswitching[1] = dummy_NRswitching[1] + sum(t_switch_DTG_F[j] * ypos[9, j, , g])
          dummy_NRswitching[2] = dummy_NRswitching[2] + sum(t_switch_DTG_F[j] * ypos[9, j, , g]) * ifelse(g %in% NRTI_res_pos, 1, 0)
        }
      }
    } else{
      dy_dt_pos[10, j, , ] = dy_dt_pos[10, j, , ] + t_switch_DTG_F[j] * ypos[9, j, , ]
      dy_dt_pos[ 9, j, , ] = dy_dt_pos[ 9, j, , ] - t_switch_DTG_F[j] * ypos[9, j, , ]
      dummy_NRswitching[1] = dummy_NRswitching[1] + sum(t_switch_DTG_F[j] * ypos[9, j, , ])
      dummy_NRswitching[2] = dummy_NRswitching[2] + sum(t_switch_DTG_F[j] * ypos[9, j, , NRTI_res_pos])
    }
    
    dummy_rates_F1NNRTI[4] = t_switch_DTG_F[j]
  }
  
  #Treatment stages and CD4, mortality
  #Infected and diagnosed
  
  dy_dt_pos[1:3, 1, , ] = dy_dt_pos[1:3, 1, , ] - cd4_untreated[1] * ypos[1:3, 1, , ];
  dy_dt_pos[1:3, 2, , ] = dy_dt_pos[1:3, 2, , ] + cd4_untreated[1] * ypos[1:3, 1, , ] - cd4_untreated[2] * ypos[1:3, 2, , ];
  dy_dt_pos[1:3, 3, , ] = dy_dt_pos[1:3, 3, , ] + cd4_untreated[2] * ypos[1:3, 2, , ] - cd4_untreated[3] * ypos[1:3, 3, , ];
  dy_dt_pos[1:3, 4, , ] = dy_dt_pos[1:3, 4, , ] + cd4_untreated[3] * ypos[1:3, 3, , ];
  
  #Mortality
  for(j in 1:4){ # CD4 stages
    dy_dt_pos[1:3, j, , ] = dy_dt_pos[1:3, j, , ] - mu[1, j] * ypos[1:3, j, , ];
    dummy_epidemic[3] = dummy_epidemic[3] + sum(mu[1, j] * ypos[1:3, j, , ])
  }
  
  #On ART
  art_line_1st =   c(1, 1, 1, 0)  # used in dummy compartment, first 3 correspond to NNRTI inel, NNRTI el, DTG, and the forth to PI -> using as factor to only get fail on first
  art_line_nnrti = c(1, 1, 0, 0)  # -> using as factor to only get fail on NNRTI
  art_line_dtg =   c(0, 0, 1, 0)  # -> using as factor to only get fail on DTG
  DTG_delay_F_to_S_by_ART =  c(1, 1, DTG_delay_F_to_S, 1)  # not DTG compartments (used to apply DTG_delay_F_to_S to DTG only)
  fordummy_genotype_nnrti_res = c(1:dim_genotype) %in% NNRTI_res_pos %>% as.numeric()
  fordummy_genotype_nrti_res = c(1:dim_genotype) %in% NRTI_res_pos %>% as.numeric()
  fordummy_genotype_dtg_res = c(1:dim_genotype) %in% DTG_res_pos %>% as.numeric()
  fordummy_genotype_dtg_res_midhigh = c(1:dim_genotype) %in% DTG_res_pos_midandhigh %>% as.numeric()
  
  progress_f = c(0, 0, 1, 0)
  modify_pi_fail_rate = c(1, failratePImodif)
  for(art in 1:4){
    art_j = art_eq[art];
    modify_pi_fail_rate_j = modify_pi_fail_rate[art_j]
    for(j in 1:12){cd4_j[j] = cd4[art_j,j];}
    
    treat_j <- array(0, c(25, dim_genotype))  # not sure the 25 is correct here, might go back to 24
    #treat_j[1:4] = rep(0, 4);
    #T to F
    for(i in 5:8){treat_j[i, ] =  modify_pi_fail_rate_j * treat[art_j, i] / alpha_g[art] * E[, art];}
    #T to S
    for(i in 1:4){treat_j[i, ] =  1.0/3.0 - treat_j[i + 4, ];}
    #S to F
    for(i in 9:12){treat_j[i, ] =  treat[art_j,i] * E2[, art] / alpha_g[art];}  
    #F to S
    for(i in 13:16){treat_j[i, ] =  DTG_delay_F_to_S_by_ART[art] * (treat[art_j, i] / E2[, art] * alpha_g[art]);}
    #Frecent to Fintermediate
    for(i in 17:20){treat_j[i, ] = 1.0/avg_months_Frecent;}
    #Finterediate to Flong
    for(i in 21:24){treat_j[i, ] =  1.0/avg_months_Finterm;}
    
    if(art == 3){
      stopifnot(any(is.na(treat_j)) == F)
      dummy_rates_F1DTG[5] <- treat_j[13, 1] # CD4 1, in wt
      dummy_rates_F1DTG[6] <- treat_j[14, 1] # CD4 2, in wt
      dummy_rates_F1DTG[7] <- treat_j[15, 1] # CD4 3, in wt
      dummy_rates_F1DTG[8] <- treat_j[16, 1] # CD4 4, in wt
    }
    if(art == 2){dummy_rates_F1NNRTI[2] = treat_j[13, 4]}
    
    for(k in 1:2){
      #Start
      dy_dt_pos[3*art+1, 1, k, ] = dy_dt_pos[3*art+1, 1, k, ] + cd4_j[4] * ypos[3*art+1, 2, k, ]                                     - cd4_j[1] * ypos[3*art+1, 1, k, ]             - (treat_j[1, ] + treat_j[5, ]) * ypos[3*art+1, 1, k, ];
      dy_dt_pos[3*art+1, 2, k, ] = dy_dt_pos[3*art+1, 2, k, ] + cd4_j[1] * ypos[3*art+1, 1, k, ] + cd4_j[5] * ypos[3*art+1, 3, k, ] - (cd4_j[4] + cd4_j[2]) * ypos[3*art+1, 2, k, ] - (treat_j[2, ] + treat_j[6, ]) * ypos[3*art+1, 2, k, ];
      dy_dt_pos[3*art+1, 3, k, ] = dy_dt_pos[3*art+1, 3, k, ] + cd4_j[2] * ypos[3*art+1, 2, k, ] + cd4_j[6] * ypos[3*art+1, 4, k, ] - (cd4_j[5] + cd4_j[3]) * ypos[3*art+1, 3, k, ] - (treat_j[3, ] + treat_j[7, ]) * ypos[3*art+1, 3, k, ];
      dy_dt_pos[3*art+1, 4, k, ] = dy_dt_pos[3*art+1, 4, k, ] + cd4_j[3] * ypos[3*art+1, 3, k, ]                                     - cd4_j[6] * ypos[3*art+1, 4, k, ]             - (treat_j[4, ] + treat_j[8, ]) * ypos[3*art+1, 4, k, ];
      
      #Suppressed
      dy_dt_pos[3*art+2, 1, k, ] = dy_dt_pos[3*art+2, 1, k, ] + cd4_j[7] * ypos[3*art+2, 2, k, ]                                    + treat_j[1, ] * ypos[3*art+1, 1, k, ] + treat_j[13, ] * ypos[3*art+3, 1, k, ] - treat_j[9, ]  * ypos[3*art+2, 1, k, ] + treat_j[13, ] * ypos[2*art+14, 1, k, ] + treat_j[13, ] * ypos[2*art+15, 1, k, ];
      dy_dt_pos[3*art+2, 2, k, ] = dy_dt_pos[3*art+2, 2, k, ] + cd4_j[8] * ypos[3*art+2, 3, k, ] - cd4_j[7] * ypos[3*art+2, 2, k, ] + treat_j[2, ] * ypos[3*art+1, 2, k, ] + treat_j[14, ] * ypos[3*art+3, 2, k, ] - treat_j[10, ] * ypos[3*art+2, 2, k, ] + treat_j[14, ] * ypos[2*art+14, 2, k, ] + treat_j[14, ] * ypos[2*art+15, 2, k, ];
      dy_dt_pos[3*art+2, 3, k, ] = dy_dt_pos[3*art+2, 3, k, ] + cd4_j[9] * ypos[3*art+2, 4, k, ] - cd4_j[8] * ypos[3*art+2, 3, k, ] + treat_j[3, ] * ypos[3*art+1, 3, k, ] + treat_j[15, ] * ypos[3*art+3, 3, k, ] - treat_j[11, ] * ypos[3*art+2, 3, k, ] + treat_j[15, ] * ypos[2*art+14, 3, k, ] + treat_j[15, ] * ypos[2*art+15, 3, k, ];
      dy_dt_pos[3*art+2, 4, k, ] = dy_dt_pos[3*art+2, 4, k, ]                                    - cd4_j[9] * ypos[3*art+2, 4, k, ] + treat_j[4, ] * ypos[3*art+1, 4, k, ] + treat_j[16, ] * ypos[3*art+3, 4, k, ] - treat_j[12, ] * ypos[3*art+2, 4, k, ] + treat_j[16, ] * ypos[2*art+14, 4, k, ] + treat_j[16, ] * ypos[2*art+15, 4, k, ];
      
      #Failed
      dy_dt_pos[3*art+3, 1, k, ] = dy_dt_pos[3*art+3, 1, k, ]                                     - cd4_j[10] * ypos[3*art+3, 1, k, ] + treat_j[5, ] * ypos[3*art+1, 1, k, ] + treat_j[9, ]  * ypos[3*art+2, 1, k, ] - treat_j[13, ] * ypos[3*art+3, 1, k, ] - progress_f[art] * treat_j[17, ] * ypos[3*art+3, 1, k, ];
      dy_dt_pos[3*art+3, 2, k, ] = dy_dt_pos[3*art+3, 2, k, ] + cd4_j[10] * ypos[3*art+3, 1, k, ] - cd4_j[11] * ypos[3*art+3, 2, k, ] + treat_j[6, ] * ypos[3*art+1, 2, k, ] + treat_j[10, ] * ypos[3*art+2, 2, k, ] - treat_j[14, ] * ypos[3*art+3, 2, k, ] - progress_f[art] * treat_j[18, ] * ypos[3*art+3, 2, k, ];
      dy_dt_pos[3*art+3, 3, k, ] = dy_dt_pos[3*art+3, 3, k, ] + cd4_j[11] * ypos[3*art+3, 2, k, ] - cd4_j[12] * ypos[3*art+3, 3, k, ] + treat_j[7, ] * ypos[3*art+1, 3, k, ] + treat_j[11, ] * ypos[3*art+2, 3, k, ] - treat_j[15, ] * ypos[3*art+3, 3, k, ] - progress_f[art] * treat_j[19, ] * ypos[3*art+3, 3, k, ];
      dy_dt_pos[3*art+3, 4, k, ] = dy_dt_pos[3*art+3, 4, k, ] + cd4_j[12] * ypos[3*art+3, 3, k, ]                                     + treat_j[8, ] * ypos[3*art+1, 4, k, ] + treat_j[12, ] * ypos[3*art+2, 4, k, ] - treat_j[16, ] * ypos[3*art+3, 4, k, ] - progress_f[art] * treat_j[20, ] * ypos[3*art+3, 4, k, ];
      
      #Fail-progression
      dy_dt_pos[2*art+14, 1, k, ] = dy_dt_pos[2*art+14, 1, k, ] + progress_f[art] * treat_j[17, ] * ypos[3*art+3, 1, k, ]                                      - cd4_j[10] * ypos[2*art+14, 1, k, ] - treat_j[13, ] * ypos[2*art+14, 1, k, ] - progress_f[art] * treat_j[21, ] * ypos[2*art+14, 1, k, ]
      dy_dt_pos[2*art+14, 2, k, ] = dy_dt_pos[2*art+14, 2, k, ] + progress_f[art] * treat_j[18, ] * ypos[3*art+3, 2, k, ] + cd4_j[10] * ypos[2*art+14, 1, k, ] - cd4_j[11] * ypos[2*art+14, 2, k, ] - treat_j[14, ] * ypos[2*art+14, 2, k, ] - progress_f[art] * treat_j[22, ] * ypos[2*art+14, 2, k, ]
      dy_dt_pos[2*art+14, 3, k, ] = dy_dt_pos[2*art+14, 3, k, ] + progress_f[art] * treat_j[19, ] * ypos[3*art+3, 3, k, ] + cd4_j[11] * ypos[2*art+14, 2, k, ] - cd4_j[12] * ypos[2*art+14, 3, k, ] - treat_j[15, ] * ypos[2*art+14, 3, k, ] - progress_f[art] * treat_j[23, ] * ypos[2*art+14, 3, k, ]
      dy_dt_pos[2*art+14, 4, k, ] = dy_dt_pos[2*art+14, 4, k, ] + progress_f[art] * treat_j[20, ] * ypos[3*art+3, 4, k, ] + cd4_j[12] * ypos[2*art+14, 3, k, ]                                      - treat_j[16, ] * ypos[2*art+14, 4, k, ] - progress_f[art] * treat_j[24, ] * ypos[2*art+14, 4, k, ]
      
      dy_dt_pos[2*art+15, 1, k, ] = dy_dt_pos[2*art+15, 1, k, ] + progress_f[art] * treat_j[21, ] * ypos[2*art+14, 1, k, ]                                      - cd4_j[10] * ypos[2*art+15, 1, k, ] - treat_j[13, ] * ypos[2*art+15, 1, k, ]
      dy_dt_pos[2*art+15, 2, k, ] = dy_dt_pos[2*art+15, 2, k, ] + progress_f[art] * treat_j[22, ] * ypos[2*art+14, 2, k, ] + cd4_j[10] * ypos[2*art+15, 1, k, ] - cd4_j[11] * ypos[2*art+15, 2, k, ] - treat_j[14, ] * ypos[2*art+15, 2, k, ]
      dy_dt_pos[2*art+15, 3, k, ] = dy_dt_pos[2*art+15, 3, k, ] + progress_f[art] * treat_j[23, ] * ypos[2*art+14, 3, k, ] + cd4_j[11] * ypos[2*art+15, 2, k, ] - cd4_j[12] * ypos[2*art+15, 3, k, ] - treat_j[15, ] * ypos[2*art+15, 3, k, ]
      dy_dt_pos[2*art+15, 4, k, ] = dy_dt_pos[2*art+15, 4, k, ] + progress_f[art] * treat_j[24, ] * ypos[2*art+14, 4, k, ] + cd4_j[12] * ypos[2*art+15, 3, k, ]                                      - treat_j[16, ] * ypos[2*art+15, 4, k, ]
      

      
      #Mortality
      for(j in 1:4){
        dy_dt_pos[3*art+1, j, k, ] = dy_dt_pos[3*art+1, j, k, ] - mu[2, j] * ypos[3*art+1, j, k, ];
        dy_dt_pos[3*art+2, j, k, ] = dy_dt_pos[3*art+2, j, k, ] - mu[3, j] * ypos[3*art+2, j, k, ];
        dy_dt_pos[3*art+3, j, k, ] = dy_dt_pos[3*art+3, j, k, ] - mu[4, j] * ypos[3*art+3, j, k, ];
        dy_dt_pos[2*art+14, j, k, ] = dy_dt_pos[2*art+14, j, k, ] - mu[4, j] * ypos[2*art+14, j, k, ];
        dy_dt_pos[2*art+15, j, k, ] = dy_dt_pos[2*art+15, j, k, ] - mu[4, j] * ypos[2*art+15, j, k, ];
        
        dummy_epidemic[3] = dummy_epidemic[3] + 
          sum(mu[2, j] * ypos[3*art+1, j, k, ] + 
                mu[3, j] * ypos[3*art+2, j, k, ] +
                mu[4, j] * (ypos[3*art+3, j, k, ] + ypos[2*art+14, j, k, ] + ypos[2*art+15, j, k, ]))
      }
    }
  }
  dummy_rates_F1NNRTI[3] = mu[4, 4]
  dummy_rates_F1DTG[ 9] = mu[4, 1]
  dummy_rates_F1DTG[10] = mu[4, 2]
  dummy_rates_F1DTG[11] = mu[4, 3]
  dummy_rates_F1DTG[12] = mu[4, 4]
  
  
  
  # LC DYNAMICS
  # LOSS-FROM-CARE (F on DTG go to Diagnosed)
  for(f in c(12, 20, 21)){
    dy_dt_pos[ f, , , ] = dy_dt_pos[f, , , ] - LC_rate * ypos[f, , , ]
    dy_dt_pos[24, , , ] = dy_dt_pos[24, , , ] + LC_rate * ypos[f, , , ]}
  # RESUMING CARE
  dy_dt_pos[24, , , ] = dy_dt_pos[24, , , ] - LCtoCare * ypos[24, , , ]
  dy_dt_pos[10, , , ] = dy_dt_pos[10, , , ] + LCtoCare * ypos[24, , , ]
  #CD4 PROGRESSION
  dy_dt_pos[24, 1, , ] = dy_dt_pos[24, 1, , ]                               - cd4[1,10] * ypos[24, 1, , ];
  dy_dt_pos[24, 2, , ] = dy_dt_pos[24, 2, , ] + cd4[1,10] * ypos[24, 1, , ] - cd4[1,11] * ypos[24, 2, , ];
  dy_dt_pos[24, 3, , ] = dy_dt_pos[24, 3, , ] + cd4[1,11] * ypos[24, 2, , ] - cd4[1,12] * ypos[24, 3, , ];
  dy_dt_pos[24, 4, , ] = dy_dt_pos[24, 4, , ] + cd4[1,12] * ypos[24, 3, , ];
  #Mortality
  for(j in 1:4){
    dy_dt_pos[24, j, , ] = dy_dt_pos[24, j, , ] - mu[4, j] * ypos[24, j, , ];
    dummy_epidemic[3] = dummy_epidemic[3] + sum(mu[4, j] * ypos[24, j, , ])
  } # CD4 stages
  
  
  # REVIEWER COMMENT: Rescue Therapy for those on failing PI
  dy_dt_pos[25, , , ] = dy_dt_pos[25, , , ] + rescue_rate * (ypos[15, , , ] + ypos[22, , , ] + ypos[23, , , ])
  dy_dt_pos[15, , , ] = dy_dt_pos[15, , , ] - rescue_rate * (ypos[15, , , ])
  dy_dt_pos[22, , , ] = dy_dt_pos[22, , , ] - rescue_rate * (ypos[22, , , ])
  dy_dt_pos[23, , , ] = dy_dt_pos[23, , , ] - rescue_rate * (ypos[23, , , ])
  # Add CD4 movement on Rescue Therapy
  dy_dt_pos[25, 1, , ] = dy_dt_pos[25, 1, , ]                              + cd4[1,7] * ypos[25, 2, , ]; # CD4 movement as in suppressed
  dy_dt_pos[25, 2, , ] = dy_dt_pos[25, 2, , ] - cd4[1,7] * ypos[25, 2, , ] + cd4[1,8] * ypos[25, 3, , ];
  dy_dt_pos[25, 3, , ] = dy_dt_pos[25, 3, , ] - cd4[1,8] * ypos[25, 3, , ] + cd4[1,9] * ypos[25, 4, , ];
  dy_dt_pos[25, 4, , ] = dy_dt_pos[25, 4, , ] - cd4[1,9] * ypos[25, 4, , ];
  # Add Mortality on Rescue Therapy
  for(j in 1:4){
    dy_dt_pos[25, j, , ] = dy_dt_pos[25, j, , ] - mu[3, j] * ypos[25, j, , ]; # Mortality as in suppressed
    dummy_epidemic[3] = dummy_epidemic[3] + sum(mu[3, j] * ypos[25, j, , ])
  }
  
  
  
  A_DTG = A * P_DTG
  A_NNRTI = A * P_NNRTI
  A_NRTI = A * P_NRTI
  
  R_DTG = R * P_DTG
  R_NNRTI = R * P_NNRTI
  R_NRTI = R * P_NRTI
  
  for(j in 1:4){  # CD4 levels
    for(k in 1:2){  # gender
  
      # Acquisition on Failing NNRTI -> NNRTI & NRTI resistance
      dy_dt_pos[6, j, k, ] = dy_dt_pos[6, j, k, ] + t(A_NRTI + A_NNRTI) %*% ypos[6, j, k, ] - rowSums(A_NRTI + A_NNRTI) * ypos[6, j, k, ];
      dy_dt_pos[9, j, k, ] = dy_dt_pos[9, j, k, ] + t(A_NRTI + A_NNRTI) %*% ypos[9, j, k, ] - rowSums(A_NRTI + A_NNRTI) * ypos[9, j, k, ];
      
      # Acquisition on Failing DTG -> DTG & NRTI resistance
      dy_dt_pos[12, j, k, ] = dy_dt_pos[12, j, k, ] + t(A_NRTI + A_DTG) %*% ypos[12, j, k, ] - rowSums(A_NRTI + A_DTG) * ypos[12, j, k, ];
      dy_dt_pos[20, j, k, ] = dy_dt_pos[20, j, k, ] + t(A_NRTI + A_DTG) %*% ypos[20, j, k, ] - rowSums(A_NRTI + A_DTG) * ypos[20, j, k, ];
      dy_dt_pos[21, j, k, ] = dy_dt_pos[21, j, k, ] + t(A_NRTI + A_DTG) %*% ypos[21, j, k, ] - rowSums(A_NRTI + A_DTG) * ypos[21, j, k, ];
      
      
      # Reversion (NRTI & NNRTI mutations)
      # Reversion in Undiagnosed
      dy_dt_pos[1, j, k, ]  = dy_dt_pos[ 1, j, k, ] + t(R_DTG + R_NNRTI + R_NRTI) %*% ypos[ 1, j, k, ] - rowSums(R_DTG + R_NNRTI + R_NRTI) * ypos[ 1, j, k, ];
      
      # Reversion in Diagnosed
      dy_dt_pos[2, j, k, ] = dy_dt_pos[2, j, k, ] + t(R_DTG + R_NNRTI + R_NRTI) %*% ypos[2, j, k, ] - rowSums(R_DTG + R_NNRTI + R_NRTI) * ypos[2, j, k, ];
      dy_dt_pos[3, j, k, ] = dy_dt_pos[3, j, k, ] + t(R_DTG + R_NNRTI + R_NRTI) %*% ypos[3, j, k, ] - rowSums(R_DTG + R_NNRTI + R_NRTI) * ypos[3, j, k, ];
      
      # Reversion in Failing NNRTI
      dy_dt_pos[6, j, k, ] = dy_dt_pos[6, j, k, ] + t(R_DTG) %*% ypos[6, j, k, ] - rowSums(R_DTG) * ypos[6, j, k, ];
      dy_dt_pos[9, j, k, ] = dy_dt_pos[9, j, k, ] + t(R_DTG) %*% ypos[9, j, k, ] - rowSums(R_DTG) * ypos[9, j, k, ];
      
      # Reversion in Failing DTG
      dy_dt_pos[12, j, k, ] = dy_dt_pos[12, j, k, ] + t(R_NNRTI) %*% ypos[12, j, k, ] - rowSums(R_NNRTI) * ypos[12, j, k, ];
      dy_dt_pos[20, j, k, ] = dy_dt_pos[20, j, k, ] + t(R_NNRTI) %*% ypos[20, j, k, ] - rowSums(R_NNRTI) * ypos[20, j, k, ];
      dy_dt_pos[21, j, k, ] = dy_dt_pos[21, j, k, ] + t(R_NNRTI) %*% ypos[21, j, k, ] - rowSums(R_NNRTI) * ypos[21, j, k, ];
      
      # Reversion in Failing PI
      dy_dt_pos[15, j, k, ] = dy_dt_pos[15, j, k, ] + t(R_NNRTI + R_DTG) %*% ypos[15, j, k, ] - rowSums(R_NNRTI + R_DTG) * ypos[15, j, k, ];
      
      # Reversion in Loss-from-care
      dy_dt_pos[24, j, k, ] = dy_dt_pos[24, j, k, ] + t(R_DTG + R_NNRTI + R_NRTI) %*% ypos[24, j, k, ] - rowSums(R_DTG + R_NNRTI + R_NRTI) * ypos[24, j, k, ];
    }
  }
  
  return(list(c(0, 0, array(dy_dt_pos), c(dummy_rates_F1NNRTI, dummy_rates_F1DTG), dummy_switch_to_PI, dummy_NRswitching, dummy_epidemic)))
}
