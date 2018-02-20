
spill_power <- function(delta = NULL, # proportion SAR increase
                        Yp = NULL, # vector of prospective years
                        cj = NULL, # sample size mulitplier
                        retro_df = NULL, # data frame of retrospective data
                        spread_factors = NULL, # list of cohort-spread factors
                        iter = 10){
 
    stopifnot(!is.null(delta),
           !is.null(Yp),
           !is.null(cj),
           !is.null(retro_df),
           !is.null(spread_factors)) 
  
    bootn <- iter
    tmp_delta <- delta
    tmp_cj <- cj
    xi_list <- spread_factors
  
  #------------------------------------------------------------------------------
  # calculate summary statistics for log transformed juveniles 
  #------------------------------------------------------------------------------
  lnJ <- retro_df %>%
    group_by(period) %>%
    summarise(lnJ = mean(log(J), na.rm = TRUE)) %>%
    dplyr::select(lnJ) %>%
    pull()
  
  Sigma_lnJ <- retro_df %>%
    mutate(lnJ = log(J)) %>%
    dplyr::select(my, period, lnJ) %>%
    spread(period, lnJ) %>%
    dplyr::select(-my) %>%
    na.omit() %>%
    as.matrix() %>%
    var(na.rm = TRUE)
  
  # vcov differs slightly from Steve's - the numbers below are for chinook H+W
  #vcov_J <- c(.69, .32, -.19, -.07, .32, .46, -.11, .12, -.19, -.11, .24, -.03, -.07, .12, -.03, .32)
  #Sigma_lnJ <- matrix(vcov_J, nrow = 4, ncol = 4) 
  
  
  #------------------------------------------------------------------------------  
  # fit retro model to get mean paramater estimates for simulation
  #------------------------------------------------------------------------------
  
  ret_mod <- glmer(cbind(A1, A0) ~ period + (1|my) + (1|obs), data = retro_df, family = binomial())
  
  # extract fixed effects terms and vcov
  beta <- fixef(ret_mod)
  Sigma_beta <- vcov(ret_mod)
  
  # extract random effects
  s2_y <-VarCorr(ret_mod)$my[1,1]
  s2_obs <-VarCorr(ret_mod)$obs[1,1]
  
  # alternate block design; first year 0011, 1100, 0011, 1100, ...
  block <- rep(c(1, 2), max(Yp))
  tmp_spill <- rep(c('0011', '1100'), max(Yp))
  
  # collect power analysis data
  tmp_p_block <- as.data.frame(matrix(NA, nrow = bootn, ncol = length(Yp)))
  tmp_p_baci <- as.data.frame(matrix(NA, nrow = bootn, ncol = length(Yp)))
  
  tmp_sar_block <- as.data.frame(matrix(NA, nrow = bootn, ncol = length(Yp)))
  tmp_sar_baci <- as.data.frame(matrix(NA, nrow = bootn, ncol = length(Yp)))

  for(t in 1:length(Yp)){ # number of years - study sample size
    
    tmp_Yp <- Yp[t]
    
    #print(paste0('Effect size ', tmp_delta, ', juvenile sample size ', tmp_cj, ' and study year ',tmp_Yp, ' at ', Sys.time(),'.'))
    
    for(b in 1:bootn){ # bootstrap iterations  

  # Step 5
  # randomly generate Yp sets of cohort specific SAR values via a - i
  
  # Step 5a - randomly generate a set of period effects
  beta_tilde <- mvrnorm(n = 1, mu = beta, Sigma = Sigma_beta)
  
  # Step 5b and 5c
  tmp_yr <- (max(retro_df$my)+1):(max(retro_df$my)+tmp_Yp)
  
  df <- tibble(year = tmp_yr,
               my = tmp_yr,
               spill = tmp_spill[1:tmp_Yp],
               y_j = rnorm(tmp_Yp, 0, sqrt(s2_y)), # generate Yp year effects
               p1 = y_j + beta_tilde[1], # create logit values for each cohort
               p2 = y_j + beta_tilde[1] + beta_tilde[2],
               p3 = y_j + beta_tilde[1] + beta_tilde[3],
               p4 = y_j + beta_tilde[1] + beta_tilde[3])
  
  # Step 5d - transform intermediate values to SAR scale
  df <- df %>%
    mutate(s1 = exp(p1) / (1 + exp(p1)),
           s2 = exp(p2) / (1 + exp(p2)),
           s3 = exp(p3) / (1 + exp(p3)),
           s4 = exp(p4) / (1 + exp(p4)))
  
  # Step 5e - for each prospective year randomly sample with replacement from
  # 1998-2014 one set of "cohort-spread factors" (Appendix C)
  
  xi_tilde <- matrix(NA, nrow = tmp_Yp, ncol = 4)
  
  p_years <- 1:tmp_Yp
  
  for(k in 1:tmp_Yp) {
    spill_pat <- block[p_years[k]] # 1 = 0011, 2 = 1100
    xi_tilde[k,] <- as.matrix(xi_list[[spill_pat]][sample(1:17, 1, replace = TRUE),][,-1])
  }
  
  colnames(xi_tilde) <- c('xi1', 'xi2', 'xi3', 'xi4')
  

  block_df <- df %>%
    bind_cols(xi_tilde %>%
                as.tibble()) %>%
    mutate(tmp_delta = tmp_delta)
  
  # Step 5f and 5g - multiply intermediate SAR value by cohort-specific effects of Gas-cap spill exposure
  block_df <- block_df %>%
    mutate(delta1_prime = (1 + (tmp_delta - 1) * xi1),
           delta2_prime = (1 + (tmp_delta - 1) * xi2),
           delta3_prime = (1 + (tmp_delta - 1) * xi3),
           delta4_prime = (1 + (tmp_delta - 1) * xi4),
           s1_tilde = s1 * delta1_prime,
           s2_tilde = s2 * delta2_prime,
           s3_tilde = s3 * delta3_prime,
           s4_tilde = s4 * delta4_prime,
           period1 = log(s1_tilde/(1 - s1_tilde)),
           period2 = log(s2_tilde/(1 - s2_tilde)),
           period3 = log(s3_tilde/(1 - s3_tilde)),
           period4 = log(s4_tilde/(1 - s4_tilde)))
  
  baci_df <- df %>% 
    mutate(tmp_delta = tmp_delta,
           s1_tilde = s1 * tmp_delta,
           s2_tilde = s2 * tmp_delta,
           s3_tilde = s3 * tmp_delta,
           s4_tilde = s4 * tmp_delta,
           period1 = log(s1_tilde/(1 - s1_tilde)),
           period2 = log(s2_tilde/(1 - s2_tilde)),
           period3 = log(s3_tilde/(1 - s3_tilde)),
           period4 = log(s4_tilde/(1 - s4_tilde)))
  
  # Convert data frame to long format to sample for obs error and create final simulated SAR estimates.
  
  block_pros_df <- block_df %>%
    dplyr::select(year, my, spill, period1:period4) %>%
    gather(key = cohort, value = block_lij_tilde, period1:period4) %>%
    arrange(my, cohort)
  
  baci_pros_df <- baci_df %>%
    dplyr::select(year, my, spill, period1:period4) %>%
    gather(key = cohort, value = baci_lij_tilde, period1:period4) %>%
    arrange(my, cohort)
  
  pros_df <- full_join(block_pros_df, baci_pros_df, by = c('year', 'my', 'spill', 'cohort'))
  
  # Step 5h - generate observation error values for each cohort and year
  # Step 5i - transform to final SAR value to each cohort
  
  pros_df <- pros_df %>%
    mutate(eij = rnorm(1:n(), 0,sqrt(s2_obs)),
           block_sij_tilde = exp(block_lij_tilde + eij) / (1 + exp(block_lij_tilde + eij)),
           baci_sij_tilde = exp(baci_lij_tilde + eij) / (1 + exp(baci_lij_tilde + eij)))
  
  # Step 6 - generate cohort sample sizes with multiplier
  
  Jmat <- as.tibble(matrix(tmp_cj * exp(mvrnorm(n = tmp_Yp, mu = lnJ, Sigma = Sigma_lnJ)),
                           ncol = 4)) %>%
    setNames(gsub('V', 'period', names(.))) %>%
    mutate(year = 2015:(2014+tmp_Yp)) %>%
    gather(key = cohort, value = Jij_tilde, period1:period4)
  
  # Step 7 - generate returning adults and non-returning adults
  
  pros_df <- full_join(pros_df, Jmat, by = c('year', 'cohort')) %>%
    mutate(Jij_tilde = round(Jij_tilde),
           period = str_extract(cohort, '[[:digit:]]+'),
           G = str_sub(spill,period, period)) %>%                 
    rowwise() %>%
    mutate(block_A1_tilde = rbinom(1, Jij_tilde, block_sij_tilde),
           block_A0_tilde = Jij_tilde - block_A1_tilde,
           baci_A1_tilde = rbinom(1, Jij_tilde, baci_sij_tilde),
           baci_A0_tilde = Jij_tilde - baci_A1_tilde) %>%
    dplyr::select(year, my, period, gas_cap = G, J = Jij_tilde, block_A1_tilde:baci_A0_tilde) %>%
    gather(key, value, block_A1_tilde:baci_A0_tilde) %>%
    separate(key, c('design', 'adults', 'tilde'), sep = '_') %>%
    spread(adults, value) %>%
    dplyr::select(-tilde) %>%
    mutate(gas_cap = ifelse(design == 'baci', 1, gas_cap)) %>%
    arrange(design, year, period)
  
  # save some summary stats for 
  b_sum <- pros_df %>%
    group_by(my, design) %>% # remove period
    summarize(J = sum(J),
              A1 = sum(A1),
              SAR = A1/J) %>%
    ungroup() %>%
    group_by(design) %>%
    summarize(SAR = mean(SAR))

  tmp_sar_baci[b,t] <- b_sum %>% filter(design == 'baci') %>% pull(SAR)
  tmp_sar_block[b,t] <- b_sum %>% filter(design == 'block') %>% pull(SAR)
  
  # Model results
  mod_df <- pros_df %>%
    dplyr::select(-my) %>%
    bind_rows(retro_df %>%
                mutate(year = my,
                       design = rep('before',n()),
                       period = as.character(period),
                       gas_cap = as.character(gas_cap)) %>%
                dplyr::select(year, period, gas_cap, J, design, A0, A1)) %>%
    mutate(my = as.factor(year),
           period = as.factor(period),
           gas_cap = as.factor(gas_cap),
           design = as.factor(design),
           obs = 1:n())

  # Step 8 - FIT MODELS - combine simulated prospective data with fixed retrospective data
  # what does this assume - all retrospective cohorts come from BiOp "0000"
  # which allows the model to estimate BiOP effects more precisely given the
  # additional years, but does the retropsective years come from the same distribution
  # of the prosepective years?  I would argue they do not, at some level all fish
  # now recieve a benefit from spill, when none in retrospective years do. In the modeling sense,
  # all years come from a distribution with an equal Var and mean of 0.  The different years come from different distributions.
  
  # calculate p-value for block spill
  
  block_mod <- mod_df %>% filter(design %in% c('before', 'block'))
  
  full <- glmer(cbind(A1, A0) ~ period + gas_cap + (1|my) + (1|obs), data = block_mod, family = binomial())
  reduced <- glmer(cbind(A1, A0) ~ period + (1|my) + (1|obs), data = block_mod, family = binomial())
  
  lik <- anova(full, reduced, test = 'Chi')
  tau <- fixef(full)[5]
  
  if(tau > 0){
    p_one <- lik[2,8]/2
  } else {
    p_one <- 1 - (lik[2,8]/2)
  }
  
  tmp_p_block[b,t] <- p_one
  
  # now for before/after design
  
  baci_mod <- mod_df %>% filter(design %in% c('before', 'baci'))        
  
  full_baci <- glmer(cbind(A1,A0) ~ period + gas_cap + (1|my) + (1|obs), data = baci_mod, family = binomial())
  reduced_baci <- glmer(cbind(A1,A0) ~ period + (1|my) + (1|obs), data = baci_mod, family = binomial())      
  
  lik_baci <- anova(full_baci, reduced_baci, test = 'Chi')
  tau_baci <- fixef(full_baci)[5]
  
  if(tau_baci > 0){
    p1 <- lik_baci[2,8]/2
  } else {
    p1 <- 1 - (lik_baci[2,8]/2)
  }
  
  tmp_p_baci[b,t] <- p1
  
  
    } # ends b loop for boot
  
  } # ends t loop
  
  names(tmp_p_block) <- Yp
  names(tmp_p_baci) <- Yp
  names(tmp_sar_block) <- Yp
  names(tmp_sar_baci) <- Yp
  
  return(list(block_pvalue = tmp_p_block, baci_pvalue = tmp_p_baci,
              block_sar = tmp_sar_block, baci_sar = tmp_sar_baci))
  
} # ends function
