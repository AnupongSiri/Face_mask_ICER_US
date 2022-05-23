################################################################################ 
# This Shiny app creat the system dynamics model and calculate ICER using      #
# mean/median of parameters.                                                   #
# Depends on:                                                                  #
#    02_ui.r                                                                   #
# To run in R using command (online by Gist.github):                           #
#    library(shiny); runGist("cd7daf9661d390dd24640840b8737723")               #    
# Author:                                                                      #
#     - Anupong Sirirungreung, <anusiri@g.ucla.edu>                            # 
################################################################################ 

rm(list = ls())

#### Load libraries ####
library(shiny)
library(tidyverse)
library(deSolve)
library(ggplot2)
require(gridExtra)
library(scales)

#########################################################
#### Cohort stratified by age group                  ####
#### A.0-17, B.18-44, C.45-64, D.65-84, E.>=85 years ####
#########################################################

NUM_COHORTS <- 5
NUM_STATES <- 7

#### Set up functions ####
# Set up stocks #
Get_stocks <- function(input){
  stocks <- c(S_A=73039150*(1-input$vc_A_value)-input$Is_A_value, 
              S_B=117818671*(1-input$vc_B_value)-input$Is_B_value, 
              S_C=83323439*(1-input$vc_C_value)-input$Is_C_value, 
              S_D=47453305*(1-input$vc_D_value)-input$Is_D_value, 
              S_E=6604958*(1-input$vc_E_value)-input$Is_E_value, 
              Sv_A=73039150*(input$vc_A_value), 
              Sv_B=117818671*(input$vc_B_value),         
              Sv_C=83323439*(input$vc_C_value), 
              Sv_D=47453305*(input$vc_D_value), 
              Sv_E=6604958*(input$vc_E_value),
              E_A=0,         E_B=0,          E_C=0,           E_D=0,         E_E=0,
              Ev_A=0,        Ev_B=0,         Ev_C=0,          Ev_D=0,        Ev_E=0,
              Is_A=input$Is_A_value,        
              Is_B=input$Is_B_value,         
              Is_C=input$Is_C_value,        
              Is_D=input$Is_D_value,        
              Is_E=input$Is_E_value,
              Ia_A=0,        Ia_B=0,         Ia_C=0,          Ia_D=0,        Ia_E=0,
              R_A=0,         R_B=0,          R_C=0,           R_D=0,         R_E=0)
  return(stocks)
}

# Set up time #
Get_time <- function(input){
  START = 0
  FINISH = 600
  STEP = 1
  simtime <- seq(START, FINISH, by=STEP)
  return(simtime)
}

Get_parms <- function(input){
  # Setup parameters #
  CohortPopulatons <- c(73039150, 117818671, 83323439, 47453305, 6604958)
  vf <- rep(input$vf_value,5) # Need to be identified
    
  # Country specific contact rate
  #CE <- as.matrix(read.csv("D:/Documents/UCLA_PhD_Epi/2023_2Spring/Simmulation_EPIDEM206/EPIDEM206_project/Mask_ICER/Mask_ICER_data/Thailand_CE_5x5.csv", header = F))
  
  # Vaccine efficacy
  VE <- rep(input$ve_value,5)  
  
  # Latent period
  D_lat  <- c(A=input$D_lat_value, B=input$D_lat_value, C=input$D_lat_value, D=input$D_lat_value, E=input$D_lat_value)
  sigma  <- 1/D_lat
  
  # Probability of symptoms
  # May need to define for vaccinated group
  # q_symp is asymptomatic infection 
  q_symp <- c(A=input$q_symp_value, B=input$q_symp_value, C=input$q_symp_value, D=input$q_symp_value, E=input$q_symp_value)
  p_symp <- 1-q_symp
  
  # Probability of developing immunity after infection
  p_imm  <- c(A=input$p_imm_value, B=input$p_imm_value, C=input$p_imm_value, D=input$p_imm_value, E=input$p_imm_value)
  
  # Infectious period duration
  D_inf  <- c(A=input$D_inf_value, B=input$D_inf_value, C=input$D_inf_value, D=input$D_inf_value, E=input$D_inf_value)
  v      <- (1/D_inf) 
  
  # Calculate beta0 to get from contact rate matrix 
  #CohortPopulatons <- c(12860865, 25105265, 18805172, 7359568, 826013)
  #beta0   <- CE/CohortPopulatons
  
  # Uniform Rt but not corespondent to age-stratified contact rate
  # May need to allow seasonal change 
  Rt     <- matrix(rep(input$R0_value, 25), nrow=5)
  beta1  <- v*(Rt/CohortPopulatons)
  
  # Apply R0 to beta0
  #R0 <- matrix(rep(input$R0_value, 25), nrow=5)
  #beta2 <- (beta0 * R0)/2
  
  # Face mask effectiveness
  ME <- c(input$ME_A_value, input$ME_B_value, input$ME_C_value, input$ME_D_value, input$ME_E_value)
 
  return(list(vf=vf, beta1=beta1, VE=VE, sigma=sigma,
           p_symp=p_symp, q_symp=q_symp, p_imm=p_imm,
           v=v, ME=ME))
}

model1 <- function(time, stocks, parms){
  with(as.list(c(stocks, parms)),{ 
    #convert the stocks vector to a matrix
    states<-matrix(stocks,nrow=NUM_COHORTS,ncol=NUM_STATES, byrow = F)
    
    S  <- states[,1] # S  - Susceptible 
    Sv <- states[,2] # Sv - Susceptible and vaccinated
    E  <- states[,3] # E  - Exposed
    Ev <- states[,4] # Ev - Exposed and vaccinated
    Is <- states[,5] # Is - Infected and symptomatic
    Ia <- states[,6] # Ia - Infected and asymptomatic
    R  <- states[,7] # R  - Recovery/Immune
    
    SSv     <- S * vf 
    SE      <- ((beta1 %*% Ia) * S * (1-ME)) + ((beta1 %*% Is) * S * (1-ME))
    SvEv    <- ((beta1 %*% Ia) * Sv * (1-VE) * (1-ME)) + ((beta1 %*% Is) * Sv * (1-VE) * (1-ME))
    EIs     <- sigma * E * p_symp
    EvIs    <- sigma * Ev * p_symp
    EIa     <- sigma * E * q_symp
    EvIa    <- sigma * Ev * q_symp
    IR      <- v * p_imm * (Ia + Is) 
    IS      <- v * (1-p_imm) * (Ia + Is) 
    
    dS_dt   <- IS - SSv - SE
    dSv_dt  <- SSv - SvEv
    dE_dt   <- SE - EIs - EIa
    dEv_dt  <- SvEv - EvIs - EvIa
    dIs_dt  <- EIs + EvIs - IR*p_symp - IS*p_symp
    dIa_dt  <- EIa + EvIa - IR*q_symp - IS*q_symp
    dR_dt   <- IR
    
    return (list(c(dS_dt, dSv_dt, dE_dt, dEv_dt, dIs_dt, dIa_dt, dR_dt)))  
  })
}

# model0: Set mask effectiveness = 0
model0 <- function(time, stocks, parms){
  with(as.list(c(stocks, parms)),{ 
    #convert the stocks vector to a matrix
    states<-matrix(stocks,nrow=NUM_COHORTS,ncol=NUM_STATES, byrow = F)
    
    S  <- states[,1] # S  - Susceptible 
    Sv <- states[,2] # Sv - Susceptible and vaccinated
    E  <- states[,3] # E  - Exposed
    Ev <- states[,4] # Ev - Exposed and vaccinated
    Is <- states[,5] # Is - Infected and symptomatic
    Ia <- states[,6] # Ia - Infected and asymptomatic
    R  <- states[,7] # R  - Recovery/Immune
    
    SSv     <- S * vf 
    SE      <- ((beta1 %*% Ia) * S ) + ((beta1 %*% Is) * S)
    SvEv    <- ((beta1 %*% Ia) * Sv * (1-VE)) + ((beta1 %*% Is) * Sv * (1-VE))
    EIs     <- sigma * E * p_symp
    EvIs    <- sigma * Ev * p_symp
    EIa     <- sigma * E * q_symp
    EvIa    <- sigma * Ev * q_symp
    IR      <- v * p_imm * (Ia + Is) 
    IS      <- v * (1-p_imm) * (Ia + Is) 
    
    dS_dt   <- IS - SSv - SE
    dSv_dt  <- SSv - SvEv
    dE_dt   <- SE - EIs - EIa
    dEv_dt  <- SvEv - EvIs - EvIa
    dIs_dt  <- EIs + EvIs - IR*p_symp - IS*p_symp
    dIa_dt  <- EIa + EvIa - IR*q_symp - IS*q_symp
    dR_dt   <- IR
    
    return (list(c(dS_dt, dSv_dt, dE_dt, dEv_dt, dIs_dt, dIa_dt, dR_dt)))  
  })
}

SimSEIR <- function(input){
  simtime <- Get_time(input)
  parms   <- Get_parms(input)
  stocks  <- Get_stocks(input)
  # Run ode
  o1<-data.frame(ode(y=stocks, times=simtime, func = model1,
                     parms=parms, method="euler"))
  
  o0<-data.frame(ode(y=stocks, times=simtime, func = model0,
                     parms=parms, method="euler"))
  
  # Transform o1 to long format to plot
  o1_long <- o1 %>% 
    mutate(S_T  = S_A  + S_B  + S_C  + S_D  + S_E,
           Sv_T = Sv_A + Sv_B + Sv_C + Sv_D + Sv_E,
           E_T  = E_A  + E_B  + E_C  + E_D  + E_E,
           Ev_T = Ev_A + Ev_B + Ev_C + Sv_D + Sv_E,
           Is_T = Is_A + Is_B + Is_C + Is_D + Is_E,
           Ia_T = Ia_A + Ia_B + Ia_C + Ia_D + Ia_E,
           R_T  = R_A  + R_B  + R_C  + R_D  + R_E) %>%
    pivot_longer(S_A:R_T, 
                 names_to = "SEIR", 
                 values_to = "n") %>%
    separate(SEIR, c("SEIR", "Age_group"))
  
  o0_long <- o0 %>% 
    mutate(S_T  = S_A  + S_B  + S_C  + S_D  + S_E,
           Sv_T = Sv_A + Sv_B + Sv_C + Sv_D + Sv_E,
           E_T  = E_A  + E_B  + E_C  + E_D  + E_E,
           Ev_T = Ev_A + Ev_B + Ev_C + Sv_D + Sv_E,
           Is_T = Is_A + Is_B + Is_C + Is_D + Is_E,
           Ia_T = Ia_A + Ia_B + Ia_C + Ia_D + Ia_E,
           R_T  = R_A  + R_B  + R_C  + R_D  + R_E) %>%
    pivot_longer(S_A:R_T, 
                 names_to = "SEIR", 
                 values_to = "n") %>%
    separate(SEIR, c("SEIR", "Age_group"))
  
  return(list(o1_long=o1_long,o0_long=o0_long))
}

Cost_QALY_cal <- function(df=dt){
  # Collapse age group 65-84 and >=85 years and get number of symptomatic infection
  n_Is <- df %>% filter(SEIR=="Is", Age_group!="Total") %>% select(n)
  n_Is <- c(n_Is$n[1:3], n_Is$n[4]+n_Is$n[5])
  
  n_S <- df %>% filter(SEIR=="S", Age_group!="Total") %>% select(n)
  n_S <- c(n_S$n[1:3], n_S$n[4]+n_S$n[5])
  
  n_Sv <- df %>% filter(SEIR=="Sv", Age_group!="Total") %>% select(n)
  n_Sv <- c(n_Sv$n[1:3], n_Sv$n[4]+n_Sv$n[5])
  
  #### Probability parameters ####
  # Ambulatory care P(ambu)
  p_ambu = 0.15
  
  # Hospitalization given infection P(hosp|Is) 
  p_hosp_A = 0.0092
  p_hosp_B = 0.0081
  p_hosp_C = 0.0826
  p_hosp_D = 0.2570
  p_hosp <- c(p_hosp_A,p_hosp_B,p_hosp_C,p_hosp_D)
  
  # ICU admission, given hospitalization P(icu|hosp)
  p_icu_A = 0.1710
  p_icu_B = 0.2380
  p_icu_C = 0.3610
  p_icu_D = 0.3530
  p_icu <- c(p_icu_A,p_icu_B,p_icu_C,p_icu_D)
  
  # Mortality, given hospitalization P(die|hosp)  
  p_die_A = 0.0061
  p_die_B = 0.0890
  p_die_C = 0.0580
  p_die_D = 0.1550
  p_die <- c(p_die_A, p_die_B, p_die_C, p_die_D)
  
  # Pneumonia, given hospitalization P(pneu|hosp)
  p_pneu = 0.790
  
  # ARDS, requiring ventilator use in ICU P(ards|icu)
  p_ards = 0.771
  
  # Hospitalize given ambu P(hosp|ambu)
  p_hosp1_ambu1 <- p_ambu*p_hosp
  
  # Hospitalize given non-ambu P(hosp|ambu)
  p_hosp1_ambu0 <- (1-p_ambu)*p_hosp
  
  # Stay symptomatic infected (IS)
  p_Is <- 1-(p_hosp1_ambu0+p_ambu)
  
  # Stay symptomatic infected (IS)
  p_Is <- 1-(p_hosp1_ambu0+p_ambu)
  
  # Using surgical masks
  p_sx_m  = 0.020
  p_n95_m = 0.345
  p_cl_m  = 0.455
  
  # Side effect due to vaccination: Severe-myocarditis/pericarditis
  p_card = 0.000023
  
  # Side effect due to vaccination: Severe-allergic reaction/anaphyraxis
  p_anap = (0.000003 + 0.000011)/2
  
  #### Estimate number of each stage ####
  # Number of ambulatory care 
  n_ambu1 = n_Is*p_ambu
  
  # Number of non-ambulatory care 
  n_ambu0 = n_Is*(1-p_ambu)
  
  # Number of hospitalization 
  n_hosp = n_Is*p_hosp 
  
  # Number of severe pneumonia, given hospitalization
  n_pneu1 = n_Is*p_hosp*p_pneu
  
  # Number of severe non-Pneumonia, given hospitalization
  n_pneu0 = n_Is*p_hosp*(1-p_pneu)
  
  # Number of ICU, given hospitalization
  n_icu = n_Is*p_hosp*p_icu
  
  # Number of ARDS, given hospitalization and ICU admission
  n_ards = n_Is*p_hosp*p_icu*p_ards
  
  # Number of sepsis, given hospitalization and ICU admission
  n_seps = n_Is*p_hosp*p_icu*(1-p_ards)
  
  # Number of deaths given hospitalization 
  n_die = n_Is*p_hosp*p_die
  
  # Number of survive given hospitalization 
  n_surv = n_Is*p_hosp*(1-p_die)
  
  # Number of Severe-myocarditis/pericarditis
  n_myoc = n_Sv*p_card/2
  n_peri = n_Sv*p_card/2 
  
  # Number of Severe-allergic reaction/anaphyraxis
  n_anap = n_Sv*p_anap 
  
  #### Costs parameters ####
  # Ambulatory care visit USD/day 
  c_ambu = (110.43 + 148.33)/2
  
  # Over counter medications USD/day
  c_med_A = 3.990 # apply to <18 years
  c_med_B = 0.470 # apply to >=18 years
  c_med <- c(c_med_A, c_med_B, c_med_B, c_med_B)
  
  # Hospitalization for severe pneumonia USD/day
  c_pneu1_A = 12877.37
  c_pneu1_B = 10945.96
  c_pneu1_C = 14129.68
  c_pneu1_D = (12632.32+11312.21)/2 # average between 65-84 and >=85 years
  c_pneu1 <- c(c_pneu1_A, c_pneu1_B, c_pneu1_C, c_pneu1_D)
  
  # Hospitalization for severe non-pneumonia USD/day
  c_pneu0 = 7093.13
  
  # Hospitalization for sepsis USD/day
  c_seps_A = 23375.13
  c_seps_B = 45091.74
  c_seps_C = 39896.27
  c_seps_D = (31217.54 + 23375.13)/2 # average between 65-84 and >=85 years
  c_seps <- c(c_seps_A, c_seps_B, c_seps_C, c_seps_D)
  
  # Hospitalization for ARDS USD/day
  c_ards_A = 43621.1
  c_ards_B = 26997.29
  c_ards_C = 20459.9
  c_ards_D = (19280.11 + 17056.54)/2 # average between 65-84 and >=85 years
  c_ards <- c(c_ards_A, c_ards_B, c_ards_C, c_ards_D)
  
  # Vaccination cost
  c_vacc     = 20
  c_vacc_adm = 40
  
  # Surgical face mask
  c_sx_m  = 0.08
  c_n95_m = 0.5
  c_cl_m  = 2.5
  c_cl_m_w = 0.007
  
  # Hospitalization for myocarditis
  c_myoc = 35289.6
  
  # Hospitalization for pericarditis
  c_peri = 16002.76
  
  # Hospitalization for allergic reaction/anaphylaxis
  c_anap = 7753.38
  
  #### Calculate costs ####
  ##!!! Need to reconsider about time unit !!!##
  # Ambulatory care
  d_ambu = 7
  C_ambu = sum((n_ambu1*c_ambu) + (n_ambu1*c_med))
  
  # Hospitalization, not admitted to ICU
  d_hosp_A = 3
  d_hosp_C = 4
  d_hosp_D = 6
  d_hosp <- c(d_hosp_A, d_hosp_A, d_hosp_C, d_hosp_D)
  C_hosp = (n_pneu1 * c_pneu1) + (n_pneu0 * c_pneu0)
  
  # Hospitalization, admitted to ICU (ARDS + SEPSIS)
  d_icu = 9
  C_icu = (n_seps * c_seps) + (n_ards * c_ards)
  
  # Total treatment costs 
  C_tx_total = sum(C_ambu, C_hosp, C_icu)
  
  # Hospitalization, myocarditis
  d_myoc = 5.9
  C_myoc = n_myoc*c_myoc
  
  # Hospitalization, pericarditis
  d_peri = 4.8
  C_peri = n_peri*c_peri
  
  # Hospitalization, allergic reaction/anaphylaxis
  d_anap = 2.3
  C_anap = n_anap*c_anap
  
  # Vaccination cost
  C_vacc = sum((c_vacc+c_vacc_adm)*n_Sv*3) + sum(C_myoc, C_peri, C_anap)
  
  # Face mask cost
  # May need to reflex face mask policy for each age-group
  n_cl_m = 2
  n_dis_m = 0.515
  C_fm = sum((n_S + n_Sv)*p_sx_m*c_sx_m/n_dis_m + 
               (n_S + n_Sv)*p_n95_m*c_n95_m/n_dis_m +
               (n_S + n_Sv)*p_cl_m*n_cl_m*c_cl_m + (n_S + n_Sv)*p_cl_m*n_cl_m*c_cl_m_w)  
  
  # Total medical cost (vaccine + treatment)
  Cost_med = C_vacc + C_tx_total 
  
  #### Utility weights ####
  # Healthy QALY weights
  qaly_A = 1
  qaly_B = 0.92
  qaly_D = 0.84
  qaly <- c(qaly_A, qaly_B, qaly_B, qaly_D)
  
  # Mild non-specific symptoms QALY weight
  qaly_Is = 0.648
  
  # Hospitalized non-pneumonia QALY weight
  qaly_pneu0 = 0.514
  
  # Pneumonia QALY weight
  qaly_pneu1 = 0.496
  
  # Sepsis QALY weight
  qaly_seps = 0.467
  
  # ARDS QALY weight
  qaly_ards = 0.100
  
  #### Calculate QALY ####
  # Healthy individuals
  QALY_S = qaly*(n_S + n_Sv)/365.5
  
  # Mild non-specific symptoms
  QALY_Is = qaly_Is*(n_Is-n_hosp)/365.5
  
  # Hospitalized non-pneumonia 
  QALY_pneu0 = (qaly_pneu0*n_pneu0)/365.5
  
  # Pneumonia QALY weight
  QALY_pneu1 = (qaly_pneu1*n_pneu1)/365.5
  
  # Sepsis QALY weight
  QALY_seps = (qaly_seps*n_seps)/365.5
  
  # ARDS QALY weight
  QALY_ards = (qaly_ards*n_ards)/365.5
  
  # Total QALY
  QALY_total = sum(QALY_S,QALY_Is,QALY_pneu0,QALY_pneu1,QALY_seps,QALY_ards)
  
  return(list(QALY_total=QALY_total, Cost_med=Cost_med, C_fm=C_fm))
}

make_icer_tbl <- function(costs0, costs1, qalys0, qalys1){
  library(kableExtra)
  library(gridExtra)
  # Format
  format_costs <- function(x) {
    formatC(x/1000000, format = "f", digits = 2)
  }
  format_qalys <- function(x) {
    formatC(x/1000000, format = "f", digits = 2)
  }
  # Computations
  total_costs0 <- sum(costs0)
  total_costs1 <- sum(costs1)
  total_qalys0 <- sum(qalys0)
  total_qalys1 <- sum(qalys1)
  incr_total_costs <- total_costs1 - total_costs0
  inc_total_qalys <- total_qalys1 - total_qalys0
  icer <- incr_total_costs/inc_total_qalys
  
  # Make table
  tibble(
    `Strategy` = c("No face mask use", "Face mask use"),
    `Costs per mil` = c(total_costs0, total_costs1) %>%
      format_costs(), 
    `QALYs per mil` = c(total_qalys0, total_qalys1) %>%
      format_qalys(),
    `Incremental costs per mil` = c("--", incr_total_costs %>% 
                                      format_costs()),
    `Incremental QALYs per mil` = c("--", inc_total_qalys %>% 
                                      format_qalys()),
    `ICER` = c("--", round(icer,2))
  ) %>%
    kable() %>%
    kable_styling() %>%
    footnote(general = "Preliminary results of Thailand COVID-19 SEIR for ICER, (USD)",
             footnote_as_chunk = TRUE)
}

#### Define server logic ####
shinyServer(function(input, output) {

    output$distPlot <- renderPlot({

      p1_total <- ggplot(data=SimSEIR(input)$o1_long %>% filter(Age_group=="T"),  
                         aes(x=time, y=n, color=SEIR)) +
        geom_line(size=1)+
        coord_cartesian(xlim = c(input$xlim[1],input$xlim[2]), 
                        ylim = c(input$ylim[1],input$ylim[2])) +
        scale_y_continuous(labels = comma) +
        scale_color_brewer(type="Qualitative",
                           palette = "Paired") +
        labs(y="Number of population") +
        theme_bw()
      
      p1_total

    })
    
    output$distPlot2 <- renderPlot({
      
      p1 <- ggplot(data=SimSEIR(input)$o1_long,  
                   aes(x=time, y=n, color=SEIR)) +
        geom_line(size=1)+
        coord_cartesian(xlim = c(input$xlim[1],input$xlim[2]), 
                        ylim = c(input$ylim[1],input$ylim[2])) +
        scale_y_continuous(labels = comma) +
        scale_color_brewer(type="Qualitative",
                           palette = "Paired") +
        labs(y="Number of population") +
        facet_wrap(vars(Age_group), ncol = 2, 
                   labeller = as_labeller(c('A'="0-17 years", 
                                            'B'="18-44 years",
                                            'C'="45-64 years",
                                            'D'="65-84 years",
                                            'E'=">=85 years",
                                            'T'="Total"))) +
        theme_bw()
      
      p1
      
    })
    
    datasetInput <- reactive({
      dat <- SimSEIR(input)$o1_long
      dat <- dat %>% 
        filter(time>=input$time_value_min, 
               time<=input$time_value_max, 
               SEIR%in%input$SEIR_select) %>%
        mutate(n=round(n,0),
               Age_group = factor(Age_group, 
                                  labels = c(
                                    "0-17 years", 
                                    "18-44 years",
                                    "45-64 years",
                                    "65-84 years",
                                    ">=85 years",
                                    "Total"))
               
        )
      dat
    })
    
    output$Table1 <- renderDataTable({
      datasetInput()
    })
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste(input$dataset, ".csv", sep = "")
      },
      content = function(file) {
        write.csv(datasetInput(), file, row.names = FALSE)
      })
    
    ICER_dat <- reactive({
      df1 <- SimSEIR(input)$o1_long
      df0 <- SimSEIR(input)$o0_long 
      
      ICER <- NULL
      Incre_costs <- NULL
      Incre_qaly <- NULL
      df0_temp <- NULL
      df1_temp <- NULL
      for (i in seq(10,600,10)) {
        df1_temp <- df1%>%filter(time==i)
        df0_temp <- df0%>%filter(time==i)
        Cost1 = (Cost_QALY_cal(df1_temp)$Cost_med + Cost_QALY_cal(df1_temp)$C_fm) 
        Cost0 = Cost_QALY_cal(df0_temp)$Cost_med
        QALY1 = Cost_QALY_cal(df1_temp)$QALY_total 
        QALY0 = Cost_QALY_cal(df0_temp)$QALY_total
        Incre_costs[i/10] = Cost1 - Cost0
        Incre_qaly[i/10]  = QALY1-QALY0
        ICER[i/10]        = (Cost1-Cost0)/(QALY1-QALY0)
      }
      
      ICER_dat <- tibble(ICER=ICER, 
                         Incre_costs=Incre_costs,
                         Incre_qaly=Incre_qaly,
                         time=seq(10,600,10))
      
      ICER_dat
    })
    
    output$Table2 <- renderDataTable({
      ICER_dat()
    })
    
    output$distPlot3 <- renderPlot({
      
      p3 <- ggplot(data=ICER_dat(),  
                   aes(x=time, y=ICER)) +
        geom_line(size=1, color="#253494")+
        geom_point(size=4, fill="#f03b20", shape = 21) +
        geom_hline(yintercept = 50000, size=1, linetype=2) +
        coord_cartesian(xlim = c(input$xlim[1],input$xlim[2]),
                        ylim = c(-200000,200000)) +
        scale_y_continuous(labels = comma) +
        labs(title = "ICER", 
             y     = "USD/QALY") +
        theme_bw()
      
      p3
    })
    
    output$distPlot4 <- renderPlot({
      
      p4 <- ggplot(data=ICER_dat(),  
                   aes(x=time, y=Incre_costs)) +
        geom_line(size=1, color="#253494")+
        geom_point(size=4, fill="#f03b20", shape = 21) +
        coord_cartesian(xlim = c(input$xlim[1],input$xlim[2])) +
        scale_y_continuous(labels = comma) +
        labs(title = "Incremental costs",
             y     = "USD") +
        theme_bw()
      
      p4
    })
    
    output$distPlot5 <- renderPlot({
      
      p5 <- ggplot(data=ICER_dat(),  
                   aes(x=time, y=Incre_qaly)) +
        geom_line(size=1, color="#253494")+
        geom_point(size=4, fill="#f03b20", shape = 21) +
        coord_cartesian(xlim = c(input$xlim[1],input$xlim[2])) +
        scale_y_continuous(labels = comma) +
        labs(title = "Incremental QALYs",
             y     = "QALYs") +
        theme_bw()
      
      p5
    })

})
