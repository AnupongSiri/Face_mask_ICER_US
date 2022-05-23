#### Moving time for loop in which improve efficiency####
rm(list = ls())
library(tidyverse)
library(triangulr)
library(mc2d)
library(ggplot2)
library(scales)
library(data.table)

#### Functions ####
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

estGammaParams <- function(mu, std) {
  shape <- (mu^2)/(std^2)
  scale <- (std^2)/mu
  return(params = list(shape = shape, scale = scale))
}

#### Cost and QALY calculate function
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
  p_ambu = rtri(n=1,min=0.06,max=0.26,mode=0.15)
  
  # Hospitalization given infection P(hosp|Is) 
  p_hosp_A = rpert(1,min=0.0081,max=0.0101,mode=0.0092)
  p_hosp_B = rpert(1,min=0.0073,max=0.0089,mode=0.0081)
  p_hosp_C = rpert(1,min=0.0744,max=0.9090,mode=0.0826)
  p_hosp_D = rpert(1,min=0.2314,max=0.2828,mode=0.2570)
  p_hosp <- c(p_hosp_A,p_hosp_B,p_hosp_C,p_hosp_D)
  
  # ICU admission, given hospitalization P(icu|hosp)
  p_icu_A = rpert(1,min=0.1540,max=0.1881,mode=0.1710)
  p_icu_B = rpert(1,min=0.2140,max=0.2620,mode=0.2380)
  p_icu_C = rpert(1,min=0.3250,max=0.3970,mode=0.3610)
  p_icu_D = rpert(1,min=0.3180,max=0.3880,mode=0.3530)
  p_icu <- c(p_icu_A,p_icu_B,p_icu_C,p_icu_D)
  
  # Mortality, given hospitalization P(die|hosp)  
  p_die_A = rpert(1,min=0.0055,max=0.0067,mode=0.0061)
  p_die_B = rpert(1,min=0.0801,max=0.0979,mode=0.0890)
  p_die_C = rpert(1,min=0.0520,max=0.0635,mode=0.0580)
  p_die_D = rpert(1,min=0.1392,max=0.1702,mode=0.1550)
  p_die <- c(p_die_A, p_die_B, p_die_C, p_die_D)
  
  # Pneumonia, given hospitalization P(pneu|hosp)
  #p_pneu = rpert(1,min=0.7110,max=0.8690,mode=0.7900)
  p_pneu = rbeta(1,37.4825,9.963703)
  
  # ARDS, requiring ventilator use in ICU P(ards|icu)
  #p_ards = rpert(1,min=0.771-0.053,max=0.771+0.053,mode=0.771)
  p_ards = rbeta(1,47.69001,14.16474)
  
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
  p_card = runif(1,min=0.0000156,max=0.0000270)
  
  # Side effect due to vaccination: Severe-allergic reaction/anaphyraxis
  p_anap = runif(1,min=0.0000030,max=0.0000110)
  
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
  c_ambu = runif(1,min=110.43,max=148.33)
  
  # Over counter medications USD/day
  c_med_A = rgamma(1,3.99) # apply to <18 years
  c_med_B = rgamma(1,0.47) # apply to >=18 years
  c_med <- c(c_med_A, c_med_B, c_med_B, c_med_B)
  
  # Hospitalization for severe pneumonia USD/day
  c_pneu1_A = rgamma(1,shape=72.91697,scale=176.6032)
  c_pneu1_B = rgamma(1,shape=109.7047,scale=99.77658)
  c_pneu1_C = rgamma(1,shape=130.1039,scale=108.6030)
  c_pneu1_D = rgamma(1,shape=577.155,scale=20.74359) # average between 65-84 and >=85 years
  c_pneu1 <- c(c_pneu1_A, c_pneu1_B, c_pneu1_C, c_pneu1_D)
  
  # Hospitalization for severe non-pneumonia USD/day
  c_pneu0 = rgamma(1,shape=35.95122,scale=197.2987)
  
  # Hospitalization for sepsis USD/day
  c_seps_A = rgamma(1,shape=157.7108,scale=148.2152)
  c_seps_B = rgamma(1,shape=70.18464,scale=642.4731)
  c_seps_C = rgamma(1,shape=214.3383,scale=186.1369)
  c_seps_D = rgamma(1,shape=285.8038,scale=95.50724) # average between 65-84 and >=85 years
  c_seps <- c(c_seps_A, c_seps_B, c_seps_C, c_seps_D)
  
  # Hospitalization for ARDS USD/day
  c_ards_A = rgamma(1,shape=107.9214,scale=404.1931)
  c_ards_B = rgamma(1,shape=300.0305,scale=89.98181)
  c_ards_C = rgamma(1,shape=2031.648,scale=10.07059)
  c_ards_D = rgamma(1,shape=1111.702,scale=16.34281) # average between 65-84 and >=85 years
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
  c_myoc = rgamma(1,shape=252.1578,scale=139.9504)
  
  # Hospitalization for pericarditis
  c_peri = rgamma(1,shape=3020.831,scale=5.29747)
  
  # Hospitalization for allergic reaction/anaphylaxis
  c_anap = rtri(n=1,min=6774.80,max=8280.31,mode=7753.38)
  
  #### Calculate costs ####
  # Ambulatory care
  d_ambu = 7
  C_ambu = sum((n_ambu1*c_ambu) + (n_ambu1*c_med))
  
  # Hospitalization, not admitted to ICU
  d_hosp_A = rpert(1,min=2,max=5,mode=3)
  d_hosp_C = rpert(1,min=2,max=7,mode=4)
  d_hosp_D = rpert(1,min=3,max=10,mode=6)
  d_hosp <- c(d_hosp_A, d_hosp_A, d_hosp_C, d_hosp_D)
  C_hosp = (n_pneu1 * c_pneu1) + (n_pneu0 * c_pneu0)
  
  # Hospitalization, admitted to ICU (ARDS + SEPSIS)
  d_icu = rpert(1,min=4,max=17,mode=9)
  C_icu = (n_seps * c_seps) + (n_ards * c_ards)
  
  # Total treatment costs 
  C_tx_total = sum(C_ambu, C_hosp, C_icu)
  
  # Hospitalization, myocarditis
  d_myoc = rgamma(1,shape=444.0051,scale=0.01328814)
  C_myoc = n_myoc*c_myoc
  
  # Hospitalization, pericarditis
  d_peri = rgamma(1,shape=6400,scale=0.00075)
  C_peri = n_peri*c_peri
  
  # Hospitalization, allergic reaction/anaphylaxis
  d_anap = rgamma(1,shape=241.0256,scale=0.009542554)
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
  qaly_Is = rbeta(1,13.28415,7.216083)
  
  # Hospitalized non-pneumonia QALY weight
  qaly_pneu0 = rbeta(1,15.69598,14.84095)
  
  # Pneumonia QALY weight
  qaly_pneu1 = rbeta(1,3.794383,3.855583)
  
  # Sepsis QALY weight
  qaly_seps = rbeta(1,3.120699,3.56174)
  
  # ARDS QALY weight
  qaly_ards = rtri(1,min=0.08,max=0.15,mode=0.10)
  
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


ICER_MonteC <- function(dt1=dt1,dt0=dt0){
  ICER <- NULL
  for(i in 1:200){
    Cost1 = (Cost_QALY_cal(dt1)$Cost_med + Cost_QALY_cal(dt1)$C_fm) 
    Cost0 = Cost_QALY_cal(dt0)$Cost_med
    QALY1 = Cost_QALY_cal(dt1)$QALY_total 
    QALY0 = Cost_QALY_cal(dt0)$QALY_total
    ICER[i] <- (Cost1-Cost0)/(QALY1-QALY0)
  }
  return(ICER)
}




#### Comparison dataset ####
df1<-read.csv("./Mask_ICER_output/Export_Is_S_Sv_numbers_wFM_US_r1350vf0210Marchcases.csv")
df0<-read.csv("./Mask_ICER_output/Export_Is_S_Sv_numbers_woFM_US_r1350vf0210Marchcases.csv")

df1 <- data.table(df1)
df0 <- data.table(df0)

ICER_t <- NULL
Incre_costs_t <- NULL
Incre_qaly_t <- NULL

start_time <- Sys.time()
for (j in 1:1000) {
  ICER <- NULL
  Incre_costs <- NULL
  Incre_qaly <- NULL
  temp1 <- NULL
  temp0 <- NULL
  
  for (i in seq(10,600,10)) {
    temp1 <- Cost_QALY_cal(df1[time==i])
    temp0 <- Cost_QALY_cal(df0[time==i])
    Cost1 = temp1$Cost_med + temp1$C_fm
    Cost0 = temp0$Cost_med
    QALY1 = temp1$QALY_total 
    QALY0 = temp0$QALY_total
    ICER[i/10] <- (Cost1-Cost0)/(QALY1-QALY0)
    Incre_costs[i/10] <- Cost1 - Cost0
    Incre_qaly[i/10] <- QALY1-QALY0
    
  
  }
  
  ICER_t <- c(ICER_t,ICER)
  Incre_costs_t <- c(Incre_costs_t,Incre_costs)
  Incre_qaly_t <- c(Incre_qaly_t,Incre_qaly)
  
}

ICER_m   <- colMeans(matrix(ICER_t, ncol=60, byrow = T))
ICER_l <- apply(matrix(ICER_t, ncol=60, byrow = T), 2, FUN=quantile , probs = 0.025)
ICER_h <- apply(matrix(ICER_t, ncol=60, byrow = T), 2, FUN=quantile , probs = 0.975)

Incre_costs_m   <- colMeans(matrix(Incre_costs_t, ncol=60, byrow = T))
Incre_costs_l <- apply(matrix(Incre_costs_t, ncol=60, byrow = T), 2, FUN=quantile , probs = 0.025)
Incre_costs_h <- apply(matrix(Incre_costs_t, ncol=60, byrow = T), 2, FUN=quantile , probs = 0.975)

Incre_qaly_m   <- colMeans(matrix(Incre_qaly_t, ncol=60, byrow = T))
Incre_qaly_l <- apply(matrix(Incre_qaly_t, ncol=60, byrow = T), 2, FUN=quantile , probs = 0.025)
Incre_qaly_h <- apply(matrix(Incre_qaly_t, ncol=60, byrow = T), 2, FUN=quantile , probs = 0.975)

end_time <- Sys.time()
end_time - start_time #Took 18.27295 mins for n=200; 1.440377 hours for n=1000; improved for loops get 35.84205 mins for n=1000
#Tool 1.15 hours for 600 time steps
ICER_dat <- tibble(time=seq(10,600,10),
                   ICER_m=ICER_m, ICER_l=ICER_l, ICER_h=ICER_h,
                   Incre_costs_m=Incre_costs_m, Incre_costs_l=Incre_costs_l, Incre_costs_h=Incre_costs_h,
                   Incre_qaly_m=Incre_qaly_m, Incre_qaly_l=Incre_qaly_l, Incre_qaly_h=Incre_qaly_h) 

write.csv(ICER_dat,"./Mask_ICER_output/ICER_Monte_1000_US_r1350vf0210Marchcases.csv", row.names = F)

p3 <- ggplot(data=ICER_dat,  
             aes(x=time, y=ICER_m)) +
  geom_ribbon(aes(ymin = ICER_l, ymax = ICER_h), fill="#a6cee3", )+
  geom_line(size=1, color="#253494")+
  geom_point(size=2, fill="#f03b20", shape = 21) +
  geom_hline(yintercept = 50000, size=1, linetype=2) +
  coord_cartesian(xlim = c(0,600),
                  ylim = c(-200000,200000)) +
  scale_y_continuous(labels = comma,
                     breaks = c(-200000,-100000,0,50000,100000,200000)) +
  labs(#title = "ICERs", 
       y = "USD/QALY",
       x = "Days",
       caption = "Monte Carlo simulation based on parameter distributions (n = 1000 times)") +
  theme_bw()

p3

p4 <- ggplot(data=ICER_dat,  
             aes(x=time, y=Incre_costs_m)) +
  geom_ribbon(aes(ymin = Incre_costs_l, ymax = Incre_costs_h), fill="#a6cee3", )+
  geom_line(size=1, color="#253494")+
  geom_point(size=2, fill="#f03b20", shape = 21) +
  geom_hline(yintercept = 0, size=1, linetype=2) +
  coord_cartesian(xlim = c(0,600)) +
  scale_y_continuous(labels = comma) +
  labs(title = "Incremental costs", 
       y     = "USD",
       x = "Days",
       caption = "Monte Carlo simulation based on parameter distributions (n = 1000 times)") +
  theme_bw()

p4

p5 <- ggplot(data=ICER_dat,  
             aes(x=time, y=Incre_qaly_m)) +
  geom_ribbon(aes(ymin = Incre_qaly_l, ymax = Incre_qaly_h), fill="#a6cee3", )+
  geom_line(size=1, color="#253494")+
  geom_point(size=2, fill="#f03b20", shape = 21) +
  geom_hline(yintercept = 0, size=1, linetype=2) +
  coord_cartesian(xlim = c(0,600)) +
  scale_y_continuous(labels = comma) +
  labs(title = "Incremental QALYs", 
       y     = "QALYs",
       x = "Days",
       caption = "Monte Carlo simulation based on parameter distributions (n = 1000 times)") +
  theme_bw()

p5


SSv <- df0 %>% filter(SEIR=="S"|SEIR=="Sv") %>%
  pivot_wider(names_from = SEIR, values_from = n) %>%
  mutate(SvS_ratio = Sv/S,
         SvSs_ratio = Sv/(S+Sv))

p_SSv <- ggplot(SSv%>%filter(Age_group=="Total"),aes(x=time, y =SvSs_ratio, color=Age_group)) +
  geom_line()+
  ylim(c(0.6,1))
p_SSv

p_SvS <- ggplot(SSv,aes(x=time, y =SvS_ratio, color=Age_group)) +
  geom_line()
p_SvS

Total_pop <- c(73039150, 117818671, 83323439, 47453305, 6604958, sum(c(73039150, 117818671, 83323439, 47453305, 6604958)))

SSv$Total_pop <- rep(Total_pop,601)

SSv <- SSv %>% mutate(SvT = Sv/Total_pop)

pVc <- ggplot(SSv,aes(x=time, y =SvT, color=Age_group)) +
  geom_line()+
  ylim(c(0,1))
pVc

#### Symptomatic infection plot between intervention ####
df0$FaceMask <- "No"
df1$FaceMask <- "Yes"

df01 <- rbind(df0,df1)
  
pIs <- ggplot(df01%>%filter(SEIR=="Is"), 
              aes(x=time, y=n, 
                  color=Age_group))+
  geom_line(aes(linetype=FaceMask),size = 1) + 
  scale_color_brewer(type="Qualitative",
                     palette = "Accent") +
  scale_y_continuous(labels=comma) +
  scale_x_continuous(limits = c(0,600)) +
  scale_linetype_manual(values=c(4,1)) +
  labs(#title = "Symptomatic infection cases by face mask intervention",
       y="Number of population",
       x="Days",
       color="Age group",
       linetype="Face mask use") +
  theme_bw()

pIs 

