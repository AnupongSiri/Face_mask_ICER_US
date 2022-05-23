################################################################################ 
# This Shiny app creat the system dynamics model and calculate ICER using      #
# mean/median of parameters.                                                   #
# Depends on:                                                                  #
#    01_server.r                                                               # 
# To run in R using command (online by Gist.github):                           #
#    library(shiny); runGist("cd7daf9661d390dd24640840b8737723")               #    
# Author:                                                                      #
#     - Anupong Sirirungreung, <anusiri@g.ucla.edu>                            # 
################################################################################ 

#### Load libraries ####
library(shiny)
library(shinyWidgets)
library(plotly)

#### Define Shiny user interface logic ####
shinyUI(fluidPage(
  
    titlePanel("US COVID-19 SEIR"),

    sidebarLayout(
        sidebarPanel(
          fluidRow(
            column(width=4,
                   setSliderColor(c(rep("#b2df8a", 3)), sliderId=c(8,9,10)),
                   h4(div(HTML("<em>Set clinical parameters...</em>"))),
                   
                   sliderInput("vf_value", "Vaccination force (???): ", 0, 0.08, 0.021, step=0.001, post = " units"),
                   sliderInput("ve_value", "Vaccination effectiveness (0.82-0.89): ", 0, 1, 0.86, step=0.01),
                   sliderInput("q_symp_value", "Probability of asymptomatic infection (0.3-0.5): ", 0, 1, 0.45, step=0.01),
                   sliderInput("p_imm_value", "Probability of developing immunity after infection: ", 0, 1, 0.91, step=0.01),
                   sliderInput("D_lat_value", "Latent period (4.1-7): ", 1, 14, 5.2, step=0.1, post = " days"),
                   sliderInput("D_inf_value", "Infectious period (4-15): ", 1, 30, 7, step=0.1, post = " days"),
                   sliderInput("R0_value", "Reproductive number: ", 0, 2.5, 1.66, step=0.01),
                   
                   hr(),
                   
            ),
            column(width=4,
                   
                   # Initial stocks, number of population stratified by age group
                   h4(div(HTML("<em>No. of symptomatic infected...</em>"))),
                   numericInput("Is_A_value", "Among 0-17 years: ", round(((0.2+0.2+0.3+0.6)/4)*73039150/100000), min = 0, max = 12860865),
                   numericInput("Is_B_value", "Among 18-44 years: ", round(((1.4+2.4+3.3)/3)*117818671/100000), min = 0, max = 25105265),
                   numericInput("Is_C_value", "Among 45-64 years: ", round(((3.3+4.0)/2)*83323439/100000), min = 0, max = 18805172),
                   numericInput("Is_D_value", "Among 65-84 years: ", round((3.8)*47453305/100000), min = 0, max = 7359568),
                   numericInput("Is_E_value", "Among >=85 years: ", round((3.8)*6604958/100000), min = 0, max = 826013),
                   
                   # Initial vaccince coverage proportions by age group
                   h4(div(HTML("<em>Vaccination coverage...</em>"))),
                   sliderInput("vc_A_value", "Among 0-17 years: ", 0, 1, 0, step=0.01),
                   sliderInput("vc_B_value", "Among 18-44 years: ", 0, 1, 0, step=0.01),
                   sliderInput("vc_C_value", "Among 45-64 years: ", 0, 1, 0, step=0.01),
                   sliderInput("vc_D_value", "Among 65-84 years: ", 0, 1, 0, step=0.01),
                   sliderInput("vc_E_value", "Among >=85 years: ", 0, 1, 0, step=0.01),
                   hr(),
                   
            ),
            column(width=4,
                   h4(div(HTML("<em>Graphic parameters...</em>"))),
                   
                   # Graph parameters
                   sliderInput("xlim", "X-axis limit: ", 0, 600, c(0,600), step=1, post=" days"),
                   sliderInput("ylim", "Y-axis limit: ", 0, 400000000, c(0,400000000), step=100000, post=" pop."),
                   hr(),
                   
                   # Min-max time (1 day step) to simulate
                   h4(div(HTML("<em>Time point need to get estimated number...</em>"))),
                   numericInput("time_value_min", "Time min:",  0, min = 0, max = 600),
                   numericInput("time_value_max", "Time max:",  600, min = 0, max = 600),
                   hr(),
                   
                   # Face mask effectiveness (ME) by age group; no face mask use (ME==0)
                   h4(div(HTML("<em>Face mask effectiveness comparison...</em>"))),
                   sliderInput("ME_A_value", "Among 0-17 years: ", 0, 0.4, 0.18, step=0.01),
                   sliderInput("ME_B_value", "Among 18-44 years: ", 0, 0.4, 0.18, step=0.01),
                   sliderInput("ME_C_value", "Among 45-64 years: ", 0, 0.4, 0.18, step=0.01),
                   sliderInput("ME_D_value", "Among 65-84 years: ", 0, 0.4, 0.18, step=0.01),
                   sliderInput("ME_E_value", "Among >=85 years: ", 0, 0.4, 0.18, step=0.01),
                   hr(),
                   
                   selectInput("SEIR_select", "SEIR status:",
                               c("Susceptible (S)" = "S",
                                 "Susceptible and vaccinated (Sv)" = "Sv",
                                 "Exposed (E)" = "E",
                                 "Exposed and vacinated (Ev)" = "E",
                                 "Symptomatic infected (Is)" = "Is",
                                 "Asymptomatic infected (Ia)" = "Ia",
                                 "Recovery (R)" = "R"),
                               selected = c("Is","S","Sv"),
                               multiple = T),
                   downloadButton("downloadData", "Download")
            )
          )
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(type = "tabs",
                      
                      tabPanel("Total population",
                               column(width=12,plotOutput("distPlot"))),
                      
                      tabPanel("Stratified by age group", 
                               plotOutput("distPlot2",
                                          height = "800px")),
                      
                      tabPanel("Download data",
                               dataTableOutput("Table1")),
                      
                      tabPanel("ICER data",
                               fluidRow(column(4,
                                               plotOutput("distPlot3")),
                                        column(4,
                                               plotOutput("distPlot4")),
                                        column(4,
                                               plotOutput("distPlot5"))),
                               dataTableOutput("Table2"))
                      )
            
        )
    )
))
