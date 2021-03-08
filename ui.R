#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(shiny)
library(tidyverse)
library(magrittr)
library(shinythemes)
library(shinydashboard)
library(rhandsontable)
library(shinyjs)
library(shinyWidgets)
library(bsplus)
library(htmltools)
library(shinyglide)
library(DT)

source("helpers.R") # Load all the code needed to show feedback on a button click

##################################################################################
## header
##################################################################################


header <- dashboardHeader(
  title = "MMDIT 1.0"#, 
  # tags$li(class = "dropdown", tags$img(class= "hsc-logo",
  #                                       src='hsc_white.png',
  #                                       width = "50",
  #                                       height = "50",
  #                                       style = "margin-right: 20px;"
  #                                       ))
                    )
                            


##################################################################################
## sidebar
##################################################################################
sidebar <- dashboardSidebar(
    width = 250,
    sidebarMenu(id= "tabs",
        menuItem(text = "Home",icon = icon("home"), tabName = "home_tab_ID"),
        menuItem(text = "Data import",icon = icon("file-import"), tabName = "data_import_ID"),
        menuItem(text = "Mixture Deconvolution",icon = icon("bars"), tabName = "continuous_analysis_method_ID"),
        menuItem(text = "Mixture Analysis", icon = icon("bars"), tabName = "semi-continuous_analysis_method_ID"),
        menuItem(text = "User guide", icon = icon("book"), tabName = "user_guide_tab_ID"),
        menuItem(text = "Github code repository", href="https://github.com/SammedMandape/MMDIT_UI", icon = icon("github")), #tabName = "code_repo_tab_ID"),
        tags$br(),tags$br(),
        actionButton("start_over_ID", "Start over", styleclass = "primary", size = "small")
    )
    
)

    
##################################################################################
## body
##################################################################################   
body <- dashboardBody(
  useShinyjs(),
  tags$head(
    tags$link(href="mydatastyles.css", rel="stylesheet", type="text/css")
  ),
  tags$head(includeHTML(("ggl_ana.html"))),
    tabItems(
      #########################################################################################################
      ## Home Tab
      #########################################################################################################
      tabItem(tabName = "home_tab_ID",
              fluidRow(
                box(width = 12, height = "auto",
                            tags$div(class="landing-wrapper",
                                    
                                    # child element 1: images
                                    tags$div(class="landing-block background-content",
                                             # images - top -> bottom, left -> right
                                             tags$img(src="MMDIT.png"))),
                             
                            tags$div(style=("width:400px;margin:auto;"),       
                                             shiny::actionButton("lets_begin_ID","Begin Analysis"),
                                     ),
                            tags$div(class="publication-details",
                                    tags$div(class="paper-text", 
                                    tags$p("Mitochondrial mixture database and interpretation tool (MMDIT) is an open-source, interactive software
                                    for the probabilistic genotyping of mitochondrial DNA mixtures based on complete mitochondrial genomes (mtGenomes). 
                                    MMDIT can perform both 'Mixture Deconvolution' and 'Mixture Analysis'. This tool is described in detail in the paper",
                                           tags$a("MMDIT tool.",href="www.google.com"), "Funding: This work was supported by Award Number 
                                           2017-DN-BX-0134 by the National Institute of Justice Office of Justice Programs, U.S. 
                                           Department of Justice. The opinions, findings, and conclusions or recommendations expressed 
                                           are those of the authors and do not necessarily reflect those of the U.S.  Department of Justice."), 
                                    tags$p("If you use this tool, please consider citing:"),
                                    tags$li("Mandape et al: MMDIT paper"), tags$li("Smart, U.; Cihlar, J.C.;
                                                                                    Mandape, S.N.; Muenzler, M.; King,
                                                                                    J.L.; Budowle, B.;Woerner, A.E. A
                                                                                    Continuous Statistical Phasing
                                                                                    Framework for the Analysis of
                                                                                    Forensic Mitochondrial DNA
                                                                                    Mixtures. Genes 2021, 12, 128.
                                                                                    https://doi.org/10.3390/
                                                                                    genes12020128"), 
                                    tags$li("Crysup, B.;Woerner, A.E; King, J.L; Budowle, B. Graph Algorithms for Mixture Interpretation. Genes 2021, 12, 185. https://doi.org/
                                              10.3390/genes12020185")
                                    )),
                    fluidRow(
                      tags$hr(class="line-break"),
                      column(width = 6, 
                             tags$h3("Download Sample Data", style = "text-align:right;margin-right:10%"),
                             tags$p("Download the sample data and follow user guide for instructions on running the tool", 
                                    style = "text-align:right;margin-right:10%"),
                             downloadButton("download_data_ID", label = "Sample Data"),
                             tags$br(),tags$br(),tags$br(),
                             actionButton("user_guide_ID", label = "User Guide"),
                      
                    ),
                    column(
                      width = 6, style = "border-left:1px solid gray",
                      tags$h3("Contact", style = "margin-left:10%"),
                      tags$div( style = "margin-left:10%",
                      tags$p("For questions, comments or suggestions about the tool, please contact,"),
                      tags$p(
                        tags$strong("Sammed Mandape"),
                        tags$br(), "Bioinformatician",
                        tags$br(), tags$a(href = "mailto: sammed.mandape@unthsc.edu", "sammed.mandape@unthsc.edu")
                      ),
                      tags$p(
                        tags$strong("Center for Human Identification"),
                        tags$br(), "3500 Camp Bowie Blvd.",
                        tags$br(), "Fort Worth, Tx 76107"
                      )
                      
                    ))
                    
                  ), ### fluidRow for columns end here
                                    
                tags$div(id = "footers",
                  tags$footer(class ="green-footer",
                    tags$a(
                      href = "https://www.unthsc.edu/graduate-school-of-biomedical-sciences/laboratory-faculty-and-staff/",
                      target = "_blank",
                      tags$img(
                        src = "UNTCHI_RDU.png",
                        width = "75",
                        height = "75",
                        align = "left",
                        title = "UNTCHI R&D Unit",
                        style = "margin: 10px 10px"
                      )
                    ), 
                    
                   tags$a(
                      href = "https://www.untchi.org/",
                      target = "_blank",
                      tags$img(
                        src = "CHI_2018.png",
                        width = "75",
                        height = "75",
                        align = "right",
                        title = "Center for Human Identification Website",
                        style = "margin: 10px 10px"
                      )
                    ),
                    
                    tags$div(class = "social-media-unthsc",
                      
                      tags$br(),
                      
                      tags$a(
                        href ="https://twitter.com/UNTCHI_RDU",
                        target = "_blank",
                        tags$img(src = "twitter.png",
                                 title = "UNTCHI RDU on Twitter",
                                 width = "30",
                                 height = "30",
                                 class = "displayed",
                                 style = "margin: 10px 10px"
                        )
                      ),
                      
                      tags$a(
                        href ="https://www.instagram.com/untchi_rdu/",
                        target = "_blank",
                        tags$img(src = "instagram.png",
                                 title = "UNTCHI RDU on Instagram",
                                 width = "30",
                                 height = "30",
                                 class = "displayed",
                                 style = "margin: 10px 10px"
                        )
                      ),
                      
                      tags$a(
                        href = "https://www.facebook.com/UNTCenterForHumanIdentification",
                        target = "_blank",
                        tags$img(src = "facebook.png",
                                 title = "UNTCHI RDU on Facebook",
                                 width = "30",
                                 height = "30",
                                 class = "displayed",
                                 style = "margin: 10px 10px"
                        )
                      )
                    )
                  ),
                  
                  ##############################################
                  #Apply secondary footer    
                  tags$footer(
                    tags$br(),
                    
                    tags$div(
                      id = "secondary-footer-docs"
                    ),
                    
                    tags$br(), 
                    
                    style = "background-color: #253746;
                    color: #FFFFFF;"
                  )
      ) ##div-footer ends here
      ) ##box ends here
      ) ##fluid row ends here
                ), ##tabitem ends here

      
      #########################################################################################################
      ## Data import tab
      #########################################################################################################
        tabItem(tabName = "data_import_ID",
                fluidRow(
                    box(title = tags$b("STEP 1 : File Import"), width = 12, height = "auto",
                        column(width = 6,
                               fileInput("select_ss_data_ID",
                                         label = "Upload one or more single source data",
                                         multiple = TRUE,
                                         accept = "text/plain"),
                               fileInput("select_mix_data_ID",
                                         label = "Upload mixture data *",
                                         accept = "text/plain"),
                               #fluidRow 
                                tags$div(
                               actionButton("clear_ID", label = "Clear selection",class="btn-info"),
                               )
                          ),
                        column(width=6,
                               verbatimTextOutput("files_selected", placeholder = T),
                               tags$head(tags$style("#files_selected{font-weight:bold;
                                                             overflow-y:scroll;
                                                             max-height: 200px;
                                                             color:blue; margin-top:24px;
                                                             white-space: pre-wrap;
                                                  }"))

                        )
                       ) 
                ),
                
                fluidRow(
                    box(title = tags$b("STEP 2 : Preprocessing"), width = 12, height = "auto",
                        column(width = 6,
                        tags$div(actionButton(
                          "analyze_data_ID", label = "Analyze data",
                        ),
                      
                        )
                    )
                )),
                
                
                fluidRow(
                  box(title = tags$b("STEP 3 : Preprocessing Results"), width = 12, height = "auto", 
                    tabsetPanel(id = "primary_analysis_ID",
                           tabPanel("SNV analysis",
                                    actionButton(inputId = "next_to_indel_analysis_ID", label = "Next"), #),
                                    tags$br(),
                                    rHandsontableOutput("empop_variant_input_snp"),
                                    tags$head(tags$style(
                                      "#empop_variant_input_snp{float:center;
                                      text-align:center;
                                      margin-left:auto;
                                      margin-right:auto;
                                      }"
                                    ))
                                    ),
                           tabPanel("Indel analysis",
                                    withBusyIndicatorUI(
                                      actionButton(inputId = "save_df_ID", label = "Save", class="btn-info")
                                      ),
                                    rHandsontableOutput("empop_variant_input_indel"),
                                    tags$head(tags$style(
                                      "#empop_variant_input_indel{float:center;
                                      text-align:center;
                                      margin-left:auto;
                                      margin-right:auto;
                                      }"
                                    ))

                                    ),
                           
                           tabPanel("Select Inclusion/Exclusion list",
                                    column(8, style = "padding-left: 0px;
                                                      margin-top: 10px;",
                                    fluidRow(style = "margin-right: -10px;
                                             margin-left: -10px;", 
                                      column(6,style="margin-top:10px;",
                                          withBusyIndicatorUI(actionButton(
                                                            "load_MMDIT_ID",
                                                            "Load MMDIT database",
                                                            style="margin-top:0px;margin-right:0px;margin-bottom:10px;"
                                                            ))),
                                      column(6,style="margin-top:10px;",
                                          withBusyIndicatorUI(actionButton("done_run_backend_MMDIT_ID",
                                                                           "Finalize list", style="margin-top:0px;margin-right:20px;")))
                                      
                                    ),
                                   fluidRow(box(title = "Choose Inclusion/Exclusion list", status = "warning",
                                        div(style = "text-align: left;",
                                            radioGroupButtons(
                                              inputId = "Id072",
                                              label = "",
                                              choices = c("Choose one option to input inclusion / exclusion list",
                                                          "Choose from Precision ID Kit",
                                                          "Upload bed file of regions to include / exclude"),
                                              checkIcon = list(
                                                yes = tags$i(class = "fa fa-check-square", 
                                                             style = "color: steelblue"),
                                                no = tags$i(class = "fa fa-square-o", 
                                                            style = "color: steelblue")),
                                              direction = "vertical",
                                              justified = TRUE
                                            ))
                                        ),
                                        box(title = "Select populations", status = "warning",
                                            pickerInput(
                                              inputId = "myPicker_accor_ID", 
                                              label = "", 
                                              choices = "", 
                                              options = list(
                                                `actions-box` = TRUE, 
                                                size = 4,
                                                `selected-text-format` = "count > 6"
                                              ), 
                                              multiple = TRUE
                                            ) %>% 
                                              bsplus::shinyInput_label_embed(
                                                shiny::icon("info-circle") %>%
                                                  bs_embed_tooltip("Choose at least 2 populations. All populations
                                                                   are selected by default.")
                                              )
                                            
                                            )
                                        
                                        )
                                   
                                    )
                                    )
                    )
                )
                ),
                fluidRow(
                  box(title = tags$b("STEP 4 : Analyze"), tags$p("Choose between one of the following methods to analyze data."), width = 12, height = "auto", style = "text-align:left;",
                      column(width = 12,
                             tags$div(actionButton(
                               "continuous_btn_ID", label = "Mixture Deconvolution",
                             ),
                             actionButton(
                               "semi_continuous_btn_ID", label = "Mixture Analysis",
                             ),
                             
                             )
                      )
                  ))
                ),
      #########################################################################################################
      ## Mixture Analysis
      #########################################################################################################
        tabItem(tabName = "semi-continuous_analysis_method_ID",
                fluidRow(
                  box(title = tags$b("Mixture Analysis"), width = 12, height = "auto",
                    column(
                    selectInput(
                      "select_mixture_type_ID",
                      label = "Enter the number of contributors in mixture",
                      choices = c('2-persons mixture', '3-persons mixture' 
                      )
                    ),
                    width = 6
                  ),
                  tags$br(),
                  withBusyIndicatorUI(actionButton("generate_mixture_statistics_ID", "Generate Mixture Statistics"))
                  )
                ),
                fluidRow(
                  box(title = tags$b("RESULTS"), width = 12, height = "auto",
                      
                      div(
                        tags$b("Random Man Not Excluded (RMNE) Statistics"),
                        rHandsontableOutput("mixStat_rmne_rhotOut")
                      ),
                      tags$br(),
                      div(
                        tags$b("Likelihood Ratios (LR) Statistics"),
                        rHandsontableOutput("mixStat_lr_rhotOut")
                      )
                      
                      
                      
                      )
                )
                
          
        ),
      #########################################################################################################
      ## Mixture Deconvolution
      #########################################################################################################
        tabItem(tabName = "continuous_analysis_method_ID",
                fluidRow(
                  box(title = tags$b("Mixture Deconvolution"), width = 12, height = "auto",
                      fluidRow(
                        column(width = 6,
                              fileInput("select_var_excel_quantition_ID",
                               label = "Upload excel file with quantitative data",
                               accept = ".xlsx"
                               ))),
                      fluidRow(
                        column(width = 3,
                               textInput("cont_edit_dist_input_ID", "Graph edit distance (GED)", value = 4)),
                        column(width = 3,
                               textInput("cont_numMCMC_input_ID","Number of Markov Chain Monte Carlo (MCMC)", value = 3000)),
                        column(width = 3, # TODO add litte info icon to mention about what happens if less than 15
                               textInput("cont_panel_size_input_ID", "Panel size", value = 25)),
                        column(width = 3,
                               numericInput("cont_numPerson_inMix_input_ID","Number of contributors",value = 2, min = 1,max = 5)),
                        ),
                      fluidRow(
                        column(width = 3,
                               textInput("cont_readCount_norm_ID","Read count normalization value", value = 100)),
                        column(width = 3,
                               textInput("cont_miscopying_rate_input_ID","Miscopying rate", value = 0.01))
                      ),
                      withBusyIndicatorUI(actionButton("genereate_cont_run_deploid", label = "Analyze",class="btn-info")),
                      tags$style("#genereate_cont_run_deploid{float:left;width:auto}")
                      )
                ),
                fluidRow(
                  box(title = tags$b("RESULTS"), width = 12, height = "auto",
                      tabsetPanel(id="cont_results_box_ID",
                                  tabPanel("Statistics",
                                           tags$br(),
                                           div(
                                             tags$b("Proportion of Mixtures"),
                                             rHandsontableOutput("cont_deploid_mixProp_rhotOut")
                                           ),
                                           tags$br(),
                                           div(
                                             tags$b("Estimated Haplotypes"),
                                             rHandsontableOutput("cont_deploid_estimatedHap_rhotOut")
                                           )
                                           
                                           
                                           
                                  ),
                                  tabPanel("Trace plot",
                                           plotOutput("cont_tace_plot_ID")
                                           ),
                                  tabPanel("Calculate Random Match Probability",
                                           tags$br(),
                                           rHandsontableOutput("cont_rmp_rhot_ID")
                                           )
                      )
                  )
                )
                ),
      ###################################################################################################################
      ## USER Guide Tab
      ###################################################################################################################
      tabItem(tabName = "user_guide_tab_ID",
              fluidRow(
                box(width=12,
                    uiOutput("samplepdf")
              )))
    ), ###Tabitems end here
  use_bs_tooltip(),
    
)


dashboardPage(
  
  header,
  sidebar,
  body
  
)

