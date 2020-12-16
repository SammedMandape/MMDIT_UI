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
#library(shinysky)
library(shinyjs)
#library(shinymaterial)
library(shinyWidgets)
library(bsplus)
library(htmltools)
library(shinyglide)
library(DT)

source("helpers.R") # Load all the code needed to show feedback on a button click

tags$head(
  tags$link(href="mydatastyles.css", rel="stylesheet", type="text/css")
) 
##################################################################################
## sidebar
##################################################################################
sidebar <- dashboardSidebar(
    sidebarMenu(
        menuItem(text = "Data import", tabName = "data_import_ID"),
        menuItem(text = "Continuous analysis method", tabName = "continuous_analysis_method_ID"),
        menuItem(text = "Semi-continuous analysis method", tabName = "semi-continuous_analysis_method_ID"),
        tags$br(),tags$br(),
        actionButton("start_over_ID", "Start over", styleclass = "primary", size = "small")
    )
    
)

    
##################################################################################
## body
##################################################################################   
body <- dashboardBody(
  useShinyjs(),
    tabItems(
        tabItem(tabName = "data_import_ID",
                fluidRow(
                #     tags$div(class="divInsteadbox",
                    # wellPanel(class="divInsteadbox",fluidRow(
                    #   tags$div(
                    box(title = tags$b("File Import"), width = 12, height = "auto",
                        column(width = 6,
                               fileInput("select_ss_data_ID",
                                         label = "Upload one or more single source data",
                                         multiple = TRUE,
                                         accept = "text/plain"),
                               fileInput("select_mix_data_ID",
                                         #"empop_variant_file_id",
                                         label = "Upload mixture data",
                                         accept = "text/plain"),
                               #fluidRow 
                               tags$div(actionButton(
                                   "analyze_data_ID", label = "Analyze data",class="btn-info"
                                ),
                               tags$style("#analyze_data_ID{float:right}"),
                               actionButton("clear_ID", label = "Clear selection",class="btn-info"),
                               tags$style("#clear_ID{float:left;width:auto}")
                               #<button type="button">Button</button>
                               )
                          ),
                        column(width=6,
                               #h5("Files uploaded"),
                               # select2Input("select2input1", label = "Files uploaded",
                               #              type= "select", multiple = TRUE
                               #              )
                               #htmlOutput("files_selected")
                               verbatimTextOutput("files_selected", placeholder = T),
                               tags$head(tags$style("#files_selected{font-weight:bold;
                                                             overflow-y:scroll;
                                                             max-height: 200px;
                                                             color:blue; margin-top:24px;
                                                             white-space: pre-wrap;
                                                  }"))

                        )
                       ) 
                ),#),tags$style(".divInsteadbox{border10px solid black;
                  #                   background-color: white;}"),
                fluidRow(
                  #conditionalPanel(condition = "(input.analyze_data_id + 1) % 2 == 0",
                  box(width = 12, height = 600, 
                      #div(style = "height:600px;"),
                    # fluidRow(
                    #   column(width = 6,
                    #     selectInput("select_data_for_analysis_ID","Choose dataset",choices = c("All")),
                    #          # tags$head(tags$style("#select_data_for_analysis_ID{float:center;box-sizing: border-box;
                    #          #     border: 10px solid black;color:red;text-align:center;margin-left:auto;
                    #          #          margin-right:auto;}"))
                    #     )
                    #   ),
                    tabsetPanel(id = "primary_analysis_ID",
                    #tabBox(width = 12, id = "primary_analysis_ID",
                           tabPanel("SNP analysis",
                                    # div(style = "display: inline-block;
                                    #     margin: 10px;
                                    #     box-sizing: border-box;
                                    #     width:auto;",
                                    actionButton(inputId = "next_to_indel_analysis_ID", label = "Next", class="btn-info"), #),
                                    tags$style("#next_to_indel_analysis_ID{
                                      margin:10px;
                                      box-sizing:border-box;
                                    }"),
                                    tags$br(),
                                    rHandsontableOutput("empop_variant_input_snp"),
                                    tags$head(tags$style(
                                      "#empop_variant_input_snp{float:center;
                                      text-align:center;
                                      margin-left:auto;
                                      margin-right:auto;
                                      }"
                                    ))
                                    #DTOutput("empop_variant_input_table")
                                    ),
                           tabPanel("Indel analysis",
                                    # div(style = "display: inline-block; margin: 10px; border: 10px Solid black;",
                                    withBusyIndicatorUI(
                                      actionButton(inputId = "save_df_ID", label = "Save", class="btn-info")
                                      ),#),
                                    tags$style("#save_df_ID{
                                      margin:10px;
                                      box-sizing:border-box;
                                    }"),
                                    rHandsontableOutput("empop_variant_input_indel"),
                                    tags$head(tags$style(
                                      "#empop_variant_input_indel{float:center;
                                      text-align:center;
                                      margin-left:auto;
                                      margin-right:auto;
                                      }"
                                    ))

                                    ),
                           
                           # tags$script("$(function() {
                           #                  $('#incl_excl_accor_ID').on('click', function(x) {
                           #                    Shiny.onInputChange('selected_tab_accordi_ID', x.target.innerText)
                           #                  });
                           #                });"),
                           tabPanel("Select Inclusion/Exclusion list",column(2),
                                    column(8, algin="center",
                                    #div(
                                    fluidRow(style = "margin-right: -10px;
                                             margin-left: -10px;",
                                      div(style = "float:left;
                                      margin: 10px; 
                                      box-sizing: border-box;",
                                          withBusyIndicatorUI(actionButton(
                                                            "load_MMDIT_ID",
                                                            "Load MMDIT database",
                                                            class='btn-info'))
                                          ),
                                      div(style = "float:right;
                                      margin: 10px;
                                      box-sizing: border-box;",
                                          withBusyIndicatorUI(actionButton("done_run_backend_MMDIT_ID",
                                                     "Finalize list", class="btn-info"))
                                      ),
                                    ),
                                    bs_accordion(id = "incl_excl_accor_ID") %>%
                                      bs_set_opts(panel_type = "success", use_heading_link = TRUE) %>%
                                      bs_append(title = "Choose Incl Exc list", 
                                                content = div(style = "text-align: left;",
                                                radioGroupButtons(
                                                  inputId = "Id072",
                                                  label = "",
                                                  choices = c("Choose one option to input inclusion / exclusion list",
                                                                "Choose from Precision ID Kit",
                                                                #"Manually input genomic coordinate intervals",
                                                                "Upload bed file of regions to include / exclude"),
                                                  checkIcon = list(
                                                    yes = tags$i(class = "fa fa-check-square", 
                                                                 style = "color: steelblue"),
                                                    no = tags$i(class = "fa fa-square-o", 
                                                                style = "color: steelblue")),
                                                  direction = "vertical",
                                                  justified = TRUE
                                                ))
                                                  ) %>%
                                      bs_set_opts(panel_type = "success", use_heading_link = TRUE) %>%
                                      bs_append(title= "Select populations",
                                                content = pickerInput(
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
                                                    #icon("info") %>%
                                                    shiny::icon("info-circle") %>%
                                                      bs_embed_tooltip("Select more than one
                                                                       population. Make sure 
                                                                       to load MMDIT database 
                                                                       first. Default is all 
                                                                       populations selected.")
                                                  ))
                                    )
                                    )#,
                    # tabPanel("Testing shiny slige",
                    #          fluidPage(id = "glideTest",
                    #           # tags$style("#glideTest{
                    #           #            background-color: #007BA7;
                    #           # }"),  
                    #          # fluidRow(column(6,align="center",
                    #           glide(
                    #            screen(
                    #              p("This is a very simple shinyglide application."),
                    #              p("Please click on Next to go to the next screen.")
                    #            ),
                    #            screen(
                    #              p("Please choose a value."),
                    #              numericInput("n", "n", value = 10, min = 10)
                    #            ),
                    #            screen(
                    #              p("And here is the result.")#,
                    #              #plotOutput("plot")
                    #            )
                    #          ))#))
                    #          )
                    )
                )
                )
                ),
        tabItem(tabName = "semi-continuous_analysis_method_ID",
                fluidRow(
                  box(title = tags$b("Semi-continuous Analysis"), width = 12, height = "auto",
                    column(
                    selectInput(
                      "select_mixture_type_ID",
                      label = "mtDNA source type",
                      choices = c('Single source', '2-persons mixture', '3-persons mixture' 
                      )
                    ),
                    width = 6
                  ),
                  column(width = 6,
                         selectInput("select_files_SC_ID", "Choose input files",c(Choose=''),multiple = TRUE, selectize = TRUE)
                    
                  ),
                  withBusyIndicatorUI(actionButton("generate_mixture_statistics_ID", "Generate Mixture Statistics",class="btn-info"))
                  )
                ),
                fluidRow(
                  box(title = tags$b("RESULTS"), width = 12, height = "auto",
                      div(
                        tags$b("RMNE Statistics"),
                        rHandsontableOutput("mixStat_rmne_rhotOut")
                      ),
                      tags$br(),
                      div(
                        tags$b("LR Statistics"),
                        rHandsontableOutput("mixStat_lr_rhotOut")
                      )
                      
                      )
                )
                
          
        ),
        tabItem(tabName = "continuous_analysis_method_ID",
                fluidRow(
                  box(title = tags$b("Continuous Analysis"), width = 12, height = "auto",
                      fluidRow(#style = "display: block;",
                        column(width = 6,
                              fileInput("select_var_excel_quantition_ID",
                               label = "Upload excel file with quantitative data",
                               accept = ".xlsx"
                               ))),
                      fluidRow(
                        #style = "display: block;", 
                        column(width = 3,
                               textInput("cont_edit_dist_input_ID", "Enter edit distance to be used", value = 4)),
                        column(width = 3,
                               textInput("cont_numMCMC_input_ID","Enter the number of MCMC steps to be used", value = 800)),
                        column(width = 3, # TODO add litte info icon to mention about what happens if less than 15
                               textInput("cont_panel_size_input_ID", "Enter panel size", value = 25)),
                        column(width = 3,
                               numericInput("cont_numPerson_inMix_input_ID","Enter the number of persons in mixture",value = 2, min = 1,max = 5)),
                        ),
                      fluidRow(
                        column(width = 3,
                               textInput("cont_recombRate_input_ID","Recombination rate", value = 0.0)),
                        column(width = 3,
                               textInput("cont_miscopying_rate_input_ID","Miscopying rate", value = 0.01))
                      ),
                      withBusyIndicatorUI(actionButton("genereate_cont_run_deploid", label = "Run DEploid",class="btn-info")),
                      tags$style("#genereate_cont_run_deploid{float:left;width:auto}")
                      )
                ),
                fluidRow(
                  box(title = tags$b("RESULTS"), width = 12, height = "auto",
                      div(
                        tags$b("Proportion of Mixtures"),
                        rHandsontableOutput("cont_deploid_mixProp_rhotOut")
                      ),
                      tags$br(),
                      div(
                        tags$b("Estimated Haplotypes"),
                        rHandsontableOutput("cont_deploid_estimatedHap_rhotOut")
                        )

                      
                  )
                )
                
                )
        
    ),
  use_bs_tooltip(),
  tags$script("$(function() {
                                            $('#incl_excl_accor_ID').on('click', function(x) {
                                              Shiny.onInputChange('selected_tab_accordi_ID', x.target.innerText)
                                            });
                                          });"),
  HTML("<input id='selected_tab_accordi_ID' type='text' style='display: none;'>")
    
    
)



dashboardPage(
  dashboardHeader(title = "MMDIT analysis tool"),
  sidebar,
  body
  
)

