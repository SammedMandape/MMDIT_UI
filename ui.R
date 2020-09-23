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
                    box(title = "File import", width = 12, height = "auto",
                        column(width = 6,
                               fileInput("select_ss_data_ID",
                                         label = "Upload one or multiple single source data",
                                         multiple = TRUE,
                                         accept = "text/plain"),
                               fileInput("select_mix_data_ID",
                                         #"empop_variant_file_id",
                                         label = "Upload mixture data",
                                         accept = "text/plain"),
                               fluidRow(actionButton(
                                   "analyze_data_id", label = "Analyze data",class="btn-info"
                                ),
                               tags$style("#analyze_data_id{float:right}"),
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
                                                             max-height: 500px;
                                                             color:red;
                                                             white-space: pre-wrap;
                                                  }"))

                        )
                       )
                ),
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
                                    actionButton(inputId = "save_df_ID", label = "Save", class="btn-info"),#),
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
                                      actionButton("load_MMDIT_ID","Load MMDIT database", class='btn-info'),
                                        # style="margin: 10px; 
                                        # box-sizing: border-box;
                                        # "),
                                    tags$style("#load_MMDIT_ID{margin: 10px; 
                                        box-sizing: border-box;
                                    }"),
                                    bs_accordion(id = "incl_excl_accor_ID") %>%
                                      bs_set_opts(panel_type = "success", use_heading_link = TRUE) %>%
                                      bs_append(title = "Choose Incl Exc list", 
                                                content = #div(class = "radioselect_inclExcl",
                                                            fluidRow(#column(1,
                                                                     radioGroupButtons(
                                                                       inputId = "ID070", label = "", choices = c(
                                                                         "Inclusion" = "incl_choice_ID",
                                                                         "Exclusion" = "excl_choice_ID"
                                                                       ), justified = TRUE, checkIcon = list(yes = icon("ok", 
                                                                                                                        lib = "glyphicon"))
                                                                     ),
                                                              #),
                                                              #column(4,
                                                                     textInput("text_incl_ID","", value = "") %>% 
                                                                       bsplus::shinyInput_label_embed(
                                                                         #icon("info") %>%
                                                                         shiny::icon("info-circle") %>%
                                                                           bs_embed_tooltip("Input semicolon separated coordinates (0-based start and 1-based stop). 
                                                                                            Default is inclusion of whole mitochondrial genome.")
                                                                       )
                                                                     #textInput("text_Excl_ID","", value = "")
                                                                     # bs_button("I'm a button") %>%
                                                                     #   bs_embed_tooltip(title = "I'm a tooltip")
                                                              #)
                                                              )
                                                              #)
                                                  
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
                                                )),
                                    # div(class = "radioselect_inclExcl",
                                    # column(1,
                                    # radioButtons("incl_excl_ID", "Choose", 
                                    #              choices = c(
                                    #                          "Inclusion" = "inclusion_ID",
                                    #                          "Exclusion" = "exclusion_ID"),selected = NULL)
                                    # )),
                                    # column(2, textInput("text_incl_ID","", value = "1-16569") %>% 
                                    #          bsplus::shinyInput_label_embed(
                                    #            #icon("info") %>%
                                    #              shiny::icon("info-circle") %>%
                                    #              bs_embed_tooltip("Input semicolon separated coordinates (0-based start and 1-based stop)")
                                    #          ),
                                    #        textInput("text_Excl_ID","", value = "")
                                    #        # bs_button("I'm a button") %>%
                                    #        #   bs_embed_tooltip(title = "I'm a tooltip")
                                    #        ),
                                    # tags$style(".radioselect_inclExcl{
                                    #            font-weight: bold;
                                    #            line-height: 2.0;
                                    #            margin: 10px;
                                    #            }"),
                                    # column(3,
                                    #        bs_accordion(id = "select_population_acco_ID") %>%
                                    #          bs_set_opts(panel_type = "success", use_heading_link = TRUE) %>%
                                    #          bs_append(title= "Select populations",
                                    #                    content = pickerInput(
                                    #          inputId = "myPicker", 
                                    #          label = "Select/deselect all + format selected", 
                                    #          choices = LETTERS, 
                                    #          options = list(
                                    #            `actions-box` = TRUE, 
                                    #            size = 4,
                                    #            `selected-text-format` = "count > 4"
                                    #          ), 
                                    #          multiple = TRUE
                                    #        ))
                                    #        )

                                    ),
                                    )
                    
                    
                           #tabPanel("")

                    )
                )
                #)
                )
                ),
        tabItem(tabName = "semi-continuous_analysis_method_ID",
                fluidRow(
                  box(title = "Semi-continuous analysis", width = 12, height = "auto",
                    column(
                    selectInput(
                      "select_mixture_type_id",
                      label = "Mixture type",
                      choices = c('Single source', '2-person mixture (2-knowns)', '2-person mixture (2-unknowns)',
                                  '2-person mixture (1-known, 1-unknown)','3-person mixture (3-unknowns)',
                                  '3-person mixture (1-known, 2-unknown)','3-person mixture (2-known, 1-unknown)'
                      )
                    ),
                    width = 6
                  ),
                  column(width = 6,
                         selectInput("select_files_SC_ID", "Choose input files",c(Choose=''),multiple = TRUE, selectize = TRUE)
                    
                  )
                  )
                )
                
          
        )
        
    ),
  use_bs_tooltip()
    
    
)



dashboardPage(
  dashboardHeader(title = "MMDIT analysis tool"),
  sidebar,
  body
  
)

