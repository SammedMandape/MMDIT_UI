#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(magrittr)
library(shinydashboard)
library(rhandsontable)
library(shinythemes)
#library(shinysky)
library(shinyjs)
library(MMDIT)
library(shinyWidgets)
library(DT)

source("Ufunctions_updated_08122020.R")


shinyServer(function(input, output, session) {
  
    values <- reactiveValues(
      x_ss_name = NULL,
      y_mix_name = NULL,
      x_ss_data = NULL,
      y_mix_data = NULL,
      foo6 = NULL,
      mydata_db = NULL,
      mydata_pops = NULL,
      mydata_amps = NULL
    )
    
    observeEvent(input$select_ss_data_ID,{
      values$x_ss_name<-input$select_ss_data_ID[['name']]
      values$x_ss_datapath <- input$select_ss_data_ID[['datapath']]
      
    })
    
    observeEvent(input$select_mix_data_ID,{
      values$y_mix_name <- input$select_mix_data_ID[['name']]
      values$y_mix_datapath <- input$select_mix_data_ID[['datapath']]
    })
    
    x_val_name <- reactive({
      if(is.null(values$x_ss_name)){
        return(NULL)
      }else{
        return(append("\tSingle source data:",values$x_ss_name))
      }
    })
    
    y_val_name <- reactive({
      if(is.null(values$y_mix_name)){
        return(NULL)
      }else{
        return(append("\n\tMixture data:",values$y_mix_name))
      }
    })
    
    output$files_selected <- renderText({
      paste(c("Uploaded file[s]\n",x_val_name(),y_val_name()),collapse = " ")
      
      
    })
    
    
    observeEvent(input$clear_ID,{
      values$x_ss_name <- NULL
      values$y_mix_name <- NULL
      reset('select_ss_data_ID')
      reset('select_mix_data_ID')
    }
    )
    
    observe({
      if(is.null(values$x_ss_name) || is.null(values$y_mix_name)){
        shinyjs::disable("analyze_data_id")
      }else{
        shinyjs::enable("analyze_data_id")
      }
    })

    
    observeEvent(input$analyze_data_id,{
    values$x_ss_data <- NULL
    for (i in 1:length(values$x_ss_datapath)){
      mydata_ss <- MMDIT::Empop2variant(values$x_ss_datapath[i]) %>% 
          pull(Variant) %>%
          MMDIT::Variant2snp() %$%
          UnfoldSNP(Pos, Allele, Type) %>%
          dplyr::mutate(FileID = values$x_ss_name[i],
                        Source = "Single")
      values$x_ss_data <- bind_rows(values$x_ss_data,mydata_ss)
    }
    
    values$y_mix_data <- MMDIT::Empop2variant(values$y_mix_datapath) %>% 
      pull(Variant) %>%
      MMDIT::Variant2snp() %$%
      UnfoldSNP(Pos, Allele, Type) %>%
      dplyr::mutate(FileID = values$y_mix_name, Source = "Mixture")
    })
    
    output$empop_variant_input_snp <- rhandsontable::renderRHandsontable({
      if(is.null(values$y_mix_data) & is.null(values$y_mix_data)){
        return(NULL)
      }
      rhot_x_ss_data <- values$x_ss_data %>% 
        dplyr::filter(Type == "Substitution") %>% 
        dplyr::mutate(Start = as.integer(Pos-1), Stop=Pos) %>%
        dplyr::select(FileID, Start, Stop, Allele, Type, Source) 
      
      rhot_y_mix_data <- values$y_mix_data %>%
        dplyr::filter(Type == "Substitution") %>% 
        dplyr::mutate(Start = as.integer(Pos-1), Stop=Pos) %>%
        dplyr::select(FileID, Start, Stop, Allele, Type, Source)
      
      rhot_2gether <- bind_rows(rhot_x_ss_data,rhot_y_mix_data)
      rhandsontable(rhot_2gether, height = 400, width = 700)
    })
    
    
    
    output$empop_variant_input_indel <- rhandsontable::renderRHandsontable({
      if(is.null(values$y_mix_data) & is.null(values$y_mix_data)){
        return(NULL)
      }
      rhot_x_ss_data_indel <- values$x_ss_data %>%
        filter(Type=="Insertion"|Type=="Deletion") %>%
        mutate(Start = ifelse(Type=="Insertion",
                              Pos,
                              as.integer(Pos - 1)),
               Stop=ifelse(Type=="Insertion",
                           Pos,
                           Pos)) %>%
        select(FileID, Start, Stop, Allele, Type, Source)
      
      rhot_y_mix_data_indel <- values$y_mix_data %>%
        filter(Type=="Insertion"|Type=="Deletion") %>%
        mutate(Start = ifelse(Type=="Insertion",
                              Pos,
                              as.integer(Pos - 1)),
               Stop=ifelse(Type=="Insertion",
                           Pos,
                           Pos)) %>%
        select(FileID, Start, Stop, Allele, Type, Source)
      
      rhot_2gether_indel <- bind_rows(rhot_x_ss_data_indel, rhot_y_mix_data_indel)
      
      rhandsontable(rhot_2gether_indel, height = 500, width = 700)
        
      
    })
    
    
    observeEvent(input$next_to_indel_analysis_ID, {
        updateTabItems(session, "primary_analysis_ID", selected = "Indel analysis") 
    })
    
    observeEvent(input$start_over_ID,{
      session$reload()
    })
    
    
    # observeEvent(input$getmitogen_ID,{
    #   output$tmp_get_mitoseq <- rhandsontable::renderRHandsontable({
    #     foo7 <- values$foo6
    #     browser()
    #     db<-loadMMDIT()
    #     getMitoGenomes(db, pop = c("AF","AM"),blk = foo7) -> get_mit_gen_tab
    #     rhandsontable(get_mit_gen_tab, height = 400, width = 500)
    #   })
    # })
    
    
    observeEvent(input$load_MMDIT_ID,{
      withBusyIndicatorServer("load_MMDIT_ID",{
        Sys.sleep(1)
        tryCatch({
        values[['mydata_db']] <- loadMMDIT()
        values$mydata_pops <- getPops(values$mydata_db)
        values$mydata_amps <- getAmpCoordinates(values$mydata_db)},
        warning = function(w){
          stop("Warning! Something went wrong, try again!")
        }, error = function(e){
          stop("MMDIT didn't load correctly, please check mmdit.sqlite3 database exists!")
        }
        )
        })
    })
    
    observe({
      # TODO is tab 'Select populations' is selected and is.null(db) then alert to load db first. if(is.null())
      if(!is.null(values$mydata_pops)){
        pops <- values$mydata_pops %>% pull()
      updatePickerInput(session, inputId = "myPicker_accor_ID", choices = pops)
      }
    })
    
    getselect <- function(x){
      if(x == "Select all amplicons"){
        return(list(selected = c(1:162)))
      }
      if(x == "Select even rows of amplicons" ){
        return(list(selected = seq(2,162,2)))
      }
      if(x == "Select odd rows of amplicons"){
        return(list(selected = seq(1,162,2)))
      }
      if(x == "Manually select amplicons (default)"){
        return('multiple')
      }
    }
    output$mydata_amps_sel_dtID <- DT::renderDT({
      mydata_amps_sel <- values$mydata_amps
      },selection = getselect(input$Id081)
      )
    output$mydata_amps_sel_textout_ID <- renderPrint({cat('Selected amplicons (row numbers):\n\n')
                              input$mydata_amps_sel_dtID_rows_selected})
    
    observe({
      if(input$Id072 == "Choose from Precision ID Kit" & is.null(values[['mydata_db']])){
        showModal(modalDialog(title = "Precision kitID amplicon table",
                              div(tags$b("Invalid! Please load MMDIT database", style = "color: red;"))
        ))
      }
      if(input$Id072 == "Choose from Precision ID Kit" & (!is.null(values$mydata_db))){
        showModal(modalDialog(
          title = "Precision kitID amplicon table (Please select amplicons to include)",
          DT::dataTableOutput("mydata_amps_sel_dtID"), 
          size = "l", 
          style = "height:500px; overflow-y: scroll;",
          footer = tagList(div(style="display:inline-block; float: left;",
                           pickerInput(
                             inputId = "Id081",
                             label = NULL, 
                             choices = c("Manually select amplicons (default)",
                                         "Select all amplicons", 
                                         "Select even rows of amplicons",
                                         "Select odd rows of amplicons"),
                             options = list(
                               style = "btn-primary")
                           )),
                           div(style="display:inline-block;",actionButton("cont_precisionID_ID","Continue")),
                           div(style="display:inline-block;margin:0px 0px 0px 7px;",modalButton("Dismiss"))
                           )
          
          ))
      }
    })
    
    
    observeEvent(input$cont_precisionID_ID,{
      if(length(input$mydata_amps_sel_dtID_rows_selected) == 0){
        showModal(modalDialog(
          div(tags$b(style="color: red;", "Please select amplicons to include!"))
        ))
      }else{
      showModal(modalDialog(
        title = "Input semicolon separated genomic coordinates to exclude in bed format",
        div(style = "max-width: 100%;",textInput("text_excl_cont_preIDkit_ID","", value = "") %>% 
                           bsplus::shinyInput_label_embed(
                             shiny::icon("info-circle") %>%
                               bs_embed_tooltip("0-based start and 1-based stop. For example, 298-300;320-325")
                           )),
        div(verbatimTextOutput("mydata_amps_sel_textout_ID")),
        footer = tagList(actionButton("done_cont_perc_kitID_id","Done"))
      ))
    }
    })
    
    observeEvent(input$done_cont_perc_kitID_id,{
      removeModal()
      browser()
    })
    
    
    
    # observeEvent(input$Id081,{
    #   if(input$Id081 == "Select even rows of amplicons"){
    #     showModal(modalDialog(
    #       title = "Precision kitID amplicon table (Please select the amplicons to include)",
    #       DT::dataTableOutput("mydata_amps_sel_dtID"), 
    #       size = "l", 
    #       style = "height:500px; overflow-y: scroll;"
    #     ))
    #   }
    # })
        
    
############################################################################################################################    
    # # There is a bug in rhandsontable in shinymodal: the table is loaded correctly the first time it is called
    # # but the second time it is called in shinymodal it truncates after first column and will only display all
    # # the columns if scrolled the table or clicked on the table. The following is a workaround to that bug,
    # # where Rhandsontable is forcefully rendered each time modal is called. This is done by indirectly passing
    # # the amplicons generated via the following reactive element.
    #  mydata_amps_tmp <- reactive({
    #    if(is.null(values[['mydata_amps_select']])){
    #      if(is.null(input$mydata_amps_rhot)){
    #          values$mydata_amps %>% mutate("Select amps" = rep(TRUE,1))
    #      } else{
    #        #values$mydata_amps %>% mutate("Select amps" = rep(FALSE,1))
    #        # TODO bug: when user continues and if everything is deselected -> user dismisses modal -> 
    #        # user goes back to that modal everything is deselected by default -> rather it should be
    #        # selected by default.
    #        hot_to_r(input$mydata_amps_rhot)
    #      }
    #    }else if(values[['mydata_amps_select']] == 2){
    #       return(values$mydata_amps %>% mutate("Select amps" = rep(FALSE,1)))
    #    }
    #  })
    #  
    # # instead of directly passing values$mydata_amps to renderRHandsontable it is indirectly passed, by bringing in 
    # # reactive element mydata_amps_tmp in between. 
    # observe({#values$mydata_amps,{
    #   if(is.null(values[['mydata_amps_select']])){
    #   output$mydata_amps_rhot <- rhandsontable::renderRHandsontable({
    #     rhandsontable(mydata_amps_tmp(), width = "100%", height = 200)
    #   })
    #   }else if(values[['mydata_amps_select']] == 2){
    #   output$mydata_amps_rhot_selected <- rhandsontable::renderRHandsontable({
    #     rhandsontable(mydata_amps_tmp(), width = "100%", height = 200)
    #   })
    #   }
    # })
    # 
    # # observe({
    # #   if(input$Id072 == "Choose from Precision ID Kit" & is.null(values[['mydata_db']])){
    # #     showModal(modalDialog(title = "Precision kitID amplicon table",
    # #       div(tags$b("Invalid! Please load MMDIT database", style = "color: red;"))
    # #     ))
    # #   }
    # #   if(input$Id072 == "Choose from Precision ID Kit" & (!is.null(values$mydata_db))){
    # #     values[['mydata_amps_select']] <- NULL
    # #     showModal(modalDialog(
    # #       title = "Precision kitID amplicon table (Please select the amplicons to include)",
    # #       rHandsontableOutput("mydata_amps_rhot"),
    # #       footer = tagList(
    # #         actionButton("deselct_all_ID", "Deselect all"),
    # #         actionButton("continue_prec_ID_kit_ID", "Continue"),tags$style("#deselct_all_ID{display:inline-block;float:left; color:red;}")
    # #         ), easyClose = FALSE
    # #     
    # #       # TODO: selectall button, change exit to escape or click anywhere
    # #       
    # #       
    # #     )
    # #     )
    # #   }
    # #   
    # #   if(input$Id072 == "Manually input genomic coordinate intervals" & is.null(values[['mydata_db']])){
    # #     showModal(modalDialog(
    # #       title = "Input genomic coordinate intervals in bed format",
    # #       div(tags$b("Invalid! Please load MMDIT database", style = "color: red;"))
    # #     ))
    # #   }else if(input$Id072 == "Manually input genomic coordinate intervals" & (!is.null(values$mydata_db))){
    # #     showModal(modalDialog(
    # #       title = "Input genomic coordinate intervals to include in bed format",
    # #       textInput("text_incl_manual_ID","", value = "") %>% bsplus::shinyInput_label_embed(
    # #         shiny::icon("info-circle") %>%
    # #           bs_embed_tooltip("0-based start and 1-based stop")
    # #       ),
    # #       footer = tagList(actionButton("continue_manual_input_ID","Continue"))
    # #     ))
    # #   }
    # #   
    # # })
    # 
    # observeEvent(input$deselct_all_ID,{
    #   values[['mydata_amps_select']] <- 2
    #   showModal(modalDialog(
    #     title = "Precision kitID amplicon table (Please select the amplicons to include)",
    #     rHandsontableOutput("mydata_amps_rhot_selected"),
    #     footer = tagList(
    #       actionButton("continue_prec_ID_kit_ID", "Continue")
    #                      )
    #   ))
    # })
    # 
    # 
    # observeEvent(input$continue_prec_ID_kit_ID,{
    #   showModal(modalDialog(
    #     title = "Input semicolon separated genomic coordinates to exclude in bed format",
    #     textInput("text_excl_cont_preIDkit_ID","", value = "") %>% 
    #                  bsplus::shinyInput_label_embed(
    #                    shiny::icon("info-circle") %>%
    #                      bs_embed_tooltip("0-based start and 1-based stop")
    #                  ),
    #     footer = tagList(actionButton("save_modal_ID", "Save"), tags$style("#save_modal_ID{display:inline-block;float:left; color:red;}"),
    #                      modalButton("Dismiss"))
    #   ))
    # })
    # 
    # observeEvent(input$continue_manual_input_ID,{
    #   showModal(modalDialog(
    #     title = "Input semicolon separated genomic coordinates to exclude in bed format",
    #     textInput("text_excl_cont_manual_ID","",value = "") %>% 
    #       bsplus::shinyInput_label_embed(
    #         shiny::icon("info-circle") %>%
    #           bs_embed_tooltip("0-based start and 1-based stop")
    #       )
    #   ))
    # })
############################################################################################################################    
    

})
