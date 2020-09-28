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


source("Ufunctions_updated_08122020.R")


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    ###################################
    # values <- reactiveValues(
    #     upload_state = NULL
    #     #browser()
    #     #x = NULL,
    #     #y = NULL
    # )
    # 
    # observeEvent({input$select_ss_data_ID
    #     input$select_mix_data_ID},{
    #     values$upload_state <- 'uploaded'
    #         #browser()
    #         #values$x
    #         #values$y
    #     })
    # 
    # observeEvent(input$clear_ID,{
    #     #browser()
    #     values$upload_state <- 'reset'
    # })
    # 
    # file_input <- reactive({
    #     #browser()
    #     if(is.null(values$upload_state)){
    #         #browser()
    #         return(NULL)
    #     }else if (values$upload_state == 'uploaded'){
    #         as_tibble(input$select_ss_data_ID) -> tmp_ss
    #         as_tibble(input$select_mix_data_ID) -> tmp_m
    #         bind_rows(tmp_ss, tmp_m) -> tmp_ss_m
    #         #browser()
    #         return(tmp_ss_m)
    #     }else if(values$upload_state == 'reset'){
    #         #browser()
    #         return(NULL)
    #     }
    # })
    # 
    # output$files_selected <- renderText({
    #     return(paste(c("Uploaded file: ", file_input()$name), collapse = " "))
    # })
    ##########################################################
  
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
      # if(!is.null(x_val())){
      #   browser()
      # }
      paste(c("Uploaded file[s]\n",x_val_name(),y_val_name()),collapse = " ")
      
      
    })
    
    # observe({
    # 
    #     values$x<-input$select_ss_data_ID[['name']]
    #     values$y<-input$select_mix_data_ID[['name']]
    #     # updateSelectInput(session, "select2input1",
    #     #                   choices = x,
    #     #                   selected = x)
    # 
    #     output$files_selected <- renderText({
    #         paste(c("Uploaded file:",x,y),collapse = " ")
    #         })
    # 
    #     #updateSelectInput(session, "select_files_SC_ID", choices = c(x,y))
    #     })
    
    observeEvent(input$clear_ID,{
      values$x_ss_name <- NULL
      values$y_mix_name <- NULL
      reset('select_ss_data_ID')
      reset('select_mix_data_ID')
      #browser()
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
      # if(!is.null(values$x_ss_datapath))
      #   browser()
    #mydata_ss <- vector(mode = "list",length = length(values$x_ss_datapath))
    values$x_ss_data <- NULL
    for (i in 1:length(values$x_ss_datapath)){
      mydata_ss <- MMDIT::Empop2variant(values$x_ss_datapath[i]) %>% 
          pull(Variant) %>%
          MMDIT::Variant2snp() %$%
          UnfoldSNP(Pos, Allele, Type) %>%
          dplyr::mutate(FileID = values$x_ss_name[i],
                        Source = "Single")
      values$x_ss_data <- bind_rows(values$x_ss_data,mydata_ss)
      #browser()
    }
    
    values$y_mix_data <- MMDIT::Empop2variant(values$y_mix_datapath) %>% 
      pull(Variant) %>%
      MMDIT::Variant2snp() %$%
      UnfoldSNP(Pos, Allele, Type) %>%
      dplyr::mutate(FileID = values$y_mix_name, Source = "Mixture")
    # values$x_ss_data <- MMDIT::Empop2variant(values$x_ss_datapath) %>% 
    #   pull(Variant) %>% 
    #   MMDIT::Variant2snp() %$% 
    #   UnfoldSNP(Pos, Allele, Type) %>% 
    #   dplyr::mutate(FileID = values$x_ss_name) 
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
      browser()
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
      
      rhandsontable(rhot_2gether_indel, height = 500, width = 500)
        
      
    })
    #output$empop_variant_input_snp_ss
    # empop_variant_file_id <- eventReactive(c(input$empop_variant_file_id),
    #                                        {
    #                                            empop_variant_file_path <- input$empop_variant_file_id
    #                                        }
    # )
    # 
    # 
    # 
    # output$empop_variant_input_snp <- rhandsontable::renderRHandsontable({
    #     #browser()
    #     empop_variant_file <- empop_variant_file_id()$datapath
    #     #print(input$analyze_data_id)
    #     #browser()
    #     mydata_empop <- Empop2variant(empop_variant_file)
    #     mydata_empop_1 <- mydata_empop %>% mutate(FileID = empop_variant_file_id()$name)
    #     if(input$analyze_data_id == 0){
    #         return(rhandsontable(mydata_empop_1, width = 700, height = 400)
    #                %>% hot_table(highlightRow = TRUE, highlightCol = TRUE))
    #     }
    # 
    #     mydata_variant_unfolded_snp <- Variant2snp(mydata_empop$Variant) %$%
    #         UnfoldSNP(Pos, Allele, Type) %>%
    #         mutate(FileID = empop_variant_file_id()$name) %>%
    #         filter(Type == "Substitution") %>%
    #         mutate(Start = as.integer(Pos-1), Stop=Pos) %>%
    #         select(FileID, Start, Stop, Allele, Type)
    #         return(rhandsontable(mydata_variant_unfolded_snp, width = 700, height = 400) %>%
    #                hot_table(highlightRow = TRUE, highlightCol = TRUE))
    # })
    # 
    # 
    # output$empop_variant_input_indel <- rhandsontable::renderRHandsontable({
    #     empop_variant_file <- empop_variant_file_id()$datapath
    #     mydata_empop <- Empop2variant(empop_variant_file)
    #     mydata_variant_unfolded_indel <- Variant2snp(mydata_empop$Variant) %$%
    #         UnfoldSNP(Pos, Allele, Type) %>%
    #         mutate(FileID = empop_variant_file_id()$name) %>%
    #         filter(Type=="Insertion"|Type=="Deletion") %>%
    #         mutate(Start = ifelse(Type=="Insertion",
    #                               Pos,
    #                               as.integer(Pos - 1)),
    #                Stop=ifelse(Type=="Insertion",
    #                            Pos,
    #                            Pos)) %>%
    #         select(FileID, Start, Stop, Allele, Type)
    #     return(rhandsontable(mydata_variant_unfolded_indel, width = 700, height = 400) %>%
    #            hot_table(highlightRow = TRUE, highlightCol = TRUE))
    # 
    # })
    
    observeEvent(input$next_to_indel_analysis_ID, {
        updateTabItems(session, "primary_analysis_ID", selected = "Indel analysis") 
    })
    #observeEvent(input)
    
    observeEvent(input$start_over_ID,{
      session$reload()
    })
    
    # observe({
    #   
    #   if(input$ID070 == "excl_choice_ID"){
    #     foo <- input$text_incl_ID
    #     foo1 <- as_tibble(foo)
    #     foo2 <- foo1 %>% separate_rows(value,sep = ";")
    #     foo1 %>% separate_rows(value,sep = ";") %>% separate(value, into = c("Start","Stop"), convert = T) -> foo3
    #     foo3 %>% mutate(len = Stop - Start) -> foo4
    #     sequence(foo4$len) + rep(foo4$Start, foo4$len) -> foo5
    #     unique(foo5) -> values$foo6
    #     #browser()
    #   }
    #   
    #   # shinyjs::toggleState("text_incl_ID",input$incl_excl_ID == "inclusion_ID")
    #   # shinyjs::toggleState("text_Excl_ID",input$incl_excl_ID == "exclusion_ID")
    #   # if(input$incl_excl_ID == "exclusion_ID"){
    #   # updateTextInput(session, "text_incl_ID", value = "")
    #   # }else{
    #   #   updateTextInput(session, "text_incl_ID", value = "1-16569")
    #   # }
    # })
    
    observeEvent(input$getmitogen_ID,{
      output$tmp_get_mitoseq <- rhandsontable::renderRHandsontable({
        foo7 <- values$foo6
        browser()
        db<-loadMMDIT()
        getMitoGenomes(db, pop = c("AF","AM"),blk = foo7) -> get_mit_gen_tab
        rhandsontable(get_mit_gen_tab, height = 400, width = 500)
        
        #browser()
        
      })
    })
    
    # observe({
    #   if(is.null(input$selected_tab_accordi_ID)){
    #     return(NULL)
    #   }else if(!is.null(input$selected_tab_accordi_ID) & (input$selected_tab_accordi_ID == "Select populations")){
    #     db <- loadMMDIT()
    #     if(is.null(db)){
    #       stop("Database is not loaded")
    #     }
    #     values$mydata_pops <- MMDIT::getPops(db) %>% pull()
    #     updatePickerInput(session, inputId = "myPicker_accor_ID", choices = values$mydata_pops)
    #   } 
    # })
    
    # observe({
    #   if(input$primary_analysis_ID == "Inclusion/Exclusion"){
    #     db <- loadMMDIT()
    #     values$mydata_pops <- MMDIT::getPops(db) %>% pull()
    #     updatePickerInput(session, inputId = "myPicker_accor_ID", choices = values$mydata_pops)
    #   }
    # })
    
    observeEvent(input$load_MMDIT_ID,{
      values[['mydata_db']] <- loadMMDIT()
      values$mydata_pops <- getPops(values$mydata_db)
      values$mydata_amps <- getAmpCoordinates(values$mydata_db)
      values$mydata_amps <- values$mydata_amps %>% mutate("Select amps" = rep(FALSE,1))
    })
    
    observe({
      # TODO is tab 'Select populations' is selected and is.null(db) then alert to load db first. if(is.null())
      if(!is.null(values$mydata_pops)){
        pops <- values$mydata_pops %>% pull()
      updatePickerInput(session, inputId = "myPicker_accor_ID", choices = pops)
      }
    })
    
    # There is a bug in rhandsontable in shinymodal: the table is loaded correctly the first time it is called
    # but the second time it is called in shinymodal it truncates after first column and will only display all
    # the columns if scrolled on the table or clicked on the table. The following is a workaround to that bug,
    # where Rhandsontable is forcefully rendered each time modal is called. This is done by indirectly passing
    # the amplicons generated via the following reactive element.
     mydata_amps_tmp <- reactive({
       if(is.null(input$mydata_amps_rhot)){
         values$mydata_amps
       } else{
         hot_to_r(input$mydata_amps_rhot)
       }
     })
     
     
    observeEvent(values$mydata_amps,{
      output$mydata_amps_rhot <- rhandsontable::renderRHandsontable({
        rhandsontable(mydata_amps_tmp(), width = "100%", height = 200)
      })
    })
    observe({
      if(input$Id072 == "Choose from Precision ID Kit" & is.null(values[['mydata_db']])){
        showModal(modalDialog("Precision kitID amplicon table",
          div(tags$b("Invalid! Please load MMDIT database", style = "color: red;"))
        ))
      }
      if(input$Id072 == "Choose from Precision ID Kit" & (!is.null(values$mydata_db))){
        #browser()
        showModal(modalDialog(
          title = "Precision kitID amplicon table",
          #"This is an important message!",
          rHandsontableOutput("mydata_amps_rhot"),
          footer = tagList(actionButton("go2", "Continue"),
                           modalButton("Stop")), easyClose = TRUE
          # TODO: selectall button, change exit to escape or click anywhere
          
          
        ))
      }
    })
    #})

})
