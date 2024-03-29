# @Author: Sammed N Mandape
# Title: Bioinformatician
# Email: sammed.mandape@unthsc.edu
# Center for Human Identification
# University of North Texas Health Science Center, Fort Worth, TX, USA
#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(MMDIT)

# Need to set options for publishing to shinyapps.io
#options(repos = BiocManager::repositories()) 

# MMDIT needs to be updaed to include continuous' updated functions instead
# of sourcing part of if from "AllCodes5.R"
source("AllCodes5.R")


shinyServer(function(input, output, session) {
  
    values <- reactiveValues(
      x_ss_name = NULL,
      y_mix_name = NULL,
      x_ss_data = NULL,
      y_mix_data = NULL,
      mydata_db = NULL,
      mydata_pops = NULL,
      mydata_amps = NULL
      
    )
    
    # Homepage
    ####### NEW ADDITIONS #####################################################
    
    showModal(modalDialog(
      title = "Somewhat important message!",
      "This app is for research use only. We collect your usage data within the MMDIT app to perform analytics of usage and improve our app. 
      For details of policy please check ",tags$a("Google Privacy & Terms.",href="https://policies.google.com/privacy"),
      easyClose = FALSE,
      footer = tagList(div(style="display:inline-block;margin:0px 0px 0px 7px;",modalButton("Accept"))),
    ))
    
    observeEvent(input$continuous_btn_ID, {
      updateNavlistPanel(session, "tabs",
                         selected = "continuous_analysis_method_ID")
    })
    
    observeEvent(input$semi_continuous_btn_ID, {
      updateNavlistPanel(session, "tabs",
                         selected = "semi-continuous_analysis_method_ID")
    })
    
    observeEvent(input$lets_begin_ID, {
      updateNavlistPanel(session, "tabs",
                         selected = "data_import_ID")
    })
    
    observeEvent(input$user_guide_ID, {
      updateNavlistPanel(session, "tabs",
                         selected = "user_guide_tab_ID")
    })
    
    output$samplepdf <- renderUI({tags$iframe(`title`="UserGuide",style = "height: 800px; width: 100%; scrolling=yes", src = "UserGuide.pdf", lang="en")})
    
    
    #############################################################################
    
    
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
        return() #returns NULL
      }else{
        return(append("\tSingle source data:",values$x_ss_name))
      }
    })
    
    y_val_name <- reactive({
      if(is.null(values$y_mix_name)){
        return() # returns NULL
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
      values$x_ss_datapath <- NULL
      reset('select_ss_data_ID')
      reset('select_mix_data_ID')
    }
    )
    
    observe({
      if(is.null(values$y_mix_name)){
        shinyjs::disable("analyze_data_ID")
      }else{
        shinyjs::enable("analyze_data_ID")
      }
    })
    
    # Checking for heteroplasmy
    all_bases <- c("","A","G","T","C","a","g","t","c")
    
    observeEvent(input$analyze_data_ID,{
      values$x_ss_data <- NULL
      if(!is.null(values$x_ss_datapath)){
        ###########################################################################################
        # Single source samples
        for (i in 1:length(values$x_ss_datapath)){
          
          mydata_variant2snp <- empop2variant(values$x_ss_datapath[i]) %>% 
            pull(Variant) %>%
            Variant2snp() %>%
            mutate(Pos = as.integer(Pos))
          
          
          sample_alleles <- mydata_variant2snp %>% pull(Allele)
          
          if(!all(sample_alleles %in% all_bases)){
            showModal(modalDialog(title = "Please check!",
                                  "There is a/ are heteroplasmic site[s] in known /
                                  proposed contributors' data and will be removed 
                                  or you can change the input and reload the files." ))
            mydata_variant2snp <- mydata_variant2snp[(sample_alleles %in% all_bases),]
          }
          
          mydata_ss <- mydata_variant2snp %$%
              UnfoldSNP(Pos, Allele, Type) %>%
              dplyr::mutate(FileID = values$x_ss_name[i],
                            Source = "Single")
          values$x_ss_data <- bind_rows(values$x_ss_data,mydata_ss)
        }
        }
      ###########################################################################################
      
      
      ###########################################################################################
      # Mixture sample
      mydata_mix_variant2snp <- empop2variant(values$y_mix_datapath) %>% 
          pull(Variant) %>%
          variant2snp() %>% 
          mutate(Pos = as.integer(Pos)) 
      
      values$y_mix_data <- mydata_mix_variant2snp %$%
          UnfoldSNP(Pos, Allele, Type) %>%
          dplyr::mutate(FileID = values$y_mix_name, Source = "Mixture")
      ###########################################################################################
    })
    
    output$empop_variant_input_snp <- rhandsontable::renderRHandsontable({
      if(is.null(values$y_mix_data) && is.null(values$y_mix_data)){
        return(NULL)
      }
      rhot_x_ss_data <- NULL
      if(!is.null(values$x_ss_data)){
        rhot_x_ss_data <- values$x_ss_data %>% 
          dplyr::filter(Type == "Substitution") %>% 
          dplyr::mutate(Start = as.integer(Pos-1), Stop=Pos) %>%
          dplyr::select(FileID, Start, Stop, Allele, Type, Source) 
      }
      
      values[['rhot_y_mix_data']] <- values$y_mix_data %>%
        dplyr::filter(Type == "Substitution") %>% 
        dplyr::mutate(Start = as.integer(Pos-1), Stop=Pos) %>%
        dplyr::select(FileID, Start, Stop, Allele, Type, Source)
      
      rhot_2gether <- bind_rows(rhot_x_ss_data,values[['rhot_y_mix_data']])
      rhandsontable(rhot_2gether, height = 400, width = 700)
    })
    
    
    
    output$empop_variant_input_indel <- rhandsontable::renderRHandsontable({
      if(is.null(values$y_mix_data) && is.null(values$y_mix_data)){
        return(NULL)
      }
      rhot_x_ss_data_indel <- NULL
      if(!is.null(values$x_ss_data)){
        rhot_x_ss_data_indel <- values$x_ss_data %>%
          filter(Type=="Insertion"|Type=="Deletion") %>%
          mutate(Start = ifelse(Type=="Insertion",
                                Pos,
                                as.integer(Pos - 1)),
                 Stop=ifelse(Type=="Insertion",
                             Pos,
                             Pos)) %>%
          select(FileID, Start, Stop, Allele, Type, Source)
      }
      
      values[['rhot_y_mix_data_indel']] <- values$y_mix_data %>%
        filter(Type=="Insertion"|Type=="Deletion") %>%
        mutate(Start = ifelse(Type=="Insertion",
                              Pos,
                              as.integer(Pos - 1)),
               Stop=ifelse(Type=="Insertion",
                           Pos,
                           Pos)) %>%
        select(FileID, Start, Stop, Allele, Type, Source)
      
      rhot_2gether_indel <- bind_rows(rhot_x_ss_data_indel, values[['rhot_y_mix_data_indel']])
      
      # Looks like rhandsontable has bug and it won't take source if the length of vector is 1,
      # following is a work around to replicate single element
      uniqFileID <- distinct(rhot_2gether_indel, FileID) %>% pull()
      if(length(uniqFileID)==1) uniqFileID <- rep(uniqFileID,2)
      
      rhandsontable(rhot_2gether_indel, height = 500, width = 700) %>%
        rhandsontable::hot_col(col = "FileID", type = "dropdown", source = uniqFileID, 
                               strict = TRUE)
      
    })
    
    
    # Save user edited rhandsontable for indel mixture data and get
    # the mixture calls by converting hot to r object and binding
    observeEvent(input$save_df_ID,{
      withBusyIndicatorServer("save_df_ID",{
        Sys.sleep(1)
        values[['my_mix_data']] <- as_tibble(bind_rows(hot_to_r(input$empop_variant_input_indel) %>%
                                                         filter(Source == "Mixture"),values[['rhot_y_mix_data']]))
      })
      
    })
    
    
    observeEvent(input$next_to_indel_analysis_ID, {
      updateTabItems(session, "primary_analysis_ID", selected = "Indel analysis") 
    })
    
    observeEvent(input$start_over_ID,{
      session$reload()
    })
    
    
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
      if(!is.null(values$mydata_pops)){
        pops <- values$mydata_pops %>% pull()
        
        if("AF" %in% pops) pops[pops %in% c("AF")] <- "African"
        if("AM" %in% pops) pops[pops %in% c("AM")] <- "American"
        if("AS" %in% pops) pops[pops %in% c("AS")] <- "Asian"
        if("EU" %in% pops) pops[pops %in% c("EU")] <- "European"
        if("OC" %in% pops) pops[pops %in% c("OC")] <- "Oceanian"
      updatePickerInput(session, inputId = "myPicker_accor_ID", choices = pops, selected = pops)
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
      if(input$Id072 == "Choose from Precision ID Kit" && is.null(values[['mydata_db']])){
        showModal(modalDialog(title = "Error!",
                              div(tags$b("Invalid! Please load MMDIT database", style = "color: red;"))
        ))
      }
      if(input$Id072 == "Choose from Precision ID Kit" && (!is.null(values[['mydata_db']]))){
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
    
    observe({
      if(input$Id072 == "Upload bed file of regions to include / exclude" && is.null(values[['mydata_db']])){
        showModal(modalDialog(title = "Error!",
                              div(tags$b("Invalid! Please load MMDIT database", style = "color: red;"))
        ))
      }
      if(input$Id072 == "Upload bed file of regions to include / exclude" && (!is.null(values[['mydata_db']]))){
        showModal(modalDialog(
          title = "Please upload a bed file of regions to include",
          fileInput("bed_file_incl_ID",
                    label = "Upload bed file"),
          footer = tagList(
                    div(style="display:inline-block;",actionButton("cont_bed_file_incl_ID","Continue")),
                    div(style="display:inline-block;margin:0px 0px 0px 7px;",modalButton("Dismiss"))
            
          )
        ))
        
      }
      
      
      
      
    })
    
    
    observeEvent(input$cont_precisionID_ID,{
      if(length(input$mydata_amps_sel_dtID_rows_selected) == 0){
        showModal(modalDialog(title = "Error!",
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
        footer = modalButton("Done")
      ))
    }
    })
    
    
    observeEvent(input$done_run_backend_MMDIT_ID,{
      withBusyIndicatorServer("done_run_backend_MMDIT_ID",{
        
      # excluded regions
      if(input$text_excl_cont_preIDkit_ID==""){
        values[['data_excl_input_pkid_final']] <- 0
        }else {
      as_tibble(input$text_excl_cont_preIDkit_ID) -> mydata_excl_input_pkid
      mydata_excl_input_pkid %>%
        separate_rows(value,sep = ";") %>% 
        separate(value, into = c("Start","Stop"), convert = T) %>%
        mutate(len = Stop - Start) -> mydata_excl_input_pkid_1
      sequence(mydata_excl_input_pkid_1[['len']]) +
        rep(mydata_excl_input_pkid_1[['Start']],
            mydata_excl_input_pkid_1[['len']]) -> mydata_excl_input_pkid_2
      unique(mydata_excl_input_pkid_2) -> values[['data_excl_input_pkid_final']]
      }
      
     
      # included regions
      values$mydata_amps[input$mydata_amps_sel_dtID_rows_selected,] %>%
        unite("temp",start,stop, remove = FALSE, sep = "-") -> mydata_incl_input_pkid_tmp
      
      if("16541-16649" %in% mydata_incl_input_pkid_tmp$temp){
        mydata_incl_input_pkid_tmp %>% # to take care of circular DNA, maybe there 
          filter(start!="16541" & stop != "16649") %>% # is a better way to do this?
          add_row(start=c(16541,0),stop=c(16569, 80)) -> mydata_incl_input_pkid_tmp 
        }
      
      mydata_incl_input_pkid_tmp %>% 
        select(-temp) %>%
        mutate(len = stop - start) -> mydata_incl_input_pkid
      sequence(mydata_incl_input_pkid[['len']]) +
        rep(mydata_incl_input_pkid[['start']], 
            mydata_incl_input_pkid[['len']]) -> mydata_incl_input_pkid_1
      unique(mydata_incl_input_pkid_1) %>% 
        sort() -> values[['data_incl_input_pkid_final']]
       
      
      
      # getting all exclusions
      mtDNALen <- getMtgenomeLength(values[['mydata_db']])
      mtDNA <- seq(mtDNALen)
      intersect(values[['data_excl_input_pkid_final']], values[['data_incl_input_pkid_final']]) -> mydata_excl_intersect_pkid_1
      setdiff(mtDNA, values[['data_incl_input_pkid_final']]) -> mydata_excl_intersect_pkid_2
      union(mydata_excl_intersect_pkid_1, mydata_excl_intersect_pkid_2) -> values[['data_excl_total_pkid']]
      
      # taking all the arguments - inclusion / exclusion list, and population 
      # generating Mitogenomes strings from database used as proxy for 
      # random mitochondrial sequences.
      if(length(input$myPicker_accor_ID) >= 2){
        values[['population_selected']] <- input$myPicker_accor_ID
        if("African" %in% values[['population_selected']]) values[['population_selected']][values[['population_selected']] %in% c("African")] <- "AF"
        if("American" %in% values[['population_selected']]) values[['population_selected']][values[['population_selected']] %in% c("American")] <- "AM"
        if("Asian" %in% values[['population_selected']]) values[['population_selected']][values[['population_selected']] %in% c("Asian")] <- "AS"
        if("European" %in% values[['population_selected']]) values[['population_selected']][values[['population_selected']] %in% c("European")] <- "EU"
        if("Oceanian" %in% values[['population_selected']]) values[['population_selected']][values[['population_selected']] %in% c("Oceanian")] <- "OC"
        
      }else{
        showModal(modalDialog(
          title = "Error!",
          div(tags$b(style="color: red;", "Please choose at least 2 populations!"))
          
          
        ))
      }

      # get the mito reference genome that will be used in both methods
      values[['rcrs']] <- MMDIT::getMtgenomeSequence(values[['mydata_db']], double=FALSE)
      
      print("Done getting mitochondrial genomes from the database")
      
      # converting the final list of exclusion to tibble that is used to remove sites
      # from input data
      mydata_excl_final_tib <- as_tibble_col(values[['data_excl_total_pkid']], column_name = "Pos")
      
      # filter excluded sites from knowns (single source) data by anti-join
      if(!is.null(values$x_ss_data)){
        
        values[['ss_W_sites_excl']] <- anti_join(values$x_ss_data, 
                                                mydata_excl_final_tib,
                                                by = "Pos")
      }
      
      # exclude sites from mixture data (snps OR substitution)
      values[['my_mix_data_snps']] <- values[['my_mix_data']] %>% 
        filter(Type=="Substitution") %>% 
        mutate(Pos = Stop) %>% 
        anti_join(mydata_excl_final_tib, by="Pos") %>% select(-Pos)
      
      # **exclude sites from mixture data (insertion). This is only for insertions
      # defined for at a single positions, if there are insertions with a range
      # in positions this will not work**
      values[['my_mix_data_ins']] <- values[['my_mix_data']] %>% filter(Type=="Insertion") %>% mutate(Pos=Stop) %>%
        anti_join(mydata_excl_final_tib, by="Pos") %>% select(-Pos)
      
      # exclude sites from mixture data (deletions).
      my_mix_data_del_temp <- values[['my_mix_data']] %>% filter(Type=="Deletion") %>% mutate(len=Stop-Start)
      
      # if loop is to check if there is a range defined for deletions
      if(nrow(my_mix_data_del_temp %>% filter(len > 1)) > 0){
        # range means more than one pos for the event (eg: deletion in range 524-528)
        # filter to only include len > 1 and check if any pos in the deletion range is in
        # the excluded sites, if yes then the modal shows up, 
        # else the anti-join is done after removing the range of pos and 
        # adding back the events' range of pos
        my_mix_data_del_temp_1 <- my_mix_data_del_temp %>% filter(len > 1) 
        sequence(my_mix_data_del_temp_1[['len']]) + rep(my_mix_data_del_temp_1[['Start']],my_mix_data_del_temp_1[['len']]) -> my_mix_data_del_temp_2
        if(nrow(mydata_excl_final_tib[mydata_excl_final_tib$Pos %in% my_mix_data_del_temp_2,])){
          showModal(modalDialog(title = "Warning!",
                                div(tags$b(style="color: red;", "The excluded positions are in the range defined under 'Indel analysis' tab. 
                              Please either modify the deletion range or modify the inclusion / exclusion list.
                                         "))
          ))
        }else{
          values[['my_mix_data_del']] <- my_mix_data_del_temp %>% filter(Type=="Deletion") %>% filter(len==1) %>% mutate(Pos=Stop) %>%
            anti_join(mydata_excl_final_tib, by="Pos") %>% bind_rows(my_mix_data_del_temp %>% filter(Type=="Deletion") %>% filter(len!=1)) %>%
            select(-c(len, Pos))
        }
        
      }else{
        values[['my_mix_data_del']] <- my_mix_data_del_temp %>% mutate(Pos=Stop) %>%
          anti_join(mydata_excl_final_tib, by="Pos") %>% select(-c(len, Pos))
        
      }
      
      # combine subs, ins, del of mixture data after excluding sites
      values[['my_mix_data_final']] <- bind_rows(values[['my_mix_data_snps']], 
                                                 values[['my_mix_data_ins']], 
                                                 values[['my_mix_data_del']] ) %>%
                                      mutate(Allele = if_else(Allele == '-', ' ', Allele))
      
      })
      
    })
    
    #######################################################################
    ## Estimation of Match Statistics
    #######################################################################
    
    # function of generating string to be applied / map to each group of knowns
    my_seqdiffs2seq <- function(x_data, y_key){
      my_tmp_sd2s <- seqdiffs2seq(values[['rcrs']], x_data$Pos, x_data$evenType, x_data$Allele) %>%
        as_tibble()
      return(my_tmp_sd2s)
      
    }
    
    
    observeEvent(input$generate_mixture_statistics_ID,{
      withBusyIndicatorServer("generate_mixture_statistics_ID",{
        
        genomes_semiCont <- MMDIT::getMitoGenomes(values[['mydata_db']], 
                                                pop = values[['population_selected']], 
                                                blk = values[['data_excl_total_pkid']],
                                           ignoreIndels = FALSE)
        genomes_semiCont$sampleid <- as.character(genomes_semiCont$sampleid)
      
        # Preprocess genome into distinct haplotypes and counts 
        mygenomes <- MMDIT::preprocessMitoGenomes(genomes_semiCont)
        genomes <- mygenomes[[1]]
        genCount <- mygenomes[[2]]
        values[['foo_genomes_semiCont']] <- genomes
        
        
        #####################################################################################################
        # adding integer event types
        values[['my_ss_2_hap_final']] <- c()
        if(!is.null(values$x_ss_data)){
          
          my_ss_2_hap_W_sites_excl <- values[['ss_W_sites_excl']] %>% 
            group_by(FileID) %>% 
            mutate(evenType = as.integer(case_when(Type == "Substitution" ~ 0, 
                                                   Type == "Deletion" ~ 1, 
                                                   Type == "Insertion" ~ 2))) 
          
          
          # calculate seqdiffs2seq for each group, if more than one known present
          values[['my_ss_2_hap_final']] <- group_map(my_ss_2_hap_W_sites_excl, my_seqdiffs2seq) %>%
            unlist()
        }
        
        
        nInMix_val <- case_when(input$select_mixture_type_ID == "2-persons mixture" ~ 2,
                  input$select_mixture_type_ID == "3-persons mixture" ~ 3)
        
        # semi-continuous wrapper
        inter <- semicontinuousWrapper(genomes, 
                                       genCount, 
                                       values[['rcrs']],
                                       values[['my_mix_data_final']]$Start,
                                       values[['my_mix_data_final']]$Stop, 
                                       values[['my_mix_data_final']]$Allele,
                                       knownHaps = values[['my_ss_2_hap_final']],
                                       nInMix = nInMix_val)
        rmneStats <- inter[[1]]
        lrStats <- inter[[2]] %>% arrange(pop)
        
        
          
        output$mixStat_rmne_rhotOut <- rhandsontable::renderRHandsontable({
          rhandsontable(rmneStats)
          
        })
        
        
          
        output$mixStat_lr_rhotOut <- rhandsontable::renderRHandsontable({
          rhandsontable(lrStats)
          
        })
      })
      
      })
    
    #######################################################################
    ## Mixture Deconvolution
    #######################################################################
    observe({
      if(is.null(input$select_var_excel_quantition_ID)) shinyjs::disable("genereate_cont_run_deploid") else shinyjs::enable("genereate_cont_run_deploid")
    })
    
    
    observeEvent(input$select_var_excel_quantition_ID,{
      values[['altref_file']] <- input$select_var_excel_quantition_ID[['datapath']]
    })
    
    
    # Trace plot function showing convergence of MCMC steps
    plot.llk <- function (dEploid.run, title = "", cex.lab = 1, cex.main = 1, cex.axis = 1 ){
      llk = dEploid.run$llks
      llkEvent = dEploid.run$llksStates
      llk_sd = sd(llk)
      llk_range = range(llk)
      #windowsFonts(A = windowsFont("Times New Roman"))
      plot(llk, lty=2, type="l", col="black", xlab="Iteration", ylab="Log-likelihood", main=title,
          cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis)
      # plot(llk, lty=2, type="l", col="black", xlab="Iteration", ylab="Log-likelihood", main=title,
      #       cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, family="A", ylim=c(-2725,-2695), cex.lab=1.2) 
      updateSingleAt = which(llkEvent == 1)
      updateBothAt = which(llkEvent == 2)
      updatePropAt = which(llkEvent == 0)
      index = c(1:length(llk))
      # points(index[updateSingleAt], llk[updateSingleAt], cex = 0.6, col="red")
      # points(index[updateBothAt], llk[updateBothAt], cex = 0.6, col="blue")
      # points(index[updatePropAt], llk[updatePropAt], cex = 0.6, col="green")
    }
    
    
    observeEvent(input$genereate_cont_run_deploid,{
      # getting user inputs
      withBusyIndicatorServer("genereate_cont_run_deploid",{
      user_input_ED = as.integer(input$cont_edit_dist_input_ID)
      user_input_numMCMC = as.integer(input$cont_numMCMC_input_ID)
      user_input_panel_size = as.integer(input$cont_panel_size_input_ID)
      user_input_nInMix = input$cont_numPerson_inMix_input_ID
      user_input_readCntNorm = as.integer(input$cont_readCount_norm_ID)
      user_input_recombRate = as.numeric(input$cont_recombRate_input_ID)
      user_input_miscopyingRate = as.numeric(input$cont_miscopying_rate_input_ID)
      
      
      # check if panel size is atleast 15
      if(user_input_panel_size <15){
        showModal(modalDialog(title = "Error!",
                              div(tags$b(style="color: red;", "Panel size should be atleast 15. Setting it to less than 
                              15 might give erroneous results!
                                         "))))
      }
      
      cont_mix_data <- values[['my_mix_data_final']] %>% select(-c(FileID,Source)) %>% filter(Type == "Substitution") %>%
        mutate(Pos=Stop)
      altref <- Empop2AltRef(cont_mix_data, values[['altref_file']], normVal = user_input_readCntNorm)
      if(values$y_mix_name == "TheRealEMPOP_07908-034.txt" || values$y_mix_name == "TheRealEMPOP_025-07908.txt"){
        altref<-filter(altref, !is.na(NormalizedCount)) ##### not needed..if using JK's realempop
      }
      
      mtGenomes_cont <- MMDIT::getMitoGenomes(values[['mydata_db']], pop = values[['population_selected']], 
                                              blk = values[['data_excl_total_pkid']],
                                              ignoreIndels = TRUE)
      allDiffs <- getSeqdiffs(values[['mydata_db']], pop = values[['population_selected']], ignoreIndels = TRUE)
      nearN <- getNeNe(mPos = altref$Pos, mAllele = altref$Allele, 
                    Count = altref$NormalizedCount, IsRef = altref$isRef, 
                    ref = values[['rcrs']], mtGenomes = mtGenomes_cont, 
                    db = values[['mydata_db']], ED = user_input_ED, l = values[['data_excl_total_pkid']])
      nn<-nearN
      nn<-nn[1:user_input_panel_size] # reduce panel size while testing
      if(length(nn)<15)  {
        stop("This might give erroneous results, try increasing the panel size!")
      } 
      allDiffs[allDiffs$sampleid%in%nn,]%>%dplyr::select(-event)->longdf
      
      #######################################################################
      # TODO - Give user ability to choose seed
      #set.seed(5) # all odd amplicons for 07908_034
      #set.seed(3)
      #######################################################################
      
      dEploid.run<-rundEploid(l=values[['data_excl_total_pkid']], mPos = altref$Pos, mAllele = altref$Allele, Count = altref$NormalizedCount, IsRef = altref$isRef, 
                              SampleID = longdf$sampleid, rPos = longdf$position, rAllele =longdf$basecall, NumMCMC=user_input_numMCMC, exportPostProb =TRUE, 
                              recomb= 0.0, k= user_input_nInMix )
      
      
      # Trace plot
      output$cont_tace_plot_ID <- renderPlot(
        plot.llk(dEploid.run)
      )
      
      #run function to get mixture proportions from runDeploid object
      MixProp<-getMixProps(dEploid.run) %>%
        as_tibble_col(column_name = "MixProp")
      
      output$cont_deploid_mixProp_rhotOut <- rhandsontable::renderRHandsontable({
        rhandsontable(MixProp)
      })
      
      #run function to get estimated haplotypes from runDeploid object
      Haps<-getdEploidHaps(dEploid.run, mPos = altref$Pos, 
                           mAllele = altref$Allele, 
                           Count = altref$NormalizedCount, 
                           IsRef = altref$isRef) %>%
        mutate(Pos=as.integer(Pos))
      print("donegetdeploid")
      output$cont_deploid_estimatedHap_rhotOut <- rhandsontable::renderRHandsontable({
        rhandsontable(Haps, height = 500, width = 300)
      })
      print(user_input_numMCMC)
      print(user_input_panel_size)
      
      
      ########################################################################
      ##get hamming distance - only for internal use
      ########################################################################
      # v1<-Empop2variant("047_empop.txt")
      # v2<-Empop2variant("025_empop.txt")
      # s1<-Variant2snp(v1$Variant)
      # s2<-Variant2snp(v2$Variant)
      # u1<-UnfoldSNP(s1$Pos, s1$Allele,s1$Type)%>%mutate(SampleID=1)%>%filter(Type!="Insertion" & Type!="Deletion")
      # u2<-UnfoldSNP(s2$Pos, s2$Allele,s2$Type)%>%mutate(SampleID=2)%>%filter(Type!="Insertion" & Type!="Deletion")
      # Truehaps<-bind_rows(u1, u2)%>%dplyr::select(SampleID , Pos)#
      # dplyr::filter(Truehaps, !Truehaps$Pos%in%dEploid.run$ExSites)->Truehaps1#remove sites that violate the infinite sites model
      # tSampleID<-Truehaps1$SampleID
      # tPos<-Truehaps1$Pos
      # HammingDist <- getHapDis(dSampleID = Haps$Samples,  dPos = Haps$Pos, tSampleID, tPos)
      
      
      ##################
      # RMP calculation
      ##################
      # update sites to be excluded
      deploid_exSites <- dEploid.run[['ExSites']]
      updated_total_exSites <- unique(c(deploid_exSites, values[['data_excl_total_pkid']]))
      updated_total_exSites_tib <- as_tibble_col(updated_total_exSites, column_name = "Pos")
      
      # exclude all the sites from estimated haplotypes from deploid
      estHap_W_sites_excl <- anti_join(Haps, updated_total_exSites_tib, by="Pos")
      
      # generate vector of strings from haplotypes
      estHap_W_sites_excl_final <- estHap_W_sites_excl %>% 
        group_by(Samples) %>% 
        mutate(evenType = 0, Allele = Nuc) %>%
        select(-c(Nuc))

      estHap_final <- group_map(estHap_W_sites_excl_final, my_seqdiffs2seq) %>%
        unlist()
    
      # recalculate getMitoGenomes using the updated exclusion sites
      genomes_cont_4EstHap <- MMDIT::getMitoGenomes(values[['mydata_db']],
                                             pop = values[['population_selected']],
                                             blk = updated_total_exSites, 
                                             ignoreIndels = FALSE
                                             )
      
      genomes_cont_4EstHap$sampleid <- as.character(genomes_cont_4EstHap$sampleid)
      mygenomes_cont_4EstHap_1 <- MMDIT::preprocessMitoGenomes(genomes_cont_4EstHap)
      genomes_cont <- mygenomes_cont_4EstHap_1[[1]]
      
      # rmp
      estHap_rmp_database <- MMDIT::rmpWrapper(estHap_final, genomes = genomes_cont)
      
      
      output$cont_rmp_rhot_ID <- rhandsontable::renderRHandsontable(
        rhandsontable(estHap_rmp_database) %>% 
          hot_col("RmpFrequentist", format = "0.00000") %>%
          hot_col("RmpTheta", format = "0.00000") %>%
          hot_col("Theta", format="0.00000")
        
      )
      
      
      })
    })
    
    output$download_data_ID <- downloadHandler(
      filename <- function() {
        paste("SampleData", "zip", sep=".")
      },
      
      content <- function(file) {
        file.copy("SampleData.zip", file)
      },
      contentType = "application/zip"
    )

})
