library(shiny)

drugnames_unsorted <- read.table(paste0("data/results_memory.txt"),
        stringsAsFactors = F,
        sep = "\t",
        h = T)[,"results.drug"]

data <- read.table(paste0("data/results_memory_sorted.txt"),
        stringsAsFactors = F,
        sep = "\t",
        h = T)

data <- data[!is.na(data$results.ROR),]
p_crit <- 0.05/(nrow(data)*16)

all_drugs <- unique(data$results.drug)

load("data/assignments_list.RData")
load("data/drugbank_list_processedList.RData")

all_genes <- sort(unique(unlist(sapply(info_list,function(x){
   unique(c(x$targets$target_gene, x$targets$target_synonyms, x$targets$name))
}))))

all_genes <- all_genes[!all_genes %in% c("","\n    ","1.-.-.-")]

shinyServer(function(input, output, session) {

    newData_allCases <- reactive({
        #-------------------------- Read in sorted data for all cases
        data <- read.table(paste0("data/results_",
                input$result_name,"_sorted.txt"),
            stringsAsFactors = F,
            sep = "\t",
            h = T)
        data <- data[!is.na(data$results.ROR),]
        data
    })
    
    newData_allCases_noNA <- reactive({
        #-------------------------- Read in sorted data for all cases
        data <- read.table(paste0("data/results_",
                input$result_name,"_sorted.txt"),
            stringsAsFactors = F,
            sep = "\t",
            h = T)
        
        data <- data[!is.na(data$results.ROR),]

        for(i_col in grep(".pval",names(data), fixed = T)){
            data[data[,i_col] == 0 & !is.na(data[,i_col]), i_col] <- exp(-500)
            data[is.na(data[,i_col]),i_col] <- 1
        }
        
        for(i_col in grep(".ROR",names(data), fixed = T)){
            data[data[,i_col] == 0 & !is.na(data[,i_col]), i_col] <- 0.0000000000001
            data[is.na(data[,i_col]),i_col] <- 1
        }
        data
    })

    newRange <- reactive({
        #-------------------------- Get range from sliders
        range_out <- round(input$range)
        range_out
    })
    
    newRangeY <- reactive({
        #-------------------------- Get range from sliders
        range_out <- round(input$rangeY)
        range_out
    })
    
    newIndication_list <- reactive({
        load(paste0("data/result_lists_",input$result_name,".RData"))
        indication_list
    })
    
    newIndication_drug_list <- reactive({
        load(paste0("data/result_lists_",input$result_name,".RData"))

        names(indication_drug_list) <- drugnames_unsorted
        indication_drug_list
    })
    
    observe({
        if(!is.null(input$pvalPlot_click$x)){
            click_use <- round(input$pvalPlot_click$x)
            click_use <- ifelse(click_use < 1,1, click_use)

            data <- newData_allCases_noNA() # is sorted
 
            if(input$plotReference != "Indications"){
                columnname_plot_use <- gsub("Order p ","",input$plotReference)
            
                if(columnname_plot_use == "all"){
                    columnname_plot_use <- "results"
                }else{
                    columnname_plot_use <- paste0("results_",columnname_plot_use)
                }
                
                logPval_plot_use <- (-log10(data[,paste0(columnname_plot_use,".pval")])*ifelse(data[,paste0(columnname_plot_use,".ROR")] < 1,-1,1))
                order_use <- order(logPval_plot_use)
                #order_use <- c(1:length(all_drugs))
            }else{
                order_use <- c(1:length(all_drugs))
            }

            updateSelectInput(session, inputId = "drug",
                selected = all_drugs[order_use][click_use])
        }else if(!is.null(input$ORPlot_click$x)){
            click_use <- round(input$ORPlot_click$x)
            click_use <- ifelse(click_use < 1, 1, click_use)
            
            data <- newData_allCases_noNA() # is sorted

            if(input$plotReference != "Indications"){
                columnname_plot_use <- gsub("Order p ","",input$plotReference)

                if(columnname_plot_use == "all"){
                    columnname_plot_use <- "results"
                }else{
                    columnname_plot_use <- paste0("results_",columnname_plot_use)
                }

                logPval_plot_use <- (-log10(data[,paste0(columnname_plot_use,".pval")])*ifelse(data[,paste0(columnname_plot_use,".ROR")] < 1,-1,1))
                order_user <- order(logPval_plot_use)
                #order_use <- c(1:length(all_drugs))
            }else{
                order_use <- c(1:length(all_drugs))
            }
            updateSelectInput(session, inputId = "drug",
                selected = all_drugs[order_use][click_use])
        }
    })
    
    observe({
        if(!is.null(input$pvalPlot_dblclick$x)){
            updateSliderInput(session, inputId = "range",
                value = c(1,nrow(data)))
        }else if(!is.null(input$ORPlot_dblclick$x)){
            updateSliderInput(session, inputId = "range",
                value = c(1,nrow(data)))
        }
    })

    observe({
        if(!is.null(input$pvalPlot_brush$xmin)){
            updateSliderInput(session, inputId = "range",
                value = c(round(input$pvalPlot_brush$xmin),
                    round(input$pvalPlot_brush$xmax)))
        }else if(!is.null(input$ORPlot_brush$xmin)){
            updateSliderInput(session, inputId = "range",
                value = c(round(input$ORPlot_brush$xmin),
                    round(input$ORPlot_brush$xmax)))
        }
    })
    
    output$ORPlot <- renderPlot({
        #-------------------------- Plot the count of adverse events (y-axis) per active ingredient (x-axis)
        
        if(input$result_name == "anorexia"){
            reac_print <- "eating disorder"
        }else if(input$result_name == "anxiety"){
            reac_print <- "anxiety"
        }else if(input$result_name == "attention"){
            reac_print <- "disturbed attention"
        }else if(input$result_name == "dementia"){
            reac_print <- "dementia"
        }else if(input$result_name == "depression"){
            reac_print <- "depression"
        }else if(input$result_name == "emotional"){
            reac_print <- "disturbed emotions"
        }else if(input$result_name == "mania"){
            reac_print <- "bipolar/mania"
        }else if(input$result_name == "memory"){
            reac_print <- "impaired memory"
        }else if(input$result_name == "panic_attack"){
            reac_print <- "panic attack"
        }else if(input$result_name == "paranoia"){
            reac_print <- "paranoia"
        }else if(input$result_name == "psychotic"){
            reac_print <- "psychotic disorder"
        }else if(input$result_name == "suicide"){
            reac_print <- "suicide"
        }else if(input$result_name == "anaemia"){
            reac_print <- "anaemia"
        }else if(input$result_name == "arthralgia"){
            reac_print <- "arthralgia"
        }else if(input$result_name == "nasopharyngitis"){
            reac_print <- "nasopharyngitis"
        }else if(input$result_name == "pneumonia"){
            reac_print <- "pneumonia"
        }
        
        data <- newData_allCases_noNA() # is sorted

        range_use <- newRange()
        
        if(input$type == "all"){
            columnname_use <- "results"
        }else{
            columnname_use <- paste0("results_",input$type)
        }
        
        if(input$plotReference != "Indications"){
            columnname_plot_use <- gsub("Order p ","",input$plotReference)

            if(columnname_plot_use == "all"){
                columnname_plot_use <- "results"
            }else{
                columnname_plot_use <- paste0("results_",columnname_plot_use)
            }

            logPval_plot_use <- -log10(data[,paste0(columnname_plot_use,".pval")])*ifelse(data[,paste0(columnname_plot_use,".ROR")] < 1,-1,1)
            data <- data[order(logPval_plot_use),]
        }

        data[,paste0(columnname_use,".ROR")] <- log(data[,paste0(columnname_use,".ROR")])
        
        plot(data[,paste0(columnname_use,".ROR")],
            bty = "n",
            main = paste0("log(Odds ratio) for ",reac_print," and ",input$drug," intake"),
            ylab = "log(Odds ratio)",
            type = "n",
            cex = 0.75,
            ylim = c(-5,5),
            xlim = range_use,
            xaxt = "n",
            xlab = "Active Ingredients")
        if(input$plotReference == "Indications"){
            abline(v = data$count_change, col = "lightgrey", xpd = F)
        }

        points(x = 1:length(data[,paste0(columnname_use,".ROR")]),
            y = data[,paste0(columnname_use,".ROR")],
            pch = 20,
            cex = 0.75,
            col = "grey")
        
        points(x = which(data[,paste0(columnname_use,".pval")] < p_crit),
            y = data[data[,paste0(columnname_use,".pval")] < p_crit &
                        !is.na(data[,paste0(columnname_use,".pval")]),
                    paste0(columnname_use,".ROR")],
            pch = 20,
            col = "blue",
            cex = 0.75)
        
        points(x = which(data$results.drug == input$drug),
            y = data[,paste0(columnname_use,".ROR")][which(data$results.drug == input$drug)],
            pch = 20,
            col = "red")
        
        lines(x = rep(which(data$results.drug == input$drug),2),
            y = c(0,data[,paste0(columnname_use,".ROR")][which(data$results.drug == input$drug)]),
            col = "red")
        abline(h = 0, col = "grey")
        
        axis(1)
    })
    
    output$pvalPlot <- renderPlot({
        #-------------------------- Plot -log10(p) (y-axis) per active ingredient (x-axis)
        
        data <- newData_allCases_noNA()
        
        range_use <- newRange()
        rangeY_use <- newRangeY()
        
        if(input$type == "all"){
            columnname_use <- "results"
        }else{
            columnname_use <- paste0("results_",input$type)
        }
        
        if(input$plotReference != "Indications"){
            columnname_plot_use <- gsub("Order p ","",input$plotReference)

            if(columnname_plot_use == "all"){
                columnname_plot_use <- "results"
            }else{
                columnname_plot_use <- paste0("results_",columnname_plot_use)
            }

            logPval_plot_use <- -log10(data[,paste0(columnname_plot_use,".pval")])*ifelse(data[,paste0(columnname_plot_use,".ROR")] < 1,-1,1)
            data <- data[order(logPval_plot_use),]
        }

        logPval_use <- -log10(data[,paste0(columnname_use,".pval")])*ifelse(data[,paste0(columnname_use,".ROR")] < 1,-1,1)
        
        plot(logPval_use, # $results.pval
            pch = 20,
            col = "grey",
            ylim = rangeY_use,
            xlim = range_use,
            bty = "n",
            main = paste0("-log10(p) Fisher's exact test (",input$result_name," within ",columnname_use,")"),
            ylab = "-log10(p)",
            cex = 0.75,
            xaxt = "n",
            xlab = "Active ingredients")
        
        if(input$plotReference == "Indications"){
            abline(v = data$count_change, col = "lightgrey", xpd = F)
        }

        points(x = which(data[,paste0(columnname_use,".pval")] < p_crit),
            y = logPval_use[data[,paste0(columnname_use,".pval")] < p_crit &
                    !is.na(data[,paste0(columnname_use,".pval")])],
            pch = 20,
            col = "blue",
            cex = 0.75)
        
        points(x = which(data$results.drug == input$drug),
               y = logPval_use[which(data$results.drug == input$drug)],
               pch = 20,
               col = "red")
        
        lines(x = rep(which(data$results.drug == input$drug),2),
            y = c(0,logPval_use[which(data$results.drug == input$drug)]),
            col = "red")
        
        abline(h = -log10(p_crit), col = "red")
        abline(h = log10(p_crit), col = "red")
        abline(h = 0, col = "grey")
        axis(1)
    })

    output$drugPlot <- renderPlot({
        #-------------------------- Plot XXXXX plots for all odds ratios
        data <- newData_allCases()
        data_drug <- data[data$results.drug == input$drug,]
        
        data_drug$results.ROR[data_drug$results.ROR == 0] <- 0.0000000000001
        data_drug$results_indication.ROR[data_drug$results_indication.ROR == 0] <- 0.0000000000001
        data_drug$results_age.ROR[data_drug$results_age.ROR == 0] <- 0.0000000000001
        data_drug$results_female.ROR[data_drug$results_female.ROR == 0] <- 0.0000000000001
        data_drug$results_male.ROR[data_drug$results_male.ROR == 0] <- 0.0000000000001
        
        data_drug$results.CI_low[data_drug$results.CI_low == 0 | data_drug$results.CI_low == -Inf] <- 0.0000000000001
        data_drug$results_indication.CI_low[data_drug$results_indication.CI_low == 0 | data_drug$results_indication.CI_low == -Inf] <- 0.0000000000001
        data_drug$results_age.CI_low[data_drug$results_age.CI_low == 0 | data_drug$results_age.CI_low == -Inf] <- 0.0000000000001
        data_drug$results_female.CI_low[data_drug$results_female.CI_low == 0 | data_drug$results_female.CI_low == -Inf] <- 0.0000000000001
        data_drug$results_male.CI_low[data_drug$results_male.CI_low == 0 | data_drug$results_male.CI_low == -Inf] <- 0.0000000000001

        data_drug$results.CI_up[data_drug$results.CI_up == Inf] <- 100
        data_drug$results_indication.CI_up[data_drug$results_indication.CI_up == Inf] <- 100
        data_drug$results_age.CI_up[data_drug$results_age.CI_up == Inf] <- 100
        data_drug$results_female.CI_up[data_drug$results_female.CI_up == Inf] <- 100
        data_drug$results_male.CI_up[data_drug$results_male.CI_up == Inf] <- 100

        oldpar <- par()
        par(mar = c(5.1,10,0.1,0))
        
        plot(x = c(log(data_drug$results.ROR),
            log(data_drug$results_indication.ROR),
            log(data_drug$results_age.ROR),
            log(data_drug$results_female.ROR),
            log(data_drug$results_male.ROR)),
            y = c(5:1),
            pch = 20,
            bty = "n",
            xlim = c(-5,5),
            main = "",
            ylab = "",
            xlab = "log(Odds ratios)",
            yaxt = "n")

        abline(v = 0, col = "grey")

        arrows(y0 = 5, y1 = 5,
            x0 = log(data_drug$results.CI_low),
            x1 = log(data_drug$results.CI_up),
            angle = 90,
            length = 0.1,
            code = 3)

        arrows(y0 = 4, y1 = 4,
            x0 = log(data_drug$results_indication.CI_low),
            x1 = log(data_drug$results_indication.CI_up),
            angle = 90,
            length = 0.1,
            code = 3)

        arrows(y0 = 3, y1 = 3,
            x0 = log(data_drug$results_age.CI_low),
            x1 = log(data_drug$results_age.CI_up),
            angle = 90,
            length = 0.1,
            code = 3)

        arrows(y0 = 2, y1 = 2,
            x0 = log(data_drug$results_female.CI_low),
            x1 = log(data_drug$results_female.CI_up),
            angle = 90,
            length = 0.1,
            code = 3)

        arrows(y0 = 1, y1 = 1,
            x0 = log(data_drug$results_male.CI_low),
            x1 = log(data_drug$results_male.CI_up),
            angle = 90,
            length = 0.1,
            code = 3)

        axis(2, at = c(5:1),
            labels = c("Across all",
                paste0("Within ",data_drug$results_indication.indication),
                paste0("... & aged ",data_drug$results_age.min_age,"-",data_drug$results_age.max_age),
                "... & female",
                "... & male"),
            las = 2)
    })
    
    output$textResults <- renderUI({
        #-------------------------- Print short description text
        data <- newData_allCases()
        data_drug <- data[data$results.drug == input$drug,]
        N_total <- data_drug$results.DrugReac + data_drug$results.DrugNoReac

        N_indication <- data_drug$results_indication.DrugReac + data_drug$results_indication.DrugNoReac

        text_print <- paste0("The active ingredient <b>",data_drug$results.drug,"</b> was administered in ",
        	N_total," reports. The most frequent indication for ",data_drug$results.drug,
        	" was <b>",data_drug$results_indication.indication,"</b> (", N_indication," cases).\n",
            "Results within indication, age, and separate genders were compared to other entries with indication ",
            data_drug$results_indication.indication,".")
        HTML(text_print)
    })
    
    output$fishersTable <- renderTable({
        #-------------------------- Print cross table drug ~ reaction
        data <- newData_allCases()
        data_drug <- data[data$results.drug == input$drug,,drop = F]

        if(input$type == "all"){
            columnname_use <- "results"
        }else{
            columnname_use <- paste0("results_",input$type)
        }
        
        #data_drug[data_drug[,paste0(columnname_use,".pval")] == 0,
        #    paste0(columnname_use,".pval")] <- exp(-700)
        
        pval <- data_drug[,paste0(columnname_use,".pval")]
        OR <- data_drug[,paste0(columnname_use,".ROR")]
        CI_low <- data_drug[,paste0(columnname_use,".CI_low")]
        CI_up <- data_drug[,paste0(columnname_use,".CI_up")]
        CI <- paste0(CI_low,"-",CI_up)
        log_OR <- ifelse(OR == 0, log(0.0000000000001), log(OR))

        pval <- sapply(c(pval),function(x){
            if(x > 0.01){
                return(as.character(round(x,2)))
            }else if(x > 0.001){
                return(as.character(round(x,3)))
            }else if(x > 0.0001){
                return(as.character(round(x,4)))
            }else{
                return(as.character(signif(x,3)))
            }
        })

        data_mat <- matrix(c(OR,CI,round(log_OR,2),pval),ncol = 1)
        
        colnames(data_mat) <- c("value")
        rownames(data_mat) <- c("OR","95% CI","log(OR)","p")
        data_mat
    },
    include.rownames = T)
    
    output$crossTable <- renderTable({
        #-------------------------- Print cross table drug ~ reaction
        data <- newData_allCases()
        data_drug <- data[data$results.drug == input$drug,]

        if(input$type == "all"){
            data_drug <- data_drug[,grepl("results.Drug|results.NoDrug",names(data_drug)),drop = F]
        }else if(input$type == "indication"){
            data_drug <- data_drug[,grepl("results_indication.Drug|results_indication.NoDrug",names(data_drug)),drop = F]
        }else if(input$type == "age"){
            data_drug <- data_drug[,grepl("results_age.Drug|results_age.NoDrug",names(data_drug)),drop = F]
        }else if(input$type == "female"){
            data_drug <- data_drug[,grepl("results_female.Drug|results_female.NoDrug",names(data_drug)),drop = F]
        }else if(input$type == "male"){
            data_drug <- data_drug[,grepl("results_male.Drug|results_male.NoDrug",names(data_drug)),drop = F]
        }
        data_mat <- matrix(unlist(data_drug),ncol = 2, byrow = T)

        perc_reacDrug <- data_mat[1,1]/(sum(as.numeric(data_mat[1,])))
        perc_reacDrug <- as.character(round(100*perc_reacDrug,2))
        perc_reacNoDrug <- data_mat[2,1]/(sum(as.numeric(data_mat[2,])))
        perc_reacNoDrug <- as.character(round(100*perc_reacNoDrug,2))

        perc_drugReac <- data_mat[1,1]/(sum(as.numeric(data_mat[,1])))
        perc_drugReac <- as.character(round(100*perc_drugReac,2))
        perc_drugNoReac <- data_mat[1,2]/(sum(as.numeric(data_mat[,2])))
        perc_drugNoReac <- as.character(round(100*perc_drugNoReac,2))

        data_mat <- cbind(data_mat,c(perc_reacDrug,perc_reacNoDrug))
        data_mat <- rbind(data_mat,c(perc_drugReac,perc_drugNoReac,""))

        colnames(data_mat) <- c("reac","no reac","% reac")
        rownames(data_mat) <- c("drug","no drug","% drug")
        data_mat
    },
    include.rownames = T)
    
    output$logregResults <- renderTable({
        #-------------------------- Print result table for logistic regression
        data <- newData_allCases()
        data_drug <- data[data$results.drug == input$drug,]
        data_drug <- data_drug[,grepl("logReg",names(data_drug)),drop = F]
        if(!is.na(data_drug[,4])){
        
            N_used <- as.character(data_drug[,14])
            N_used <- substr(N_used,2,nchar(N_used)-1)
            data_drug <- data_drug[,-c(1:3,5,14)]
            data_drug <- sapply(data_drug,function(x){
                if(x > 0.01){
                    return(as.character(round(x,2)))
                }else if(x > 0.001){
                    return(as.character(round(x,3)))
                }else if(x > 0.0001){
                    return(as.character(round(x,4)))
                }else{
                    return(as.character(signif(x,3)))
                }
            })
            data_drug <- matrix(data_drug,ncol = 3, byrow = T)
            data_drug <- cbind(data_drug,rep(N_used,3))

        }else{
            data_drug <- rep(NA,12)
            data_drug <- matrix(data_drug,ncol = 4, byrow = T)
        }
        colnames(data_drug) <- c("log(OR)","SE","p","N")
        rownames(data_drug) <- c("drug","age","sex")
        data_drug
    }, include.rownames = T)
    
    output$drugActiveIngredients <- renderDataTable({
        #-------------------------- Print table of all drugnames assigned to the active ingredient
        drugActiveIngredients_table <- assignments_list[[input$drug]]

        drugActiveIngredients_df <- data.frame(Drugname = names(drugActiveIngredients_table),
            Count = as.numeric(drugActiveIngredients_table),
            stringsAsFactors = F)
        drugActiveIngredients_df
    })
    
    output$drugIndications <- renderDataTable({
        #-------------------------- Print table of the 10 top indications of the active ingredient
        indication_list <- newIndication_list()
        
        data <- newData_allCases()
        N_use <- data[data$results.drug == input$drug,"results.N"]
        
        indication_vect <- indication_list[[N_use]]
        indication_df <- data.frame(Indications = names(indication_vect),
            Count = as.numeric(indication_vect),
            stringsAsFactors = F)
        indication_df
    })
    
    output$topIndicationsDrug <- renderDataTable({
        #-------------------------- Print table of the 10 most administered drugnames for the
        #                           top indication
        indication_drug_list <- newIndication_drug_list()
        
        data <- newData_allCases()
        N_use <- data[data$results.drug == input$drug,"results.N"]
        
        drug_vect <- indication_drug_list[[N_use]]
        
        drugActiveIngredients<- names(assignments_list[[input$drug]])
        
        names(drug_vect)[names(drug_vect) %in% drugActiveIngredients] <- paste0(names(drug_vect)[names(drug_vect) %in% drugActiveIngredients],"*")
        
        drug_df <- data.frame(Drugname = names(drug_vect),
            Count = unname(drug_vect),
            stringsAsFactors = F)
        drug_df
    })
    
    output$infoDrug <- renderTable({
        if(input$drug %in% names(info_list)){
            use_name <- input$drug
            df_use <- info_list[[use_name]]$targets
            df_use[,names(df_use) != "go_classifiers"]
        }else if(any(names(assignments_list[[input$drug]]) %in% names(info_list))){
            name_ind <- names(assignments_list[[input$drug]]) %in% names(info_list)
            use_name <- names(assignments_list[[input$drug]])[name_ind]
            df_use <- info_list[[use_name]]$targets
            df_use[,names(df_use) != "go_classifiers"]
        }else{
            data.frame(name = "",
                actions = "",
                target_function = "",
                cellular_location = "",
                target_gene = "",
                target_synonyms = "",
                stringsAsFactors = F)
        }
    })
    
    output$infoDrugDescription <- renderUI({
        if(input$drug %in% names(info_list)){
            use_name <- input$drug
            text_print <- info_list[[use_name]]$description
        }else if(any(names(assignments_list[[input$drug]]) %in% names(info_list))){
            name_ind <- names(assignments_list[[input$drug]]) %in% names(info_list)
            use_name <- names(assignments_list[[input$drug]])[name_ind]
            text_print <- info_list[[use_name]]$description
        }else{
            text_print <- ""
            HTML(text_print)
        }
    })
    
    output$infoDrugPharmacodynamic <- renderUI({
        if(input$drug %in% names(info_list)){
            use_name <- input$drug
            text_print <- info_list[[use_name]]$pharmacodynamic
        }else if(any(names(assignments_list[[input$drug]]) %in% names(info_list))){
            name_ind <- names(assignments_list[[input$drug]]) %in% names(info_list)
            use_name <- names(assignments_list[[input$drug]])[name_ind]
            text_print <- info_list[[use_name]]$pharmacodynamic
        }else{
            text_print <- ""
            HTML(text_print)
        }
    })
    
    output$infoDrugIndication <- renderUI({
        if(input$drug %in% names(info_list)){
            use_name <- input$drug
            text_print <- info_list[[use_name]]$indication
        }else if(any(names(assignments_list[[input$drug]]) %in% names(info_list))){
            name_ind <- names(assignments_list[[input$drug]]) %in% names(info_list)
            use_name <- names(assignments_list[[input$drug]])[name_ind]
            text_print <- info_list[[use_name]]$indication
        }else{
            text_print <- ""
            HTML(text_print)
        }
    })
    
    output$GOannot <- renderDataTable({
        if(input$drug %in% names(info_list)){
            use_name <- input$drug
            df_use <- info_list[[use_name]]$targets
            data.frame(annotation = strsplit(df_use$go_classifiers, "; ",fixed = T)[[1]], stringsAsFactors = F)
        }else if(any(names(assignments_list[[input$drug]]) %in% names(info_list))){
            name_ind <- names(assignments_list[[input$drug]]) %in% names(info_list)
            use_name <- names(assignments_list[[input$drug]])[name_ind]
            df_use <- info_list[[use_name]]$targets
            data.frame(annotation = strsplit(df_use$go_classifiers, "; ",fixed = T)[[1]], stringsAsFactors = F)
        }else{
            data.frame(annotation = "",
                stringsAsFactors = F)
        }
    })
    
    output$targetGenes <- renderDataTable({
        if(input$selectGene %in% all_genes){
            target_ind <- sapply(info_list,function(x){
                targets_tmp <- unique(c(x$targets$target_gene, x$targets$target_synonyms, x$targets$name))
                input$selectGene %in% targets_tmp
            })
            
            names_targets <- names(info_list)[target_ind]

            data <- newData_allCases_noNA()

            result_list <- lapply(names_targets,function(x){
                df_tmp <- data[grepl(x,data$results.drug,ignore.case = T, fixed = T),c("results.drug",
                    "results.ROR",
                    "results.CI_low",
                    "results.CI_up",
                    "results.pval",
                    "results_indication.indication",
                    "results_indication.ROR",
                    "results_indication.CI_low",
                    "results_indication.CI_up",
                    "results_indication.pval")]
                if(nrow(df_tmp) != 0){
                    df_tmp <- as.data.frame(cbind(rep(x, nrow(df_tmp)),df_tmp),stringsAsFactors = F)
                }else{
                    df_tmp <- data.frame(x,"","","","","","","","","","")
                }

                for(i in c(3:5,8:10)){
                    df_tmp[,i] <- round(as.numeric(df_tmp[,i]),2)
                }

                for(i in c(6,11)){
                    df_tmp[,i] <- signif(as.numeric(df_tmp[,i]),2)
                }

                names(df_tmp) <- c("Ingr. DB",
                    "Ingr. FAERS",
                    "OR all",
                    "CI low all",
                    "CI up all",
                    "p all",
                    "indication",
                    "OR ind",
                    "CI low ind",
                    "CI up ind",
                    "p ind")
                df_tmp
            })

            do.call("rbind",result_list)
        }else{
            df_tmp <- data.frame(x,"","","","","","","","","","")
            names(df_tmp) <- c("Ingr. DB",
                "Ingr. FAERS",
                "OR all",
                "CI low all",
                "CI up all",
                "p all",
                "indication",
                "OR ind",
                "CI low ind",
                "CI up ind",
                "p ind")
            df_tmp
        }
    })


    output$downloadData <- downloadHandler(
        filename = function(){
            paste0(input$result_name,"_",
                gsub(" ","_",input$drug),"_FishersExactTest.txt")
        },
        content = function(file){
            data <- newData_allCases()
            data_drug <- data[data$results.drug == input$drug,,drop = F]
            names_vars <- c("ROR","CI_low","CI_up","pval","DrugReac","DrugNoReac","NoDrugReac","NoDrugNoReac")
            names_types <- c("results","results_indication","results_age","results_male","results_female")

            names_out <- sapply(names_types,function(x){
                paste0(x,".",names_vars)
            })

            data_out <- apply(names_out,2,function(x){
                mat <- data_drug[,x]
                names(mat) <- names_vars
                mat
            })
            
            data_out <- do.call("rbind",data_out)
            data_out <- cbind(c("all","within_indication","within_age_of_Drug","within_male","within_female"),data_out)
            colnames(data_out)[1] <- "analysis"

            write.table(data_out,
                file,
                quote = F,
                sep = "\t",
                row.names = F)
        }
    )

    output$downloadDataLog <- downloadHandler(filename = function(){
            paste0(input$result_name,"_",
                gsub(" ","_",input$drug),"_LogReg.txt")
        },
        content = function(file){
            data <- newData_allCases()
            data_drug <- data[data$results.drug == input$drug,]
            data_drug <- data_drug[,grepl("logReg",names(data_drug)),drop = F]
            N_used <- as.character(data_drug[,14])
            N_used <- substr(N_used,2,nchar(N_used)-1)
            data_drug <- data_drug[,-c(1:3,5,14)]
            data_drug <- sapply(data_drug,function(x){
                if(x > 0.01){
                    return(as.character(round(x,2)))
                }else if(x > 0.001){
                    return(as.character(round(x,3)))
                }else if(x > 0.0001){
                    return(as.character(round(x,4)))
                }else{
                    return(as.character(signif(x,3)))
                }
            })
            data_drug <- matrix(data_drug,ncol = 3, byrow = T)
            data_drug <- cbind(data_drug,rep(N_used,3))
            colnames(data_drug) <- c("log(OR)","SE","p","N")
            rownames(data_drug) <- c("drug","age","sex")

            write.table(data_drug,
                file,
                quote = F,
                sep = "\t")
        } 
    )
})
