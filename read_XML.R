require(XML)
data <- xmlParse("/Users/tegli/Downloads/full_database.xml")
xml_data <- xmlToList(data)

save(xml_data, file = "/Users/tegli/Desktop/FAERS/drugbank_list.RData")
length(xml_data[[1]])
str(xml_data[[1]])
names(xml_data[[1]][[1]])

names(xml_data) <- all_names
save(xml_data, file = "/Users/tegli/Desktop/FAERS/drugbank_list_named.RData")

all_names <- tolower(unname(sapply(xml_data[-length(xml_data)],function(x){
    x[["name"]]
}, USE.NAMES = F)))

all_descriptions <- sapply(xml_data[-length(xml_data)],function(x){
    x[["description"]]
}, USE.NAMES = F)

all_indications <- sapply(xml_data[-length(xml_data)],function(x){
    x[["indication"]]
}, USE.NAMES = F)

all_pharmacodynamics <- sapply(xml_data[-length(xml_data)],function(x){
    x[["pharmacodynamics"]]
}, USE.NAMES = F)


info_list <- lapply(xml_data[-length(xml_data)],function(x){
    targets_list <- lapply(x[["targets"]],function(y){
        
        if(!is.null(y$polypeptide$`go-classifiers`)){
            go_classifiers <- sort(unname(sapply(y$polypeptide$`go-classifiers`,function(go_tmp){
                if("description" %in% names(go_tmp)){
                    go_tmp$description
                }else{
                    ""
                }
            }, USE.NAMES = F)))
        }else{
            go_classifiers <- ""
        }
        
        name_tmp <- ifelse(is.null(y$name),"",y$name)
        action_tmp <- ifelse(is.null(y$actions),"",y$actions)
        general_function <- ifelse(is.null(y$polypeptide$`general-function`),
            "",y$polypeptide$`general-function`)
        cellular_location <- ifelse(is.null(y$polypeptide$`cellular-location`),
            "",y$polypeptide$`cellular-location`)
        target_gene <- ifelse(is.null(y$polypeptide$`gene-name`),
            "",y$polypeptide$`gene-name`)
        target_synonyms <- ifelse(is.null(y$polypeptide$synonyms),
            "",y$polypeptide$synonyms)
            
        data.frame(name = name_tmp,
            actions = paste(action_tmp, collapse = "; "),
            target_function = general_function,
            cellular_location = cellular_location,
            target_gene = target_gene,
            target_synonyms = paste(target_synonyms, collapse = "; "),
            go_classifiers = paste(go_classifiers, collapse = "; "),
            stringsAsFactors = F)
    })
    do.call("rbind",targets_list)
})

save(info_list, file = "/Users/tegli/Desktop/FAERS/drugbank_list_processedList.RData")
