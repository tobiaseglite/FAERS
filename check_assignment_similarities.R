data <- read.table(paste0("results_memory.txt"),
        stringsAsFactors = F,
        sep = "\t",
        h = T)

drugnames_unsorted <- data[,"results.drug"]

load("assignments_list.RData")

assignments_list <- assignments_list[data$results.N]

result_list <- lapply(1:length(assignments_list),function(x){
    data_tmp <- assignments_list[[x]]
    name_tmp <- names(assignments_list)[x]
    data.frame(activeIngredient = rep(name_tmp,length(data_tmp)),
        drugname = names(data_tmp),
        count = as.numeric(data_tmp),
        stringsAsFactors = F)
})
result_df <- do.call("rbind",result_list)
result_df <- result_df[result_df$count > 1000,]

sorensen_dice <- function(x1,x2){
    (2*length(intersect(x1,x2)))/(length(x1) + length(x2))
}

ingredients <- unique(result_df$activeIngredient)
n_ingredients <- length(ingredients)
result_mat <- matrix(rep(NA,3*n_ingredients^2), ncol = 3)

counter <- 0
for(x in 1:n_ingredients){
    for(y in 1:n_ingredients){
        if(y <= x){
            next
        }
        counter <- counter+1
        print(counter)
        result_mat[counter,1] <- x
        result_mat[counter,2] <- y
        result_mat[counter,3] <- round(sorensen_dice(result_df[result_df$activeIngredient == ingredients[x], "drugname"],
            result_df[result_df$activeIngredient == ingredients[y], "drugname"]),4)
    }
}
result_mat <- result_mat[!is.na(result_mat[,1]),]
save(result_mat, file = "result_lists_emotional_drugSimilarity_all.RData")
result_mat <- result_mat[result_mat[,3] != 0,]
colnames(result_mat) <- c("from","to","weight")

nodes_df <- data.frame(id = c(1:n_ingredients), ingredient = ingredients,
    stringsAsFactors = F)

result_mat_plot <- result_mat

result_mat_plot <- merge(x = result_mat_plot, y = nodes_df,
    by.x = "from", by.y = "id", all.x = T)
result_mat_plot <- merge(x = result_mat_plot, y = nodes_df,
    by.x = "to", by.y = "id", all.x = T)
result_mat_plot <- result_mat_plot[order(result_mat_plot[,1],result_mat_plot[,2]),]

head(result_mat_plot)

assignments_list[["niacinamide"]][assignments_list[["niacinamide"]] > 1000]
assignments_list[["pyridoxine"]][assignments_list[["pyridoxine"]] > 1000]



library(igraph)
net <- graph.data.frame(result_mat_plot, nodes_df, directed=T)
net

net <- simplify(net, remove.multiple = F, remove.loops = T)

quartz()
plot(net,
    vertex.shape="none",
    vertex.label=V(net)$ingredient,
    vertex.label.font=2,
    vertex.label.color="gray40",
    vertex.label.cex=.5,
    edge.color="gray85")
