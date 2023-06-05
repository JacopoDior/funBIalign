# CASE STUDY: Food_Infl - WORLD -------
library(readr)
library(fda)
library(tidyverse)
library(KernSmooth)
library(Rcpp)
library(dendextend)
library(Rcpp)

# Source functons and upload data
setwd("~/Desktop/funBIalign")
source("funBIalign_functions.R")

# world data
long_data <- read_csv("case_studies/data/food_price_inflation_world.csv")

View(long_data)

# PREPROCESSING -----
## Preprocessing: From Long to Wide ----
long_data <- long_data %>%
  dplyr::select(Area, Year, Months, Value) %>%
  mutate(Time = paste(Year, Months))

all_time <- lapply(unique(long_data$Year), function(x){paste(x,unique(long_data$Months))}) %>% unlist()
all_time <- cbind.data.frame(all_time, 1:length(all_time))
colnames(all_time) <- c("Time", "Order")

long_data <- merge(long_data, all_time, by="Time")
long_data <- long_data[order(long_data$Order),] %>%
  dplyr::select(Area, Order, Value)

wide_data <- pivot_wider(long_data, names_from = Order, values_from = Value)
area_names <- wide_data$Area
wide_data <- wide_data %>%
  dplyr::select(-Area)
rownames(wide_data) <- area_names
wide_data <- as.matrix(wide_data)
#save(wide_data, file = "food_price_inflation_world_wide.rds")


## Preprocessing: Smoothing ----
data <- t(wide_data)

# Local smoothing
foodinfl_loc <- c()
for(k in 1:dim(data)[2]){
  temp <- locpoly(1:dim(data)[1], data[,k], bandwidth = 1.5)
  temp <- temp$y[seq(1,length(temp$y), len = 258)]
  foodinfl_loc <- cbind(foodinfl_loc, temp)
}

# FUNBIALIGN ----
portion_len <- 12
min_card  <- 4

foodinfl_loc <- t(foodinfl_loc)
matplot(t(foodinfl_loc), type='l')

## STEP 1 and 2 -----
# step 1
window_data      <- create_lots(foodinfl_loc, portion_len)

# step 2: compute fMRS-based dissimilarity matrix
D_fmsr  <- create_distance(foodinfl_loc, window_data, portion_len)

## STEP 3 -----
# step 3: get the sub-trees (tree_s)
minidend         <- get_minidend(D_fmsr, window_data)

# step 3: identify seeds and corresponding families
all_paths   <- get_path_complete(minidend, window_data, min_card = min_card)
all_paths   <- all_paths[!(lapply(all_paths, is.null) %>% unlist())] 

# step 3: get recommended nodes and their info (cardinality and score)
# collect all recommended nodes (as an array of portion ids)
all_recommended_labels <- lapply(all_paths, function(x){x$recommended_node_labels})
list_of_recommendations <- lapply(rapply(
  all_recommended_labels, enquote, how='unlist'),
  eval)
# get recommended node cardinality
vec_of_card <- lapply(list_of_recommendations, length) %>% unlist()

# get vector of adjusted fMSR
vec_of_scores <- lapply(all_paths, function(x){x$recommended_node_scores}) %>% unlist()

## STEP 4 -----
# STEP 4: post-processing and rearranging results using different criteria

### CRITERIUM: adjusted fMSR ----
best_order <- vec_of_scores %>% order()
vec_of_scores_ordered <- vec_of_scores[best_order]
# ordered list_of_recommendations
list_of_recommendations_ordered <- list_of_recommendations[best_order]

#for every recommended motif, compute all its accolites
all_accolites <- lapply(list_of_recommendations_ordered,
                        function(x){
                          lapply(x, get_accolites, window_data, portion_len, FALSE) %>% unlist()
                        })

# 
delete <- c()
for(i in 2:length(list_of_recommendations_ordered)){ # compare every motif (from the second one)
  node_1 <- list_of_recommendations_ordered[[i]]
  for(j in 1:(i-1)){ # to the other higher ranked nodes
    node_2 <- list_of_recommendations_ordered[[j]]
    if(length(node_1) <= length(node_2)){ 
      accolites_2      <- all_accolites[[j]]
      common_elements  <- sum(node_1 %in% accolites_2)
      if((common_elements == length(node_1)) && vec_of_scores_ordered[i]>vec_of_scores_ordered[j]){
        delete <- c(delete,i)
        break()
      }
    }
  }
}

# perc_of_appearance
list_of_recommendations_ordered <- list_of_recommendations_ordered[-delete]
vec_of_scores_ordered <- vec_of_scores_ordered[-delete]


# Plot in the line
library(scales)
for(q in 1:length(list_of_recommendations_ordered)){
  temp_motif <- list_of_recommendations_ordered[[q]]
  lots_in_motif <- lapply(temp_motif,
                          function(x){
                            strsplit(x, '_') %>% 
                              unlist() %>%
                              as.numeric()}) %>% 
    unlist() %>%
    matrix(ncol=3, byrow=T)
  
  title   <- paste0("Number of instances: ", dim(plot_me)[2],
                    " - adj fMSR:  ", vec_of_scores_ordered[q] %>% round(3))
  
  matplot(t(foodinfl_loc), col = "grey40", ylab='', xlab='', axes='FALSE',
          lty = 1, type = 'l', main = title)
  for(k in 1:nrow(lots_in_motif)){
    rect(lots_in_motif[k,2], min(foodinfl_loc)-10, lots_in_motif[k,3], max(foodinfl_loc) +10,
         border = alpha("firebrick3", 0.2), col = alpha("firebrick3", 0.2))
    matplot(lots_in_motif[k,2]:lots_in_motif[k,3],
            foodinfl_loc[lots_in_motif[k,1],lots_in_motif[k,2]:lots_in_motif[k,3]], type='l', add=T, col='firebrick3', lwd = 2)
    box(col="grey40",lwd = 2)
    axis(1, at=seq(1, length(foodinfl_loc), by=12), labels = all_time$Time[seq(1, length(foodinfl_loc), by=12)], col="grey40", col.ticks="grey40", col.axis="grey60", cex.axis=1.5)
    axis(2, col="grey40", col.ticks="grey40", col.axis="grey60", cex.axis=1.5)
  }
}

# CRITERIUM: variance ----
motif_var <- lapply(list_of_recommendations_ordered, 
                    function(x){
                      temp <- window_data[x,]
                      temp_mean <- temp %>% colMeans()
                      ((temp - temp_mean)^2 %>% sum())/(nrow(temp))
                    }) %>% unlist()

var_order <- motif_var %>% order(decreasing = TRUE)
list_of_recommendations_ordered <- list_of_recommendations_ordered[var_order]
vec_of_scores_ordered <- vec_of_scores_ordered[var_order]











