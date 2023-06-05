# CASE STUDY: Temperature Change - 19 REGIONS -------
library(tidyverse)
library(fda)
library(scales)
library(magic)

# #FUNCTIONS AND DATA UPLOAD -----
# perform GCV b-spline smoothing
smooth_curves <- function(data, degree = 3, weights = NULL, breaks = seq_len(nrow(data)), 
                          daytime = seq_len(nrow(data))){
  # create basis functions
  spline_basis = create.bspline.basis(norder = degree+1, breaks = breaks)
  # chose lambda by GCV
  lambda = 10^seq(-4, 4, by = 0.2)
  rough_fd = vector("list", length(lambda))
  GCV = matrix(NA, nrow = length(lambda), ncol = ncol(data))
  for(i in 1:length(lambda)){
    rough_i = smooth.basis(daytime, data, fdParobj = fdPar(fdobj = spline_basis, Lfdobj = 2, lambda = lambda[i]), wtvec=weights)
    rough_fd[[i]] = rough_i$fd
    GCV[i,] = rough_i$gcv
    
    # re-smooth the curves with NA
    index_na = which(is.na(GCV[i,]))
    for(j in index_na){
      times_not_na = !is.na(data[,j])
      rough_i = smooth.basis(daytime[times_not_na], data[times_not_na,j], fdParobj = fdPar(fdobj = spline_basis, Lfdobj = 2, lambda = lambda[i]))
      rough_fd[[i]]$coefs[,j] = rough_i$fd$coefs
      GCV[i,j] = rough_i$gcv
    }
  }
  
  # GCV_all_curves
  meanGCV = rowMeans(GCV, na.rm = T)
  index_selected = which.min(meanGCV)
  lambda_selected = lambda[index_selected]
  GCV_all_curves = rough_fd[[index_selected]]
  
  plot(log10(lambda), meanGCV, type = 'l', xlab = 'log10(lambda)', ylab = 'Mean GCV', main = 'GCV all curves')
  abline(v = log10(lambda_selected), lty = 2)
  
  return(GCV_all_curves)
}

# Source functions and upload data
setwd("~/Desktop/funBIalign")
source("funBIalign_functions.R")

# 19 regions data
mydata <- read_csv("case_studies/data/FAOSTAT_tempchange_region.csv")

# PREPROCESSING -----
## Preprocessing: From Long to Wide ----

# select only interesting columns
long_data <- mydata[, c(4,8,10,12)] 

long_data <- long_data %>%
  dplyr::select(Area, Year, Months, Value) %>%
  mutate(Time = paste(Year,Months))

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

## Preprocessing: b-splines smoothing -------
# Smooth the data or...
temp <- smooth_curves(t(wide_data))
data_smooth <- t(eval.fd(1:(dim(wide_data)[2]), temp)) #(it might be time consuming)
data_smooth[which(is.na(wide_data)),] <- NA
colnames(data_smooth) <- 1:ncol(data_smooth)

# alternatively, load the smoothed data directly
load('case_studies/data/FAOSTAT_tempchange_region_smooth.RData')

## Plotting the data -----
matplot(t(data_smooth), col = alpha("grey40", 0.45), ylab='', xlab='', axes='FALSE', type = 'l')
box(col="grey40",lwd = 2)
axis(1, at=seq(1,ncol(data_smooth), by = 24), labels = all_time$Time[seq(1, ncol(data_smooth), by = 24)], col="grey50", col.ticks="grey50", col.axis="grey60", cex.axis=1)
axis(2, col="grey40", col.ticks="grey50", col.axis="grey50", cex.axis=1)

# funBIalign ----
portion_len <- 12*10 #3 years period
min_card  <- 3 

## STEP 1 and 2 -----
# step 1
window_data      <- create_lots(data_smooth, portion_len)

# step 2: compute fMRS-based dissimilarity matrix
D_fmsr  <- create_distance(data_smooth, window_data, portion_len)

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

### CRITERION: adjusted fMSR ----
best_order <- vec_of_scores %>% order()
vec_of_scores_ordered <- vec_of_scores[best_order]
# ordered list_of_recommendations
list_of_recommendations_ordered <- list_of_recommendations[best_order]

#for every recommended motif, compute all its accolites
all_accolites <- lapply(list_of_recommendations_ordered,
                        function(x){
                          lapply(x, get_accolites, window_data, portion_len, FALSE) %>% unlist()
                        })

#Starting from the top, we compare each motif to those with higher rank. If all portions of
# are acolytes to portions of an higher ranking motif, we filter it out; 
# otherwise we retain it.

# we identify the nodes to delete
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

# we delete the recommended nodes and we order the remaining ones
list_of_recommendations_ordered <- list_of_recommendations_ordered[-delete]
vec_of_scores_ordered <- vec_of_scores_ordered[-delete]

### CRITERION: Rank Sum ----
best_order <- (rank(-vec_of_card) + rank(vec_of_scores)) %>% order()

vec_of_scores_ordered <- vec_of_scores[best_order]
# ordered list_of_recommendations
list_of_recommendations_ordered <- list_of_recommendations[best_order]

#for every recommended motif, compute all its accolites
all_accolites <- lapply(list_of_recommendations_ordered,
                        function(x){
                          lapply(x, get_accolites, window_data, portion_len, FALSE) %>% unlist()
                        })

#Starting from the top, we compare each motif to those with higher rank. If all portions of
# are acolytes to portions of an higher ranking motif, we filter it out; 
# otherwise we retain it.
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

# we delete the recommended nodes and we order the remaining ones
list_of_recommendations_ordered <- list_of_recommendations_ordered[-delete]
vec_of_scores_ordered <- vec_of_scores_ordered[-delete]

### CRITERION: Variance 
# Here used just to order and not to delete nodes
motif_var <- lapply(list_of_recommendations_ordered, 
                    function(x){
                      temp <- window_data[x,]
                      temp_mean <- temp %>% colMeans()
                      ((temp - temp_mean)^2 %>% sum())/(nrow(temp))
                    }) %>% unlist()

var_order <- motif_var %>% order(decreasing = TRUE)
list_of_recommendations_ordered <- list_of_recommendations_ordered[var_order]
vec_of_scores_ordered <- vec_of_scores_ordered[var_order]

# Creating some plot ----
## Plot the data and highlight the motif occurrences in red -----
for(q in 1:length(list_of_recommendations_ordered)){
  temp_motif <- list_of_recommendations_ordered[[q]]
  lots_in_motif <- lapply(temp_motif,
                          function(x){
                            strsplit(x, '_') %>% 
                              unlist() %>%
                              as.numeric()}) %>% 
    unlist() %>%
    matrix(ncol=3, byrow=T)
  
  title   <- paste0("Number of instances: ", dim(temp_motif)[2],
                    " - adj fMSR:  ", vec_of_scores_ordered[q] %>% round(3))
  
  matplot(t(data_smooth), col = "grey40", ylab='', xlab='', axes='FALSE',
          lty = 1, type = 'l', main = title)
  for(k in 1:nrow(lots_in_motif)){
    rect(lots_in_motif[k,2], min(data_smooth)-10, lots_in_motif[k,3], max(data_smooth) +10,
         border = alpha("firebrick3", 0.2), col = alpha("firebrick3", 0.2))
    matplot(lots_in_motif[k,2]:lots_in_motif[k,3],
            data_smooth[lots_in_motif[k,1],lots_in_motif[k,2]:lots_in_motif[k,3]], type='l', add=T, col='firebrick3', lwd = 2)
    box(col="grey40",lwd = 2)
    axis(1, at=seq(1, length(data_smooth), by=12), labels = all_time$Time[seq(1, length(data_smooth), by=12)], col="grey40", col.ticks="grey40", col.axis="grey60", cex.axis=1.5)
    axis(2, col="grey40", col.ticks="grey40", col.axis="grey60", cex.axis=1.5)
  }
}

## Plot only the motif -----
for(q in 1:length(list_of_recommendations_ordered)){
  
  temp_motif <- list_of_recommendations_ordered[[q]]
  plot_me <- t(window_data[temp_motif,])
  title   <- paste0("Number of occurrences: ", dim(plot_me)[2],
                    " - adj fMSR:  ", vec_of_scores_ordered[q] %>% round(3))
  
  matplot(plot_me, col = "#648a60", ylab='', xlab='', axes='FALSE',
          lty = 1, type = 'l', lwd = 4, main = title)
  box(col="grey40", lwd = 2)
  axis(1, at=1:portion_len, labels = 1:portion_len, col="grey40", col.ticks="grey40", col.axis="grey60", cex.axis=1.5)
  axis(2, col="grey40", col.ticks="grey40", col.axis="grey60", cex.axis=1.5)
}

