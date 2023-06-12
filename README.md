# funBIalign: a hierachical algorithm for functional motif discovery based on mean squared residue scores
Welcome to the repository containing the functions and the scripts for the paper **funBIalign: a hierachical algorithm for functional motif discovery based on mean squared residue scores**!
This folder contains:
- **funBIalign_functions.R:** a script with all the funBIalign functions;
- **case_study:** a folder with everything needed to reproduce the two case studies presented in the paper;

The **case_study** folder contains:
- **data:** a folder with all the data needed for the food price inflation case study (**food_price_inflation_world.csv** and **food_price_inflation_region.csv**) and the temperature change one (**FAOSTAT_tempchange_region.csv** and **FAOSTAT_tempchange_region_smooth.RData**);
- **regional_food_price_inflation.R:** a script to reproduce the world food price inflation case study;
- **world_food_price_inflation.R:** a script to reproduce the food price inflation case study for the 19 geographical regions;
- **temperature_change.R:** a script to reproduce the temperature change case study for the 19 geographical regions;

# How to use funBIalign?
funBIalign can be applied to any data set having functional observations (one or more) as rows. All the functions required are inside the **funBIalign_functions.R:** script file.
Before using funBIalign to your own dataset, we suggest to check one of the three case study script files, where the method is applied step by step from smoothing the data to plot the resulting motifs.

Here we review the three fundamental steps of the algorithm.
## STEP 1 - Portion creation and alignment
Select the length (as number of data points) `portion_len` and the minimum number of portions `min_card` characterizing the resulting motifs and perform Step 1
```
# Select the length (as number of data points) and the minimum number of portions
portion_len <- 12
min_card <- 3 

# Creates and align all the portions
window_data <- create_lots(mydata, portion_len)
```

## STEP 2 -  Hierarchical clustering based on the fMRS dissimilarity
Compute fMRS-based dissimilarity matrix, perform the agglomerative hierarchical clustering and get the sub-trees.
```
#Compute fMRS-based dissimilarity matrix
D_fmsr <- create_distance(foodinfl_loc, window_data, portion_len)

#Get the sub-trees (tree_s)
minidend <- get_minidend(D_fmsr, window_data)
```

## STEP 3 - Collection of candidate functional motifs
Identify seeds node and their families. Get for every family a recommended node (candidate motif).
```
# Identify seeds and corresponding families
all_paths <- get_path_complete(minidend, window_data, min_card = min_card)
all_paths <- all_paths[!(lapply(all_paths, is.null) %>% unlist())] 

# Get recommended nodes and their info (cardinality and score)
# Collect all recommended nodes (as an array of portion ids)
all_recommended_labels <- lapply(all_paths, function(x){x$recommended_node_labels})
list_of_recommendations <- lapply(rapply(
  all_recommended_labels, enquote, how='unlist'),
  eval)
  
# get recommended node cardinality
vec_of_card <- lapply(list_of_recommendations, length) %>% unlist()

# get recommended node adjusted fMSR
vec_of_scores <- lapply(all_paths, function(x){x$recommended_node_scores}) %>% unlist()
```
For the Step 4 - Post-processing we invite the user to refer to the case studies scripts where commented examples with different criteria are shown

# Simulations
Following the paper **funBIalign: a hierachical algorithm for functional motif discovery based on mean squared residue scores**, we present here the 100 different simulated data set we used to assess the performance of our method. We consider 10 alternative sets of 4 motifs; 2 alternative numbers of occurrences (8 or 10) which, for simplicity, are the same for all motifs; and 4 alternative levels of noise, expressed by σ = 0.1, 0.5, 1 or 2.

The first 80 simulations are stored in the **same_noise** folder that contains motifs sharing the same σ, and we use all σ values in turn.

The last 20 simulations are in the **different_noise** folder that counts simulations having each motif with a different σ.

Every dataset is stored in a `.rds` in the form of a list: the actual data are in the `line` field, while the `motifs` one keeps information about the motifs embedded (curve id, starting data point, and ending data point).
