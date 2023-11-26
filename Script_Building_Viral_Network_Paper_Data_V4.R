#Cleaning
rm(list=ls())

#Libraries
library(readxl)  #For reading Excel files
library(dplyr)   #For data manipulation
library(RCy3)    #For interfacing with Cytoscape
library(igraph)  #For graph data structures and algorithms



####Testing Ground########
#Setting path,
setwd("C:/Users/Shade/Desktop/Master/Network Biology/Project_Results/Viral_Pertubation_Results")
####Individual paths
#Directory paths for saving the networks
ebv_dir <- "C:/Users/Shade/Desktop/Master/Network Biology/Project_Results/Viral_Pertubation_Results/EBV_Networks"
hpv_dir <- "C:/Users/Shade/Desktop/Master/Network Biology/Project_Results/Viral_Pertubation_Results/HPV16_Networks"
#File
excel_file_path <- "C:/Users/Shade/Desktop/Master/Network Biology/Project_Data/Viral_perturbation_study/Dataset_S1.xls" 

#Reading dataset, and filtering on type of interaction
dataset <- read_excel(excel_file_path)


#Extract the first column to identify interaction types
interaction_types <- dataset[[1]]

#Create a character vector of unique interaction types
unique_interaction_types <- unique(interaction_types)

#Print out the unique interaction types for inspection
print(unique_interaction_types)


###########INFO RELATED TO DATASASAT, this is imporant for later on variabel
###OUTPUT EBV:
#[1] "VH-PPI"            "VH-PDI"            "PPI"               "PDI"              
#[5] "MCI_KEGG"          "MCI_FCA|KEGG"      "MCI_FCA"           "MCI_BiGG|KEGG"    
#[9] "MCI_BiGG|FCA|KEGG" "MCI_BiGG"
###OUTPUT HPV16:
#[1] "VH-PPI"                "PPI"                   "PDI"                  
#[4] "MCI_KEGG"              "MCI_FluxCoupling|KEGG" "MCI_BiGG|KEGG"        
#[7] "MCI_BiGG"
###########



#Read all sheet names from the Excel file
sheet_names <- excel_sheets(excel_file_path)

#Function to create edges for a specific type of interaction within the dataset
create_edges <- function(df, interaction_type, debug = FALSE) {
  #Print the current interaction type if debugging is enabled
  if(debug) {
    print(paste("Creating edges for", interaction_type))
    print("Original column names:")
    print(colnames(df))
  }
  
  #Sanitize the column names to ensure consistency
  colnames(df) <- gsub("\n", " ", colnames(df))  #Replace newline characters with spaces
  colnames(df) <- gsub(" +", " ", colnames(df))  #Replace multiple spaces with a single space
  colnames(df) <- trimws(colnames(df))           #Trim leading and trailing whitespace
  
  #Rename columns to match expected names
  df <- rename(df, `Type of interaction` = colnames(df)[1],
               `Viral protein` = colnames(df)[2],
               `Intermediate human protein` = colnames(df)[3],
               `Human protein` = colnames(df)[4],
               `Disease Id` = colnames(df)[5],
               `ICD-9` = colnames(df)[6],
               `OMIM Description` = colnames(df)[7])
  
  ####OLD CODE, CAN BE IGNORED####
  #if(debug) print("Updated column names to standard names.")
  ##Check if expected columns exist in the dataframe, if not, stop the function and print an error message
  #expected_cols <- c("Type of interaction", "Viral protein", "Intermediate human protein", "Human protein", "Disease Id")
  #missing_cols <- setdiff(expected_cols, colnames(df))
  #if(length(missing_cols) > 0) {
  #  stop("The following expected columns are missing in the dataframe: ", paste(missing_cols, collapse=", "), call. = FALSE)
  #}
  
  #Filter the dataframe for the specified interaction type
  df_filtered <- df %>% filter(`Type of interaction` == interaction_type)
  
  #Create edges considering the presence of an intermediate protein and disease association
  edges <- df_filtered %>%
    select(
      Source = `Viral protein`, 
      Intermediate = `Intermediate human protein`, 
      Target = `Human protein`,
      DiseaseId = `Disease Id`
    ) %>%
    mutate(
      #Create a flag for intermediate protein presence
      IntermediateFlag = ifelse(Intermediate == "-", FALSE, TRUE),
      #Create a flag for disease association presence
      DiseaseAssociated = ifelse(DiseaseId == "-", FALSE, TRUE)
    ) %>%
    filter(Source != "-", Target != "-") %>%
    distinct() %>%
    mutate(
      #Modify the Source to include the intermediate protein if present
      Source = ifelse(IntermediateFlag, paste0(Source, "-", Intermediate), Source),
      #Modify the Target to indicate disease association if present
      Target = ifelse(DiseaseAssociated, paste0(Target, "*"), Target)
    ) %>%
    select(Source, Target) %>%
    distinct() #Removing double edges, and such
  
  #Return the edge list dataframe
  return(edges)
}


#Function to build networks for each interaction type and combine them
build_networks <- function(sheet_index, interaction_types, debug = FALSE, Virus_Type) {
  #Read data for the specified sheet in the Excel file
  df <- read_excel(excel_file_path, sheet = sheet_names[sheet_index])
  if(debug) print(paste("Processing sheet:", sheet_names[sheet_index]))
  
  #Initialize an empty list to store networks keyed by interaction type
  networks <- list()
  
  #Iterate over each interaction type and create a network
  for(interaction_type in interaction_types) {
    if(interaction_type == "MCI_FluxCoupling|KEGG" & Virus_Type == "EBV"){
      next
    }
  
    if(debug) print(paste("Creating network for interaction type:", interaction_type))
    edges <- create_edges(df, interaction_type, debug)
    
    #Check if there are any edges to create a network
    if(nrow(edges) > 0) {
      network <- graph_from_data_frame(edges, directed = TRUE)
      
      #Ensure that all vertices have names
      V(network)$name <- ifelse(is.na(V(network)$name), "", V(network)$name)
      
      networks[[interaction_type]] <- network
    } else {
      if(debug) print(paste("No edges for interaction type:", interaction_type))
    }
  }
  
  #Check if we have at least one network to combine
  if(length(networks) > 0) {
    #Get the names of the first network as base for the combined network
    combined_network <- networks[[1]]
    
    #Combine the networks
    for(i in 2:length(networks)) {
      combined_network <- graph.union(combined_network, networks[[i]], byname = TRUE)
    }
  } else {
    combined_network <- make_empty_graph(directed = TRUE)
  }
  
  return(list("individual" = networks, "combined" = combined_network))
}

#Define the interaction types of interest
interaction_types_of_interest <- c("VH-PPI","VH-PDI" ,"PPI", "PDI", "MCI_KEGG","MCI_FCA|KEGG", "MCI_FCA","MCI_BiGG|KEGG","MCI_BiGG|FCA|KEGG", "MCI_BiGG","MCI_FluxCoupling|KEGG") #Fixed, this


#Build networks for each virus
ebv_networks <- build_networks(1, interaction_types_of_interest, debug = TRUE, "EBV")
hpv_networks <- build_networks(2, interaction_types_of_interest, debug = TRUE,"HPV16")




###Making save functions, to be sure
#Function to save individual networks to files
save_networks <- function(network_list, dir_path) {
  #Create directory if it doesn't exist
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  #Iterate over the list of networks and save each network
  for (interaction_type in names(network_list$individual)) {
    network <- network_list$individual[[interaction_type]]
    #Replace invalid characters in file names
    sanitized_interaction_type <- gsub("\\|", "_", interaction_type)
    file_name <- paste0(dir_path, "/", sanitized_interaction_type, "_network.rds")
    saveRDS(network, file = file_name)
  }
}


#Function to save the combined network to a file
save_combined_network <- function(network, dir_path, virus_name) {
  #Define the file name for the combined network
  file_name <- paste0(dir_path, "/", virus_name, "_combined_network.rds")
  
  #Save the combined network
  saveRDS(network, file = file_name)
}


###SAVING Networks
#Save the individual networks and the combined network for EBV
save_networks(ebv_networks, ebv_dir)
save_combined_network(ebv_networks$combined, ebv_dir, "EBV")

#Save the individual networks and the combined network for HPV
save_networks(hpv_networks, hpv_dir)
save_combined_network(hpv_networks$combined, hpv_dir, "HPV16")




#Print summaries of the networks
print(summary(ebv_networks$combined))
print(summary(hpv_networks$combined))

