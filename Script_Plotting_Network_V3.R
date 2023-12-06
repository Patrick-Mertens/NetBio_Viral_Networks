#Cleaning
rm(list=ls())

#Libaries
library(RCy3)
library(igraph)

#Paths
#Setting path
setwd("C:/Users/Shade/Desktop/Master/Network Biology/Project_Results/Viral_Pertubation_Results")
#Individual paths
#Directory paths for saving the networks
ebv_dir <- "C:/Users/Shade/Desktop/Master/Network Biology/Project_Results/Viral_Pertubation_Results/EBV_Networks_V2"
hpv_dir <- "C:/Users/Shade/Desktop/Master/Network Biology/Project_Results/Viral_Pertubation_Results/HPV16_Networks_V2"
#File
excel_file_path <- "C:/Users/Shade/Desktop/Master/Network Biology/Project_Data/Viral_perturbation_study/Dataset_S1.xls" 

##Functions
create_network_from_rds <- function(rds_path, style_name, df_protein_types,virus_type) {
  #Load the network from the RDS file
  network <- readRDS(rds_path)
  
  #Count the number of edges before simplification
  original_edge_count <- ecount(network)
  
  #Simplify the network by removing duplicate edges and self-loops
  network <- simplify(network, remove.multiple = TRUE, remove.loops = TRUE)
  
  #Count the number of edges after simplification
  simplified_edge_count <- ecount(network)
  
  
  ##Controle
  #Print out the original and simplified edge counts
  cat("Network from", rds_path, "originally had", original_edge_count, "edges.\n")
  cat("After simplification, it has", simplified_edge_count, "edges.\n")
  cat("Removed", original_edge_count - simplified_edge_count, "duplicate/loop edges.\n\n")
  ##Output: NO DOUBLE EDGES!!!!
  
  #Set the 'type' attribute based on protein lists
  V(network)$type <- sapply(V(network)$name, function(n) {
    if (grepl("^[^-]+--[^-]+$", n)) { #Regex to ensure "--" is between characters
      return("complex") #Nodes with a dash between characters are intermediates forming a complex
    } else if (grepl("\\*\\*$", n)) {
      return("disease")
    #Nodes with an asterisk are associated with a disease
    } else if (n %in% df_protein_types$viral_proteins) {
      return("viral")
    } else {
      return("human") #Default to human if none of the above conditions are met
    }
  })
  
  #Extract the interaction type from the file name
  interaction_type <- gsub(".rds$", "", basename(rds_path))
  title <- paste(virus_type, interaction_type, "Network", sep = "_")
  
  #Create network in Cytoscape with the updated title
  network_suid <- createNetworkFromIgraph(network, title = title)
  
  #Create network in Cytoscape with the title as the basename of the RDS file
  #network_suid <- createNetworkFromIgraph(network, title = basename(rds_path))
  
  #Set the visual style for the network
  setVisualStyle(style_name, network = network_suid)
  
  #Apply the layout
  layoutNetwork("force-directed", network = network_suid)
  
  #Return the simplified network for further manipulation
  return(network)
}


#Function to apply the visual style and layout to a network in Cytoscape
apply_style_and_layout <- function(network_suid) {
  RCy3::setVisualStyle(style.name = style_name, network = network_suid)
  RCy3::layoutNetwork("force-directed", network = network_suid)
}


#Function to check if Cytoscape is running and clear any existing sessions
check_and_clear_cytoscape <- function() {
  #Try pinging Cytoscape and catch any potential errors
  ping_result <- tryCatch({
    cytoscapePing()
    TRUE  #If cytoscapePing() is successful, return TRUE
  }, error = function(e) {
    FALSE  #If there is an error, return FALSE
  })
  
  if (ping_result) {
    #Clear the current session if Cytoscape is running
    deleteAllNetworks()
  } else {
    #Start Cytoscape if not running or if ping failed
    #print("This doesnt work") #This works, for some reason
    cytoscapeOpen()
    Sys.sleep(360) #To be safe
    deleteAllNetworks() #Remvogin, Networks if there are still in them
  }
}  




######Execution groind####################
#Loading data
df_ebv <- read_excel(excel_file_path, sheet = 1)
df_HPV16 <- read_excel(excel_file_path, sheet = 2)

##Colanmes of excel is weird, resetting it to something more normal
#df_ebv reseting
df_ebv <- rename(df_ebv, `Type of interaction` = colnames(df_ebv)[1],
             `Viral protein` = colnames(df_ebv)[2],
             `Intermediate human protein` = colnames(df_ebv)[3],
             `Human protein` = colnames(df_ebv)[4],
             `Disease Id` = colnames(df_ebv)[5],
             `ICD-9` = colnames(df_ebv)[6],
             `OMIM Description` = colnames(df_ebv)[7])
#df_HPV16 Reseting
df_HPV16 <- rename(df_HPV16, `Type of interaction` = colnames(df_HPV16)[1],
             `Viral protein` = colnames(df_HPV16)[2],
             `Intermediate human protein` = colnames(df_HPV16)[3],
             `Human protein` = colnames(df_HPV16)[4],
             `Disease Id` = colnames(df_HPV16)[5],
             `ICD-9` = colnames(df_HPV16)[6],
             `OMIM Description` = colnames(df_HPV16)[7])



#protein list
df_protein_types_ebv <- list(
  viral_proteins = unique(df_ebv$`Viral protein`),
  intermediate_proteins = unique(df_ebv$`Intermediate human protein` != "-"),
  disease_associated_proteins = unique(df_ebv$`Human protein`[df_ebv$`Disease Id` != "-"])
)

#Define protein lists for HPV16 based on your dataset
df_protein_types_hpv <- list(
  viral_proteins = unique(df_HPV16$`Viral protein`),
  intermediate_proteins = unique(df_HPV16$`Intermediate human protein`[df_HPV16$`Intermediate human protein` != "-"]),
  disease_associated_proteins = unique(df_HPV16$`Human protein`[df_HPV16$`Disease Id` != "-"])
)




#Define the colors for each type of node
node_colors <- c("viral"="#FF0000", "complex"="#FFA500", "human"="#00FF00", "disease"="#0000FF", "default"="#CCCCCC")

style_name <- 'Viral_Networks'
#Check if the style already exists
existing_styles <- RCy3::getVisualStyleNames()
if (!(style_name %in% existing_styles)) {
  #If the style doesn't exist, create it
  RCy3::createVisualStyle(style.name = style_name)
}

#Now set the default node color
RCy3::setNodeColorDefault(new.color = node_colors["default"], style.name = style_name)

#Map the 'type' attribute to node colors

#RCy3::setNodeColorMapping(table.column = "type", 
#                          table.column.values = c("viral", "complex", "human", "disease"), 
#                          colors = c(node_colors["viral"], node_colors["complex"], node_colors["human"], node_colors["disease"]), 
#                          style.name = style_name)



##COLOURING, #THis works some times, really weird
colouring <- FALSE
if(colouring){
  networks_list<- RCy3::getNetworkList()
  Whole_networks <- c("EBV_EBV_combined_network_Network","HPV_HPV16_combined_network_Network")
  
  for (net in Whole_networks) {
    network_name <- Whole_networks[net]
    ID <- RCy3::getNetworkSuid(title = network_name)
    RCy3::setNodeColorMapping("type",
                              c("viral", "complex", "human", "disease"),
                              c(node_colors["viral"], node_colors["complex"], node_colors["human"], node_colors["disease"]),
                              style.name = style_name, mapping.type = "d",network = ID)
    
  }
  # RCy3::setNodeColorMapping("type",
  #                     c("viral", "complex", "human", "disease"),
  #                     c(node_colors["viral"], node_colors["complex"], node_colors["human"], node_colors["disease"]),
  #                     style.name = style_name, mapping.type = "d")
}
#Get the list of RDS files for EBV
ebv_files <- list.files(ebv_dir, pattern = "\\.rds$", full.names = TRUE)

#Get the list of RDS files for HPV16
hpv_files <- list.files(hpv_dir, pattern = "\\.rds$", full.names = TRUE)



##Cytoscape
check_and_clear_cytoscape()

#Creating
#Pass the df_protein_types list to the function for EBV
for (file_path in ebv_files) {
  create_network_from_rds(file_path, style_name, df_protein_types_ebv, "EBV")
}

#Pass the df_protein_types list to the function for HPV16
for (file_path in hpv_files) {
  create_network_from_rds(file_path, style_name, df_protein_types_hpv,"HPV")
}