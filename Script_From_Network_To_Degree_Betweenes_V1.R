#Cleaning
rm(list=ls())

if (!requireNamespace("leiden", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
  remotes::install_github("TomKellyGenetics/leiden")
}
if (!requireNamespace("markovchain", quietly = TRUE)) {
  install.packages("markovchain")
}




##Load the RCy3 library if not already loaded
if (!requireNamespace("RCy3", quietly = TRUE)) {
  install.packages("RCy3")
}
library(RCy3)
library(igraph)
library(ggplot2)
library(scales)
library(markovchain)
library(leiden)


##Connect to Cytoscape
#PING
cytoscapePing()

###Degrees
#Get Network Names  (This is runned after the network plot creation)
network_list <- getNetworkList()

#Initialize a list to store igraph objects for each network
igraph_list <- list()

for (net in network_list) {
  ##Get network data as a dataframe
  network_data <- getNetworkSuid(net)
  edge_data <- getTableColumns(table = 'edge', columns = c('source', 'target'), network = network_data)
  node_data <- getTableColumns(table = 'node', columns = c('name', 'type'), network = network_data)
  
  ##Merge node type information with edge data
  edges_with_type <- merge(edge_data, node_data, by.x = 'source', by.y = 'name')
  edges_with_type <- merge(edges_with_type, node_data, by.x = 'target', by.y = 'name', suffixes = c("_source", "_target"))
  
  ##Convert to igraph object
  g <- graph_from_data_frame(edges_with_type, directed = TRUE, vertices = node_data)
  igraph_list[[net]] <- g
}


#empty list
network_degrees <- list()

##Initialize a list to store the overall degree distribution for each network
overall_degree_dist_list <- list()

##Initialize a list to store the degree distribution for each node type within each network
degree_dist_by_type_list <- list()
#degree_dist_list <- list()


##This does not work optimal,

for (net_name in names(igraph_list)) {
  g <- igraph_list[[net_name]]
  
  ##Get degrees for each node
  deg <- degree(g)
  
  ##Calculate overall degree distribution of the network
  overall_deg_dist <- table(deg) / length(deg)
  
  ##Store the overall degree distribution
  overall_degree_dist_list[[net_name]] <- overall_deg_dist
  
  ##Get the 'type' attribute for each node
  types <- V(g)$type
  
  ##Create a dataframe of degrees and types
  degree_type_df <- data.frame(deg, types)
  
  ##Calculate degree distribution for each type
  degree_dist_by_type <- aggregate(deg ~ types, degree_type_df, function(x) table(x) / length(x))
  
  ##Store the result
  degree_dist_by_type_list[[net_name]] <- degree_dist_by_type
}

Whole_list <- list()
Whole_list[["EBV_EBV_combined_network_Network"]] <- igraph_list$EBV_EBV_combined_network_Network
Whole_list[["HPV_HPV16_combined_network_Network"]] <- igraph_list$HPV_HPV16_combined_network_Network

#c(igraph_list$EBV_EBV_combined_network_Network, igraph_list$HPV_HPV16_combined_network_Network)


Combo_overall_degree_dist_list <- list()
Combo_degree_dist_by_type_list <- list()
#Only combo
for (net_name in names(Whole_list)) {
  g <- igraph_list[[net_name]]
  
  ##Get degrees for each node
  deg <- degree(g)
  
  ##Calculate overall degree distribution of the network
  overall_deg_dist <- table(deg) / length(deg)
  
  ##Store the overall degree distribution
  Combo_overall_degree_dist_list[[net_name]] <- overall_deg_dist
  
  ##Get the 'type' attribute for each node
  types <- V(g)$type
  
  ##Create a dataframe of degrees and types
  degree_type_df <- data.frame(deg, types)
  
  ##Calculate degree distribution for each type
  degree_dist_by_type <- aggregate(deg ~ types, degree_type_df, function(x) table(x) / length(x))
  
  ##Store the result
  Combo_degree_dist_by_type_list[[net_name]] <- degree_dist_by_type
}



##Lets plot Combo_overall_degree_dist_list, and caculate the average degree per type



##Initialize a list to store the plot objects
plot_list <- list()

##Plotting each network's degree distribution from Combo_overall_degree_dist_list
for (net_name in names(Combo_overall_degree_dist_list)) {
  ##Assuming each item in the list is a dataframe with columns "Degree" and "Frequency"
  plot_data <- as.data.frame(Combo_overall_degree_dist_list[[net_name]])
  names(plot_data) <- c("Degree", "Frequency")
  
  ##Create the plot and store it in the list
  plot_list[[net_name]] <- ggplot(plot_data, aes(x = Degree, y = Frequency)) +
    geom_bar(stat = "identity") +
    ggtitle(paste("Degree Distribution for Network:", net_name)) +
    xlab("Degree") + ylab("Frequency")
}

##You can now access and print each plot individually
##For example, to print the plot for the first network
print(plot_list[[1]]) #This works


##Initialize a list to store the average degree by type for each network
average_degree_by_type_list <- list()

for (net_name in names(Whole_list)) {
  g <- igraph_list[[net_name]]
  
  ##Get degrees for each node
  deg <- degree(g)
  
  ##Get the 'type' attribute for each node
  types <- V(g)$type
  
  ##Create a dataframe of degrees and types
  degree_type_df <- data.frame(deg, types)
  
  ##Calculate the average degree for each type
  average_degree_by_type <- aggregate(deg ~ types, degree_type_df, mean)
  
  ##Store the result
  average_degree_by_type_list[[net_name]] <- average_degree_by_type
}



##Inspecting the results for a specific networ
print(average_degree_by_type_list$EBV_EBV_combined_network_Network)
print(average_degree_by_type_list$HPV_HPV16_combined_network_Network)

##Assuming igraph_list contains your network data
frequency_table_1 <- list()  ##Frequency of each type attached to "disease"
frequency_table_2 <- list()  ##Frequency of types attached to "human" nodes connected to "disease"

for (net_name in names(Whole_list)) {
  g <- igraph_list[[net_name]]
  
  ##Identify disease nodes
  disease_nodes <- V(g)[V(g)$type == "disease"]
  
  ##Table 1: Frequency of each type attached to "disease"
  neighbors_of_disease <- unlist(lapply(disease_nodes, function(node) neighbors(g, node)$type))
  frequency_table_1[[net_name]] <- table(neighbors_of_disease)
  
  ##Table 2: Frequency of types attached to "human" nodes connected to "disease"
  human_neighbors_data <- unlist(lapply(disease_nodes, function(disease_node) {
    human_neighbors <- V(g)[neighbors(g, disease_node)$type == "human"]
    return(unlist(lapply(human_neighbors, function(human_node) {
      ##Get neighbors of the human node excluding the original disease node
      other_neighbors <- setdiff(neighbors(g, human_node), disease_node)
      return(other_neighbors$type)
    })))
  }))
  frequency_table_2[[net_name]] <- table(human_neighbors_data)
}





########BETWEEENENENS
##Assuming igraph_list contains your network data
Combo_betweenness_by_type_list <- list()
Combo_betweenness_dist_plots <- list()

for (net_name in names(Whole_list)) {
  g <- igraph_list[[net_name]]
  
  ##Isolate the largest component
  components <- components(g)
  largest_component <- which.max(components$csize)
  subgraph_g <<- induced_subgraph(g, which(components$membership == largest_component))
  
  ##Calculate betweenness for each node in the largest component
  btw <- betweenness(subgraph_g, normalized = TRUE)
  
  ##Get the 'type' attribute for each node in the largest component
  types <- V(subgraph_g)$type
  
  ##Create a dataframe of betweenness and types for the largest component
  betweenness_type_df <- data.frame(btw, types)
  
  ##Calculate average betweenness for each type
  avg_btw_by_type <- aggregate(btw ~ types, betweenness_type_df, mean)
  Combo_betweenness_by_type_list[[net_name]] <- avg_btw_by_type
  
  ##Create the betweenness distribution plot
  betweenness_dist_df <- data.frame(Betweenness = btw)
  p <- ggplot(betweenness_dist_df, aes(x = Betweenness)) +
    geom_histogram(binwidth = max(betweenness_dist_df$Betweenness) / 30, fill = "blue", alpha = 0.7) +
    ggtitle(paste("Betweenness Distribution for Largest Component of Network:", net_name)) +
    xlab("Normalized Betweenness") + ylab("Frequency")
  
  Combo_betweenness_dist_plots[[net_name]] <- p
}

##Example to view the results for a specific network
print(Combo_betweenness_by_type_list$EBV_EBV_combined_network_Network)
print(Combo_betweenness_dist_plots$HPV_HPV16_combined_network_Network)





####SHortestPATH

long_comp <- FALSE
if(long_comp){
  print("long version path length")
  ##Assuming igraph_list contains your network data
  path_length_distributions <- list()
  path_length_plots <- list()
  
  ##ATTENTION,THIS TAKE A LOT OF TIME
  for (net_name in names(Whole_list)) {
    g <- igraph_list[[net_name]]
    
    ##Isolate the largest component
    components <- components(g)
    largest_component <- which.max(components$csize)
    subgraph_g <- induced_subgraph(g, which(components$membership == largest_component))
    
    ##Get nodes of type 'viral', 'complex', 'human', and 'disease'
    viral_nodes <- V(subgraph_g)[V(subgraph_g)$type == "viral"]
    complex_nodes <- V(subgraph_g)[V(subgraph_g)$type == "complex"]
    human_nodes <- V(subgraph_g)[V(subgraph_g)$type == "human"]
    disease_nodes <- V(subgraph_g)[V(subgraph_g)$type == "disease"]
    
    ##Function to calculate shortest path lengths from a set of source nodes to target nodes
    get_path_lengths <- function(source_nodes, target_nodes, graph) {
      path_lengths <- numeric()
      for (source in source_nodes) {
        for (target in target_nodes) {
          paths <- shortest_paths(graph, source, target, mode = "out")
          path_lengths <- c(path_lengths, lengths(paths$vpath))
        }
      }
      return(path_lengths)
    }
    
    ##Calculate path length distributions
    viral_to_disease_lengths <- get_path_lengths(viral_nodes, disease_nodes, subgraph_g)
    complex_to_disease_lengths <- get_path_lengths(complex_nodes, disease_nodes, subgraph_g)
    human_to_disease_lengths <- get_path_lengths(human_nodes, disease_nodes, subgraph_g)
    
    ##Store the results in a list
    path_length_distributions[[net_name]] <- list(
      viral_to_disease = viral_to_disease_lengths,
      complex_to_disease = complex_to_disease_lengths,
      human_to_disease = human_to_disease_lengths
    )
    
    ##Store the path length data for plotting
    path_length_data <- data.frame(
      length = c(viral_to_disease_lengths, complex_to_disease_lengths, human_to_disease_lengths),
      type = rep(c("viral_to_disease", "complex_to_disease", "human_to_disease"),
                 c(length(viral_to_disease_lengths), length(complex_to_disease_lengths), length(human_to_disease_lengths)))
    )
    
    ##Plot the path length distribution for each type
    p <- ggplot(path_length_data, aes(x = length, fill = type)) +
      geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
      scale_fill_brewer(palette = "Set1") +
      ggtitle(paste("Path Length Distribution from Types to Disease for Network:", net_name)) +
      xlab("Path Length") + ylab("Frequency") +
      theme_minimal()
    
    ##Store the plot in a list
    path_length_plots[[net_name]] <- p
  }
  long_path_length_distributions <- path_length_distributions
  long_path_length_plots <- path_length_plots
} else{
  print("Reduced version path length")
  #THIS IS REDUCED VERSION, if path longer than 5,it is cutof
  ##Assuming igraph_list contains your network data
  path_length_distributions <- list()
  path_length_plots <- list()
  log_path_length_plots <- list()
  
  for (net_name in names(Whole_list)) {
    g <- igraph_list[[net_name]]
    
    ##Isolate the largest component
    components <- components(g)
    largest_component <- which.max(components$csize)
    subgraph_g <- induced_subgraph(g, which(components$membership == largest_component))
    
    ##Get nodes of type 'viral', 'complex', 'human', and 'disease'
    viral_nodes <- V(subgraph_g)[V(subgraph_g)$type == "viral"]
    complex_nodes <- V(subgraph_g)[V(subgraph_g)$type == "complex"]
    human_nodes <- V(subgraph_g)[V(subgraph_g)$type == "human"]
    disease_nodes <- V(subgraph_g)[V(subgraph_g)$type == "disease"]
    
    ##Calculate shortest path lengths and cap them at 5
    get_capped_path_lengths <- function(source_nodes, target_nodes, graph) {
      path_lengths <- distances(graph, v = source_nodes, to = target_nodes)
      path_lengths[path_lengths > 5] <- 5
      path_lengths <- as.vector(path_lengths)  ##Flatten the matrix to a vector
      return(path_lengths)
    }
    
    ##Calculate path length distributions
    viral_to_disease_lengths <- get_capped_path_lengths(viral_nodes, disease_nodes, subgraph_g)
    complex_to_disease_lengths <- get_capped_path_lengths(complex_nodes, disease_nodes, subgraph_g)
    human_to_disease_lengths <- get_capped_path_lengths(human_nodes, disease_nodes, subgraph_g)
    
    ##Bin the path lengths, treating lengths greater than 5 as '5+'
    bin_path_lengths <- function(lengths) {
      lengths <- as.factor(ifelse(lengths > 5, "5+", lengths))
      return(lengths)
    }
    
    ##Store the binned path lengths for plotting
    path_length_data <- data.frame(
      length = bin_path_lengths(c(viral_to_disease_lengths, complex_to_disease_lengths, human_to_disease_lengths)),
      type = rep(c("viral_to_disease", "complex_to_disease", "human_to_disease"),
                 c(length(viral_to_disease_lengths), length(complex_to_disease_lengths), length(human_to_disease_lengths)))
    )
    
    ##Plot the binned path length distribution for each type
    p <- ggplot(path_length_data, aes(x = length, fill = type)) +
      geom_bar(position = "dodge") +
      scale_fill_brewer(palette = "Set1") +
      ggtitle(paste("Binned Path Length Distribution from Types to Disease for Network:", net_name)) +
      xlab("Path Length") + ylab("Count") +
      theme_minimal()
    
    #Saving the plots
    path_length_plots[[net_name]] <- p
    
    #Log
    ##Create the plot with a logarithmic y-axis
    p_log <- ggplot(path_length_data, aes(x = length, fill = type)) +
      geom_bar(position = "dodge", stat = "count") +
      scale_fill_brewer(palette = "Set1") +
      ggtitle(paste("Binned Path Length Distribution from Types to Disease for Network:", net_name)) +
      xlab("Path Length") + ylab("Count") +
      theme_minimal() +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
      theme(axis.text.y = element_text(angle = 45, hjust = 1))
    
    ##Store the plot in the new list
    log_path_length_plots[[net_name]] <- p_log
  }
}


print(path_length_plots$EBV_EBV_combined_network_Network)
print(log_path_length_plots$EBV_EBV_combined_network_Network)
print(log_path_length_plots$HPV_HPV16_combined_network_Network)


print(path_length_plots$EBV_EBV_combined_network_Network)




#############PLOT FOR RAPPORTS
##Extract data from the log_path_length_plots list
ebv_path_length_data <- log_path_length_plots$EBV_EBV_combined_network_Network$data
hpv_path_length_data <- log_path_length_plots$HPV_HPV16_combined_network_Network$data

##Add 'Network' column to distinguish between the two networks
ebv_path_length_data$Network <- 'EBV'
hpv_path_length_data$Network <- 'HPV'

##Clean the 'type' names by removing underscores and adding a label for 'Type Nodes'
ebv_path_length_data$type <- gsub("_to_", " to ", ebv_path_length_data$type)
hpv_path_length_data$type <- gsub("_to_", " to ", hpv_path_length_data$type)

##Combine the data from both networks
combined_path_length_data <- rbind(ebv_path_length_data, hpv_path_length_data)

##Plot the binned path length distribution for each type with combined networks
combined_log_plot <- ggplot(combined_path_length_data, aes(x = length, fill = type)) +
  geom_bar(stat = "count", position = "dodge") +
  facet_wrap(~ Network, scales = 'free_y') + ##Separate plots for each network
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Path Length to Disease",
       x = "Path Length",
       y = "Count (log)",
       fill = "Type Nodes") + ##Clean legend title
  theme_minimal() +
  scale_y_log10() + ##Log scale for y-axis
  theme(legend.title = element_text(size = 13), ##Custom legend title size
        legend.position = "bottom", ##Move legend to bottom
        axis.text.y = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) ##Rotate y-axis text

##Print the combined log plot
print(combined_log_plot)






##Geting viral nodes for python script |THIS IS NOT NEEEDED FOR ASSIGNMENET
##Assuming 'g' is your igraph graph object and it has a 'type' attribute for nodes
EBV_viral_nodes <- V(Whole_list$EBV_EBV_combined_network_Network)[V(Whole_list$EBV_EBV_combined_network_Network)$type == "viral"]
HPV16_viral_nodes <- V(Whole_list$HPV_HPV16_combined_network_Network)[V(Whole_list$HPV_HPV16_combined_network_Network)$type == "viral"]
##To get the names of these viral nodes
EBV_viral_node_names <- V(Whole_list$EBV_EBV_combined_network_Network)$name[EBV_viral_nodes]
HPV16_viral_node_names <- V(Whole_list$HPV_HPV16_combined_network_Network)$name[HPV16_viral_nodes]


print(EBV_viral_node_names)
print(HPV16_viral_node_names)
