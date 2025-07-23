####### Phylogenetic Clustering and Monophyly Analysis #######
# Author: Shinnam Yoo
# Description: This script processes phylogenetic tree clusters to solve monophyly issues iteratively.

# Load necessary libraries ---------------------------------------------------
library(data.table)
library(stringr)
library(tibble)
library(dplyr)
library(MonoPhy)
library(tidyr)
library(ape)
library(treeio)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(viridis)
library(extrafont)
library(ggstar)
library(readr)
library(treestats)
library(gridExtra)
library(pastecs)
library(phangorn)

# User-defined Parameters --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 1 && length(args) >= idx + 1) {
    return(args[idx + 1])
  }
  return(default)
}

input <- get_arg("--in")
path <- get_arg("--path")
db_fasta <- get_arg("--DB")

setwd(paste0(path,"/cluster"))

# 1. Solve monophyly ------------------------------------------------------

## 1) Load Cluster and Tree Data --------------------------------------------------

# Get list of files matching the pattern
files <- list.files(pattern = "raxml.expanded.nwk.rooted_") # parameter

# Extract numeric identifiers from filenames, sort, and convert to character
num <- gsub(".*?([0-9]+).*", "\\1", files) %>%  # parameter
  as.numeric() %>% 
  sort() %>% 
  as.character()

# Create empty lists to store results
tree_ls_ini <- list()    # For storing tree objects
sol_ls_ini <- list()     # For storing monophyly solutions
clst_ls_ini <- list()    # For storing cluster data


## 2) Solve monophyly for each subtree ----------------------------------------

for (i in num) {
  
  # Import cluster data
  t <- fread(paste0('cluster_',i,'.tsv'), sep = '\t', header = T)
  t$ClusterNumber <- t$ClusterNumber %>% as.character()  # Convert ClusterNumber to character
  t$cluster <- paste0(t$ClusterNumber, '_', t$SequenceName)  # Create combined cluster identifier
  
  # Store the cluster data in the list
  clst_ls_ini[[i]] <- t
  
  # Import corresponding subtree (tree) data
  tr <- read.tree(paste0('raxml.expanded.nwk.rooted_',i,'.nwk'))
  
  # Rename taxa in the tree based on the cluster data
  tr <- rename_taxa(tr, clst_ls_ini[[i]], SequenceName, cluster)
  
  # Collapse polytomies
  tr <- di2multi(tr, tol = 1e-05) # parameter
  
  # Store the processed tree in the list
  tree_ls_ini[[i]] <- tr
  
  # Assess monophyly for the tree and store the result
  sol <- AssessMonophyly(tr, outlierlevel = 0.99)  # parameter # if outlierlevel > 0.5, this iterates until only monophyletic groups are left
  
  # Store the monophyly solution in the list
  sol_ls_ini[[i]] <- sol
}

# Combine all data frames from sol_ls_ini and distinguish by list index
TreeCluster_stats <- do.call(rbind, 
                             lapply(names(sol_ls_ini), function(i) {
                               # Extract summary of monophyly for each subtree
                               s <- GetSummaryMonophyly(sol_ls_ini[[i]])$Genera
                               # Convert the summary to a data frame
                               df <- as.data.frame(do.call(cbind, unclass(s)))
                               # Add Category and Subtree information
                               df$Category <- attr(s, "row.names")
                               df$Subtree <- i
                               return(df)
                             }))

# Remove row names
rownames(TreeCluster_stats) <- NULL

# Reorder columns
TreeCluster_stats <- TreeCluster_stats %>% select("Subtree", "Category", "Taxa", "Tips")

# Write the combined data frame to a CSV file
write.csv(TreeCluster_stats, paste0(path,'/taxonomy/TreeCluster_stats.csv'), quote = F, row.names = F)

## 3) Iteratively solve monophyly for each subtree ----------------------------

# Clone initial lists
tree_ls <- tree_ls_ini    # tree objects
sol_ls <- sol_ls_ini     # monophyly solutions
clst_ls <- clst_ls_ini    # cluster data

iteration <- 1

# Repeat until no outliers are left
repeat {
  
  start_time <- Sys.time()  # Record the start time
  cat(sprintf("Iteration %d is running...\n", iteration))
  
  # Solve outliers
  for (i in num) {
    tryCatch({
      
      # Get list of outlier tips from the previous solution
      ls <- GetOutlierTips(sol_ls[[i]])$Genera 
      
      # Get vector of outlier tips (singletons are marked with -1)
      out <- ls[!grepl("-1", names(ls))] %>% unlist() %>% as.vector()
      
      # Skip iteration if no outliers are left
      if (is.null(out) || length(out) == 0) {
        next
      }
      
      # Assign new groups to outliers
      clst_ls[[i]][clst_ls[[i]]$cluster %in% out,]$ClusterNumber <- paste0(str_replace_all(clst_ls[[i]][clst_ls[[i]]$cluster %in% out,]$ClusterNumber, "\\..*",""), ".", iteration + 1)
      clst_ls[[i]]$cluster <- paste0(clst_ls[[i]]$ClusterNumber, '_', clst_ls[[i]]$SequenceName)
      
      # Solve monophyly
      
      # Import subtree and append updated cluster information
      tr <- read.tree(paste0('raxml.expanded.nwk.rooted_',i,'.nwk'))
      tr <- rename_taxa(tr, clst_ls[[i]], SequenceName, cluster)
      
      # Collapse polytomies
      tr <- di2multi(tr, tol = 1e-05)
      
      # Store updated tree
      tree_ls[[i]] <- tr
      
      # Solve monophyly for the updated tree
      sol <- AssessMonophyly(tr, outlierlevel = 0.99) 
      
      # Store the solution
      sol_ls[[i]] <- sol
      
    }, error = function(e) {
      # Handle exceptions and print error messages
      print(paste("An error occurred for", i, ":"))
      print(e)
    })
  }
  
  # Solve results
  results <- sapply(sol_ls, function(x) {
    o <- GetOutlierTips(x)$Genera[!grepl("-1", names(GetOutlierTips(x)$Genera))]
    length(o) == 0
  }) %>% unlist()  
  
  # If no outliers are left, break the loop
  if (all(results)) {
    end_time <- Sys.time()  # Record end time
    cat(sprintf("Iteration %d completed in %f seconds\n", iteration, as.numeric(end_time - start_time, units = "secs")))
    break
  }
  
  end_time <- Sys.time()  # Record end time
  cat(sprintf("Iteration %d completed in %f seconds\n", iteration, as.numeric(end_time - start_time, units = "secs")))
  
  # Increment iteration
  iteration <- iteration + 1
  
}

# Combine all data frames from sol_ls and distinguish by list index

Monophyly_stats <- do.call(rbind, 
                           lapply(names(sol_ls), function(i) {
                             s <- GetSummaryMonophyly(sol_ls[[i]])$Genera
                             df <- as.data.frame(do.call(cbind, unclass(s)))
                             df$Category <- attr(s, "row.names")
                             df$Subtree <- i
                             return(df)
                           })
)

# Remove row names
rownames(Monophyly_stats) <- NULL

# Reorder columns
Monophyly_stats <- Monophyly_stats %>% select("Subtree", "Category", "Taxa", "Tips")

# Write the combined data frame to a CSV file
write.csv(Monophyly_stats, paste0(path,'/taxonomy/Monophyly_stats.csv'), quote = F, row.names = F)

# Get summary
Monophyly_summary <- do.call(rbind, 
                             lapply(names(sol_ls), function(i) {
                               s <- GetResultMonophyly(sol_ls[[i]])$Genera
                               s <- s %>% rownames_to_column(var = "Subcluster")
                               df <- as.data.frame(do.call(cbind, unclass(s)))
                               df$Subtree <- i
                               return(df)
                             })
)

# Remove row names
rownames(Monophyly_summary) <- NULL

# Write the combined data frame to a CSV file
write.csv(Monophyly_summary, paste0(path,'/taxonomy/Monophyly_summary.csv'), quote = F, row.names = F)

# Get final cluster result
pt <- data.frame()

for (i in names(sol_ls)) {
  
  df <- sol_ls[[i]]$Genera$TipStates
  
  # Add subtree column
  df$Subtree <- i 
  
  # Bind the current dataframe to the result dataframe
  pt <- rbind(pt, df)
  
}

# Recover original tip name 
pt$Tip <- str_replace(pt$Tip, paste0(pt$Taxon, "_"), "")

# Replace -1 to Monotypic
pt[pt$Taxon == -1, ]$Taxon <- "Monotypic"

# Make cluster column
pt$Cluster <- paste0(pt$Subtree, ".", pt$Taxon)

# Identify unique cluster values
unique_taxa <- names(which(table(pt$Cluster) == 1))

# Replace unique cluster values with "Monotypic"
pt <- pt %>%
  mutate(Cluster = ifelse(Cluster %in% unique_taxa, paste0(Subtree,".Monotypic"), Cluster))

# Synchronize Taxon to Cluster
pt[grepl("Monotypic",pt$Cluster),]$Taxon <- "Monotypic"

# Add Source column
pt$Source <- NA

headers <- readLines(db_fasta)
vec <- grep("^>", headers, value = TRUE)
vec <- sub("^>", "", vec)

pt$Source <- ifelse(pt$Tip %in% vec, "DB", pt$Source)

# Tips that are not DB are query
pt[is.na(pt$Source), ]$Source <- "query"

# Remove row names
rownames(pt) <- NULL

write.csv(pt, paste0(path,'/taxonomy/phylotype.csv'), quote = F, row.names = F)


# 2. LCA-based taxonomic assignment ------------------------------------------

## 1) Define functions -----------------------------------------------------

# Define a function to extract common parts from each position
common_strings <- function(strings) {
  if (length(strings) == 0) return(character(0))
  Reduce(function(x, y) {
    common <- which(x == y)
    if (length(common) == 0) return(character(0))
    x[common]
  }, strings)
}

# Define a function to extract common substring parts
# input(taxonomy) should be a vector whose element is ; separated taxonomy
consensus_tax <- function(taxonomy) {
  # Remove NA values from the vector
  taxonomy_clean <- na.omit(taxonomy)
  
  # Split each element into a list based on semicolons
  split_taxonomy <- str_split(taxonomy_clean, ";")
  
  # Extract common parts from each position
  common_parts <- common_strings(split_taxonomy)
  
  # Join common parts with semicolons and return
  paste(common_parts, collapse = ";")
}

# Define the iterative_lca function to process LCA taxonomy
# This function may overclassify long-branched query tips, but this behavior might be relevant in some contexts
iterative_lca <- function(tips, tree_table_with_taxonomy, tree) {
  
  # Get the LCA (Lowest Common Ancestor) node of the clade
  lca <- ifelse(length(tips) == 1, 
                tips, 
                getMRCA(tree %>% as.phylo(), tips))
  
  # Get the cluster associated with the tips, excluding monotypic taxon
  cluster <- tree_table_with_taxonomy %>% 
    filter(node %in% tips) %>% 
    filter(Taxon != 'Monotypic') %>% 
    pull(Cluster) %>% 
    unique()
  
  first_time <- TRUE  # Flag to track the first iteration
  
  repeat {
    # Get all descendant tips of the MRCA (Lowest Common Ancestor)
    des <- getDescendants(tree, lca)
    
    # Get information of the descendant tips
    des_df <- tree_table_with_taxonomy %>% 
      filter(node %in% des & isTip == TRUE)
    
    # Exclude long-branched DB tips that could mislead taxonomic assignment
    des_df <- des_df %>% filter(Long == "no")
    
    # Count the number of DBs in the descendants
    db_num <- grep("DB", des_df$Source) %>% length()
    
    # 1. If the number of DBs under MRCA is >= 1
    if (db_num >= 1) {
      
      # 1.1 If more than one DB is under MRCA (Only executed during the first iteration)
      if (first_time && db_num > 1) {
        
        # Identify DB tips and their parent nodes
        db_tips <- des_df %>% 
          filter(Source == "DB") %>% 
          filter(Cluster == cluster) %>% 
          pull(node)
        
        db_parents <- Ancestors(tree, db_tips, type = 'parent') %>% unique()
        
        # Rearrange DB parents from the smallest clade to avoid assigning large clades to ambiguous taxonomic levels
        dt2 <- data.table()
        
        for (db_parent in db_parents) {
          dt1 <- data.table()
          dt1$des <- getDescendants(tree, db_parent)
          dt1$parent <- db_parent
          dt2 <- rbind(dt2, dt1)
        }
        
        if (nrow(dt2) > 0) {
          db_parents <- dt2 %>% 
            group_by(parent) %>% 
            summarise(count = n()) %>% 
            arrange(count) %>% 
            pull(parent)
        }
        
        # Traverse DB parents and gather query taxonomy
        for (db_parent in db_parents) {
          
          prev_db_des_df <- NULL
          db_des_df <- NULL
          
          repeat {
            
            # Get all descendant tips of the DB parent
            db_des <- getDescendants(tree, db_parent)
            
            # Get information of the descendants
            db_des_df <- tree_table_with_taxonomy %>% 
              filter(node %in% db_des & isTip == TRUE)
            
            # Exclude long-branched DB tips
            db_des_df <- db_des_df %>% filter(Long == "no")
            
            # Count the number of DBs among the descendants
            db_num <- grep("DB", db_des_df$Source) %>% length()
            
            # If no other DB is included, continue walking up the ancestor tree
            if (db_num == 1) {
              db_parent <- Ancestors(tree, db_parent, type = 'parent') %>% unique()
              prev_db_des_df <- db_des_df  # Update the previous descendant dataframe
            } else { 
              # If DBs are included, assign consensus taxonomy to the query nodes
              if (is.null(prev_db_des_df)) {
                prev_db_des_df <- db_des_df
              }
              
              db_tax <- prev_db_des_df[prev_db_des_df$Source == "DB", ]$taxonomy %>% na.omit() %>% as.vector()
              con_tax <- consensus_tax(db_tax)
              
              # Maintain queries only within the target clade
              prev_db_des_df <- prev_db_des_df %>% filter(Source == "query")
              prev_db_des_df[prev_db_des_df$node %in% tips, ]$taxonomy <- con_tax
              prev_db_des_df <- prev_db_des_df %>% filter(!is.na(taxonomy))  # Remove incomplete taxonomies
              
              if (nrow(result) == 0) {
                result <<- prev_db_des_df
              } else {
                new_rows <- anti_join(prev_db_des_df, result, by = "label")
                result <<- rbind(result, new_rows)  # Add only non-duplicate rows
              }
              
              break  # Exit the loop
            }
          }
        }
        
        # Handle remaining queries and DBs not included in the result
        remain_df <- anti_join(des_df, result, by = "label")
        db_tax <- remain_df[remain_df$Source == "DB", ]$taxonomy %>% na.omit() %>% as.vector()
        con_tax <- consensus_tax(db_tax)
        
        remain_df <- remain_df %>% filter(Source == "query")
        remain_df[remain_df$node %in% tips, ]$taxonomy <- con_tax
        remain_df <- remain_df %>% filter(!is.na(taxonomy))
        
        result <<- rbind(result, remain_df)  # Update the global result variable
        
      } else {
        # 1.2 If DB number under MRCA is 1 during the first iteration
        # Get consensus taxonomy of the DBs
        db_info <- des_df[des_df$Source == "DB", ]
        db_tax <- db_info$taxonomy %>% na.omit() %>% as.vector()
        con_tax <- consensus_tax(db_tax)
        
        # Remove species level if necessary (depending on iteration and DB tip cluster)
        if (str_count(con_tax, ";") == 6 && (!first_time || (first_time && db_info$Cluster != cluster))) {
          con_tax <- str_split(con_tax, ";")[[1]]
          con_tax <- paste(con_tax[-length(con_tax)], collapse = ";")
        }
        
        des_df <- des_df %>% filter(Source == "query")  # Maintain only query tips
        des_df[des_df$node %in% tips, ]$taxonomy <- con_tax
        des_df <- des_df %>% filter(!is.na(taxonomy))  # Remove incomplete taxonomies
        
        result <<- rbind(result, des_df)  # Update the global result variable
      }
      break  # Exit the loop
      
    } else {
      # 2. If DB number under MRCA is less than 1 (i.e., 0), move up to the parent
      lca <- tree_table_with_taxonomy[tree_table_with_taxonomy$node == lca,]$parent
    }
    
    first_time <- F  # Set flag to F after the first iteration
  }
}

# Define the get_lca function to assign taxonomy
get_lca <- function(clade_vec, tree_table_with_taxonomy, tree) {
  
  for (clade in clade_vec) {
    
    if (clade == "Monotypic") {
      
      tips <- tree_table_with_taxonomy %>% 
        filter(Taxon == clade & Source == 'query') %>%  # Only process query tips for Monotypic taxa
        pull(node)
      
      for (tip in tips) {
        iterative_lca(tip, tree_table_with_taxonomy, tree)  # Call iterative_lca function
      }
    } else {
      
      tips <- tree_table_with_taxonomy %>% 
        filter(Taxon == clade) %>%  # Handle non-Monotypic clades
        pull(node)
      
      iterative_lca(tips, tree_table_with_taxonomy, tree)  # Call iterative_lca function
    }
  }
  
  return(result)  # Return the final result
}


## 2) Load phylotype & long-branch data ----------------------------------------

# Load phylotype results
pt <- fread(paste0(path, '/taxonomy/phylotype.csv'), sep = ',')
names(pt)[1] <- "label"

# Load long-branched tips
tips_long <- fread(paste0('long-branched_tips.tsv'), sep = '\t', header = F)

# Filter long-branched DB tips only
tips_long <- tips_long[V2 %in% vec]

names(tips_long) <- c('Subtree', 'label')
tips_long$Long <- 'yes'

## 3) Load taxonomy ------------------------------------------------------------

ref_tax <- fread(paste0(input, '/taxonomy.txt'), sep = '\t', header = F, quote = F) %>% 
  unique() %>% 
  tibble() # parameter

names(ref_tax) <- c("label", "taxonomy")
ref_tax$label <- as.character(ref_tax$label)

# Get list of files matching the pattern
files <- list.files(pattern = "raxml.expanded.nwk.rooted_")

# Extract numeric identifiers from filenames, sort, and convert to character
num <- as.numeric(gsub(".*?([0-9]+).*", "\\1", files)) %>% sort() %>% as.character()

# Create empty lists and a data frame to store results
raw_tree <- list()
raw_plot <- list()
shrunk_tree <- list()
shrunk_plot <- list()
result <- data.frame()
query_tax_df <- data.frame()
lump_df <- data.frame()
query_tax_lump_df <- data.frame()

## 4) Taxonomic assignment --------------------------------

for (i in num) {
  tryCatch({
    start_time <- Sys.time()  # Record start time
    
    ### Step 1: Load Tree ---------------------------------------
    # Load the phylogenetic tree generated by RAxML
    tr <- read.tree(paste0('raxml.expanded.nwk.rooted_', i, '.nwk'))
    
    # Collapse polytomies
    tr <- di2multi(tr, tol = 1e-05)
    raw_tree[[i]] <- tr  # Store the original tree
    
    # Visualize the tree
    p <- ggtree(tr, layout = "rectangular", size = 0.1)
    
    ## Step 2: Retrieve phylotype & long-branch annotations --------------------
    # Filter metadata for the current subtree
    pt_sub <- pt %>% filter(Subtree == i) # subset phylotype data
    long_sub <- tips_long %>% filter(Subtree == i) %>% select(label, Long) # subset long-branch data
    
    # Assign phylotype
    dat <- left_join(p$data, pt_sub, by = "label") %>% 
      left_join(long_sub, by = "label")
    
    # Assign "no" to tips without long-branch designation
    dat[is.na(dat$Long) & dat$isTip == TRUE,]$Long <- "no"
    p$data <- dat  # Apply updated data
    raw_plot[[i]] <- p  # Store the original plot
    
    # Identify nodes labeled as long-branched
    node_num <- dat[dat$Long == "yes", ]$node
    
    ## Step 3: Prune tree ------------------------------------------------------
    
    # Remove long-branched DB tips from the tree
    tr2 <- drop.tip(tr, node_num[!is.na(node_num)])
    shrunk_tree[[i]] <- tr2  # Store the pruned tree
    
    # Visualize the pruned tree
    # dat2 is a tree table object w/o long-branched DB tips. 
    # dat2 must be generated from tree object (tr2, p2). if not, MRCA calculation may be misled.
    p2 <- ggtree(tr2, layout = "rectangular", size = 0.1)
    dat2 <- left_join(p2$data, pt_sub, by = "label") %>% 
      left_join(long_sub, by = "label")
    
    # Assign "no" to tips without long-branch designation
    dat2[is.na(dat2$Long) & dat2$isTip == TRUE,]$Long <- "no"
    p2$data <- dat2  # Apply updated data
    shrunk_plot[[i]] <- p2  # Store the pruned plot
    
    ## Step 4: Assign LCA taxonomy ---------------------------------------------
    
    # Merge tree data with reference taxonomy information
    tax_df <- left_join(dat2, ref_tax, by = 'label')
    
    # Get the clades containing query tips (include monotypic query tips)
    suppressMessages({
      clade2 <- dat2 %>% filter(Source == 'query') %>% pull(Taxon) %>% unique()
    })
    
    # Initialize an empty dataframe to store LCA results
    result <- data.frame() %>% 
      mutate(label = NA_character_)
    
    # Assing LCA taxonomy using custom functions
    suppressMessages({
      get_lca(clade2, tax_df, tr2)
    })

    # Filter LCA results to retain only valid taxonomic assignments
    q_lca <- filter(result[, c('label', 'taxonomy', 'Taxon', 'Source', 'Subtree')], !is.na(taxonomy))
    
    # Remove trailing semicolons from taxonomy strings
    q_lca$taxonomy <- gsub(";+$", "", q_lca$taxonomy)
    
    # Extract the lowest assigned taxonomic rank
    q_lca$assign_max <- sapply(q_lca$taxonomy %>% str_split(';'), function(x) tail(x, 1))
    
    ## Step 5: Assign final phylotype taxonomy ---------------------------------
    q_lca$phylotax <- NA  # Initialize phylotype variable
    
    # If taxonomy has exactly six levels, retain the lowest taxonomic rank
    case1 <- q_lca[str_count(q_lca$taxonomy, ';') == 6,] %>% mutate(phylotax = assign_max)
    
    # If taxonomy has fewer than six levels, assign a new phylotype name
    case2 <- q_lca[str_count(q_lca$taxonomy, ';') != 6,] %>% arrange(desc(Taxon)) %>%
      group_by(Taxon, assign_max) %>%
      distinct() %>%
      mutate(phylotax = ifelse(Taxon == "Monotypic",
                               paste0(assign_max, "_pt_", Subtree, ".M.", row_number()),
                               paste0(assign_max, "_pt_", Subtree, ".", Taxon))) %>%
      ungroup()
    
    ## Result 1: All species-level ---------------------------------------------

    # This result distinguish among multiple species within a single phylotype
    # Useful when you want to know exact species information
    
    # Combine both cases into the final phylotype taxonomy table
    q_tax <- rbind(case1, case2)[, c('label', 'taxonomy', 'phylotax', 'Subtree')]
    
    # Append phylotype information to the taxonomy column
    q_tax <- q_tax %>%
      rowwise() %>%
      mutate(taxonomy = if_else(str_count(taxonomy, ";") == 6,
                                sub("([^;]+)$", phylotax, taxonomy),
                                paste0(taxonomy, strrep(";", 6 - str_count(taxonomy, ";")), phylotax))) %>%
      ungroup() %>% 
      arrange(taxonomy)
    
    # Append results to query_tax_df
    query_tax_df <- rbind(query_tax_df, q_tax)
    
    ## Result 2: Lump multiple species  ----------------------------------------
  
    # This result does not distinguish among multiple species within a single phylotype.
    # Useful for comparing multiple heterogeneous studies
    
    query_tax_lump <- q_lca %>% 
      mutate(taxonomy = str_replace(taxonomy, "^((?:[^;]*;){5}[^;]*);.*", "\\1")) # If species part present in the taxonomy column, remove

    # Revise assign_max column
    query_tax_lump$assign_max <- sapply(query_tax_lump$taxonomy %>% str_split(';'), function(x) 
      tail(x,1)) %>% as.vector()
    
    # pad ";" at the end of the taxonomy column
    query_tax_lump <- query_tax_lump %>%
      mutate(taxonomy = if_else(str_count(taxonomy, ";") < 6,
                                paste0(taxonomy, strrep(";", 6 - str_count(taxonomy, ";"))),
                                taxonomy))
    
    # Assign phylotype
    query_tax_lump <- query_tax_lump %>%
      group_by(Taxon, assign_max) %>%
      distinct() %>%
      mutate(phylotax = ifelse(Taxon == "Monotypic",
                               paste0(assign_max, "_pt_", Subtree, ".M.", row_number()),
                               paste0(assign_max, "_pt_", Subtree, ".", Taxon))) %>%
      ungroup() %>% 
      select('label', 'taxonomy', 'phylotax', 'Subtree')
    
    # Compare q_phylotax and query_tax_lump
    compare <- merge(query_tax_lump, q_tax, all.x = T, by = 'label') %>% 
      select(phylotax.x, phylotax.y) %>% 
      distinct() %>% 
      arrange(phylotax.x)
    
    # Divide compare table into lump and exact match
    lump <- compare %>% filter(duplicated(phylotax.x) | duplicated(phylotax.x, fromLast = TRUE)) %>% arrange(phylotax.x)
    exact <- setdiff(compare, lump)
    
    # Save lumping results
    lump_df <- rbind(lump_df, lump)
    
    # Replacing phylotax in query_tax_lump based on exact species mapping
    query_tax_lump$phylotax <- dplyr::coalesce(
      exact$phylotax.y[match(query_tax_lump$phylotax, exact$phylotax.x)], 
      query_tax_lump$phylotax
    )
    
    # Append phylotax to taxonomy
    query_tax_lump <- query_tax_lump %>% 
      mutate(taxonomy = paste0(taxonomy, phylotax)) %>% 
      arrange(taxonomy)
    
    # Save results to query_tax_lump_df
    query_tax_lump_df <- rbind(query_tax_lump_df, query_tax_lump)
    
    ####
    
    # Log Processing Time
    end_time <- Sys.time()  # Record end time
    cat(sprintf("Subtree %s completed in %f seconds\n", i, as.numeric(end_time - start_time, units = "secs")))

  }, error = function(e) {
    # Print error message if an error occurs
    cat("Error in processing subtree", i, "\n")
    print(e)
  })
}

## 5) Save Results ------------------------------------------------------------
write.csv(query_tax_df, paste0(path, '/taxonomy/query_tax.csv'), quote = F, row.names = F)
write.csv(query_tax_lump_df, paste0(path, '/taxonomy/query_tax_lump.csv'), quote = F, row.names = F)
rename(lump_df, 'lump' = 'phylotax.x', 'species' = 'phylotax.y')
write.csv(lump_df, paste0(path, '/taxonomy/lump.csv'), quote = F, row.names = F)

x <- list(raw_tree, raw_plot, shrunk_tree, shrunk_plot, query_tax_df, query_tax_lump_df)

saveRDS(x, file = paste0(path, '/taxonomy/LCA_taxonomy_results.rds'))