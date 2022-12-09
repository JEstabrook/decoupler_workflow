library(dplyr)
library(tibble)
library(reticulate)

args = commandArgs(trailingOnly=TRUE)
print(args)
expr_fname = args[1]
netw_fname = args[2]
netw_noise_out = args[3]
sec_intx_netw_noise_out = args[4]
noise_perc = args[5]
n_rep = args[6]

# Randomly sample edges from expression and randomly replace with regulons in network
expr <- readRDS(expr_fname)
netw <- readRDS(netw_fname)
edges <- rownames(expr)

# Randomly sample edges from network and regulons to replace
set.seed(42)
n_sample <- floor(nrow(netw) * as.numeric(noise_perc))
rand_edges <- sample(edges, n_sample, replace = TRUE)
rand_regs <- sample(1:nrow(netw), n_sample)

# Replace random edges
netw$target[rand_regs] <- rand_edges

# Save 
saveRDS(netw, netw_noise_out)

# Identify secondary edges for Priori
# Re-format primary interaction tibble
primary_int <- netw %>% 
  mutate("Type" = "controls-expression-of") %>% 
  select("tf", "Type", "target") %>% 
  rename(UpGene = tf, DownGene = target)

primary_int <- primary_int[order(primary_int$UpGene),]

# Initialize secondary interaction tibble
secondary_int <- tibble("UpGene" = character(), 
                        "Type" = character(),
                        "DownGene" = character())

# Loop through the regulators in the network
regulators <- unique(primary_int$UpGene)
for (reg in regulators) {
  
  # Identify primary and secondary targets
  sub_primary_targets <- unique(primary_int$DownGene[primary_int$UpGene == reg])
  sub_secondary_targets <- unique(primary_int$DownGene[primary_int$UpGene %in% sub_primary_targets])
  
  # Initialize tibble with regulator and secondary targets
  sub_secondary_int <- tibble("UpGene" = reg, 
                        "Type" = "controls-expression-of",
                        "DownGene" = sub_secondary_targets)
  sub_secondary_int <- sub_secondary_int %>% 
    distinct()
  
  # Bind to secondary interactions tibble
  secondary_int <- bind_rows(secondary_int, sub_secondary_int)
}

# Bind to primary interactions tibble
primary_secondary_int <- bind_rows(primary_int, secondary_int) %>% 
  distinct()

primary_secondary_int <- primary_secondary_int[order(primary_secondary_int$UpGene),]

# Save as pickle file
pandas <- import("pandas")
primary_secondary_int_form <- pandas$DataFrame(data = primary_secondary_int, columns = c('UpGene', 'Type', 'DownGene'))
py_save_object(primary_secondary_int, filename = sec_intx_netw_noise_out, pickle = "pickle")