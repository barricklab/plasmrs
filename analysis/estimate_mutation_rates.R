require(rsalvador)
require(dplyr)
require(ggplot2)


# mutant_plasmid_threshold is the *number* of plasmids that must be present for it to count as a
# mutant detected by the Luria-Delbrück fluctuation test
estimateMutRate <- function(population_filename, mutant_plasmid_threshold) {

    data = read.csv(filename)  

    data$key = paste(data$N_0, data$N, data$N_p, data$u_p, sep="/");
    est = data.frame()
    for (key in unique(data$key)) {
    
      key_pieces = unlist(strsplit(key, "/"))
      key_N_0 = as.numeric(key_pieces[1])
      key_N = as.numeric(key_pieces[2])
      key_N_p = as.numeric(key_pieces[3])
      key_u_p = as.numeric(key_pieces[4])
      
      key_data = data %>% filter(key==key)
      
      #Loop over replicates

      m_counts = c()
      for (on_repl in unique(key_data$repl)) {
        repl_m = 0
        repl_data = key_data %>% filter(repl==on_repl)
        print(repl_data)
        for (i in 1:nrow(repl_data)) {
          if (repl_data$p_mut[i] >= mutant_plasmid_threshold) {
            repl_m = repl_m + repl_data$cells[i]
          }
        }
        
        m_counts = c(m_counts, repl_m)
      }
      
      print(m_counts)
      
      
      key_m = newton.LD(m_counts)
      key_mu = key_m / key_N
      

      key_row = data.frame(N_0 = key_N_0, 
                           N = key_N,
                           N_p = key_N_p, 
                           u_p = key_u_p, 
                           m = key_m,
                           mu = key_mu
                           )
      est = est %>% bind_rows(key_row)
    }
    
    return(est)
}

est = estimateMutRate(population_filename = "population.csv", mutant_plasmid_threshold = 1)
