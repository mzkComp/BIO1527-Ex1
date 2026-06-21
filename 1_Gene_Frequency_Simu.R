library(ggplot2)

# ====== 1. Recursion Function ======
## Fitness formula: w11 = 1 - s, w12 = 1 - h*s, w22 = 1 (if Deleterious)
simulate_evolution <- function(p0, h, s, generations) {
  p <- numeric(generations + 1)
  p[1] <- p0
  
  for (t in 1:generations) {
    q <- 1 - p[t]
    
    if (Allele_type == 'Deleterious') {
      w11 <- 1 - s
      w12 <- 1 - h * s
      w22 <- 1
    } else {
      w11 <- 1 + s
      w12 <- 1 + h * s
      w22 <- 1
    }
    
    w_bar <- (p[t]^2 * w11) + (2 * p[t] * q * w12) + (q^2 * w22)
    p_next <- ((p[t]^2 * w11) + (p[t] * q * w12)) / w_bar
    
    # let p[t+1] in [0, 1]
    p[t + 1] <- max(0, min(1, p_next)) 
  }
  return(p)
}

# ===== 2. Initialization =====
Allele_type <- 'Deleterious' ## Allele_type: 'Deleterious' OR 'Advantageous'

if (Allele_type == 'Deleterious') {
  p0 <- 0.8
  h_values <- c(0.4, -0.4, 1.5, 1, 0)
} else {
  p0 <- 0.2
  h_values <- c(0.5, 1, 0, 0.4)
}

s_values <- c(0.05, 0.06, 0.07, 0.08, 0.09)
generations <- 1000
plot_data_list <- list() # Store data for different h-values

# ===== 3. Get Simulation Data =====
for (i in 1:length(h_values)) {
  h <- h_values[i]
  
  # Create new DataFrame
  df <- data.frame()
  for (s in s_values) {
    p_trajectory <- simulate_evolution(p0, h, s, generations)
    temp_df <- data.frame(
      Generation = 0:generations,
      Frequency = p_trajectory,
      s = factor(s) # turn into factors, as to plot
    )
    df <- rbind(df, temp_df)
  }
  
  title_text <- sprintf("1-(%d) h = %.1f", i, h) # Example: "1-(1) h = 0.4"
  
  plot_data_list[[i]] <- list(data = df, title = title_text)
}

# ===== 4. plot the data =====
for (i in 1:length(h_values)) {
  df <- plot_data_list[[i]]$data
  title_text <- plot_data_list[[i]]$title
  
  p <- ggplot(df, aes(x = Generation, y = Frequency, color = s)) +
    geom_line(linewidth = 1) +
    scale_color_viridis_d(option = "G", name = "Selection Coefficient (s)") +
    theme_minimal(base_size = 14) +
    labs(
      title = title_text,
      subtitle = "1000 Generations Simulation by Student X", 
      x = "Generation",
      y = "Allele Frequency (p)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "darkgray", size = 10),
      legend.position = "bottom"
    )
  
  ggsave(filename = sprintf('1-%d.png', i), dpi = 600, path = "./graph")
}