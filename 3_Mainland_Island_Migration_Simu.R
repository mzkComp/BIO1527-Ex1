library(ggplot2)

# ===== 1. Initialize =====
p_M <- 0.8    # Mainland allele frequency
p_0 <- 0.2    # Initial island allele frequency
m   <- 0.01   # Migration rate

generations <- 600

# ===== 2. Simulate the process and store =====
p_island <- numeric(generations + 1)
p_island[1] <- p_0

for (t in 1:generations) {
  p_island[t + 1] <- (1 - m) * p_island[t] + m * p_M
}

## Data frame for plot
df <- data.frame(
  Generation = 0:generations,
  Frequency = p_island
)

## Find the exact data point where generation == 69 to highlight it
target_point <- df[df$Generation == 69, ]

# ==== 3. Plot =====
ggplot(df, aes(x = Generation, y = Frequency)) +
  
  geom_line(color = "blue", linewidth = 1) +
  
  ## Reference lines and annotations
  geom_hline(yintercept = p_M, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = p_0, linetype = "dashed", color = "gray", linewidth = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
  
  geom_vline(xintercept = 69, linetype = "dashed", color = "darkgreen", linewidth = 0.6) +
  
  geom_point(data = target_point, 
             aes(x = Generation, y = Frequency), 
             color = "darkgreen", size = 4, shape = 18) +
  
  annotate("text", x = 500, y = p_M - 0.03, label = "Mainland p = 0.8", color = "red", size = 4.5) +
  annotate("text", x = 500, y = 0.5 - 0.03, label = "Target p = 0.5", color = "darkgreen", size = 4.5) +
  annotate("text", x = 500, y = p_0 + 0.03, label = "Initial p = 0.2", color = "gray", size = 4.5) +
  annotate("text", x = 69, y = 0.12, label = "Gen 69", color = "darkgreen", size = 4.5) +
  
  ## Aesthetics and theme
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-10, 610), expand = c(0, 0)) +
  labs(
    title = "Island Allele Frequency Change via Migration",
    subtitle = "Simulation by Student X, m=0.01",
    x = "Generation",
    y = "Allele Frequency (p)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "darkgray", size = 10),
        legend.position = "bottom")

ggsave(filename = "3.png", dpi = 600, path = "./graph")
