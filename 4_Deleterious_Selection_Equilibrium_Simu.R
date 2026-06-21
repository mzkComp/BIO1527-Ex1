# Load required library
library(ggplot2)

# ===== 1. Initialization =====
mu <- 2e-6           # Mutation rate (A -> a)
n_gen <- 2500

q <- numeric(n_gen)  # q[i] : deleterious allele frequency of i-th gen
q[1] <- 1e-8         

w11 <- 1.00
w12 <- 0.95
w22 <- 0.00

## Calculate theoretical equilibrium
## By solving quadratic: q^2 + 0.05q - 2e-6 = 0
q_eq_new <- (-0.05 + sqrt(0.05^2 + 4 * mu)) / 2
q_eq_old <- sqrt(mu) # Old equilibrium, (w11, w12, w22) = (1, 1, 0)

# ===== 2. Simulation Loop (General Discrete Generation Recursion) =====
for (i in 1:(n_gen - 1)) {
  p <- 1 - q[i]
  w_bar <- (p^2 * w11) + (2 * p * q[i] * w12) + (q[i]^2 * w22)
  
  ## Selection
  q_sel <- ((p * q[i] * w12) + (q[i]^2 * w22)) / w_bar
  
  ## Mutation
  q[i + 1] <- q_sel * (1 - mu) + (1 - q_sel) * mu
}

q_tmp <- q # Temporately store data of first simu

q <- numeric(n_gen) 
q[1] <- 1e-8 
w11 <- 1.00
w12 <- 1.00
w22 <- 0.00

for (i in 1:(n_gen - 1)) {
  p <- 1 - q[i]
  w_bar <- (p^2 * w11) + (2 * p * q[i] * w12) + (q[i]^2 * w22)
  
  ## Selection
  q_sel <- ((p * q[i] * w12) + (q[i]^2 * w22)) / w_bar
  
  ## Mutation
  q[i + 1] <- q_sel * (1 - mu) + (1 - q_sel) * mu
}

df <- data.frame(
  Generation = 1:n_gen,
  Allele_Frequency_1 = q_tmp, 
  Allele_Frequency_2 = q
) 

df |> ggplot(aes(x = Generation)) +
  geom_line(aes(y = Allele_Frequency_1, color = 'w12=0.95'), linewidth = 1) +
  geom_line(aes(y = Allele_Frequency_2, color = 'w12=1.00'), linewidth = 1) +
  scale_color_manual('w12', values=c('red', 'steelblue')) +
  
  # Theoretical equilibrium hline
  geom_hline(yintercept = q_eq_new, 
             linetype = "dashed", 
             color = "red", 
             alpha = 0.3, 
             linewidth = 0.5) +
  annotate("text", 
           x = n_gen * 0.65, y = q_eq_new * 2.0, 
           label = paste("Equilibrium (q ≈", round(q_eq_new, 6), ")"), 
           color = "red", size = 4.5) +
  geom_hline(yintercept = q_eq_old, 
             linetype = "dashed", 
             color = "red",  
             alpha = 0.3,
             linewidth = 0.5) +
  annotate("text", 
           x = n_gen * 0.65, y = q_eq_old * 1.03, 
           label = paste("Equilibrium (q ≈", round(q_eq_old, 6), ")"), 
           color = "red", size = 4.5) +
  
  # Labels and Theme
  labs(title = paste0("Mutation-Selection Balance", "\n", "Complete recessive w12 = 1.00, Heterozygote disadvantage w12=0.95"),
       subtitle = "Simulation by Student X",
       x = "Generation",
       y = "Allele Frequency (q)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "darkgray", size = 10),
        legend.position = "bottom")

ggsave(filename = sprintf('4.png'), dpi = 600, path = "./graph")