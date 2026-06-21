library(ggplot2)

# ===== 1. CCR5 Statistics Analysis =====
## O1 for Observed Count 1
O1 <- c(480, 280.8, 39.2)
names(O1) <- c("1/1", "1/Delta32", "Delta32/Delta32")

## Expected Probabilities, E1 for Expected Count 1 
p1 <- (2*O1[1] + O1[2]) / (2 * sum(O1))
q1 <- 1 - p1

E_prob1 <- c(p1^2, 2*p1*q1, q1^2)
E1 <- E_prob1 * sum(O1)

## Chi-Square Test
test1 <- chisq.test(O1, p = E_prob1, rescale.p = TRUE)
cat("---------- CCR5 ----------\n")
cat("Chi-Square:", round(test1$statistic, 4), "\n")
cat("P-value:", test1$p.value, "\n\n")

## Plot
df1 <- data.frame(
  Genotype = factor(c("1/1", "1/Δ32", "Δ32/Δ32"), levels = c("1/1", "1/Δ32", "Δ32/Δ32")),
  Count = c(O1, E1),
  Type = rep(c("Observed", "Expected"), each = 3)
)

plot1 <- ggplot(df1, aes(x = Genotype, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  geom_text(aes(label = round(Count, 1)), position = position_dodge(width = 0.6), vjust = -0.5) +
  labs(
    title = paste0("CCR5 genotypes"),
    subtitle = paste0("Chi-Square Analysis by Student X")) +
  ylim(0, max(df1$Count) * 1.15) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") + 
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "darkgray", size = 10),
    legend.position = "bottom"
  )
ggsave(filename = sprintf('2-1.png'), plot = plot1, dpi = 600, path = "./graph")

# ===== 2. Sickle-cell hemoglobin Statistics Analysis =====
O2 <- c(907.2, 290.4, 2.4)
names(O2) <- c("AA", "AS", "SS")

p2 <- (2*O2[1] + O2[2]) / (2 * sum(O2))
q2 <- 1 - p2

E_prob2 <- c(p2^2, 2*p2*q2, q2^2)
E2 <- E_prob2 * sum(O2)

test2 <- chisq.test(O2, p = E_prob2, rescale.p = TRUE)
cat("\n---------- Sickle-cell hemoglobin ----------\n")
cat("Chi-square:", round(test2$statistic, 4), "\n")
cat("P-value:", test2$p.value, "\n\n")

df2 <- data.frame(
  Genotype = factor(c("AA", "AS", "SS"), levels = c("AA", "AS", "SS")),
  Count = c(O2, E2),
  Type = rep(c("Observed", "Expected"), each = 3)
)

plot2 <- ggplot(df2, aes(x = Genotype, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  geom_text(aes(label = round(Count, 1)), position = position_dodge(width = 0.6), vjust = -0.5) +
  labs(
    title = paste0("Sickle-cell hemoglobin"),
    subtitle = paste0("Chi-Square Analysis by Student X")) +
  ylim(0, max(df2$Count) * 1.15) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") + 
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "darkgray", size = 10),
    legend.position = "bottom"
  )

ggsave(filename = sprintf('2-2.png'), plot = plot2, dpi = 600, path = "./graph")