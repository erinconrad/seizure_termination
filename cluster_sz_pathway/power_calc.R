# Define parameters
p1 <- 0.60   # Proportion stopping in stimulation group
p2 <- 0.19   # Proportion stopping in no stimulation group
alpha <- 0.05
power <- 0.80

# Run power analysis
result <- power.prop.test(p1 = p1, p2 = p2, sig.level = alpha, power = power, alternative = "two.sided")

# Round up to nearest integer for practical use
cat("Required number of seizures per group:", ceiling(result$n), "\n")

# Need 70 per group or 140 if it's 38% vs 17%
# Need 22 per group or 44 if it's 60% vs 19%