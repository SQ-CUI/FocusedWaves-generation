#Due to the limitations in DualSPHysics case files where focused waves cannot provide focusing time and only JONSWAP and P-M spectra are available. 
#The DualSPHyiscs team have provided an external .dat file for reading wave piston displacements. 
#This code is generated using R and runs in RStudio. 
#The wave generation principle applicable to this code is derived from Chen et al.(2024)https://doi.org/10.1016/j.oceaneng.2024.118516
#The size of wave tank was from DUT by Chen et al.(2024)
# ==================== Load necessary packages ====================
library(ggplot2)
library(gridExtra)

# ==================== Parameter settings ====================
d <- 0.4                # water depth (m)
A <- 0.128              # target wave amplitude (m)
x_b <- 6.3              # focusing position (m)
t_b <- 20               # focusing time (s)
f_min <- 0.5            # minimum frequency (Hz)
f_max <- 1.5            # maximum frequency (Hz)
df <- 0.02              # frequency interval (Hz)
g <- 9.81               # gravitational acceleration (m/s^2)

# Frequency sequence
f <- seq(f_min, f_max, by = df)
N <- length(f)          # 51 components
omega <- 2 * pi * f

# ==================== Solve for wave number k (dispersion relation) ====================
k <- numeric(N)
for (i in 1:N) {
  disp_rel <- function(k) omega[i]^2 - g * k * tanh(k * d)
  k[i] <- uniroot(disp_rel, interval = c(0.1, 20), tol = 1e-8)$root
}

# ==================== Component amplitudes a_i (Equation 2) ====================
sum_invk <- sum(1 / k)
a_i <- -A / (k * sum_invk)          # unit: m

# ==================== Transfer function T(f_i) (Equation 4) ====================
T_i <- 4 * (sinh(k * d))^2 / (2 * k * d + sinh(2 * k * d))

# ==================== Amplitude coefficients for each wave paddle component ====================
C_i <- a_i / T_i

# ==================== Time series (0 ~ 20 s, step 0.02 s) ====================
t_start <- 0
t_end <- 40
dt <- 0.02
t <- seq(t_start, t_end, by = dt)
n_steps <- length(t)

# ==================== Calculate wave paddle displacement S(t) (Equation 3, x=0) ====================
S <- numeric(n_steps)
for (i in 1:N) {
  S <- S + C_i[i] * sin(-k[i] * x_b - omega[i] * (t - t_b))
}

# ==================== Calculate theoretical wave surface at focusing point η(x_b, t) (Equation 1) ====================
eta_focus <- numeric(n_steps)
for (i in 1:N) {
  eta_focus <- eta_focus + a_i[i] * cos(-omega[i] * (t - t_b))
}

# ==================== Construct data frames ====================
df <- data.frame(
  time = t,
  displacement = S,
  eta_focus = eta_focus
)

# Spectral data (for Figure 3)
spectrum_df <- data.frame(
  frequency = f,
  amplitude_mm = a_i * 1000,        # convert to mm
  transfer = T_i
)

# ==================== Output DAT file (two columns, three spaces separated, no column names) ====================
# Format time and displacement as strings with 4 decimal places
time_str <- sprintf("%.4f", t)
disp_str <- sprintf("%.4f", S)
# Join the two columns with three spaces and write to file
writeLines(paste(time_str, disp_str, sep = "   "), "paddle_displacement.dat")
cat("\nWave paddle displacement data saved to: paddle_displacement.dat\n")

# ==================== Plotting ====================
# Figure 1: Wave paddle displacement
p1 <- ggplot(df, aes(x = time, y = displacement)) +
  geom_line(color = "blue", size = 0.5) +
  geom_vline(xintercept = t_start, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(xintercept = t_b, linetype = "dashed", color = "green", alpha = 0.7) +
  annotate("text", x = t_start, y = max(df$displacement), 
           label = paste("Wave start (", t_start, "s)"), 
           color = "red", hjust = -0.1, size = 3) +
  annotate("text", x = t_b, y = max(df$displacement) * 0.9, 
           label = paste("Focusing time (", t_b, "s)"), 
           color = "green", hjust = -0.1, size = 3) +
  labs(x = "Time (s)", y = "Paddle displacement (m)", 
       title = "Wave paddle displacement time series") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)

# Figure 2: Theoretical wave surface at focusing point
p2 <- ggplot(df, aes(x = time, y = -eta_focus)) +
  geom_line(color = "green", size = 0.5) +
  geom_vline(xintercept = t_b, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_hline(yintercept = A, linetype = "dotted", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = -A, linetype = "dotted", color = "blue", alpha = 0.5) +
  annotate("text", x = t_b, y = max(df$eta_focus), 
           label = paste("Focusing time (", t_b, "s)"), 
           color = "red", hjust = -0.1, size = 3) +
  annotate("text", x = max(df$time) * 0.8, y = A * 1.1, 
           label = paste("Target amplitude ±", A, "m"), 
           color = "blue", size = 3) +
  labs(x = "Time (s)", y = "Wave surface elevation (m)", 
       title = paste("Theoretical wave surface at focusing point x =", x_b, "m")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3)

# Figure 3: Amplitude spectrum and transfer function (dual Y-axis)
# Amplitude in mm
max_amp <- max(spectrum_df$amplitude_mm)
max_trans <- max(spectrum_df$transfer)
scale_factor <- max_amp / max_trans * 2   # Scaling factor to give the transfer function curve an appropriate height relative to amplitude spectrum

p3 <- ggplot(spectrum_df, aes(x = frequency)) +
  geom_line(aes(y = amplitude_mm, color = "Wave amplitude (mm)"), size = 0.8) +
  geom_point(aes(y = amplitude_mm, color = "Wave amplitude (mm)"), size = 1.5) +
  geom_line(aes(y = transfer * scale_factor, color = "Transfer function (scaled)"), 
            size = 0.8, linetype = "dashed") +
  geom_point(aes(y = transfer * scale_factor, color = "Transfer function (scaled)"), 
             size = 1.5, shape = 17) +
  geom_hline(yintercept = A * 1000, linetype = "dashed", color = "gray", alpha = 0.7) +
  annotate("text", x = max(f) * 0.8, y = A * 1000 * 1.05, 
           label = paste("Total amplitude", round(A * 1000, 1), "mm"), 
           color = "gray", size = 3) +
  scale_y_continuous(
    name = "Wave amplitude (mm)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Transfer function T")
  ) +
  scale_color_manual(
    name = "",
    values = c("Wave amplitude (mm)" = "blue", "Transfer function (scaled)" = "red"),
    labels = c("Wave amplitude (mm)", "Transfer function")
  ) +
  labs(x = "Frequency (Hz)", title = "Wave amplitude spectrum and transfer function") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# Figure 4: FFT spectrum analysis of wave paddle displacement
fs <- 1 / dt                     # Sampling frequency
Y <- fft(S)                      # FFT
P2 <- abs(Y / n_steps)           # Two-sided amplitude spectrum
P1 <- P2[1:(n_steps/2 + 1)]      # Take one-sided
P1[2:(length(P1)-1)] <- 2 * P1[2:(length(P1)-1)]
f_fft <- fs * (0:(n_steps/2)) / n_steps
fft_df <- data.frame(frequency = f_fft, magnitude = P1)

p4 <- ggplot(fft_df %>% filter(frequency <= 2.5), aes(x = frequency, y = magnitude)) +
  geom_line(color = "purple") +
  geom_vline(xintercept = c(f_min, f_max), linetype = "dashed", color = "red", alpha = 0.5) +
  annotate("text", x = f_min, y = max(fft_df$magnitude) * 0.9, 
           label = "f_min", color = "red", hjust = 1.1, size = 3) +
  annotate("text", x = f_max, y = max(fft_df$magnitude) * 0.9, 
           label = "f_max", color = "red", hjust = -0.1, size = 3) +
  labs(x = "Frequency (Hz)", y = "Amplitude (m)", title = "Spectrum analysis of wave paddle displacement") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# ==================== Combine four plots into 2×2 layout ====================
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

# Save combined plot
ggsave("paddle_analysis.png", combined_plot, width = 12, height = 10, dpi = 150)
cat("\nAnalysis plot saved to: paddle_analysis.png\n")

# ==================== Output first 10 rows of data example ====================
cat("\nFirst 10 rows of data example:\n")
cat("Time(s)\tDisplacement(m)\n")
for (i in 1:min(10, n_steps)) {
  cat(sprintf("%.3f\t%.6f\n", t[i], S[i]))
}

# ==================== Additional analysis: Check wave surface at focusing moment ====================
t_b_idx <- which.min(abs(t - t_b))
cat(sprintf("\nAnalysis at focusing moment (t = %.2f s):\n", t[t_b_idx]))
cat(sprintf("  Wave paddle displacement: %.4f m\n", S[t_b_idx]))
cat(sprintf("  Theoretical wave surface elevation: %.4f m\n", eta_focus[t_b_idx]))
cat(sprintf("  Target amplitude A: %.4f m\n", A))

cat(sprintf("  Relative error: %.2f%%\n", abs(eta_focus[t_b_idx] - A) / A * 100))
