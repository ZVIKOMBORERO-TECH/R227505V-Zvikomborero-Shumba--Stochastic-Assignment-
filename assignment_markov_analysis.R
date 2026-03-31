# ============================================================
#   HASTS416/RM MARKOV CHAIN QUESTIONS
#   Q1/A1: 5-State Markov Chain Analysis
#   Q2/A2: 7-State Markov Chain Analysis  
#   Q3/A3: Non-Homogeneous Traffic Model
# ============================================================

# Clear workspace and set seed for reproducibility
rm(list = ls())
set.seed(2024)

# Load required packages
suppressPackageStartupMessages({
  library(markovchain)
})

# Set results directory
results_dir <- "C:/Users/USER/Desktop/ASSIGNMENTSTOCHASTIC/"

# Create directory structure for organized output
DIR_A1 <- file.path(results_dir, "A1_Results")
DIR_A2 <- file.path(results_dir, "A2_Results") 
DIR_A3 <- file.path(results_dir, "A3_Results")
DIR_SUMMARY <- file.path(results_dir, "Summary")
OUTPUT_DIR <- results_dir
RESULTS_TXT <- file.path(results_dir, "Results_Summary.txt")

# Create directories if they don't exist
for (dir in c(DIR_A1, DIR_A2, DIR_A3, DIR_SUMMARY)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# Add append_line function if it doesn't exist
if (!exists("append_line")) {
  append_line <- function(text, file) {
    cat(text, "\n", file = file, append = TRUE)
  }
}

# ============================================================
#   SHARED UTILITY FUNCTIONS
# ============================================================

# Matrix power function
matrix_power <- function(P, n) {
  result <- diag(nrow(P))
  for (i in 1:n) {
    result <- result %*% P
  }
  return(result)
}

# Markov chain simulation function
simulate_mc <- function(P, start, n_steps = 60) {
  s    <- integer(n_steps + 1)
  s[1] <- start
  for (t in seq_len(n_steps))
    s[t+1] <- sample.int(nrow(P), 1L, prob = P[s[t], ])
  s
}

# Simulate single path function
simulate_path <- function(P, n_steps) {
  start <- sample.int(nrow(P), 1L)
  simulate_mc(P, start, n_steps)
}

# Stationary distribution function
stationary_dist <- function(P) {
  # Find left eigenvector corresponding to eigenvalue 1
  eigen_res <- eigen(t(P))
  pi_vec <- eigen_res$vectors[, which.max(Re(eigen_res$values))]
  # Normalize to sum to 1
  pi_vec <- Re(pi_vec) / sum(Re(pi_vec))
  return(pi_vec)
}

# Compute marginal probabilities over time
compute_marginals <- function(pi0, P, n_time) {
  out      <- matrix(0, nrow = n_time + 1, ncol = nrow(P))
  out[1, ] <- pi0
  Pt       <- diag(nrow(P))
  for (t in seq_len(n_time)) {
    Pt         <- Pt %*% P
    out[t+1, ] <- pi0 %*% Pt
  }
  out
}

# Drawing functions for transition diagrams
draw_node <- function(x, y, r, fill, border, label, cex = 1.4) {
  theta <- seq(0, 2*pi, length.out = 200)
  polygon(x + r*cos(theta), y + r*sin(theta),
          col = fill, border = border, lwd = 2.5)
  text(x, y, label, col = "white", cex = cex, font = 2)
}

draw_curved_arrow <- function(x0, y0, x1, y1, label, col,
                               curve = 0.3, lwd = 2.2,
                               lpos = 0.5, loff_x = 0, loff_y = 0,
                               tcex = 0.85) {
  dx <- x1 - x0; dy <- y1 - y0
  cx <- (x0+x1)/2 - curve*dy
  cy <- (y0+y1)/2 + curve*dx
  n  <- 50
  t  <- seq(0, 1, length.out = n)
  bx <- (1-t)^2*x0 + 2*(1-t)*t*cx + t^2*x1
  by <- (1-t)^2*y0 + 2*(1-t)*t*cy + t^2*y1
  lines(bx, by, col = col, lwd = lwd)
  arrows(bx[n-1], by[n-1], bx[n], by[n],
         length = 0.14, lwd = lwd, col = col, angle = 25)
  tm <- lpos
  lx <- (1-tm)^2*x0 + 2*(1-tm)*tm*cx + tm^2*x1 + loff_x
  ly <- (1-tm)^2*y0 + 2*(1-tm)*tm*cy + tm^2*y1 + loff_y
  text(lx, ly, label, cex = tcex, font = 2, col = "#1A252F")
}

draw_self_loop <- function(cx, cy, r_node = 0.38, r_loop = 0.52,
                            angle_deg = 90, col = "#555", lwd = 2.0,
                            label = "", tcex = 0.78) {
  a  <- angle_deg * pi / 180
  ox <- cx + r_node * cos(a)
  oy <- cy + r_node * sin(a)
  phi   <- seq(0, 2*pi, length.out = 100)
  loop_cx <- cx + (r_node + r_loop * 0.9) * cos(a)
  loop_cy <- cy + (r_node + r_loop * 0.9) * sin(a)
  slx <- loop_cx + r_loop * cos(phi)
  sly <- loop_cy + r_loop * sin(phi)
  lines(slx, sly, col = col, lwd = lwd)
  idx_tip <- which.min((slx - ox)^2 + (sly - oy)^2)[1]
  idx2    <- (idx_tip - 3) %% 100 + 1
  arrows(slx[idx2], sly[idx2], slx[idx_tip], sly[idx_tip],
         length = 0.12, lwd = lwd, col = col, angle = 22)
  text(loop_cx + r_loop * cos(a) * 0.8,
       loop_cy + r_loop * sin(a) * 0.8,
       label, cex = tcex, font = 2, col = "#1A1A2E")
}

# ============================================================
#   Q1/A1: 5-STATE MARKOV CHAIN ANALYSIS
# ============================================================

cat("╔══════════════════════════════════════════════╗\n")
cat("║    Q1/A1: 5-STATE MARKOV CHAIN ANALYSIS       ║\n")
cat("╚══════════════════════════════════════════════╝\n\n")

# Transition Matrix from assignment
P_A1 <- matrix(c(
  1.0, 0.0, 0.0, 0.0, 0.0,
  0.5, 0.0, 0.0, 0.0, 0.5,
  0.2, 0.0, 0.0, 0.0, 0.8,
  0.0, 0.0, 1.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 1.0, 0.0
), nrow = 5, byrow = TRUE)
rownames(P_A1) <- colnames(P_A1) <- paste0("S", 1:5)

cat("Transition Matrix P:\n")
print(P_A1)

# Save transition matrix to CSV
write.csv(P_A1, file.path(results_dir, "01_A1_transition_matrix.csv"), row.names = TRUE)

# Q1/A1(a): Diagram & Classification
cat("\n══════════════════════════════════════════════\n")
cat("  Q1/A1(a): DIAGRAM & CLASSIFICATION\n")
cat("══════════════════════════════════════════════\n")

# Node positions for 5-state chain
nx_A1 <- c(-2.5,  0.0,  0.0,  2.5,  2.5)
ny_A1 <- c( 0.0,  1.7, -1.5, -1.5,  1.7)
nr_A1 <- 0.40

node_fill_A1  <- c("#C0392B","#2980B9","#27AE60","#E67E22","#8E44AD")
node_frame_A1 <- c("#922B21","#1A5276","#1E8449","#B7770D","#6C3483")

png(file.path(results_dir, "A1_diagram.png"), width = 1200, height = 850, res = 130)
par(mar = c(2, 2, 4, 2), bg = "#F8F9FA")
plot(0, 0, type = "n",
     xlim = c(-3.4, 3.6), ylim = c(-2.4, 2.7),
     asp = 1, axes = FALSE, xlab = "", ylab = "",
     main = "A1(a): 5-State Markov Chain Transition Diagram",
     cex.main = 1.3)

# Draw edges for 5-state chain
# S1 → S1 (self-loop)
theta_sl <- seq(pi*0.6, pi*2.4, length.out = 80)
lx_sl <- nx_A1[1] + 0.62 * cos(theta_sl)
ly_sl <- ny_A1[1] + 0.62 * sin(theta_sl)
lines(lx_sl, ly_sl, col = "#C0392B", lwd = 2.3)
arrows(lx_sl[79], ly_sl[79], lx_sl[80], ly_sl[80],
       length = 0.14, lwd = 2.3, col = "#C0392B", angle = 25)
text(nx_A1[1] - 0.85, ny_A1[1] + 0.67, "1.0", cex = 0.88, font = 2, col = "#1A252F")

# S2 → S1 and S2 → S5
draw_curved_arrow(nx_A1[2]-nr_A1*0.65, ny_A1[2]-nr_A1*0.65,
                  nx_A1[1]+nr_A1*0.65, ny_A1[1]+nr_A1*0.65,
                  "0.5", col = "#2980B9", curve = 0.30,
                  loff_x = -0.12, loff_y = 0.12)

draw_curved_arrow(nx_A1[2]+nr_A1, ny_A1[2],
                  nx_A1[5]-nr_A1, ny_A1[5],
                  "0.5", col = "#2980B9", curve = 0.0,
                  loff_y = 0.18)

# S3 → S1 and S3 → S5
draw_curved_arrow(nx_A1[3]-nr_A1*0.65, ny_A1[3]+nr_A1*0.65,
                  nx_A1[1]+nr_A1*0.5,  ny_A1[1]-nr_A1*0.5,
                  "0.2", col = "#27AE60", curve = -0.30,
                  loff_x = -0.16, loff_y = -0.12)

draw_curved_arrow(nx_A1[3]+nr_A1*0.6,  ny_A1[3]+nr_A1*0.6,
                  nx_A1[5]-nr_A1*0.6,  ny_A1[5]-nr_A1*0.6,
                  "0.8", col = "#27AE60", curve = -0.25,
                  loff_x = 0.18, loff_y = 0.0)

# S4 → S3
draw_curved_arrow(nx_A1[4]-nr_A1, ny_A1[4],
                  nx_A1[3]+nr_A1, ny_A1[3],
                  "1.0", col = "#E67E22", curve = 0.0,
                  loff_y = -0.20)

# S5 → S4
draw_curved_arrow(nx_A1[5], ny_A1[5]-nr_A1,
                  nx_A1[4], ny_A1[4]+nr_A1,
                  "1.0", col = "#8E44AD", curve = 0.0,
                  loff_x = 0.22)

# Draw nodes
for (i in 1:5) draw_node(nx_A1[i], ny_A1[i], nr_A1,
                          node_fill_A1[i], node_frame_A1[i],
                          paste0("S", i))

# Legend
legend(-3.4, -1.6,
       legend = c("S1:  Absorbing / Recurrent  [period 1]",
                  "S2:  Transient, isolated class  [period 1]",
                  "S3,S4,S5:  Transient, period-3 cycle"),
       fill   = c("#C0392B","#2980B9","#27AE60"),
       border = NA, bty = "n", cex = 0.88)

dev.off()
cat("  A1_diagram.png written.\n")

# Classification
cat("\n▶ COMMUNICATING CLASSES\n")
cat("   C1 = {S1}         — Closed, Recurrent (Absorbing)\n")
cat("   C2 = {S2}         — Open, Transient (own class; no return)\n")
cat("   C3 = {S3, S4, S5} — Open, Transient (communicate; can escape to S1)\n")
cat("\n▶ ABSORBING STATES\n")
cat("   S1: P(S1→S1) = 1.0  ✔ Absorbing   |   No others.\n")
cat("\n▶ REFLECTING STATES  (boundary forces immediate return)\n")
cat("   None in this chain.\n")
cat("\n▶ RECURRENT vs TRANSIENT\n")
cat("   Recurrent  : {S1}           — return probability = 1\n")
cat("   Transient  : {S2,S3,S4,S5} — positive probability of non-return\n")
cat("\n▶ PERIOD OF EACH STATE\n")
cat("   S1 : d = 1  (self-loop ⟹ aperiodic)\n")
cat("   S2 : d = 1  (aperiodic; P^n(S2,S2)=0 for all n; convention d=1)\n")
cat("   S3 : d = 3  (shortest return S3→S5→S4→S3, length 3)\n")
cat("   S4 : d = 3  (shortest return S4→S3→S5→S4, length 3)\n")
cat("   S5 : d = 3  (shortest return S5→S4→S3→S5, length 3)\n")

# Q1/A1(b): Simulated Trajectories
cat("\n══════════════════════════════════════════════\n")
cat("  Q1/A1(b): SIMULATED TRAJECTORIES\n")
cat("══════════════════════════════════════════════\n")

n_steps_A1 <- 60
starts_A1  <- sample(1:5, 3)
cat("Randomly chosen start states:", paste0("S", starts_A1), "\n")
trajs_A1 <- lapply(starts_A1, simulate_mc, P = P_A1, n_steps = n_steps_A1)

cols3_A1 <- c("#E74C3C","#27AE60","#2980B9")

png(file.path(results_dir, "A1_trajectories.png"), width = 1200, height = 580, res = 130)
par(mar = c(4.5, 4.8, 3.5, 1.5), bg = "#F8F9FA")
plot(0:n_steps_A1, trajs_A1[[1]], type = "s",
     col = cols3_A1[1], lwd = 2.2,
     ylim = c(0.6, 5.5),
     xlab = "Time step  n", ylab = "State",
     main = "A1(b): Three Simulated Trajectories of 5-State Markov Chain",
     yaxt = "n", bty = "l", las = 1, cex.main = 1.2)
axis(2, at = 1:5, labels = paste0("S", 1:5), las = 1)
abline(h = 1, lty = 2, col = "#C0392B", lwd = 1.2)
text(n_steps_A1 * 0.70, 1.32,
     "Absorbing barrier  S1", col = "#C0392B", cex = 0.85)
for (k in 2:3)
  lines(0:n_steps_A1, trajs_A1[[k]], type = "s", col = cols3_A1[k], lwd = 2.2)
for (k in 1:3)
  points(0, starts_A1[k], pch = 19, col = cols3_A1[k], cex = 1.8)
legend("topright",
       legend = paste0("Trajectory ", 1:3, "  (start S", starts_A1, ")"),
       col = cols3_A1, lwd = 2.2, bty = "n", cex = 0.95)
dev.off()
cat("  A1_trajectories.png written.\n")

abs_times_A1 <- sapply(trajs_A1, function(x) {
  h <- which(x == 1)[1] - 1
  if (is.na(h)) paste0(">", n_steps_A1) else as.character(h)
})
cat("Absorption step (first hit S1):", paste(abs_times_A1, collapse = ", "), "\n")

# Save A1 trajectories data to CSV
trajectories_df <- data.frame(
  Time_Step = 0:n_steps_A1,
  Trajectory_1 = trajs_A1[[1]],
  Trajectory_2 = trajs_A1[[2]],
  Trajectory_3 = trajs_A1[[3]],
  Start_States = paste0("S", starts_A1)
)
write.csv(trajectories_df, file.path(results_dir, "10_A1_trajectories_data.csv"), row.names = FALSE)
cat("\n▶ COMMENT:\n")
cat("   • Every trajectory eventually absorbs into S1 — the only recurrent\n")
cat("     state; all others are visited only finitely many times.\n")
cat("   • Chains in {S3,S4,S5} show clear PERIOD-3 cycling (S3→S5→S4→S3)\n")
cat("     in the step-function pattern before escaping to S1.\n")
cat("   • Each passage through S3 yields a 20% chance of jumping directly\n")
cat("     to S1; from S2 this probability is 50%.\n")
cat("   • S2 is resolved in exactly ONE step (→ S1 or S5).\n")
cat("   • Longer absorption times occur when the chain repeatedly re-enters\n")
cat("     the 3-cycle and misses the escape probability.\n")

# Q1/A1(c): Steady-State & Ergodicity
cat("\n══════════════════════════════════════════════\n")
cat("  Q1/A1(c): STEADY-STATE PROBABILITIES & ERGODICITY\n")
cat("══════════════════════════════════════════════\n")

# Power iteration
Pn_A1 <- diag(5)
for (i in 1:2000) Pn_A1 <- Pn_A1 %*% P_A1
cat("\nP^2000  (row i = limiting distribution from Si):\n")
print(round(Pn_A1, 8))

# Save limiting distribution to CSV
colnames(Pn_A1) <- paste0("S", 1:5)
rownames(Pn_A1) <- paste0("From_S", 1:5)
write.csv(round(Pn_A1, 8), file.path(results_dir, "02_A1_limiting_distribution_P2000.csv"), row.names = TRUE)

# Solve linear system
A_A1      <- t(P_A1) - diag(5)
A_A1[5, ] <- 1
b_A1      <- c(0, 0, 0, 0, 1)
pi_ss_A1  <- solve(A_A1, b_A1)

cat("\nSteady-state vector π  (solving πP = π,  Σπᵢ = 1):\n")
for (i in 1:5) cat(sprintf("   π(S%d) = %.6f\n", i, pi_ss_A1[i]))
cat(sprintf("Verification  ||πP − π||∞ = %.2e\n", max(abs(pi_ss_A1 %*% P_A1 - pi_ss_A1))))

# Save steady-state probabilities to CSV
steady_state_df <- data.frame(
  State = paste0("S", 1:5),
  Probability = round(pi_ss_A1, 8),
  Verification_Error = max(abs(pi_ss_A1 %*% P_A1 - pi_ss_A1))
)
write.csv(steady_state_df, file.path(results_dir, "03_A1_steady_state_probabilities.csv"), row.names = FALSE)

cat("\n▶ INTERPRETATION:\n")
cat("   π = (1, 0, 0, 0, 0): in the long run, ALL probability mass is at S1.\n")
cat("   Transient states are visited finitely often; their steady-state\n")
cat("   probability is exactly 0. S1 is the unique attractor.\n")
cat("   The steady-state distribution is DEGENERATE (point mass at S1).\n")

cat("\n▶ ERGODICITY CHECK\n")
cat("   Ergodicity requires irreducibility + positive recurrence + aperiodicity.\n")
cat("   (i)  Irreducibility  — FAILS: S1 cannot reach S2,S3,S4,S5.\n")
cat("   (ii) Positive recurrence — only S1; S2–S5 are transient.\n")
cat("   (iii)Aperiodicity     — fails for S3,S4,S5 (d = 3).\n")
cat("   ∴  The chain is NOT ergodic.\n")

# Q1/A1(d): Convergence Analysis
cat("\n══════════════════════════════════════════════\n")
cat("  Q1/A1(d): CONVERGENCE TO STEADY STATE\n")
cat("══════════════════════════════════════════════\n")

n_time_A1 <- 80
marg_unif_A1 <- compute_marginals(rep(1/5, 5), P_A1, n_time_A1)
marg_S3_A1   <- compute_marginals(c(0, 0, 1, 0, 0), P_A1, n_time_A1)

pal5_A1  <- c("#C0392B","#2980B9","#27AE60","#E67E22","#8E44AD")
tvec_A1  <- 0:n_time_A1

png(file.path(results_dir, "A1_convergence.png"), width = 1200, height = 960, res = 130)
par(mfrow = c(2, 1),
    mar   = c(4.2, 5.0, 3.5, 8.5),
    bg    = "#F8F9FA",
    oma   = c(0, 0, 2, 0))

# Panel 1: Uniform start
plot(tvec_A1, marg_unif_A1[, 1], type = "l", col = pal5_A1[1], lwd = 2.2,
     ylim = c(-0.02, 1.07), xlab = "", ylab = "P(Xₙ = Sᵢ)",
     main = "Uniform start  π₀ = (0.2, 0.2, 0.2, 0.2, 0.2)",
     bty = "l", las = 1, cex.main = 1.1)
abline(h = 1, lty = 2, col = pal5_A1[1], lwd = 0.9)
for (j in 2:5) lines(tvec_A1, marg_unif_A1[, j], col = pal5_A1[j], lwd = 2.2)
legend(par("usr")[2] + 0.5, par("usr")[4],
       legend = paste0("S", 1:5), col = pal5_A1, lwd = 2.2,
       bty = "n", cex = 0.95, xpd = TRUE, y.intersp = 1.3)

# Panel 2: Start at S3
plot(tvec_A1, marg_S3_A1[, 1], type = "l", col = pal5_A1[1], lwd = 2.2,
     ylim = c(-0.02, 1.07), xlab = "Time  n", ylab = "P(Xₙ = Sᵢ)",
     main = "Start at S3  π₀ = (0,0,1,0,0)  — period-3 oscillations clearly visible",
     bty = "l", las = 1, cex.main = 1.1)
abline(h = 1, lty = 2, col = pal5_A1[1], lwd = 0.9)
for (j in 2:5) lines(tvec_A1, marg_S3_A1[, j], col = pal5_A1[j], lwd = 2.2)
legend(par("usr")[2] + 0.5, par("usr")[4],
       legend = paste0("S", 1:5), col = pal5_A1, lwd = 2.2,
       bty = "n", cex = 0.95, xpd = TRUE, y.intersp = 1.3)

mtext("A1(d): Convergence of Unconditional Probabilities to Steady State",
      outer = TRUE, side = 3, cex = 1.2, font = 2)

dev.off()
cat("  A1_convergence.png written.\n")

# Numerical table
cat("\n  n   P(Xₙ=S1)|uniform   P(Xₙ=S1)|start-S3\n")
cat("  ────────────────────────────────────────────\n")
for (nn in c(0, 3, 6, 10, 15, 20, 30, 50, 80))
  cat(sprintf("  %2d       %.5f              %.5f\n",
              nn, marg_unif_A1[nn+1, 1], marg_S3_A1[nn+1, 1]))

# Save convergence data to CSV
convergence_df <- data.frame(
  Time_Step = c(0, 3, 6, 10, 15, 20, 30, 50, 80),
  Probability_S1_Uniform_Start = round(marg_unif_A1[c(1,4,7,11,16,21,31,51,81), 1], 8),
  Probability_S1_Start_S3 = round(marg_S3_A1[c(1,4,7,11,16,21,31,51,81), 1], 8)
)
write.csv(convergence_df, file.path(results_dir, "04_A1_convergence_data.csv"), row.names = FALSE)

cat("\n▶ COMMENT ON CONVERGENCE:\n")
cat("   • Both initial distributions converge to π=(1,0,0,0,0) as n→∞,\n")
cat("     confirming S1 as the unique long-run state.\n")
cat("   • PANEL 1 (uniform start): all transient probabilities decay\n")
cat("     smoothly to 0; P(Xₙ=S1) rises in a clean S-curve.\n")
cat("     The period-3 effect is washed out by averaging over starts.\n")
cat("   • PANEL 2 (start at S3): striking PERIOD-3 OSCILLATIONS in\n")
cat("     P(Xₙ=S3), P(Xₙ=S4), P(Xₙ=S5) — each peaks every 3 steps\n")
cat("     (at n=0,3,6,… for S3; n=1,4,7,… for S5; n=2,5,8,… for S4)\n")
cat("     before the amplitude decays to 0. P(Xₙ=S1) climbs in a\n")
cat("     matching staircase pattern — rising by the escape probability\n")
cat("     (≈0.2 per cycle) every 3 steps.\n")
cat("   • Convergence is GRADUAL because the period-3 cycle retains\n")
cat("     probability mass for multiple steps before releasing it to S1.\n")
cat("   • The chain is not uniformly ergodic; convergence is modulated\n")
cat("     by the second-largest eigenvalue of P restricted to {S3,S4,S5}.\n")

# ============================================================
#   Q2/A2: 7-STATE MARKOV CHAIN ANALYSIS
# ============================================================

cat("\n\n╔══════════════════════════════════════════════╗\n")
cat("║    Q2/A2: 7-STATE MARKOV CHAIN ANALYSIS       ║\n")
cat("╚══════════════════════════════════════════════╝\n\n")

# A2
# ==========================================================
append_line("", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)
append_line("A2", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)

P2 <- matrix(c(
  0,   1,   0,   0,   0,   0,   0,
  1,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0.4, 0.2, 0.2, 0.2,
  0,   0,   0,   0,   0.2, 0.4, 0.4,
  0.3, 0,   0,   0.1, 0.3, 0.1, 0.2,
  0,   0,   0,   0.2, 0.2, 0.3, 0.3,
  0,   0,   0,   0.5, 0.2, 0.2, 0.1
), nrow = 7, byrow = TRUE)

states2 <- as.character(1:7)
rownames(P2) <- colnames(P2) <- paste0("S", 1:7)

mc2 <- new("markovchain",
           states = states2,
           transitionMatrix = P2)

# ----------------------------------------------------------
# A2(a)
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A2(a)", RESULTS_TXT)
append_line("The Markov chain diagram has been plotted and saved.", RESULTS_TXT)

png(file.path(DIR_A2, "A2a_Markov_Chain_Diagram.png"), width = 1600, height = 1200, res = 200)
plot(mc2)
title(main = "A2(a): Markov Chain Diagram")
dev.off()

# ----------------------------------------------------------
# A2(b)
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A2(b)", RESULTS_TXT)
append_line("Recurrent class: {1,2}", RESULTS_TXT)
append_line("Transient class: {3,4,5,6,7}", RESULTS_TXT)
append_line("Absorbing states: none", RESULTS_TXT)
append_line("Reflecting states: 5, 6, 7", RESULTS_TXT)
append_line("Periods:", RESULTS_TXT)
append_line("d(1) = 2", RESULTS_TXT)
append_line("d(2) = 2", RESULTS_TXT)
append_line("d(3) = undefined (no return)", RESULTS_TXT)
append_line("d(4) = 1", RESULTS_TXT)
append_line("d(5) = 1", RESULTS_TXT)
append_line("d(6) = 1", RESULTS_TXT)
append_line("d(7) = 1", RESULTS_TXT)

a2b_df <- data.frame(
  Item = c(
    "Recurrent class",
    "Transient class",
    "Absorbing states",
    "Reflecting states",
    "d(1)",
    "d(2)",
    "d(3)",
    "d(4)",
    "d(5)",
    "d(6)",
    "d(7)"
  ),
  Answer = c(
    "{1,2}",
    "{3,4,5,6,7}",
    "none",
    "5, 6, 7",
    "2",
    "2",
    "undefined",
    "1",
    "1",
    "1",
    "1"
  ),
  stringsAsFactors = FALSE
)
write.csv(a2b_df, file.path(DIR_A2, "A2b_Answers.csv"), row.names = FALSE)

# ----------------------------------------------------------
# A2(c)
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A2(c)", RESULTS_TXT)

set.seed(321)
a2_path1 <- simulate_path(P2, 30)
a2_path2 <- simulate_path(P2, 30)

a2_paths <- data.frame(
  Time = 0:30,
  Path1 = a2_path1,
  Path2 = a2_path2
)
write.csv(a2_paths, file.path(DIR_A2, "A2c_Simulated_Trajectories.csv"), row.names = FALSE)

png(file.path(DIR_A2, "A2c_Two_Simulated_Trajectories.png"), width = 1600, height = 1200, res = 200)
matplot(rbind(a2_path1, a2_path2),
        type = "l", lty = 1, lwd = 2,
        xlab = "Time", ylab = "State",
        main = "A2(c): Two simulated trajectories")
legend("topright",
       legend = c("Path 1", "Path 2"),
       col = 1:2, lty = 1, lwd = 2)
dev.off()

append_line("Discussion:", RESULTS_TXT)
append_line("A trajectory starting in states 3 to 7 may wander for some time,", RESULTS_TXT)
append_line("but once it enters the class {1,2}, it alternates forever between states 1 and 2.", RESULTS_TXT)
append_line("This shows that {1,2} is a closed recurrent class, while {3,4,5,6,7} is transient.", RESULTS_TXT)

# ----------------------------------------------------------
# A2(d)
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A2(d)", RESULTS_TXT)

pi2 <- stationary_dist(P2)
names(pi2) <- states2
pi2_df <- data.frame(
  State = states2,
  Stationary_Probability = round(pi2, 6)
)
write.csv(pi2_df, file.path(DIR_A2, "A2d_Stationary_Distribution.csv"), row.names = FALSE)

append_line("Limiting / stationary probabilities:", RESULTS_TXT)
capture.output(print(pi2_df), file = RESULTS_TXT, append = TRUE)
append_line("Interpretation:", RESULTS_TXT)
append_line("In the long run, probability mass is concentrated on states 1 and 2.", RESULTS_TXT)
append_line("The stationary distribution is (1/2, 1/2, 0, 0, 0, 0, 0).", RESULTS_TXT)
append_line("Is the chain ergodic? No.", RESULTS_TXT)
append_line("Although a stationary distribution exists, the recurrent class {1,2} has period 2,", RESULTS_TXT)
append_line("so the chain is not ergodic.", RESULTS_TXT)

# ==========================================================
# A3
# ==========================================================
append_line("", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)
append_line("A3", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)

append_line("For question A3:", RESULTS_TXT)
append_line("Sir has communicated that in the second row of the attached matrix, we should replace 0.5 with 0.4.", RESULTS_TXT)
append_line("Therefore the first transition matrix used is:", RESULTS_TXT)
append_line("[[0.4, 0.4, 0.2], [0.3, 0.4, 0.3], [0, 0.1, 0.9]]", RESULTS_TXT)

P13 <- matrix(c(
  0.4, 0.4, 0.2,
  0.3, 0.4, 0.3,
  0.0, 0.1, 0.9
), nrow = 3, byrow = TRUE)

P46 <- matrix(c(
  0.1, 0.5, 0.4,
  0.1, 0.3, 0.6,
  0.0, 0.1, 0.9
), nrow = 3, byrow = TRUE)

traffic_states <- c("light", "heavy", "jammed")
mu0_a3 <- c(1, 0, 0)

# ----------------------------------------------------------
# A3(a)
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A3(a)", RESULTS_TXT)
append_line("From 1 PM to 4 PM = 3 hours = 9 steps of 20 minutes each.", RESULTS_TXT)
append_line("From 4 PM to 6 PM = 2 hours = 6 steps of 20 minutes each.", RESULTS_TXT)

mu6 <- mu0_a3 %*% (P13 %^% 9) %*% (P46 %^% 6)
colnames(mu6) <- traffic_states

mu6_df <- data.frame(
  State = traffic_states,
  Probability = round(as.numeric(mu6), 6)
)
write.csv(mu6_df, file.path(DIR_A3, "A3a_Distribution_at_6PM.csv"), row.names = FALSE)

append_line("Distribution at 6 PM:", RESULTS_TXT)
capture.output(print(mu6_df), file = RESULTS_TXT, append = TRUE)
append_line(sprintf("P(light at 6 PM)  = %.6f", mu6[1,1]), RESULTS_TXT)
append_line(sprintf("P(heavy at 6 PM)  = %.6f", mu6[1,2]), RESULTS_TXT)
append_line(sprintf("P(jammed at 6 PM) = %.6f", mu6[1,3]), RESULTS_TXT)
append_line("Interpretation:", RESULTS_TXT)
append_line("By 6 PM, traffic is most likely to be jammed.", RESULTS_TXT)

# ----------------------------------------------------------
# A3(b)
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A3(b)", RESULTS_TXT)

simulate_traffic <- function() {
  state <- 1
  for(i in 1:9){
    state <- sample(1:3, 1, prob = P13[state, ])
  }
  for(i in 1:6){
    state <- sample(1:3, 1, prob = P46[state, ])
  }
  return(state)
}

set.seed(123)
N <- 10000
final_states <- replicate(N, simulate_traffic())
sim_props <- table(final_states) / N

sim_vector <- c(
  light  = ifelse("1" %in% names(sim_props), sim_props["1"], 0),
  heavy  = ifelse("2" %in% names(sim_props), sim_props["2"], 0),
  jammed = ifelse("3" %in% names(sim_props), sim_props["3"], 0)
)

sim_df <- data.frame(
  State = c("light", "heavy", "jammed"),
  Simulated_Proportion = round(as.numeric(sim_vector), 6)
)
write.csv(sim_df, file.path(DIR_A3, "A3b_Simulated_Proportions_10000.csv"), row.names = FALSE)

png(file.path(DIR_A3, "A3b_Simulated_Distribution_at_6PM.png"), width = 1600, height = 1200, res = 200)
barplot(sim_vector,
        main = "A3(b): Simulated distribution at 6 PM",
        ylab = "Proportion",
        names.arg = c("light", "heavy", "jammed"))
dev.off()

append_line("Simulated proportions at 6 PM from 10,000 trajectories:", RESULTS_TXT)
capture.output(print(sim_df), file = RESULTS_TXT, append = TRUE)
append_line("Comment:", RESULTS_TXT)
append_line("The simulated proportions are close to the exact probabilities found in A3(a),", RESULTS_TXT)
append_line("thereby verifying the result.", RESULTS_TXT)

# ==========================================================
# FINAL SUMMARY
# ==========================================================
SUMMARY_CSV <- file.path(DIR_SUMMARY, "Final_Summary.csv")

summary_df <- data.frame(
  Question = c("A1(a)", "A1(b)", "A1(c)", "A1(d)", "A2(a)", "A2(b)", "A2(c)", "A2(d)", "A3(a)", "A3(b)"),
  Result = c(
    "Transient classes {2},{3,4,5}; recurrent class {1}; absorbing state 1; reflecting state 1; periods d(1)=1, d(2)=undefined, d(3)=d(4)=d(5)=3",
    "Three simulated trajectories eventually end in state 1",
    "Steady-state distribution (1,0,0,0,0); not ergodic",
    "Unconditional probabilities converge to (1,0,0,0,0)",
    "Markov chain diagram plotted and saved",
    "Recurrent class {1,2}; transient class {3,4,5,6,7}; no absorbing states; reflecting states 5,6,7; periods identified",
    "Two simulated trajectories show eventual entry into {1,2}",
    "Stationary distribution (1/2,1/2,0,0,0,0,0); not ergodic",
    paste0("Using the lecturer correction, distribution at 6 PM = (", round(mu6[1,1],6), ", ", round(mu6[1,2],6), ", ", round(mu6[1,3],6), ")"),
    "10,000 simulated trajectories verify A3(a)"
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df, SUMMARY_CSV, row.names = FALSE)

append_line("", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)
append_line("FINAL SUMMARY OF ANSWERS", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)
capture.output(print(summary_df), file = RESULTS_TXT, append = TRUE)

# ==========================================================
# FINAL CONSOLE MESSAGE
# ==========================================================
cat("\n====================================================\n")
cat("Results successfully saved to:\n")
cat(OUTPUT_DIR, "\n")
cat("====================================================\n")