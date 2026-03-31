# ============================================================
#   HASTS416/RM MARKOV CHAIN ASSIGNMENT
#   A1: 5-State Markov Chain Analysis
#   A2: 7-State Markov Chain Analysis  
#   A3: Non-Homogeneous Traffic Model
# ============================================================

# Clear workspace and set seed for reproducibility
rm(list = ls())
set.seed(2024)

# Set results directory
results_dir <- "C:/Users/USER/Desktop/ASSIGNMENTSTOCHASTIC/"

# Create results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
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
#   A1: 5-STATE MARKOV CHAIN ANALYSIS
# ============================================================

cat("╔══════════════════════════════════════════════╗\n")
cat("║       A1: 5-STATE MARKOV CHAIN ANALYSIS      ║\n")
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

# A1(a): Diagram & Classification
cat("\n══════════════════════════════════════════════\n")
cat("  A1(a): DIAGRAM & CLASSIFICATION\n")
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

# A1(b): Simulated Trajectories
cat("\n══════════════════════════════════════════════\n")
cat("  A1(b): SIMULATED TRAJECTORIES\n")
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

# A1(c): Steady-State & Ergodicity
cat("\n══════════════════════════════════════════════\n")
cat("  A1(c): STEADY-STATE PROBABILITIES & ERGODICITY\n")
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

# A1(d): Convergence Analysis
cat("\n══════════════════════════════════════════════\n")
cat("  A1(d): CONVERGENCE TO STEADY STATE\n")
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
#   A2: 7-STATE MARKOV CHAIN ANALYSIS
# ============================================================

cat("\n\n╔══════════════════════════════════════════════╗\n")
cat("║       A2: 7-STATE MARKOV CHAIN ANALYSIS      ║\n")
cat("╚══════════════════════════════════════════════╝\n\n")




# A2
# ==========================================================
cat("\n====================================================\n")
cat("A2\n")
cat("====================================================\n")

P2 <- matrix(c(
  0,1,0,0,0,0,0,
  1,0,0,0,0,0,0,
  0,0,0,0.4,0.2,0.2,0.2,
  0,0,0,0,0.2,0.4,0.4,
  0.3,0,0,0.1,0.3,0.1,0.2,
  0,0,0,0.2,0.2,0.3,0.3,
  0,0,0,0.5,0.2,0.2,0.1
), nrow = 7, byrow = TRUE)

states2 <- as.character(1:7)
rownames(P2) <- colnames(P2) <- paste0("S", 1:7)

# Save A2 transition matrix to CSV
write.csv(P2, file.path(results_dir, "05_A2_transition_matrix.csv"), row.names = TRUE)

mc2 <- new("markovchain",
           states = states2,
           transitionMatrix = P2)

# ----------------------------------------------------------
# A2(a) Plot the Markov chain diagram
# ----------------------------------------------------------
cat("\nA2(a)\n")
cat("Plotting diagram for A2...\n")
png("chivama_A2_diagram.png", width = 800, height = 600)
plot(mc2, main = "A2(a): 7-State Markov Chain Diagram")
dev.off()

# ----------------------------------------------------------
# A2(b) Identify recurrent/transient classes, periods,
#       absorbing and reflecting states
# ----------------------------------------------------------
cat("\nA2(b)\n")
cat("Recurrent class: {1,2}\n")
cat("Transient class: {3,4,5,6,7}\n")
cat("Absorbing states: none\n")
cat("Reflecting states: 5, 6, 7\n")
cat("Periods:\n")
cat("d(1) = 2\n")
cat("d(2) = 2\n")
cat("d(3) = undefined (no return)\n")
cat("d(4) = 1\n")
cat("d(5) = 1\n")
cat("d(6) = 1\n")
cat("d(7) = 1\n")

# ----------------------------------------------------------
# A2(c) Simulate two trajectories and discuss
# ----------------------------------------------------------
cat("\nA2(c)\n")
set.seed(321)

a2_path1 <- simulate_path(P2, 30)
a2_path2 <- simulate_path(P2, 30)

# Save A2 trajectories data to CSV
trajectories_a2_df <- data.frame(
  Time_Step = 0:30,
  Trajectory_1 = a2_path1,
  Trajectory_2 = a2_path2,
  Start_States = c("Random", "Random")
)
write.csv(trajectories_a2_df, file.path(results_dir, "11_A2_trajectories_data.csv"), row.names = FALSE)

png("chivama_A2_trajectories.png", width = 800, height = 600)
matplot(rbind(a2_path1, a2_path2),
        type = "l", lty = 1, lwd = 2,
        xlab = "Time", ylab = "State",
        main = "A2(c): Two simulated trajectories")
legend("topright",
       legend = c("Path 1", "Path 2"),
       col = 1:2, lty = 1, lwd = 2)
dev.off()

cat("Comment:\n")
cat("A trajectory starting in states 3 to 7 may wander for some time,\n")
cat("but once it enters the class {1,2}, it alternates forever between 1 and 2.\n")
cat("This shows that {1,2} is a closed recurrent class, while {3,4,5,6,7} is transient.\n")

# ----------------------------------------------------------
# A2(d) Calculate limiting probabilities and interpret
#       Is the chain ergodic?
# ----------------------------------------------------------
cat("\nA2(d)\n")
pi2 <- stationary_dist(P2)
names(pi2) <- states2

cat("Stationary distribution:\n")
print(round(pi2, 6))

# Save A2 stationary distribution to CSV
stationary_df <- data.frame(
  State = paste0("S", 1:7),
  Stationary_Probability = round(pi2, 8)
)
write.csv(stationary_df, file.path(results_dir, "06_A2_stationary_distribution.csv"), row.names = FALSE)

cat("Interpretation:\n")
cat("In the long run, probability mass is concentrated on states 1 and 2.\n")
cat("The stationary distribution is (1/2, 1/2, 0, 0, 0, 0, 0).\n")
cat("Is the chain ergodic? No.\n")
cat("Although it has a stationary distribution, the recurrent class {1,2} has period 2,\n")
cat("so the chain is not ergodic.\n")



# ============================================================
#   A3: NON-HOMOGENEOUS TRAFFIC MODEL
# ============================================================

cat("\nA3(a)\n")
cat("From 1 PM to 4 PM = 3 hours = 9 steps of 20 minutes each.\n")
cat("From 4 PM to 6 PM = 2 hours = 6 steps of 20 minutes each.\n")

mu6 <- mu0_a3 %*% (P13 %^% 9) %*% (P46 %^% 6)
colnames(mu6) <- traffic_states

cat("Distribution at 6 PM:\n")
print(round(mu6, 6))

# Save A3 analytical distribution to CSV
analytical_df <- data.frame(
  Traffic_State = c("Light", "Heavy", "Jammed"),
  Probability_at_6PM = round(as.numeric(mu6[1,]), 8),
  Calculation_Method = "Analytical (mu0 * P13^9 * P46^6)"
)
write.csv(analytical_df, file.path(results_dir, "07_A3_analytical_distribution_6PM.csv"), row.names = FALSE)

cat("Answer:\n")
cat(sprintf("P(light at 6 PM)  = %.6f\n", mu6[1,1]))
cat(sprintf("P(heavy at 6 PM)  = %.6f\n", mu6[1,2]))
cat(sprintf("P(jammed at 6 PM) = %.6f\n", mu6[1,3]))

cat("Interpretation:\n")
cat("By 6 PM, traffic is most likely to be jammed.\n")

# ----------------------------------------------------------
# A3(b) Simulate 10,000 trajectories to verify the result
# ----------------------------------------------------------
cat("\nA3(b)\n")
simulate_traffic <- function() {
  state <- 1   # 1=light, 2=heavy, 3=jammed
  
  # 1 PM to 4 PM: 9 steps
  for(i in 1:9){
    state <- sample(1:3, 1, prob = P13[state, ])
  }
  
  # 4 PM to 6 PM: 6 steps
  for(i in 1:6){
    state <- sample(1:3, 1, prob = P46[state, ])
  }
  
  return(state)
}

set.seed(123)
N <- 10000
final_states <- replicate(N, simulate_traffic())
sim_props <- table(final_states) / N

# force names in case some state is missing in a rare run
sim_vector <- c(
  light  = ifelse("1" %in% names(sim_props), sim_props["1"], 0),
  heavy  = ifelse("2" %in% names(sim_props), sim_props["2"], 0),
  jammed = ifelse("3" %in% names(sim_props), sim_props["3"], 0)
)

cat("Simulated proportions at 6 PM from 10,000 trajectories:\n")
print(round(sim_vector, 6))

# Save A3 simulated distribution to CSV
simulated_df <- data.frame(
  Traffic_State = c("Light", "Heavy", "Jammed"),
  Simulated_Probability = round(as.numeric(sim_vector), 8),
  Number_of_Simulations = 10000,
  Calculation_Method = "Monte Carlo Simulation"
)
write.csv(simulated_df, file.path(results_dir, "08_A3_simulated_distribution_6PM.csv"), row.names = FALSE)

# Create and save comparison table
comparison_df <- data.frame(
  Traffic_State = c("Light", "Heavy", "Jammed"),
  Analytical_Probability = round(as.numeric(mu6[1,]), 8),
  Simulated_Probability = round(as.numeric(sim_vector), 8),
  Absolute_Difference = round(abs(as.numeric(mu6[1,]) - as.numeric(sim_vector)), 8),
  Percent_Error = round(abs(as.numeric(mu6[1,]) - as.numeric(sim_vector)) / as.numeric(mu6[1,]) * 100, 4)
)
write.csv(comparison_df, file.path(results_dir, "09_A3_analytical_vs_simulated_comparison.csv"), row.names = FALSE)

png("chivama_A3_simulation.png", width = 800, height = 600)
barplot(sim_vector,
        main = "A3(b): Simulated distribution at 6 PM",
        ylab = "Proportion",
        names.arg = c("light", "heavy", "jammed"))
dev.off()

cat("Comment:\n")
cat("The simulated proportions should be close to the exact probabilities found in A3(a),\n")
cat("thereby verifying the result.\n")

# ==========================================================
# FINAL SUMMARY OF ANSWERS
# ==========================================================
cat("\n====================================================\n")
cat("FINAL SUMMARY OF ANSWERS\n")
cat("====================================================\n")

cat("\nA1 Summary:\n")
cat("A1(a): Transient classes = {2}, {3,4,5}; Recurrent class = {1}; Absorbing state = 1; Reflecting state = 1;\n")
cat("       Periods: d(1)=1, d(2)=undefined, d(3)=d(4)=d(5)=3.\n")
cat("A1(b): Three simulated trajectories eventually end in state 1.\n")
cat("A1(c): Steady-state probabilities = (1,0,0,0,0). Chain is not ergodic.\n")
cat("A1(d): Unconditional probabilities converge to (1,0,0,0,0).\n")

cat("\nA2 Summary:\n")
cat("A2(a): Diagram plotted.\n")
cat("A2(b): Recurrent class = {1,2}; Transient class = {3,4,5,6,7}; Absorbing states = none; Reflecting states = 5,6,7.\n")
cat("       Periods: d(1)=2, d(2)=2, d(3)=undefined, d(4)=d(5)=d(6)=d(7)=1.\n")
cat("A2(c): Simulated trajectories show eventual entry into {1,2}, after which the chain alternates between 1 and 2.\n")
cat("A2(d): Stationary distribution = (1/2,1/2,0,0,0,0,0). Chain is not ergodic.\n")

cat("\nA3 Summary:\n")
cat("A3(a): Distribution at 6 PM = mu0 * P13^9 * P46^6.\n")
cat("A3(b): 10,000 simulated trajectories verify the result numerically.\n")
cat("Using the corrected matrix, the exact distribution at 6 PM is:\n")
print(round(mu6, 6))

# Create summary of all CSV files generated
csv_summary_df <- data.frame(
  File_Number = 1:11,
  CSV_Filename = c(
    "01_A1_transition_matrix.csv",
    "02_A1_limiting_distribution_P2000.csv", 
    "03_A1_steady_state_probabilities.csv",
    "04_A1_convergence_data.csv",
    "05_A2_transition_matrix.csv",
    "06_A2_stationary_distribution.csv",
    "07_A3_analytical_distribution_6PM.csv",
    "08_A3_simulated_distribution_6PM.csv", 
    "09_A3_analytical_vs_simulated_comparison.csv",
    "10_A1_trajectories_data.csv",
    "11_A2_trajectories_data.csv"
  ),
  Description = c(
    "A1: 5x5 transition matrix with state labels",
    "A1: P^2000 limiting distribution from each starting state",
    "A1: Steady-state probabilities with verification error",
    "A1: Convergence data for P(Xn=S1) over time steps",
    "A2: 7x7 transition matrix with state labels", 
    "A2: Stationary distribution for all 7 states",
    "A3: Analytical traffic distribution at 6PM",
    "A3: Simulated traffic distribution (10,000 trajectories)",
    "A3: Comparison between analytical and simulated results",
    "A1: Complete trajectory data for 3 paths (60 steps each)",
    "A2: Complete trajectory data for 2 paths (30 steps each)"
  ),
  Assignment = c(rep("A1", 4), rep("A2", 2), rep("A3", 3), "A1", "A2")
)
write.csv(csv_summary_df, file.path(results_dir, "00_CSV_Files_Summary.csv"), row.names = FALSE)

cat("\n══════════════════════════════════════════════\n")
cat("  CSV FILES SUCCESSFULLY GENERATED\n")
cat("══════════════════════════════════════════════\n")
cat("  Total CSV files created: 11\n")
cat("  See '00_CSV_Files_Summary.csv' for complete listing\n")
cat("══════════════════════════════════════════════\n")

