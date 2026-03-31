# ============================================================
#   HASTS416/RM MARKOV CHAIN QUESTIONS
#   A1: 5-State Markov Chain Analysis
#   A2: 7-State Markov Chain Analysis
#   A3: Non-Homogeneous Traffic Model
#   BASE R ONLY — no external packages required
# ============================================================

rm(list = ls())
set.seed(2024)

# ================================================================
#  OUTPUT DIRECTORY  (change this path to suit your machine)
# ================================================================
results_dir <- "STOCHASTIC_RESULTS"   # relative to working directory
for (d in c(results_dir,
            file.path(results_dir, "A1_Results"),
            file.path(results_dir, "A2_Results"),
            file.path(results_dir, "A3_Results"),
            file.path(results_dir, "Summary"))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

DIR_A1      <- file.path(results_dir, "A1_Results")
DIR_A2      <- file.path(results_dir, "A2_Results")
DIR_A3      <- file.path(results_dir, "A3_Results")
DIR_SUMMARY <- file.path(results_dir, "Summary")
RESULTS_TXT <- file.path(results_dir, "Results_Summary.txt")

# Initialise (overwrite) the text results file
cat("STOCHASTIC MARKOV CHAIN ANALYSIS — RESULTS SUMMARY\n",
    file = RESULTS_TXT, append = FALSE)
cat(paste0("Generated: ", Sys.time(), "\n\n"),
    file = RESULTS_TXT, append = TRUE)

# Convenience wrapper
append_line <- function(text, file = RESULTS_TXT) {
  cat(text, "\n", file = file, append = TRUE)
}

# ================================================================
#  SHARED UTILITY FUNCTIONS
# ================================================================

# Matrix power (replaces the %^% operator from markovchain)
matrix_power <- function(P, n) {
  result <- diag(nrow(P))
  if (n == 0) return(result)
  for (i in seq_len(n)) result <- result %*% P
  result
}

# Simulate a single Markov chain trajectory
simulate_mc <- function(P, start, n_steps = 60) {
  s    <- integer(n_steps + 1)
  s[1] <- start
  for (t in seq_len(n_steps))
    s[t + 1] <- sample.int(nrow(P), 1L, prob = P[s[t], ])
  s
}

# Simulate one path from a randomly chosen start state
simulate_path <- function(P, n_steps) {
  start <- sample.int(nrow(P), 1L)
  simulate_mc(P, start, n_steps)
}

# Stationary distribution via linear system  (P^T - I)π = 0, Σπ = 1
# Works correctly for any irreducible chain or recurrent class
stationary_dist <- function(P) {
  n  <- nrow(P)
  A  <- t(P) - diag(n)
  A[n, ] <- 1          # replace last equation with normalisation
  b  <- c(rep(0, n - 1), 1)
  tryCatch(solve(A, b), error = function(e) rep(NA, n))
}

# Compute unconditional marginal probabilities  π₀ P^t  for t = 0..n_time
compute_marginals <- function(pi0, P, n_time) {
  out      <- matrix(0, nrow = n_time + 1, ncol = nrow(P))
  out[1, ] <- pi0
  Pt       <- diag(nrow(P))
  for (t in seq_len(n_time)) {
    Pt         <- Pt %*% P
    out[t + 1, ] <- pi0 %*% Pt
  }
  out
}

# ---- Drawing helpers (base R, no packages) ----

draw_node <- function(x, y, r, fill, border, label, cex = 1.4) {
  theta <- seq(0, 2 * pi, length.out = 200)
  polygon(x + r * cos(theta), y + r * sin(theta),
          col = fill, border = border, lwd = 2.5)
  text(x, y, label, col = "white", cex = cex, font = 2)
}

draw_curved_arrow <- function(x0, y0, x1, y1, label, col,
                               curve = 0.3, lwd = 2.2,
                               lpos = 0.5, loff_x = 0, loff_y = 0,
                               tcex = 0.85) {
  dx <- x1 - x0;  dy <- y1 - y0
  cx <- (x0 + x1) / 2 - curve * dy
  cy <- (y0 + y1) / 2 + curve * dx
  tv <- seq(0, 1, length.out = 50)
  bx <- (1 - tv)^2 * x0 + 2 * (1 - tv) * tv * cx + tv^2 * x1
  by <- (1 - tv)^2 * y0 + 2 * (1 - tv) * tv * cy + tv^2 * y1
  lines(bx, by, col = col, lwd = lwd)
  arrows(bx[49], by[49], bx[50], by[50],
         length = 0.14, lwd = lwd, col = col, angle = 25)
  tm <- lpos
  lx <- (1 - tm)^2 * x0 + 2 * (1 - tm) * tm * cx + tm^2 * x1 + loff_x
  ly <- (1 - tm)^2 * y0 + 2 * (1 - tm) * tm * cy + tm^2 * y1 + loff_y
  text(lx, ly, label, cex = tcex, font = 2, col = "#1A252F")
}

draw_self_loop <- function(cx, cy, r_node = 0.38, r_loop = 0.52,
                            angle_deg = 90, col = "#555555", lwd = 2.0,
                            label = "", tcex = 0.78) {
  a       <- angle_deg * pi / 180
  ox      <- cx + r_node * cos(a)
  oy      <- cy + r_node * sin(a)
  loop_cx <- cx + (r_node + r_loop * 0.9) * cos(a)
  loop_cy <- cy + (r_node + r_loop * 0.9) * sin(a)
  phi     <- seq(0, 2 * pi, length.out = 100)
  slx     <- loop_cx + r_loop * cos(phi)
  sly     <- loop_cy + r_loop * sin(phi)
  lines(slx, sly, col = col, lwd = lwd)
  idx_tip <- which.min((slx - ox)^2 + (sly - oy)^2)[1]
  idx2    <- (idx_tip - 3) %% 100 + 1
  arrows(slx[idx2], sly[idx2], slx[idx_tip], sly[idx_tip],
         length = 0.12, lwd = lwd, col = col, angle = 22)
  text(loop_cx + r_loop * cos(a) * 0.8,
       loop_cy + r_loop * sin(a) * 0.8,
       label, cex = tcex, font = 2, col = "#1A1A2E")
}


# ================================================================
#   A1: 5-STATE MARKOV CHAIN
# ================================================================

cat("╔══════════════════════════════════════════════╗\n")
cat("║    A1: 5-STATE MARKOV CHAIN ANALYSIS         ║\n")
cat("╚══════════════════════════════════════════════╝\n\n")

P_A1 <- matrix(c(
  1.0, 0.0, 0.0, 0.0, 0.0,
  0.5, 0.0, 0.0, 0.0, 0.5,
  0.2, 0.0, 0.0, 0.0, 0.8,
  0.0, 0.0, 1.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 1.0, 0.0
), nrow = 5, byrow = TRUE)
rownames(P_A1) <- colnames(P_A1) <- paste0("S", 1:5)

cat("Transition Matrix P_A1:\n"); print(P_A1)
write.csv(P_A1, file.path(DIR_A1, "01_A1_transition_matrix.csv"))

# ------------------------------------------------------------------
# A1(a): Diagram
# ------------------------------------------------------------------
cat("\n── A1(a): Diagram & Classification ──\n")

nx_A1 <- c(-2.5,  0.0,  0.0,  2.5,  2.5)
ny_A1 <- c( 0.0,  1.7, -1.5, -1.5,  1.7)
nr_A1 <- 0.40
fill_A1  <- c("#C0392B","#2980B9","#27AE60","#E67E22","#8E44AD")
frame_A1 <- c("#922B21","#1A5276","#1E8449","#B7770D","#6C3483")

png(file.path(DIR_A1, "A1a_diagram.png"), width = 1200, height = 850, res = 130)
par(mar = c(2, 2, 4, 2), bg = "#F8F9FA")
plot(0, 0, type = "n", xlim = c(-3.4, 3.6), ylim = c(-2.4, 2.7),
     asp = 1, axes = FALSE, xlab = "", ylab = "",
     main = "A1(a): 5-State Markov Chain — Transition Diagram", cex.main = 1.3)

# S1 self-loop
theta_sl <- seq(pi * 0.6, pi * 2.4, length.out = 80)
lines(nx_A1[1] + 0.62 * cos(theta_sl), ny_A1[1] + 0.62 * sin(theta_sl),
      col = "#C0392B", lwd = 2.3)
arrows(nx_A1[1] + 0.62 * cos(theta_sl[79]),
       ny_A1[1] + 0.62 * sin(theta_sl[79]),
       nx_A1[1] + 0.62 * cos(theta_sl[80]),
       ny_A1[1] + 0.62 * sin(theta_sl[80]),
       length = 0.14, lwd = 2.3, col = "#C0392B", angle = 25)
text(nx_A1[1] - 0.85, ny_A1[1] + 0.67, "1.0", cex = 0.88, font = 2, col = "#1A252F")

draw_curved_arrow(nx_A1[2] - nr_A1*0.65, ny_A1[2] - nr_A1*0.65,
                  nx_A1[1] + nr_A1*0.65, ny_A1[1] + nr_A1*0.65,
                  "0.5", col="#2980B9", curve=0.30, loff_x=-0.12, loff_y=0.12)
draw_curved_arrow(nx_A1[2]+nr_A1, ny_A1[2], nx_A1[5]-nr_A1, ny_A1[5],
                  "0.5", col="#2980B9", curve=0.0, loff_y=0.18)
draw_curved_arrow(nx_A1[3]-nr_A1*0.65, ny_A1[3]+nr_A1*0.65,
                  nx_A1[1]+nr_A1*0.5, ny_A1[1]-nr_A1*0.5,
                  "0.2", col="#27AE60", curve=-0.30, loff_x=-0.16, loff_y=-0.12)
draw_curved_arrow(nx_A1[3]+nr_A1*0.6, ny_A1[3]+nr_A1*0.6,
                  nx_A1[5]-nr_A1*0.6, ny_A1[5]-nr_A1*0.6,
                  "0.8", col="#27AE60", curve=-0.25, loff_x=0.18)
draw_curved_arrow(nx_A1[4]-nr_A1, ny_A1[4], nx_A1[3]+nr_A1, ny_A1[3],
                  "1.0", col="#E67E22", curve=0.0, loff_y=-0.20)
draw_curved_arrow(nx_A1[5], ny_A1[5]-nr_A1, nx_A1[4], ny_A1[4]+nr_A1,
                  "1.0", col="#8E44AD", curve=0.0, loff_x=0.22)

for (i in 1:5) draw_node(nx_A1[i], ny_A1[i], nr_A1,
                          fill_A1[i], frame_A1[i], paste0("S", i))
legend(-3.4, -1.6,
       legend = c("S1: Absorbing / Recurrent  [d=1]",
                  "S2: Transient, isolated  [d=1]",
                  "S3,S4,S5: Transient, period-3 cycle"),
       fill = c("#C0392B","#2980B9","#27AE60"),
       border = NA, bty = "n", cex = 0.88)
dev.off()
cat("  A1a_diagram.png written.\n")

cat("\n▶ COMMUNICATING CLASSES\n")
cat("   C1 = {S1}         — Closed, Recurrent (Absorbing)\n")
cat("   C2 = {S2}         — Open, Transient (own class; no return)\n")
cat("   C3 = {S3, S4, S5} — Open, Transient (communicate; escape to S1)\n")
cat("▶ ABSORBING: S1 only  |  REFLECTING: None\n")
cat("▶ PERIODS: S1=1, S2=1, S3=3, S4=3, S5=3\n")
cat("▶ NOT ergodic (reducible; S3–S5 periodic d=3)\n")

# ------------------------------------------------------------------
# A1(b): Three Trajectories
# ------------------------------------------------------------------
cat("\n── A1(b): Three Simulated Trajectories ──\n")

n_steps_A1 <- 60
starts_A1  <- sample(1:5, 3)
cat("Start states:", paste0("S", starts_A1), "\n")
trajs_A1 <- lapply(starts_A1, simulate_mc, P = P_A1, n_steps = n_steps_A1)

cols3_A1 <- c("#E74C3C","#27AE60","#2980B9")

png(file.path(DIR_A1, "A1b_trajectories.png"), width = 1200, height = 580, res = 130)
par(mar = c(4.5, 4.8, 3.5, 1.5), bg = "#F8F9FA")
plot(0:n_steps_A1, trajs_A1[[1]], type = "s", col = cols3_A1[1], lwd = 2.2,
     ylim = c(0.6, 5.5), xlab = "Time step n", ylab = "State",
     main = "A1(b): Three Simulated Trajectories — 5-State Chain",
     yaxt = "n", bty = "l", las = 1, cex.main = 1.2)
axis(2, at = 1:5, labels = paste0("S", 1:5), las = 1)
abline(h = 1, lty = 2, col = "#C0392B", lwd = 1.2)
text(n_steps_A1 * 0.68, 1.32, "Absorbing barrier S1", col = "#C0392B", cex = 0.85)
for (k in 2:3) lines(0:n_steps_A1, trajs_A1[[k]], type = "s",
                     col = cols3_A1[k], lwd = 2.2)
for (k in 1:3) points(0, starts_A1[k], pch = 19, col = cols3_A1[k], cex = 1.8)
legend("topright",
       legend = paste0("Trajectory ", 1:3, "  (start S", starts_A1, ")"),
       col = cols3_A1, lwd = 2.2, bty = "n", cex = 0.95)
dev.off()

abs_A1 <- sapply(trajs_A1, function(x) {
  h <- which(x == 1)[1] - 1
  if (is.na(h)) paste0(">", n_steps_A1) else as.character(h)
})
cat("Absorption steps:", paste(abs_A1, collapse = ", "), "\n")

# Save trajectory data
write.csv(data.frame(Time = 0:n_steps_A1,
                     Traj1 = trajs_A1[[1]],
                     Traj2 = trajs_A1[[2]],
                     Traj3 = trajs_A1[[3]]),
          file.path(DIR_A1, "A1b_trajectories_data.csv"), row.names = FALSE)
cat("▶ All three paths absorb into S1. Period-3 cycling visible in {S3,S4,S5}.\n")

# ------------------------------------------------------------------
# A1(c): Steady-State & Ergodicity
# ------------------------------------------------------------------
cat("\n── A1(c): Steady-State Probabilities ──\n")

Pn_A1 <- matrix_power(P_A1, 2000)
rownames(Pn_A1) <- paste0("From_S", 1:5)
colnames(Pn_A1) <- paste0("S", 1:5)
cat("\nP^2000:\n"); print(round(Pn_A1, 8))

A_mat         <- t(P_A1) - diag(5)
A_mat[5, ]    <- 1
pi_ss_A1      <- solve(A_mat, c(0, 0, 0, 0, 1))
names(pi_ss_A1) <- paste0("S", 1:5)
cat("\nSteady-state π:\n"); print(round(pi_ss_A1, 8))
cat(sprintf("Verification ||πP-π||∞ = %.2e\n",
            max(abs(pi_ss_A1 %*% P_A1 - pi_ss_A1))))

write.csv(data.frame(State = paste0("S",1:5),
                     Pi    = round(pi_ss_A1, 8)),
          file.path(DIR_A1, "A1c_steady_state.csv"), row.names = FALSE)
write.csv(round(Pn_A1, 8),
          file.path(DIR_A1, "A1c_P2000.csv"))
cat("▶ π = (1,0,0,0,0). NOT ergodic (reducible + periodic transient class).\n")

# ------------------------------------------------------------------
# A1(d): Convergence
# ------------------------------------------------------------------
cat("\n── A1(d): Convergence of Unconditional Probabilities ──\n")

n_time_A1    <- 80
marg_unif_A1 <- compute_marginals(rep(1/5, 5), P_A1, n_time_A1)
marg_S3_A1   <- compute_marginals(c(0,0,1,0,0), P_A1, n_time_A1)
pal5_A1 <- c("#C0392B","#2980B9","#27AE60","#E67E22","#8E44AD")
tvec_A1 <- 0:n_time_A1

png(file.path(DIR_A1, "A1d_convergence.png"), width = 1200, height = 960, res = 130)
par(mfrow = c(2,1), mar = c(4.2,5.0,3.5,8.5), bg = "#F8F9FA", oma = c(0,0,2,0))

plot(tvec_A1, marg_unif_A1[,1], type="l", col=pal5_A1[1], lwd=2.2,
     ylim=c(-0.02,1.07), xlab="", ylab="P(Xn = Si)",
     main="Uniform start  π0 = (0.2, 0.2, 0.2, 0.2, 0.2)",
     bty="l", las=1, cex.main=1.1)
abline(h=1, lty=2, col=pal5_A1[1], lwd=0.9)
for (j in 2:5) lines(tvec_A1, marg_unif_A1[,j], col=pal5_A1[j], lwd=2.2)
legend(par("usr")[2]+0.5, par("usr")[4], legend=paste0("S",1:5),
       col=pal5_A1, lwd=2.2, bty="n", cex=0.95, xpd=TRUE, y.intersp=1.3)

plot(tvec_A1, marg_S3_A1[,1], type="l", col=pal5_A1[1], lwd=2.2,
     ylim=c(-0.02,1.07), xlab="Time n", ylab="P(Xn = Si)",
     main="Start at S3  — period-3 oscillations clearly visible",
     bty="l", las=1, cex.main=1.1)
abline(h=1, lty=2, col=pal5_A1[1], lwd=0.9)
for (j in 2:5) lines(tvec_A1, marg_S3_A1[,j], col=pal5_A1[j], lwd=2.2)
legend(par("usr")[2]+0.5, par("usr")[4], legend=paste0("S",1:5),
       col=pal5_A1, lwd=2.2, bty="n", cex=0.95, xpd=TRUE, y.intersp=1.3)

mtext("A1(d): Convergence of Unconditional Probabilities to Steady State",
      outer=TRUE, side=3, cex=1.2, font=2)
dev.off()
cat("  A1d_convergence.png written.\n")

conv_tab <- data.frame(
  n = c(0,3,6,10,15,20,30,50,80),
  P_S1_uniform = round(marg_unif_A1[c(1,4,7,11,16,21,31,51,81), 1], 6),
  P_S1_startS3 = round(marg_S3_A1 [c(1,4,7,11,16,21,31,51,81), 1], 6)
)
print(conv_tab)
write.csv(conv_tab, file.path(DIR_A1, "A1d_convergence_table.csv"), row.names=FALSE)
cat("▶ Period-3 oscillations visible from S3 start; smooth S-curve from uniform start.\n")


# ================================================================
#   A2: 7-STATE MARKOV CHAIN
# ================================================================

cat("\n\n╔══════════════════════════════════════════════╗\n")
cat("║    A2: 7-STATE MARKOV CHAIN ANALYSIS         ║\n")
cat("╚══════════════════════════════════════════════╝\n\n")

append_line("====================================================")
append_line("A2: 7-STATE MARKOV CHAIN ANALYSIS")
append_line("====================================================")

P2 <- matrix(c(
  0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.4, 0.2, 0.2, 0.2,
  0.0, 0.0, 0.0, 0.0, 0.2, 0.4, 0.4,
  0.3, 0.0, 0.0, 0.1, 0.3, 0.1, 0.2,
  0.0, 0.0, 0.0, 0.2, 0.2, 0.3, 0.3,
  0.0, 0.0, 0.0, 0.5, 0.2, 0.2, 0.1
), nrow = 7, byrow = TRUE)
rownames(P2) <- colnames(P2) <- paste0("S", 1:7)

cat("Transition Matrix P2:\n"); print(P2)
cat("Row sums:", round(rowSums(P2), 10), "\n")
write.csv(P2, file.path(DIR_A2, "A2_transition_matrix.csv"))

# ------------------------------------------------------------------
# A2(a): Diagram  (fully hand-drawn — no markovchain package)
# ------------------------------------------------------------------
cat("\n── A2(a): Transition Diagram ──\n")
append_line("A2(a): Transition diagram saved as A2a_diagram.png")

nx2 <- c(-2.8, -2.8,  0.0,  2.8,  2.0,  3.6,  2.8)
ny2 <- c( 0.8, -0.8,  2.2,  1.5, -0.2,  0.2, -1.5)
nr2 <- 0.38
pal2 <- c("#C0392B","#C0392B","#7D3C98","#2980B9","#2980B9","#2980B9","#2980B9")
bdr2 <- c("#922B21","#922B21","#5B2C6F","#1A5276","#1A5276","#1A5276","#1A5276")

png(file.path(DIR_A2, "A2a_diagram.png"), width = 1400, height = 950, res = 130)
par(mar = c(2,2,4.5,2), bg = "#F7F9FC")
plot(0,0,type="n", xlim=c(-3.8,4.8), ylim=c(-2.5,3.2),
     asp=1, axes=FALSE, xlab="", ylab="",
     main="A2(a): 7-State Markov Chain — Transition Diagram", cex.main=1.35)

rect(-3.6,-0.4,-1.8,1.6, col="#FADBD8", border="#E74C3C", lty=2, lwd=1.5)
text(-2.8,1.75,"Recurrent\n{S1, S2}", col="#C0392B", cex=0.82, font=3)
rect(1.2,-2.2,4.5,2.2, col="#D6EAF8", border="#2980B9", lty=2, lwd=1.5)
text(2.85,2.38,"Transient Class  {S4,S5,S6,S7}", col="#1A5276", cex=0.82, font=3)
rect(-0.75,1.65,0.75,2.85, col="#E8DAEF", border="#7D3C98", lty=2, lwd=1.5)
text(0.0,3.08,"Transient {S3}", col="#7D3C98", cex=0.82, font=3)

# S1 <-> S2
draw_curved_arrow(nx2[1], ny2[1]-nr2, nx2[2], ny2[2]+nr2,
                  "1.0","#C0392B",curve=-0.25,loff_x=-0.22)
draw_curved_arrow(nx2[2], ny2[2]+nr2*0.5, nx2[1], ny2[1]-nr2*0.5,
                  "1.0","#C0392B",curve=-0.25,loff_x=0.22)
# S3 -> S4,S5,S6,S7
draw_curved_arrow(nx2[3]+nr2*0.5, ny2[3]-nr2*0.7,
                  nx2[4]-nr2*0.3, ny2[4]+nr2*0.8, "0.4","#7D3C98",curve=0.10,loff_x=-0.22)
draw_curved_arrow(nx2[3]+nr2*0.7, ny2[3]-nr2*0.5,
                  nx2[5]-nr2*0.5, ny2[5]+nr2*0.7, "0.2","#7D3C98",curve=-0.10,loff_y=0.18)
draw_curved_arrow(nx2[3]+nr2,     ny2[3],
                  nx2[6]-nr2*0.5, ny2[6]+nr2*0.7, "0.2","#7D3C98",curve=-0.15,loff_y=0.18)
draw_curved_arrow(nx2[3]+nr2*0.6, ny2[3]-nr2*0.8,
                  nx2[7]-nr2*0.4, ny2[7]+nr2*0.7, "0.2","#7D3C98",curve=-0.18,loff_x=0.12)
# S5 -> S1 (escape)
draw_curved_arrow(nx2[5]-nr2*0.9, ny2[5]+nr2*0.4,
                  nx2[1]+nr2*0.8, ny2[1]-nr2*0.5, "0.3","#E74C3C",
                  curve=-0.2,lwd=2.4,loff_y=-0.22,tcex=0.82)
# Within transient cluster
draw_curved_arrow(nx2[4]-nr2*0.5, ny2[4]-nr2*0.8,
                  nx2[5]+nr2*0.4, ny2[5]+nr2*0.8, "0.2","#1F618D",curve=0.2,loff_x=-0.18)
draw_curved_arrow(nx2[4]+nr2,     ny2[4],
                  nx2[6]-nr2,     ny2[6],          "0.4","#1F618D",curve=0.2,loff_y=0.18)
draw_curved_arrow(nx2[4]+nr2*0.3, ny2[4]-nr2*0.9,
                  nx2[7]-nr2*0.4, ny2[7]+nr2*0.6, "0.4","#1F618D",curve=0.2,loff_x=0.18)
draw_curved_arrow(nx2[5]+nr2*0.4, ny2[5]+nr2*0.9,
                  nx2[4]-nr2*0.5, ny2[4]-nr2*0.9, "0.1","#1F618D",curve=0.2,loff_x=0.18)
draw_self_loop(nx2[5],ny2[5],r_node=nr2,r_loop=0.42,angle_deg=200,
               col="#1F618D",lwd=1.8,label="0.3")
draw_curved_arrow(nx2[5]+nr2*0.8, ny2[5]+nr2*0.5,
                  nx2[6]-nr2*0.8, ny2[6]-nr2*0.5, "0.1","#1F618D",curve=0.2,loff_y=0.18)
draw_curved_arrow(nx2[5]+nr2*0.3, ny2[5]-nr2*0.9,
                  nx2[7]-nr2*0.5, ny2[7]+nr2*0.5, "0.2","#1F618D",curve=0.2,loff_x=0.18)
draw_curved_arrow(nx2[6]-nr2,     ny2[6],
                  nx2[4]+nr2,     ny2[4],          "0.2","#1F618D",curve=0.2,loff_y=-0.18)
draw_curved_arrow(nx2[6]-nr2*0.8, ny2[6]-nr2*0.5,
                  nx2[5]+nr2*0.8, ny2[5]+nr2*0.5, "0.2","#1F618D",curve=0.2,loff_y=-0.18)
draw_self_loop(nx2[6],ny2[6],r_node=nr2,r_loop=0.42,angle_deg=20,
               col="#1F618D",lwd=1.8,label="0.3")
draw_curved_arrow(nx2[6]-nr2*0.3, ny2[6]-nr2*0.9,
                  nx2[7]+nr2*0.5, ny2[7]+nr2*0.5, "0.3","#1F618D",curve=0.2,loff_x=0.18)
draw_curved_arrow(nx2[7]-nr2*0.5, ny2[7]+nr2*0.8,
                  nx2[4]+nr2*0.3, ny2[4]-nr2*0.9, "0.5","#1F618D",curve=0.2,loff_x=-0.18)
draw_curved_arrow(nx2[7]-nr2*0.5, ny2[7]+nr2*0.6,
                  nx2[5]+nr2*0.3, ny2[5]-nr2*0.9, "0.2","#1F618D",curve=-0.2,loff_x=-0.18)
draw_curved_arrow(nx2[7]+nr2*0.5, ny2[7]+nr2*0.5,
                  nx2[6]-nr2*0.3, ny2[6]-nr2*0.9, "0.2","#1F618D",curve=0.2,loff_x=0.18)
draw_self_loop(nx2[7],ny2[7],r_node=nr2,r_loop=0.42,angle_deg=270,
               col="#1F618D",lwd=1.8,label="0.1")

for (i in 1:7) draw_node(nx2[i],ny2[i],nr2,pal2[i],bdr2[i],paste0("S",i))

legend(-3.75,-1.45,
       legend=c("Recurrent {S1,S2} — period 2, reflecting",
                "Transient {S4,S5,S6,S7} — aperiodic",
                "Transient {S3} — isolated",
                "Escape to recurrent class (S5->S1, p=0.3)"),
       col=c("#C0392B","#2980B9","#7D3C98","#E74C3C"),
       lwd=c(3,3,3,2.5), bty="n", cex=0.82, y.intersp=1.3)
dev.off()
cat("  A2a_diagram.png written.\n")

# ------------------------------------------------------------------
# A2(b): Classification
# ------------------------------------------------------------------
cat("\n── A2(b): Classification ──\n")
append_line("A2(b): Classification")
append_line("Recurrent class: {S1, S2}  — closed, period d=2, reflecting pair")
append_line("Transient class: {S3,S4,S5,S6,S7}")
append_line("Absorbing states: NONE")
append_line("Reflecting states: S1 and S2 (mutual reflectors, P(S1->S2)=P(S2->S1)=1)")
append_line("Periods: S1=2, S2=2, S3=1, S4=1, S5=1 (self-loop), S6=1 (self-loop), S7=1 (self-loop)")

a2b_df <- data.frame(
  State  = paste0("S", 1:7),
  Class  = c("Recurrent","Recurrent","Transient","Transient",
             "Transient","Transient","Transient"),
  Period = c(2, 2, 1, 1, 1, 1, 1),
  Notes  = c("Reflecting pair","Reflecting pair","Isolated singleton",
             "Aperiodic","Self-loop d=1","Self-loop d=1","Self-loop d=1"),
  stringsAsFactors = FALSE
)
print(a2b_df)
write.csv(a2b_df, file.path(DIR_A2, "A2b_classification.csv"), row.names=FALSE)

cat("\n▶ Absorbing states : NONE\n")
cat("▶ Reflecting states: S1 and S2 form a REFLECTING PAIR\n")
cat("▶ Recurrent        : {S1, S2}  — period d = 2\n")
cat("▶ Transient        : {S3, S4, S5, S6, S7}\n")

# ------------------------------------------------------------------
# A2(c): Two Simulated Trajectories
# ------------------------------------------------------------------
cat("\n── A2(c): Two Simulated Trajectories ──\n")
append_line("A2(c): Two simulated trajectories plotted.")

set.seed(321)
a2_path1 <- simulate_path(P2, 30)
a2_path2 <- simulate_path(P2, 30)
start1 <- a2_path1[1]; start2 <- a2_path2[1]
cat("Start states: S", start1, " and S", start2, "\n", sep="")

write.csv(data.frame(Time=0:30, Path1=a2_path1, Path2=a2_path2),
          file.path(DIR_A2, "A2c_trajectories.csv"), row.names=FALSE)

cols2_A2 <- c("#E74C3C","#2980B9")

png(file.path(DIR_A2, "A2c_trajectories.png"), width = 1300, height = 660, res = 130)
par(mar = c(4.8,5.0,4.0,1.8), bg = "#F7F9FC")
plot(0:30, a2_path1, type="s", col=cols2_A2[1], lwd=2.4,
     ylim=c(0.5,7.6), xlab="Time step n", ylab="State",
     main="A2(c): Two Simulated Trajectories — 7-State Chain",
     yaxt="n", bty="l", las=1, cex.main=1.2)
axis(2, at=1:7, labels=paste0("S",1:7), las=1)
rect(-1, 0.5, 32, 2.5, col="#FADBD8", border=NA)
text(26, 2.35, "Recurrent {S1,S2}", col="#C0392B", cex=0.80)
rect(-1, 2.5, 32, 7.6, col="#D6EAF8", border=NA)
text(26, 7.35, "Transient {S3-S7}", col="#1A5276", cex=0.80)
lines(0:30, a2_path1, type="s", col=cols2_A2[1], lwd=2.4)
lines(0:30, a2_path2, type="s", col=cols2_A2[2], lwd=2.4)
points(0, start1, pch=21, bg=cols2_A2[1], col="white", cex=2.2)
points(0, start2, pch=21, bg=cols2_A2[2], col="white", cex=2.2)

# Mark entry into recurrent class
abs1 <- which(a2_path1 <= 2)[1]
abs2 <- which(a2_path2 <= 2)[1]
if (!is.na(abs1)) points(abs1-1, a2_path1[abs1], pch=8, col=cols2_A2[1], cex=2.0, lwd=2)
if (!is.na(abs2)) points(abs2-1, a2_path2[abs2], pch=8, col=cols2_A2[2], cex=2.0, lwd=2)

legend("topright",
       legend=c(paste0("Path 1 (start S",start1,")"),
                paste0("Path 2 (start S",start2,")")),
       col=cols2_A2, lwd=2.4, pch=21,
       pt.bg=cols2_A2, pt.cex=1.5, bty="n", cex=0.92)
dev.off()
cat("  A2c_trajectories.png written.\n")

append_line("Discussion:")
append_line("Once either path enters {S1,S2} it alternates S1<->S2 forever (period-2).")
append_line("Time in transient region varies; S5 provides the escape with p=0.3 per visit.")
append_line("S3, if visited, is left after exactly one step and never revisited.")
cat("▶ Period-2 oscillation S1↔S2 clearly visible after absorption.\n")

# ------------------------------------------------------------------
# A2(d): Limiting Probabilities & Ergodicity
# ------------------------------------------------------------------
cat("\n── A2(d): Limiting Probabilities & Ergodicity ──\n")
append_line("A2(d): Limiting probabilities and ergodicity")

Pn2 <- matrix_power(P2, 5000)
cat("\nP^5000 (rows = starting state, columns = limiting prob):\n")
print(round(Pn2, 6))

# Cesaro average for time-average stationary distribution
N_ces <- 3000
Pk  <- diag(7); Psum <- matrix(0,7,7)
for (k in seq_len(N_ces)) { Pk <- Pk %*% P2; Psum <- Psum + Pk }
Pi_ces <- Psum / N_ces
cat("\nCesaro average row for start S5:\n")
print(round(Pi_ces[5,], 6))

pi2_df <- data.frame(
  State = paste0("S", 1:7),
  Pi_cesaro_from_S5 = round(Pi_ces[5,], 6)
)
write.csv(pi2_df, file.path(DIR_A2, "A2d_limiting_distribution.csv"), row.names=FALSE)
write.csv(round(Pn2,6), file.path(DIR_A2, "A2d_P5000.csv"))

cat("\n▶ Limiting distribution (Cesaro): π ≈ (0.5, 0.5, 0, 0, 0, 0, 0)\n")
cat("▶ In the long run chain spends HALF its time in S1, half in S2.\n")
cat("▶ NOT ergodic: reducible (S1 cannot reach S3–S7) AND recurrent\n")
cat("   class {S1,S2} is periodic (d=2). Pointwise limit P^n(S1,S1)\n")
cat("   alternates 1,0,1,0,… — Cesaro limit exists at 0.5.\n")

append_line("Limiting / stationary (Cesaro) probabilities:")
append_line("pi = (0.5, 0.5, 0, 0, 0, 0, 0)")
append_line("Interpretation: chain concentrates on {S1,S2} long-run, 50/50.")
append_line("Is the chain ergodic? NO.")
append_line("Reasons: (1) Reducible — S1 cannot reach S3-S7.")
append_line("         (2) Recurrent class {S1,S2} has period d=2 (periodic).")
append_line("         Cesaro-sense limit exists but strong ergodic theorem fails.")

# Convergence plot for A2
n_time_A2    <- 100
marg_unif_A2 <- compute_marginals(rep(1/7,7), P2, n_time_A2)
marg_S3_A2   <- compute_marginals(c(0,0,1,0,0,0,0), P2, n_time_A2)
pal7 <- c("#C0392B","#E74C3C","#7D3C98","#1A5276","#2980B9","#1F618D","#85C1E9")
tvec2 <- 0:n_time_A2

png(file.path(DIR_A2, "A2d_convergence.png"), width = 1300, height = 1000, res = 130)
par(mfrow=c(2,1), mar=c(4.2,5.2,3.5,9.0), bg="#F7F9FC", oma=c(0,0,3,0))

plot(tvec2, marg_unif_A2[,1], type="l", col=pal7[1], lwd=2.2,
     ylim=c(-0.02,0.65), xlab="", ylab="P(Xn = Si)",
     main="Uniform start  pi0 = (1/7,...,1/7)", bty="l", las=1, cex.main=1.05)
for (j in 2:7) lines(tvec2, marg_unif_A2[,j], col=pal7[j], lwd=2.2)
abline(h=0.5, lty=2, col="#C0392B", lwd=0.9)
legend(par("usr")[2]+0.5, par("usr")[4], legend=paste0("S",1:7),
       col=pal7, lwd=2.2, bty="n", cex=0.88, xpd=TRUE, y.intersp=1.3)

plot(tvec2, marg_S3_A2[,1], type="l", col=pal7[1], lwd=2.2,
     ylim=c(-0.02,1.05), xlab="Time n", ylab="P(Xn = Si)",
     main="Start at S3 — period-2 oscillations in S1/S2 visible",
     bty="l", las=1, cex.main=1.05)
for (j in 2:7) lines(tvec2, marg_S3_A2[,j], col=pal7[j], lwd=2.2)
abline(h=0.5, lty=2, col="#C0392B", lwd=0.9)
legend(par("usr")[2]+0.5, par("usr")[4], legend=paste0("S",1:7),
       col=pal7, lwd=2.2, bty="n", cex=0.88, xpd=TRUE, y.intersp=1.3)

mtext("A2(d): Convergence of Unconditional Probabilities",
      outer=TRUE, side=3, cex=1.25, font=2)
dev.off()
cat("  A2d_convergence.png written.\n")


# ================================================================
#   A3: NON-HOMOGENEOUS TRAFFIC MODEL
# ================================================================

cat("\n\n╔══════════════════════════════════════════════╗\n")
cat("║    A3: NON-HOMOGENEOUS TRAFFIC MODEL         ║\n")
cat("╚══════════════════════════════════════════════╝\n\n")

append_line("====================================================")
append_line("A3: NON-HOMOGENEOUS TRAFFIC MODEL")
append_line("====================================================")
append_line("Correction applied: P13 row-2 uses 0.4 (not 0.5).")
append_line("P13 (1 PM - 4 PM): [[0.4,0.4,0.2],[0.3,0.4,0.3],[0,0.1,0.9]]")
append_line("P46 (4 PM - 6 PM): [[0.1,0.5,0.4],[0.1,0.3,0.6],[0,0.1,0.9]]")

# 1 PM -> 4 PM matrix (corrected: row 2 uses 0.4 not 0.5)
P13 <- matrix(c(
  0.4, 0.4, 0.2,
  0.3, 0.4, 0.3,
  0.0, 0.1, 0.9
), nrow = 3, byrow = TRUE)

# 4 PM -> 6 PM matrix
P46 <- matrix(c(
  0.1, 0.5, 0.4,
  0.1, 0.3, 0.6,
  0.0, 0.1, 0.9
), nrow = 3, byrow = TRUE)

traffic_states <- c("Light", "Heavy", "Jammed")
rownames(P13) <- colnames(P13) <- traffic_states
rownames(P46) <- colnames(P46) <- traffic_states

cat("P13 (1PM-4PM, corrected):\n"); print(P13)
cat("P46 (4PM-6PM):\n"); print(P46)
cat("P13 row sums:", round(rowSums(P13),10), "\n")
cat("P46 row sums:", round(rowSums(P46),10), "\n")

write.csv(P13, file.path(DIR_A3, "A3_P13_matrix.csv"))
write.csv(P46, file.path(DIR_A3, "A3_P46_matrix.csv"))

# Initial distribution: light traffic at 1 PM
mu0_A3 <- c(1, 0, 0)

# ------------------------------------------------------------------
# A3(a): Exact distribution at 6 PM
# ------------------------------------------------------------------
cat("\n── A3(a): Distribution at 6 PM (exact) ──\n")
append_line("A3(a): Distribution at 6 PM")
append_line("1 PM to 4 PM = 3 hrs = 9 steps of 20 min each -> P13^9")
append_line("4 PM to 6 PM = 2 hrs = 6 steps of 20 min each -> P46^6")

# 9 steps of P13 then 6 steps of P46
P13_9 <- matrix_power(P13, 9)
P46_6 <- matrix_power(P46, 6)
mu6   <- mu0_A3 %*% P13_9 %*% P46_6

cat("Distribution at 6 PM:\n")
for (i in 1:3) cat(sprintf("   P(%s) = %.6f\n", traffic_states[i], mu6[i]))
cat("Most likely state at 6 PM:", traffic_states[which.max(mu6)], "\n")

mu6_df <- data.frame(State = traffic_states,
                     Probability = round(as.numeric(mu6), 6))
write.csv(mu6_df, file.path(DIR_A3, "A3a_distribution_6pm.csv"), row.names=FALSE)

append_line(sprintf("P(Light at 6PM)  = %.6f", mu6[1]))
append_line(sprintf("P(Heavy at 6PM)  = %.6f", mu6[2]))
append_line(sprintf("P(Jammed at 6PM) = %.6f", mu6[3]))
append_line(paste0("Most likely state at 6PM: ", traffic_states[which.max(mu6)]))
append_line("Interpretation: By 6PM, jammed traffic is most probable.")

# Probability path over the full day (every 20 min)
times_13 <- seq(13, 16, by = 1/3)   # 13:00 to 16:00  (9 steps)
times_46 <- seq(16, 18, by = 1/3)   # 16:00 to 18:00  (6 steps)

probs_13 <- matrix(0, nrow=10, ncol=3)
probs_46 <- matrix(0, nrow=7,  ncol=3)
probs_13[1,] <- mu0_A3
Pt_tmp <- diag(3)
for (k in 1:9) { Pt_tmp <- Pt_tmp %*% P13; probs_13[k+1,] <- mu0_A3 %*% Pt_tmp }
mu_at_4pm <- probs_13[10,]
Pt_tmp <- diag(3)
probs_46[1,] <- mu_at_4pm
for (k in 1:6) { Pt_tmp <- Pt_tmp %*% P46; probs_46[k+1,] <- mu_at_4pm %*% Pt_tmp }

all_times <- c(times_13, times_46[-1])
all_probs <- rbind(probs_13, probs_46[-1,])
colnames(all_probs) <- traffic_states

png(file.path(DIR_A3, "A3a_traffic_evolution.png"), width=1200, height=600, res=130)
par(mar=c(4.5,5.0,3.5,7.5), bg="#F7F9FC")
pal3 <- c("#27AE60","#E67E22","#C0392B")
plot(all_times, all_probs[,1], type="l", col=pal3[1], lwd=2.2,
     ylim=c(0,1), xlab="Hour of day", ylab="Probability",
     main="A3(a): Traffic State Probabilities — 1 PM to 6 PM",
     xaxt="n", bty="l", las=1, cex.main=1.2)
axis(1, at=13:18, labels=paste0(13:18,":00"))
abline(v=16, lty=2, col="#555555", lwd=1.3)
text(16.03, 0.92, "4 PM\nP changes", cex=0.80, col="#555555", adj=0)
for (j in 2:3) lines(all_times, all_probs[,j], col=pal3[j], lwd=2.2)
legend(par("usr")[2]+0.05, par("usr")[4],
       legend=traffic_states, col=pal3, lwd=2.2,
       bty="n", cex=0.92, xpd=TRUE, y.intersp=1.3)
dev.off()
cat("  A3a_traffic_evolution.png written.\n")

write.csv(data.frame(Hour=all_times, all_probs),
          file.path(DIR_A3, "A3a_full_evolution.csv"), row.names=FALSE)

# ------------------------------------------------------------------
# A3(b): Monte Carlo simulation (10,000 trajectories)
# ------------------------------------------------------------------
cat("\n── A3(b): Monte Carlo Verification (N=10,000) ──\n")
append_line("A3(b): Monte Carlo simulation with N=10,000 trajectories")

set.seed(123)
N <- 10000

simulate_traffic <- function() {
  state <- 1   # start in Light at 1 PM
  for (i in seq_len(9)) state <- sample.int(3, 1L, prob = P13[state,])
  for (i in seq_len(6)) state <- sample.int(3, 1L, prob = P46[state,])
  state
}

final_states <- replicate(N, simulate_traffic())
sim_props    <- tabulate(final_states, nbins=3) / N

cat("Simulated proportions at 6 PM:\n")
for (i in 1:3)
  cat(sprintf("   P(%s) = %.4f  (exact: %.6f)\n",
              traffic_states[i], sim_props[i], mu6[i]))

sim_df <- data.frame(
  State          = traffic_states,
  Simulated      = round(sim_props, 6),
  Exact          = round(as.numeric(mu6), 6),
  Absolute_Error = round(abs(sim_props - as.numeric(mu6)), 6)
)
print(sim_df)
write.csv(sim_df, file.path(DIR_A3, "A3b_simulation_vs_exact.csv"), row.names=FALSE)

png(file.path(DIR_A3, "A3b_simulation_barplot.png"), width=1100, height=650, res=130)
par(mar=c(4.5,5.0,3.5,1.5), bg="#F7F9FC")
bar_mat <- rbind(Simulated=sim_props, Exact=as.numeric(mu6))
bp <- barplot(bar_mat, beside=TRUE,
              col=c("#2980B9","#E74C3C"),
              names.arg=traffic_states,
              ylab="Probability", ylim=c(0, max(bar_mat)*1.25),
              main="A3(b): Simulated vs Exact Distribution at 6 PM  (N=10,000)",
              cex.main=1.2, las=1)
# Add value labels
for (i in 1:2) for (j in 1:3)
  text(bp[i,j], bar_mat[i,j]+0.012,
       sprintf("%.3f", bar_mat[i,j]), cex=0.82, font=2)
legend("topright",
       legend=c("Simulated (N=10,000)","Exact (analytical)"),
       fill=c("#2980B9","#E74C3C"), border=NA, bty="n", cex=0.95)
dev.off()
cat("  A3b_simulation_barplot.png written.\n")

append_line("Simulated proportions closely match exact values from A3(a).")
append_line("Maximum absolute error:")
append_line(sprintf("%.6f", max(abs(sim_props - as.numeric(mu6)))))
append_line("Comment: Monte Carlo verifies the analytical result. By 6 PM,")
append_line("jammed traffic is most probable, confirming peak-hour congestion.")


# ================================================================
#   FINAL SUMMARY
# ================================================================

cat("\n\n╔══════════════════════════════════════════════╗\n")
cat("║           FINAL SUMMARY OF RESULTS           ║\n")
cat("╚══════════════════════════════════════════════╝\n")

summary_df <- data.frame(
  Question = c("A1(a)","A1(b)","A1(c)","A1(d)",
               "A2(a)","A2(b)","A2(c)","A2(d)",
               "A3(a)","A3(b)"),
  Key_Result = c(
    "C1={S1} absorbing; C2={S2} transient; C3={S3,S4,S5} transient period-3; no reflecting states",
    "3 trajectories all absorb into S1; period-3 cycling visible",
    "pi=(1,0,0,0,0); NOT ergodic",
    "Convergence to pi=(1,0,0,0,0); period-3 oscillations from S3 start",
    "Diagram saved as A2a_diagram.png",
    "Recurrent {S1,S2} period-2 reflecting pair; transient {S3–S7}; no absorbing; no reflecting (other than S1,S2 pair)",
    "Both paths absorbed into {S1,S2}; S1<->S2 oscillation clearly visible",
    "Cesaro limit pi=(0.5,0.5,0,0,0,0,0); NOT ergodic (reducible + periodic)",
    sprintf("P(Light)=%.4f  P(Heavy)=%.4f  P(Jammed)=%.4f", mu6[1],mu6[2],mu6[3]),
    sprintf("MC sim (N=10000): P(Light)=%.4f P(Heavy)=%.4f P(Jammed)=%.4f",
            sim_props[1], sim_props[2], sim_props[3])
  ),
  stringsAsFactors = FALSE
)

print(summary_df)
write.csv(summary_df, file.path(DIR_SUMMARY, "Final_Summary.csv"), row.names=FALSE)

append_line("====================================================")
append_line("FINAL SUMMARY")
append_line("====================================================")
for (i in seq_len(nrow(summary_df)))
  append_line(paste0(summary_df$Question[i], ": ", summary_df$Key_Result[i]))

cat("\n====================================================\n")
cat("All results saved to folder: ", results_dir, "\n")
cat("  A1_Results/  A2_Results/  A3_Results/  Summary/\n")
cat("  Results_Summary.txt\n")
cat("====================================================\n")

