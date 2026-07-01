library(qgraph)

# ==========================================================
# 1. GENERATE THE ADJACENCY MATRIX (26 Nodes Total)
# ==========================================================
set.seed(123)
net_mat <- matrix(rnorm(26^2), 26, 26)
net_mat[lower.tri(net_mat)] <- t(net_mat)[lower.tri(net_mat)]  # simétrica
diag(net_mat) <- 0

# ==========================================================
# 2. CREATE THE BASE RING LAYOUT
# ==========================================================
inner_n <- 16
outer_n <- 10   

inner_layout <- matrix(0, nrow = inner_n, ncol = 2)

# Pushed the outer protein circle way out to 5.5 to allow a massive center expansion
outer_layout <- t(sapply(1:outer_n, function(i) {
  angle <- 2 * pi * (i - 1) / outer_n
  c(5.5 * cos(angle), 5.5 * sin(angle))
}))

layout_fixed <- rbind(inner_layout, outer_layout)

# ==========================================================
# 3. RADIAL EXPANSION & TARGETED CLUSTER DE-CLUMPING
# ==========================================================
coord_matrix <- matrix(c(
  -0.18447951,  0.12208006,  # 1
  0.16672297,  0.02699783,  # 2
  0.23182076, -0.63120213,  # 3
  -0.57402036, -0.20762068,  # 4
  -0.19963610, -0.74039055,  # 5
  0.04447803, -0.34755373,  # 6
  0.04541781, -0.09140243,  # 7
  0.03974421,  0.21725014,  # 8
  -0.18906295, -0.28760982,  # 9
  0.40579570,  0.01448752,  # 10
  -0.55537431,  0.16233900,  # 11
  -0.25683438,  0.32406681,  # 12
  -0.31776237, -0.08367255,  # 13
  -0.01704467, -0.56865497,  # 14
  -0.75447747,  0.11021163,  # 15
  -0.03996360,  0.63326899   # 16
), ncol = 2, byrow = TRUE)

# Center coordinates relative to their own spatial midpoint
center_x <- mean(coord_matrix[, 1])
center_y <- mean(coord_matrix[, 2])
centered_coords <- cbind(coord_matrix[, 1] - center_x, coord_matrix[, 2] - center_y)

# Aggressive global separation factor for large node sizes
expansion_factor <- 5.0  
expanded_center  <- centered_coords * expansion_factor

# ----------------------------------------------------------
# SURGICAL SEPARATION FOR TIGHT OVERLAPS (Tweak if needed)
# ----------------------------------------------------------

# Dispersing the Left Central Pile (1, 11, 12, 13)
expanded_center[11, 1] <- expanded_center[11, 1] - 0.6  # 11 Left
expanded_center[11, 2] <- expanded_center[11, 2] + 0.4  # 11 Up
expanded_center[12, 2] <- expanded_center[12, 2] + 0.7  # 12 Up
expanded_center[1, 2]  <- expanded_center[1, 2]  - 0.4  # 1 Down
expanded_center[13, 1] <- expanded_center[13, 1] - 0.5  # 13 Left
expanded_center[13, 2] <- expanded_center[13, 2] - 0.5  # 13 Down

# Dispersing the Right Central Pile (2, 6, 7, 8, 10)
expanded_center[8, 2]  <- expanded_center[8, 2]  + 0.6  # 8 Up
expanded_center[10, 1] <- expanded_center[10, 1] + 0.8  # 10 Right
expanded_center[2, 1]  <- expanded_center[2, 1]  + 0.5  # 2 Right
expanded_center[7, 1]  <- expanded_center[7, 1]  - 0.3  # 7 Nudge Left
expanded_center[6, 2]  <- expanded_center[6, 2]  - 0.5  # 6 Down

# Fine-tuning other close neighbors
expanded_center[15, 1] <- expanded_center[15, 1] - 0.8  # Move 15 further left
expanded_center[16, 2] <- expanded_center[16, 2] + 0.5  # Move 16 further up
expanded_center[14, 2] <- expanded_center[14, 2] - 0.4  # Move 14 further down

# Center the whole group visually
expanded_center[, 1] <- expanded_center[, 1] + 0.1
expanded_center[, 2] <- expanded_center[, 2] + 0.1

layout_fixed[1:16, ] <- expanded_center

# ==========================================================
# 4. PLOT THE GRAPH
# ==========================================================
qgraph(net_mat, layout = layout_fixed, labels = 1:26, vsize = 8)













library(qgraph)

# ==========================================================
# 1. GENERATE THE NEW ADJACENCY MATRIX (23 Nodes Total)
# ==========================================================
set.seed(123)
net_mat_2 <- matrix(rnorm(23^2), 23, 23)
net_mat_2[lower.tri(net_mat_2)] <- t(net_mat_2)[lower.tri(net_mat_2)]  # simétrica
diag(net_mat_2) <- 0

# ==========================================================
# 2. RUN THE EXACT SAME CENTER DE-CLUMPING LOGIC
# ==========================================================
coord_matrix <- matrix(c(
  -0.18447951,  0.12208006,  # 1
  0.16672297,  0.02699783,  # 2
  0.23182076, -0.63120213,  # 3
  -0.57402036, -0.20762068,  # 4
  -0.19963610, -0.74039055,  # 5
  0.04447803, -0.34755373,  # 6
  0.04541781, -0.09140243,  # 7
  0.03974421,  0.21725014,  # 8
  -0.18906295, -0.28760982,  # 9
  0.40579570,  0.01448752,  # 10
  -0.55537431,  0.16233900,  # 11
  -0.25683438,  0.32406681,  # 12
  -0.31776237, -0.08367255,  # 13
  -0.01704467, -0.56865497,  # 14
  -0.75447747,  0.11021163,  # 15
  -0.03996360,  0.63326899   # 16
), ncol = 2, byrow = TRUE)

center_x <- mean(coord_matrix[, 1])
center_y <- mean(coord_matrix[, 2])
centered_coords <- cbind(coord_matrix[, 1] - center_x, coord_matrix[, 2] - center_y)

expansion_factor <- 5.0  
expanded_center  <- centered_coords * expansion_factor

# Exact same manual adjustments to completely prevent overlaps
expanded_center[11, 1] <- expanded_center[11, 1] - 0.6
expanded_center[11, 2] <- expanded_center[11, 2] + 0.4
expanded_center[12, 2] <- expanded_center[12, 2] + 0.7
expanded_center[1, 2]  <- expanded_center[1, 2]  - 0.4
expanded_center[13, 1] <- expanded_center[13, 1] - 0.5
expanded_center[13, 2] <- expanded_center[13, 2] - 0.5
expanded_center[8, 2]  <- expanded_center[8, 2]  + 0.6
expanded_center[10, 1] <- expanded_center[10, 1] + 0.8
expanded_center[2, 1]  <- expanded_center[2, 1]  + 0.5
expanded_center[7, 1]  <- expanded_center[7, 1]  - 0.3
expanded_center[6, 2]  <- expanded_center[6, 2]  - 0.5
expanded_center[15, 1] <- expanded_center[15, 1] - 0.8
expanded_center[16, 2] <- expanded_center[16, 2] + 0.5
expanded_center[14, 2] <- expanded_center[14, 2] - 0.4
expanded_center[, 1]   <- expanded_center[, 1]   + 0.1
expanded_center[, 2]   <- expanded_center[, 2]   + 0.1

# ==========================================================
# 3. CREATE LAYOUT_FIXED_2 (16 Center Nodes + 7 Outer Nodes)
# ==========================================================
inner_n_2 <- 16
outer_n_2 <- 7  # 7 nodes on the outside instead of 10

inner_layout_2 <- matrix(0, nrow = inner_n_2, ncol = 2)

# Generate a symmetrical outer ring for exactly 7 nodes
outer_layout_2 <- t(sapply(1:outer_n_2, function(i) {
  angle <- 2 * pi * (i - 1) / outer_n_2
  c(5.5 * cos(angle), 5.5 * sin(angle))
}))

# Combine them into layout_fixed_2 (Total 23 rows)
layout_fixed_2 <- rbind(inner_layout_2, outer_layout_2)

# Overwrite the first 16 rows with our exact, polished center coordinates
layout_fixed_2[1:16, ] <- expanded_center

# ==========================================================
# 4. PLOT THE NEW GRAPH
# ==========================================================
qgraph(net_mat_2, layout = layout_fixed_2, labels = 1:23, vsize = 8)

