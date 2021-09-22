## ---- echo = FALSE------------------------------------------------------------
oldpar <- par();

## ---- echo = FALSE, fig.width = 7, fig.height = 6, fig.cap = "Figure 1: Network mapping loci to traits through an intermediate set of hidden layers in the mine_gmatrix function"----
par(mar = c(0.2, 0.2, 0, 0.2));
plot(x = 0, y = 0, type = "n", xlim = c(0, 1000), ylim = c(0, 1000), 
     xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "");
mtext(side = 1, text = "Hidden layers", cex = 2, col = "blue", line = -1);
mtext(side = 2, text = "Loci", cex = 2, col = "darkgreen", line = -2);
mtext(side = 4, text = "Traits", cex = 2, col = "red", line = -1);
Lx  <- rep(x = 100, times = 12);
Ly  <- seq(from = 100, to = 900, length = 12);
points(x = Lx, y = Ly, pch = 20, cex = 6, col = "darkgreen");
H1x <- rep(x = 300, times = 4);
H1y <- seq(from = 100, to = 900, length = 4);
points(x = H1x, y = H1y, pch = 15, cex = 6, col = "blue");
H2x <- rep(x = 500, times = 4);
H2y <- seq(from = 100, to = 900, length = 4);
points(x = H2x, y = H2y, pch = 15, cex = 6, col = "blue");
H3x <- rep(x = 700, times = 4);
H3y <- seq(from = 100, to = 900, length = 4);
points(x = H3x, y = H3y, pch = 15, cex = 6, col = "blue");
Tx  <- rep(x = 900, times = 4);
Ty  <- seq(from = 100, to = 900, length = 4);
points(x = Tx, y = Ty, pch = 18, cex = 8, col = "red");
for(i in 1:4){
    e1 <- seq(from = -36, to = 36, length = 4);
    arrows(x0 = Lx + 26, x1 = H1x[i] - 40, y0 = Ly, y1 = H1y[i] + e1, 
           length = 0.08);
    arrows(x0 = H1x + 36, x1 = H2x[i] - 38, y0 = H1y, y1 = H2y[i] + e1, 
           length = 0.08);
    arrows(x0 = H2x + 36, x1 = H3x[i] - 38, y0 = H2y, y1 = H3y[i] + e1, 
           length = 0.08);
    e2 <- seq(from = -12, to = 12, length = 4)
    arrows(x0 = H3x + 36, x1 = Tx[i] - 38, y0 = H3y, y1 = Ty[i] + e2, 
           length = 0.08);
}

## ---- echo = FALSE------------------------------------------------------------
gmt <- matrix(data = 0, nrow = 4, ncol = 4);
diag(gmt) <- 1;
gmt[1, 2] <-  -0.5;
gmt[2, 1] <-  -0.5;
gmt[1, 3] <-  0.2;
gmt[1, 4] <-  0.2;
gmt[3, 1] <-  0.2;
gmt[4, 1] <-  0.2;
gmt[2, 3] <-  0.2;
gmt[2, 4] <-  0.2;
gmt[3, 2] <-  0.2;
gmt[4, 2] <-  0.2;
gmt[3, 4] <- -0.5;
gmt[4, 3] <- -0.5;
gmt[3, 3] <-  1.0;
gmt[4, 4] <-  1.0;
print(gmt);

## ---- eval = FALSE------------------------------------------------------------
#  mg  <- mine_gmatrix(gmatrix = gmt, loci = 12, indivs = 2000, npsize = 12000,
#                      max_gen = 1200, sampleK = 1200, chooseK = 6, layers = 6,
#                      mu_pr = 0.05, mu_sd = 0.01, pr_cross = 0.05,
#                      term_cri = -5.3, sd_ini = 0.1, use_cor = FALSE);

## ---- eval = FALSE------------------------------------------------------------
#  data("gmatrices");

## ---- echo = FALSE------------------------------------------------------------
load("../data/gmatrices.rda");
mg_v1[[1]][11] <- -5.3;

## ---- echo = FALSE------------------------------------------------------------
print(mg_v1[[1]]);

## ---- echo = FALSE------------------------------------------------------------
print(mg_v1[[2]]);

## ---- echo = FALSE------------------------------------------------------------
print(mg_v1[[3]]);

## ---- echo = FALSE------------------------------------------------------------
print(mg_v1[[4]]);

## ---- echo = FALSE------------------------------------------------------------
print(mg_v1[[5]]);

## ---- echo = FALSE------------------------------------------------------------
print(mg_v1[[6]]);

## ---- echo = FALSE------------------------------------------------------------
print(mg_v1[[7]]);

## ---- echo = FALSE------------------------------------------------------------
print(mg_v1[[8]]);

## ---- eval = FALSE------------------------------------------------------------
#  sim <- run_farm_sim(mine_output = mg_v1, N = 1000, neutral_loci  = 1000,
#                      xdim = 64, ydim = 64, repro = "sexual", max_age = 1,
#                      selfing = TRUE, food_consume = c(2, 2),
#                      pesticide_consume = c(0, 0), food_needed_surv = 1,
#                      food_needed_repr = 2, reproduction_type = "food_based",
#                      pesticide_tolerated_surv = 0, pesticide_rotation_type = 2,
#                      crop_rotation_type = 2, min_age_reproduce = 0,
#                      farms = 24, time_steps = 40, mutation_pr = 0.001,
#                      crossover_pr = 0.1, net_mu_layers = 0,
#                      crop_rotation_time = 1, pesticide_rotation_time = 1,
#                      crop_per_cell = 1, pesticide_per_cell = 0.4,
#                      crop_number = 2, pesticide_number = 2, print_inds = FALSE,
#                      K_on_birth = 1000000, min_age_move = 0,
#                      age_food_threshold = 0, min_age_feed = 0,
#                      pesticide_start = 20, print_last = TRUE,
#                      mating_distance = 2, immigration_rate = 0, metabolism = 0,
#                      baseline_metabolism = 0, movement_bouts = 4,
#                      feed_while_moving = TRUE, pesticide_while_moving = TRUE);

## ---- eval = FALSE------------------------------------------------------------
#  sim <- run_farm_sim(mine_output = mg_v1, N = 1000, neutral_loci  = 1000,
#                      xdim = 64, ydim = 64, repro = "sexual", max_age = 4,
#                      selfing = TRUE, food_consume = c(2, 2),
#                      pesticide_consume = c(0, 0), food_needed_surv = 1,
#                      food_needed_repr = 2, reproduction_type = "food_based",
#                      pesticide_tolerated_surv = 0, pesticide_rotation_type = 2,
#                      crop_rotation_type = 2, min_age_reproduce = 0,
#                      max_age_feed = 1, farms = 24, time_steps = 40,
#                      mutation_pr = 0.001, crossover_pr = 0.1, net_mu_layers = 0,
#                      crop_rotation_time = 1, pesticide_rotation_time = 1,
#                      crop_per_cell = 1, pesticide_per_cell = 0.4,
#                      crop_number = 2, pesticide_number = 2, print_inds = FALSE,
#                      K_on_birth = 1000000, min_age_move = 0,
#                      age_food_threshold = 0, min_age_feed = 0,
#                      pesticide_start = 20, print_last = TRUE,
#                      mating_distance = 2, immigration_rate = 0, metabolism = 0,
#                      baseline_metabolism = 10000, movement_bouts = 4,
#                      feed_while_moving = TRUE, pesticide_while_moving = TRUE);

## ---- eval = FALSE------------------------------------------------------------
#  sim <- run_farm_sim(mine_output = mg_v1, N = 1000, neutral_loci = 1000,
#                      xdim = 64, ydim = 64, repro = "asexual", max_age = 4,
#                      selfing = TRUE, food_consume = c(2, 2),
#                      pesticide_consume = c(0, 0), food_needed_surv = 1,
#                      food_needed_repr = 1, reproduction_type = "food_based",
#                      pesticide_tolerated_surv = 0, pesticide_rotation_type  = 2,
#                      crop_rotation_type = 2, min_age_reproduce = 4,
#                      max_age_feed = 3, farms = 24, time_steps = 48,
#                      mutation_pr = 0.001, crossover_pr = 0.1, net_mu_layers = 0,
#                      crop_rotation_time = 12, pesticide_rotation_time = 12,
#                      crop_per_cell = 4, pesticide_per_cell = 0.4,
#                      crop_number = 2, pesticide_number = 2, min_age_move = 4,
#                      age_food_threshold = 3, min_age_feed = 1,
#                      pesticide_start = 100, print_last = TRUE,
#                      mating_distance = 2, metabolism = 0);

## ---- eval = FALSE------------------------------------------------------------
#  sim <- run_farm_sim(mine_output = mg_v1, N = 1000, neutral_loci = 1000,
#                      xdim = 64, ydim = 64, repro = "asexual", max_age = 4,
#                      selfing = TRUE, food_consume = c("T1", "T2"),
#                      pesticide_consume = c("T3", "T4"),
#                      food_needed_surv = 1, food_needed_repr = 1,
#                      reproduction_type = "food_based",
#                      pesticide_tolerated_surv = 0, pesticide_rotation_type  = 2,
#                      crop_rotation_type = 2, min_age_reproduce = 4,
#                      max_age_feed = 3, farms = 24, time_steps = 48,
#                      mutation_pr = 0.001, crossover_pr = 0.1, net_mu_layers = 0,
#                      crop_rotation_time = 12, pesticide_rotation_time = 12,
#                      crop_per_cell = 4, pesticide_per_cell = 0.4,
#                      crop_number = 2, pesticide_number = 2, min_age_move = 4,
#                      age_food_threshold = 3, min_age_feed = 1,
#                      pesticide_start = 100, print_last = TRUE,
#                      mating_distance = 2, metabolism = 0);

## ---- echo = FALSE------------------------------------------------------------
par(oldpar);

