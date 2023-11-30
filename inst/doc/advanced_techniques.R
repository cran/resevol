## ---- echo = FALSE------------------------------------------------------------
oldpar <- par();

## -----------------------------------------------------------------------------
library("resevol");
gmt <- matrix(data = c(1.0, -0.5, 0.2, 0.2, -0.5, 1.0, 0.2, 0.2, 0.2, 
                       0.2, 0.4, -0.6, 0.2, 0.2, -0.6, 0.4), nrow = 4);
print(gmt);

## ---- eval = FALSE------------------------------------------------------------
#  set.seed(2022);
#  mg  <- mine_gmatrix(gmatrix = gmt, loci = 12, indivs = 2000, npsize = 12000,
#                      max_gen = 5400, sampleK = 1200, chooseK = 6, layers = 4,
#                      mu_pr = 0.2, pr_cross = 0.2, mu_sd = 0.004,
#                      term_cri = -8);

## -----------------------------------------------------------------------------
load(system.file("advanced_eg.rda", package = "resevol"));

## ---- echo = FALSE------------------------------------------------------------
print(mg[[6]]);

## ---- fig.height = 6, fig.width = 6, fig.align = "left", fig.cap = "**Figure 1**. *Distribution of stress for initialised pest trait covariances in the resevol R package. Stress values are computed using genome values produced by the `mine_gmatrix` function for 1000 replicate populations of initialised pest loci values.*"----
sim_stress <- stress_test(mine_output = mg, indivs = 2000, reps = 1000);
hist(x = sim_stress, main = "", breaks = 20, xlab = "Initialised stress");

## -----------------------------------------------------------------------------
simple_terrain              <- matrix(data = 0, nrow = 12, ncol = 12)
simple_terrain[1:2, 1:3]    <- 1;
simple_terrain[3:7, 1:2]    <- 1;
simple_terrain[8:12, 1:2]   <- 2;
simple_terrain[1:2, 6:9]    <- 3;
simple_terrain[3:6, 7:9]    <- 3;
simple_terrain[1:6, 10:12]  <- 4;
simple_terrain[10:12, 7]    <- 5;
simple_terrain[7:12, 8:12]  <- 5;
simple_terrain[3:12, 3:4]   <- 6;
simple_terrain[7:9, 5]      <- 6;
simple_terrain[1:2, 4:5]    <- 7;
simple_terrain[3:6, 5:6]    <- 7;
simple_terrain[7:9, 6:7]    <- 7;
simple_terrain[10:12, 5:6]  <- 7;
print(simple_terrain)

## ---- fig.width=5, fig.height=5, fig.align = "left", fig.cap = "**Figure 2**. *Image representing a simple landscape for the resevol package created from a 12 by 12 matrix of natural numbers from 1-7.*"----
par(mar = c(0, 0, 0, 0));
image(t(simple_terrain), xaxt = "n", yaxt = "n", useRaster = TRUE);

## ---- fig.width=7, fig.height=7, fig.align = "left", fig.cap = "**Figure 3**. *A complex landscape to be used in the resevol R package, including 14 separate farms, grassland, forest, and water. The image is represented in the code by a matrix in which element numbers correspond to different terrain colours (e.g., elements corresponding to water are numbered 17).*"----
land_file <- system.file("landscape_eg.csv", package = "resevol");
land_dat  <- read.csv(file = land_file, header = FALSE);
land_eg   <- t(as.matrix(land_dat));
farm_cols <- c("#f4eadc", "#6a4b20", "#cea05f", "#e1c59d", "#a97833", "#cea05f",
               "#f2e6d6", "#6a4b20", "#cc9c59", "#dfc197", "#a27331", "#f0e3d0",
               "#5d421c", "#ca9852");
land_cols <- c(farm_cols, "#00ab41", "#234F1E", "#2832C2");
par(mar = c(0, 0, 0, 0));
image(land_eg, xaxt = "n", yaxt = "n", col = land_cols);
points(x = 0.2, y = 0.05, cex = 9, pch = 20);
text(x = 0.2, y = 0.05, labels = "1", cex = 2, col = "red");
points(x = 0.4, y = 0.1, cex = 9, pch = 20);
text(x = 0.4, y = 0.1, labels = "2", cex = 2, col = "red");
points(x = 0.4, y = 0.1, cex = 9, pch = 20);
text(x = 0.4, y = 0.1, labels = "2", cex = 2, col = "red");
points(x = 0.25, y = 0.27, cex = 9, pch = 20);
text(x = 0.25, y = 0.27, labels = "3", cex = 2, col = "red");
points(x = 0.20, y = 0.43, cex = 9, pch = 20);
text(x = 0.20, y = 0.43, labels = "4", cex = 2, col = "red");
points(x = 0.42, y = 0.48, cex = 9, pch = 20);
text(x = 0.42, y = 0.48, labels = "5", cex = 2, col = "red");
points(x = 0.28, y = 0.58, cex = 9, pch = 20);
text(x = 0.28, y = 0.58, labels = "6", cex = 2, col = "red");
points(x = 0.1, y = 0.8, cex = 9, pch = 20);
text(x = 0.1, y = 0.8, labels = "7", cex = 2, col = "red");
points(x = 0.7, y = 0.05, cex = 9, pch = 20);
text(x = 0.7, y = 0.05, labels = "8", cex = 2, col = "red");
points(x = 0.9, y = 0.2, cex = 9, pch = 20);
text(x = 0.9, y = 0.2, labels = "9", cex = 2, col = "red");
points(x = 0.85, y = 0.4, cex = 9, pch = 20);
text(x = 0.85, y = 0.4, labels = "10", cex = 2, col = "red");
points(x = 0.92, y = 0.58, cex = 9, pch = 20);
text(x = 0.92, y = 0.58, labels = "11", cex = 2, col = "red");
points(x = 0.88, y = 0.76, cex = 9, pch = 20);
text(x = 0.88, y = 0.76, labels = "12", cex = 2, col = "red");
points(x = 0.86, y = 0.93, cex = 9, pch = 20);
text(x = 0.86, y = 0.93, labels = "13", cex = 2, col = "red");
points(x = 0.52, y = 0.91, cex = 9, pch = 20);
text(x = 0.52, y = 0.91, labels = "14", cex = 2, col = "red");
points(x = 0.05, y = 0.59, cex = 7, pch = 20);
text(x = 0.05, y = 0.59, labels = "15", cex = 1.7, col = "white");
points(x = 0.46, y = 0.28, cex = 7, pch = 20);
text(x = 0.46, y = 0.28, labels = "16", cex = 1.7, col = "white");
points(x = 0.62, y = 0.36, cex = 7, pch = 20);
text(x = 0.62, y = 0.36, labels = "17", cex = 1.7, col = "white");

## -----------------------------------------------------------------------------
initial_crop <- c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3);

## -----------------------------------------------------------------------------
rotate_crop       <- matrix(data = 0, nrow = 3, ncol = 3);
rotate_crop[1, 2] <- 1;
rotate_crop[2, 1] <- 1;
rotate_crop[3, 3] <- 1;
print(rotate_crop);

## -----------------------------------------------------------------------------
initial_pesticide <- c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3);

## -----------------------------------------------------------------------------
rotate_pesticide           <- matrix(data = 0, nrow = 3, ncol = 3);
rotate_pesticide[1, 1]     <- 1;
rotate_pesticide[2, 2]     <- 1;
rotate_pesticide[3, 3]     <- 1;
print(rotate_pesticide);

## ---- eval = FALSE------------------------------------------------------------
#  set.seed(2022);
#  sim1       <- run_farm_sim(mine_output             = mg,
#                            terrain                  = land_eg,
#                            crop_init                = initial_crop,
#                            crop_rotation_type       = rotate_crop,
#                            pesticide_init           = initial_pesticide,
#                            pesticide_rotation_type  = rotate_pesticide,
#                            food_consume             = c("T1", "T2", 0),
#                            pesticide_consume        = c("T3", "T4", 0),
#                            crop_number              = 3,
#                            pesticide_number         = 3,
#                            trait_means              = c(2, 2, 0.0, 0.0),
#                            max_age                  = 6,
#                            min_age_feed             = 0,
#                            max_age_feed             = 2,
#                            min_age_move             = 3,
#                            max_age_move             = 6,
#                            min_age_metabolism       = 3,
#                            max_age_metabolism       = 6,
#                            metabolism               = 0.5,
#                            food_needed_surv         = 1,
#                            reproduction_type        = "food_based",
#                            food_needed_repr         = 2,
#                            N                        = 1000,
#                            repro                    = "biparental",
#                            mating_distance          = 4,
#                            rand_age                 = TRUE,
#                            pesticide_tolerated_surv = 2,
#                            movement_bouts           = 4,
#                            move_distance            = 2,
#                            crop_per_cell            = 10,
#                            crop_sd                  = 0,
#                            pesticide_per_cell       = 1,
#                            pesticide_sd             = 0,
#                            crop_rotation_time       = 18,
#                            pesticide_rotation_time  = 9,
#                            time_steps               = 240,
#                            pesticide_start          = 81,
#                            immigration_rate         = 100,
#                            land_edge                = "reflect",
#                            mutation_pr              = 0.01,
#                            crossover_pr             = 0.01,
#                            print_gens              = TRUE,
#                            print_last              = TRUE);

## -----------------------------------------------------------------------------
population_data_file_sim1 <- system.file("population_data_sim1.csv", 
                                         package = "resevol");
population_data_sim1      <- read.csv(file = population_data_file_sim1);

## -----------------------------------------------------------------------------
print(head(population_data_sim1));

## ---- fig.height = 7, fig.width = 7, fig.align = "left", fig.cap = "**Figure 4**. *Pest abundance over time in a population simulated using the resevol R package in which pesticides are not rotated over time. The grey shaded regions show individual crop seasons; pesticides are applied at the start and midpoint of each seasons, beginning at the time step indicated by the red vertical line.*"----
mbox <- function(x0, x1, y0, y1){
    xx <- seq(from=x0, to=x1, length.out = 100);
    yy <- seq(from=y0, to=y1, length.out = 100);
    xd <- c(rep(x0, 100), xx, rep(x1,100), rev(xx));
    yd <- c(yy, rep(y1,100), rev(yy), rep(y0, 100));
    return(list(x=xd, y=yd));
}
season          <- seq(from = 0, to = 240, by = 18);
blocks          <- length(season) - 1;
plot(x = population_data_sim1[["time_step"]], type = "n",
     y = population_data_sim1[["population_size"]], cex.lab = 1.25, 
     cex.axis = 1.25, ylab = "Pest population abundance", xlab = "Time step", 
     ylim = c(0, max(population_data_sim1[["population_size"]])));
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 21000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim1[["time_step"]], 
     y = population_data_sim1[["population_size"]], type = "l", lwd = 2);
abline(v = 81, lwd = 3, col = "red");
box();

## ---- fig.height = 7, fig.width = 7, fig.align = "left", fig.cap = "**Figure 5**. *Pest consumption over time in a population simulated using the resevol R package in which pesticides are not rotated over time. Panels show (a) mean amount of crop 1 and (b) crop 2 consumed per pest, and (c) mean amount of pesticide 1 and (d) pesticide 2 consumed per pest. Grey shaded regions show individual crop seasons; pesticides are applied at the start and midpoint of each seasons, beginning at the time step 100.*"----
par(mfrow = c(2, 2), mar = c(1, 4, 1, 1));
plot(x = population_data_sim1[["time_step"]], type = "n", ylim = c(0, 0.8),
     y = population_data_sim1[["mean_food1_consumed"]], cex.lab = 1.25, 
     ylab = "Mean crop 1 consumed", xlab = "Time step", xaxt = "n");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim1[["time_step"]], 
       y = population_data_sim1[["mean_food1_consumed"]], type = "l", lwd = 2);
text(x = 5, y = 0.79, labels = "a", cex = 2);
box();
par(mar = c(1, 4, 1, 1));
plot(x = population_data_sim1[["time_step"]], type = "n", ylim = c(0, 0.8),
     y = population_data_sim1[["mean_food2_consumed"]], cex.lab = 1.25 ,
     ylab = "Mean crop 2 consumed", xlab = "Time step", xaxt = "n");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim1[["time_step"]], 
       y = population_data_sim1[["mean_food2_consumed"]], type = "l", lwd = 2);
text(x = 5, y = 0.79, labels = "b", cex = 2);
box();
par(mar = c(4, 4, 1, 1));
plot(x = population_data_sim1[["time_step"]], type = "n", ylim = c(0, 0.8),
     y = population_data_sim1[["mean_pesticide1_consumed"]], cex.lab = 1.25, 
     ylab = "Mean pesticide 1 consumed", xlab = "Time step");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim1[["time_step"]], 
       y = population_data_sim1[["mean_pesticide1_consumed"]], type = "l", lwd = 2);
box();
text(x = 5, y = 0.78, labels = "c", cex = 2);
par(mar = c(4, 4, 1, 1));
plot(x = population_data_sim1[["time_step"]], type = "n", ylim = c(0, 0.8),
     y = population_data_sim1[["mean_pesticide2_consumed"]], cex.lab = 1.25, 
     ylab = "Mean pesticide 2 consumed", xlab = "Time step");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim1[["time_step"]], 
       y = population_data_sim1[["mean_pesticide2_consumed"]], type = "l", 
       lwd = 2);
box();
text(x = 5, y = 0.78, labels = "d", cex = 2);

## ---- fig.height = 7, fig.width = 7, fig.align = "left", fig.cap = "**Figure 6**. *Mean values of pest traits T1 (crop 1 consumption rate), T2 (crop 2 consumption rate), T3 (pesticide 1 consumption rate), and T4 (pesticide 2 consumption rate) in a population simulated using the resevol R package in which pesticides are not rotated over time. Grey shaded regions show individual crop seasons; pesticides are applied at the start and midpoint of each seasons, beginning at the time step 100.*"----
par(mfrow = c(2, 2), mar = c(1, 4, 1, 1));
plot(x = population_data_sim1[["time_step"]], type = "n", ylim = c(2, 3.6),
     y = population_data_sim1[["trait1_mean_value"]], cex.lab = 1.25, 
     ylab = "Mean crop 1 consumption trait", xlab = "Time step", xaxt = "n");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim1[["time_step"]], 
     y = population_data_sim1[["trait1_mean_value"]], type = "l", lwd = 2);
text(x = 5, y = 3.56, labels = "a", cex = 2);
box();
par(mar = c(1, 4, 1, 1));
plot(x = population_data_sim1[["time_step"]], type = "n", ylim = c(2, 3.6),
     y = population_data_sim1[["trait2_mean_value"]], cex.lab = 1.25 ,
     ylab = "Mean crop 2 consumption trait", xlab = "Time step", xaxt = "n");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim1[["time_step"]], 
       y = population_data_sim1[["trait2_mean_value"]], type = "l", lwd = 2);
text(x = 5, y = 3.56, labels = "b", cex = 2);
box();
par(mar = c(4, 4, 1, 1));
plot(x = population_data_sim1[["time_step"]], type = "n", ylim = c(-0.8, 1.4),
     y = population_data_sim1[["trait3_mean_value"]], cex.lab = 1.25, 
     ylab = "Mean pesticide 1 consumption trait", xlab = "Time step");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim1[["time_step"]], 
       y = population_data_sim1[["trait3_mean_value"]], type = "l", lwd = 2);
box();
text(x = 5, y = 1.34, labels = "c", cex = 2);
par(mar = c(4, 4, 1, 1));
plot(x = population_data_sim1[["time_step"]], type = "n", ylim = c(-0.8, 1.4),
     y = population_data_sim1[["trait4_mean_value"]], cex.lab = 1.25, 
     ylab = "Mean pesticide 2 consumption trait", xlab = "Time step");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim1[["time_step"]], 
       y = population_data_sim1[["trait4_mean_value"]], type = "l", lwd = 2);
box();
text(x = 5, y = 1.34, labels = "d", cex = 2);

## -----------------------------------------------------------------------------
last_time_step_file_sim1 <- system.file("last_time_step_sim1.csv", 
                                        package = "resevol");
pop_last_time_step_sim1  <- read.csv(file = last_time_step_file_sim1);

## ---- fig.width=7, fig.height=7, fig.align = "left", fig.cap = "**Figure 7**. *A complex landscape for simulating the ecology of pests and the evolution of pesticide resistance on farmland in the resevol R package. Terrain includes farms (brown colours), grassland (light green), forest (dark green), and water (blue). Black points show the locations of individual pests after 240 time steps for a simulation in which pesticides are not rotated over time.*"----
par(mar = c(0, 0, 0, 0));
landscape    <- land_eg;
for(i in 1:dim(pop_last_time_step_sim1)[1]){
    xloc <- pop_last_time_step_sim1[i, 1] + 1;
    yloc <- pop_last_time_step_sim1[i, 2] + 1;
    landscape[xloc, yloc] <- 18;
}
land_cols <- c(land_cols, "#000000");
image(landscape, xaxt = "n", yaxt = "n", col = land_cols);
head(pop_last_time_step_sim1);
print(land_cols);

## ---- eval = TRUE-------------------------------------------------------------
end_trait_covs <- cov(pop_last_time_step_sim1[,3:6]);
print(end_trait_covs);

## -----------------------------------------------------------------------------
rotate_pesticide           <- matrix(data = 0, nrow = 3, ncol = 3);
rotate_pesticide[1, 2]     <- 1;
rotate_pesticide[2, 1]     <- 1;
rotate_pesticide[3, 3]     <- 1;

## ---- eval = FALSE------------------------------------------------------------
#  set.seed(2022);
#  sim2       <- run_farm_sim(mine_output             = mg,
#                            terrain                  = land_eg,
#                            crop_init                = initial_crop,
#                            crop_rotation_type       = rotate_crop,
#                            pesticide_init           = initial_pesticide,
#                            pesticide_rotation_type  = rotate_pesticide,
#                            food_consume             = c("T1", "T2", 0),
#                            pesticide_consume        = c("T3", "T4", 0),
#                            crop_number              = 3,
#                            pesticide_number         = 3,
#                            trait_means              = c(2, 2, 0.0, 0.0),
#                            max_age                  = 6,
#                            min_age_feed             = 0,
#                            max_age_feed             = 2,
#                            min_age_move             = 3,
#                            max_age_move             = 6,
#                            min_age_metabolism       = 3,
#                            max_age_metabolism       = 6,
#                            metabolism               = 0.5,
#                            food_needed_surv         = 1,
#                            reproduction_type        = "food_based",
#                            food_needed_repr         = 2,
#                            N                        = 1000,
#                            repro                    = "biparental",
#                            mating_distance          = 4,
#                            rand_age                 = TRUE,
#                            pesticide_tolerated_surv = 2,
#                            movement_bouts           = 4,
#                            move_distance            = 2,
#                            crop_per_cell            = 10,
#                            crop_sd                  = 0,
#                            pesticide_per_cell       = 1,
#                            pesticide_sd             = 0,
#                            crop_rotation_time       = 18,
#                            pesticide_rotation_time  = 9,
#                            time_steps               = 240,
#                            pesticide_start          = 81,
#                            immigration_rate         = 100,
#                            land_edge                = "reflect",
#                            mutation_pr              = 0.01,
#                            crossover_pr             = 0.01,
#                            print_gens              = TRUE,
#                            print_last              = TRUE);

## -----------------------------------------------------------------------------
population_data_file_sim2 <- system.file("population_data_sim2.csv", 
                                    package = "resevol");
population_data_sim2      <- read.csv(file = population_data_file_sim2);

## ---- echo = FALSE, fig.height = 7, fig.width = 7, fig.align = "left", fig.cap = "**Figure 8**. *Pest abundance over time in a population simulated using the resevol R package in which pesticides rotates every 9 time steps. The grey shaded regions show individual crop seasons; pesticides are applied at the start and midpoint of each seasons, beginning at the time step indicated by the red vertical line.*"----
mbox <- function(x0, x1, y0, y1){
    xx <- seq(from=x0, to=x1, length.out = 100);
    yy <- seq(from=y0, to=y1, length.out = 100);
    xd <- c(rep(x0, 100), xx, rep(x1,100), rev(xx));
    yd <- c(yy, rep(y1,100), rev(yy), rep(y0, 100));
    return(list(x=xd, y=yd));
}
season          <- seq(from = 0, to = 240, by = 18);
blocks          <- length(season) - 1;
plot(x = population_data_sim2[["time_step"]], type = "n",
     y = population_data_sim2[["population_size"]], cex.lab = 1.25, 
     cex.axis = 1.25, ylab = "Pest population abundance", xlab = "Time step", 
     ylim = c(0, max(population_data_sim2[["population_size"]])));
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 21000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim2[["time_step"]], 
     y = population_data_sim2[["population_size"]], type = "l", lwd = 2);
abline(v = 81, lwd = 3, col = "red");
box();

## ---- echo = FALSE, fig.height = 7, fig.width = 7, fig.align = "left", fig.cap = "**Figure 9**. *Pest consumption over time in a population simulated using the resevol R package in which pesticides rotate every 9 time steps. Panels show (a) mean amount of crop 1 and (b) crop 2 consumed per pest, and (c) mean amount of pesticide 1 and (d) pesticide 2 consumed per pest. Grey shaded regions show individual crop seasons; pesticides are applied at the start and midpoint of each seasons, beginning at the time step 100.*"----
par(mfrow = c(2, 2), mar = c(1, 4, 1, 1));
plot(x = population_data_sim2[["time_step"]], type = "n", ylim = c(0, 0.8),
     y = population_data_sim2[["mean_food1_consumed"]], cex.lab = 1.25, 
     ylab = "Mean crop 1 consumed", xlab = "Time step", xaxt = "n");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim2[["time_step"]], 
       y = population_data_sim2[["mean_food1_consumed"]], type = "l", lwd = 2);
text(x = 5, y = 0.79, labels = "a", cex = 2);
box();
par(mar = c(1, 4, 1, 1));
plot(x = population_data_sim2[["time_step"]], type = "n", ylim = c(0, 0.8),
     y = population_data_sim2[["mean_food2_consumed"]], cex.lab = 1.25 ,
     ylab = "Mean crop 2 consumed", xlab = "Time step", xaxt = "n");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim2[["time_step"]], 
       y = population_data_sim2[["mean_food2_consumed"]], type = "l", lwd = 2);
text(x = 5, y = 0.79, labels = "b", cex = 2);
box();
par(mar = c(4, 4, 1, 1));
plot(x = population_data_sim2[["time_step"]], type = "n", ylim = c(0, 0.8),
     y = population_data_sim2[["mean_pesticide1_consumed"]], cex.lab = 1.25, 
     ylab = "Mean pesticide 1 consumed", xlab = "Time step");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim2[["time_step"]], 
       y = population_data_sim2[["mean_pesticide1_consumed"]], type = "l", lwd = 2);
box();
text(x = 5, y = 0.78, labels = "c", cex = 2);
par(mar = c(4, 4, 1, 1));
plot(x = population_data_sim2[["time_step"]], type = "n", ylim = c(0, 0.8),
     y = population_data_sim2[["mean_pesticide2_consumed"]], cex.lab = 1.25, 
     ylab = "Mean pesticide 2 consumed", xlab = "Time step");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim2[["time_step"]], 
       y = population_data_sim2[["mean_pesticide2_consumed"]], type = "l", 
       lwd = 2);
box();
text(x = 5, y = 0.78, labels = "d", cex = 2);

## ---- echo = FALSE, fig.height = 7, fig.width = 7, fig.align = "left", fig.cap = "**Figure 10**. *Mean values of pest traits T1 (crop 1 consumption rate), T2 (crop 2 consumption rate), T3 (pesticide 1 consumption rate), and T4 (pesticide 2 consumption rate) in a population simulated using the resevol R package in which pesticides rotate every 9 time steps. Grey shaded regions show individual crop seasons; pesticides are applied at the start and midpoint of each seasons, beginning at the time step 100.*"----
par(mfrow = c(2, 2), mar = c(1, 4, 1, 1));
plot(x = population_data_sim2[["time_step"]], type = "n", ylim = c(2, 3.6),
     y = population_data_sim2[["trait1_mean_value"]], cex.lab = 1.25, 
     ylab = "Mean crop 1 consumption trait", xlab = "Time step", xaxt = "n");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim2[["time_step"]], 
     y = population_data_sim2[["trait1_mean_value"]], type = "l", lwd = 2);
text(x = 5, y = 3.56, labels = "a", cex = 2);
box();
par(mar = c(1, 4, 1, 1));
plot(x = population_data_sim2[["time_step"]], type = "n", ylim = c(2, 3.6),
     y = population_data_sim2[["trait2_mean_value"]], cex.lab = 1.25 ,
     ylab = "Mean crop 2 consumption trait", xlab = "Time step", xaxt = "n");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim2[["time_step"]], 
       y = population_data_sim2[["trait2_mean_value"]], type = "l", lwd = 2);
text(x = 5, y = 3.56, labels = "b", cex = 2);
box();
par(mar = c(4, 4, 1, 1));
plot(x = population_data_sim2[["time_step"]], type = "n", ylim = c(-0.8, 1.4),
     y = population_data_sim2[["trait3_mean_value"]], cex.lab = 1.25, 
     ylab = "Mean pesticide 1 consumption trait", xlab = "Time step");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim2[["time_step"]], 
       y = population_data_sim2[["trait3_mean_value"]], type = "l", lwd = 2);
box();
text(x = 5, y = 1.34, labels = "c", cex = 2);
par(mar = c(4, 4, 1, 1));
plot(x = population_data_sim2[["time_step"]], type = "n", ylim = c(-0.8, 1.4),
     y = population_data_sim2[["trait4_mean_value"]], cex.lab = 1.25, 
     ylab = "Mean pesticide 2 consumption trait", xlab = "Time step");
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = population_data_sim2[["time_step"]], 
       y = population_data_sim2[["trait4_mean_value"]], type = "l", lwd = 2);
box();
text(x = 5, y = 1.34, labels = "d", cex = 2);

## ---- echo = FALSE, fig.width=7, fig.height=7, fig.align = "left", fig.cap = "**Figure 11**. *A complex landscape for simulating the ecology of pests and the evolution of pesticide resistance on farmland in the resevol R package. Terrain includes farms (brown colours), grassland (light green), forest (dark green), and water (blue). Black points show the locations of individual pests after 240 time steps for a simulation in which pesticides are not rotated over time.*"----
last_time_step_file_sim2 <- system.file("last_time_step_sim2.csv", 
                                        package = "resevol");
pop_last_time_step_sim2  <- read.csv(file = last_time_step_file_sim2);
par(mar = c(0, 0, 0, 0));
landscape    <- land_eg;
for(i in 1:dim(pop_last_time_step_sim2)[1]){
    xloc <- pop_last_time_step_sim2[i, 1] + 1;
    yloc <- pop_last_time_step_sim2[i, 2] + 1;
    landscape[xloc, yloc] <- 18;
}
image(landscape, xaxt = "n", yaxt = "n", col = land_cols);

## ---- eval = TRUE, echo = FALSE-----------------------------------------------
end_trait_covs2 <- cov(pop_last_time_step_sim2[,3:6]);
print(end_trait_covs2);

## ---- echo = FALSE------------------------------------------------------------
suppressWarnings(par(oldpar));

