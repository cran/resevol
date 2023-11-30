## ---- echo = FALSE------------------------------------------------------------
library("resevol");
oldpar <- par();

## ----network, echo = FALSE, fig.width = 6, fig.height = 4, fig.cap = "**Figure 1**. *Example network mapping loci (green circles) to traits (red diamonds) through an intermediate set of hidden layers (blue squares) in the `mine_gmatrix` function. Individual genomes in the resevol R package consist of standard random normal values for loci, real values for black arrows linking nodes, and real values for traits. Values shown for loci and arrows are an example for illustration.*"----
par(mar = c(0.2, 0.2, 0, 0.2));
plot(x = 0, y = 0, type = "n", xlim = c(0, 1000), ylim = c(0, 1000), 
     xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "");
mtext(side = 1, text = "Internal nodes", cex = 1.25, col = "blue", line = -3);
mtext(side = 2, text = "Loci", cex = 1.25, col = "darkgreen", line = -2);
mtext(side = 4, text = "Traits", cex = 1.25, col = "red", line = -2);
Lx  <- rep(x = 100, times = 3);
Ly  <- seq(from = 250, to = 750, length = 3);
H1x <- rep(x = 400, times = 2);
H1y <- seq(from = 250, to = 750, length = 2);
H2x <- rep(x = 600, times = 2);
H2y <- seq(from = 250, to = 750, length = 2);
Tx  <- rep(x = 900, times = 2);
Ty  <- seq(from = 250, to = 750, length = 2);
points(x = Tx, y = Ty, pch = 23, cex = 8, bg = "red", lwd = 2);
sep_L_points <- c(-50, 0, 50);
sep_N_points <- c(-30, 30);
arrows(x0 = H2x + 52, x1 = Tx[1] - 58, y0 = H2y, y1 = Ty[1] + sep_N_points, 
       length = 0.08, lwd = 3);
arrows(x0 = H2x + 52, x1 = Tx[2] - 58, y0 = H2y, y1 = Ty[2] + sep_N_points, 
       length = 0.08, lwd = 3);
arrows(x0 = Lx + 38, x1 = H1x[1] - 58, y0 = Ly, y1 = H1y[1] + sep_L_points, 
       length = 0.08, lwd = 3);
points(x = H2x, y = H2y, pch = 22, cex = 9, bg = "blue", lwd = 2);
arrows(x0 = H1x + 52, x1 = H2x[1] - 58, y0 = H1y, y1 = H2y[1] + sep_N_points, 
       length = 0.08, lwd = 3);
arrows(x0 = H1x + 52, x1 = H2x[2] - 58, y0 = H1y, y1 = H2y[2] + sep_N_points, 
       length = 0.08, lwd = 3);
points(x = H1x, y = H1y, pch = 22, cex = 9, bg = "blue", lwd = 2);
arrows(x0 = Lx + 38, x1 = H1x[2] - 58, y0 = Ly, y1 = H1y[2] + sep_L_points, 
       length = 0.08, lwd = 3);
points(x = Lx, y = Ly, pch = 21, cex = 6, bg = "darkgreen", lwd = 2);
srnm <- c(0.2, -0.1, 0.1);
text(x = Lx, y = Ly, labels = srnm, col = "white", cex = 0.9);
text(x = 200, y = 790, labels = "0.4", srt = 10, cex = 0.8);
text(x = 200, y = 660, labels = "0.1", srt = -58, cex = 0.8);
text(x = 270, y = 695, labels = "-1.2", srt = 46, cex = 0.8);
text(x = 180, y = 480, labels = "-0.6", srt = -40, cex = 0.8);
text(x = 160, y = 340, labels = "0.2", srt = 54, cex = 0.8);
text(x = 220, y = 255, labels = "2.2", srt = -12, cex = 0.8);
text(x = 490, y = 790, labels = "1.3", srt = 10, cex = 0.8);
text(x = 538, y = 410, labels = "-0.5", srt = -76, cex = 0.8);
text(x = 460, y = 390, labels = "-0.1", srt = 72, cex = 0.8);
text(x = 495, y = 260, labels = "0.0", srt = -12, cex = 0.8);
text(x = 720, y = 790, labels = "2.2", srt = 8, cex = 0.8);
text(x = 720, y = 630, labels = "0.9", srt = -58, cex = 0.8);
text(x = 690, y = 400, labels = "1.6", srt = 54, cex = 0.8);
text(x = 740, y = 260, labels = "0.8", srt = -8, cex = 0.8);
Loci  <- matrix(data = c(0.1, -0.1, 0.2), nrow = 1);
Hdn1  <- matrix(data = c(0.4, -1.2, 0.2, 0.1, -0.6, 2.2), nrow = 3, ncol = 2);
Hdn2  <- matrix(data = c(1.3, -0.1, -0.5, 0.0), nrow = 2, ncol = 2);
Hdn3  <- matrix(data = c(2.2, 1.6, 0.9, 0.8), nrow = 2, ncol = 2);
Tvals <- Loci %*% Hdn1 %*% Hdn2 %*% Hdn3;
text(x = Tx, y = Ty, labels = rev(Tvals), cex = 0.8);

## ---- eval = FALSE------------------------------------------------------------
#  trait_covs  <- matrix(data = c(1, -0.4, -0.4, 1), nrow = 2, ncol = 2);
#  new_network <- mine_gmatrix(loci = 3, layers = 2, gmatrix = trait_covs,
#                              max_gen = 1000, term_cri = -6.0, prnt_out = FALSE);

## ---- echo = FALSE------------------------------------------------------------
trait_covs  <- matrix(data = c(1, -0.4, -0.4, 1), nrow = 2, ncol = 2);
load(system.file("new_network.rda", package = "resevol"));

## ---- echo = FALSE------------------------------------------------------------
print(new_network[[6]]);

## ----flowchart, echo=FALSE, fig.height=3, fig.width=8, fig.cap = "**Figure 2**. *Overview of simulated events in the resevol R package. Note that metabolism, feeding, pesticide consumption, movement, and reproduction are all subject to a minimum and maximum pest age. Consequently, simulation order might not reflect the order of events from the perspective of a focal pest (e.g., pests might move from ages 1-2, but only feed from ages 2-4). Crops and pesticides are also not necessarily rotated in each time step (see Landscape). Statistics collected within a time step are printed to a CSV file.*"----
mbox <- function(x0, x1, y0, y1){
    xx <- seq(from=x0, to=x1, length.out = 100);
    yy <- seq(from=y0, to=y1, length.out = 100);
    xd <- c(rep(x0, 100), xx, rep(x1,100), rev(xx));
    yd <- c(yy, rep(y1,100), rev(yy), rep(y0, 100));
    return(list(x=xd, y=yd));
}
par(mar=c(0,0,0,0));
# ===============================================================
plot(x=0, y=0, type="n", xlim=c(0,100), ylim=c(40,95), xaxt="n", yaxt="n",
     xlab="",ylab="");
ibox <- mbox(x0 = 0,  x1 = 14, y0 = 63, y1 = 47);
polygon(x = ibox$x, y = ibox$y, lwd = 5, border = "black", col = "white");
cbox <- mbox(x0 = 0, x1 = 10, y0 = 88, y1 = 72);
polygon(x = cbox$x, y = cbox$y, lwd = 3, border = "black", col = "white");
ubox <- mbox(x0 = 15, x1 = 25, y0 = 88, y1 = 72);
polygon(x = ubox$x, y = ubox$y, lwd = 3, border = "black", col = "white");
nbox <- mbox(x0 = 30, x1 = 40, y0 = 88, y1 = 72);
polygon(x = nbox$x, y = nbox$y, lwd = 3, border = "black", col = "white");
fbox <- mbox(x0 = 45, x1 = 55, y0 = 88, y1 = 72);
polygon(x = fbox$x, y = fbox$y, lwd = 3, border = "black", col = "white");
tbox <- mbox(x0 = 60, x1 = 70, y0 = 88, y1 = 72);
polygon(x = tbox$x, y = tbox$y, lwd = 3, border = "black", col = "white");
tbox <- mbox(x0 = 75, x1 = 85, y0 = 88, y1 = 72);
polygon(x = tbox$x, y = tbox$y, lwd = 3, border = "black", col = "white");
tbox <- mbox(x0 = 90, x1 = 100, y0 = 88, y1 = 72);
polygon(x = tbox$x, y = tbox$y, lwd = 3, border = "black", col = "white");
ibox <- mbox(x0 = 90,  x1 = 100, y0 = 63, y1 = 47);
polygon(x = ibox$x, y = ibox$y, lwd = 3, border = "black", col = "white");
tbox <- mbox(x0 = 75, x1 = 85, y0 = 63, y1 = 47);
polygon(x = tbox$x, y = tbox$y, lwd = 3, border = "black", col = "white");
text(x=7,  y=59, labels="Initialise", col="black", cex = 0.8);
text(x=7,  y=55, labels="pests and", col="black", cex = 0.8);
text(x=7,  y=51, labels="landscape", col="black", cex = 0.8);
text(x=5, y=84, labels="Crops and", col="black", cex = 0.8);
text(x=5, y=80, labels="pesticides", col="black", cex = 0.8);
text(x=5, y=76, labels="rotated?", col="black", cex = 0.8);
text(x=20, y=82, labels="Pests", col="black", cex = 0.8);
text(x=20, y=78, labels="metabolise", col="black", cex = 0.8);
text(x=35, y=82, labels="Pests", col="black", cex = 0.8);
text(x=35, y=78, labels="feed", col="black", cex = 0.8);
text(x=50, y=82, labels="Pesticide", col="black", cex = 0.8);
text(x=50, y=78, labels="consumed", col="black", cex = 0.8);
text(x=65, y=82, labels="Pests", col="black", cex = 0.8);
text(x=65, y=78, labels="move", col="black", cex = 0.8);
text(x=80, y=82, labels="Pests", col="black", cex = 0.8);
text(x=80, y=78, labels="reproduce", col="black", cex = 0.8);
text(x=95, y=82, labels="Pest", col="black", cex = 0.8);
text(x=95, y=78, labels="mortality", col="black", cex = 0.8);
text(x=95,  y=57, labels="Pests", col="black", cex = 0.8);
text(x=95,  y=53, labels="immigrate", col="black", cex = 0.8);
text(x=80,  y=59, labels="Collect", col="black", cex = 0.8);
text(x=80,  y=55, labels="statistics", col="black", cex = 0.8);
text(x=80,  y=51, labels="as CSV", col="black", cex = 0.8);
arrows(x0=2, x1=2,   y0=63, y1=71.5, lwd=2, length=0.10);
arrows(x0=10, x1=14.5, y0=80, y1=80, lwd=2, length=0.10);
arrows(x0=25, x1=29.5, y0=80, y1=80, lwd=2, length=0.10);
arrows(x0=40, x1=44.5, y0=80, y1=80, lwd=2, length=0.10);
arrows(x0=55, x1=59.5, y0=80, y1=80, lwd=2, length=0.10);
arrows(x0=70, x1=74.5, y0=80, y1=80, lwd=2, length=0.10);
arrows(x0=85, x1=89.5, y0=80, y1=80, lwd=2, length=0.10);
arrows(x0=95, x1=95,   y0=72, y1=63.5, lwd=2, length=0.10);
arrows(x0=90, x1=85.5, y0=55, y1=55, lwd=2, length=0.10);
arrows(x0=75, x1=70.5, y0=60, y1=60, lwd=2, length=0.10);
arrows(x0=60, x1=55.5, y0=60, y1=60, lwd=2, length=0.10);
rbox <- mbox(x0 = 60, x1 = 70, y0 = 57, y1 = 63);
polygon(x = rbox$x, y = rbox$y, lwd = 3, border = "black", col = "black");
text(x=65, y=60, labels="Extinction?", col="white", cex = 0.75);
rbox <- mbox(x0 = 45, x1 = 55, y0 = 57, y1 = 63);
polygon(x = rbox$x, y = rbox$y, lwd = 3, border = "black", col = "black");
text(x=50, y=60, labels="Max. time?", col="white", cex = 0.75);
text(x = 58.5, y = 58.5, labels = "No", col = "black", cex = 0.75);
text(x = 43.5, y = 58.5, labels = "No", col = "black", cex = 0.75);
text(x = 68, y = 54, labels = "Yes", col = "black", cex = 0.75);
text(x = 53, y = 54, labels = "Yes", col = "black", cex = 0.75);
arrows(x0=66, x1=66, y0=57, y1=50, lwd=2, length=0.0);
arrows(x0=51, x1=51, y0=57, y1=50, lwd=2, length=0.0);
rbox <- mbox(x0 = 30, x1 = 40, y0 = 47, y1 = 53);
polygon(x = rbox$x, y = rbox$y, lwd = 3, border = "black", col = "black");
text(x=35, y=50, labels="STOP", col="white", cex = 1);
rbox <- mbox(x0 = 14, x1 = 20, y0 = 47, y1 = 63);
polygon(x = rbox$x, y = rbox$y, lwd = 5, border = "black", col = "black");
text(x=17, y=55, labels="START", col="white", cex = 1, srt = 90);
arrows(x0=66, x1=40.5, y0=50, y1=50, lwd=2, length=0.10);
arrows(x0=45, x1=35, y0=60, y1=60, lwd=2, length=0.0);
arrows(x0=35, x1=35, y0=60, y1=67, lwd=2, length=0.0);
arrows(x0=35, x1=7, y0=67, y1=67, lwd=2, length=0.0);
arrows(x0=7, x1=7, y0=67, y1=71.5, lwd=2, length=0.1);

## ---- eval = FALSE------------------------------------------------------------
#  sim <- run_farm_sim(mine_output = new_network, repro = "asexual",
#                      pesticide_number = 2, pesticide_init = "random",
#                      pesticide_consume = c("T1", "T2"), farms = 9,
#                      pesticide_rotation_time = 16, pesticide_rotation_type = 3,
#                      pesticide_tolerated_surv = 0, pesticide_per_cell = 1,
#                      crop_rotation_time = 16, crop_number = 1, crop_per_cell = 4,
#                      food_consume = 1, reproduction_type = "food_based",
#                      food_needed_surv = 1, food_needed_repr = 1, max_age = 4,
#                      min_age_feed = 0, max_age_feed = 2, min_age_move = 3,
#                      max_age_move = 4, min_age_reproduce = 4, print_gens = FALSE,
#                      max_age_reproduce = 4, age_pesticide_threshold = 2,
#                      rand_age = TRUE, move_distance = 2, immigration_rate = 10,
#                      time_steps = 160, print_last = TRUE, xdim = 64, ydim = 64,
#                      trait_means = c(0.1, 0.1), land_edge = "torus");

## ----overtime, echo=FALSE, fig.width=8, fig.cap = "**Figure 3**. *Agricultural pest ecological and evolutionary dynamics over 160 time steps from an individual-based simulation using the resevol R package. Panels show (a) pesticide abundance change, (b) mean realised amount of pesticides 1 and 2 uptaken per pest, (c) mean food consumed per pest, and (d) mean value of evolving traits 1 and 2 underlying pest uptake over time. Note that only pests with positive values for Traits 1 or 2 can uptake Pesticide 1 or 2, respectively. Pests with negative trait values will be unaffected by corresponding pesticides (i.e., pesticide consumption is a threshold trait), hence the difference between realised pesticide consumption (b) and the traits underlying it (d). White and grey vertical stripes indicate seasons of a single crop and pesticide application.*"----
mbox <- function(x0, x1, y0, y1){
    xx <- seq(from=x0, to=x1, length.out = 100);
    yy <- seq(from=y0, to=y1, length.out = 100);
    xd <- c(rep(x0, 100), xx, rep(x1,100), rev(xx));
    yd <- c(yy, rep(y1,100), rev(yy), rep(y0, 100));
    return(list(x=xd, y=yd));
}
pop_dat_file    <- system.file("population_data.csv", package = "resevol");
population_data <- read.csv(pop_dat_file);
row_number      <- 161; # dim(population_data)[1];
time_step       <- population_data[["time_step"]][2:row_number];
population_size <- population_data[["population_size"]][2:row_number];
food_consumed   <- population_data[["mean_food_consumed"]][2:row_number];
p1_consumed     <- population_data[["mean_pesticide1_consumed"]][2:row_number];
p2_consumed     <- population_data[["mean_pesticide2_consumed"]][2:row_number];
trait1          <- population_data[["trait1_mean_value"]][2:row_number];
trait2          <- population_data[["trait2_mean_value"]][2:row_number];
season          <- seq(from = 0, to = 160, by = 16);
blocks          <- length(season) - 1;
par(mfrow = c(2, 2), mar = c(1, 5, 3, 1));
# Pest abundance plot
plot(x = time_step, y = population_size, type = "l", lwd = 2, xaxt = "n",
     xlab = "", ylab = "Pest abundance", cex.lab = 1, ylim = c(0, 4600),
     cex.axis = 1);
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1000, y1 = 6000);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
text(x = 5, y = 4300, labels = "a", cex = 2);
points(x = time_step, y = population_size, type = "l", lwd = 2);
box();
# Pesticide consumed
plot(x = time_step, y = p1_consumed, type = "n", lwd = 2, xlab = "", xaxt = "n",
     ylab = "Mean pesticide consumed", cex.lab = 1, ylim = c(0, 0.4),
     cex.axis = 1);
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1, y1 = 3);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = time_step, y = p1_consumed, type = "l", lwd = 2);
points(x = time_step, y = p2_consumed, type = "l", lwd = 2, col = "red",
       lty = "dashed");
legend(x = 158, y = 0.4, col = c("black", "red"), lty = c("solid", "dashed"),
       legend = c("Pesticide 1", "Pesticide 2"), cex = 1, lwd = 2, bg = "white",
       xjust = 1);
text(x = 5, y = 0.37, labels = "b", cex = 2);
box();
par(mar = c(4, 5, 0, 1));
# Pest food consumption plot
plot(x = time_step, y = food_consumed, type = "n", lwd = 2, xlab = "Time step",
     ylab = "Mean food consumed", cex.lab = 1, ylim = c(0, 2),
     yaxt = "n");
axis(side = 2, at = c(0, 0.5, 1, 1.5), cex.axis = 1);
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1, y1 = 3);
    if(i %% 2 == 0){
      polygon(x = rbox$x, y = rbox$y, lwd = 2, border = NA, 
              col = "grey90");
    }
}
points(x = time_step, y = food_consumed, type = "l", lwd = 2);
text(x = 5, y = 0.2, labels = "c", cex = 2);
box();
# Pesticide traits
plot(x = time_step, y = trait1, type = "n", lwd = 2, xlab = "Time step",
     ylab = "Mean trait value", cex.lab = 1, ylim = c(-0.9, 0.1),
     yaxt = "n");
axis(side = 2, at = c(-0.8, -0.6, -0.4, -0.2, 0.0), cex.axis = 1);
for(i in 1:blocks){
    rbox <- mbox(x0 = season[i], x1 = season[i + 1], y0 = -1, y1 = 3);
    if(i %% 2 == 0){
        polygon(x = rbox$x, y = rbox$y, lwd = 3, border = NA, 
                col = "grey90");
    }
}
points(x = time_step, y = trait1, type = "l", lwd = 2);
points(x = time_step, y = trait2, type = "l", lwd = 2, col = "red",
       lty = "dashed");
legend(x = 158, y = 0.1, col = c("black", "red"), lty = c("solid", "dashed"),
       legend = c("Trait 1", "Trait 2"), cex = 1, lwd = 2, xjust = 1, 
       bg = "white");
text(x = 5, y = -0.82, labels = "d", cex = 2);
box();

## ----landscape, echo=FALSE, fig.width=6.4, fig.height=3, fig.cap = "**Figure 4**. *Locations of pests (black) across a landscape that includes nine farms (coloured blocks) in the last time step of a simulation using the resevol R package (left panel). The right panel shows which farms apply pesticide 1 (dark grey) and 2 (light grey).*"----
load(system.file("sim.rda", package = "resevol"));
last_time_file <- system.file("last_time_step.csv", package = "resevol");
last_time_step <- read.csv(last_time_file, header = FALSE);
landscape      <- sim[[2]][,,1];
landscape_vals <- unique(as.vector(landscape));
reord_land     <- sample(landscape_vals, size = length(landscape_vals));
for(i in 1:dim(landscape)[1]){
    for(j in 1:dim(landscape)[2]){
        old_val         <- landscape[i, j];
        new_val         <- reord_land[old_val];
        landscape[i, j] <- new_val;
    }
}
for(i in 1:dim(last_time_step)[1]){
    xloc <- last_time_step[i, 1] + 1;
    yloc <- last_time_step[i, 2] + 1;
    landscape[xloc, yloc] <- 10;
}
land_cols      <- c(hcl.colors(9, "YlOrRd", rev = TRUE), "#000000");
par(mfrow = c(1, 2), mar = c(0, 0.1, 0, 0.1));
image(landscape, xaxt = "n", yaxt = "n", col = land_cols);
box();
pesticide_position <- sim[[2]][,,12];
pesticide_cols     <- c("grey20", "grey60");
image(pesticide_position, xaxt = "n", yaxt = "n", col = pesticide_cols);
box();

## ---- echo = FALSE------------------------------------------------------------
suppressWarnings(par(oldpar));

