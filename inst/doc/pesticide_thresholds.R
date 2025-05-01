## ----echo = FALSE-------------------------------------------------------------
oldpar <- par();

## ----eval = FALSE-------------------------------------------------------------
# data(gmatrices);

## ----echo = FALSE-------------------------------------------------------------
mg_dat <- c(1.069133287, -0.072243449,  0.064732772, -0.056487174, 
            -0.072243449,  0.874066557,  0.007063041, -0.007362562, 0.064732772,  
            0.007063041,  0.948800593, -0.009399430, -0.056487174, -0.007362562, 
            -0.009399430,  1.001389459);
mg_n1  <- matrix(data = mg_dat, nrow = 4, ncol = 4, byrow = TRUE);
print(mg_n1);

## ----eval = FALSE-------------------------------------------------------------
# sim <- run_farm_sim(mine_output = mg_n1, repro = "asexual",
#                     pesticide_number = 1, pesticide_init = "random",
#                     pesticide_consume = c("T1"), farms = 5,
#                     pesticide_rotation_time = 16, pesticide_rotation_type = 3,
#                     pesticide_tolerated_surv = 0, pesticide_per_cell = 1,
#                     crop_rotation_time = 4, crop_number = 1, crop_per_cell = 8,
#                     food_consume = 1, reproduction_type = "food_based",
#                     food_needed_surv = 1, food_needed_repr = 1, max_age = 4,
#                     min_age_feed = 0, max_age_feed = 2, min_age_move = 3,
#                     max_age_move = 4, min_age_reproduce = 4, print_gens = FALSE,
#                     max_age_reproduce = 4, age_pesticide_threshold = 2,
#                     rand_age = TRUE, move_distance = 2, immigration_rate = 10,
#                     time_steps = 50, print_last = FALSE, xdim = 18, ydim = 18,
#                     trait_means = c(1, 1, 1, 1), land_edge = "torus",
#                     pesticide_threshold = c(10000, 0, 0, 4, 10000),
#                     pesticide_delay = c(1, 1, 10, 10, 10));

## ----echo = FALSE-------------------------------------------------------------
suppressWarnings(par(oldpar));

