#' Mine G-matrices
#'
#' Mine networks for establishing the link between genome and g-matrix. The 
#' output from this function is required to run individual-based simulations
#' in the rest of the package. The key input to this function, 'gmatrix', is a
#' (square) covariance matrix, with each row and column representing a trait for
#' the individual-based model. This function will run an evolutionary algorithm
#' to try to find a network that produces traits with the covariance structure
#' of gmatrix from a set of random standard normal values. The network from loci
#' values to trait values goes through a number of linked nodes to achieve this, 
#' and each generation tests the stress of the resulting network in terms of 
#' expected squared deviation of trait covariances from the input gmatrix. 
#' Simulations can take minutes to hours or longer, depending on parameters 
#' chosen and the number of traits. See vignettes for a more comprehensive
#' explanation for what this function is doing.
#'
#'@param loci The number of loci that individuals in the model will have
#'@param layers The number of layers in the network from loci to traits
#'@param indivs The number of individuals to test the covariance matrix
#'@param npsize The size of the network population in the evolutionary algorithm
#'@param mu_pr The probability of a network value to mutate
#'@param mu_sd The standard deviation of mutation effect size
#'@param max_gen The maximum number of generations of the evolutionary algorithm
#'@param pr_cross The probability of a crossover occurring for a network
#'@param sampleK Number of networks sampled to take part in tournament selection
#'@param chooseK Number of winners in tournament selection
#'@param term_cri Stress criteria (ln) for evolutionary algorithm terminating
#'@param sd_ini StDev of initialised networked values
#'@param use_cor Compare correlation matrix rather than the covariance matrix
#'@param prnt_out Print out progress showing stress for each generation
#'@param gmatrix G-matrix that the evolutionary algorithm will match
#'@return A list of eight elements that includes the following: (1) A vector of
#'input parameters, (2) the pre-specified covariance matrix, (3) matrix defining
#'the effects of loci values on the first layer of the network, (4) a three
#'dimensional array link the first network layer to trait values, (5) a matrix
#'of the marginal effect of each locus on each trait, (6) the mined covariance
#'structure, (7) all network values to be inserted into individual genomes, and
#'(8) the log stress of the mined matrix against the pre-specified matrix.
#'@examples
#'gmt       <- matrix(data = 0, nrow = 3, ncol = 3);
#'diag(gmt) <- 1;
#'mg        <- mine_gmatrix(gmatrix = gmt, loci = 6, layers = 3, indivs = 100, 
#'                          npsize = 100, max_gen = 2, prnt_out = FALSE);
#'@useDynLib resevol
#'@importFrom stats rnorm rpois runif rbinom
#'@importFrom utils read.csv write.csv
#'@export
mine_gmatrix <- function(loci     = 18, 
                         layers   = 6, 
                         indivs   = 1000,
                         npsize   = 2000, 
                         mu_pr    = 0.05, 
                         mu_sd    = 0.01,
                         max_gen  = 1000,
                         pr_cross = 0.05,
                         sampleK  = 40,
                         chooseK  = 4,
                         term_cri = -5.3,
                         sd_ini   = 0.1,
                         use_cor  = FALSE,
                         prnt_out = TRUE,
                         gmatrix){
    
    paras <- c(loci, layers, indivs, npsize, mu_pr, mu_sd, max_gen, pr_cross,
               sampleK, chooseK, term_cri, sd_ini, use_cor, as.numeric(prnt_out)
               );
    
    if(loci < 2 | loci %% 1 > 0){
        stop("ERROR: 'loci' needs to be an integer value above 1.")
    }
    if(layers < 1 | layers %% 1 > 0){
        stop("ERROR: 'layers' needs to be an integer value above 0.")
    }
    if(indivs < 10 | indivs %% 1 > 0){
        stop("ERROR: 'indivs' needs to be an integer value above 9.")
    }
    if(npsize < 10 | npsize %% 1 > 0){
        stop("ERROR: 'npsize' needs to be an integer value above 9.")
    }
    if(mu_pr > 1 | mu_pr < 0){
        stop("ERROR: 'mu_pr' must be a value between 0 and 1.")
    }
    if(mu_sd <= 0){
        stop("ERROR: 'mu_sd' must be a value greater than 0.")
    }
    if(max_gen < 1 | max_gen %% 1 > 0){
        stop("ERROR: 'max_gen' needs to be an integer value above 0.")
    }
    if(pr_cross > 1 | pr_cross < 0){
        stop("ERROR: 'pr_cross' must be a value between 0 and 1.")
    }
    if(sampleK < 1 | sampleK %% 1 > 0){
        stop("ERROR: 'sampleK' needs to be an integer value above 0.")
    }
    if(chooseK < 1 | chooseK %% 1 > 0){
        stop("ERROR: 'chooseK' needs to be an integer value above 0.")
    }
    if(sd_ini <= 0){
        stop("ERROR: 'sd_ini' must be a value greater than 0.")
    }
    if(dim(gmatrix)[1] != dim(gmatrix)[2]){
        stop("ERROR: 'gmatrix' needs to be a square matrix");
    }
    GOUT     <- run_mine_gmatrix(PARAS = paras, GMATRIX_c = gmatrix);
    return(GOUT); # Returns the values to produce the desired g-matrix
}

run_mine_gmatrix <- function(PARAS, GMATRIX_c){
    .Call("mine_gmatrix", PARAS, GMATRIX_c);
}