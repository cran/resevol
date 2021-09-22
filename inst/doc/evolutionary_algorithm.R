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

## -----------------------------------------------------------------------------
arrows_1_dat <- rnorm(n = 4 * 12, mean = 0, sd = 0.1);
arrows_1_mat <- matrix(data = arrows_1_dat, nrow = 12, ncol = 4);
print(arrows_1_mat);

## -----------------------------------------------------------------------------
loci <- rnorm(n = 12, mean = 0, sd = 1);

## -----------------------------------------------------------------------------
print(loci %*% arrows_1_mat);

## -----------------------------------------------------------------------------
arrows_2_dat <- rnorm(n = 4 * 4, mean = 0, sd = 0.1);
arrows_2_mat <- matrix(data = arrows_2_dat, nrow = 4, ncol = 4);
print(arrows_2_mat);

## -----------------------------------------------------------------------------
print(loci %*% arrows_1_mat %*% arrows_2_mat);

## ---- echo=FALSE, fig.height=1.75, fig.width=6, fig.cap = "Figure 2: Conceptual overview of the evolutionary algorithm used in the resevol package."----
mbox <- function(x0, x1, y0, y1){
    xx <- seq(from=x0, to=x1, length.out = 100);
    yy <- seq(from=y0, to=y1, length.out = 100);
    xd <- c(rep(x0, 100), xx, rep(x1,100), rev(xx));
    yd <- c(yy, rep(y1,100), rev(yy), rep(y0, 100));
    return(list(x=xd, y=yd));
}
par(mar=c(0,0,0,0));
# ===============================================================
plot(x=0, y=0, type="n", xlim=c(0,85), ylim=c(50,95), xaxt="n", yaxt="n",
     xlab="",ylab="");
ibox <- mbox(x0 = 0,  x1 = 10, y0 = 90, y1 = 70);
polygon(x = ibox$x, y = ibox$y, lwd = 3, border = "black", col = "white");
cbox <- mbox(x0 = 15, x1 = 25, y0 = 90, y1 = 70);
polygon(x = cbox$x, y = cbox$y, lwd = 3, border = "black", col = "white");
ubox <- mbox(x0 = 30, x1 = 40, y0 = 90, y1 = 70);
polygon(x = ubox$x, y = ubox$y, lwd = 3, border = "black", col = "white");
nbox <- mbox(x0 = 45, x1 = 55, y0 = 90, y1 = 70);
polygon(x = nbox$x, y = nbox$y, lwd = 3, border = "black", col = "white");
fbox <- mbox(x0 = 60, x1 = 70, y0 = 90, y1 = 70);
polygon(x = fbox$x, y = fbox$y, lwd = 3, border = "black", col = "white");
tbox <- mbox(x0 = 75, x1 = 85, y0 = 90, y1 = 70);
polygon(x = tbox$x, y = tbox$y, lwd = 3, border = "black", col = "white");
text(x=5,  y=80, labels="Initialisation", col="black", cex = 0.5);
text(x=20, y=80, labels="Crossover", col="black", cex = 0.5);
text(x=35, y=80, labels="Mutation", col="black", cex = 0.5);
text(x=50, y=82, labels="Fitness", col="black", cex = 0.5);
text(x=50, y=78, labels="evaluation", col="black", cex = 0.5);
text(x=65, y=82, labels="Tournament", col="black", cex = 0.5);
text(x=65, y=78, labels="selection", col="black", cex = 0.5);
text(x=80, y=80, labels="Replacement", col="black", cex = 0.5);
arrows(x0=10, x1=14, y0=80, y1=80, lwd=2, length=0.10);
arrows(x0=25, x1=29, y0=80, y1=80, lwd=2, length=0.10);
arrows(x0=40, x1=44, y0=80, y1=80, lwd=2, length=0.10);
arrows(x0=55, x1=59, y0=80, y1=80, lwd=2, length=0.10);
arrows(x0=70, x1=74, y0=80, y1=80, lwd=2, length=0.10);
arrows(x0 = 80, x1 = 80, y0 = 70, y1 = 68, lwd = 2, length = 0);
arrows(x0 = 80, x1 = 20, y0 = 65, y1 = 65, lwd = 2, length = 0);
arrows(x0 = 20, x1 = 20, y0 = 65, y1 = 70, lwd = 2, length = 0.05);
text(x=73, y=63, labels="No", col="black", cex = 0.5);
rbox <- mbox(x0 = 75, x1 = 85, y0 = 62, y1 = 68);
polygon(x = rbox$x, y = rbox$y, lwd = 3, border = "black", col = "black");
text(x=80, y=65, labels="Termination?", col="white", cex = 0.5);
fbox <- mbox(x0 = 60, x1 = 70, y0 = 60, y1 = 50);
polygon(x = fbox$x, y = fbox$y, lwd = 3, border = "black", col = "white");
text(x=65, y=57, labels="Individual", col="black", cex = 0.5);
text(x=65, y=53, labels="network", col="black", cex = 0.5);
arrows(x0 = 80, x1 = 80, y0 = 62, y1 = 55, lwd = 2, length = 0);
arrows(x0 = 80, x1 = 70, y0 = 55, y1 = 55, lwd = 2, length = 0.05);
text(x=73, y=53, labels="Yes", col="black", cex = 0.5);

## ---- echo = FALSE------------------------------------------------------------
par(oldpar);

