### Directories:
work_path <- getwd()
source_path <- paste0(work_path,"/sources")
figure_path <- paste0(work_path, "/figures")
### Input the functions needed:
source(paste0(source_path,"/loads.R"))
source(paste0(source_path,"/functions.R"))
library(latex2exp)
font_size <- 1.3
axis_size <- 1.3


################################################################################# 
################################################################################# 
##################### Show an variogram cor(5,t): ###############################
################################################################################# 
################################################################################# 
p = 1
s <- 5
t_vec <- seq(0,15, by = 1)
cov_st <- c(0)
appr_cov_st1 <- c(0)
appr_cov_st2 <- c(0)
appr_cov_st3 <- c(0)
appr_cov_st4 <- c(0)
k <- c(5,10,30,100)
for (i in 2:length(t_vec)) {
  t <- t_vec[i]
  x <- sort(c(s,t))
  Sigma_IWP <- compute_Wp_cov(x, p = p)
  cov_st[i] <- Sigma_IWP[1,2]/sqrt(Sigma_IWP[1,1]*Sigma_IWP[2,2])
  knots <- seq(0,15, length.out = (k[1]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st1[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
  knots <- seq(0,15, length.out = (k[2]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st2[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
  knots <- seq(0,15, length.out = (k[3]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st3[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
  knots <- seq(0,15, length.out = (k[4]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st4[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
}
pdf(file = paste0(figure_path, "/vari_p1.pdf"), width = 5, height = 5)
par(mar = c(2.5, 3, 0.1, 0.1), mgp = c(1.4, 0.5, 0))
plot(cov_st~t_vec, type = "p", ylab = TeX(r'($\rho(5,x)$)'), xlab = "x", cex=2, ylim = c(0,1), cex.lab=font_size, cex.axis=axis_size,
)
lines(appr_cov_st1~t_vec, lty = 'dotted', col = "green", lwd = 4)
lines(appr_cov_st2~t_vec, lty = 'dashed', col = "brown", lwd = 2)
lines(appr_cov_st3~t_vec, lty = 'dotdash', col = "red", lwd = 2)
lines(appr_cov_st4~t_vec, lty = 'solid', col = "purple", lwd = 2)
legend(8, 0.6, legend=c("True", "k = 5", "k = 10", "k = 30", "k = 100"),
       col=c("black", "green", "brown", "red", "purple"), lty= c(NA, "dotted", "dashed", "dotdash", "solid"), 
       cex=1.5, box.lty = 0, lwd = 2,
       pch = c(1, NA, NA, NA, NA)
)
dev.off()

p = 2
s <- 5
t_vec <- seq(0,15, by = 1)
cov_st <- c(0)
appr_cov_st1 <- c(0)
appr_cov_st2 <- c(0)
appr_cov_st3 <- c(0)
appr_cov_st4 <- c(0)
k <- c(5,10,30,100)
for (i in 2:length(t_vec)) {
  t <- t_vec[i]
  x <- sort(c(s,t))
  Sigma_IWP <- compute_Wp_cov(x, p = p)
  cov_st[i] <- Sigma_IWP[1,2]/sqrt(Sigma_IWP[1,1]*Sigma_IWP[2,2])
  knots <- seq(0,15, length.out = (k[1]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st1[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
  knots <- seq(0,15, length.out = (k[2]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st2[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
  knots <- seq(0,15, length.out = (k[3]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st3[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
  knots <- seq(0,15, length.out = (k[4]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st4[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
}
pdf(paste0(figure_path, "/vari_p2.pdf"), width = 5, height = 5)
par(mar = c(2.5, 3, 0.1, 0.1), mgp = c(1.4, 0.5, 0))
plot(cov_st~t_vec, type = "p", ylab = TeX(r'($\rho(5,x)$)'), xlab = "x", cex = 2, ylim = c(0,1), xlim = c(0,15), cex.lab=font_size, cex.axis=axis_size,
)
lines(appr_cov_st1~t_vec, lty = 'dotted', col = "green", lwd = 4)
lines(appr_cov_st2~t_vec, lty = 'dashed', col = "brown", lwd = 2)
lines(appr_cov_st3~t_vec, lty = 'dotdash', col = "red", lwd = 2)
lines(appr_cov_st4~t_vec, lty = 'solid', col = "purple", lwd = 2)
dev.off()

p = 3
s <- 5
t_vec <- seq(0,15, by = 1)
cov_st <- c(0)
appr_cov_st1 <- c(0)
appr_cov_st2 <- c(0)
appr_cov_st3 <- c(0)
appr_cov_st4 <- c(0)
k <- c(5,10,30,100)
for (i in 2:length(t_vec)) {
  t <- t_vec[i]
  x <- sort(c(s,t))
  Sigma_IWP <- compute_Wp_cov(x, p = p)
  cov_st[i] <- Sigma_IWP[1,2]/sqrt(Sigma_IWP[1,1]*Sigma_IWP[2,2])
  knots <- seq(0,15, length.out = (k[1]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st1[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
  knots <- seq(0,15, length.out = (k[2]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st2[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
  knots <- seq(0,15, length.out = (k[3]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st3[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
  knots <- seq(0,15, length.out = (k[4]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st4[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
}
pdf(paste0(figure_path,"/vari_p3.pdf"), width = 5, height = 5)
par(mar = c(2.5, 3, 0.1, 0.1), mgp = c(1.4, 0.5, 0))
plot(cov_st~t_vec, type = "p", ylab = TeX(r'($\rho(5,x)$)'), xlab = "x", cex = 2, ylim = c(0,1), cex.lab=font_size, cex.axis=axis_size,
)
lines(appr_cov_st1~t_vec, lty = 'dotted', col = "green", lwd = 4)
lines(appr_cov_st2~t_vec, lty = 'dashed', col = "brown", lwd = 2)
lines(appr_cov_st3~t_vec, lty = 'dotdash', col = "red", lwd = 2)
lines(appr_cov_st4~t_vec, lty = 'solid', col = "purple", lwd = 2)
dev.off()

p = 4
s <- 5
t_vec <- seq(0,15, by = 1)
cov_st <- c(0)
appr_cov_st1 <- c(0)
appr_cov_st2 <- c(0)
appr_cov_st3 <- c(0)
appr_cov_st4 <- c(0)
k <- c(5,10,30,100)
for (i in 2:length(t_vec)) {
  t <- t_vec[i]
  x <- sort(c(s,t))
  Sigma_IWP <- compute_Wp_cov(x, p = p)
  cov_st[i] <- Sigma_IWP[1,2]/sqrt(Sigma_IWP[1,1]*Sigma_IWP[2,2])
  knots <- seq(0,15, length.out = (k[1]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st1[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
  knots <- seq(0,15, length.out = (k[2]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st2[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
  knots <- seq(0,15, length.out = (k[3]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st3[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
  knots <- seq(0,15, length.out = (k[4]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cov_st4[i] <- Sigma_OS[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS[2,2])
}
pdf(paste0(figure_path, "/vari_p4.pdf"), width = 5, height = 5)
par(mar = c(2.5, 3, 0.1, 0.1), mgp = c(1.4, 0.5, 0))
plot(cov_st~t_vec, type = "p", ylab = TeX(r'($\rho(5,x)$)'), xlab = "x", cex = 2, ylim = c(0,1), cex.lab=font_size, cex.axis=axis_size,
)
lines(appr_cov_st1~t_vec, lty = 'dotted', col = "green", lwd = 4)
lines(appr_cov_st2~t_vec, lty = 'dashed', col = "brown", lwd = 2)
lines(appr_cov_st3~t_vec, lty = 'dotdash', col = "red", lwd = 2)
lines(appr_cov_st4~t_vec, lty = 'solid', col = "purple", lwd = 2)
dev.off()





################################################################################# 
################################################################################# 
##################### Do the same for cross-correlation computation: ############
################################################################################# 
################################################################################# 

### Compute the correlation bewteen Wp(y) and Wq(x)

### Fixed the location y = 5:
y = 5

### Compute the exact cross-correlation
p = 2
q <- p-1
x = 2
C <- solve(Compute_Aug_Wp_Prec(svec = sort(c(x,y)), p = 2))
exact_cross_cov <- C[1,(2*p)]/sqrt(C[(2*p),(2*p)] * C[1,1])

k <- 300
knots <- seq(0,15, length.out = k)
OS_design <- local_poly(x = knots, refined_x = sort(c(x,y)), p = p)
OS_deriv_design <- local_poly(x = knots, refined_x = sort(c(x,y)), p = q)
OS_prec <- compute_weights_precision(knots)
Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)

approxi_cross_cov <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
mar.default <- c(5.1, 4.1, 4.1, 2.1)

### Plotting:
p = 2
q <- p - 1
s <- 5
t_vec <- seq(0,15, by = 1)
t_vec <- t_vec[t_vec!=5]
cross_cov_st <- c(0)
appr_cross_cov_st1 <- c(0)
appr_cross_cov_st2 <- c(0)
appr_cross_cov_st3 <- c(0)
appr_cross_cov_st4 <- c(0)
k <- c(5,10,30,100)
for (i in 2:length(t_vec)) {
  t <- t_vec[i]
  x <- sort(c(s,t))
  C <- solve(Compute_Aug_Wp_Prec(svec = x, p = p))
  cross_cov_st[i] <- C[1,(p+2)]/sqrt(C[1,1] * C[(p+2), (p+2)])
  knots <- seq(0,15, length.out = (k[1]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st1[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
  knots <- seq(0,15, length.out = (k[2]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st2[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
  knots <- seq(0,15, length.out = (k[3]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st3[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
  knots <- seq(0,15, length.out = (k[4]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st4[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
}
pdf(file = paste0(figure_path, "/cross_p2.pdf"), width = 5, height = 5)
par(mar = c(2.5, 3, 0.1, 0.1), mgp = c(1.4, 0.5, 0))
plot(cross_cov_st~t_vec, type = "p", ylab = TeX(r'($\rho^{(0,1)}(5,x)$)'), xlab = "x", cex = 2, ylim = c(0,1), cex.lab=font_size, cex.axis=axis_size,
)
lines(appr_cross_cov_st1~t_vec, lty = 'dotted', col = "green", lwd = 4)
lines(appr_cross_cov_st2~t_vec, lty = 'dashed', col = "brown", lwd = 2)
lines(appr_cross_cov_st3~t_vec, lty = 'dotdash', col = "red", lwd = 2)
lines(appr_cross_cov_st4~t_vec, lty = 'solid', col = "purple", lwd = 2)
legend(8, 0.4, legend=c("True","k = 5", "k = 10", "k = 30", "k = 100"),
       col=c("black", "green", "brown", "red", "purple"), lty= c(NA, "dotted", "dashed", "dotdash", "solid"), 
       cex=1.5, box.lty = 0, lwd = 2,
       pch = c(1, NA, NA, NA, NA))
dev.off()

### Plotting:
p = 3
q <- p - 1
s <- 5
t_vec <- seq(0,15, by = 1)
t_vec <- t_vec[t_vec!=5]
cross_cov_st <- c(0)
appr_cross_cov_st1 <- c(0)
appr_cross_cov_st2 <- c(0)
appr_cross_cov_st3 <- c(0)
appr_cross_cov_st4 <- c(0)

k <- c(5,10,30,100)
for (i in 2:length(t_vec)) {
  t <- t_vec[i]
  x <- sort(c(s,t))
  C <- solve(Compute_Aug_Wp_Prec(svec = x, p = p))
  cross_cov_st[i] <- C[1,(p+2)]/sqrt(C[1,1] * C[(p+2), (p+2)])
  knots <- seq(0,15, length.out = (k[1]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st1[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
  knots <- seq(0,15, length.out = (k[2]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st2[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
  knots <- seq(0,15, length.out = (k[3]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st3[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
  knots <- seq(0,15, length.out = (k[4]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st4[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
}
pdf(paste0(figure_path, "/cross_p3.pdf"), width = 5, height = 5)
par(mar = c(2.5, 3, 0.1, 0.1), mgp = c(1.4, 0.5, 0))
plot(cross_cov_st~t_vec, type = "p", ylab = TeX(r'($\rho^{(0,1)}(5,x)$)'), xlab = "x", cex = 2, ylim = c(0,1), cex.lab=font_size, cex.axis=axis_size,
)
lines(appr_cross_cov_st1~t_vec, lty = 'dotted', col = "green", lwd = 4)
lines(appr_cross_cov_st2~t_vec, lty = 'dashed', col = "brown", lwd = 2)
lines(appr_cross_cov_st3~t_vec, lty = 'dotdash', col = "red", lwd = 2)
lines(appr_cross_cov_st4~t_vec, lty = 'solid', col = "purple", lwd = 2)
dev.off()

### Plotting:
p = 4
q <- p - 1
s <- 5
t_vec <- seq(0,15, by = 1)
t_vec <- t_vec[t_vec!=5]
cross_cov_st <- c(0)
appr_cross_cov_st1 <- c(0)
appr_cross_cov_st2 <- c(0)
appr_cross_cov_st3 <- c(0)
appr_cross_cov_st4 <- c(0)
k <- c(5,10,30,100)
for (i in 2:length(t_vec)) {
  t <- t_vec[i]
  x <- sort(c(s,t))
  C <- solve(Compute_Aug_Wp_Prec(svec = x, p = p))
  cross_cov_st[i] <- C[1,(p+2)]/sqrt(C[1,1] * C[(p+2), (p+2)])
  knots <- seq(0,15, length.out = (k[1]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st1[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
  knots <- seq(0,15, length.out = (k[2]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st2[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
  knots <- seq(0,15, length.out = (k[3]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st3[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
  knots <- seq(0,15, length.out = (k[4]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st4[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
}
pdf(paste0(figure_path, "/cross_p4.pdf"), width = 5, height = 5)
par(mar = c(2.5, 3, 0.1, 0.1), mgp = c(1.4, 0.5, 0))
plot(cross_cov_st~t_vec, type = "p", ylab = TeX(r'($\rho^{(0,1)}(5,x)$)'), xlab = "x", cex = 2, ylim = c(0,1), cex.lab=font_size, cex.axis=axis_size,
)
lines(appr_cross_cov_st1~t_vec, lty = 'dotted', col = "green", lwd = 4)
lines(appr_cross_cov_st2~t_vec, lty = 'dashed', col = "brown", lwd = 2)
lines(appr_cross_cov_st3~t_vec, lty = 'dotdash', col = "red", lwd = 2)
lines(appr_cross_cov_st4~t_vec, lty = 'solid', col = "purple", lwd = 2)
dev.off()

### Plotting:
p = 4
q <- p - 2
s <- 5
t_vec <- seq(0,15, by = 1)
t_vec <- t_vec[t_vec!=5]
cross_cov_st <- c(0)
appr_cross_cov_st1 <- c(0)
appr_cross_cov_st2 <- c(0)
appr_cross_cov_st3 <- c(0)
appr_cross_cov_st4 <- c(0)
k <- c(5,10,30,100)
for (i in 2:length(t_vec)) {
  t <- t_vec[i]
  x <- sort(c(s,t))
  C <- solve(Compute_Aug_Wp_Prec(svec = x, p = p))
  cross_cov_st[i] <- C[1,(p+3)]/sqrt(C[1,1] * C[(p+3), (p+3)])
  knots <- seq(0,15, length.out = (k[1]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st1[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
  knots <- seq(0,15, length.out = (k[2]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st2[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
  knots <- seq(0,15, length.out = (k[3]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st3[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
  knots <- seq(0,15, length.out = (k[4]))
  OS_design <- local_poly(x = knots, refined_x = x, p = p)
  OS_deriv_design <- local_poly(x = knots, refined_x = x, p = q)
  OS_prec <- compute_weights_precision(knots)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_cross <- OS_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS_deriv <- OS_deriv_design %*% solve(OS_prec) %*% t(OS_deriv_design)
  Sigma_OS <- OS_design %*% solve(OS_prec) %*% t(OS_design)
  appr_cross_cov_st4[i] <- Sigma_OS_cross[1,2]/sqrt(Sigma_OS[1,1]*Sigma_OS_deriv[2,2])
}
pdf(paste0(figure_path, "/cross_p4_2driv.pdf"), width = 5, height = 5)
par(mar = c(2.5, 3, 0.1, 0.1), mgp = c(1.4, 0.5, 0))
plot(cross_cov_st~t_vec, type = "p", ylab = TeX(r'($\rho^{(0,2)}(5,x)$)'), xlab = "x", cex = 2, ylim = c(0,1), cex.lab=font_size, cex.axis=axis_size,
)
lines(appr_cross_cov_st1~t_vec, lty = 'dotted', col = "green", lwd = 4)
lines(appr_cross_cov_st2~t_vec, lty = 'dashed', col = "brown", lwd = 2)
lines(appr_cross_cov_st3~t_vec, lty = 'dotdash', col = "red", lwd = 2)
lines(appr_cross_cov_st4~t_vec, lty = 'solid', col = "purple", lwd = 2)
dev.off()

