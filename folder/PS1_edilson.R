########################################### Question 6 ###########################################
# Importing data
PS1_data_edilson <- read.csv("~/Documentos/EPGE/Econometria_1/Listas/PS1_data_edilson.csv")
attach(PS1_data_edilson)

# Defining X matrix
n <- dim(PS1_data_edilson)[1]
X <- matrix(c(X1,rep(0,n),rep(0,n),X2),
            nrow=2*n, ncol=2)
Y <- c(Y1,Y2)

# Estimating beta_hat_mme
beta_hat_MME <- solve(t(X)%*%X)%*%t(X)%*%Y

########################################### Question 7 ###########################################
# Defining the objective function
objective_func <- function(b,y1,y2,x1,x2,M){
  v = c(sum((x1*(y1-x1*b[1]))), 
        sum((x2*(y2-x2*b[2]))),
        sum((x1*((y1-x1*b[1])^3))),
        sum((x2*((y2-x2*b[2])^3))))
  
  t(v)%*%M%*%v
}

# Optimizing
opt <- optim(par=c(1,1),
             objective_func,
             method = "BFGS",
             y1=Y1, y2=Y2, 
             x1=X1,x2=X2, M=diag(4))
beta_hat_GMM <- opt$par
beta_hat_GMM

########################################### Question 8 ###########################################
# Estimating G_0
jacobian_h <- array(dim = c(4,2,n))
G_hat <- matrix(0,nrow = 4,ncol = 2)
for (i in 1:n) {
  jacobian_h[,,i] <- -matrix(c(X1[i]^2,0,
                               0,X2[i]^2,
                               3*X1[i]^2*((Y1[i] -X1[i]*beta_hat_GMM[1])^2), 0,
                               0, 3*X2[i]^2*((Y2[i] -X2[i]*beta_hat_GMM[2])^2)),
                             byrow=TRUE,
                             ncol=2)
  
  G_hat <- (jacobian_h[,,i]+G_hat)
}
G_hat <- (1/n)*G_hat

# Estimating S_0
h_beta_GMM <- matrix(nrow=n,ncol=4)
S_hat <- matrix(0,nrow=4,ncol=4)
for (i in 1:n) {
  h_beta_GMM[i,] <- c(X1[i]*(Y1[i] - X1[i]*beta_hat_GMM[1]),
                      X2[i]*(Y2[i] - X2[i]*beta_hat_GMM[2]),
                      X1[i]*((Y1[i] - X1[i]*beta_hat_GMM[1])^3),
                      X2[i]*((Y2[i] - X2[i]*beta_hat_GMM[2])^3))
  S_hat <- S_hat + (h_beta_GMM[i,]%*%t(h_beta_GMM[i,]))
}
S_hat <- (1/n)*S_hat

# Calculating the asymptotic variance of beta_hat_GMM
V_hat <- (1/n)*solve(t(G_hat)%*%G_hat)%*%t(G_hat)%*%S_hat%*%G_hat%*%solve(t(G_hat)%*%G_hat)
V_hat

########################################### Question 9 ###########################################
# Finding the optimal weighting matrix
S_hat_inv <- solve(S_hat)
S_hat_inv

# Estimating beta_hat_OGMM
opt1 <- optim(par=c(1,1),
              objective_func,
              method = "BFGS",
              y1=Y1, y2=Y2, 
              x1=X1,x2=X2, M=S_hat_inv)
beta_hat_OGMM <- opt1$par
beta_hat_OGMM

########################################### Question 9 ###########################################
# Calculating G_tilde
jacobian_h_tilde <- array(dim = c(4,2,n))
G_tilde <- matrix(0,nrow = 4,ncol = 2)
for (i in 1:n) {
  jacobian_h_tilde[,,i] <- -matrix(c(X1[i]^2,0,
                                     0,X2[i]^2,
                                     3*X1[i]^2*((Y1[i] -X1[i]*beta_hat_OGMM[1])^2), 0,
                                     0, 3*X2[i]^2*((Y2[i] -X2[i]*beta_hat_OGMM[2])^2)),
                                   byrow=TRUE,
                                   ncol=2)
  
  G_tilde <- (jacobian_h_tilde[,,i]+G_tilde)
}
G_tilde <- (1/n)*G_tilde

# Estimating S_0
h_beta_OGMM <- matrix(nrow=n,ncol=4)
S_tilde <- matrix(0,nrow=4,ncol=4)
for (i in 1:n) {
  h_beta_OGMM[i,] <- c(X1[i]*(Y1[i] - X1[i]*beta_hat_OGMM[1]),
                      X2[i]*(Y2[i] - X2[i]*beta_hat_OGMM[2]),
                      X1[i]*((Y1[i] - X1[i]*beta_hat_OGMM[1])^3),
                      X2[i]*((Y2[i] - X2[i]*beta_hat_OGMM[2])^3))
  S_tilde <- S_tilde + (h_beta_OGMM[i,]%*%t(h_beta_OGMM[i,]))
}
S_tilde <- (1/n)*S_tilde

# Calculating the asymptotic variance of beta_hat_GMM
V_tilde <- (1/n)*solve(t(G_tilde)%*%solve(S_tilde)%*%G_tilde)
V_tilde
