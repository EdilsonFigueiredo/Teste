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




