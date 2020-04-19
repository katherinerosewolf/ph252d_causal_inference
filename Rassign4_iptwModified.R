set.seed(1)

# The true value of Qbar0(A,W) = E_0[Y|A,W]
Qbar0 = function(A,W){1000 + plogis(W*A)}
# The true value of g0(1|W) = Pr(A=1|W)
g0 = function(W){0.2 + 0.6*W}

# A function which returns a data frame with n i.i.d. observations from P0
gen.data = function(n){
	W = rbinom(n,1,1/2)
	A = rbinom(n,1,g0(W))
	Y = 1000+rbinom(n,1,Qbar0(A,W)-1000)
	return(data.frame(W=W,A=A,Y=Y))
}

# samples size
n= 1000
# Number of Monte Carlo draws
R = 2000
# Matrix of estimates from IPTW, modified Horvitz-Thompson, and my.est
est = matrix(NA, 
             nrow = R,
             ncol = 5)
colnames(est) = c('IPTW', 'Modifed HT', 'simple.mean', 'my.est', 'my.est.2')
for(r in 1:R){
	# Generate data with sample size
	ObsData = gen.data(n)
	W = ObsData$W
	A = ObsData$A
	Y = ObsData$Y
	# IPTW estimate
	IPTW.est = mean(A * Y/g0(W))
	# Modified Horvitz-Thompson estimate
	HT.est = mean(A * Y/g0(W))/mean(A/g0(W))
	# You should replace the NA below with your own estimate
	simple.mean = mean(Y)
	# scale outcome
	l <- min(Y)
	u <- max(Y)
	l.w <- min(A/g0(W))
	u.w <- max(A/g0(W))
	scale.w <- (A/g0(W) - l.w)/(u.w - l.w)
	scale_Y <- (Y - l)/(u - l)
	iptw <- A * Y/g0(W)
	scale_iptw <- (A * Y/g0(W) - min(A * Y/g0(W)))/(max(A * Y/g0(W)) - min(A * Y/g0(W)))
	my.est = mean(A * (Y-mean(Y))/g0(W) + mean(Y))
	my.est.2 = mean(scale_iptw) + min(Y)
	# Put the estimates into the est matrix
	est[r,] = c(IPTW.est,HT.est,simple.mean,my.est,my.est.2)
}

# Calculate the true value of EE[Y|A=1,W]
truth = 1/2*(Qbar0(1,0) + Qbar0(1,1))
# note: we know P_0(W=1)=0.5
truth

# Calculate the estimated bias, variance, and MSE
est.bias = colMeans(est) - truth
est.var = apply(est,2,var)
est.mse = est.bias^2 + est.var

# Only can report estimated bias/variance/MSE because only took finitely many Monte Carlo draws (2000)
print('The estimators have (estimated) bias:')
print(est.bias)
print('The estimators have (estimated) variance:')
print(est.var)
print('The estimators have (estimated) MSE:')
print(est.mse)

