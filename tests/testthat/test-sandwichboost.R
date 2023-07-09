
set.seed(111)

I = 50
N = 3
n_obs = N*I
g_0 = function(x) tanh(x[,1])
m_0 = function(x) cos(x[,1])
f_0 = function(x) 2 + 1*tanh(x[,1])
beta = 1
X = as.matrix(runif(n_obs,-5, 5))
V = as.matrix(rnorm(n_obs,0,1))
U = as.matrix(rnorm(n_obs,0,1))
epsilon = as.matrix(rep(0,n_obs))
xi = as.matrix(rep(0,n_obs))

for (i in seq_len(I)) {
  x = as.matrix(X[(N*(i-1)+1):(N*i),])
  Sigma0 = matrix(0,N,N)
  f_0x <- f_0(x)
  for (j in seq_len(N)) {
    for (jj in seq_len(N)) {
      Sigma0[j,jj] <- 0.8^(abs(j-jj))*f_0x[j]*f_0x[jj]
    }
  }
  epsilon[(N*(i-1)+1):(N*i),] = ((expm::sqrtm(Sigma0)))%*%U[(N*(i-1)+1):(N*i),]
}

D = m_0(X) + xi
Y = beta*D + g_0(X) + epsilon
id =rep(1:I, each=N)

rdf <- data.frame(X=X, D=D, Y=Y, id=id)

test_that("sandwich boosting works for 1 dimensional X & gams", {
  tst <- PLR_sandboost(l_formula = Y ~ s(X,bs="cr"), l_learner="gam", m_formula = D ~ s(X,bs="cr"), m_learner="gam", s_formula = u_s ~ s(X,bs="cr"), s_learner="gam", data=rdf, proxyCCF="equicorr")
  expect_type(tst, "list")
})

test_that("sandwich boosting works for 1 dimensional X & random forest fits", {
  tst <- PLR_sandboost(l_formula = Y ~ X, l_learner="randomforest", m_formula = D ~ X, m_learner="lm", s_formula = u_s ~ X, s_learner="randomforest", data=rdf, K=2L, proxyCCF="equicorr", m_stop=20)
  expect_type(tst, "list")
})

test_that("sandwich boosting works for 1 dimensional X & a mix of lm/gam", {
  tst <- PLR_sandboost(l_formula = Y ~ X, l_learner="lm", m_formula = D ~ X, m_learner="lm", s_formula = u_s ~ s(X,bs="cr"), s_learner="gam", data=rdf, K=2L, proxyCCF="autoreg", init_CCF="mem", m_stop=20)
  expect_type(tst, "list")
})


g_0 = function(x) tanh(x[,1]) + cosh(x[,3])
m_0 = function(x) cos(x[,1]) + 2*x[,2]
f_0 = function(x) 2 + 1*tanh(x[,1])
beta = 1
num_X = 10
X = matrix(runif(n_obs*num_X,-5, 5), nrow=n_obs, ncol=num_X)
V = as.matrix(rnorm(n_obs,0,1))
U = as.matrix(rnorm(n_obs,0,1))
epsilon = as.matrix(rep(0,n_obs))
xi = as.matrix(rep(0,n_obs))

for (i in seq_len(I)) {
  x = as.matrix(X[(N*(i-1)+1):(N*i),])
  Sigma0 = matrix(0,N,N)
  f_0x <- f_0(x)
  for (j in seq_len(N)) {
    for (jj in seq_len(N)) {
      Sigma0[j,jj] <- 0.8^(abs(j-jj))*f_0x[j]*f_0x[jj]
    }
  }
  epsilon[(N*(i-1)+1):(N*i),] = ((expm::sqrtm(Sigma0)))%*%U[(N*(i-1)+1):(N*i),]
}

D = m_0(X) + xi
Y = beta*D + g_0(X) + epsilon
id =rep(1:I, each=N)

rdf <- data.frame(X=X, D=D, Y=Y, id=id)

test_that("sandwich boosting works for multidimensional X & gams/randomforests", {
  tst <- PLR_sandboost(l_formula = Y ~ s(X.1,bs="cr"), l_learner="gam", m_formula = D ~ s(X.2,bs="cr"), m_learner="gam", s_formula = u_s ~ X.3 + X.5, s_learner="randomforest", data=rdf, m_stop=20)
  expect_type(tst, "list")
})

test_that("sandwich boosting works for multidimensional X & a mix of lasso/ridge/gam", {
  tst <- PLR_sandboost(l_formula = Y ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, l_learner="lasso", m_formula = D ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, m_learner="ridge", s_formula = u_s ~ s(X.1,bs="cr") + s(X.2,bs="cr"), s_learner="gam", data=rdf, variable_steps = TRUE)
  expect_type(tst, "list")
})
