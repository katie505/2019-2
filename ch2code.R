#10월1일 통계계산
# Chapter 2

Bisection = function(x0, x1, epsilon = 1e-5) #이분법
{
  fx0 = f(x0)
  fx1 = f(x1)
  if (fx0 * fx1 >0)  
    return("wrong initial values")
  error = abs(x1 - x0)
  N = 1
  while (error > epsilon)
  {
    N = N + 1
    error = error / 2
    x2 = (x0 + x1) / 2
    fx2 = f(x2)
    if (fx0 * fx2 < 0)
    {
      x1 = x2; fx1 = fx2
    } else
    {
      x0 = x2; fx0 = fx2
    }
  }
  
  return(list(x = x2, n = N))
}


Newton = function(x0, epsilon = 1e-5, n = 100)
{
  e = 1
  N = 1
  d = epsilon
  while (e > epsilon)
  {
    N = N + 1
    if (N > n) 
      return("not converge after 100 iterations")
    x1 = x0 - f(x0) * d / (f(x0 + d) - f(x0))
    e = abs(x1 - x0)
    x0 = x1
  }
  
  return(list(x = x1, n = N))
}


Integral = function(a, b, n)
{
  integral = 0
  h = (b - a) / n
  for (i in 1:n)
    integral = integral + h * f(a + (i-1/2) * h)
  
  return(integral)
}

Trapezoid = function(a, b, n = 50)
{
  h = (b - a) / n
  integral = (f(a) + f(b)) / 2
  
  x = a
  n1 = n - 1
  for (i in 1:n1)
  {
    x = x + h
    integral = integral + f(x)
  }
  integral = integral * h
  
  return(integral)
}

Simpson = function(a, b, n = 12)
{
  h = (b - a) / n
  integral = f(a) + f(b)
  x2 = a
  x3 = a + h
  even = 0
  odd = f(x3)
  h2 = 2 * h
  n1 = n / 2 - 1
  for (i in 1:n1)
  {
    x2 = x2 + h2
    x3 = x3 + h2
    even = even + f(x2)
    odd = odd + f(x3)
  }
  integral = (integral + 4 * odd + 2 * even) * h / 3
  
  return(integral)
}



# Example 2-1
f = function(x) {x^2-3}
result = Bisection(1,2)
result

Newton(1)


result2 = Newton(1)
abs(sqrt(3)-result2$x) /*참값과의 차이*/

# Example 2-2 #정규분포를 -1,1 구간에서 적분
f = function(x) dnorm(x)
Trapezoid(-1,1) #사다리꼴적분
Simpson(-1,1)
2*(pnorm(1)-0.5)

# Example 2-3

Trapezoid(3,4)
Simpson(3,4)
pnorm(4)-pnorm(3)

Trapezoid(3, 4, n=100)
Simpson(3, 4, n=24)


# Example 2-4
zq = function(p, x0=0, epsilon = 1e-5, n=100) {
  f = function(x) dnorm(x)
  F = function(x){Simpson(-4, x, n=24)-p}
  e=1
  N=1
  while (e > epsilon) {
    N = N+1
    if (N>n) return("not converge after 100 iterations")
    x1 = x0 - F(x0) / f(x0)
    e = abs(x1-x0)
    x0 = x1
  }
  return(list(x1, N))
}

zq(0.9)
qnorm(0.9)


## Optimization

f = function(x)
{
  f = (x[1] - 1)^2 + (x[2] - 1)^2 - x[1] * x[2]
}

df = function(x)
{
  df1 = 2 * (x[1] - 1) - x[2]
  df2 = 2 * (x[2] - 1) - x[1]
  df = c(df1, df2)
  return(df)
}

Norm = function(u)
{
  return(sqrt(sum(u^2)))
}

# Steepest decscent method
m = 100
par(mfrow=c(1,2), pty="s")
x1 = x2 = seq(-10.5, 10.5, length=m)
xg = expand.grid(x1, x2)
z = matrix(apply(xg, 1, f), m, m)
xh = NULL; fh = NULL
x0 = c(-10, -3); fx0 = f(x0); ni = 0
for (i in 1:10)
{  
  xh = rbind(xh, x0); fh = c(fh, fx0); ni = ni+1
  cat("iteration=", round(i,2))
  cat("  x0=", round(x0,2), "  f(x0)=", round(f(x0),3), "\n")
  d = df(x0)
  for (iters in 1:20)
  {
    x = x0 - d; fx = f(x)
    if (fx < fx0) break
    d = d / 2
  }
  x0 = x; fx0 = fx
}
contour(x1, x2, z, levels=round(fh, 2))
for (i in 1:(ni-1))
{
  points(xh[i,1], xh[i,2], pch=as.character(i))
  x1=xh[i,1]; y1=xh[i,2]; x2=xh[i+1,1]; y2=xh[i+1,2]
  arrows(x1, y1, x2, y2, length=0.1, col="red", lwd=0.5)
}
points(xh[ni,1], xh[ni,2], pch=as.character(ni))


# Newton-Raphson method
x1 = x2 = seq(-10.5, 10.5, length=m)
xg = expand.grid(x1, x2)
z = matrix(apply(xg, 1, f), m, m)
xh = NULL; fh = NULL
x0 = c(-10, -3); fx0 = f(x0); ni = 0
df2 = matrix(c(2,-1,-1,2),2,2); v = solve(df2)
for (i in 1:10)
{  
  xh = rbind(xh, as.vector(x0)); fh = c(fh, fx0); ni = ni+1
  cat("iteration=", round(i,2))
  cat("  x0=", round(x0,2), "  f(x0)=", round(f(x0),3), "\n")
  #   d = df(x0)
  d = v %*% df(x0)
  for (iters in 1:20)
  {
    x = x0 - d; fx = f(x)
    if (fx < fx0) break
    d = d / 2
  }
  if (abs(fx-fx0) < 1e-5) break
  x0 = x; fx0 = fx
}
contour(x1, x2, z, levels=round(fh, 2))
for (i in 1:(ni-1))
{
  points(xh[i,1], xh[i,2], pch=as.character(i))
  x1=xh[i,1]; y1=xh[i,2]; x2=xh[i+1,1]; y2=xh[i+1,2]
  arrows(x1, y1, x2, y2, length=0.1, col="red", lwd=0.5)
}
points(xh[ni,1], xh[ni,2], pch=as.character(ni))

## lp
library(lpSolve)
eg.lp = lp(objective.in=c(5,8),
           const.mat=matrix(c(1,1,1,2),nrow=2),
           const.rhs=c(2,3),
           const.dir=c("<=","="), direction="max")
eg.lp$solution

# degeneracy
degen.lp = lp(objective.in = c(3,1),
              const.mat = matrix(c(1,1,1,4,1,2,3,1), nrow=4),
              const.rhs = c(2,3,4,4), const.dir = rep(">=",4))
degen.lp
degen.lp$solution	# still finds the optimum

# infeasibility
eg.lp = lp(objective.in=c(5,8),
           const.mat=matrix(c(1,1,1,1),nrow=2),
           const.rhs=c(2,1),
           const.dir=c(">=","<="))
eg.lp

# unboundedness
eg.lp = lp(objective.in=c(5,8),
           const.mat=matrix(c(1,1,1,2),nrow=2),
           const.rhs=c(2,3),
           const.dir=c(">=",">="), direction="max")
eg.lp


## QP
library(quadprog)
x = c(.45, .08, -1.08, .92, 1.65, .53, .52, -2.15, -2.20,
      -.32, -1.87, -.16, -.19, -.98, -.20, .67, .08, .38,
      .76, -.78)
y = c(1.26, .58, -1, 1.07, 1.28, -.33, .68, -2.22, -1.82,
      -1.17, -1.54, .35, -.23, -1.53, .16, .91, .22, .44,
      .98, -.98)
X = cbind(rep(1,20), x)
XX = t(X) %*% X
Xy = t(X) %*% y
A = matrix(c(0, 1), ncol = 1)
b = 1
solve.QP(Dmat = XX, dvec = Xy, Amat = A, bvec = b)

# optimal portfolio for investing in 3 stocks
# beta_i : ratio of investment, nonnegative and sum to 1
# x: daily return for stocks: 0.002, 0.005, 0.01
# D: variability in the returns (covariance)
A = cbind(rep(1,3), diag(rep(1,3)))
D = matrix(c(.01,.002,.002,.002,.01,.002,.002,.002,.01), nrow=3)
x = c(.002,.005,.01)
b = c(1,0,0,0)
solve.QP(2*D, x, A, b, meq=1)

# optimal strategy: invest 10.4%, 29.2%, 60.4% for stocks 1,2,3
# optimal value is 0.002



### R functions(실제로 r코드에서 이분법이 쓰일 때)
## finding zeros of a function (based on bisection)
f = function(x, a) x^(1/3) * sin(5*x) - a*x^(1/2)
curve(f(x, a=0.5), 0, 5)
abline(h=0, lty=3)

uniroot(f, c(.1,1), a=0.5) # uniroot: Brent's metods, region: (0.1, 1) #비선형 근 찾기(이분법에 기초)
uniroot(f, c(.1,.5), a=0.5) # region: (0.1, 0.5)

library(rootSolve)
zpts = uniroot.all(f, c(0,5), a=0.5) # find all zeros in the region
zpts
yz = rep(0, length(zpts))
points(zpts, yz)

ff = function(x) sin(x)+1
uniroot(ff, c(-pi,0)) # error when touches but not cross
uniroot(ff, c(-pi, -pi/2))

## finding zeros of a function (based on Newton-Raphson)
library(rootSolve)
fs = function(s) s^3 - 3*s^2 + 4*rho
rho = 0.96
curve(fs(x), 0, 3); abline(h=0)

multiroot(fs, c(1.5, 2.5))

model = function(x) c(F1 = 10*x[1]+3*x[2]^2-3, 
                      F2 = x[1]^2 - exp(x[2])-2)
(ss = multiroot(model, c(1,1))) # simultaneous equations without an explicit Jacobian

derivs = function(x) matrix(c(10, 6*x[2], 2*x[1], -exp(x[2])), 
                            norw=2, byrow=T)
(ssJ = multiroot(model, c(0,0), jacfunc=derivs)) # providing a Jacobian

# nleqslv function using a Broyden (Newton with a linesearch...)

# BBsolve in BB package: finds a local optimum (useful for high-dim)

## Numerical differentiation: numDeriv package
library(numDeriv)

f = function(u) {
  x = u[1]; y = u[2]; z = u[3]
  return(2*x + 3*y^2 - sin(z))
}
# gradient of f at (1,1,0)
grad(f, c(1,1,0)) 

# Jacobian of F at (2,1)
F = function(x) c(x[1]^2 + 2*x[2]^2-3, cos(pi*x[1]/2) - 5*x[2]^3)
jacobian(F, c(2,1))

# Hessian of f at (1,1,0)
hessian(f, c(1,1,0))

zapsmall(hessian(f, c(1,1,0)))

## numerical differentiation using pracma package

## integrate: basic integration function
f = function(x) exp(-x) * cos(x)
(q = integrate(f, 0, pi))
str(q)

v = 0.5*(1+exp(-pi))
abs(q$value-v) # error

# trapz function in pracma 
library(pracma)

xs = seq(0, pi, length.out=101)
ys = f(xs)
trapz(xs, ys)

# Simpson's rule
simp = integral(f, 0, pi, method="Simp")
abs(simp-v) # error

# Gaussian quadrature: library(gaussquad)

## integrals in higher dimension
# numerical integration does not work for dimension higher than 3
# More accurate: SparseGrid package


### Optimization 

## 1-dim optimization 
# optimize function
f = function(x) x * sin(4*x)
optimize(f, c(0,3))
optimize(f, c(1,3), maximum=T)

## multi-dim optimization
# optim function
f = function(x) {
  x1 = x[1]; x2 = x[2]
  return(1/x1 + 1/x2 + (1-x2)/(x2*(1-x1)) +
           1/((1-x1)*(1-x2)))
}

optim(c(.5,.5), f) # default: Nelder-Mead

# nlm function, BB packages...

## constrained optimization with constrOptim function


### optimization with local maxima

## simulated annealing: GenSA package

## genetic algorithm: DEoptim package