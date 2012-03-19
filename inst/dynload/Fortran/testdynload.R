require(odesolve)
dyn.load("mymod.so")
parms <- c(k1 = 0.04, k2 = 1e4, k3=3e7)
my.atol <- c(1e-6, 1e-10, 1e-6)
times <- c(0, 4*10^(-1:10))

print(system.time(
                  out1 <- lsoda(c(1.0,0.0,0.0),times,"myderivs",
                                parms,
                                rtol=1e-4,atol=my.atol,jacfunc="myjac",
                                dllname="mymod")
                  )
      )
print(out1)

lsexamp <- function(t, y, p)
  {
    yd1 <- -p["k1"] * y[1] + p["k2"] * y[2]*y[3]
    yd3 <- p["k3"] * y[2]^2
    list(c(yd1,-yd1-yd3,yd3),c(massbalance=sum(y)))
  }
exampjac <- function(t, y, p)
  {
    c(-p["k1"],	 p["k1"],  0,

        p["k2"]*y[3],
      - p["k2"]*y[3] - 2*p["k3"]*y[2],
                       2*p["k3"]*y[2],

      p["k2"]*y[2],  -p["k2"]*y[2],  0
      )
  }
print(system.time(
out2 <- lsoda(c(1,0,0),times,lsexamp, parms, rtol=1e-4, atol= my.atol,
              jac = exampjac)
))

print(all.equal(as.vector(out1),as.vector(out2[,1:4])))
