N <- 10

x <- rnorm(2*N)
y <- rnorm(2*N)
x[(N+1):(2*N)] <- x[(N+1):(2*N)] + 0.1
xy <- complex(real=x,imaginary=y)
group <- rep(1:2,each=N)
participant <- rep(1:N,times=2)

data <- data.frame(xy,group)
group=NULL
participant=NULL
results <- analysecplx(data)

data <- data.frame(xy,group,participant)
analysecplx(data)


x <- rnorm(2*N)
y <- rnorm(2*N)
x <- x + 0.5*y
xy <- complex(real=x,imaginary=y)
group <- rep(1:2,each=N)
participant <- rep(1:N,times=2)

input <- data.frame(xy,group)
analysecplx(input)

input <- data.frame(xy,group,participant)
analysecplx(input)


N <- 20
x <- rnorm(3*N)
y <- rnorm(3*N)
x[(N+1):(2*N)] <- x[(N+1):(2*N)] + 0.1
xy <- complex(real=x,imaginary=y)
group <- rep(1:3,each=N)
participant <- rep(1:N,times=3)

data <- data.frame(xy,group)
analysecplx(input)

input <- data.frame(xy,group,participant)
analysecplx(input)



N <- 20
x <- rnorm(3*N)
y <- rnorm(3*N) + 0.5*x
x[(N+1):(2*N)] <- x[(N+1):(2*N)] + 0.5
xy <- complex(real=x,imaginary=y)
group <- rep(1:3,each=N)
participant <- rep(1:N,times=3)

input <- data.frame(xy,group)
analysecplx(input)

input <- data.frame(xy,group,participant)
analysecplx(input)





N <- 10

x <- rnorm(N)
y <- rnorm(N)
xy1 <- complex(real=x,imaginary=y)
x <- rnorm(N)
y <- rnorm(N)
xy2 <- complex(real=x,imaginary=y)

tsqh.test(xy1,xy2)
h <- Hotelling::hotelling.test(x=data.frame(Re(xy1),Im(xy1)),y=data.frame(Re(xy2),Im(xy2)))
