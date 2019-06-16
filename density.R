d0 <- DOC1[, 'class'] == 'control'
d1 <- DOC1[, 'class'] == 'disease'
base.rate.d1 <- sum(d1) / (sum(d1) + sum(d0))

D1G.density    <- density (DOC1$D1G)
D1G.d0.density <- density (DOC1$D1G[d0])
D1G.d1.density <- density (DOC1$D1G[d1])

D1G.d0.f <- approxfun(D1G.d0.density$x, D1G.d0.density$y)
D1G.d1.f <- approxfun(D1G.d1.density$x, D1G.d1.density$y)

p.d.given.D1G <- function(D1G, base.rate.d1)
{
  p1 <- D1G.d1.f(DOC1$D1G) * base.rate.d1
  p0 <- D1G.d0.f(DOC1$D1G) * (1 - base.rate.d1)
  p1 / (p0 + p1)
}

x <- 1:4000000
y <- p.d.given.D1G (x, base.rate.d1)
plot(x, y, type='l', col='red', xlab='D1G', ylab='estimated p(disease|D1G)')


plot(density(DOC1$D1G), col='blue', xlab='D1G', ylab='estimate p(D1G), p(D1G|disease), p(D1G|control)', main=NA)lines(density(DOC1$D1G), col='red')

fy.x <- npcdens(class~D1G, nmulti=2, data=DOC1)

Pima.eval <- data.frame(type=factor("control"),
                        D1G=seq(min(DOC1$D1G), max(DOC1$D1G), length=4000000))

plot(x, y, type='l', lty=2, col='red', xlab='glu',
     ylab='estimated p(diabetes|glu)')
lines(Pima.eval$D1G, predict(fy.x, newdata=Pima.eval), col="blue")
legend(0, 1, c("Unconditional bandwidth", "Conditional bandwidth"),
       col=c("red", "blue"), lty=c(2, 1))


cdens <- cdplot(as.factor(class)~D1G, data= DOC1, plot = FALSE)
plot(I(as.numeric((as.factor(class)) - 1) ~ jitter(D1G, factor = 2))
     spineplot(as.factor(class)~D1G, data= DOC1, breaks = 6)
     spineplot(fail ~ temperature, breaks = quantile(temperature)))


