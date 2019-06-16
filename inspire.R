xs <- xcmsSet(file,method="centWave",ppm=20,peakwidth=c(10,50))
xset2 <- group(xset2, bw = 10)
xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden")
an <- xsAnnotate(xs1, polarity="negative", nSlaves=2)