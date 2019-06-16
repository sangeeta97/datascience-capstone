g<- Book2$Label
cyl.f <- factor(g, levels= c("control","disease"),
                                      labels = c("DM", "DN")) 
sm.density.compare(y, cyl.f, xlab="peak intensity value")
title(main="Cer(d18:0/16:0) distribution")
colfill<-c(2:(2+length(levels(cyl.f)))) 
legend(locator(1), levels(cyl.f), fill=colfill)
y1<- Book2$`619.4746`
sm.density.compare(y1, cyl.f, xlab="peak intensity value")
title(main="DG(16:1n7/0:0/22:6n3) distribution")
colfill<-c(2:(2+length(levels(cyl.f))))
legend(locator(1), levels(cyl.f), fill=colfill)
sm.density.compare(y1, cyl.f, xlab="peak intensity value", model = "equal")
cyl.f <- factor(g, levels= c("control","disease"),
                +                 labels = c("CN", "DM"))
