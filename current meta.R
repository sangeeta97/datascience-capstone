View(Book2)
> g<- Book2$Label
> y1<- Book2$`619.4746`
> sm.density.compare(y1, cyl.f, xlab="peak intensity value")
Error: could not find function "sm.density.compare"
> library(sm)
Package 'sm', version 2.2-5.4: type help(sm) for summary information
> sm.density.compare(y1, cyl.f, xlab="peak intensity value")
Error in factor(group) : object 'cyl.f' not found
> cyl.f <- factor(g, levels= c("control","disease"),
                  +                                       labels = c("DM", "DN"))
> sm.density.compare(y, cyl.f, xlab="peak intensity value")
Error in is.vector(x) : object 'y' not found
> sm.density.compare(y1, cyl.f, xlab="peak intensity value")
> sm.density.compare(y1, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> View(Book2)
> y2<- Book2$1093`1093.5107`
Error: unexpected numeric constant in "y2<- Book2$1093"
> y2<- Book2$`1093.5107`
> sm.density.compare(y2, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.28 
> colfill<-c(2:(2+length(levels(cyl.f))))
> legend(locator(1), levels(cyl.f), fill=colfill)
> title(main="PIP2(20:1(11Z)/18:2(9Z,12Z))distribution")
> y3<- Book2$`1095.5129`
> sm.density.compare(y3, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.3 
> colfill<-c(2:(2+length(levels(cyl.f))))
> legend(locator(1), levels(cyl.f), fill=colfill)
> title(main="PIP2(22:2(13Z,16Z)/16:0) distribution")
> y1<- Book2$`567.4413`
> cyl.f <- factor(g, levels= c("control","disease"),
                  +                 labels = c("CN", "DM"))
> sm.density.compare(y1, cyl.f, xlab="peak intensity value")
> sm.density.compare(y1, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="DG(14:0/0:0/20:5n3) distribution")
> colfill<-c(2:(2+length(levels(cyl.f))))
> legend(locator(1), levels(cyl.f), fill=colfill)
> y2<- Book2$`574.4997`
> sm.density.compare(y2, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="Cer(d18:0/16:0) distribution")
> colfill<-c(2:(2+length(levels(cyl.f))))
> legend(locator(1), levels(cyl.f), fill=colfill)
> y3<- Book2$`583.436`
> sm.density.compare(y3, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="DG(14:1n5/0:0/20:5n3) distribution")
> colfill<-c(2:(2+length(levels(cyl.f))))
> legend(locator(1), levels(cyl.f), fill=colfill)
> y4<- Book2$`609.5596`
> sm.density.compare(y4, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.02 
> colfill<-c(2:(2+length(levels(cyl.f))))
> legend(locator(1), levels(cyl.f), fill=colfill)
> y5<- Book2$`631.5674`
> sm.density.compare(y5, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="DG(14:0/24:1/(15Z)/0:0) distribution")
> colfill<-c(2:(2+length(levels(cyl.f))))
> legend(locator(1), levels(cyl.f), fill=colfill)
> y6<- Book2$`657.5952`
> sm.density.compare(y6, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="CE(20:2/6Z,9Z)) distribution")
> colfill<-c(2:(2+length(levels(cyl.f))))
> legend(locator(1), levels(cyl.f), fill=colfill)
> Y7<- Book2$`658.5942`
> sm.density.compare(y6, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> sm.density.compare(y7, cyl.f, xlab="peak intensity value", model = "equal")
Error in is.vector(x) : object 'y7' not found
> sm.density.compare(Y7, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="Cer(d18:0/22:0) distribution")
> colfill<-c(2:(2+length(levels(cyl.f))))
> legend(locator(1), levels(cyl.f), fill=colfill)
> legend(locator(1), levels(cyl.f), fill=colfill)
> y8<- Book2$`671.4488`
> sm.density.compare(y8, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.03 
> title(main="CL(HMDB57478) distribution")
> colfill<-c(2:(2+length(levels(cyl.f))))
> legend(locator(1), levels(cyl.f), fill=colfill)
> Book2$`673.5268`
[1]  58173.22  53487.54  44851.91  38553.00  37387.91  30108.54  35995.31
[8]  30718.20  91642.96  47011.27  66193.10  38418.74  56094.33  41117.88
[15] 198537.63 201451.34 218467.27 253962.07  50313.22  93269.01
> y9<- Book2$`673.5268`
> sm.density.compare(y9, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.02 
> title(main="SM(d18:1/14:0) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y10<- Book2$`675.5383`
> sm.density.compare(y10, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="DG(20:1n9/0:0/22:6n3) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y11<- Book2$`704.5829`
> sm.density.compare(y11, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.01 
> title(main="SM(d18:0/16:0) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y12<- Book2$'766.5411'
> sm.density.compare(y12, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.17 
> title(main="PE(18:1/20:3) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y13<- Book2$'791.6904'
> sm.density.compare(y13, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.65 
> title(main="TG(14:1/18:4/o-18:0) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y14<- Book2$'820.4921'
> sm.density.compare(y14, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.3 
> title(main="PS(18:1/18:2) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y15<- Book2$`843.7782`
> sm.density.compare(y15, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.57 
> title(main="TG(18:1/14:0/20:0) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y16<- Book2$'957.5902'
> sm.density.compare(y16, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.01 
> title(main="LysoPE(0:0/18:1) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y17<- Book2$'997.5302'
> sm.density.compare(y17, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="LysoPE(0:0/20:5) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y18<- Book2$`1005.592`
> sm.density.compare(y18, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="LysoPE(0:0/20:3) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y19<- Book2$.`1029.5925`
Error: unexpected symbol in "y19<- Book2$.`1029.5925`"
> y19<- Book2$`1029.5925`
> sm.density.compare(y19, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="LysoPC(18:4) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y20<- Book2$`1133.6981`
> sm.density.compare(y20, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.04 
> title(main="GM3(d18:1/16:0) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y21<- Book2$`1135.7124`
> sm.density.compare(y21, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.01 
> title(main="GM3(d18:0/16:0) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y22<- Book2$`1159.7011`
> sm.density.compare(y22, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="GM3(d18:1/18:1) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y23<- Book2$'657.5807'
> sm.density.compare(y23, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.1 
> title(main="DG(18:0/0:0/22:2n6) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y24<- Book2$`684.61`
> sm.density.compare(y24, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.02 
> title(main="Cer(d18:0/24:1) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y25<- Book2$'685.6138'
> sm.density.compare(y25, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.03 
> title(main="DG(20:0/0:0/22:2n6) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y26<- Book2$'685.6248'
Warning message:
  Unknown column '685.6248' 
> y26<- Book2$`686.6248`
> sm.density.compare(y26, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.19 
> title(main="Cer(d18:0/24:0) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y27<- Book2$'700.5294'
> sm.density.compare(y27, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.26 
> title(main="PE(16:0/P-18:1) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y28<- Book2$'712.6413'
> sm.density.compare(y28, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.06 
> title(main="Cer(d18:0/26:0) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y29<- Book2$'714.6563'
> sm.density.compare(y29, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.08 
> y30<- Book2$'716.5259'
> sm.density.compare(y30, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.27 
> title(main="PE(14:0/20:1) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y31<- Book2$'719.6187'
> sm.density.compare(y31, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.16 
> title(main="TG(14:0/14:0/14:1) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y32<- Book2$'722.516'
> sm.density.compare(y32, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.03 
> 
  > title(main="PE(P-18:1/18:3) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y33<- Book2$'729.64'
> sm.density.compare(y33, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.16 
> title(main="TG(14:0/16:0/14:1) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y34<- Book2$'731.6172'
> sm.density.compare(y34, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="TG(15:0/14:1/14:1) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y35<- Book2$'731.6535'
> sm.density.compare(y35, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.29 
> title(main="DG(22:0/0:0/22:2n6) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y36<- Book2$'745.6349'
> sm.density.compare(y36, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.33 
> title(main="TG(14:0/14:1/16:1) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y37<- Book2$'747.6503'
> sm.density.compare(y37, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.22 
> title(main="TG(16:0/14:0/14:1) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y38<- Book2$'750.5473'
> sm.density.compare(y38, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.36 
> title(main="PE(P-18:1/20:3) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y39<- Book2$'759.6855'
> sm.density.compare(y39, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.14 
> title(main="TG(14:1/14:1/o-18:1) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y40<- Book2$'773.6658'
> sm.density.compare(y40, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.42 
> title(main="TG(14:1/14:0/18:1) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y41<- Book2$'774.5776'
> sm.density.compare(y41, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.63 
> title(main="PC(P-18:0/20:4) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y42<- Book2$'775.6815'
> sm.density.compare(y42, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.1 
> title(main="TG(14:0/18:1/14:0) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y43<- Book2$'787.7168'
> sm.density.compare(y43, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.14 
> title(main="TG(14:1/16:1/o-18:0) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y44<- Book2$'792.4604'
> sm.density.compare(y44, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.04 
> title(main="PS(14:0/20:3) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y45<- Book2$'802.6093'
> sm.density.compare(y45, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.54 
> title(main="PC(22:4/P-18:0) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y46<- Book2$'822.6351'
> sm.density.compare(y46, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.17 
> title(main="PC(O-20:0/20:4) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y47<- Book2$'865.6723'
> sm.density.compare(y47, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.31 
> 
  > title(main="TG(15:0/18:4/22:6) distribution")
> 
  > legend(locator(1), levels(cyl.f), fill=colfill)
> y48<- Book2$'909.7045'
> sm.density.compare(y48, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.57 
> title(main="TG(18:4/20:2n6/20:5) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y49<- Book2$'945.7976'
> sm.density.compare(y49, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="TG(22:2/20:2n6/18:3) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y50<- Book2$'742.5443'
> sm.density.compare(y50, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="PC(15:0/18:2) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y51<- Book2$'773.5365'
> sm.density.compare(y51, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="PG(18:1/18:1) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y52<- Book2$'774.5425'
> sm.density.compare(y52, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="PC(15:0/22:5) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> y53<- Book2$'784.5901'
> sm.density.compare(y53, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="PC(14:0/22:5) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> sm.density.compare(y53, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="PC(14:0/22:2) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> 
  > y55<- Book2$'1189.761'
> sm.density.compare(y55, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.04 
> title(main="GM3(d18:1/20:0) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> sm.density.compare(y53, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0 
> title(main="PC(14:0/22:2) distribution")
> legend(locator(1), levels(cyl.f), fill=colfill)
> sm.density.compare(y23, cyl.f, xlab="peak intensity value", model = "equal")

Test of equal densities:  p-value =  0.11 
> 