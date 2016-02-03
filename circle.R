# Not everything that looks non-linear is non-linear 
radius <- 10

theta<-c(0:365)

y <- radius * sin(theta) + rnorm(366,mean=0,sd=1)
x <- radius * cos(theta)

plot(x,y)

#OLS
#bestfit <- lm(y~sin(theta))
#yest <- bestfit$coefficients[1]+bestfit$coefficients[2]*sin(theta)
# beta0 <-bestfit$coefficients[1]
# beta <- bestfit$coefficients[2]

# Orthogonal regression
v <- prcomp(cbind(x,y))$rotation
beta <- v[2,1]/v[1,1]
beta0 <- mean(y) - beta*mean(x)

yest <- beta0 + beta * sin(theta)



par(new=T)
plot(x,yest,yaxt="n",ylab="",pch="-")
text(0,2," y=0.4712516 - 12.40906 X")
text(0,0,"Die gerade linie is gottlos und unmoralisch.")
text(0,-2,"-Friedensreich Hundertwasser, 1958")

#summary(bestfit)