#The largest plant genome is 150GBP (150e9 bases)
#The probability of having any number of incorrect basecalls for a genome of this size
#	is given by a binomial distribtion with parameters 150e9 and p,
#	where p is the probability that a base is called incorrectly.
#	To find the probability of calling each base such that there is a 95%
#	probability that there are no errors in the genome, we use the pdf of the
#	binomial to find F(X = 0).

x <- seq(from = -14, to = -10, length.out = 100)
p_0_incorrect = function(x){ dbinom(0,150e9,x) }
y<-sapply(10^x, FUN=p_0_incorrect)
ggplot() + geom_line(aes(x,y)) + ylab('P(0 incorrect calls)') + xlab('Probability of an incorrect call, exponent') + ggtitle("Probability of 0 incorrect calls in a genome of length 150GB")
abline(a = .95, b = 0)

#Approximately 95% is ~ x = 10^-12.466
#Analytically the solution is 1-exp(log(.95)/150e9) (appx 3.419487e-13)
#Human: 3.3e9
#Amoeba:670e9

##
## plot a 3d space of to find appropriate
##
p_incorrect = function(cov, percent, error){ sum(dbinom( floor(percent * cov ):cov, cov, error  )) }

#full plot
x <- 1:50
y<-seq(.01,1,by=.01)
z<- outer(x,y,Vectorize(p_incorrect), error=.01)
p <- list(x=x,y=y,z=z)
image2D(p,xlab="Coverage",ylab="Consensus Level",main="Probability of an Incorrect Call",contour=list(levels=c(.05,.001,1e-7,1e-11),col="white",labcex=.8,lwd=4),lighting=TRUE,rasterImage=TRUE)

#Note that as the amount of coverage increases, the level of consensus required to call a base at a particular level of certainty decreases exponentially

#This also means that at a fixed consensus level, the probability of an incorrect call decreases exponentially with coverage
x <- 2:50
y <- sapply(x, FUN=p_incorrect, percent = .4, error = .01)
ggplot() + geom_line(aes(x,y)) + scale_y_log10() + ylab('P(incorrect) under binomial, log10 scale') + xlab('Coverage') + ggtitle("Probability incorrect call with >40% consensus, 1% error")

#We expect the runtime of tools like GATK to increase exponentially with coverage as well, but this requires data collection



p_incorrect_cutoff = function(cov, cutoff, error){sum(dbinom(cutoff:cov,cov,error))}
minimum_cutoff = function(coverages, p){
	coverage = 1:coverages
	min(which(sapply(coverage, FUN=p_incorrect_cutoff, cov = max(coverage), error = .01) < p))
	} #gives first cutoff where p_incorrect < p
coverage <- 1:50
cutoffs <- sapply(coverage, FUN=minimum_cutoff, p=1e-11) #seems some kind of number pattern here, possibly associated with the binomial
#For p=1e-11, the number of times the next cutoff appears equals the number of times the current cutoff appears + 3 if the current number is even, +2 if the current number is odd
ggplot() + geom_line(aes(coverage,cutoffs)) + ggtitle("Cutoff if P(incorrect call) = 1e-11, P(error) = .01")




####EM MODEL
pn <- function(zn,z,mu){zn/z * (1 - 3 * mu) + (z - zn)/z * mu}
like <- function(xn,x,pn){sum(log(dbinom(xn,x,pn)))}
multilike <- function(xn,pn){prod(pn^xn)}
binlike <- function(xn,pn){dbinom(4,sum(xn),pn^xn}
likemu <- function(xn,x,zn,z,mu){sum(log(dbinom(xn,x,pn(zn,z,mu))))}
#A, T, G, C
xn <- c(3,0,1,1)
x <- c(5,5,5,5)
zn <- c(2,0,1,0)
z <- c(3,3,3,3)
zz <- z - zn

mu <- (xn / (x + 1) - zn / z) * (3*z)/(zz - 3*zn) #mu to maximize probability given N
mu<-seq(0,1,by=.01) #mu for graphing
vals <- sapply(mu,FUN = function(y) likemu(xn,x,zn,z,y))
p <- sapply(mu, FUN=function(y) pn(zn,z,y))
phat <- sapply(mu, FUN=function(y) (prod(pn(zn,z,y)^(xn)))^(1/sum(xn)))
which(vals==max(vals)) #find index of max mu

test <- data.frame(xn, x, zn, z, zz)
apply(X = test,MARGIN = 1,FUN = function(x) likemu(x['xn'],x['x'],x['zn'],x['z']))
sapply(mu,FUN = function(y) likemu(xn[1],x[1],zn[1],z[1],y))



ggplot() + geom_line(aes(mu, vals))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

sapply(X=seq(1,length(mu)), FUN = function(y) gm_mean(c(pa[y],pa[y],pa[y],pg[y])))





