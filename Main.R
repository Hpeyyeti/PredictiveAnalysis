heartval <- read.table("heart_val.txt", header = TRUE)
library(leaps)
names(heartval) <- c( "cost", "gender","age",  "interv", "drugs", "ervisit", "complic", "comorb", "dur")

head(heartval)
summary(heart$cost)
summary(heart$age)
summary(heart$dur)
plot(heart)

sub1.heart <- subset(heart, cost > 0 & cost< 10000)
head(sub1.heart)
dim(sub1.heart)

fit1 <- lm(cost ~ age + gender + interv + drugs + ervisit + complic + comorb + dur, data=sub1.heart)

summary(fit1)
qqnorm(rstudent(fit1))

abline(0,1)
plot(fitted(fit1), rstudent(fit1),xlab = "fitted values", ylab = "Studentized residuals", main = "Residual plot")
abline(0,0)

heart_box <- read.table("heart_box.txt", header = TRUE)
head(heart_box)
names(heart_box) <- c("id", "cost", "gender", "age", "interv", "drugs", "ervisit", "complic", "comorb", "dur","claims")

boxplot(interv~claims, data=heart_box ,main="Number of interventions", font.main=3, cex.main=1.2, xlab="Total cost of claims", ylab="interventions", font.lab=3)

boxplot(drugs~claims, data=heart_box ,main="Number of Drugs", font.main=3, cex.main=1.2, xlab="Total cost of claims", ylab="# of drugs", font.lab=3)

boxplot(ervisit~claims, data=heart_box ,main="Number of ERvists", font.main=3, cex.main=1.2, xlab="Total cost of claims", ylab="# of ERvistis", font.lab=3)

boxplot(dur~claims, data=heart_box ,main="Duration of Treatment", font.main=3, cex.main=1.2, xlab="Total cost of claims", ylab="# of days", font.lab=3)

boxplot(complic~claims, data=heart_box ,main="Number of complications", font.main=3, cex.main=1.2, xlab="Total cost of claims", ylab="# of complications", font.lab=3)

boxplot(comorb~claims, data=heart_box ,main="Number of comorbidities", font.main=3, cex.main=1.2, xlab="Total cost of claims", ylab="# of comorbidities", font.lab=3)
axis(side=2,at=seq(0, 3, by=1))
###################################################Check transformation required###########################################
library(MASS)
fit_box= lm(cost ~ age + gender + interv + drugs + ervisit + complic + comorb + dur, data = sub1.heart)
boxcox(fit_box)
boxcox(fit_box, plot=F)
attributes(BC)
BC <- boxcox(fit_box)
BC$x[BC$y==max(BC$y)]
# confidence interval
S <- max(BC$y) - 0.5*qchisq(0.95,1)
S
BC$x[BC$y>S]


########################################Exhaustive search for best subset model###############################################
sub <- regsubsets(x=cbind(sub1.heart$age,sub1.heart$drugs,sub1.heart$interv,sub1.heart$ervisit,sub1.heart$complic,sub1.heart$comorb,sub1.heart$dur), y=(sub1.heart$cost^0.1),  method = "exhaustive", all.best = FALSE, nbest = 3, data=sub1.heart)
summary(sub)
dim(sub1.heart)
Cp <- summary(sub)$cp
AdjR2 <- summary(sub)$adjr2
SSRes <- summary(sub)$rss
R2 <- summary(sub)$rsq
Matrix <- summary(sub)$which

# We want this for our table
p <- apply(Matrix,1, sum)
MSE <- SSRes/(731-p) # n = 785, change if necessary

# Create a table
output <- cbind(p, Matrix, SSRes, R2, AdjR2, MSE, Cp)
colnames(output)[3:9] <- c("age","drugs", "interv","ervisit", "complic", "comorb", "dur") 
output

cor(sub1.heart[3:10])

#########################Modedl selection using cost^0.10101##########################

# Forward variable selection
#use alpha-to-enter=0.1
fit.0 <- lm(cost^0.10101~1, data = sub1.heart)
add1(fit.0, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, test = "F")
fit.1 <- lm(cost^0.10101~interv, data = sub1.heart)
add1(fit.1, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, test = "F")
fit.2 <- lm(cost^0.10101~comorb+interv, data = sub1.heart)
add1(fit.2, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, test = "F")
fit.3 <- lm(cost^0.10101~comorb+interv+dur, data = sub1.heart)
add1(fit.3, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, test = "F")
fit.4 <- lm(cost^0.10101~comorb+interv+dur+complic, data = sub1.heart)
add1(fit.4, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, test = "F")
fit.5 <- lm(cost^0.10101~ervisit+comorb+interv+dur+complic, data = sub1.heart)
add1(fit.5, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, test = "F")
fit.6 <- lm(cost^0.10101~ervisit+comorb+interv+dur+complic+age, data = sub1.heart)
add1(fit.6, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, test = "F")
summary(fit.6)
####conclusion:except gender drugs Rsq=0.6271
#######################################################
# Backward variable selection 
#use alpha-to-leave=0.1
fit.4 <- lm(cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, data=sub1.heart)
drop1(fit.4, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, test = "F")	
fit.3 <- lm(cost^0.10101~age + interv + drugs + ervisit + complic + comorb + dur, data = sub1.heart)
drop1(fit.3, cost^0.10101~age + interv + drugs + ervisit + complic + comorb + dur, test = "F")
fit.2 <- lm(cost^0.10101~age + interv + ervisit + comorb + complic + dur, data = sub1.heart)
drop1(fit.2, cost^0.10101~age + interv + ervisit + comorb + complic + dur, test = "F")
summary(fit.2)
####conclusion:except gender drugs Rsq=0.6271
#######################################################

########### Stepwise variable selection 
#alpha-to-leave=alpha-to-enter=0.1
fit.0 <- lm(cost^0.10101~1, data = sub1.heart)
add1(fit.0, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, data=sub1.heart, test = "F")
fit.1 <- lm(cost^0.10101~interv, data = sub1.heart)
drop1(fit.1, cost^0.10101~interv, test = "F")
add1(fit.1, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, data=sub1.heart, test = "F")
fit.2 <- lm(cost^0.10101~comorb+interv, data = sub1.heart)
drop1(fit.2, cost^0.10101~comorb+interv, test = "F")
add1(fit.2, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, data=sub1.heart, test = "F")
fit.3 <- lm(cost^0.10101~comorb+interv+dur, data = sub1.heart)
drop1(fit.3, cost^0.10101~comorb+interv+dur, test = "F")
add1(fit.3, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, data=sub1.heart, test = "F")
fit.4 <- lm(cost^0.10101~comorb+interv+dur+complic, data = sub1.heart)
drop1(fit.4, cost^0.10101~comorb+interv+dur+complic, test = "F")
add1(fit.4, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, data=sub1.heart, test = "F")
fit.5 <- lm(cost^0.10101~ervisit+comorb+interv+dur+complic, data = sub1.heart)
drop1(fit.5, cost^0.10101~ervisit+comorb+interv+dur+complic, test = "F")
add1(fit.5, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, data=sub1.heart, test = "F")
fit.6 <- lm(cost^0.10101~ervisit+comorb+interv+dur+complic+age, data = sub1.heart)
drop1(fit.6, cost^0.10101~ervisit+comorb+interv+dur+complic+age, test = "F")
add1(fit.6, cost^0.10101~age + gender + interv + drugs + ervisit + complic + comorb + dur, data=sub1.heart, test = "F")
summary(fit.6)
####conclusion:except gender drugs Rsq=0.6271

#################################Check model assumptions#####################
qqnorm(rstandard(fit.6))
abline(0,1)

plot(fitted(fit.6), rstudent(fit.6),xlab = "fitted values", ylab = "Studentized residuals", main = "Residual plot")
abline(0,0)

plot(fitted(fit.6), rstandard(fit.6))
abline(0,0)


###########################Outlier Analysis#################################
#standardized residuals
mse=summary(fit.6)$sigma^2
d=resid(fit.6)/sqrt(mse)
rabs <- abs(d)
a <- cbind(sub1.heart, rabs)
asorted <- a[order(-rabs), ]
asorted[1:5, ]

rstud <- pres(fit.6)
b <- cbind(sub1.heart, rstud)
bsorted <- b[order(-rstud), ]
bsorted[1:10, ]

Press.res=resid(fit.6)/(1-h.ii)
Press.res[5]
####check leverage
X <- as.matrix(cbind(1,sub1.heart$age,sub1.heart$interv,sub1.heart$ervisit,sub1.heart$complic,sub1.heart$comorb,sub1.heart$dur))
h.ii = hat(X)
heart.hii <- cbind(a,h.ii)
sort.hii <- heart.hii[order(-h.ii), ]
sort.hii[1:10,]
sort.rabs <- heart.hii[order(-rabs), ]
sort.rabs[1:5,]
h.ii[5]
p=7
n=785
2*p/n

cook <- cooks.distance(fit.6)
sort.cook <- cook[order(-cook)]
sort.cook
qf(0.25,7,778)
#cook's distance(highest is .46) < 0.607 i.e remoing obs i doesnt moves coeff outside 70% ci 
##45 79 183

1-3*p*1/n
1+3*p*1/n
2/sqrt(n)
2*sqrt(p/n)
cooks.distance(fit.6)[5]
Press.res=resid(fit.6)/(1-h.ii)
Press.res[45]
dfbetas(fit.6)[5,]
dffits(fit.6)[5]
covratio(fit.6)[5]
fit45 <- update(fit.6, subset = -5)
summary(fit45)
anova(fit45)
anova(fit.6)
qt(0.025,778)

fit.7 <- lm(cost^0.10101~ervisit+comorb+interv+dur+complic, data = sub1.heart)
summary(fit.7)
anova(fit.7)

qqnorm(rstandard(fit.7))
abline(0,1)

plot(fitted(fit.7), rstudent(fit.7),xlab = "fitted values", ylab = "Studentized residuals", main = "Residual plot")
abline(0,0)

fit.7 <- lm(cost^0.1~comorb+interv+dur, data = sub1.heart)
summary(fit.7)
anova(fit.7)
PRESS(fit.7)
1- (31.39155/65.1157)

########################Duration model###################


heart <- read.table("heart_costs.txt", header = FALSE)
library(leaps)
names(heart) <- c("id", "cost", "age", "gender", "interv", "drugs", "ervisit", "complic", "comorb", "dur")

head(heart)

plot(heart)

sub1.heart <- subset(heart, cost > 0)

dim(sub1.heart)
cost2=log(heart$cost+1)
model2 <- lm(dur ~ age + gender + interv +drugs + ervisit + complic + comorb + cost, data=sub1.heart)
plot(model2)
fit1 <- lm(I(log(sub1.heart$cost)) ~ sub1.heart$age + sub1.heart$gender + sub1.heart$interv + sub1.heart$drugs + sub1.heart$ervisit + sub1.heart$complic + sub1.heart$comorb + sub1.heart$dur, data=sub1.heart)

summary(model2)
qqnorm(rstudent(fit2))
qqline(rstudent(fit2))
abline(0,1)
plot(fitted(fit2), rstudent(fit2))
abline(0,0)

summary(fit2)


sample.val <- sample(c(1:731),300)
sample.est <- (1:731)[-sample(c(1:731),300)]
heart.est <- heartval[sample.est,]
heart.val <- heartval[sample.val,]
head(pinot.val)
fit.c <- lm(I(cost^.1)~ interv + dur + comorb, data=heartval)
coefficients(fit.c)
y_hat <- predict(fit.c, heart.val[,c(7:9)])
y.new=heart.val[,1]
pred_error <- y.new - y_hat

rsq_pred = 1- sum(pred_error^2)/sum((y.new - mean(y.new))^2)

# Repeat very often:
beta0 <- numeric() 
beta1 <- numeric()
beta2 <- numeric()
rsq_pred <- numeric()
MSP <- numeric()
for (i in 1:10){
	sample.val <- sample(c(1:731),220)
	sample.est <- (1:731)[-sample(c(1:731),220)]
	heart.est <- heartval[sample.est,]
	heart.val <- heartval[sample.val,]

fit.c <- lm(I(cost^.1)~ interv + dur + comorb, data=heartval)
coefficients(fit.c)
y_hat <- predict(fit.c, heart.val[,c(7:9)])
y.new=heart.val[,1]
pred_error <- y.new - y_hat
	t=(y.new - mean(y.new))
		
	beta0[i] <- coef(fit)[1]
	beta1[i] <- coef(fit)[2]
	beta2[i] <- coef(fit)[3]
	MSP[i] <- sum(pred_error^2)/220
	rsq_pred[i] = 1- sum(pred_error^2)/sum(t^2)
}

hist(rsq_pred)
summary(rsq_pred)


X <- cbind( pinot$Flavor, pinot$Oakiness)
solve(t(X)%*%X)

solve(cor(pinot[4:5]))
