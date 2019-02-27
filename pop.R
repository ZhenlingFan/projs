rm(list = ls())
library(spdep)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(rgdal)
library(MASS)
library(sp)
library(ncf)
#--------Exploratory Data Analysis--------
popdf <- read.table("D:/UWM3/STAT 679/Project/pop.csv", header = TRUE, sep=",")
nrow(popdf)
str(popdf)
head(popdf)
summary(popdf)
#sample correlation
cor(popdf$pop, popdf$NTL)
cor(popdf$pop, popdf$GDP)
cor(popdf$pop, popdf$emp)
cor(popdf$pop, popdf$invest)
cor(popdf$pop, popdf$income)
cor(popdf$pop, popdf$retail)
cor(popdf$pop, popdf$edu)
cor(popdf$pop, popdf$health)
#Visualization
p1 <- ggplot(popdf, aes(x = NTL, y = pop)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth() + 
  #ggtitle("NTL") + 
  theme_bw()
p2 <- ggplot(popdf, aes(x = GDP, y = pop)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth() + 
  #ggtitle("GDP") + 
  theme_bw()
p3 <- ggplot(popdf, aes(x = emp, y = pop)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth() + 
  #ggtitle("emp") + 
  theme_bw()
p4 <- ggplot(popdf, aes(x = invest, y = pop)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth() + 
  #ggtitle("invest") + 
  theme_bw()
p5 <- ggplot(popdf, aes(x = income, y = pop)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth() + 
  #ggtitle("income") + 
  theme_bw()
p6 <- ggplot(popdf, aes(x = retail, y = pop)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth() + 
  #ggtitle("retail") + 
  theme_bw()
p7 <- ggplot(popdf, aes(x = edu, y = pop)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth() + 
  #ggtitle("edu") + 
  theme_bw()
p8 <- ggplot(popdf, aes(x = health, y = pop)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth() + 
  #ggtitle("health") + 
  theme_bw()
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol=4,
             bottom = "Relationships of Socio-economic variables and population")
#--------Spatial Correlation--------
#read in shape file
pop_shp <- readOGR(dsn="D:/UWM3/STAT 679/Project", layer="GD_pop", verbose = FALSE)
#class(pop_shp)
names(pop_shp)
pop_nb = poly2nb(pop_shp, row.names = pop_shp$cvalue)
#W - default row-standardized weights
listw_popW = nb2listw(pop_nb, style = "W", zero.policy = TRUE)
#B - binary weights
listw_popB = nb2listw(pop_nb, style = "B", zero.policy = TRUE)
#Visualize the weight
coords = coordinates(pop_shp)
plot(listw_popB, coords, col="blue", cex=0.1)
#--------Exploratory Spatial Data Analysis--------
at = seq(from=0, to=100, length.out=5)
p0 = spplot(pop_shp, "pop", cuts=4, main="Pop")
p1 = spplot(pop_shp, "NTL", cuts=4, main="NTL")
p2 = spplot(pop_shp, "GDP", cuts=4, main="GDP")
p3 = spplot(pop_shp, "emp", cuts=4, main="emp")
p4 = spplot(pop_shp, "invest", cuts=4, main="invest")
p5 = spplot(pop_shp, "income", cuts=4, main="income")
p6 = spplot(pop_shp, "retail", cuts=4, main="retail")
p7 = spplot(pop_shp, "edu", cuts=4, main="edu")
p8 = spplot(pop_shp, "health", cuts=4, main="health")
grid.arrange(p0,p1,p2,p3,p4,p5,p6,p7,p8, ncol=3,
             bottom = "Visualization of variables ")

#ptest = spplot(pop_shp, "GDP", cuts=4,main="GDP")
#ptest
#http://www.nickeubank.com/wp-content/uploads/2015/10/RGIS3_MakingMaps_part1_mappingVectorData.html
#--------Moran's I and Geary's C--------
#Moran's I test based on randomization
set.seed(1)
moran.test(pop_shp$pop, listw_popB, zero.policy = TRUE, alternative="two.sided")
#Moran's I test based on Monte Carlo
set.seed(2)
moran.mc(pop_shp$pop, listw_popB, zero.policy = TRUE, nsim=999)
#Geary's C test
set.seed(3)
geary.test(pop_shp$pop, listw_popB, zero.policy = TRUE, alternative = "two.sided")
set.seed(4)
geary.mc(pop_shp$pop, listw_popB, zero.policy = TRUE, nsim=999)

#--------Ordinary Least Squares--------
data = pop_shp@data
m1 <- lm(pop~NTL + GDP + emp + invest + income + retail + edu + health,
         data = popdf)
summary(m1)
#Confidence intervals
cbind(coefest = coef(m1), confint(m1))

#Variable selection
#Backward elimination based on AIC
m2 = step(m1, direction = "back")
summary(m2)
#Backward elimination based on BIC
n = nrow(popdf) #n is the sample size
m3 = step(m1, k=log(n))
summary(m3)

stdresids <- stdres((m3))
plot(x=fitted(m3),y=stdresids,main="Standardized Residuals vs. Fitted",
     ylim=c(-3.5, 3.5),
     xlab="Fitted values",
     ylab="Standardized residuals")
abline(h = c(-3,0,3), col=c("red", "black", "red"), lty=2)

#Model Diagnostics
par(mfrow=c(1,2))
plot(m3,which=c(1,2), cex=0.1)

m3.log = lm(log(pop) ~ GDP + emp + invest + edu + health,
            data=popdf)
par(mfrow=c(1,2))
plot(m3.log, which=c(1,2), cex=0.1)
#Moran's I test
moran.test(resid(m3), listw_popW, randomisation = FALSE, zero.policy = TRUE)

#------Spatial Linear Models--------
#SAR
m3_sar = spautolm(pop ~ GDP + emp + invest + edu + health,
                  data = popdf,
                  listw = listw_popW, zero.policy = TRUE, family = "SAR")
summary(m3_sar)
set.seed(5)
moran.mc(resid(m3_sar), listw_popW, zero.policy = TRUE, nsim=999)
#CAR
set.seed(6)
m3_car = spautolm(pop ~ GDP + emp + invest + edu + health,
                  data = popdf,
                  listw = listw_popW, zero.policy = TRUE, family = "CAR")
summary(m3_car)
set.seed(6)
moran.mc(resid(m3_car), listw_popW, zero.policy = TRUE, nsim=999)

#
m2_sar = spautolm(pop ~ GDP + emp + invest + income + edu + health,
                  data = popdf,
                  listw = listw_popW, zero.policy = TRUE, family = "SAR")
summary(m2_sar)
set.seed(7)
moran.mc(resid(m2_sar), listw_popW, zero.policy = TRUE, nsim=999)
#
m2_car = spautolm(pop ~ GDP + emp + invest + income + edu + health,
                  data = popdf,
                  listw = listw_popW, zero.policy = TRUE, family = "CAR")
summary(m2_car)
set.seed(8)
moran.mc(resid(m2_car), listw_popW, zero.policy = TRUE, nsim=999)
