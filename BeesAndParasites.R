#####################################################################
#Section 1.4
#####################################################################
#Set the working directory and load the data
wd <- "D:/InterestingStuff/statistics/A beginners guide to GLM and GLMM with R/GLM Bees and parasites"
wd <- "C:/Documenten/QRM/Rscripts/"
setwd(wd)
Bees <- read.table(file = 'workerbees.txt', 
                   header = TRUE)
names(Bees)



##################################################################
#House keeping
#Convert the data to absence/presence data
Bees$Parasites[Bees$Parasites>0] <- 1
##################################################################



##################################################################
table(Bees$Parasites)

#Plot the data and add a linear regression line
#Figure 1.17
points(x = Bees$CellSize,
     y = Bees$Parasites,
     xlab = "Cell size",
     ylab = "Absence or presence of parasites",
     pch = 19)

M0 <- lm(Parasites ~ CellSize, 
         data = Bees)
abline(M0, lwd = 3)


#Boxplot of cell size conditional on absence/presence of parasites
boxplot(CellSize ~ Parasites,
        data = Bees)
##################################################################






##################################################################
M1 <- glm(Parasites ~ CellSize, data = Bees, family = binomial)
summary(M1)


#Figure 1.18
plot(x=Bees$CellSize, y = Bees$Parasites,
     main = "Logit model on the original data")
range(Bees$CellSize)
MyData <- data.frame(CellSize = seq(0.35,0.66, length = 50) )
prediction.with.se <- predict(M1, newdata=MyData, type = "response", se.fit = TRUE)
MyData$P <- prediction.with.se$fit

# confidence bands as 1.96 times the standard error
CIhigh <- MyData$P + 1.96 * prediction.with.se$se.fit
CIlow <- MyData$P - 1.96 * prediction.with.se$se.fit
polygon(x = c(MyData$CellSize, rev(MyData$CellSize)),
        y = c(CIhigh, rev(CIlow)), col = "light grey",
        border = FALSE)

# plot the model
lines(MyData$CellSize, MyData$P)

# CI for predictions
# MyData$P <- predict(M1, newdata=MyData, type = "response")

#For Figure 1.19: just change the link in the family argument,
#and rerun



#############################################################
M2 <- glm(Parasites ~ CellSize, data = Bees, family = quasibinomial)
summary(M2)
#See Hilbe (2009)
#############################################################


##################### BART ##################################
# what would a kernel density approach give?

# First with a histogram of all cells and one of only the infected ones
# Deviding the latter by the first
hist.all <- hist(Bees$CellSize, breaks = seq(0.35,0.67, by = 0.02))
plot(hist.all)
hist.para <- hist(Bees$CellSize[Bees$Parasites == 1], breaks = seq(0.35,0.67, by = 0.02))
plot(hist.para, add = TRUE, col = "red")

plot(x = hist.all$mids, y = hist.para$counts / hist.all$counts)
lines(MyData$CellSize, MyData$P)

# with a smoothed kernel
jpeg("BeesAndParasites_Impact_of_kernel_bandwidth.jpg")
plot(MyData$CellSize, MyData$P, col = "red", lwd = 3, type = "l",
     ylim = c(0, 1), xlab = "Cell size", ylab = "Probability on parasites",
     main = "Impact of kernel bandwith")
points(x=Bees$CellSize, y = Bees$Parasites)
CellSize.vec <- seq(0.35,0.67, by = 0.001)
bw.list <- c(0.005, 0.02, 0.05)
bw.colors <- c("blue", "green", "orange")

for (i in 1:length(bw.list)) {
  bw <- bw.list[i]
  bw.color <- bw.colors[i]
  pd.all <- rep(0, times = length(CellSize.vec))
  pd.parasites <- rep(0, times = length(CellSize.vec))
  for (i in 1:NROW(Bees)) {
    smoothed.point <- dnorm(CellSize.vec, mean = Bees$CellSize[i], sd = bw)
    pd.all <- pd.all + smoothed.point
    if (Bees$Parasites[i] == 1) {
      pd.parasites <- pd.parasites + smoothed.point
    }
  }
  # add the line to the graph
  lines(CellSize.vec, pd.parasites / pd.all, col = bw.color)
}
legend("bottomright", 
       legend = c("data", "logit GLM", paste0("bw = ", bw.list)),
       col = c("black", "red", bw.colors), lty=c(NA, 1, 1, 1, 1),
       pch = c(1, rep(NA, 4)))
dev.off()

# confidence intervals
confint(M1)

###################### half the data ##############################
Bees.half <- Bees[rbinom(NROW(Bees), size = 1, prob = 0.5) == 1,]
table(Bees.half$Parasites)
M1.half <- glm(Parasites ~ CellSize, data = Bees.half, family = binomial)
summary(M1.half)


#Figure 1.18
plot(x=Bees.half$CellSize, y = Bees.half$Parasites,
     main = "Logit model on half of the data")
range(Bees$CellSize)
MyData <- data.frame(CellSize = seq(0.35,0.66, length = 50) )
prediction.with.se.half <- predict(M1.half, newdata=MyData, type = "response", se.fit = TRUE)
MyData$P <- prediction.with.se.half$fit

# confidence bands as 1.96 times the standard error
CIhigh <- MyData$P + 1.96 * prediction.with.se.half$se.fit
CIlow <- MyData$P - 1.96 * prediction.with.se.half$se.fit
polygon(x = c(MyData$CellSize, rev(MyData$CellSize)),
        y = c(CIhigh, rev(CIlow)), col = "light grey",
        border = FALSE)

# plot the model
lines(MyData$CellSize, MyData$P)


###################### double the data ##############################
Bees.double <- rbind(Bees, Bees)
table(Bees.double$Parasites)
M1.double <- glm(Parasites ~ CellSize, data = Bees.double, family = binomial)
summary(M1.double)


#Figure 1.18
plot(x=Bees.double$CellSize, y = Bees.double$Parasites,
     main = "Logit model on double of the data")
range(Bees.double$CellSize)
MyData <- data.frame(CellSize = seq(0.35,0.66, length = 50) )
prediction.with.se.double <- predict(M1.double, newdata=MyData, type = "response", se.fit = TRUE)
MyData$P <- prediction.with.se.double$fit

# confidence bands as 1.96 times the standard error
CIhigh <- MyData$P + 1.96 * prediction.with.se.double$se.fit
CIlow <- MyData$P - 1.96 * prediction.with.se.double$se.fit
polygon(x = c(MyData$CellSize, rev(MyData$CellSize)),
        y = c(CIhigh, rev(CIlow)), col = "light grey",
        border = FALSE)

# plot the model
lines(MyData$CellSize, MyData$P)

