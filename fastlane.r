#-----------------------------------------------------------------chapters 1-9
# -------------------------------------------------------------------------------------------------------------------------------------------
# mean(Nile)
print(Nile)
print(length(Nile))
print(mean(Nile))
hist(Nile) # creates pdf file in the same directory as the program file
#print(hist(Nile)) # не работает в терминале
# ?hist
print(Nile[2])
print(Nile[c(2,5,6)]) # c function concatenates vectors
print(Nile[2:4])
print(mean(Nile[81:100]))
n81100 <- Nile[81:100] # variable assignment
print(mean(n81100))
print(sd(n81100))
print(sum(c(5,12,13)))
# --------------------------------------
x <- c(5,12,13)
x > 8
sum(x > 8)
# --------------------------------------
# way 1
print(sum(Nile > 1200)) # how many years had a flow above 1200
which(Nile > 1200) # which years had a flow above 1200 - shows indeces, not years
# way 2
which1200 <- which(Nile > 1200)
which1200
length(which1200)
Nile[which1200] # what were the river flows in those 7 years
# or
Nile[Nile > 1200]
# ---------------------------------------
x <- c(5,12,13,8)
x[-1] # Here we are asking for all of x except for x[1]
x[c(-1,-4)]
# ---------------------------------------
# dataframes
# ?ToothGrowth                                                     
head(ToothGrowth)
tg <- ToothGrowth
mean(tg$len) # Dollar signs are used to denote the individual columns, e.g. ToothGrowth$dose for the dose column
mean(tg[,1]) # Some data frames don't have column names, but that is no obstacle. We can use column numbers
tg[3,1] # To get the element in row 3, column 1
z <- tg[2:5,c(1,3)] # rows 2 through 5, and columns 1 and 3, assigning the result to z
z
nrow(ToothGrowth) # nrow function tells us the number of rows in any data frame
head(ToothGrowth$len) # head function works on vectors too
head(ToothGrowth$len,10) # second argument, specifying how many elements to print
# creating a dataframem-----------------------------------------------------------------------------
x <- c(5,12,13)
y <- c('abc','de','z')
d <- data.frame(x,y)
d
# ------------------------------------------------
z <- tg[,-2] # last 2 columns
z
# -------------------------------------------------------------------------------------------------------
# classes
class(tg)
class(tg$supp)
# R factors are used when we have categorical variables
levels(tg$supp) # list of categories for tg$supp as follows
whichOJ <- which(tg$supp == 'OJ')
whichVC <- which(tg$supp == 'VC')
mean(tg[whichOJ,1])
mean(tg[whichVC,1])
# or
tgoj <- tg[tg$supp == 'OJ',]
tgvc <- tg[tg$supp == 'VC',]
mean(tgoj$len)
mean(tgvc$len)
#
tg[tg$supp=='OJ' & tg$len < 8.8,]
tg[tg$len > 28 | tg$dose == 1.0,]
endose <- tg[,c(1,3)]
endose
lendose <- tg[,c('len','dose')]
head(lendose)
exts <- Nile[Nile < 800 | Nile > 1300]
head(exts)
length(Nile[Nile < 800 | Nile > 1300])
tapply(tg$len,tg$supp,mean) # Split the vector tg$len into two groups, according to the value of tg$supp, then apply mean to each group
z <- tapply(tg$len,tg$supp,mean)
z[1]
# ---------------------------------------------------------------
x <- c(8,5,12,13)
g <- c('M',"F",'M','M')
tapply(x,g,mean) # Split x into two piles, according to the corresponding elements of g, and then find the mean in each pile
tapply(PlantGrowth$weight,PlantGrowth$group,length)
#
head(mtcars)
row.names(mtcars)
row.names(mtcars)[7]
row.names(mtcars)[7] <- 'Dustpan'
row.names(mtcars)[7]

# ------------------------------------------------------------chapters 10-18
#------------------------------------------------------------------------------------------------------------------------------------
# pima <- read.csv('http://heather.cs.ucdavis.edu/FasteR/data/Pima.csv', header=TRUE)
pima <- read.csv('Pima.csv', header=TRUE)
head(pima)
dim(pima) # dim function tells us that there are 768 people in the study, 9 variables measured on each

table(pima$glucose)
# the first, third, fifth and so on lines are the glucose values
# second, fourth, sixth and so on lines are the counts of women having those values
pg <- pima$glucose
pg1 <- pg[pg > 0]
length(pg1)
mean(pg)
mean(pg1)

pima$glucose[pima$glucose == 0] <- NA # recode zeros as NAs
# break things up into smaller steps
glc <- pima$glucose # first line just makes a copy of the original vector, to avoid clutter in the code
z <- glc == 0 # determines which elements of glc are 0s, resulting in z being a vector of TRUEs and FALSEs
glc[z] <- NA # assigns NA to those elements in glc corresponding to the TRUEs
pima$glucose <- glc # we need to have the changes in the original data, so we copy glc to it

sum(is.na(pima$glucose)) # let's verify that we now have 5 NAs in the glucose variable
mean(pima$glucose, na.rm = TRUE) # instruct the function to skip the NAs

mtmpg <- mtcars$mpg
mtmpg
mtcars$cyl == 4
mt4 <- mtmpg[mtcars$cyl == 4]
# a cleaner way
mtl <- split(mtmpg,mtcars$cyl)
mtl
mtl[[1]]
head(mtcars$cyl)

l <- list(a = c(2, 5), b = "sky")
l

names(mtl) <- c('four', 'six', 'eight')
mtl
mtl[[2]][3]
mtl$six[3]

plot(Nile)
which(Nile < 600)
plot(mtcars$wt, mtcars$mpg)
is.numeric(Nile)
class(Nile)
Nile[1 + 1925 - 1871]

# load(url('https://github.com/matloff/fasteR/blob/master/data/prgeng.RData?raw=true'))
load("prgeng.RData")
head(prgeng)
plot(prgeng$age,prgeng$wageinc)

# indxs <- sample(1:nrow(prgeng),2500) # nrow() function returns the number of rows in the argument
# prgeng2500 <- prgeng[indxs,]
# plot(prgeng2500$age,prgeng2500$wageinc)
# plot(prgeng2500$age,prgeng2500$wageinc,col=prgeng2500$sex)
# plot(pe2500$age,pe2500$wageinc,col=as.factor(pe2500$sex),xlab='age',ylab='wage',cex=0.6) # 'cex = 0.6' means "Draw the dots at 60% of default size

wageByGender <- split(prgeng$wageinc,prgeng$sex)
dm <- density(wageByGender[[1]])
dw <- density(wageByGender[[2]])
plot(dw,col='red')
points(dm,cex=0.2)

gt1200 <- which(Nile > 1200)
nileSubsetGT1200 <- Nile[gt1200]
mean(nileSubsetGT1200)
mean(Nile[Nile > 1200])

mgd <- function(x,d) mean(x[x > d])
mgd(Nile,1200)
save(mgd,file='mean_greater_than_d')
# load('mean_greater_than_d')

rng <- function(y) max(y) - min(y)
rng(Nile)

for (i in 1:9) print(sum(pima[,i] == 0))
colnames(pima)

for (i in 2:6) {
    zeroIndices <- which(pima[,i] == 0)
    pima[zeroIndices,i] <- NA
}

# ------------------------------------------------------------chapters 19-27
#------------------------------------------------------------------------------------------------------------------------------------
nile <- ifelse(Nile > 1150,3,2)
nile <- ifelse(Nile < 800,1,nile)
table(nile)

load("mlb.RData")
head(mlb)
age <- round(mlb$Age)
table(age)
taout <- tapply(mlb$Weight,age,mean)
taout
names(taout)
plot(23:35,taout[3:15])

lm(Weight ~ Age,data=mlb) # linear regression
abline(181.4366,0.6936)

lm(Weight ~ Height + Age, data=mlb)

table(mlb$PosCategory)

rownums <- split(1:nrow(mlb),mlb$PosCategory)
str(rownums)

posNames <- c('Catcher','Infielder','Outfielder','Pitcher')
m <- data.frame()
for (pos in posNames) {
   posRows <- rownums[[pos]]
   lmo <- lm(Weight ~ Age, data = mlb[posRows,])
   newrow <- lmo$coefficients
   m <- rbind(m,newrow)
}
m
row.names(m) <- posNames
names(m) <- c('intercept','slope')
m

plot(mlb$Age,mlb$Weight,col=mlb$PosCategory)
with(mlb,plot(Age,Weight,col=PosCategory,cex=0.6))

library(ggplot2)
p <- ggplot(mlb) # make an empty plot, based on the data frame mlb
p + geom_point(aes(x = Age, y = Weight, col = PosCategory),cex=0.6) # geom_point() is a ggplot2 function. Its task is to produce scatter plots

# ------------------------------------------------------------chapters 28-36
#------------------------------------------------------------------------------------------------------------------------------------
zlm <- function(rws) lm(Weight ~ Age, data=mlb[rws,])$coefficients
w <- lapply(rownums,zlm)
w

apply(pima,2,max)

abt <- readLines('aboutR.txt')
str(abt)
head(abt)
paste('abc','987')
abt1 <- paste(abt,collapse=' ') 
str(abt1)
substr(abt1,288,336)
y <- strsplit(abt1,' ')
str(y)
y1 <- y[[1]]
y2 <- y1[y1 != '']

day <- read.csv('day.csv',header=TRUE)
head(day) 
day1 <- day[,c(8,10,12:14)]
head(day1)
lmout <- lm(casual ~ .,data=day1) # "casual ~ ." means, "regress casual against all the other variables in this dataset
lmout

summary(lmout)
newCase <- data.frame(workingday=1, temp=0.26, hum=0.55, windspeed=0.18)
predict(lmout,newCase)

hds <- day$dteday[day$holiday == 1]
hds
hd <- as.Date(hds)
hp <- as.POSIXlt(hd)
hp

# logistic model
glout <- glm(test ~ bmi + age, data=pima, family=binomial)
summary(glout)

# The '&&' operator stands for "and"
