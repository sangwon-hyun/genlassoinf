set.seed(0)
n = 7
rep = c(4,3,4,2,3,3,4)
lev = c(100,93,86,88,91,83,86)+rnorm(n,0,1)
y = rep(lev,rep)
x = 1:length(y)

w = 5; h = 5
outputdir = "output/"
mar = c(4.5,4.5,0.5,0.5)

# Plot, repeated IC case
pdf(file=file.path(outputdir, "ic1.pdf"),width=w,height=h)
par(mar=mar)
plot(x, y, xlab="Step", ylab="IC", main="")
i = c(1,which(diff(y)!=0)+1)
points(x[i], y[i], col="blue", pch=19, cex=1.1)
points(x[8], y[8], col="red", cex=2)
legend("topright", col=c("black","blue","red"),
       pch=c(21,19,21), 
       legend=c("All steps","Candidate steps","IC-selected step"))
graphics.off()

set.seed(0)
n = 25
x = 1:n
y = (x - 18)^2
y = (y - min(y))/max(y)*20 + 80
y = y + rnorm(n,0,1)
plot(x,y)

# Plot, non-repeated IC case
pdf(file=file.path(outputdir, "ic2.pdf"),width=w,height=h)
par(mar=mar)
plot(x, y, xlab="Step", ylab="IC", main="")
i = 1:n
points(x[i], y[i], col="blue", pch=19, cex=1.1)
j = which(diff(y)>0)
m = j[min(which(diff(j)==1))]
points(x[m], y[m], col="red", cex=2)
legend("topright", col=c("black","blue","red"),
       pch=c(21,19,21), 
       legend=c("All steps","Candidate steps","IC-selected step"))
graphics.off()
