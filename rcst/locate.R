# Use R --slave --args name axes [legend] < locate.R

args = commandArgs()

name = args[4]

x = 6.2
y = 3

data <- read.csv(file = paste(name, ".csv", sep = ""), head = FALSE, sep = ";", dec = ".", check.names = FALSE)
pdf(file = paste(name, ".pdf", sep = ""), width = x, height = y, paper = "special",
  family = "Helvetica", pointsize = 11)
par(mar=c(4, 4, 1, 1))

xrange = c(0, 10)
xscale = c(0, 2, 4, 6, 8, 10)
xtitle = ""
xlabs = F

yrange = c(100, 30000)
yscale = c(100, 1000, 10000)
ytitle = ""
ylabs = F

if(grepl("x", args[5]))
{
  xtitle = "Size (bpc)"
  xlabs = xscale
}
if(grepl("y", args[5]))
{
  ytitle = "Time (s)"
  ylabs = yscale
}

plot(c(1),
  c(1),
  type = "n",
  axes = F,
  main = "",
  xlab = xtitle,
  ylab = ytitle,
  xlim = xrange,
  ylim = yrange,
  log = "y")

axis(1, at = xscale, lab = xlabs, cex.axis = 0.8)
axis(2, at = yscale, lab = ylabs, cex.axis = 0.8)
box()

nc = ncol(data)
labs = c(1, 3)
points(data[1, 2:nc], data[2, 2:nc], type = "b", pch = 03, cex = 1.5)  # SSA
points(data[3, 2:nc], data[4, 2:nc], type = "b", pch = 04, cex = 1.5)  # RFM

if(length(args) > 5)
{
  legend(max(xrange), max(yrange), xjust = 1, yjust = 1, data[labs, 1],
  pch = c(03, 04), pt.cex = 1.5)
}

dev.off()
q()
