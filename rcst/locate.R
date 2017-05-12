
# Use R --slave --args name axes [legend] < locate.R

args = commandArgs()

name = args[4]

x = 3.4
y = 3

data <- read.csv(file = paste(name, ".csv", sep = ""), head = FALSE, sep = ";", dec = ".", check.names = FALSE)
pdf(file = paste(name, ".pdf", sep = ""), width = x, height = y, paper = "special",
  family = "Helvetica", pointsize = 11)
par(mar = c(4, 4, 1, 1))

xrange = c(0, 10)
xscale = c(0, 2, 4, 6, 8, 10)
xtitle = ""
xlabs = F

yrange = c(0.3, 300)
yscale = c(1, 10, 100)
ytitle = ""
ylabs = F

if(grepl("x", args[5]))
{
  xtitle = "Size (bpc)"
  xlabs = xscale
}
if(grepl("y", args[5]))
{
  ytitle = "Time (µs / occurrence)"
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

axis(1, at = xscale, lab = xlabs)
axis(2, at = yscale, lab = ylabs)
box()

nc = ncol(data)
labs = c(2, 4, 6)
time_per_occ = 1000000.0 / 254997642

points(data[2, 2:nc], data[3, 2:nc] * time_per_occ, type = "b", pch = 01, cex = 1.5)  # SSA
points(data[4, 2:nc], data[5, 2:nc] * time_per_occ, type = "b", pch = 02, cex = 1.5)  # SSA-RRR
points(data[6, 2:nc], data[7, 2:nc] * time_per_occ, type = "b", pch = 00, cex = 1.5)  # RFM

if(length(args) > 5)
{
  legend(max(xrange), max(yrange), xjust = 1, yjust = 1, data[labs, 1],
  pch = c(01, 02, 00), cex = 0.8, pt.cex = 1.2)
}

dev.off()
q()
