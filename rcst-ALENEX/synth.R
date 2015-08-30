# Use R --slave --args name axes [legend] < synth.R

args = commandArgs()

name = args[4]

x = 3.3
y = 3

data <- read.csv(file = paste(name, ".csv", sep = ""), head = FALSE, sep = ";", dec = ".", check.names = FALSE)
pdf(file = paste(name, ".pdf", sep = ""), width = x, height = y, paper = "special",
  family = "Helvetica", pointsize = 11)
par(mar=c(4, 4, 1, 1))

xrange = c(0.0001, 0.1)
xscale = c("0.0001", "0.001", "0.01", "0.1")
xtitle = ""
xlabs = xscale

yrange = c(0, 15)
yscale = c(0, 3, 6, 9, 12, 15)
ytitle = ""
ylabs = yscale

if(grepl("x", args[5]))
{
  xtitle = "Mutation rate"
}
if(grepl("y", args[5]))
{
  ytitle = "Size (bpc)"
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
  log = "x")

axis(1, at = xscale, lab = xlabs, cex.axis = 0.8)
axis(2, at = yscale, lab = ylabs, cex.axis = 0.8)
box()

nr = nrow(data)
nc = ncol(data)
xpos = c(0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1)
symbols = c(04, 03, 08)
if(nr < 3)
{
  symbols = c(08, 03)
}

points(xpos, data[1, 2:nc], type = "b", pch = symbols[1], cex = 1.5)
points(xpos, data[2, 2:nc], type = "b", pch = symbols[2], cex = 1.5)

if(nr >= 3)
{
  points(xpos, data[1, 2:nc] + data[2, 2:nc], type = "b", pch = symbols[3], cex = 1.5)  # RLCP
}

if(length(args) > 5)
{
  legend(min(xrange), max(yrange), xjust = 0, yjust = 1, data[1:nr, 1],
  pch = symbols, pt.cex = 1.5)
}

dev.off()
q()
