
# Use R --slave --args name axes [legend] < synth.R

args = commandArgs()

name = args[4]

x = 3.4
y = 3

data <- read.csv(file = paste(name, ".csv", sep = ""), head = FALSE, sep = ";", dec = ".", check.names = FALSE)
pdf(file = paste(name, ".pdf", sep = ""), width = x, height = y, paper = "special",
  family = "Helvetica", pointsize = 11)
par(mar = c(4, 4, 1, 1))

xrange = c(0.0001, 0.1)
xscale = c("0.0001", "0.001", "0.01", "0.1")
xtitle = ""
xlabs = FALSE

yrange = c(0, 15)
yscale = c(0, 3, 6, 9, 12, 15)
ytitle = ""
ylabs = FALSE

if(grepl("x", args[5]))
{
  xtitle = "Mutation rate"
  xlabs = xscale
}
if(grepl("y", args[5]))
{
  ytitle = "Size (bpc)"
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
  log = "x")

axis(1, at = xscale, lab = xlabs)
axis(2, at = yscale, lab = ylabs)
box()

nr = nrow(data)
nc = ncol(data)
xpos = c(0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1)

symbols = c(00, 01, 02, 07)
if(nr == 2)
{
  symbols = c(07, 12)
}

for(i in c(1:nr))
{
  points(xpos, data[i, 2:nc], type = "b", pch = symbols[i], cex = 1.5)
}

if(length(args) > 5)
{
  legend(min(xrange), max(yrange), xjust = 0, yjust = 1, data[1:nr, 1],
  pch = symbols, pt.cex = 1.5)
}

dev.off()
q()
