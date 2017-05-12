
# Use R --slave --args name axes [legend] < locate.R

args = commandArgs()

name = args[4]

x = 3.4
y = 3

data <- read.csv(file = paste(name, ".csv", sep = ""), head = FALSE, sep = ";", dec = ".", check.names = FALSE)
pdf(file = paste(name, "_bars.pdf", sep = ""), width = x, height = y, paper = "special",
  family = "Helvetica", pointsize = 11)
par(mar = c(4, 4, 1, 1))

nc = ncol(data)

xtitle = NULL

yrange = c(0.3, 300)
yscale = c(1, 10, 100)
ytitle = ""
ylabs = FALSE

time_per_occ = 1000000.0 / 254997642

if(grepl("x", args[5]))
{
  xtitle = "Sample interval"
}
if(grepl("y", args[5]))
{
  ytitle = "Time (µs / occurrence)"
  ylabs = yscale
}

labs = c(3, 5, 7)
times = cbind(t(data[labs, 2:nc])) * time_per_occ
colnames(times) = data[labs, 1]
rownames(times) = data[1, 2:nc]
times

colors = gray.colors(3)
bp = barplot(t(times),
  beside = TRUE,
  col = colors,
  axes = FALSE,
  xlab = xtitle,
  ylab = ytitle,
  ylim = yrange,
  log = "y")

axis(side = 2, at = ylabs, lab = ylabs)
box()

if(length(args) > 5)
{
  legend(1, 0.75 * max(yrange), xjust = 0, yjust = 1, colnames(times), fill = colors, cex = 0.8)
}

dev.off()
q()
