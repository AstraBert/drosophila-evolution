library(vcfR)
vcf <- read.vcfR('DEST.four_pops.neutral.vcf.gz')
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)

library(reshape2)
library(ggplot2) 
library(cowplot)

# Melt our matrix into a long form data.frame.
dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
dpf <- dpf[ dpf$Depth > 0,]

# Create a row designator.
# You may want to adjust this
#samps_per_row <- 20
samps_per_row <- 4
myRows <- ceiling(length(levels(dpf$Sample))/samps_per_row)
myList <- vector(mode = "list", length = myRows)

for(i in 1:myRows){
  myIndex <- c(i*samps_per_row - samps_per_row + 1):c(i*samps_per_row)
  myIndex <- myIndex[myIndex <= length(levels(dpf$Sample))]
  myLevels <- levels(dpf$Sample)[myIndex]
  myRegex <- paste(myLevels, collapse = "$|^")
  myRegex <- paste("^", myRegex, "$", sep = "")
  myList[[i]] <- dpf[grep(myRegex, dpf$Sample),]
  myList[[i]]$Sample <- factor(myList[[i]]$Sample)
}

# Create the plot.
myPlots <- vector(mode = "list", length = myRows)
for(i in 1:myRows){
  myPlots[[i]] <- ggplot(myList[[i]], aes(x=Sample, y=Depth)) + 
                  geom_violin(fill="#8dd3c7", adjust=1.0, scale = "count", trim=TRUE)

  myPlots[[i]] <- myPlots[[i]] + theme_bw()
  myPlots[[i]] <- myPlots[[i]] + theme(axis.title.x = element_blank(), 
                  axis.text.x = element_text(angle = 60, hjust = 1))
  myPlots[[i]] <- myPlots[[i]] + scale_y_continuous(breaks=c(1, 10, 100, 800),
                  minor_breaks=c(1:10, 2:10*10, 2:8*100))
  myPlots[[i]] <- myPlots[[i]] + theme( panel.grid.major.y=element_line(color = "#A9A9A9", size=0.6) )
  myPlots[[i]] <- myPlots[[i]] + theme( panel.grid.minor.y=element_line(color = "#C0C0C0", size=0.2) )
}

# Plot the plot.
plot_grid(plotlist = myPlots, nrow = myRows)
