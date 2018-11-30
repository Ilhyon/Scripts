myData <- aggregate(mtcars$mpg,
                    by = list(cyl = mtcars$cyl, gears = mtcars$gear),
                    FUN = function(x) c(mean = mean(x), sd = sd(x),
                                        n = length(x)))

myData <- do.call(data.frame, myData)
myData$se <- myData$x.sd / sqrt(myData$x.n)
colnames(myData) <- c("cyl", "gears", "mean", "sd", "n", "se")
myData$names <- c(paste(myData$cyl, "cyl /",
                        myData$gears, " gear"))
########################################
# Recreate ugly non-grouped barplot
# using ggplot2
########################################

# load in ggplot2 library
#install.packages("ggplot2")
library(ggplot2)

# Set spacing between bars
dodge <- position_dodge(width = 0.9)

# Set y limits
limits <- aes(ymax = myData$mean + myData$se,
              ymin = myData$mean - myData$se)

# Construct ggplot object
p <- ggplot(data = myData, aes(x = names, y = mean, fill = names))

# Add in bars, error bars, and fiddle with
# axis names and ticks and titles
p + geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  labs(x = "No. Cylinders and Gears", y = "Miles Per Gallon") +
  ggtitle("Mileage by No. Cylinders\nand No. Gears")


########################################
# Group our bar graphs in ggplot2
########################################

# Construct ggplot object specifying cylinders
# and gears as factors to make sure they group
p <- ggplot(data = myData, aes(x = factor(cyl), y = mean,
                               fill = factor(gears)))

# Same stuff as last time: add in bars and
# error bars and titles and labels
p + geom_bar(stat = "identity",
             position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.25) +
  labs(x = "No. Cylinders", y = "Miles Per Gallon") +
  ggtitle("Mileage by No. Cylinders\nand No. Gears") +
  scale_fill_discrete(name = "No. Gears")

