#This function is an alternative to imageWithTextColorbar
#It has all the same arguments as colors.from.values, plus
#arguments about customizing the position and look of the
#colorbar.
#the bar can be positioned using the legend position words,
#such as "bottomright". The margin.buffer (distance from
#plot edge) and label.buffer (text shift) can be adjusted
#to fine-tune the position. These are both percentages of
#the total plot width or height.
#num.scale.pts determines how many numbers will be displayed
#on the scale.
#orientation determines whether the bar will be displayed
#vertically or horizontally. cex is used to adjust the size
#of the bar. With a cex of 1, the legend will be limited
#to one quadrant of the plot, but this can be contracted
#or expanded by adjusting cex up or down.
#as with legend, this function can only be used on an 
#existing plot.

add_color_bar <- function(values, position = "bottomright", num.scale.pts = 5, 
    margin.buffer = 0.05, label.buffer = 0.01, legend.title = "",
    orientation = c("v", "h"), split.at.vals = FALSE, split.points = 0, 
    col.scale = c("green", "purple", "orange", "blue", "brown", "gray"), light.dark = "f", 
    grad.dir = c("high", "low", "middle", "ends"), color.fun = c("linear", "exponential"), 
    exp.steepness = 1, global.color.scale = FALSE, global.min = NULL, global.max = NULL, 
    use.pheatmap.colors = FALSE, na.col = "lightgray", custom.colors = NULL){

            #set up some reasonable defaults
            orientation = orientation[1]
            num.lines = 100

            range.min <- min(values, na.rm = TRUE)
            range.max <- max(values, na.rm = TRUE)

            val.range <- segment.region(range.min, range.max, num.lines, "ends")

            range.col <- colors.from.values(val.range, split.at.vals = split.at.vals, 
                split.points = split.points, col.scale = col.scale, light.dark = light.dark, 
                grad.dir = grad.dir, color.fun = color.fun, exp.steepness = exp.steepness, 
                global.color.scale = global.color.scale, global.min = global.min, 
                global.max = global.max, use.pheatmap.colors = use.pheatmap.colors, 
                na.col = na.col, custom.colors = custom.colors)

            plot.dim <- par("usr")
            plot.height <- plot.dim[4] - plot.dim[3]
            plot.width <- plot.dim[2] - plot.dim[1]

            #setup default positions for the bar based on 
            #the position words
            on_left <- as.logical(length(grep("left", position)))
            on_bottom <- as.logical(length(grep("bottom", position)))

            if(on_bottom){ #bottom half; y position relative to minimum of plot
                ymin <- plot.dim[3] + (plot.height*margin.buffer)
                ymax <- mean(plot.dim[c(3:4)]) - (plot.height*margin.buffer)
            }else{ #top half; 
                ymin <- mean(plot.dim[c(3:4)]) + (plot.height*margin.buffer)
                ymax <- plot.dim[4] - (plot.height*margin.buffer)
            }

            if(on_left){ #left half; x position relative to minimum of plot
                xmin <- plot.dim[1] + (plot.width*margin.buffer)
                xmax <- mean(plot.dim[c(1:2)]) - (plot.width*margin.buffer)
            }else{ #on right
                xmin <- mean(plot.dim[c(1:2)]) + (plot.height*margin.buffer)
                xmax <- plot.dim[2] - (plot.width*margin.buffer)
            }
    

            #if we have a vertical orientation, the x position will 
            #be the min x position if on the left and the max if on 
            #the right
            #the y position will be a sequence from the min to the max
            #the axis labels will be to the left of the points
            if(orientation == "v"){
                if(on_left){
                    xpos <- rep(xmin, length(range.col))
                }else{
                    xpos <- rep(xmax, length(range.col))
                }
                ypos <- segment.region(ymin, ymax, length(range.col), "ends")
                scalex <- xpos[1:num.scale.pts] - (plot.width*label.buffer)
                scaley <- segment.region(ymin, ymax, num.scale.pts, "ends")
                text(x = xpos[1], y = max(ypos) + (plot.height*label.buffer*2), 
                    labels = legend.title)
            }else{
                #if we have a horizontal orientation, the y position will 
                #be the min y position if on the bottom and the max if on 
                #the top
                #the x position will be a sequence from the min to the max
                #the axis labels will be below the points
                if(on_bottom){
                    ypos <- rep(ymin, length(range.col))
                }else{
                    ypos <- rep(ymax, length(range.col))
                }
                xpos <- segment.region(xmin, xmax, length(range.col), "ends")
                scaley <- ypos[1:num.scale.pts] - (plot.width*label.buffer)
                scalex <- segment.region(xmin, xmax, num.scale.pts, "ends")
                text(x = mean(xpos), y = ypos[1] + (plot.width*label.buffer*2), 
                    labels = legend.title)
            }
          
            points(xpos, ypos, col = range.col, pch = "-", cex = 2)
            scale.labels <- signif(segment.region(range.min, range.max, num.scale.pts, "ends"), 2)
            text(scalex, scaley, labels = scale.labels, adj = 1, cex = 0.7)

}