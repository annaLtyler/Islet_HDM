#This function plots a network separating out defined modules
#set vertex.names to NA to remove the vertex labels

plot.modular.net <- function(net, modules, module.names = NULL, 
micro.layout = layout_nicely, macro.layout = layout_nicely, 
macro.layout.fun = "mean", cluster.layout.matrix = NULL, 
shiftx = 10, shifty = 10, vertex.col = NULL, vertex.size = 1, 
vertex.names = NA, edge.color = NULL, edge.split.vals = TRUE, edge.color.split = 0,
edge.colorscale = c("blue", "brown"), edge.grad.dir = c("ends"), 
vertex.label.dist = 1, vertex.label.cex = 1, edge.arrow.size = 1, 
edge.width = 1, edge.alpha = 0.5, highlight.nodes = NULL){

	if(is.null(V(net)$name)){stop("Network vertices must have names.")}

	require(RColorBrewer)
	mod.cols <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Set1"))
	
	u_modules <- sort(unique(modules))
	if(is.null(module.names)){
		mod.names <- u_modules
		}else{
		mod.names <- module.names	
		}
		
		
	if(is.null(cluster.layout.matrix)){
		layout.net <- module.weights(net, modules, summ.fun = macro.layout.fun)
		cluster.layout.matrix <- macro.layout(layout.net)
		# cluster.pane.matrix <- get.layout.mat(length(u_modules))
		}

	module.locale <- lapply(u_modules, function(x) which(modules == x))


	if(is.null(vertex.col)){
		node.col <- rep(mod.cols[1], vcount(net))
		for(m in 2:length(u_modules)){
			node.col[module.locale[[m]]] <- mod.cols[m]
			}
		}else{
		node.col <- vertex.col	
		}
		
		
	if(!is.null(highlight.nodes)){
		highlight.pos <- match(highlight.nodes, V(net)$name)
		node.col[highlight.pos] <- "red"
		}


	mod.layout <- vector(mode = "list", length = length(u_modules))
	names(mod.layout) <- u_modules
	for(i in 1:length(u_modules)){
		subnet <- induced_subgraph(net, v = module.locale[[i]])
		mod.layout[[i]]  <- micro.layout(subnet)
		rownames(mod.layout[[i]]) <- V(subnet)$name
		# pane.x.y <- which(cluster.layout.matrix == i, arr.ind = TRUE)
		pane.x.y <- cluster.layout.matrix[i,,drop=FALSE]
		mod.layout[[i]][,1] <- mod.layout[[i]][,1] + (shiftx * pane.x.y[,1])
		mod.layout[[i]][,2] <- mod.layout[[i]][,2] + (shifty * pane.x.y[,2])
		if(any(which(!is.finite(mod.layout[[i]])))){stop()}
		}	

	full.layout <- Reduce("rbind", mod.layout)
	full.layout.order <- match(V(net)$name, rownames(full.layout))
	full.layout <- full.layout[full.layout.order,]
	
	if(is.null(edge.color)){
		eweight <- E(net)$weight
		if(length(unique(eweight)) == 1){
			edge.col = rep("gray", ecount(net))
		}else{
			edge.col <- colors.from.values(vals = eweight, 
			col.scale = edge.color.scale, split.points = edge.color.split,
			grad.dir = edge.grad.dir, split.at.vals = edge.split.vals)
		}
	}else{
		edge.col <- edge.color
	}
	

	par(mar = c(0, 2, 2, 2))
	layout(matrix(c(1,2,3,2), nrow = 2, byrow = TRUE), heights = c(1,0.2), widths = c(1,0.2))
	# layout.show(3)

	#1: network	
	plot(net, layout = full.layout, vertex.color = node.col, 
	edge.color = edge.col, vertex.size = vertex.size, 
	vertex.label = vertex.names, vertex.label.dist = vertex.label.dist,
	vertex.label.cex = vertex.label.cex, edge.arrow.size = edge.arrow.size,
	edge.width = edge.width)

	#2: edge weight legend
	par(mar = c(3,3,3,3))
	if(is.null(edge.color)){
		imageWithTextColorbar(matrix(E(net)$weight, ncol = 1), col.scale = edge.color.scale, 
		axis.line = 0, cex = 1, grad.dir = edge.grad.dir, split.at.vals = edge.split.vals,
		split.points = edge.color.split)	
	}

	par(mar = c(2,0,0,0))
	plot.new()
	plot.window(xlim = c(0,1), ylim = c(0,1))
	legend(x = 0.2, y = 1, fill = mod.cols[1:length(u_modules)], legend = mod.names)	

	net <- set_edge_attr(net, "edge.color", value = edge.col)
	net <- set_vertex_attr(net, "vertex.color", value = node.col)
	net <- set_vertex_attr(net, "x", value = full.layout[,1])
	net <- set_vertex_attr(net, "y", value = full.layout[,2])
	
	invisible(net)
	
	}