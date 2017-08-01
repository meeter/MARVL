######star
mystar <- function(coords, v=NULL, params) {
       vertex.color <- params("vertex", "color")
       if (length(vertex.color) != 1 && !is.null(v)) {
         vertex.color <- vertex.color[v]
       }
       vertex.size  <- 1/200 * params("vertex", "size")
       if (length(vertex.size) != 1 && !is.null(v)) {
         vertex.size <- vertex.size[v]
       }
			 vertex.frame.color <- params("vertex", "frame.color")
  		 if (length(vertex.frame.color) != 1 && !is.null(v)) {
    	 vertex.frame.color <- vertex.frame.color[v]
  		 }
			 vertex.frame.width <- params("vertex", "frame.width")
  		 if (length(vertex.frame.width) != 1 && !is.null(v)) {
    	 vertex.frame.width <- vertex.frame.width[v]
  		 }
       norays <- params("vertex", "norays")
       if (length(norays) != 1 && !is.null(v)) {
         norays <- norays[v]
       }
       symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color, lwd=vertex.frame.width,
          stars=matrix(c(vertex.size, vertex.size/2), nrow=1, ncol=norays*2),
          add=TRUE, inches=FALSE)
       }
add_shape("star", clip=shapes("circle")$clip, plot=mystar, parameters=list(vertex.norays=6, vertex.frame.width=3))

####triangle
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
	vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
	symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color, lwd=vertex.frame.width,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip, plot=mytriangle)

######fcircle
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color, lwd=vertex.frame.width,
          circles=vertex.size, add=TRUE, inches=FALSE)
}
add_shape("fcircle", clip=shapes("circle")$clip, plot=mycircle)


####frectangle
myrectangle <- function(coords, v=NULL, params) {

  vertex.color       <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.size        <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.size <- rep(vertex.size, length=nrow(coords))   
  vertex.size2       <- 1/200 * params("vertex", "size2")
  if (length(vertex.size2) != 1 && !is.null(v)) {
    vertex.size2 <- vertex.size2[v]
  }
  vertex.size <- cbind(vertex.size2, vertex.size)
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color, lwd=vertex.frame.width,
          rectangles=2*vertex.size, add=TRUE, inches=FALSE)
}
add_shape("frectangle", clip=igraph.shape.noclip, plot=myrectangle)
