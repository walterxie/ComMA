# Graph Internal Functions
# Author: Walter Xie
# Accessed on 21 Apr 2016


# df, x.id, y.id are compulsory, 
# fill.id, group.id, colour.id are optional
ggInit <- function(df, x.id, y.id, fill.id=NULL, group.id=NULL, colour.id=NULL, verbose=TRUE) {
  col.names <- colnames(df)
  if (!is.element(tolower(x.id), tolower(col.names)))
    stop("Data frame do NOT have column name \"", x.id, "\" !")
  if (!is.element(tolower(y.id), tolower(col.names)))
    stop("Data frame do NOT have column name \"", y.id, "\" !")
  
  suppressMessages(require(ggplot2))
  aes.string <- paste0("aes(x=", x.id, ", y=", y.id)
  if (! is.null(fill.id)) {
    if (!is.element(tolower(fill.id), tolower(col.names)))
      stop("Data frame do NOT have column name \"", fill.id, "\" !")
    
    aes.string <- paste0(aes.string, ", fill=", fill.id)
  } 
  if (! is.null(group.id)) {
    if (!is.element(tolower(group.id), tolower(col.names)))
      stop("Data frame do NOT have column name \"", group.id, "\" !")
    
    aes.string <- paste0(aes.string, ", group=", group.id)
  } 
  if (! is.null(colour.id))  {
    if (!is.element(tolower(colour.id), tolower(col.names)))
      stop("Data frame do NOT have column name \"", colour.id, "\" !")
    
    aes.string <- paste0(aes.string, ", colour=", colour.id)
  } 
  aes.string <- paste0(aes.string, ")")
  
  if (verbose)
    cat("ggplot aesthetic mappings are ", aes.string, ".\n")
  
  p <- ggplot(df, eval(parse(text = aes.string)))
  return(p)
}

ggOptFacetGrid <- function(p, col.names, x.facet.id=NULL, y.facet.id=NULL) {
  if (! is.null(x.facet.id) || ! is.null(y.facet.id)) {
    fac1 <- "."
    fac2 <- "."
    if (! is.null(x.facet.id)) {
      if (!is.element(tolower(x.facet.id), tolower(col.names)))
        stop("Data frame do NOT have column name \"", x.facet.id, "\" !")
      
      fac1 <- x.facet.id
    }
    if (! is.null(y.facet.id)) {
      if (!is.element(tolower(y.facet.id), tolower(col.names)))
        stop("Data frame do NOT have column name \"", y.facet.id, "\" !")
      fac2 <- y.facet.id
    }
    p <- p + facet_grid(reformulate(fac2,fac1))
  }
  return(p)
}

# if shape.id not NULL, no shape for points
ggOptPointAndShape <- function(p, col.names, shape.id=NULL, data=NULL, shapes=NULL, point.size=3) {
  if (! is.null(shape.id)) {
    if (!is.element(tolower(shape.id), tolower(col.names)))
      stop("Data frame do NOT have column name \"", shape.id, "\" !")
    
    p <- p + geom_point(data=data, aes_string(shape=shape.id), size = point.size) 
    if (! is.null(shapes)) 
      p <- p + scale_shape_manual(values=shapes) 
  } else {
    p <- p + geom_point(data=data, size = point.size)
  }
  return(p)
}

# normally ellipsed.id == colour.id to show clusters
ggOptEllipse <- function(p, col.names, ellipsed.id=NULL) {
  if (! is.null(ellipsed.id)) {
    if (!is.element(tolower(ellipsed.id), tolower(col.names)))
      stop("Data frame do NOT have column name \"", ellipsed.id, "\" !")
    p <- p + stat_ellipse(aes_string(mapping = ellipsed.id), type = "t", linetype = 2)
  }
  return(p)
}

# scale_*_brewer, * has "colour", "fill" 
ggOptPalette <- function(p, scale.to="colour", palette=NULL) {
  if (! is.null(palette)) {
    if (length(palette) == 1)
      scale.string <- paste0("scale_", scale.to, "_brewer(palette=palette)")
    else if (length(palette) <= 3)
      scale.string <- paste0("scale_", scale.to, "_gradientn(colours=palette)")
    else 
      scale.string <- paste0("scale_", scale.to, "_manual(values=palette)")
    
    p <- p + eval(parse(text = scale.string))
  }
  return(p)
}

# add text.id before ggOptText, such as df$row.names <- rownames(df)
ggOptText <- function(p, col.names, text.id=NULL, text.data=NULL, colour.id=NULL, text.size=3, 
                      text.hjust=-0.1, text.vjust=-0.2, text.alpha=0.5) {
  if (! is.null(text.id)) {
    if (!is.element(tolower(text.id), tolower(col.names)))
      stop("Data frame do NOT have column name \"", text.id, "\" !")
    
    aes.string <- paste0("aes(label=", text.id)
    if (! is.null(colour.id)) 
      aes.string <- paste0(aes.string, ", colour=", colour.id)
    aes.string <- paste0(aes.string, ")") 
    p <- p + geom_text(data=text.data, eval(parse(text = aes.string)), size=text.size, 
                       hjust=text.hjust, vjust=text.vjust, alpha=text.alpha)
  }
  return(p)
}

# 1) auto.scale.x.max = max(df[,x.id]) to determine positive breaks automatically,
# breaks=c(start, 10, 100, ...), default start=c(0.1, 1) in get_breaks_positive_values.
# 2) breaks.interval > 0, seq(interval.min, interval.max, breaks.interval)
ggOptScaleAxis <- function(p, axis="y", scale="continuous", trans="identity", 
                           expand=c(0,0), breaks=waiver(), labels=waiver(), 
                           auto.scale.max=NULL, breaks.start=c(0), verbose=TRUE) {
  if ( ! is.element(axis, c("x", "y")) )
    stop("Incorrect axis ", axis, ", use x or y !")
  if ( ! is.element(scale, c("continuous", "discrete")) )
    stop("Incorrect scale ", scale, ", use continuous or discrete !")
  
  # scale_y_continuous(trans=trans, expand=expand, breaks=breaks, labels=labels)
  scale.string <- paste0("scale_", axis, "_", scale, "(expand=expand, breaks=breaks")
  if (scale=="continuous")
    scale.string <- paste0(scale.string, ", trans=trans") 
  if (! is.null(auto.scale.max)) {
    breaks <- ComMA::get_breaks_positive_values(auto.scale.max, start=breaks.start)
    scale.string <- paste0(scale.string, ", labels=ComMA::scientific_10)") 
  } else {
    if (trans=="per") {
      suppressMessages(require(scales))
      trans="identity"
      scale.string <- paste0(scale.string, ", labels=percent_format())") 
    } else {
      scale.string <- paste0(scale.string, ")") 
    }
  }
  
  if (verbose) {
    cat("Axis", axis, "scale :", scale.string, ".\n")
    cat("Axis", axis, "breaks :", paste(breaks, collapse = ","), ".\n")
  }
  
  p <- p + eval(parse(text = scale.string))
  return(p)
}


# Setting limits on the coordinate system will zoom the plot, 
# but will not change the underlying data like setting limits on a scale will
ggOptCoordCartesian <- function(p, df, x.id, y.id, x.lim.cart=NULL, y.lim.cart=NULL) {
  if (! is.null(x.lim.cart)) {
    if (length(x.lim.cart) != 2)
      stop("Invalid format, use x.lim.cart = c(1000,NA) !")
    if (which(is.na(x.lim.cart))==1) {
      x.lim.cart[1] = min(df[,x.id])
    } else if (which(is.na(x.lim.cart))==2) {
      x.lim.cart[2] = max(df[,x.id])
    }
    p <- p + coord_cartesian(xlim=x.lim.cart)
  } 
  if (! is.null(y.lim.cart)) {
    if (length(y.lim.cart) != 2)
      stop("Invalid format, use y.lim.cart = c(1000,NA) !")
    if (which(is.na(y.lim.cart))==1) {
      y.lim.cart[1] = min(df[,y.id])
    } else if (which(is.na(y.lim.cart))==2) {
      y.lim.cart[2] = max(df[,y.id])
    }
    p <- p + coord_cartesian(ylim=y.lim.cart)
  } 
  return(p)
}

ggOptLegend <- function(p, legend.title=NULL, legend.col=1, legend.nrow=0) {
  if (!is.null(legend.title))
    p <- p + labs(fill=legend.title)
  
  if (legend.col > 1 && legend.nrow > 0)
    warning("Cannot change legend.col and legend.nrow at the same time ! Skip both changes !")
  
  if (legend.col > 1 && legend.nrow == 0)
    p <- p + guides(fill=guide_legend(ncol=legend.col))
  if (legend.nrow > 0 && legend.col == 1)
    p <- p + guides(fill=guide_legend(nrow=legend.nrow,byrow=TRUE))
  return(p)
}

ggLabTitle <- function(p, x.id, y.id, title, x.lab="x.id", y.lab="y.id") {
  if (x.lab=="x.id") 
    x.lab = x.id
  if (y.lab=="y.id") 
    y.lab = y.id
  
  p <- p + theme_bw() + xlab(x.lab) + ylab(y.lab) + ggtitle(title)
  
  return(p)
}

ggThemeRotateXText <- function(p, x.text.angle=0) {
  if (x.text.angle > 0) 
    p <- p + theme(axis.text.x = element_text(angle = x.text.angle, hjust = 1, vjust = 0.3))
  return(p)
}

ggThemePanelBorder <- function(p, title.size=10) {
  p <- p + theme(plot.title = element_text(size = title.size), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank()) 
  return(p)
}

ggThemeAxis <- function(p, title.size=10) {
  p <- p + theme(plot.title = element_text(size = title.size), panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(), panel.background = element_blank(),
                 axis.line = element_line(size = 3, colour = "black"), panel.border = element_blank()) 
  return(p)
}


