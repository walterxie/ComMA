# Graph Internal Functions
# Author: Walter Xie
# Accessed on 21 Apr 2016


# df, x.id, y.id are compulsory, 
# fill.id, group.id, colour.id are optional
ggInit <- function(df, x.id, y.id=NULL, fill.id=NULL, group.id=NULL, colour.id=NULL, verbose=TRUE) {
  col.names <- colnames(df)
  if (!is.element(tolower(x.id), tolower(col.names)))
    stop("Data frame do NOT have column name \"", x.id, "\" !")

  suppressMessages(suppressWarnings(require(ggplot2)))
  aes.string <- paste0("aes(x=", x.id)
  if (! is.null(y.id)) {
    if (substring(y.id, 0, 2) != ".." && !is.element(tolower(y.id), tolower(col.names)))
      stop("Data frame do NOT have column name \"", y.id, "\" !")
    aes.string <- paste0(aes.string, ", y=", y.id)
  } 
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

# reformulate(fac2,fac1) => fac1 ~ fac2
ggOptFacetGrid <- function(p, col.names, x.facet.id=NULL, y.facet.id=NULL, 
                           scales = "fixed", space = "fixed", shrink = TRUE, 
                           drop = TRUE, verbose=TRUE) {
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
    if (verbose)
      cat("facet_grid(", deparse(reformulate(fac2,fac1)), ").\n")
    
    p <- p + facet_grid(reformulate(fac2,fac1), scales = scales, 
                        space = space, shrink = shrink, drop = drop)
  }
  return(p)
}

# if shape.id not NULL, no shape for points
# If NULL, the default, the data is inherited from the plot data.
ggOptPointAndShape <- function(p, col.names, shape.id=NULL, data=NULL, shapes=NULL, 
                               point.alpha=1, point.size=3) {
  if (! is.null(shape.id)) {
    if (!is.element(tolower(shape.id), tolower(col.names)))
      stop("Data frame do NOT have column name \"", shape.id, "\" !")
    
    p <- p + geom_point(data=data, aes_string(shape=shape.id), 
                        size = point.size, alpha=point.alpha) 
    if (! is.null(shapes)) 
      p <- p + scale_shape_manual(values=shapes) 
  } else {
    p <- p + geom_point(data=data, size = point.size, alpha=point.alpha)
  }
  return(p)
}

# add text.id before ggOptText, such as df$row.names <- rownames(df)
# If check_overlap=TRUE, text that overlaps previous text in the same layer will not be plotted
ggOptText <- function(p, col.names, text.id=NULL, text.data=NULL, colour.id=NULL, 
                      text.size=3, text.hjust=-0.1, text.vjust=-0.2, 
                      text.alpha=0.5, text.avoid.overlap=TRUE, verbose=TRUE) {
  if (! is.null(text.id)) {
    if (!is.element(tolower(text.id), tolower(col.names)))
      stop("Data frame do NOT have column name \"", text.id, "\" !")
    
    aes.string <- paste0("aes(label=", text.id)
    if (! is.null(colour.id)) 
      aes.string <- paste0(aes.string, ", colour=", colour.id)
    
    if (verbose)
      cat("geom_text(", aes.string, ", ...\n")
    
    if (length(text.size) > 1) {
      # prevent Error: Aesthetics must be either length 1 or the same as the data (5): size
      if (is.null(text.data))
        stop("Data for text cannot be NULL if given size vector !")
      # geom_text(data=df, aes(size = text.size))
      aes.string <- paste0(aes.string, ", size=text.size)") 
      p <- p + geom_text(eval(parse(text = aes.string)), data=text.data, 
                         hjust=text.hjust, vjust=text.vjust, alpha=text.alpha, 
                         check_overlap = text.avoid.overlap)
    } else {
      aes.string <- paste0(aes.string, ")") 
      p <- p + geom_text(eval(parse(text = aes.string)), data=text.data, size=text.size, 
                         hjust=text.hjust, vjust=text.vjust, alpha=text.alpha, 
                         check_overlap = text.avoid.overlap)
    }
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
ggOptPalette <- function(p, scale.to="colour", palette=NULL, verbose=TRUE) {
  if (! is.null(palette)) {
    if (length(palette) == 1)
      scale.string <- paste0("scale_", scale.to, "_brewer(palette=palette)")
    else if (length(palette) <= 3)
      scale.string <- paste0("scale_", scale.to, "_gradientn(colours=palette)")
    else 
      scale.string <- paste0("scale_", scale.to, "_manual(values=palette)")
    
    if (verbose)
      cat("colour scale : ", scale.string, "\n")
    
    p <- p + eval(parse(text = scale.string))
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
  if (scale=="continuous") {
    if (trans=="log")
      # show 1 at log scale
      scale.string <- paste0(scale.string, ", trans=ComMA::mylog_trans(base=10, from=-0.3)")  
    else
      scale.string <- paste0(scale.string, ", trans=trans")     
  }
  
  if (! is.null(auto.scale.max)) {
    breaks <- ComMA::get_breaks_positive_values(auto.scale.max, start=breaks.start)
    scale.string <- paste0(scale.string, ", labels=ComMA::scientific_10)") 
  } else {
    if (trans=="per") {
      suppressMessages(suppressWarnings(require(scales)))
      trans="identity"
      scale.string <- paste0(scale.string, ", labels=percent_format())") 
    } else {
      scale.string <- paste0(scale.string, ")") 
    }
  }
  
  if (verbose) {
    cat("Axis", axis, "scale :", scale.string, "\n")
    cat("Axis", axis, "breaks :", paste(breaks, collapse = ","), "\n")
  }
  
  p <- p + eval(parse(text = scale.string))
  return(p)
}


# Setting limits on the coordinate system will zoom the plot, 
# but will not change the underlying data like setting limits on a scale will
ggOptCoordCartesian <- function(p, df, x.id, y.id, x.lim.cart=NULL, y.lim.cart=NULL, 
                                coord.flip=FALSE, verbose=TRUE) {
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
  
  if (coord.flip) {
    p <- p + coord_flip()
    if (verbose)
      cat("Flip cartesian coordinates.\n")
  }
  return(p)
}

ggOptLegend <- function(p, legend.title.fill=NULL, legend.title.colour=NULL, 
                        legend.title.shape=NULL, legend.title.group=NULL,
                        legend.title.size=NULL, legend.col=1, legend.row=0, 
                        no.legend=NULL, verbose=TRUE) {
  if (!is.null(legend.title.fill))
    p <- p + labs(fill=legend.title.fill)
  if (!is.null(legend.title.colour))
    p <- p + labs(colour=legend.title.colour)
  if (!is.null(legend.title.shape))
    p <- p + labs(shape=legend.title.shape)
  if (!is.null(legend.title.group))
    p <- p + labs(group=legend.title.group, 
                  size=legend.title.size)
  if (!is.null(legend.title.size))
    p <- p + labs(size=legend.title.size)
  
  if (legend.col > 1 && legend.row > 0)
    warning("Cannot change legend.col and legend.row at the same time ! Skip both changes !")
  
  if (legend.col > 1 && legend.row == 0)
    p <- p + guides(fill=guide_legend(ncol=legend.col))
  if (legend.row > 0 && legend.col == 1) {
    p <- p + guides(fill=guide_legend(nrow=legend.row,byrow=TRUE))
  }
  
  if (!is.null(no.legend)) {
    guides.string <- paste0("guides(", no.legend, "=FALSE)")
    p <- p + eval(parse(text = guides.string))
    if (verbose) 
      cat("no.legend = ", no.legend, "\n")
  }
    
  return(p)
}

ggLabTitle <- function(p, x.id, y.id, title, x.lab=NULL, y.lab=NULL) {
  p <- p + theme_bw() + ggtitle(title)
  
  if (! is.null(x.lab)) 
    p <- p + xlab(x.lab) 
  
  if (! is.null(y.lab)) 
    p <- p + ylab(y.lab)
  
  return(p)
}

ggThemeOthers <- function(p, x.text=TRUE, y.text=TRUE, x.ticks=TRUE, y.ticks=TRUE, 
                          x.text.angle=0, verbose=TRUE, plot.margin.cm=NULL,
                          legend.position="right", legend.direction="vertical") {
  theme.string <- "theme(legend.position=legend.position, legend.direction=legend.direction"
  # hide x or y axis labels
  if (!x.text)
    theme.string <- paste(theme.string, "axis.text.x = element_blank()", sep = ",")
  if (!y.text)
    theme.string <- paste(theme.string, "axis.text.y = element_blank()", sep = ",")
  # rotate x labels
  if (x.text.angle > 0) 
    theme.string <- paste(theme.string, "axis.text.x = element_text(angle = x.text.angle, hjust = 1, vjust = 0.3)", sep = ",")
  # hide x or y axis ticks
  if (!x.ticks)
    theme.string <- paste(theme.string, "axis.ticks.x = element_blank()", sep = ",")
  if (!y.ticks)
    theme.string <- paste(theme.string, "axis.ticks.y = element_blank()", sep = ",")

  theme.string <- paste0(theme.string, ")")
  
  if (verbose)
    cat("theme : ", theme.string, "\n")
  
  p <- p + eval(parse(text = theme.string))
  
  if (! is.null(plot.margin.cm)) {
    suppressMessages(suppressWarnings(require(grid)))
    p <- p + theme(plot.margin=unit(plot.margin.cm,"cm"))
  }
  return(p)
}

ggThemePanelBorder <- function(p, title.size=10) {
  p <- p + theme(plot.title = element_text(size = title.size), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank()) 
  return(p)
}

# This funcion is broken in ggplot2, bug in axis.line
ggThemeAxis <- function(p, title.size=10) {
  p <- p + theme(plot.title = element_text(size = title.size), panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(), panel.background = element_blank(), 
                 panel.border = element_blank())
  #p <- p + theme(axis.line = element_line(size = 3, colour = "black")) 
  p <- p + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
  return(p)
}


