
circ_mean <- function(x, # bearings -- in degrees, not radians
                      w=NULL,
                      return = "angle",
                      ...){
        # a custom weighted version of the "mean resultant vector" from circular stats
        # see: https://doi.org/10.3389/fpsyg.2018.02040 and help(CircStats::circ.disp)
        # my implementation: convert the bearings into XY, then take the mean of the XYs, weight by wind speed
        
        require(geosphere)
        if(is.null(w)) w <- rep(1, length(x))
        x <- x / 180 * pi
        xy <- cbind(x = sin(x), y = cos(x))
        xy <- apply(xy, 2, weighted.mean, w=w, ...)  # mean resultant vector
        rbar <- sqrt(sum(xy^2)) # mean resultant vector length
        if(return == "rbar") return(rbar)
        iso <- sqrt(1-rbar) # circular standard deviation
        if(return == "iso") return(iso)
        angle <- bearing(c(0,0), xy)
        if(return == "angle") return(angle)
}

cw2d <- function(data, colors = c("black", "yellow", "green", 
                                  "cyan", "blue", "magenta", "red"), 
                 origin = NULL, xyratio = NULL, kernel = NULL){
        require(colormap)
        #data = select(x, est5c, est6c)
        result <- rep(NA, nrow(data))
        a <- which(!is.na(apply(data, 1, sum)))
        data <- na.omit(data)
        if(is.null(origin)){
                origin <- c(sum(range(data[, 1], na.rm = T))/2, sum(range(data[, 
                                                                               2], na.rm = T))/2)}
        xrange <- range(data[, 1])
        yrange <- range(data[, 2])
        xmag <- plyr::round_any(max(abs(xrange)), (xrange[2] - xrange[1])/20, 
                                ceiling)
        ymag <- plyr::round_any(max(abs(yrange)), (yrange[2] - yrange[1])/20, 
                                ceiling)
        if(is.null(xyratio)) xyratio <- xmag/ymag
        pdata <- as.data.frame(polarize(data, xyratio = xyratio, 
                                        xorigin = origin[1], yorigin = origin[2]))
        names(pdata) <- c("distance", "angle")
        if(!is.null(kernel)) pdata$distance <- kernel(pdata$distance)
        pdata$angle <- pdata$angle/360
        n <- length(colors) - 1
        pdata$cl <- ceiling(pdata$angle * n) + 1
        pdata$fl <- floor(pdata$angle * n) + 1
        col <- matrix(NA, length(pdata$angle), 3)
        mx <- max(pdata$distance)
        colors <- col2rgb(colors)
        pal <- colors[, c(2:ncol(colors), 2)]/255
        center <- colors[, 1]/255
        center <- as.vector(center)
        getcol <- function(x) {
                interp <- x[2] * n - x[4] + 1
                col_angle <- (as.vector(pal[, x[3]]) * interp + as.vector(pal[, 
                                                                              x[4]]) * (1 - interp))
                col_angle * x[1]/mx + center * (1 - x[1]/mx)
        }
        col <- t(apply(pdata, 1, getcol))
        col[pdata$distance == 0, ] <- center
        result[a] <- rgb(col)
        return(result)
}

inv_logit <- function(x) exp(x)/(1+exp(x))

log10inc <- function(x){ # change in log10 precip
        y <- rnorm(1)
        z <- rnorm(1)
        ((y*10^(z+x)) - (y*10^z)) / (y*10^z)
}

inv_log10inc <- function(x){ # proportional change
        y <- runif(1)
        y2 = (x * y) + y
        log10(y2) - log10(y)
}

style <- theme_bw() + 
        theme(strip.text = element_text(color = "white"),
              strip.background = element_rect(fill = "black", color = "black"),
              legend.position = "bottom")


# build basis splines
tensor_splines <- function(x, y, xbounds = NULL, ybounds = NULL,
                         knots = 2, degree = 3){
        k = seq(0, 1, length.out = knots)
        k = k[2:(length(k)-1)]
        if(knots == 2) k <- NULL
        
        library(splines)
        if(is.null(xbounds)) xbounds <- range(x)
        if(is.null(ybounds)) xbounds <- range(y)
        B1 <- bs(x, knots = k, Boundary.knots = xbounds, degree=degree, intercept = TRUE)
        B2 <- bs(y, knots = k, Boundary.knots = xbounds, degree=degree, intercept = TRUE)
        B <- rep(1, nrow(B1))
        for(i in 1:ncol(B1)) for(j in 1:ncol(B2)) B <- cbind(B, B1[,i] * B2[,j])
        B <- B[,2:ncol(B)]
        B
}

# indices of adjacent spline bases
tensor_adj <- function(s, d = 1){
        b <- matrix(1:ncol(s), sqrt(ncol(s)))
        f <- function(x, d){
                n <- length(x)
                y <- lapply(1:(d+1), function(i) x[i:(n+i-1)])
                y <- do.call("paste", y)
                y[!grepl("NA", y)]
        }
        bx <- as.vector(apply(b, 1, f, d = d))
        by <- as.vector(apply(b, 2, f, d = d))
        tibble(x = c(bx, by)) %>% 
                separate(x, paste0("i", 1:(d+1)), convert = T) %>% 
                as.matrix()
}


define <- function(var, vars, define = T){
        def <- function(var){
                list(bio5 = "summer max temp", 
                     bio6 = "winter min temp", 
                     bio12 = "precipitation",
                     bio1 = "tmean")[[var]]
        }
        recode(var,
               "1" = def(vars[1]),
               "2" = def(vars[2]),
               "3" = def(vars[3]))
        y <- recode(var,
                    "1" = vars[1],
                    "2" = vars[2],
                    "3" = vars[3])
        if(define) y <- recode(var,
                               "1" = def(vars[1]),
                               "2" = def(vars[2]),
                               "3" = def(vars[3]))
        return(y)
}

label_vars <- function(x){
        data.frame(z = recode(as.character(x[,1]),
                              "windward" = "windward exposure",
                              "tpi" = "elevational position",
                              "tpis" = "prominence",
                              "bio12" = "moisture",
                              "bio5" = "high temperature",
                              "bio6" = "low temperature"))
}

label_vars2 <- function(x){
        y <- data.frame(z = recode(as.character(x[,1]),
                                   "windward" = "windward\nexposure",
                                   "tpi" = "elevational\nposition",
                                   "tpis" = "prominence",
                                   "bio12" = "moisture",
                                   "bio5" = "high temperature",
                                   "bio6" = "low temperature"))
        if(ncol(x)>1) y$z <- paste0(y$z, "\n", recode(as.character(x[,2]),
                                                      "predicted" = "[predicted]",
                                                      "observed" = "[observed]"))
        y
}

set_order <- function(x) x %>%
        mutate(clim_var = factor(clim_var, 
                                 levels = c("bio12", "bio5", "bio6", "MEAN")),
               topo_var = factor(topo_var, 
                                 levels = c("northness", "eastness", "windward", "tpi", "MEAN")))



GeomArrow2 <- ggplot2::ggproto("GeomArrow2", Geom,
                               required_aes = c("x", "y", "mag", "angle"),
                               default_aes = ggplot2::aes(color = "black", fill = "white", scale = 1),
                               draw_key = draw_key_polygon,
                               draw_panel = function(data, panel_scales, coord, scale = 1, 
                                                     theta = 15, length = .5, alpha = 1) {
                                       coords <- coord$transform(data, panel_scales)
                                       
                                       Mmag <- max(coords$mag)
                                       coords$mag <- with(coords, mag/Mmag*coords$scale)
                                       
                                       coords$dx <- with(coords, cos(angle)*mag)*scale
                                       coords$dy <- with(coords, sin(angle)*mag)*scale
                                       
                                       xx <- unit.c(unit(coords$x, "npc"), 
                                                    unit(coords$x, "npc") + unit(coords$dx, "snpc"))
                                       yy <- unit.c(unit(coords$y, "npc"), 
                                                    unit(coords$y, "npc") + unit(coords$dy, "snpc"))
                                       pol <- grid::polylineGrob(x = xx, y = yy,
                                                                 default.units = "npc",
                                                                 arrow = grid::arrow(angle = theta, type = "closed", 
                                                                                     length = unit(length, "lines")),
                                                                 gp = grid::gpar(col = coords$colour,
                                                                                 fill = coords$fill,
                                                                                 alpha = alpha),
                                                                 id = rep(seq(nrow(coords)), 2))
                                       pol
                                       
                               })

geom_arrow2 <- function(mapping = NULL, data = NULL, stat = "identity",
                        position = "identity", na.rm = FALSE, show.legend = NA,
                        inherit.aes = TRUE, scale = 1, ...) {
        layer(geom = GeomArrow2,
              mapping = mapping,
              data = data,
              stat = stat,
              position = position,
              show.legend = show.legend,
              inherit.aes = inherit.aes,
              params = list(na.rm = na.rm, scale = scale, ...)
        )
}

weighted.cor <- function(a, b, w) cov.wt(cbind(a, b), wt = w, cor = TRUE)$cor[1,2]
decile <- function(x) floor(rank(x) / (length(x)+1) * 10)
logit <- function(p) log(p / (1-p))
