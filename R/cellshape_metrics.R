



#This function converts an alphashape output into a points path suitable to create a polygon in R
#THe base convex hull function doesnt require this
ashape2poly <- function(ashape){
    # Convert node numbers into characters
    ashape$edges[,1] <- as.character(ashape$edges[,1])
    ashape_graph <- igraph::graph_from_edgelist(ashape$edges[,1:2], directed = FALSE)
    if (!igraph::is.connected(ashape_graph)) {
        stop("Graph not connected")
    }
    if (any(igraph::degree(ashape_graph) != 2)) {
        stop("Graph not circular")
    }
    if (igraph::clusters(ashape_graph)$no > 1) {
        stop("Graph composed of more than one circle")
    }
    # Delete one edge to create a chain
    cut_graph <- ashape_graph - igraph::E(ashape_graph)[1]
    # Find chain end points
    ends <- names(which(igraph::degree(cut_graph) == 1))
    path <- igraph::get.shortest.paths(cut_graph, ends[1], ends[2])$vpath[[1]]
    # this is an index into the points
    pathX <- as.numeric(igraph::V(ashape_graph)[path]$name)
    # join the ends
    pathX <- c(pathX, pathX[1])
    return(pathX)
}



ashape_to_poly2 <- function(coord, alpha_var_range = 1:100){


    for (alpha_var in alpha_var_range){
        alpha_shape <- tryCatch(alphahull::ashape(coord, alpha = alpha_var),
                                error = function(e) e)

        if (!inherits(alpha_shape, "error")) {
            pth <- tryCatch(ashape2poly(alpha_shape),
                            error = function(e) e)

            if ( !inherits(pth, "error")){

                hull_points <-  coord[pth,]
                poly <- sp::Polygon(hull_points)
                poly_list <- list(poly)

                # Create spatial polygon object
                sp_poly <- sp::SpatialPolygons(list(sp::Polygons(poly_list,
                                                                 ID = 1)))


                return(list(sp_poly = sp_poly, alpha = alpha_var))
            }
        }

    }
}


#this function creates a convex hull, this is necessary for many metrics
chull_to_poly2 <- function(coord){



    hull <- grDevices::chull(coord)

    hull_points <- coord[hull, ]

    # Create polygon object
    poly <- sp::Polygon(hull_points)

    # Create list of polygon objects
    poly_list <- list(poly)

    # Create spatial polygon object
    sp_poly <- sp::SpatialPolygons(list(sp::Polygons(poly_list, ID = 1)))

    return(sp_poly)
}





# for functions requiring the bounding box we just need lines rather than sp_polygon,
# you can also get the bounding box from polygon but this is more computationally expensive
get_ashape_lines <- function(coord){
    coord <- as.matrix(coord)
    colnames(coord) <- c("x", "y")
    #this function creates a polygon from an alpha hull from points
    for (alpha_var in 1:100){
        alpha_shape <- alphahull::ashape(coord, alpha = alpha_var)
        pth <- tryCatch(ashape2poly(alpha_shape), error=function(e) e)

        if ( !inherits(pth, "error")){
            hull_points <-  coord[pth,]
            return(hull_points)
        }
    }
}

#
# ashape_to_poly <- function(x,y){
#     points <- cbind(x,y)
#
#     for (alpha_var in 1:100){
#         alpha_shape <- alphahull::ashape(points, alpha = alpha_var)
#         pth <- tryCatch(ashape2poly(alpha_shape),error=function(e) e)
#
#         if ( !inherits(pth, "error")){
#
#             hull_points <-  points[pth,]
#             poly <- sp::Polygon(hull_points)
#             poly_list <- list(poly)
#
#             # Create spatial polygon object
#             sp_poly <- sp::SpatialPolygons(list(sp::Polygons(poly_list, ID = 1)))
#
#
#             return(sp_poly)
#         }
#     }
# }
#

get_ashape_chull_poly <- function(x) {

    polygon <- try(ashape_to_poly2(x), silent = TRUE)
    if (methods::is(polygon, "try-error")) {
        polygon <- list(NA, NA)
    }

    ch_polygon <- try(chull_to_poly2(x), silent = TRUE)
    if (methods::is(ch_polygon, "try-error")) {
        ch_polygon <- NA
    }

    return(list(ashape_poly = polygon[[1]],
                chull_poly = ch_polygon,
                ashape_alpha = polygon[[2]]))
}

elongation <- function(x) {

    if (methods::is(x, "matrix") | methods::is(x, "data.frame")) {
        # Calculate the alpha shape
        lines <- get_ashape_lines(x)
    }

    if (methods::is(x, "SpatialPolygons")) {
        lines <- x@polygons[[1]]@Polygons[[1]]@coords
    }

    #nice function to extract bbox that is not a flat square
    bb <- shotGroups::getMinBBox(lines)

    #pretty sure that width here is opposite of what is on slides
    #http://www.cyto.purdue.edu/cdroms/micro2/content/education/wirth10.pdf page 25
    return(bb$height/bb$width)

}


eccentricity <- function(x){



    if (methods::is(x, "matrix") | methods::is(x, "data.frame")) {
        # Calculate the alpha shape
        lines <- get_ashape_lines(x)
    }

    if (methods::is(x, "SpatialPolygons")) {
        lines <- x@polygons[[1]]@Polygons[[1]]@coords
    }


    #nice function to extract bbox that is not a flat square
    bb <- shotGroups::getMinBBox(lines)


    #this is super similar to elongation but basically we take the max and min, for many cells we get the same result as elongation
    major_axis <- max(bb$height, bb$width)
    minor_axis <- min(bb$height, bb$width)

    return(minor_axis/major_axis)

}



sphericity <- function(x){

    if (methods::is(x, "matrix") | methods::is(x, "data.frame")) {
        # Calculate the alpha shape
        res <- ashape_to_poly2(x)[[1]]
    }

    if (methods::is(x, "SpatialPolygons")) {
        res <- x
    }


    #spatstat needs a different type of spatial object
    pol <- maptools::as.owin.SpatialPolygons(res)
    #if (is(pol, "try-error")) { return(NULL) }
    #calculate outer circle
    outer_radius <- spatstat.geom::boundingradius(pol)[[1]]
    #get inner circle
    inner_radius <- spatstat.geom::incircle(pol)
    #perform calculation of radius of inner circle divided by radius of outer circle
    return(inner_radius$r/outer_radius)
}


solidity <- function(x, chull_polygon = NULL){

    if (methods::is(x, "matrix") | methods::is(x, "data.frame")) {
        polygon <- ashape_to_poly2(x)[[1]]
        ch_polygon <- chull_to_poly2(x)
    }

    if (methods::is(x, "SpatialPolygons") & methods::is(chull_polygon, "SpatialPolygons")) {
        polygon <- x
        ch_polygon <- chull_polygon
    }


    area <- as.numeric(rgeos::gArea(polygon))

    ch_area <- as.numeric(rgeos::gArea(ch_polygon))

    return(area/ch_area)

}


convexity <- function(x, chull_polygon = NULL){


    if (methods::is(x, "matrix") | methods::is(x, "data.frame")) {
        polygon <- ashape_to_poly2(x)[[1]]
        ch_polygon <- chull_to_poly2(x)
    }

    if (methods::is(x, "SpatialPolygons") & methods::is(chull_polygon, "SpatialPolygons")) {
        polygon <- x
        ch_polygon <- chull_polygon
    }



    # Calculate area of polygon object
    perimeter <- as.numeric(rgeos::gLength(polygon))

    # Calculate perimeter of convex hull object
    ch_perimeter <- as.numeric(rgeos::gLength(ch_polygon))

    return(ch_perimeter/perimeter)

}


circularity <- function(x, chull_polygon = NULL){
    #This one kinda sucks for us because we dont really have a ground truth polygon, its basically just alpha convex hull(tighter polygon)
    # versus a normal convex hull


    if (methods::is(x, "matrix") | methods::is(x, "data.frame")) {
        polygon <- ashape_to_poly2(x)[[1]]
        ch_polygon <- chull_to_poly2(x)
    }

    if (methods::is(x, "SpatialPolygons") & methods::is(chull_polygon, "SpatialPolygons")) {
        polygon <- x
        ch_polygon <- chull_polygon
    }


    # Calculate area of polygon object
    area <- as.numeric(rgeos::gArea(polygon))

    # Calculate perimeter of convex hull object
    convex_perimeter <- as.numeric(rgeos::gLength(ch_polygon))

    #circularity calculation
    return((4*pi*area)/ ((convex_perimeter)^2))

}



compactness <- function(x) {

    if (methods::is(x, "matrix") | methods::is(x, "data.frame")) {
        # Calculate the alpha shape
        polygon <- ashape_to_poly2(x)[[1]]
    }

    if (methods::is(x, "SpatialPolygons")) {
        polygon <- x
    }

    # Calculate area of polygon object
    area <- as.numeric(rgeos::gArea(polygon))

    # Calculate perimeter of convex hull object
    perimeter <- as.numeric(rgeos::gLength(polygon))

    #Compactness calculation
    return((4*pi*area)/ ((perimeter)^2))
}




cell_shape_fn <- function(x, method) {

    if (is.null(x$ashape_poly)) {
        return(NA)
    }

    if (method == "elongation") {
        return(elongation(x$ashape_poly))
    }

    if (method == "compactness") {
        return(compactness(x$ashape_poly))
    }

    if (method == "eccentricity") {
        return(eccentricity(x$ashape_poly))
    }

    if (method == "sphericity") {
        return(sphericity(x$ashape_poly))
    }

    if (method == "convexity") {
        if (is.null(x$chull_poly)) {
            return(NA)
        }
        return(convexity(x$ashape_poly, x$chull_poly))
    }

    if (method == "solidity") {
        if (is.null(x$chull_poly)) {
            return(NA)
        }
        return(solidity(x$ashape_poly, x$chull_poly))
    }

    if (method == "circularity") {
        if (is.null(x$chull_poly)) {
            return(NA)
        }
        return(circularity(x$ashape_poly, x$chull_poly))
    }

}

.run_cell_shape_metrics <- function(spe, poly_objects, method,
                                    use_BPPARAM = use_BPPARAM,
                                    verbose = verbose) {


    if (verbose) {
        print(paste("Calculating", method))
    }

    res <- BiocParallel::bplapply(poly_objects, function(x) {
        res_tmp <- try(cell_shape_fn(x, method), silent = TRUE)
        if (methods::is(res_tmp, "try-error")) {
            return(NA)
        } else {
            return(res_tmp)
        }
    }, BPPARAM = use_BPPARAM)
    res <- unlist(res)
    colData(spe)[, method] <- res[as.character(spe$cell_id)]
    spe <- .add_metrics(spe, method, "cell_level")
    return(spe)
}


