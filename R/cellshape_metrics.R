



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



ashape_to_poly2 <- function(coord, alpha_var_range = 1:1000){


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



get_ashape_chull_poly <- function(x) {

    poly_res <- try(ashape_to_poly2(x), silent = TRUE)
    if (methods::is(poly_res, "try-error")) {
        poly_res <- list(NA, NA)
    }

    ch_poly_res <- try(chull_to_poly2(x), silent = TRUE)
    if (methods::is(ch_poly_res, "try-error")) {
        ch_poly_res <- NA
    }

    return(list(ashape_poly = poly_res[[1]],
                chull_poly = ch_poly_res,
                ashape_alpha = poly_res[[2]]))
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
        sp_poly <- ashape_to_poly2(x)[[1]]
        sf_poly <- sf::st_as_sf(sp_poly)
    }

    if (methods::is(x, "SpatialPolygons")) {
        sf_poly <- sf::st_as_sf(x)
    }

    major_axis <- find_major_axis(sf_poly)
    minor_axis <- find_minor_axis(sf_poly,major_axis$major_axis_endpoints)

    return(as.numeric(minor_axis$minor_axis_length/major_axis$major_axis_length))

}

sphericity <- function(x){

    if (methods::is(x, "matrix") | methods::is(x, "data.frame")) {
        # Calculate the alpha shape
        res <- ashape_to_poly2(x)[[1]]
    }

    if (methods::is(x, "SpatialPolygons")) {
        res <- x
    }

    bbox <- res@bbox
    pol <- spatstat.geom::owin(xrange = c(bbox["x", "min"], bbox["x", "max"]), yrange = c(bbox["y", "min"], bbox["y", "max"]))
    #calculate outer circle
    outer_radius <- spatstat.geom::boundingradius(pol)[[1]]
    #get inner circle
    inner_radius <- spatstat.geom::incircle(pol)

    #perform calculation of radius of inner circle divided by radius of outer circle
    return(inner_radius$r/outer_radius)
}

solidity <- function(x, chull_polygon = NULL) {
    if (methods::is(x, "matrix") | methods::is(x, "data.frame")) {
        poly_res <- ashape_to_poly2(x)[[1]]
        ch_polygon <- chull_to_poly2(x)
    }

    if (methods::is(x, "SpatialPolygons") & methods::is(chull_polygon, "SpatialPolygons")) {
        poly_res <- x
        ch_polygon <- chull_polygon
    }

    area_p <- as.numeric(raster::area(poly_res))
    ch_area <- as.numeric(raster::area(ch_polygon))

    return(area_p / ch_area)
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
    perimeter <- as.numeric(lwgeom::st_perimeter_2d(sf::st_as_sf(polygon)))

    # Calculate perimeter of convex hull object
    ch_perimeter <- as.numeric(lwgeom::st_perimeter_2d(sf::st_as_sf(ch_polygon)))

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
    area <- as.numeric(raster::area(polygon))

    # Calculate perimeter of convex hull object
    ch_perimeter <- as.numeric(lwgeom::st_perimeter_2d(sf::st_as_sf(ch_polygon)))

    #circularity calculation
    return((4*pi*area)/ ((ch_perimeter)^2))

}

compactness <- function(x) {

    if (methods::is(x, "matrix") | methods::is(x, "data.frame")) {
        # Calculate the alpha shape
        poly_res <- ashape_to_poly2(x)[[1]]
    }

    if (methods::is(x, "SpatialPolygons")) {
        poly_res <- x
    }

    # Calculate area of polygon object
    area <- as.numeric(raster::area(poly_res))

    # Calculate perimeter of hull object
    perimeter <- as.numeric(lwgeom::st_perimeter_2d(sf::st_as_sf(poly_res)))


    return((4*pi*area)/ ((perimeter)^2))
}

# support functions for eccentricity to find major and minor axis, this function essentially breaks a line into #num samples number of points to find the longest line through the polygon that is within the polygon
# the more samples the more accurate but also computational expensive, angle threshhold allows for slight leeway in minor not being perfectly perpendicular to major, again a computationally expensive trade off for better accuracy
find_major_axis <- function(poly_res, num_samples = 25) {

    #we get the boundary poly_res and then create lots of samples along its perimeter
    boundary <- sf::st_boundary(poly_res)
    boundary_length <- sf::st_length(boundary)[[1]]
    segment_length <- boundary_length / num_samples
    boundary_points <- sf::st_coordinates(sf::st_segmentize(boundary, segment_length))

    max_distance <- -Inf
    max_indices <- c(0, 0)
    num_boundary_points <- nrow(boundary_points)
    #loop through sampled x,y boundary points to find lines, keep biggest
    for (i in 1:(num_boundary_points - 1)) {
        for (j in (i + 1):num_boundary_points) {
            line_segment <- sf::st_linestring(matrix(c(boundary_points[i,], boundary_points[j,]), nrow = 2, byrow = TRUE))
            contains_result <- sf::st_contains(poly_res, sf::st_sfc(line_segment))

            if (length(contains_result[[1]]) > 0) {  # line segment is contained inside the poly_res or part of its perimeter
                current_distance <- sqrt((boundary_points[i, 1] - boundary_points[j, 1])^2 +
                                             (boundary_points[i, 2] - boundary_points[j, 2])^2)

                if (current_distance > max_distance) {
                    max_distance <- current_distance
                    max_indices <- c(i, j)
                }
            }
        }
    }

    major_axis_endpoints <- boundary_points[max_indices, ]
    major_axis_length <- max_distance
    #here we return the endpoints and length, ultimately eccentricity just is minor axis/major axis but this way we can do plots
    result <- list(
        major_axis_endpoints = major_axis_endpoints,
        major_axis_length = major_axis_length
    )

    return(result)
}
find_minor_axis <- function(poly_res, major_axis_endpoints, num_samples = 25, angle_threshold = 2) {
    boundary <- sf::st_boundary(poly_res)
    boundary_length <- sf::st_length(boundary)[[1]]
    segment_length <- boundary_length / num_samples
    boundary_points <- sf::st_coordinates(sf::st_segmentize(boundary, segment_length))

    major_axis_vector <- major_axis_endpoints[2,] - major_axis_endpoints[1,]
    major_axis_angle <- atan2(major_axis_vector[2], major_axis_vector[1])
    major_axis_line <- sf::st_linestring(major_axis_endpoints)

    max_distance <- -Inf
    max_indices <- c(0, 0)
    num_boundary_points <- nrow(boundary_points)

    for (i in 1:(num_boundary_points - 1)) {
        for (j in (i + 1):num_boundary_points) {
            line_segment <- sf::st_linestring(matrix(c(boundary_points[i,], boundary_points[j,]), nrow = 2, byrow = TRUE))
            contains_result <- sf::st_contains(poly_res, sf::st_sfc(line_segment))

            if (length(contains_result[[1]]) > 0) {  # line segment is contained inside the poly_res or part of its perimeter
                # Check if the line segment intersects the major axis
                intersects_result <- sf::st_intersects(major_axis_line, sf::st_sfc(line_segment))

                if (length(intersects_result[[1]]) > 0) {
                    line_segment_vector <- boundary_points[j,] - boundary_points[i,]
                    line_segment_angle <- atan2(line_segment_vector[2], line_segment_vector[1])

                    angle_difference <- abs(major_axis_angle - line_segment_angle) * (180 / pi)
                    angle_difference <- min(angle_difference, 360 - angle_difference)

                    if (angle_difference >= (90 - angle_threshold) && angle_difference <= (90 + angle_threshold)) {
                        current_distance <- sqrt((boundary_points[i, 1] - boundary_points[j, 1])^2 +
                                                     (boundary_points[i, 2] - boundary_points[j, 2])^2)

                        if (current_distance > max_distance) {
                            max_distance <- current_distance
                            max_indices <- c(i, j)
                        }
                    }
                }
            }
        }
    }

    minor_axis_endpoints <- boundary_points[max_indices, ]
    minor_axis_length <- max_distance

    result <- list(
        minor_axis_endpoints = minor_axis_endpoints,
        minor_axis_length = minor_axis_length
    )

    return(result)
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


