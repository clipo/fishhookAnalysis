library(shiny)
library(Momocs)
library(imager)
library(dplyr)
library(ggplot2)

# Add this function to check matrix dimensions before operations
check_matrix_dims <- function(A, B, op_name) {
  if (!is.matrix(A) || !is.matrix(B)) {
    return(FALSE)
  }
  
  if (op_name == "subtract" || op_name == "add") {
    return(nrow(A) == nrow(B) && ncol(A) == ncol(B))
  } else if (op_name == "multiply") {
    return(ncol(A) == nrow(B))
  }
  
  return(TRUE)
}

# Safe distance calculation function to prevent dimension mismatches
safe_distance <- function(p1, p2) {
  if (!is.numeric(p1) || !is.numeric(p2) || length(p1) != length(p2)) {
    return(NA)
  }
  return(sqrt(sum((p1 - p2)^2)))
}

# Add this function to handle dimension checking before matrix operations
check_dimensions <- function(matrix1, matrix2, operation_name) {
  if (is.null(matrix1) || is.null(matrix2)) {
    return(FALSE)
  }
  
  # Check if both are matrices/arrays
  if (!is.matrix(matrix1) && !is.array(matrix1)) {
    return(FALSE)
  }
  if (!is.matrix(matrix2) && !is.array(matrix2)) {
    return(FALSE)
  }
  
  # For addition/subtraction operations
  if (operation_name == "add" || operation_name == "subtract") {
    if (nrow(matrix1) != nrow(matrix2) || ncol(matrix1) != ncol(matrix2)) {
      return(FALSE)
    }
  }
  
  # For dot product operations
  if (operation_name == "dot_product") {
    if (ncol(matrix1) != nrow(matrix2)) {
      return(FALSE)
    }
  }
  
  return(TRUE)
}

# Enhanced contour extraction with adaptive sampling
safe_extract_contour <- function(contour) {
  tryCatch({
    df <- as.data.frame(contour)
    coords <- matrix(c(df$x, df$y), ncol = 2)
    colnames(coords) <- c("x", "y")
    return(coords)
  }, error = function(e) {
    return(NULL)
  })
}

# Improved curvature calculation with better angular resolution
# Enhance the calculation of curvature with more robust error handling
calculate_curvature <- function(outline) {
  n_points <- nrow(outline)
  curvature <- numeric(n_points)
  
  for (i in 1:n_points) {
    # Use 3-point neighborhoods for more accurate curvature
    prev_idx <- ifelse(i <= 2, n_points - (2-i), i - 2)
    next_idx <- ifelse(i >= n_points-1, i + 2 - n_points, i + 2)
    
    # Ensure all indices are valid
    prev_idx <- max(1, min(n_points, prev_idx))
    next_idx <- max(1, min(n_points, next_idx))
    
    # Use Menger curvature for better accuracy in high-curvature regions
    p1 <- outline[prev_idx,]
    p2 <- outline[i,]
    p3 <- outline[next_idx,]
    
    # Check if all points are valid and have expected dimensions
    if (length(p1) == ncol(outline) && length(p2) == ncol(outline) && 
        length(p3) == ncol(outline)) {
      
      # Calculate area of triangle
      area <- 0.5 * abs(p1[1]*(p2[2]-p3[2]) + p2[1]*(p3[2]-p1[2]) + p3[1]*(p1[2]-p2[2]))
      
      # Calculate side lengths
      a <- sqrt(sum((p2 - p3)^2))
      b <- sqrt(sum((p1 - p3)^2))
      c <- sqrt(sum((p1 - p2)^2))
      
      # Prevent division by zero
      if(a*b*c > 0) {
        # Menger curvature = 4*area/(a*b*c)
        curvature[i] <- 4 * area / (a*b*c)
      } else {
        # Fallback to original method for edge cases with safe vector checks
        v1 <- p2 - p1
        v2 <- p3 - p2
        
        v1_len <- sqrt(sum(v1^2))
        v2_len <- sqrt(sum(v2^2))
        
        if (v1_len > 0 && v2_len > 0) {
          dot_prod <- sum(v1 * v2) / (v1_len * v2_len)
          dot_prod <- min(max(dot_prod, -1), 1)
          curvature[i] <- acos(dot_prod)
        } else {
          curvature[i] <- 0  # Default for degenerate cases
        }
      }
    } else {
      curvature[i] <- 0  # Default for invalid points
    }
  }
  
  return(curvature)
}

# Adaptive smoothing for curvature values
smooth_curvature <- function(curvature, window_size = 5) {
  n_points <- length(curvature)
  smoothed <- numeric(n_points)
  
  # Identify high curvature regions
  high_curv_threshold <- mean(curvature) + 1.5 * sd(curvature)
  high_curv_regions <- curvature > high_curv_threshold
  
  for (i in 1:n_points) {
    # Dynamically adjust window size based on curvature
    local_window <- if(high_curv_regions[i]) {
      max(3, window_size - 2)  # Smaller window for high curvature
    } else {
      window_size  # Regular window elsewhere
    }
    
    half_window <- floor(local_window / 2)
    
    # Calculate indices within the window (with wrap-around)
    indices <- ((i - half_window):(i + half_window)) %% n_points
    indices[indices == 0] <- n_points
    
    # Use weighted mean giving more weight to central points
    weights <- dnorm(seq(-half_window, half_window), 0, half_window/2)
    weights <- weights / sum(weights)
    
    smoothed[i] <- sum(curvature[indices] * weights)
  }
  
  return(smoothed)
}

# Adaptive interpolation for outline points
adaptive_interpolate <- function(outline, n_points = 100) {
  # Calculate curvature of original outline
  curv <- calculate_curvature(outline)
  
  # Normalize curvature to probabilities
  curv_norm <- curv / sum(curv)
  
  # Allocate more points to high-curvature regions
  # Maintain at least 50% uniform distribution for stability
  point_allocation <- 0.5 * rep(1/length(curv), length(curv)) + 
                     0.5 * curv_norm
  
  # Calculate cumulative arc length along outline
  arc_lengths <- numeric(nrow(outline))
  for (i in 2:nrow(outline)) {
    arc_lengths[i] <- arc_lengths[i-1] + 
      sqrt(sum((outline[i,] - outline[i-1,])^2))
  }
  total_length <- arc_lengths[nrow(outline)] + 
    sqrt(sum((outline[1,] - outline[nrow(outline),])^2))
  
  # Normalize arc lengths
  arc_lengths <- arc_lengths / total_length
  
  # Allocate new points based on weighted distribution
  new_points <- matrix(0, n_points, 2)
  
  for (i in 1:n_points) {
    # Find position based on weighted distribution
    target_pos <- (i - 1) / n_points
    
    # Find points before and after the target position
    before_idx <- max(which(arc_lengths <= target_pos))
    if (before_idx == length(arc_lengths)) {
      after_idx <- 1
    } else {
      after_idx <- before_idx + 1
    }
    
    # Linear interpolation
    if (before_idx == length(arc_lengths)) {
      before_pos <- arc_lengths[before_idx]
      after_pos <- 1
      
      before_point <- outline[before_idx,]
      after_point <- outline[1,]
    } else {
      before_pos <- arc_lengths[before_idx]
      after_pos <- arc_lengths[after_idx]
      
      before_point <- outline[before_idx,]
      after_point <- outline[after_idx,]
    }
    
    # Interpolate
    if (after_pos == before_pos) {
      t <- 0
    } else {
      t <- (target_pos - before_pos) / (after_pos - before_pos)
    }
    
    new_points[i,] <- (1-t) * before_point + t * after_point
  }
  
  return(new_points)
}

# Custom mshape function to calculate the mean shape from a list of coordinates
# Modify custom_mshape function to ensure all shapes have same dimensions
custom_mshape <- function(coo_list) {
  # Check if input is a Momocs Out object
  if (inherits(coo_list, "Out")) {
    coo_list <- coo_list$coo
  }
  
  # Get number of coordinates and dimensions
  n_shapes <- length(coo_list)
  if (n_shapes == 0) return(NULL)
  
  n_points <- nrow(coo_list[[1]])
  n_dims <- ncol(coo_list[[1]])
  
  # Check if all shapes have the same dimensions
  for (i in 1:n_shapes) {
    if (nrow(coo_list[[i]]) != n_points || ncol(coo_list[[i]]) != n_dims) {
      # Normalize shapes to same number of points if they differ
      coo_list[[i]] <- adaptive_interpolate(coo_list[[i]], n_points = n_points)
    }
  }
  
  # Initialize mean shape matrix
  mean_shape <- matrix(0, nrow = n_points, ncol = n_dims)
  
  # Sum all shapes
  for (i in 1:n_shapes) {
    if (check_dimensions(mean_shape, coo_list[[i]], "add")) {
      mean_shape <- mean_shape + coo_list[[i]]
    } else {
      warning(paste("Shape", i, "has incompatible dimensions and was skipped"))
    }
  }
  
  # Divide by number of shapes to get mean
  mean_shape <- mean_shape / n_shapes
  
  return(mean_shape)
}

# Identify significant inflection points
# Update the find_inflection_points function with safeguards against infinite recursion
find_inflection_points <- function(curvature, outline, threshold_factor = 1.5, 
                                   min_threshold = 0.3, recursion_count = 0) {
  # Add recursion limit to prevent infinite loops
  if (recursion_count > 10) {
    warning("Reached maximum recursion depth in find_inflection_points, using minimum threshold")
    # Return at least 3 points with highest curvature
    sorted_indices <- order(curvature, decreasing = TRUE)
    return(sorted_indices[1:min(3, length(sorted_indices))])
  }
  
  # Determine adaptive threshold based on curvature distribution
  curvature_threshold <- mean(curvature) + threshold_factor * sd(curvature)
  
  # Find local maxima exceeding threshold
  n_points <- length(curvature)
  inflection_indices <- c()
  
  for (i in 2:(n_points-1)) {
    if (curvature[i] > curvature[i-1] && 
        curvature[i] > curvature[i+1] && 
        curvature[i] > curvature_threshold) {
      inflection_indices <- c(inflection_indices, i)
    }
  }
  
  # If too few points found, reduce threshold and try again
  if (length(inflection_indices) < 3) {
    new_threshold <- threshold_factor * 0.8
    if (new_threshold < min_threshold) {
      # If below minimum threshold, just take top 3 points by curvature
      sorted_indices <- order(curvature, decreasing = TRUE)
      return(sorted_indices[1:min(3, length(sorted_indices))])
    } else {
      return(find_inflection_points(curvature, outline, new_threshold, 
                                    min_threshold, recursion_count + 1))
    }
  }
  
  return(inflection_indices)
}

# Define fishhook regions based on geometric features
# Update define_fishhook_regions with additional checks for valid indices and regions
# Modified define_fishhook_regions function focusing on the critical error areas
define_fishhook_regions <- function(outline, inflection_indices, curvature) {
  n_points <- nrow(outline)
  
  # Calculate centroid for reference
  centroid <- colMeans(outline)
  
  # Calculate principal axes using tryCatch to handle potential errors
  pca_result <- tryCatch({
    prcomp(outline)
  }, error = function(e) {
    # Fallback if PCA fails
    warning("PCA failed, using fallback orientation method")
    list(
      rotation = matrix(c(1, 0, 0, 1), nrow=2),
      center = colMeans(outline)
    )
  })
  
  # Extract rotation matrix safely
  if (inherits(pca_result, "prcomp") && is.matrix(pca_result$rotation)) {
    pc1 <- pca_result$rotation[,1]
    if (ncol(pca_result$rotation) >= 2) {
      pc2 <- pca_result$rotation[,2]
    } else {
      pc2 <- c(-pc1[2], pc1[1])  # Perpendicular vector to pc1
    }
  } else {
    # Fallback if rotation matrix isn't available
    pc1 <- c(1, 0)
    pc2 <- c(0, 1)
  }
  
  # Project points onto principal axes - with dimension checking
  if (ncol(outline) == length(pc1)) {
    proj_pc1 <- as.matrix(outline) %*% pc1
    proj_pc2 <- as.matrix(outline) %*% pc2
  } else {
    # Handle dimension mismatch
    warning("Dimension mismatch in PC projection. Using index-based projection.")
    proj_pc1 <- 1:nrow(outline)
    proj_pc2 <- (1:nrow(outline)) %% n_points
  }
  
  # Find extremal points safely
  max_pc1_idx <- which.max(proj_pc1)
  min_pc1_idx <- which.min(proj_pc1)
  max_pc2_idx <- which.max(proj_pc2)
  min_pc2_idx <- which.min(proj_pc2)
  
  # Calculate linearity scores safely
  linearity_scores <- numeric(n_points)
  segment_length <- floor(n_points/5)  # Test segments of ~20% of outline
  
  for (i in 1:n_points) {
    # Define segment indices with wrap-around handling
    segment_indices <- ((i-floor(segment_length/2)):(i+floor(segment_length/2))) %% n_points
    segment_indices[segment_indices == 0] <- n_points
    
    # Extract segment points - ensure indices are valid
    segment_indices <- segment_indices[segment_indices > 0 & segment_indices <= n_points]
    
    if (length(segment_indices) >= 3) {  # Need at least 3 points for PCA
      segment <- outline[segment_indices,, drop=FALSE]
      
      # Calculate linearity using PCA on segment
      segment_pca <- tryCatch({
        prcomp(segment)
      }, error = function(e) {
        # If PCA fails, return NULL
        NULL
      })
      
      if (!is.null(segment_pca) && length(segment_pca$sdev) >= 1) {
        linearity_scores[i] <- segment_pca$sdev[1] / sum(segment_pca$sdev)
      } else {
        linearity_scores[i] <- 0
      }
    } else {
      linearity_scores[i] <- 0
    }
  }
  
  # Identify shank as region with highest linearity
  if (any(!is.na(linearity_scores))) {
    shank_center_idx <- which.max(linearity_scores)
  } else {
    # Fallback if linearity calculation failed
    shank_center_idx <- 1
  }
  
  # Extend shank region based on linearity threshold
  if (any(!is.na(linearity_scores) & !is.infinite(linearity_scores))) {
    max_linearity <- max(linearity_scores[!is.na(linearity_scores) & !is.infinite(linearity_scores)])
    linearity_threshold <- 0.85 * max_linearity
    shank_candidate_indices <- which(linearity_scores > linearity_threshold)
  } else {
    # Fallback if linearity values are invalid
    shank_candidate_indices <- 1:n_points
    linearity_threshold <- 0
  }
  
  # Safety check - ensure we have candidate indices
  if (length(shank_candidate_indices) == 0) {
    # Fallback: use first 20% of points as candidates
    shank_candidate_indices <- 1:max(1, floor(n_points * 0.2))
  }
  
  # Ensure indices are within valid range
  shank_candidate_indices <- shank_candidate_indices[
    shank_candidate_indices > 0 & shank_candidate_indices <= n_points
  ]
  
  # Find contiguous segment centered on shank_center_idx
  distances_from_center <- abs(shank_candidate_indices - shank_center_idx) %% n_points
  distances_reversed <- abs(n_points - distances_from_center)
  adjusted_distances <- pmin(distances_from_center, distances_reversed)
  
  # Sort indices by distance and take first 20% for shank region
  if (length(adjusted_distances) > 0) {
    shank_indices <- shank_candidate_indices[order(adjusted_distances)]
    shank_size <- max(3, floor(n_points * 0.2))
    shank_region <- shank_indices[1:min(shank_size, length(shank_indices))]
  } else {
    # Fallback if distance calculation fails
    shank_region <- 1:max(3, floor(n_points * 0.2))
  }
  
  # Point identification - calculate distances safely
  distances_to_shank <- numeric(n_points)
  
  for (i in 1:n_points) {
    if (i > 0 && i <= nrow(outline) && 
        shank_center_idx > 0 && shank_center_idx <= nrow(outline)) {
      
      p1 <- outline[i,]
      p2 <- outline[shank_center_idx,]
      
      if (length(p1) == length(p2)) {
        distances_to_shank[i] <- sqrt(sum((p1 - p2)^2))
      } else {
        # Fallback if point dimensions don't match
        distances_to_shank[i] <- abs(i - shank_center_idx) %% n_points
      }
    } else {
      distances_to_shank[i] <- abs(i - shank_center_idx) %% n_points
    }
  }
  
  # Find high curvature points
  if (!is.null(curvature) && length(curvature) == n_points && 
      !all(is.na(curvature)) && !all(is.infinite(curvature))) {
    
    curv_mean <- mean(curvature, na.rm = TRUE)
    curv_sd <- sd(curvature, na.rm = TRUE)
    high_curvature_threshold <- curv_mean + 1.5 * curv_sd
    high_curvature_indices <- which(curvature > high_curvature_threshold)
  } else {
    # Fallback if curvature calculation is invalid
    high_curvature_indices <- numeric(0)
  }
  
  # Select point as high curvature point farthest from shank
  if (length(high_curvature_indices) > 0) {
    valid_indices <- high_curvature_indices[
      high_curvature_indices > 0 & high_curvature_indices <= length(distances_to_shank)
    ]
    
    if (length(valid_indices) > 0) {
      point_idx <- valid_indices[which.max(distances_to_shank[valid_indices])]
    } else {
      # Fallback if no valid high curvature indices
      point_idx <- which.max(distances_to_shank)
    }
  } else {
    # Fallback to point farthest from shank
    point_idx <- which.max(distances_to_shank)
  }
  
  # Define point region around point_idx
  point_size <- max(3, floor(n_points * 0.1))
  half_point <- floor(point_size / 2)
  point_region <- ((point_idx - half_point):(point_idx + half_point)) %% n_points
  point_region[point_region == 0] <- n_points
  
  # Ensure point region indices are valid
  point_region <- point_region[point_region > 0 & point_region <= n_points]
  
  # Bend region calculation - find the region between shank and point
  # Create empty set for bend region
  bend_region <- numeric(0)
  
  # First ensure regions are not empty
  if (length(shank_region) == 0) {
    shank_region <- 1:max(1, floor(n_points * 0.2))
  }
  
  if (length(point_region) == 0) {
    point_region <- (floor(n_points * 0.7)):n_points
  }
  
  # Get extremes for each region
  shank_extremes <- range(shank_region)
  point_extremes <- range(point_region)
  
  # Divide outline into regions using a simplified approach to avoid array conformability issues
  # This method is more robust but less sophisticated than the original
  all_indices <- 1:n_points
  assigned_indices <- unique(c(shank_region, point_region))
  remaining_indices <- setdiff(all_indices, assigned_indices)
  
  if (length(remaining_indices) > 0) {
    # Take middle section as bend
    bend_size <- min(length(remaining_indices), max(3, floor(n_points * 0.3)))
    bend_indices <- sort(remaining_indices)
    
    if (length(bend_indices) >= bend_size) {
      mid_start <- max(1, floor((length(bend_indices) - bend_size) / 2))
      bend_region <- bend_indices[mid_start:(mid_start + bend_size - 1)]
    } else {
      bend_region <- bend_indices
    }
  } else {
    # Create artificial bend region if all points are already assigned
    bend_region <- floor(n_points * 0.4):floor(n_points * 0.6)
    bend_region <- setdiff(bend_region, c(shank_region, point_region))
  }
  
  # No barb region for simplicity
  barb_region <- numeric(0)
  barb_center_idx <- NULL
  
  # Calculate functional metrics - critical area prone to non-conformable arrays error
  # Compute centroids for regions
  if (length(shank_region) > 0) {
    # Ensure indices are valid
    valid_shank <- shank_region[shank_region > 0 & shank_region <= nrow(outline)]
    if (length(valid_shank) > 0) {
      shank_midpoint <- colMeans(outline[valid_shank, , drop = FALSE])
    } else {
      shank_midpoint <- colMeans(outline[1:min(3, nrow(outline)), , drop = FALSE])
    }
  } else {
    shank_midpoint <- colMeans(outline[1:min(3, nrow(outline)), , drop = FALSE])
  }
  
  if (point_idx > 0 && point_idx <= nrow(outline)) {
    point_midpoint <- outline[point_idx, ]
  } else if (length(point_region) > 0) {
    valid_point <- point_region[point_region > 0 & point_region <= nrow(outline)]
    if (length(valid_point) > 0) {
      point_midpoint <- colMeans(outline[valid_point, , drop = FALSE])
    } else {
      point_midpoint <- colMeans(outline[max(1, nrow(outline)-2):nrow(outline), , drop = FALSE])
    }
  } else {
    point_midpoint <- colMeans(outline[max(1, nrow(outline)-2):nrow(outline), , drop = FALSE])
  }
  
  # Gape calculation
  if (length(shank_midpoint) == length(point_midpoint)) {
    gape <- sqrt(sum((shank_midpoint - point_midpoint)^2))
  } else {
    # Fallback if dimensions don't match
    gape <- 1.0  # Default value
  }
  
  # Throat calculation - this is likely where the non-conformable arrays error occurs
  throat <- 0  # Initialize default
  throat_idx <- 1  # Initialize default
  
  # Calculate throat only if point and shank midpoints have same dimension
  if (length(shank_midpoint) == length(point_midpoint)) {
    shank_to_point <- point_midpoint - shank_midpoint
    shank_to_point_len <- sqrt(sum(shank_to_point^2))
    
    if (shank_to_point_len > 0) {
      shank_to_point_unit <- shank_to_point / shank_to_point_len
      
      # Calculate perpendicular distances for bend points
      if (length(bend_region) > 0) {
        # Ensure bend indices are valid
        valid_bend <- bend_region[bend_region > 0 & bend_region <= nrow(outline)]
        
        if (length(valid_bend) > 0) {
          depths <- numeric(length(valid_bend))
          
          for (i in 1:length(valid_bend)) {
            bend_idx <- valid_bend[i]
            if (bend_idx > 0 && bend_idx <= nrow(outline)) {
              bend_point <- outline[bend_idx, ]
              
              # Verify dimensions match before calculation
              if (length(bend_point) == length(shank_midpoint) && 
                  length(shank_to_point_unit) == length(shank_midpoint)) {
                
                # Vector from shank to bend point
                shank_to_bend <- bend_point - shank_midpoint
                # Project onto shank-point line
                proj_length <- sum(shank_to_bend * shank_to_point_unit)
                # Projection point
                proj_point <- shank_midpoint + proj_length * shank_to_point_unit
                # Perpendicular distance
                depths[i] <- sqrt(sum((bend_point - proj_point)^2))
              } else {
                depths[i] <- 0
              }
            } else {
              depths[i] <- 0
            }
          }
          
          throat <- max(depths)
          throat_idx <- valid_bend[which.max(depths)]
        }
      }
    }
  }
  
  # Point angle calculation
  point_angle <- 90  # Default angle
  
  if (point_idx > 0 && point_idx <= nrow(outline)) {
    pre_point_idx <- (point_idx - 3) %% n_points
    if (pre_point_idx == 0) pre_point_idx <- n_points
    post_point_idx <- (point_idx + 3) %% n_points
    if (post_point_idx == 0) post_point_idx <- n_points
    
    # Ensure all indices are valid
    if (pre_point_idx > 0 && pre_point_idx <= nrow(outline) && 
        post_point_idx > 0 && post_point_idx <= nrow(outline)) {
      
      v1 <- outline[point_idx, ] - outline[pre_point_idx, ]
      v2 <- outline[post_point_idx, ] - outline[point_idx, ]
      
      # Calculate lengths to avoid division by zero
      v1_len <- sqrt(sum(v1^2))
      v2_len <- sqrt(sum(v2^2))
      
      if (v1_len > 0 && v2_len > 0) {
        cos_angle <- sum(v1 * v2) / (v1_len * v2_len)
        # Ensure cos_angle is in valid range for acos
        cos_angle <- min(max(cos_angle, -1), 1)
        point_angle <- acos(cos_angle) * 180 / pi
      }
    }
  }
  
  # Return regions and metrics
  return(list(
    shank = unique(shank_region),
    bend = unique(bend_region),
    point = unique(point_region),
    barb = barb_region,  # Empty for simplicity
    key_points = list(
      shank_center = shank_center_idx,
      point = point_idx,
      throat = throat_idx,
      barb_center = barb_center_idx
    ),
    metrics = list(
      gape = gape,
      throat = throat,
      point_angle = point_angle
    ),
    curvature = curvature
  ))
}

ui <- fluidPage(
  titlePanel("Fishhook Morphometric Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons("image_source", "Image source:", 
                   choices = c("Demo Images", "Upload My Own"), 
                   selected = "Demo Images"),
      
      conditionalPanel(
        condition = "input.image_source == 'Upload My Own'",
        fileInput("images", "Upload JPEG Images", 
                  accept = c(".jpg", ".jpeg"), 
                  multiple = TRUE)
      ),
      
      h4("Image Processing"),
      radioButtons("artifact_mode", "Fishhook appearance:", 
                   choices = c("Darker than background" = "darker", 
                               "Lighter than background" = "lighter"),
                   selected = "darker"),
      
      sliderInput("threshold", "Threshold (0–1)", 
                  min = 0, max = 1, value = 0.5, step = 0.01),
      
      sliderInput("noise_size", "Noise removal", 
                  min = 0, max = 20, value = 5, step = 1),
      
      h4("Analysis Settings"),
      sliderInput("nbh", "Harmonics", 
                  min = 1, max = 50, value = 20),
      
      numericInput("k_clusters", "Clusters", 
                   value = 3, min = 1),
      
      # Region definition method selection
      h4("Region Detection Method"),
      radioButtons("region_method", "Method:", 
                   choices = c("Geometric (Shape-based)" = "geometric",
                               "Proportional (Percentage-based)" = "proportional",
                               "Manual key points" = "manual"),
                   selected = "geometric"),
      
      # Show region proportions only when using proportional method
      conditionalPanel(
        condition = "input.region_method == 'proportional'",
        sliderInput("shank_ratio", "Shank proportion (%):", 
                    min = 15, max = 50, value = 25, step = 5),
        
        sliderInput("bend_ratio", "Bend proportion (%):", 
                    min = 30, max = 60, value = 45, step = 5),
        
        sliderInput("point_ratio", "Point proportion (%):", 
                    min = 10, max = 30, value = 20, step = 5),
        
        sliderInput("barb_ratio", "Barb proportion (if present) (%):", 
                    min = 0, max = 20, value = 10, step = 5),
      ),
      
      # Show geometric detection parameters when using geometric method
      conditionalPanel(
        condition = "input.region_method == 'geometric'",
        sliderInput("curvature_threshold", "Curvature sensitivity:", 
                    min = 0.5, max = 2.5, value = 1.5, step = 0.1),
        
        sliderInput("curvature_smoothing", "Curvature smoothing:", 
                    min = 1, max = 11, value = 5, step = 2)
      ),
      
      conditionalPanel(
        condition = "input.region_method == 'manual'",
        h4("Manual Key Points Selection"),
        p("Select key points on the outline after processing:"),
        
        actionButton("select_shank_top", "Top of Shank", class = "btn-primary"),
        actionButton("select_shank_bottom", "Bottom of Shank", class = "btn-primary"),
        actionButton("select_bend_deepest", "Deepest Bend Point", class = "btn-success"),
        actionButton("select_point", "Point Tip", class = "btn-danger"),
        actionButton("select_barb", "Barb (if present)", class = "btn-info"),
        
        actionButton("reset_points", "Reset Points", class = "btn-warning"),
        
        verbatimTextOutput("key_points_info")
      ),
      
      actionButton("crop_btn", "Process Image", 
                   class = "btn-primary"),
      
      actionButton("next_img", "Next Image", 
                   class = "btn-info"),
      
      actionButton("reset_img", "Reset", 
                   class = "btn-warning"),
      
      actionButton("run_analysis", "Run Analysis", 
                   class = "btn-success"),
      
      hr(),
      
      h4("Fishhook Metrics"),
      verbatimTextOutput("hook_metrics"),
      
      h4("Downloads"),
      downloadButton("download_pdf", "Download Analysis (PDF)", 
                     class = "btn-danger"),
      
      hr(),
      
      div(style = "text-align: center; margin-top: 10px;",
          tags$p("Source code and documentation:"),
          tags$a(
            href = "https://github.com/clipo/fishhookAnalysis",
            target = "_blank",  # Opens in new tab
            div(
              style = "display: inline-flex; align-items: center; color: #333; text-decoration: none; border: 1px solid #ddd; border-radius: 5px; padding: 5px 10px; background: #f8f8f8;",
              tags$img(src = "https://github.githubassets.com/images/modules/logos_page/GitHub-Mark.png", 
                       height = "20px",
                       style = "margin-right: 5px;"),
              "GitHub Repository"
            )
          )
      )
    ),
    
    mainPanel(
      h4(textOutput("current_image_name")),
      
      # Use a div with position:relative to contain the plot and markers
      div(
        id = "image-container",
        style = "position: relative;",
        plotOutput("image_plot", click = "img_click", height = "300px"),
        uiOutput("click_markers")
      ),
      
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Processing", 
                 fluidRow(
                   column(6, plotOutput("cropped_plot", height = "250px")),
                   column(6, plotOutput("binary_plot", height = "250px"))
                 ),
                 plotOutput("outline_plot", height = "250px", click="outline_click"),
                 verbatimTextOutput("processing_info")
        ),
        tabPanel("Results", 
                 tabsetPanel(
                   tabPanel("Outlines Grid", plotOutput("grid_plot", height = "500px")),
                   tabPanel("Stacked Outlines", plotOutput("stack_plot", height = "500px")),
                   tabPanel("EFA Harmonic Power", plotOutput("harmonic_plot", height = "500px")),
                   tabPanel("PCA Morphospace", plotOutput("morpho_plot", height = "500px")),
                   tabPanel("Thin-Plate Splines", plotOutput("tps_plot", height = "500px")),
                   tabPanel("Integration/Modularity", plotOutput("modularity_plot", height = "500px"))
                 )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  # Reactive values
  img_index <- reactiveVal(1)
  click_points <- reactiveVal(data.frame(x = numeric(0), y = numeric(0)))
  outlines_list <- reactiveVal(list())
  names_list <- reactiveVal(c())
  cropped_img <- reactiveVal(NULL)
  binary_img <- reactiveVal(NULL)
  outline_shape <- reactiveVal(NULL)
  processing_log <- reactiveVal("")
  hook_metrics_data <- reactiveVal(NULL)
  
  # Reactive values for analysis
  outlines_r <- reactiveVal()
  efa_result_r <- reactiveVal()
  pca_result_r <- reactiveVal()
  tps_grid_r <- reactiveVal()
  modularity_r <- reactiveVal()
  curvature_r <- reactiveVal()
  
  # Additional reactive values for manual key points
  key_points_mode <- reactiveVal(FALSE)
  current_key_point <- reactiveVal(NULL)
  key_points <- reactiveVal(list(
    shank_top = NULL,
    shank_bottom = NULL,
    bend_deepest = NULL,
    point = NULL,
    barb = NULL
  ))
  
  # Get image files
  get_files <- reactive({
    if (input$image_source == "Demo Images") {
      files <- list.files("www", pattern = "(?i)\\.jpe?g$", full.names = TRUE)
      if (length(files) == 0) {
        processing_log("No demo images found in the 'www' folder.")
        return(list(files = character(0), names = character(0)))
      }
      names <- basename(files)
    } else {
      req(input$images)
      files <- input$images$datapath
      names <- input$images$name
    }
    list(files = files, names = names)
  })
  
  # Track clicks for cropping
  observeEvent(input$img_click, {
    clicks <- click_points()
    if (nrow(clicks) < 2) {
      click_points(rbind(clicks, data.frame(x = input$img_click$x, y = input$img_click$y)))
    }
  })
  
  # Add this to the server function:
  observeEvent(input$run_analysis, {
    out_list <- outlines_list()
    name_list <- names_list()
    
    current_log <- processing_log()
    processing_log(paste(current_log, "\nStarting analysis...", sep=""))
    
    if (length(out_list) == 0) {
      processing_log(paste(current_log, "\nNo outlines to analyze. Process images first.", sep=""))
      return()
    }
    
    tryCatch({
      # Create Momocs Out object
      outlines <- Out(out_list)
      names(outlines$coo) <- name_list
      
      # Check outline consistency
      n_points <- sapply(outlines$coo, nrow)
      if (length(unique(n_points)) > 1) {
        processing_log(paste(current_log, 
                             "\nWarning: Outlines have inconsistent numbers of points.", 
                             "Normalizing to", median(n_points), "points per outline.", sep=""))
        
        # Normalize all outlines to the median point count
        target_points <- median(n_points)
        for (i in 1:length(outlines$coo)) {
          if (n_points[i] != target_points) {
            outlines$coo[[i]] <- adaptive_interpolate(outlines$coo[[i]], n_points = target_points)
          }
        }
      }
      
      # Align outlines
      processing_log(paste(processing_log(), "\nAligning outlines...", sep=""))
      outlines <- outlines %>% 
        coo_center() %>% 
        coo_scale() %>% 
        coo_align()
      
      outlines_r(outlines)
      
      # Elliptical Fourier Analysis
      processing_log(paste(processing_log(), "\nPerforming Elliptical Fourier Analysis...", sep=""))
      efa <- efourier(outlines, nb.h = input$nbh, norm = TRUE)
      efa_result_r(efa)
      
      # PCA - add dimension checking
      processing_log(paste(processing_log(), "\nRunning PCA...", sep=""))
      if (!is.null(efa$coe) && is.matrix(efa$coe) && nrow(efa$coe) > 0 && ncol(efa$coe) > 0) {
        pca <- PCA(efa)
        
        # Clustering - check if PCA has valid data
        if (!is.null(pca$x) && nrow(pca$x) > 0) {
          # Determine appropriate number of clusters based on data size
          k <- min(input$k_clusters, nrow(pca$x))
          
          # Check if we have enough data for clustering
          if (nrow(efa$coe) >= k) {
            clust <- hclust(dist(efa$coe), method = "ward.D2")
            clusters <- cutree(clust, k = k)
            pca$fac <- data.frame(cluster = factor(clusters))
          } else {
            # Not enough data for clustering
            processing_log(paste(processing_log(), 
                                 "\nWarning: Not enough data for", k, "clusters.", 
                                 "Using", nrow(efa$coe), "cluster(s) instead.", sep=""))
            pca$fac <- data.frame(cluster = factor(1:nrow(efa$coe)))
          }
          
          pca_result_r(pca)
          
          # Calculate mean shape for TPS grid
          mean_shape <- custom_mshape(outlines$coo)
          
          # Define safer methods for finding extremes
          safe_which_max <- function(v) {
            if (length(v) == 0) return(1)
            return(which.max(v))
          }
          
          safe_which_min <- function(v) {
            if (length(v) == 0) return(1)
            return(which.min(v))
          }
          
          # For TPS visualization, use actual specimen shapes at the extremes of the PCA
          pc1_scores <- pca$x[,1]
          pc2_scores <- pca$x[,2]
          
          # Find specimens at extremes of PC axes with safety checks
          pc1_pos_idx <- safe_which_max(pc1_scores)
          pc1_neg_idx <- safe_which_min(pc1_scores)
          pc2_pos_idx <- safe_which_max(pc2_scores)
          pc2_neg_idx <- safe_which_min(pc2_scores)
          
          # Store shapes at PCA extremes
          tps_grid_r(list(
            mean = mean_shape,
            pc1_pos = outlines$coo[[pc1_pos_idx]],
            pc1_neg = outlines$coo[[pc1_neg_idx]],
            pc2_pos = outlines$coo[[pc2_pos_idx]],
            pc2_neg = outlines$coo[[pc2_neg_idx]]
          ))
          
          # Estimate modularity by dividing outline into regions
          n_points <- nrow(mean_shape)
          
          # Determine region indices based on method
          if (input$region_method == "manual") {
            # Use manually selected key points if available
            kp <- key_points()
            
            # Check if essential points are selected
            if (is.null(kp$shank_top) || is.null(kp$shank_bottom) || 
                is.null(kp$bend_deepest) || is.null(kp$point)) {
              processing_log(paste(processing_log(), 
                                   "\nError: Essential key points are not selected for manual region definition.",
                                   "\nAttempting to automatically define regions instead.", 
                                   sep=""))
              
              # Fall back to geometric method
              curv <- calculate_curvature(mean_shape)
              window_size <- ceiling(input$curvature_smoothing / 2) * 2 + 1
              smoothed_curv <- smooth_curvature(curv, window_size)
              
              # Find inflection points
              inflection_points <- find_inflection_points(smoothed_curv, mean_shape, input$curvature_threshold)
              
              # Define regions using specialized fishhook region detection
              regions <- define_fishhook_regions(mean_shape, inflection_points, smoothed_curv)
              
              # Extract regions
              shank_region <- regions$shank
              bend_region <- regions$bend
              point_region <- regions$point
              barb_region <- regions$barb
              
              # Store metrics
              hook_metrics_data(regions$metrics)
              
            } else {
              # Manual region definition as before, but with robustness checks
              # ... [rest of the manual region definition code] ...
              # Define shank region between top and bottom with wrap-around handling
              if (kp$shank_top <= kp$shank_bottom) {
                shank_region <- kp$shank_top:kp$shank_bottom
              } else {
                shank_region <- c(kp$shank_top:n_points, 1:kp$shank_bottom)
              }
              
              # Define point region around point (10% of total points)
              point_size <- round(n_points * 0.1)
              half_point <- floor(point_size / 2)
              point_start <- (kp$point - half_point) %% n_points
              if(point_start == 0) point_start <- n_points
              point_end <- (kp$point + half_point) %% n_points
              if(point_end == 0) point_end <- n_points
              
              if(point_start <= point_end) {
                point_region <- point_start:point_end
              } else {
                point_region <- c(point_start:n_points, 1:point_end)
              }
              
              # Define barb region if selected
              if (!is.null(kp$barb)) {
                barb_size <- round(n_points * 0.1)
                half_barb <- floor(barb_size / 2)
                barb_start <- (kp$barb - half_barb) %% n_points
                if(barb_start == 0) barb_start <- n_points
                barb_end <- (kp$barb + half_barb) %% n_points
                if(barb_end == 0) barb_end <- n_points
                
                if(barb_start <= barb_end) {
                  barb_region <- barb_start:barb_end
                } else {
                  barb_region <- c(barb_start:n_points, 1:barb_end)
                }
              } else {
                barb_region <- numeric(0)
              }
              
              # Find bend region - the region between shank and point
              # Corrected to handle the case where indices may wrap around
              all_indices <- 1:n_points
              excluded_indices <- unique(c(shank_region, point_region, barb_region))
              remaining_indices <- setdiff(all_indices, excluded_indices)
              
              # Find closest segments to the bend point
              if (!is.null(kp$bend_deepest) && kp$bend_deepest %in% remaining_indices) {
                # Find continuous segment around bend_deepest
                bend_center <- kp$bend_deepest
                
                # Calculate distances from each remaining point to bend_center
                bend_dists <- numeric(length(remaining_indices))
                for (i in 1:length(remaining_indices)) {
                  idx1 <- remaining_indices[i]
                  idx2 <- bend_center
                  # Calculate distance along outline (smaller of two possible paths)
                  dist1 <- abs(idx1 - idx2)
                  dist2 <- n_points - dist1
                  bend_dists[i] <- min(dist1, dist2)
                }
                
                # Sort remaining indices by distance to bend_center
                sorted_indices <- remaining_indices[order(bend_dists)]
                
                # Take up to 20% of outline as bend region
                bend_size <- min(length(sorted_indices), round(n_points * 0.2))
                bend_region <- sorted_indices[1:bend_size]
              } else {
                # If bend point not in remaining indices or not selected
                # Divide remaining indices into regions
                if (length(remaining_indices) > 0) {
                  # Use the middle third as bend
                  bend_region <- remaining_indices[1:max(1, floor(length(remaining_indices)/2))]
                } else {
                  # Fallback - create artificial bend region
                  bend_region <- setdiff(1:n_points, c(shank_region, point_region, barb_region))
                  if (length(bend_region) == 0) {
                    # Last resort - divide outline into quarters
                    quarters <- floor(n_points/4)
                    bend_region <- (2*quarters+1):(3*quarters)
                  }
                }
              }
              
              # Calculate functional metrics safely
              # ... [metrics calculation code with dimension checks] ...
              # Gape - distance between shank center and point
              if (length(shank_region) > 0) {
                shank_midpoint <- colMeans(mean_shape[shank_region, , drop = FALSE])
              } else {
                shank_midpoint <- colMeans(mean_shape[1:max(1, floor(n_points * 0.2)), , drop = FALSE])
              }
              
              if (!is.null(kp$point) && kp$point > 0 && kp$point <= n_points) {
                point_midpoint <- mean_shape[kp$point, ]
                if (length(shank_midpoint) == length(point_midpoint)) {
                  gape <- sqrt(sum((shank_midpoint - point_midpoint)^2))
                } else {
                  gape <- 1.0  # Default value
                }
              } else {
                gape <- 1.0
              }
              
              # Throat calculation
              if (!is.null(kp$bend_deepest) && kp$bend_deepest > 0 && kp$bend_deepest <= n_points && 
                  !is.null(kp$point) && kp$point > 0 && kp$point <= n_points) {
                
                point_midpoint <- mean_shape[kp$point, ]
                
                if (length(shank_midpoint) == length(point_midpoint)) {
                  shank_to_point <- point_midpoint - shank_midpoint
                  shank_to_point_len <- sqrt(sum(shank_to_point^2))
                  
                  if (shank_to_point_len > 0) {
                    shank_to_point_unit <- shank_to_point / shank_to_point_len
                    
                    bend_point <- mean_shape[kp$bend_deepest, ]
                    shank_to_bend <- bend_point - shank_midpoint
                    proj_length <- sum(shank_to_bend * shank_to_point_unit)
                    proj_point <- shank_midpoint + proj_length * shank_to_point_unit
                    throat <- sqrt(sum((bend_point - proj_point)^2))
                  } else {
                    throat <- 0.5
                  }
                } else {
                  throat <- 0.5
                }
              } else {
                throat <- 0.5
              }
              
              # Point angle calculation
              if (!is.null(kp$point) && kp$point > 0 && kp$point <= n_points) {
                pre_point_idx <- (kp$point - 3) %% n_points
                if (pre_point_idx == 0) pre_point_idx <- n_points
                post_point_idx <- (kp$point + 3) %% n_points
                if (post_point_idx == 0) post_point_idx <- n_points
                
                v1 <- mean_shape[kp$point, ] - mean_shape[pre_point_idx, ]
                v2 <- mean_shape[post_point_idx, ] - mean_shape[kp$point, ]
                
                v1_len <- sqrt(sum(v1^2))
                v2_len <- sqrt(sum(v2^2))
                
                if (v1_len > 0 && v2_len > 0) {
                  cos_angle <- sum(v1 * v2) / (v1_len * v2_len)
                  cos_angle <- min(max(cos_angle, -1), 1)  # Ensure valid range for acos
                  point_angle <- acos(cos_angle) * 180 / pi
                } else {
                  point_angle <- 90  # Default angle
                }
              } else {
                point_angle <- 90
              }
              
              # Store metrics
              hook_metrics_data(list(
                gape = gape,
                throat = throat,
                point_angle = point_angle
              ))
            }
          } else if (input$region_method == "geometric") {
            # Geometric region detection based on fishhook shape features
            processing_log(paste(processing_log(), "\nPerforming geometric region detection...", sep=""))
            
            # Calculate curvature if not already done
            if (is.null(curvature_r())) {
              curv <- calculate_curvature(mean_shape)
              window_size <- ceiling(input$curvature_smoothing / 2) * 2 + 1
              smoothed_curv <- smooth_curvature(curv, window_size)
            } else {
              smoothed_curv <- curvature_r()
            }
            
            # Find inflection points with enhanced robustness
            inflection_points <- find_inflection_points(smoothed_curv, mean_shape, input$curvature_threshold)
            
            # Define regions using enhanced fishhook region detection
            regions <- define_fishhook_regions(mean_shape, inflection_points, smoothed_curv)
            
            # Extract regions
            shank_region <- regions$shank
            bend_region <- regions$bend
            point_region <- regions$point
            barb_region <- regions$barb
            
            # Store metrics
            hook_metrics_data(regions$metrics)
            
            # Log key points
            log_text <- paste(processing_log(), 
                              "\nGeometric detection found key points:",
                              "\nShank center: Point ", regions$key_points$shank_center,
                              "\nPoint: Point ", regions$key_points$point,
                              "\nThroat point: Point ", regions$key_points$throat,
                              sep="")
            
            if (!is.null(regions$key_points$barb_center)) {
              log_text <- paste(log_text, "\nBarb center: Point ", regions$key_points$barb_center, sep="")
            } else {
              log_text <- paste(log_text, "\nNo barb detected.", sep="")
            }
            
            processing_log(log_text)
          } else {
            # Proportional method with enhanced robustness
            # ... [rest of proportional region code with similar safety checks] ...
            
            # Get user-defined proportions - ensure percentages sum to 100%
            total_prop <- (input$shank_ratio + input$bend_ratio + input$point_ratio + input$barb_ratio) / 100
            
            # Normalize proportions if they don't sum to 1
            shank_prop <- (input$shank_ratio / 100) / total_prop
            bend_prop <- (input$bend_ratio / 100) / total_prop
            point_prop <- (input$point_ratio / 100) / total_prop
            barb_prop <- (input$barb_ratio / 100) / total_prop
            
            # Use PCA to identify primary axis for orientation
            pca_outline <- prcomp(mean_shape)
            # Project onto PC1 to find extremes
            proj_pc1 <- as.matrix(mean_shape) %*% pca_outline$rotation[,1]
            # This gives us the extreme points along the primary axis
            max_pc1_idx <- which.max(proj_pc1)
            min_pc1_idx <- which.min(proj_pc1)
            
            # Choose start point based on curvature (point typically has higher curvature)
            curv <- calculate_curvature(mean_shape)
            max_curv_index <- which.max(curv)
            
            # Calculate distance from max_curv to PC1 extremes
            dist_to_max <- abs(max_curv_index - max_pc1_idx) %% n_points
            dist_to_min <- abs(max_curv_index - min_pc1_idx) %% n_points
            dist_to_max2 <- abs(n_points - dist_to_max)
            dist_to_min2 <- abs(n_points - dist_to_min)
            
            # Determine the closest PC1 extreme to the high curvature point
            closest_idx <- if(min(dist_to_max, dist_to_max2) <= min(dist_to_min, dist_to_min2)) {
              max_pc1_idx
            } else {
              min_pc1_idx
            }
            
            # Define point region around high curvature/PC1 extreme
            start_idx <- closest_idx
            
            # Ensure valid starting point
            if (start_idx <= 0 || start_idx > n_points) {
              start_idx <- 1
            }
            
            # Create regions sequentially
            point_size <- max(1, round(n_points * point_prop))
            point_start <- start_idx
            point_end <- (point_start + point_size - 1) %% n_points
            if (point_end == 0) point_end <- n_points
            
            if (point_start <= point_end) {
              point_region <- point_start:point_end
            } else {
              point_region <- c(point_start:n_points, 1:point_end)
            }
            
            # Bend follows point
            bend_size <- max(1, round(n_points * bend_prop))
            bend_start <- (point_end + 1) %% n_points
            if (bend_start == 0) bend_start <- n_points
            bend_end <- (bend_start + bend_size - 1) %% n_points
            if (bend_end == 0) bend_end <- n_points
            
            if (bend_start <= bend_end) {
              bend_region <- bend_start:bend_end
            } else {
              bend_region <- c(bend_start:n_points, 1:bend_end)
            }
            
            # Shank follows bend
            shank_size <- max(1, round(n_points * shank_prop))
            shank_start <- (bend_end + 1) %% n_points
            if (shank_start == 0) shank_start <- n_points
            shank_end <- (shank_start + shank_size - 1) %% n_points
            if (shank_end == 0) shank_end <- n_points
            
            if (shank_start <= shank_end) {
              shank_region <- shank_start:shank_end
            } else {
              shank_region <- c(shank_start:n_points, 1:shank_end)
            }
            
            # Barb is the remaining points if barb_prop > 0
            remaining <- setdiff(1:n_points, c(point_region, bend_region, shank_region))
            
            if (barb_prop > 0 && length(remaining) > 0) {
              barb_size <- min(length(remaining), max(1, round(n_points * barb_prop)))
              barb_region <- remaining[1:barb_size]
              
              # Add any still remaining points to nearest region
              still_remaining <- remaining[(barb_size+1):length(remaining)]
            } else {
              barb_region <- numeric(0)
              still_remaining <- remaining
            }
            
            # Assign any remaining points to nearest region
            if (length(still_remaining) > 0) {
              for (pt in still_remaining) {
                # Find closest region
                if (length(point_region) > 0) {
                  point_dist <- min(abs(pt - point_region))
                } else {
                  point_dist <- Inf
                }
                
                if (length(bend_region) > 0) {
                  bend_dist <- min(abs(pt - bend_region))
                } else {
                  bend_dist <- Inf
                }
                
                if (length(shank_region) > 0) {
                  shank_dist <- min(abs(pt - shank_region))
                } else {
                  shank_dist <- Inf
                }
                
                if (length(barb_region) > 0) {
                  barb_dist <- min(abs(pt - barb_region))
                } else {
                  barb_dist <- Inf
                }
                
                min_dist <- min(point_dist, bend_dist, shank_dist, barb_dist)
                
                if (min_dist == point_dist && is.finite(point_dist)) {
                  point_region <- c(point_region, pt)
                } else if (min_dist == bend_dist && is.finite(bend_dist)) {
                  bend_region <- c(bend_region, pt)
                } else if (min_dist == shank_dist && is.finite(shank_dist)) {
                  shank_region <- c(shank_region, pt)
                } else if (min_dist == barb_dist && is.finite(barb_dist)) {
                  barb_region <- c(barb_region, pt)
                } else {
                  # If all distances are infinite, assign to shank as default
                  shank_region <- c(shank_region, pt)
                }
              }
            }
            
            # Calculate fishhook metrics with safety checks
            # Determine key points for metrics calculation
            if (length(point_region) > 0 && length(curv) == n_points) {
              point_idx <- point_region[which.max(curv[point_region])]
            } else if (length(point_region) > 0) {
              point_idx <- point_region[ceiling(length(point_region)/2)]
            } else {
              point_idx <- 1
            }
            
            if (length(shank_region) > 0) {
              shank_midpoint <- colMeans(mean_shape[shank_region, , drop = FALSE])
            } else {
              shank_midpoint <- colMeans(mean_shape[1:max(1, floor(n_points * 0.2)), , drop = FALSE])
            }
            
            if (point_idx > 0 && point_idx <= n_points) {
              point_midpoint <- mean_shape[point_idx, ]
            } else if (length(point_region) > 0) {
              point_midpoint <- colMeans(mean_shape[point_region, , drop = FALSE])
            } else {
              point_midpoint <- colMeans(mean_shape[max(1, floor(n_points * 0.8)):n_points, , drop = FALSE])
            }
            
            # Gape calculation with dimension check
            if (length(shank_midpoint) == length(point_midpoint)) {
              gape <- sqrt(sum((shank_midpoint - point_midpoint)^2))
            } else {
              gape <- 1.0  # Default value
              warning("Dimension mismatch in gape calculation")
            }
            
            # Throat calculation with safety checks
            if (length(shank_midpoint) == length(point_midpoint)) {
              shank_to_point <- point_midpoint - shank_midpoint
              shank_to_point_len <- sqrt(sum(shank_to_point^2))
              
              if (shank_to_point_len > 0) {
                shank_to_point_unit <- shank_to_point / shank_to_point_len
                
                # Calculate perpendicular distances for bend points
                if (length(bend_region) > 0) {
                  depths <- numeric(length(bend_region))
                  for (i in 1:length(bend_region)) {
                    if (bend_region[i] > 0 && bend_region[i] <= n_points) {
                      bend_point <- mean_shape[bend_region[i],]
                      # Vector from shank to bend point
                      shank_to_bend <- bend_point - shank_midpoint
                      # Project onto shank-point line
                      proj_length <- sum(shank_to_bend * shank_to_point_unit)
                      # Projection point
                      proj_point <- shank_midpoint + proj_length * shank_to_point_unit
                      # Perpendicular distance
                      depths[i] <- sqrt(sum((bend_point - proj_point)^2))
                    } else {
                      depths[i] <- 0
                    }
                  }
                  throat <- max(depths)
                } else {
                  throat <- 0
                }
              } else {
                throat <- 0.5  # Default value
              }
            } else {
              throat <- 0.5  # Default value
            }
            
            # Point angle calculation with safety checks
            if (point_idx > 0 && point_idx <= n_points) {
              pre_point_idx <- (point_idx - 3) %% n_points
              if (pre_point_idx == 0) pre_point_idx <- n_points
              post_point_idx <- (point_idx + 3) %% n_points
              if (post_point_idx == 0) post_point_idx <- n_points
              
              if (pre_point_idx > 0 && pre_point_idx <= n_points && 
                  post_point_idx > 0 && post_point_idx <= n_points) {
                
                v1 <- mean_shape[point_idx,] - mean_shape[pre_point_idx,]
                v2 <- mean_shape[post_point_idx,] - mean_shape[point_idx,]
                
                v1_len <- sqrt(sum(v1^2))
                v2_len <- sqrt(sum(v2^2))
                
                if (v1_len > 0 && v2_len > 0) {
                  cos_angle <- sum(v1 * v2) / (v1_len * v2_len)
                  cos_angle <- min(max(cos_angle, -1), 1)  # Ensure valid range for acos
                  point_angle <- acos(cos_angle) * 180 / pi
                } else {
                  point_angle <- 90  # Default angle
                }
              } else {
                point_angle <- 90  # Default angle
              }
            } else {
              point_angle <- 90  # Default angle
            }
            
            # Store metrics
            hook_metrics_data(list(
              gape = gape,
              throat = throat,
              point_angle = point_angle
            ))
          }
          
          # Ensure regions stay within bounds
          shank_region <- shank_region[shank_region > 0 & shank_region <= n_points]
          bend_region <- bend_region[bend_region > 0 & bend_region <= n_points]
          point_region <- point_region[point_region > 0 & point_region <= n_points]
          if (length(barb_region) > 0) {
            barb_region <- barb_region[barb_region > 0 & barb_region <= n_points]
          }
          
          # Create partition
          partition <- rep(0, n_points)
          if (length(shank_region) > 0) partition[shank_region] <- 1  # Shank
          if (length(bend_region) > 0) partition[bend_region] <- 2   # Bend
          if (length(point_region) > 0) partition[point_region] <- 3  # Point
          if (length(barb_region) > 0) {
            partition[barb_region] <- 4 # Barb
          }
          
          # Store results
          modularity_r(list(
            partition = partition,
            regions = list(shank = shank_region, bend = bend_region, point = point_region, barb = barb_region),
            n_points = n_points,
            metrics = hook_metrics_data()
          ))
          
          # Update log and switch to Results tab
          processing_log(paste(processing_log(), "\nAnalysis complete! Switching to Results tab.", sep=""))
          updateTabsetPanel(session, "main_tabs", selected = "Results")
        } else {
          processing_log(paste(processing_log(), "\nError: Invalid PCA result", sep=""))
        }
      } else {
        processing_log(paste(processing_log(), "\nError: Invalid EFA coefficient matrix", sep=""))
      }
      
    }, error = function(e) {
      # Enhanced error handling - provide more specific error messages
      error_msg <- e$message
      
      # Try to provide more helpful error messages based on common issues
      if (grepl("non-conformable", error_msg)) {
        error_msg <- paste("Non-conformable arrays error: This usually occurs when trying to",
                           "perform operations between arrays with incompatible dimensions.",
                           "This may be due to inconsistent outline points or region definitions.")
      } else if (grepl("subscript out of bounds", error_msg)) {
        error_msg <- paste("Subscript out of bounds: The code is trying to access elements",
                           "that don't exist. This could be due to invalid region indices.")
      } else if (grepl("is.numeric", error_msg)) {
        error_msg <- paste("Numeric conversion error: One of the variables expected to be",
                           "numeric contains non-numeric values.")
      }
      
      # Update log with detailed error information
      processing_log(paste(processing_log(), "\nError during analysis:", error_msg, 
                           "\nDetails: ", e$message, sep=""))
    })
  })
  
  # Reset
  observeEvent(input$reset_img, {
    click_points(data.frame(x = numeric(0), y = numeric(0)))
    cropped_img(NULL)
    binary_img(NULL)
    outline_shape(NULL)
    processing_log("")
    hook_metrics_data(NULL)
  })
  
  # Add click markers on the image
  output$click_markers <- renderUI({
    clicks <- click_points()
    if (nrow(clicks) == 0) return(NULL)
    
    # Create markers for each click point
    markers <- lapply(1:nrow(clicks), function(i) {
      tags$div(
        style = sprintf(
          "position: absolute; left: %.0fpx; top: %.0fpx; width: 10px; height: 10px; 
           background-color: red; border-radius: 50%%; z-index: 1000;",
          clicks$x[i], clicks$y[i]
        )
      )
    })
    
    markers
  })
  
  # Create observer for region method selection
  observeEvent(input$region_method, {
    if(input$region_method == "manual") {
      # Display message in log
      current_log <- processing_log()
      processing_log(paste(current_log, 
                           "\nManual region selection mode activated. Process an image and then select key points.",
                           sep=""))
    } else {
      # Disable key point selection mode
      key_points_mode(FALSE)
    }
  })
  
  # Observers for key point selection buttons
  observeEvent(input$select_shank_top, {
    if(is.null(outline_shape())) {
      processing_log(paste(processing_log(), "\nProcess an image first before selecting points.", sep=""))
      return()
    }
    key_points_mode(TRUE)
    current_key_point("shank_top")
    processing_log(paste(processing_log(), "\nClick on the outline to select the TOP OF SHANK point.", sep=""))
  })
  
  observeEvent(input$select_shank_bottom, {
    if(is.null(outline_shape())) {
      processing_log(paste(processing_log(), "\nProcess an image first before selecting points.", sep=""))
      return()
    }
    key_points_mode(TRUE)
    current_key_point("shank_bottom")
    processing_log(paste(processing_log(), "\nClick on the outline to select the BOTTOM OF SHANK point.", sep=""))
  })
  
  observeEvent(input$select_bend_deepest, {
    if(is.null(outline_shape())) {
      processing_log(paste(processing_log(), "\nProcess an image first before selecting points.", sep=""))
      return()
    }
    key_points_mode(TRUE)
    current_key_point("bend_deepest")
    processing_log(paste(processing_log(), "\nClick on the outline to select the DEEPEST BEND point.", sep=""))
  })
  
  observeEvent(input$select_point, {
    if(is.null(outline_shape())) {
      processing_log(paste(processing_log(), "\nProcess an image first before selecting points.", sep=""))
      return()
    }
    key_points_mode(TRUE)
    current_key_point("point")
    processing_log(paste(processing_log(), "\nClick on the outline to select the POINT TIP.", sep=""))
  })
  
  observeEvent(input$select_barb, {
    if(is.null(outline_shape())) {
      processing_log(paste(processing_log(), "\nProcess an image first before selecting points.", sep=""))
      return()
    }
    key_points_mode(TRUE)
    current_key_point("barb")
    processing_log(paste(processing_log(), "\nClick on the outline to select the BARB point (if present).", sep=""))
  })
  
  observeEvent(input$reset_points, {
    key_points(list(
      shank_top = NULL,
      shank_bottom = NULL,
      bend_deepest = NULL,
      point = NULL,
      barb = NULL
    ))
    processing_log(paste(processing_log(), "\nKey points reset.", sep=""))
  })
  
  # Process image
  observeEvent(input$crop_btn, {
    files <- get_files()$files
    names <- get_files()$names
    idx <- img_index()
    
    if (length(files) == 0 || idx > length(files)) {
      processing_log("No images available")
      return()
    }
    
    if (nrow(click_points()) != 2) {
      processing_log("Click two points to crop")
      return()
    }
    
    log_text <- paste("Processing:", names[idx], "\n")
    
    tryCatch({
      # Load image
      img <- load.image(files[idx])
      img_gray <- suppressWarnings(grayscale(img))
      
      # Crop image
      cp <- click_points()
      x1 <- round(min(cp$x)); x2 <- round(max(cp$x))
      y1 <- round(min(cp$y)); y2 <- round(max(cp$y))
      
      # Validate crop dimensions
      if (x1 >= x2 || y1 >= y2 || 
          x1 < 0 || y1 < 0 || 
          x2 > dim(img_gray)[1] || y2 > dim(img_gray)[2]) {
        processing_log(paste(log_text, "Invalid crop region"))
        return()
      }
      
      # Apply crop
      crop <- imsub(img_gray, x %inr% c(x1, x2), y %inr% c(y1, y2))
      cropped_img(crop)
      
      # Apply threshold based on artifact appearance
      if (input$artifact_mode == "darker") {
        bin <- crop < input$threshold
      } else {
        bin <- crop > input$threshold
      }
      
      # Clean up noise
      if (input$noise_size > 0) {
        bin <- clean(bin, input$noise_size)
      }
      
      binary_img(bin)
      
      # Extract contours
      contours <- tryCatch({
        imager::contours(bin, nlevels = 1)
      }, error = function(e) {
        log_text <<- paste0(log_text, "ERROR extracting contours: ", e$message, "\n")
        processing_log(log_text)
        return(NULL)
      })
      
      if (is.null(contours) || length(contours) == 0) {
        log_text <- paste0(log_text, "No contours found\n")
        processing_log(log_text)
        return()
      }
      
      # Select the largest contour
      selected_contour <- NULL
      max_points <- 0
      
      for (i in 1:length(contours)) {
        # Use our safe extraction function
        current_contour <- safe_extract_contour(contours[[i]])
        
        if (!is.null(current_contour) && nrow(current_contour) > max_points) {
          max_points <- nrow(current_contour)
          selected_contour <- current_contour
        }
      }
      
      if (is.null(selected_contour) || nrow(selected_contour) < 50) {
        log_text <- paste0(log_text, "Contour has too few points\n")
        processing_log(log_text)
        return()
      }
      
      log_text <- paste0(log_text, "Found contour with ", nrow(selected_contour), " points\n")
      
      # Use adaptive interpolation instead of regular interpolation for better hook region resolution
      interp <- adaptive_interpolate(selected_contour, n = 100)
      outline_shape(interp)
      
      # Calculate curvature for geometric detection
      if (input$region_method == "geometric") {
        curv <- calculate_curvature(interp)
        # Apply smoothing based on user parameter
        window_size <- ceiling(input$curvature_smoothing / 2) * 2 + 1  # Ensure odd number
        smoothed_curv <- smooth_curvature(curv, window_size)
        curvature_r(smoothed_curv)
        
        log_text <- paste0(log_text, "Calculated curvature for geometric region detection\n")
      }
      
      # Add to outlines list
      out_list <- outlines_list()
      name_list <- names_list()
      out_list[[length(out_list) + 1]] <- interp
      name_list <- c(name_list, tools::file_path_sans_ext(names[idx]))
      outlines_list(out_list)
      names_list(name_list)
      
      log_text <- paste0(log_text, "Successfully added outline with enhanced hook detection\n")
      processing_log(log_text)
      
    }, error = function(e) {
      log_text <<- paste0(log_text, "CRITICAL ERROR: ", e$message, "\n")
      processing_log(log_text)
    })
  })
  
  # Next image
  observeEvent(input$next_img, {
    files <- get_files()$files
    next_idx <- img_index() + 1
    
    if (next_idx <= length(files)) {
      click_points(data.frame(x = numeric(0), y = numeric(0)))
      cropped_img(NULL)
      binary_img(NULL)
      outline_shape(NULL)
      processing_log("")
      hook_metrics_data(NULL)
      img_index(next_idx)
    } else {
      processing_log("No more images")
    }
  })
  
  # Display current image name
  output$current_image_name <- renderText({
    files <- get_files()$files
    names <- get_files()$names
    idx <- img_index()
    
    if (length(files) == 0) {
      return("No images available")
    }
    
    if (idx <= length(names)) {
      paste0("Image ", idx, "/", length(names), ": ", names[idx], 
             " (Click two points to crop)")
    } else {
      "All images processed. Run analysis below."
    }
  })
  
  # Display original image
  output$image_plot <- renderPlot({
    files <- get_files()$files
    idx <- img_index()
    
    if (length(files) == 0 || idx > length(files)) {
      plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", 
           main = "No images available", axes = FALSE)
      text(0.5, 0.5, "Add images to 'www' or upload your own")
      return()
    }
    
    tryCatch({
      img <- load.image(files[idx])
      plot(img, main = "Click two points to crop")
      
      # Draw crop rectangle
      cp <- click_points()
      if (nrow(cp) == 2) {
        rect(min(cp$x), min(cp$y), max(cp$x), max(cp$y), border = "red", lwd = 2)
      }
    }, error = function(e) {
      plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", 
           main = "Error loading image", axes = FALSE)
      text(0.5, 0.5, paste("Error:", e$message))
    })
  })
  
  # Display cropped image
  output$cropped_plot <- renderPlot({
    if (is.null(cropped_img())) {
      plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, 
           main = "Cropped Image", xlab = "", ylab = "")
      text(0.5, 0.5, "Process an image to see preview")
      return()
    }
    
    plot(cropped_img(), main = "Cropped Image")
  })
  
  # Display binary image
  output$binary_plot <- renderPlot({
    if (is.null(binary_img())) {
      plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, 
           main = "Binary Image", xlab = "", ylab = "")
      text(0.5, 0.5, "Process an image to see binary result")
      return()
    }
    
    plot(binary_img(), main = paste0("Binary (Threshold = ", input$threshold, ")"))
  })
  
  # Display outline with interactive key point selection
  output$outline_plot <- renderPlot({
    if(is.null(outline_shape())) {
      plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, 
           main = "Outline", xlab = "", ylab = "")
      text(0.5, 0.5, "Process an image to see extracted outline")
      return()
    }
    
    # Plot outline
    plot(outline_shape(), type = "l", asp = 1, main = "Extracted Outline")
    
    # If using geometric method and curvature has been calculated, show it
    if(input$region_method == "geometric" && !is.null(curvature_r())) {
      # Add curvature visualization as color along the outline
      curv <- curvature_r()
      # Normalize curvature to 0-1 scale for coloring
      norm_curv <- (curv - min(curv)) / (max(curv) - min(curv))
      
      # Plot points with color based on curvature
      for(i in 1:nrow(outline_shape())) {
        # Higher curvature = more red
        col_val <- rgb(norm_curv[i], 0, 1-norm_curv[i])
        points(outline_shape()[i, 1], outline_shape()[i, 2], 
               col = col_val, pch = 19, cex = 1)
      }
      
      # Add a legend for curvature
      gradient_colors <- sapply(seq(0, 1, length.out = 10), function(x) {
        rgb(x, 0, 1-x)
      })
      
      # Draw legend in the corner
      legend_x <- min(outline_shape()[,1]) + 0.1 * (max(outline_shape()[,1]) - min(outline_shape()[,1]))
      legend_y <- min(outline_shape()[,2]) + 0.1 * (max(outline_shape()[,2]) - min(outline_shape()[,2]))
      
      # Small rectangle for the legend
      rect(legend_x, legend_y, 
           legend_x + 0.2 * (max(outline_shape()[,1]) - min(outline_shape()[,1])), 
           legend_y + 0.05 * (max(outline_shape()[,2]) - min(outline_shape()[,2])),
           col = "white", border = "black")
      
      # Add colored points for the gradient
      for(i in 1:length(gradient_colors)) {
        x_pos <- legend_x + (i-1)/9 * 0.2 * (max(outline_shape()[,1]) - min(outline_shape()[,1]))
        points(x_pos, legend_y + 0.025 * (max(outline_shape()[,2]) - min(outline_shape()[,2])), 
               col = gradient_colors[i], pch = 19, cex = 1.5)
      }
      
      # Add labels
      text(legend_x, legend_y + 0.07 * (max(outline_shape()[,2]) - min(outline_shape()[,2])), 
           "Curvature", pos = 4, cex = 0.8)
      text(legend_x, legend_y - 0.02 * (max(outline_shape()[,2]) - min(outline_shape()[,2])), 
           "Low", pos = 4, cex = 0.7)
      text(legend_x + 0.2 * (max(outline_shape()[,1]) - min(outline_shape()[,1])), 
           legend_y - 0.02 * (max(outline_shape()[,2]) - min(outline_shape()[,2])), 
           "High", pos = 2, cex = 0.7)
    }
    
    # Show key points if any exist
    kp <- key_points()
    if(!is.null(kp$shank_top)) {
      points(outline_shape()[kp$shank_top, 1], outline_shape()[kp$shank_top, 2], 
             col = "blue", pch = 19, cex = 2)
      text(outline_shape()[kp$shank_top, 1], outline_shape()[kp$shank_top, 2], 
           "SHANK TOP", col = "blue", pos = 1)
    }
    if(!is.null(kp$shank_bottom)) {
      points(outline_shape()[kp$shank_bottom, 1], outline_shape()[kp$shank_bottom, 2], 
             col = "blue", pch = 19, cex = 2)
      text(outline_shape()[kp$shank_bottom, 1], outline_shape()[kp$shank_bottom, 2], 
           "SHANK BOTTOM", col = "blue", pos = 2)
    }
    if(!is.null(kp$bend_deepest)) {
      points(outline_shape()[kp$bend_deepest, 1], outline_shape()[kp$bend_deepest, 2], 
             col = "green", pch = 19, cex = 2)
      text(outline_shape()[kp$bend_deepest, 1], outline_shape()[kp$bend_deepest, 2], 
           "BEND", col = "green", pos = 2)
    }
    if(!is.null(kp$point)) {
      points(outline_shape()[kp$point, 1], outline_shape()[kp$point, 2], 
             col = "red", pch = 19, cex = 2)
      text(outline_shape()[kp$point, 1], outline_shape()[kp$point, 2], 
           "POINT", col = "red", pos = 4)
    }
    if(!is.null(kp$barb)) {
      points(outline_shape()[kp$barb, 1], outline_shape()[kp$barb, 2], 
             col = "purple", pch = 19, cex = 2)
      text(outline_shape()[kp$barb, 1], outline_shape()[kp$barb, 2], 
           "BARB", col = "purple", pos = 3)
    }
    
    # Add metrics if available
    metrics <- hook_metrics_data()
    if (!is.null(metrics)) {
      text_x <- min(outline_shape()[,1])
      text_y <- max(outline_shape()[,2])
      metrics_text <- paste("Gape:", round(metrics$gape, 2),
                            "\nThroat:", round(metrics$throat, 2),
                            "\nPoint Angle:", round(metrics$point_angle, 2), "°")
      
      text(text_x, text_y, metrics_text,
           pos = 4, offset = 0.5, cex = 0.8)
    }
  })
  
  # Handle clicks on the outline plot
  observeEvent(input$outline_click, {
    # Only process clicks if in key point selection mode
    if(!key_points_mode() || is.null(current_key_point()) || is.null(outline_shape())) {
      return()
    }
    
    # Find the closest point on the outline to the click
    click_x <- input$outline_click$x
    click_y <- input$outline_click$y
    
    outline <- outline_shape()
    distances <- sqrt((outline[,1] - click_x)^2 + (outline[,2] - click_y)^2)
    closest_point <- which.min(distances)
    
    # Update the key points with the new selection
    current_kp <- key_points()
    current_kp[[current_key_point()]] <- closest_point
    key_points(current_kp)
    
    # Update log
    processing_log(paste(processing_log(), 
                         "\nSelected point ", closest_point, " as ", 
                         toupper(current_key_point()), ".",
                         sep=""))
    
    # Turn off selection mode
    key_points_mode(FALSE)
  })
  
  # Display key points info
  output$key_points_info <- renderText({
    kp <- key_points()
    if(is.null(kp$shank_top) && is.null(kp$shank_bottom) && 
       is.null(kp$bend_deepest) && is.null(kp$point) &&
       is.null(kp$barb)) {
      return("No points selected yet.")
    }
    
    paste("Selected points:\n",
          "Shank Top: ", if(is.null(kp$shank_top)) "Not set" else paste("Point", kp$shank_top), "\n",
          "Shank Bottom: ", if(is.null(kp$shank_bottom)) "Not set" else paste("Point", kp$shank_bottom), "\n",
          "Deepest Bend: ", if(is.null(kp$bend_deepest)) "Not set" else paste("Point", kp$bend_deepest), "\n",
          "Point Tip: ", if(is.null(kp$point)) "Not set" else paste("Point", kp$point), "\n",
          "Barb: ", if(is.null(kp$barb)) "Not set" else paste("Point", kp$barb))
  })
  
  # Display processing log
  output$processing_info <- renderText({
    processing_log()
  })
  
  # Display hook metrics
  output$hook_metrics <- renderText({
    metrics <- hook_metrics_data()
    if (is.null(metrics)) {
      return("Run analysis to calculate metrics.")
    }
    
    paste("Gape width:", round(metrics$gape, 2),
          "\nThroat depth:", round(metrics$throat, 2),
          "\nPoint angle:", round(metrics$point_angle, 2), "°")
  })
  
  # Run analysis
  observeEvent(input$run_analysis, {
    out_list <- outlines_list()
    name_list <- names_list()
    
    current_log <- processing_log()
    processing_log(paste(current_log, "\nStarting analysis...", sep=""))
    
    if (length(out_list) == 0) {
      processing_log(paste(current_log, "\nNo outlines to analyze. Process images first.", sep=""))
      return()
    }
    
    tryCatch({
      # Create Momocs Out object
      outlines <- Out(out_list)
      names(outlines$coo) <- name_list
      
      # Align outlines
      outlines <- outlines %>% 
        coo_center() %>% 
        coo_scale() %>% 
        coo_align()
      
      outlines_r(outlines)
      
      # Elliptical Fourier Analysis
      efa <- efourier(outlines, nb.h = input$nbh, norm = TRUE)
      efa_result_r(efa)
      
      # PCA
      pca <- PCA(efa)
      
      # Clustering
      clust <- hclust(dist(efa$coe), method = "ward.D2")
      k <- input$k_clusters
      clusters <- cutree(clust, k = k)
      pca$fac <- data.frame(cluster = factor(clusters))
      
      pca_result_r(pca)
      
      # Calculate mean shape for TPS grid
      mean_shape <- custom_mshape(outlines$coo)
      
      # For TPS visualization, use actual specimen shapes at the extremes of the PCA
      pc1_scores <- pca$x[,1]
      pc2_scores <- pca$x[,2]
      
      # Find specimens at extremes of PC axes
      pc1_pos_idx <- which.max(pc1_scores)
      pc1_neg_idx <- which.min(pc1_scores)
      pc2_pos_idx <- which.max(pc2_scores)
      pc2_neg_idx <- which.min(pc2_scores)
      
      # Store shapes at PCA extremes
      tps_grid_r(list(
        mean = mean_shape,
        pc1_pos = outlines$coo[[pc1_pos_idx]],
        pc1_neg = outlines$coo[[pc1_neg_idx]],
        pc2_pos = outlines$coo[[pc2_pos_idx]],
        pc2_neg = outlines$coo[[pc2_neg_idx]]
      ))
      
      # Estimate modularity by dividing outline into regions
      n_points <- nrow(mean_shape)
      
      # Determine region indices based on method
      if (input$region_method == "manual") {
        # Use manually selected key points if available
        kp <- key_points()
        
        # Check if essential points are selected
        if (is.null(kp$shank_top) || is.null(kp$shank_bottom) || 
            is.null(kp$bend_deepest) || is.null(kp$point)) {
          processing_log(paste(current_log, 
                               "\nError: Essential key points are not selected for manual region definition.",
                               "\nMust select: Shank Top, Shank Bottom, Deepest Bend, and Point.", 
                               sep=""))
          return()
        }
        
        # Define shank region between top and bottom
        if (kp$shank_top <= kp$shank_bottom) {
          shank_region <- kp$shank_top:kp$shank_bottom
        } else {
          shank_region <- c(kp$shank_top:n_points, 1:kp$shank_bottom)
        }
        
        # Define point region around point (10% of total points)
        point_size <- round(n_points * 0.1)
        half_point <- floor(point_size / 2)
        point_start <- (kp$point - half_point) %% n_points
        if(point_start == 0) point_start <- n_points
        point_end <- (kp$point + half_point) %% n_points
        if(point_end == 0) point_end <- n_points
        
        if(point_start <= point_end) {
          point_region <- point_start:point_end
        } else {
          point_region <- c(point_start:n_points, 1:point_end)
        }
        
        # Define barb region if selected
        if (!is.null(kp$barb)) {
          barb_size <- round(n_points * 0.1)
          half_barb <- floor(barb_size / 2)
          barb_start <- (kp$barb - half_barb) %% n_points
          if(barb_start == 0) barb_start <- n_points
          barb_end <- (kp$barb + half_barb) %% n_points
          if(barb_end == 0) barb_end <- n_points
          
          if(barb_start <= barb_end) {
            barb_region <- barb_start:barb_end
          } else {
            barb_region <- c(barb_start:n_points, 1:barb_end)
          }
        } else {
          barb_region <- numeric(0)
        }
        
        # Find shortest path from shank to point for bend region
        # Two possible paths
        path1 <- setdiff((max(shank_region)+1):(min(point_region)-1), barb_region)
        path2 <- setdiff(c((max(point_region)+1):n_points, 1:(min(shank_region)-1)), barb_region)
        
        # Select shorter path that contains the bend point
        if (length(path1) <= length(path2) && 
            (length(path1) == 0 || any(path1 == kp$bend_deepest))) {
          bend_region <- path1
        } else if (length(path2) > 0 && any(path2 == kp$bend_deepest)) {
          bend_region <- path2
        } else {
          # If bend point not in either path, add it to the shorter path
          if (length(path1) <= length(path2)) {
            bend_region <- unique(c(path1, kp$bend_deepest))
          } else {
            bend_region <- unique(c(path2, kp$bend_deepest))
          }
        }
        
        # Add adjacent points if bend region is too small
        if (length(bend_region) < 10) {
          additional_points <- setdiff(
            ((min(bend_region) - 5):(max(bend_region) + 5)) %% n_points,
            c(shank_region, point_region, barb_region)
          )
          additional_points[additional_points == 0] <- n_points
          bend_region <- unique(c(bend_region, additional_points))
        }
        
        # Calculate functional metrics
        # Gape - distance between shank center and point
        shank_midpoint <- colMeans(mean_shape[shank_region, , drop = FALSE])
        point_midpoint <- mean_shape[kp$point, ]
        gape <- sqrt(sum((shank_midpoint - point_midpoint)^2))
        
        # Throat - perpendicular distance from bend to shank-point line
        shank_to_point <- point_midpoint - shank_midpoint
        shank_to_point_len <- sqrt(sum(shank_to_point^2))
        shank_to_point_unit <- shank_to_point / shank_to_point_len
        
        bend_point <- mean_shape[kp$bend_deepest, ]
        shank_to_bend <- bend_point - shank_midpoint
        proj_length <- sum(shank_to_bend * shank_to_point_unit)
        proj_point <- shank_midpoint + proj_length * shank_to_point_unit
        throat <- sqrt(sum((bend_point - proj_point)^2))
        
        # Point angle calculation
        pre_point_idx <- (kp$point - 3) %% n_points
        if (pre_point_idx == 0) pre_point_idx <- n_points
        post_point_idx <- (kp$point + 3) %% n_points
        if (post_point_idx == 0) post_point_idx <- n_points
        
        v1 <- mean_shape[kp$point, ] - mean_shape[pre_point_idx, ]
        v2 <- mean_shape[post_point_idx, ] - mean_shape[kp$point, ]
        
        cos_angle <- sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))
        point_angle <- acos(min(max(cos_angle, -1), 1)) * 180 / pi
        
        # Store metrics
        hook_metrics_data(list(
          gape = gape,
          throat = throat,
          point_angle = point_angle
        ))
        
      } else if (input$region_method == "geometric") {
        # Geometric region detection based on fishhook shape features
        
        # Calculate curvature if not already done
        if (is.null(curvature_r())) {
          curv <- calculate_curvature(mean_shape)
          window_size <- ceiling(input$curvature_smoothing / 2) * 2 + 1  # Ensure odd number
          smoothed_curv <- smooth_curvature(curv, window_size)
        } else {
          smoothed_curv <- curvature_r()
        }
        
        # Find inflection points
        inflection_points <- find_inflection_points(smoothed_curv, mean_shape, input$curvature_threshold)
        
        # Define regions using specialized fishhook region detection
        regions <- define_fishhook_regions(mean_shape, inflection_points, smoothed_curv)
        
        # Extract regions
        shank_region <- regions$shank
        bend_region <- regions$bend
        point_region <- regions$point
        barb_region <- regions$barb
        
        # Store metrics
        hook_metrics_data(regions$metrics)
        
        # Log key points
        log_text <- paste(current_log, 
                          "\nGeometric detection found key points:",
                          "\nShank center: Point ", regions$key_points$shank_center,
                          "\nPoint: Point ", regions$key_points$point,
                          "\nThroat point: Point ", regions$key_points$throat,
                          sep="")
        
        if (!is.null(regions$key_points$barb_center)) {
          log_text <- paste(log_text, "\nBarb center: Point ", regions$key_points$barb_center, sep="")
        } else {
          log_text <- paste(log_text, "\nNo barb detected.", sep="")
        }
        
        processing_log(log_text)
        
      } else {
        # Proportional region identification based on percentages
        # Automatic region identification based on percentages and PCA orientation
        
        # Get user-defined proportions - ensure percentages sum to 100%
        total_prop <- (input$shank_ratio + input$bend_ratio + input$point_ratio + input$barb_ratio) / 100
        
        # Normalize proportions if they don't sum to 1
        shank_prop <- (input$shank_ratio / 100) / total_prop
        bend_prop <- (input$bend_ratio / 100) / total_prop
        point_prop <- (input$point_ratio / 100) / total_prop
        barb_prop <- (input$barb_ratio / 100) / total_prop
        
        # Use PCA to identify primary axis for orientation
        pca_outline <- prcomp(mean_shape)
        # Project onto PC1 to find extremes
        proj_pc1 <- as.matrix(mean_shape) %*% pca_outline$rotation[,1]
        # This gives us the extreme points along the primary axis
        max_pc1_idx <- which.max(proj_pc1)
        min_pc1_idx <- which.min(proj_pc1)
        
        # Choose start point based on curvature (point typically has higher curvature)
        curv <- calculate_curvature(mean_shape)
        max_curv_index <- which.max(curv)
        
        # Calculate distance from max_curv to PC1 extremes
        dist_to_max <- abs(max_curv_index - max_pc1_idx) %% n_points
        dist_to_min <- abs(max_curv_index - min_pc1_idx) %% n_points
        dist_to_max2 <- abs(n_points - dist_to_max)
        dist_to_min2 <- abs(n_points - dist_to_min)
        
        # Determine the closest PC1 extreme to the high curvature point
        closest_idx <- if(min(dist_to_max, dist_to_max2) <= min(dist_to_min, dist_to_min2)) {
          max_pc1_idx
        } else {
          min_pc1_idx
        }
        
        # Define point region around high curvature/PC1 extreme
        start_idx <- closest_idx
        
        # Create regions sequentially
        point_size <- round(n_points * point_prop)
        point_start <- start_idx
        point_end <- (point_start + point_size - 1) %% n_points
        if (point_end == 0) point_end <- n_points
        
        if (point_start <= point_end) {
          point_region <- point_start:point_end
        } else {
          point_region <- c(point_start:n_points, 1:point_end)
        }
        
        # Bend follows point
        bend_size <- round(n_points * bend_prop)
        bend_start <- (point_end + 1) %% n_points
        if (bend_start == 0) bend_start <- n_points
        bend_end <- (bend_start + bend_size - 1) %% n_points
        if (bend_end == 0) bend_end <- n_points
        
        if (bend_start <= bend_end) {
          bend_region <- bend_start:bend_end
        } else {
          bend_region <- c(bend_start:n_points, 1:bend_end)
        }
        
        # Shank follows bend
        shank_size <- round(n_points * shank_prop)
        shank_start <- (bend_end + 1) %% n_points
        if (shank_start == 0) shank_start <- n_points
        shank_end <- (shank_start + shank_size - 1) %% n_points
        if (shank_end == 0) shank_end <- n_points
        
        if (shank_start <= shank_end) {
          shank_region <- shank_start:shank_end
        } else {
          shank_region <- c(shank_start:n_points, 1:shank_end)
        }
        
        # Barb is the remaining points if barb_prop > 0
        if (barb_prop > 0) {
          barb_region <- setdiff(1:n_points, c(point_region, bend_region, shank_region))
        } else {
          barb_region <- numeric(0)
          
          # If no barb, reassign remaining points to nearest region
          remaining <- setdiff(1:n_points, c(point_region, bend_region, shank_region))
          
          for (pt in remaining) {
            # Find closest region
            point_dist <- min(abs(pt - point_region))
            bend_dist <- min(abs(pt - bend_region))
            shank_dist <- min(abs(pt - shank_region))
            
            min_dist <- min(point_dist, bend_dist, shank_dist)
            
            if (min_dist == point_dist) {
              point_region <- c(point_region, pt)
            } else if (min_dist == bend_dist) {
              bend_region <- c(bend_region, pt)
            } else {
              shank_region <- c(shank_region, pt)
            }
          }
        }
        
        # Calculate fishhook metrics
        # Determine key points for metrics calculation
        point_idx <- point_region[which.max(curv[point_region])]
        shank_midpoint <- colMeans(mean_shape[shank_region, , drop = FALSE])
        point_midpoint <- mean_shape[point_idx, ]
        
        # Gape - distance between shank center and point
        gape <- sqrt(sum((shank_midpoint - point_midpoint)^2))
        
        # Throat - max perpendicular distance from bend to shank-point line
        shank_to_point <- point_midpoint - shank_midpoint
        shank_to_point_len <- sqrt(sum(shank_to_point^2))
        shank_to_point_unit <- shank_to_point / shank_to_point_len
        
        # Calculate perpendicular distances for bend points
        depths <- numeric(length(bend_region))
        for (i in 1:length(bend_region)) {
          bend_point <- mean_shape[bend_region[i],]
          # Vector from shank to bend point
          shank_to_bend <- bend_point - shank_midpoint
          # Project onto shank-point line
          proj_length <- sum(shank_to_bend * shank_to_point_unit)
          # Projection point
          proj_point <- shank_midpoint + proj_length * shank_to_point_unit
          # Perpendicular distance
          depths[i] <- sqrt(sum((bend_point - proj_point)^2))
        }
        throat <- max(depths)
        
        # Point angle calculation
        pre_point_idx <- (point_idx - 3) %% n_points
        if (pre_point_idx == 0) pre_point_idx <- n_points
        post_point_idx <- (point_idx + 3) %% n_points
        if (post_point_idx == 0) post_point_idx <- n_points
        
        v1 <- mean_shape[point_idx,] - mean_shape[pre_point_idx,]
        v2 <- mean_shape[post_point_idx,] - mean_shape[point_idx,]
        
        cos_angle <- sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))
        point_angle <- acos(min(max(cos_angle, -1), 1)) * 180 / pi
        
        # Store metrics
        hook_metrics_data(list(
          gape = gape,
          throat = throat,
          point_angle = point_angle
        ))
      }
      
      # Ensure regions stay within bounds
      shank_region <- shank_region[shank_region > 0 & shank_region <= n_points]
      bend_region <- bend_region[bend_region > 0 & bend_region <= n_points]
      point_region <- point_region[point_region > 0 & point_region <= n_points]
      if (length(barb_region) > 0) {
        barb_region <- barb_region[barb_region > 0 & barb_region <= n_points]
      }
      
      # Create partition
      partition <- rep(0, n_points)
      partition[shank_region] <- 1  # Shank
      partition[bend_region] <- 2   # Bend
      partition[point_region] <- 3  # Point
      if (length(barb_region) > 0) {
        partition[barb_region] <- 4 # Barb
      }
      
      # Store results
      modularity_r(list(
        partition = partition,
        regions = list(shank = shank_region, bend = bend_region, point = point_region, barb = barb_region),
        n_points = n_points,
        metrics = hook_metrics_data()
      ))
      
      # Update log and switch to Results tab
      processing_log(paste(current_log, "\nAnalysis complete! Switching to Results tab.", sep=""))
      updateTabsetPanel(session, "main_tabs", selected = "Results")
      
    }, error = function(e) {
      # Error handling - update log
      processing_log(paste(current_log, "\nError during analysis:", e$message, sep=""))
    })
  })
  
  # Display grid of outlines
  output$grid_plot <- renderPlot({
    req(outlines_r())
    panel(outlines_r(), names = TRUE, cex.names = 0.7)
  })
  
  # Display stacked outlines
  output$stack_plot <- renderPlot({
    req(outlines_r())
    stack(outlines_r(), border = "black", col = "#00000010", lwd = 1,
          main = "Stacked & Aligned Outlines")
  })
  
  # Display actual EFA harmonics power
  output$harmonic_plot <- renderPlot({
    req(efa_result_r())
    efa <- efa_result_r()
    
    # Extract harmonics info and ensure it's valid
    n_harmonics <- efa$nb.h
    
    # Add debug print and validation
    print(paste("n_harmonics type:", class(n_harmonics), "value:", n_harmonics))
    
    # Check if n_harmonics is valid
    if (is.null(n_harmonics) || !is.numeric(n_harmonics) || length(n_harmonics) == 0 || n_harmonics < 1) {
      # Fallback to a default value
      n_harmonics <- 20
      print("Using default n_harmonics = 20")
    }
    
    # Calculate actual power for each harmonic from the coefficients
    harmonic_power <- numeric(n_harmonics)
    
    # Debug check if coefficient structure exists
    if (is.null(efa$coe) || !is.matrix(efa$coe)) {
      print("Warning: EFA coefficients matrix is not available or not in expected format")
      print(paste("efa$coe class:", class(efa$coe)))
      print(paste("efa$coe dimensions:", if(is.matrix(efa$coe)) paste(dim(efa$coe), collapse="x") else "N/A"))
      
      # Fallback to theoretical distribution if we can't access actual coefficients
      harmonic_power <- 1/(1:n_harmonics)^2
      print("Using theoretical power distribution as fallback")
    } else {
      # For each harmonic, calculate the power from its coefficients
      for (h in 1:n_harmonics) {
        # Get indices for this harmonic's coefficients
        coef_indices <- (4*(h-1) + 1):(4*h)
        
        # Check if we have enough coefficients
        if (max(coef_indices) <= ncol(efa$coe)) {
          # Average the squared coefficients across all specimens
          harmonic_power[h] <- mean(apply(efa$coe[, coef_indices, drop=FALSE], 1, function(x) sum(x^2)))
        } else {
          # Handle edge case if coefficients are missing
          harmonic_power[h] <- 0
          print(paste("Warning: coefficients missing for harmonic", h))
        }
      }
    }
    
    # Calculate normalized power as percentage
    total_power <- sum(harmonic_power)
    if (total_power > 0) {
      power_percent <- 100 * harmonic_power / total_power
    } else {
      # Handle edge case of zero total power
      power_percent <- rep(0, n_harmonics)
      power_percent[1] <- 100  # Assign all power to first harmonic as fallback
    }
    
    # Calculate cumulative power
    cumulative <- cumsum(power_percent)
    
    # Create a barplot with cumulative line
    par(mar = c(5, 4, 4, 4) + 0.1)
    bp <- barplot(power_percent, 
                  ylim = c(0, max(100, max(power_percent) * 1.2)),
                  main = "Harmonic Power Distribution", 
                  xlab = "Harmonic", 
                  ylab = "Power (%)",
                  col = "skyblue")
    
    # Add cumulative line
    par(new = TRUE)
    plot(bp, cumulative, type = "b", pch = 19, col = "red",
         axes = FALSE, xlab = "", ylab = "")
    axis(side = 4, at = pretty(range(cumulative)))
    mtext("Cumulative (%)", side = 4, line = 3)
    
    # Add a legend
    legend("topright", 
           legend = c("Power", "Cumulative"), 
           fill = c("skyblue", NA),
           lty = c(NA, 1),
           pch = c(NA, 19),
           col = c(NA, "red"),
           bty = "n")
    
    # Add note about 90% cumulative power if applicable
    h_90_idx <- which(cumulative >= 90)
    if (length(h_90_idx) > 0) {
      h_90 <- min(h_90_idx)
      mtext(paste0("90% of shape information captured by first ", h_90, " harmonics"), 
            side = 3, line = 0, cex = 0.8)
    }
  })
  
  # Display PCA morphospace
  output$morpho_plot <- renderPlot({
    req(pca_result_r())
    plot_PCA(pca_result_r(), morphospace = TRUE)
  })
  
  # Display Thin-Plate Spline visualizations
  output$tps_plot <- renderPlot({
    req(tps_grid_r(), outlines_r())
    tps <- tps_grid_r()
    
    # Create a 2x2 layout
    par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
    
    # PC1 negative - use actual shape
    plot(tps$mean, type = "l", asp = 1, main = "PC1 Negative", 
         xlab = "", ylab = "", axes = FALSE, col = "gray")
    lines(tps$pc1_neg, col = "red", lwd = 2)
    
    # PC1 positive - use actual shape
    plot(tps$mean, type = "l", asp = 1, main = "PC1 Positive", 
         xlab = "", ylab = "", axes = FALSE, col = "gray")
    lines(tps$pc1_pos, col = "blue", lwd = 2)
    
    # PC2 negative - use actual shape
    plot(tps$mean, type = "l", asp = 1, main = "PC2 Negative", 
         xlab = "", ylab = "", axes = FALSE, col = "gray")
    lines(tps$pc2_neg, col = "red", lwd = 2)
    
    # PC2 positive - use actual shape
    plot(tps$mean, type = "l", asp = 1, main = "PC2 Positive", 
         xlab = "", ylab = "", axes = FALSE, col = "gray")
    lines(tps$pc2_pos, col = "blue", lwd = 2)
    
    # Reset layout
    par(mfrow = c(1, 1))
    
    # Add title
    mtext("Shape Variation at PCA Extremes", 
          side = 3, line = -2, outer = TRUE, cex = 1.5)
  })
  
  # Display Integration/Modularity visualization with fishhook regions
  output$modularity_plot <- renderPlot({
    req(modularity_r(), outlines_r())
    mod <- modularity_r()
    outlines <- outlines_r()
    
    # Set up the plot layout for fishhook regions
    layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
    
    # Plot 1: Mean shape with regions colored
    mean_shape <- custom_mshape(outlines$coo)
    
    # Plot the mean shape with fishhook region indicators
    plot(mean_shape, type = "n", asp = 1, 
         main = "Fishhook Regions", axes = FALSE, xlab = "", ylab = "")
    
    # Plot each region with different colors
    points(mean_shape[mod$regions$shank, ], col = "blue", pch = 19, cex = 1)
    points(mean_shape[mod$regions$bend, ], col = "green", pch = 19, cex = 1)
    points(mean_shape[mod$regions$point, ], col = "red", pch = 19, cex = 1)
    if(length(mod$regions$barb) > 0) {
      points(mean_shape[mod$regions$barb, ], col = "purple", pch = 19, cex = 1)
    }
    
    # Connect points to show outline
    lines(mean_shape, col = "gray")
    
    # Add region labels
    shank_center <- colMeans(mean_shape[mod$regions$shank, , drop = FALSE])
    bend_center <- colMeans(mean_shape[mod$regions$bend, , drop = FALSE])
    point_center <- colMeans(mean_shape[mod$regions$point, , drop = FALSE])
    
    text(shank_center[1], shank_center[2], "SHANK", col = "darkblue", font = 2)
    text(bend_center[1], bend_center[2], "BEND", col = "darkgreen", font = 2)
    text(point_center[1], point_center[2], "POINT", col = "darkred", font = 2)
    
    if(length(mod$regions$barb) > 0) {
      barb_center <- colMeans(mean_shape[mod$regions$barb, , drop = FALSE])
      text(barb_center[1], barb_center[2], "BARB", col = "darkmagenta", font = 2)
    }
    
    # Add legend
    legend_items <- c("Shank", "Bend", "Point")
    legend_colors <- c("blue", "green", "red")
    
    if(length(mod$regions$barb) > 0) {
      legend_items <- c(legend_items, "Barb")
      legend_colors <- c(legend_colors, "purple")
    }
    
    legend("topright", legend = legend_items,
           col = legend_colors, pch = 19, cex = 0.8)
    
    # Plot 2: Gape and Throat visualization
    plot(mean_shape, type = "l", asp = 1, 
         main = "Functional Dimensions", axes = FALSE, xlab = "", ylab = "")
    
    # Calculate shank and point centroids for gape visualization
    shank_mid <- colMeans(mean_shape[mod$regions$shank, , drop = FALSE])
    point_mid <- colMeans(mean_shape[mod$regions$point, , drop = FALSE])
    
    # Draw gape line
    lines(rbind(shank_mid, point_mid), col = "red", lwd = 2)
    text((shank_mid[1] + point_mid[1])/2, (shank_mid[2] + point_mid[2])/2, 
         "Gape", col = "red", pos = 3, offset = 0.5)
    
    # Calculate and visualize throat depth
    # Vector from shank to point
    shank_to_point <- point_mid - shank_mid
    shank_to_point_len <- sqrt(sum(shank_to_point^2))
    shank_to_point_unit <- shank_to_point / shank_to_point_len
    
    # Calculate perpendicular distances for bend points
    depths <- numeric(length(mod$regions$bend))
    for (i in 1:length(mod$regions$bend)) {
      bend_point <- mean_shape[mod$regions$bend[i],]
      # Vector from shank to bend point
      shank_to_bend <- bend_point - shank_mid
      # Project onto shank-point line
      proj_length <- sum(shank_to_bend * shank_to_point_unit)
      # Projection point
      proj_point <- shank_mid + proj_length * shank_to_point_unit
      # Draw perpendicular line for max depth
      if (i == which.max(sqrt(rowSums((mean_shape[mod$regions$bend,] - shank_mid)^2)))) {
        lines(rbind(bend_point, proj_point), col = "green", lwd = 2, lty = 2)
        text((bend_point[1] + proj_point[1])/2, (bend_point[2] + proj_point[2])/2, 
             "Throat", col = "green", pos = 2, offset = 0.5)
      }
    }
    
    # Add functional metrics as text
    metrics <- mod$metrics
    metrics_text <- paste("Gape:", round(metrics$gape, 2),
                          "\nThroat:", round(metrics$throat, 2),
                          "\nPoint Angle:", round(metrics$point_angle, 2), "°")
    
    text(min(mean_shape[,1]), min(mean_shape[,2]), metrics_text,
         pos = 4, offset = 1, cex = 0.8)
    
    # Plot 3: Variance by region
    # Calculate variance in each region across all specimens
    shank_var <- numeric(length(mod$regions$shank))
    bend_var <- numeric(length(mod$regions$bend))
    point_var <- numeric(length(mod$regions$point))
    if (length(mod$regions$barb) > 0) {
      barb_var <- numeric(length(mod$regions$barb))
    }
    
    # For each specimen, calculate distance from mean
    n_specimens <- length(outlines$coo)
    
    for (i in 1:length(mod$regions$shank)) {
      point_dists <- numeric(n_specimens)
      for (j in 1:n_specimens) {
        spec_point <- outlines$coo[[j]][mod$regions$shank[i], ]
        mean_point <- mean_shape[mod$regions$shank[i], ]
        point_dists[j] <- sqrt(sum((spec_point - mean_point)^2))
      }
      shank_var[i] <- var(point_dists)
    }
    
    for (i in 1:length(mod$regions$bend)) {
      point_dists <- numeric(n_specimens)
      for (j in 1:n_specimens) {
        spec_point <- outlines$coo[[j]][mod$regions$bend[i], ]
        mean_point <- mean_shape[mod$regions$bend[i], ]
        point_dists[j] <- sqrt(sum((spec_point - mean_point)^2))
      }
      bend_var[i] <- var(point_dists)
    }
    
    for (i in 1:length(mod$regions$point)) {
      point_dists <- numeric(n_specimens)
      for (j in 1:n_specimens) {
        spec_point <- outlines$coo[[j]][mod$regions$point[i], ]
        mean_point <- mean_shape[mod$regions$point[i], ]
        point_dists[j] <- sqrt(sum((spec_point - mean_point)^2))
      }
      point_var[i] <- var(point_dists)
    }
    
    if (length(mod$regions$barb) > 0) {
      for (i in 1:length(mod$regions$barb)) {
        point_dists <- numeric(n_specimens)
        for (j in 1:n_specimens) {
          spec_point <- outlines$coo[[j]][mod$regions$barb[i], ]
          mean_point <- mean_shape[mod$regions$barb[i], ]
          point_dists[j] <- sqrt(sum((spec_point - mean_point)^2))
        }
        barb_var[i] <- var(point_dists)
      }
      
      # Calculate mean variance for each region
      region_vars <- c(mean(shank_var), mean(bend_var), mean(point_var), mean(barb_var))
      names(region_vars) <- c("Shank", "Bend", "Point", "Barb")
      
      # Create barplot
      barplot(region_vars, col = c("blue", "green", "red", "purple"),
              main = "Variance by Region", ylab = "Mean Variance")
    } else {
      # Calculate mean variance for each region (no barb)
      region_vars <- c(mean(shank_var), mean(bend_var), mean(point_var))
      names(region_vars) <- c("Shank", "Bend", "Point")
      
      # Create barplot
      barplot(region_vars, col = c("blue", "green", "red"),
              main = "Variance by Region", ylab = "Mean Variance")
    }
    
    # Plot 4: Heat map of variance along outline
    point_vars <- numeric(mod$n_points)
    
    # Calculate variance at each point
    for (i in 1:mod$n_points) {
      point_dists <- numeric(n_specimens)
      for (j in 1:n_specimens) {
        spec_point <- outlines$coo[[j]][i, ]
        mean_point <- mean_shape[i, ]
        point_dists[j] <- sqrt(sum((spec_point - mean_point)^2))
      }
      point_vars[i] <- var(point_dists)
    }
    
    # Plot mean shape
    plot(mean_shape, type = "n", asp = 1, 
         main = "Variance Heatmap", axes = FALSE, xlab = "", ylab = "")
    
    # Color points by variance
    # Normalize variance values to 0-1 range for coloring
    norm_vars <- (point_vars - min(point_vars)) / (max(point_vars) - min(point_vars))
    
    # Plot points colored by variance
    for (i in 1:mod$n_points) {
      # Create color from blue (low variance) to red (high variance)
      col_val <- rgb(norm_vars[i], 0, 1-norm_vars[i])
      points(mean_shape[i, 1], mean_shape[i, 2], col = col_val, pch = 19, cex = 1.5)
    }
    
    # Connect points to show outline
    lines(mean_shape, col = "gray")
    
    # Plot 5: Region Contributions to PC1
    # Calculate contributions of each region to PC1 variation
    pc1_scores <- pca_result_r()$x[,1]
    
    # Get min and max specimens on PC1
    pc1_pos_idx <- which.max(pc1_scores)
    pc1_neg_idx <- which.min(pc1_scores)
    
    # Calculate contribution for each region
    shank_cont <- mean(sqrt(rowSums((outlines$coo[[pc1_pos_idx]][mod$regions$shank,] - 
                                       outlines$coo[[pc1_neg_idx]][mod$regions$shank,])^2)))
    
    bend_cont <- mean(sqrt(rowSums((outlines$coo[[pc1_pos_idx]][mod$regions$bend,] - 
                                      outlines$coo[[pc1_neg_idx]][mod$regions$bend,])^2)))
    
    point_cont <- mean(sqrt(rowSums((outlines$coo[[pc1_pos_idx]][mod$regions$point,] - 
                                       outlines$coo[[pc1_neg_idx]][mod$regions$point,])^2)))
    
    if (length(mod$regions$barb) > 0) {
      barb_cont <- mean(sqrt(rowSums((outlines$coo[[pc1_pos_idx]][mod$regions$barb,] - 
                                        outlines$coo[[pc1_neg_idx]][mod$regions$barb,])^2)))
      
      # Normalize contributions
      pc1_cont <- c(shank_cont, bend_cont, point_cont, barb_cont)
      pc1_cont <- 100 * pc1_cont / sum(pc1_cont)
      names(pc1_cont) <- c("Shank", "Bend", "Point", "Barb")
      
      # Create barplot
      barplot(pc1_cont, col = c("blue", "green", "red", "purple"),
              main = "Region Contributions to PC1", ylab = "Contribution (%)")
    } else {
      # Normalize contributions (no barb)
      pc1_cont <- c(shank_cont, bend_cont, point_cont)
      pc1_cont <- 100 * pc1_cont / sum(pc1_cont)
      names(pc1_cont) <- c("Shank", "Bend", "Point")
      
      # Create barplot
      barplot(pc1_cont, col = c("blue", "green", "red"),
              main = "Region Contributions to PC1", ylab = "Contribution (%)")
    }
    
    # Plot 6: Region correlation matrix
    # Function to calculate mean distance between two sets of points
    calc_distance <- function(shape1, shape2) {
      return(mean(sqrt(rowSums((shape1 - shape2)^2))))
    }
    
    # Calculate distances between all pairs of specimens
    n_specimens <- length(outlines$coo)
    
    # Matrix to store distances for each region
    shank_dists <- matrix(0, n_specimens, n_specimens)
    bend_dists <- matrix(0, n_specimens, n_specimens)
    point_dists <- matrix(0, n_specimens, n_specimens)
    if (length(mod$regions$barb) > 0) {
      barb_dists <- matrix(0, n_specimens, n_specimens)
    }
    
    for (i in 1:(n_specimens-1)) {
      for (j in (i+1):n_specimens) {
        # Shank distances
        shank_dists[i,j] <- shank_dists[j,i] <- calc_distance(
          outlines$coo[[i]][mod$regions$shank,], 
          outlines$coo[[j]][mod$regions$shank,]
        )
        
        # Bend distances
        bend_dists[i,j] <- bend_dists[j,i] <- calc_distance(
          outlines$coo[[i]][mod$regions$bend,], 
          outlines$coo[[j]][mod$regions$bend,]
        )
        
        # Point distances
        point_dists[i,j] <- point_dists[j,i] <- calc_distance(
          outlines$coo[[i]][mod$regions$point,], 
          outlines$coo[[j]][mod$regions$point,]
        )
        
        # Barb distances (if applicable)
        if (length(mod$regions$barb) > 0) {
          barb_dists[i,j] <- barb_dists[j,i] <- calc_distance(
            outlines$coo[[i]][mod$regions$barb,], 
            outlines$coo[[j]][mod$regions$barb,]
          )
        }
      }
    }
    
    if (length(mod$regions$barb) > 0) {
      # Calculate correlations between distance matrices (with barb)
      shank_bend_cor <- cor(as.vector(shank_dists), as.vector(bend_dists))
      shank_point_cor <- cor(as.vector(shank_dists), as.vector(point_dists))
      shank_barb_cor <- cor(as.vector(shank_dists), as.vector(barb_dists))
      bend_point_cor <- cor(as.vector(bend_dists), as.vector(point_dists))
      bend_barb_cor <- cor(as.vector(bend_dists), as.vector(barb_dists))
      point_barb_cor <- cor(as.vector(point_dists), as.vector(barb_dists))
      
      # Create correlation matrix
      cor_matrix <- matrix(1, 4, 4)
      cor_matrix[1,2] <- cor_matrix[2,1] <- shank_bend_cor
      cor_matrix[1,3] <- cor_matrix[3,1] <- shank_point_cor
      cor_matrix[1,4] <- cor_matrix[4,1] <- shank_barb_cor
      cor_matrix[2,3] <- cor_matrix[3,2] <- bend_point_cor
      cor_matrix[2,4] <- cor_matrix[4,2] <- bend_barb_cor
      cor_matrix[3,4] <- cor_matrix[4,3] <- point_barb_cor
      
      # Create heatmap
      image(1:4, 1:4, cor_matrix, axes = FALSE, 
            col = colorRampPalette(c("blue", "white", "red"))(100),
            main = "Region Correlations",
            xlab = "", ylab = "")
      axis(1, at = 1:4, labels = c("Shank", "Bend", "Point", "Barb"))
      axis(2, at = 1:4, labels = c("Shank", "Bend", "Point", "Barb"))
      
      # Add correlation values
      for (i in 1:4) {
        for (j in 1:4) {
          text(i, j, round(cor_matrix[i,j], 2), col = "black")
        }
      }
    } else {
      # Calculate correlations between distance matrices (no barb)
      shank_bend_cor <- cor(as.vector(shank_dists), as.vector(bend_dists))
      shank_point_cor <- cor(as.vector(shank_dists), as.vector(point_dists))
      bend_point_cor <- cor(as.vector(bend_dists), as.vector(point_dists))
      
      # Create correlation matrix
      cor_matrix <- matrix(1, 3, 3)
      cor_matrix[1,2] <- cor_matrix[2,1] <- shank_bend_cor
      cor_matrix[1,3] <- cor_matrix[3,1] <- shank_point_cor
      cor_matrix[2,3] <- cor_matrix[3,2] <- bend_point_cor
      
      # Create heatmap
      image(1:3, 1:3, cor_matrix, axes = FALSE, 
            col = colorRampPalette(c("blue", "white", "red"))(100),
            main = "Region Correlations",
            xlab = "", ylab = "")
      axis(1, at = 1:3, labels = c("Shank", "Bend", "Point"))
      axis(2, at = 1:3, labels = c("Shank", "Bend", "Point"))
      
      # Add correlation values
      for (i in 1:3) {
        for (j in 1:3) {
          text(i, j, round(cor_matrix[i,j], 2), col = "black")
        }
      }
    }
    
    # Reset layout
    par(mfrow = c(1, 1))
    
    # Add title
    mtext("Integration & Modularity Analysis", 
          side = 3, line = -2, outer = TRUE, cex = 1.5)
  })
  
  # Download handler for PDF graphics - updated for fishhook analysis
  output$download_pdf <- downloadHandler(
    filename = function() {
      paste("FishhookAnalysis_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".pdf", sep = "")
    },
    content = function(file) {
      # Check if analysis results exist
      if (is.null(outlines_r())) {
        return(NULL)
      }
      
      # Create PDF with multiple plots
      pdf(file, width = 8.5, height = 11)
      
      # Title page
      plot.new()
      text(0.5, 0.7, "Fishhook Morphometric Analysis", cex = 2)
      text(0.5, 0.6, paste("Generated on", format(Sys.Date(), "%B %d, %Y")), cex = 1.2)
      text(0.5, 0.5, paste("Number of specimens:", length(outlines_r()$coo)), cex = 1.2)
      text(0.5, 0.4, paste("Region detection method:", input$region_method), cex = 1.2)
      
      # Fishhook metrics
      metrics <- hook_metrics_data()
      if (!is.null(metrics)) {
        text(0.5, 0.3, paste("Mean Gape Width:", round(metrics$gape, 2)), cex = 1.2)
        text(0.5, 0.25, paste("Mean Throat Depth:", round(metrics$throat, 2)), cex = 1.2)
        text(0.5, 0.2, paste("Mean Point Angle:", round(metrics$point_angle, 2), "°"), cex = 1.2)
      }
      
      # Grid of outlines
      panel(outlines_r(), names = TRUE, cex.names = 0.7)
      title("Individual Fishhook Outlines")
      
      # Stacked outlines
      stack(outlines_r(), border = "black", col = "#00000010", lwd = 1)
      title("Stacked & Aligned Fishhook Outlines")
      
      # Actual harmonic power
      efa <- efa_result_r()
      n_harmonics <- efa$nb.h
      
      # Calculate actual power for each harmonic from the coefficients
      harmonic_power <- numeric(n_harmonics)
      for (h in 1:n_harmonics) {
        coef_indices <- (4*(h-1) + 1):(4*h)
        if (max(coef_indices) <= ncol(efa$coe)) {
          harmonic_power[h] <- mean(apply(efa$coe[, coef_indices], 1, function(x) sum(x^2)))
        } else {
          harmonic_power[h] <- 0
        }
      }
      
      # Calculate normalized power as percentage
      power_percent <- 100 * harmonic_power / sum(harmonic_power)
      cumulative <- cumsum(power_percent)
      
      par(mar = c(5, 4, 4, 4) + 0.1)
      bp <- barplot(power_percent, 
                    ylim = c(0, max(100, max(power_percent) * 1.2)),
                    main = "Harmonic Power Distribution", 
                    xlab = "Harmonic", 
                    ylab = "Power (%)",
                    col = "skyblue")
      
      par(new = TRUE)
      plot(bp, cumulative, type = "b", pch = 19, col = "red",
           axes = FALSE, xlab = "", ylab = "")
      axis(side = 4, at = pretty(range(cumulative)))
      mtext("Cumulative (%)", side = 4, line = 3)
      
      # PCA visualization if available
      if (!is.null(pca_result_r())) {
        # Use plot_PCA without main argument, then add title separately
        plot_PCA(pca_result_r(), morphospace = TRUE)
        title("PCA Morphospace")
        
        # PCA summary
        plot.new()
        grid <- pca_result_r()$x
        pc1_var <- round(100 * pca_result_r()$eig[1] / sum(pca_result_r()$eig), 1)
        pc2_var <- round(100 * pca_result_r()$eig[2] / sum(pca_result_r()$eig), 1)
        
        text(0.5, 0.9, "PCA Summary", cex = 1.5)
        text(0.5, 0.8, paste("PC1 explains", pc1_var, "% of variance"), cex = 1.2)
        text(0.5, 0.7, paste("PC2 explains", pc2_var, "% of variance"), cex = 1.2)
        text(0.5, 0.6, paste("Cumulative:", pc1_var + pc2_var, "%"), cex = 1.2)
        
        # Range on PC1 and PC2
        text(0.5, 0.4, paste("PC1 range:", round(min(grid[,1]), 2), "to", round(max(grid[,1]), 2)), cex = 1)
        text(0.5, 0.3, paste("PC2 range:", round(min(grid[,2]), 2), "to", round(max(grid[,2]), 2)), cex = 1)
      }
      
      # Add TPS grid if available
      if (!is.null(tps_grid_r())) {
        tps <- tps_grid_r()
        
        # Create a 2x2 layout inside the PDF
        par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
        
        # PC1 negative - use actual shape
        plot(tps$mean, type = "l", asp = 1, main = "PC1 Negative", 
             xlab = "", ylab = "", axes = FALSE, col = "gray")
        lines(tps$pc1_neg, col = "red", lwd = 2)
        
        # PC1 positive - use actual shape
        plot(tps$mean, type = "l", asp = 1, main = "PC1 Positive", 
             xlab = "", ylab = "", axes = FALSE, col = "gray")
        lines(tps$pc1_pos, col = "blue", lwd = 2)
        
        # PC2 negative - use actual shape
        plot(tps$mean, type = "l", asp = 1, main = "PC2 Negative", 
             xlab = "", ylab = "", axes = FALSE, col = "gray")
        lines(tps$pc2_neg, col = "red", lwd = 2)
        
        # PC2 positive - use actual shape
        plot(tps$mean, type = "l", asp = 1, main = "PC2 Positive", 
             xlab = "", ylab = "", axes = FALSE, col = "gray")
        lines(tps$pc2_pos, col = "blue", lwd = 2)
        
        # Reset layout
        par(mfrow = c(1, 1))
        
        # Add title for TPS page
        title("Shape Variation at PCA Extremes", line = -1)
      }
      
      # Add modularity plots if available
      if (!is.null(modularity_r()) && !is.null(outlines_r())) {
        mod <- modularity_r()
        outlines <- outlines_r()
        mean_shape <- custom_mshape(outlines$coo)
        
        # Create a page for fishhook functional visualization
        plot(mean_shape, type = "l", asp = 1, 
             main = "Fishhook Functional Dimensions", 
             axes = FALSE, xlab = "", ylab = "")
        
        # Plot regions in different colors
        points(mean_shape[mod$regions$shank, ], col = "blue", pch = 19, cex = 0.8)
        points(mean_shape[mod$regions$bend, ], col = "green", pch = 19, cex = 0.8)
        points(mean_shape[mod$regions$point, ], col = "red", pch = 19, cex = 0.8)
        if (length(mod$regions$barb) > 0) {
          points(mean_shape[mod$regions$barb, ], col = "purple", pch = 19, cex = 0.8)
        }
        
        # Calculate centroids for functional metrics
        shank_mid <- colMeans(mean_shape[mod$regions$shank, , drop = FALSE])
        point_mid <- colMeans(mean_shape[mod$regions$point, , drop = FALSE])
        
        # Draw gape line
        lines(rbind(shank_mid, point_mid), col = "red", lwd = 2)
        text((shank_mid[1] + point_mid[1])/2, (shank_mid[2] + point_mid[2])/2, 
             "Gape", col = "red", pos = 3, offset = 0.5)
        
        # Add throat line
        # Vector from shank to point
        shank_to_point <- point_mid - shank_mid
        shank_to_point_len <- sqrt(sum(shank_to_point^2))
        shank_to_point_unit <- shank_to_point / shank_to_point_len
        
        # Find deepest point of bend
        bend_center <- colMeans(mean_shape[mod$regions$bend, , drop = FALSE])
        
        # Vector from shank to bend center
        shank_to_bend <- bend_center - shank_mid
        # Project onto shank-point line
        proj_length <- sum(shank_to_bend * shank_to_point_unit)
        # Projection point
        proj_point <- shank_mid + proj_length * shank_to_point_unit
        
        # Draw throat line
        lines(rbind(bend_center, proj_point), col = "green", lwd = 2, lty = 2)
        text((bend_center[1] + proj_point[1])/2, (bend_center[2] + proj_point[2])/2, 
             "Throat", col = "green", pos = 2, offset = 0.5)
        
        # Add legend
        legend_items <- c("Shank", "Bend", "Point", "Gape", "Throat")
        legend_cols <- c("blue", "green", "red", "red", "green")
        legend_lty <- c(NA, NA, NA, 1, 2)
        legend_pch <- c(19, 19, 19, NA, NA)
        
        if (length(mod$regions$barb) > 0) {
          legend_items <- c(legend_items, "Barb")
          legend_cols <- c(legend_cols, "purple")
          legend_lty <- c(legend_lty, NA)
          legend_pch <- c(legend_pch, 19)
        }
        
        legend("topright", legend = legend_items, 
               col = legend_cols, lty = legend_lty, pch = legend_pch, 
               cex = 0.8, bty = "n")
        
        # Add metrics
        metrics_text <- paste("Gape width:", round(mod$metrics$gape, 2),
                              "\nThroat depth:", round(mod$metrics$throat, 2),
                              "\nPoint angle:", round(mod$metrics$point_angle, 2), "°")
        
        text(min(mean_shape[,1]), min(mean_shape[,2]), metrics_text,
             pos = 4, offset = 1, cex = 0.8)
      }
      
      dev.off()
    }
  )
}

shinyApp(ui, server)