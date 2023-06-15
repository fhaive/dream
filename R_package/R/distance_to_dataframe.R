#' Convert a distance matrix to a data frame.
#'
#' This is an internal function that converts a distance matrix to a data frame format. It considers the distance matrix to be symmetric and only uses the upper triangle of the matrix.
#'
#' @param distance_matrix A symmetric distance matrix.
#' @return A data frame with three columns: "drug1", "drug2", and "distance", representing pairs of drugs and their corresponding distances.
#' @export
distance_to_dataframe <- function(distance_matrix) {
    # we consider the distance matrix symmetric, 
    # and we only use the upper triangle
    distance_matrix[lower.tri(distance_matrix, diag=TRUE)] <- NA
    dist_df = na.omit(data.frame(as.table(distance_matrix)))
    colnames(dist_df) <- c("drug1", "drug2", "distance")
    dist_df
}