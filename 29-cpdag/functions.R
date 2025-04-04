#' Convert a graphNEL object to Mermaid diagram code
#'
#' This function takes a \code{graphNEL} object and converts it to a Mermaid diagram
#' code string with a left-to-right ("graph LR") orientation. For each pair of nodes,
#' it prints an undirected edge if both directions exist; otherwise, it prints the
#' directed edge that exists.
#'
#' @param graph A \code{graphNEL} object representing the graph.
#'
#' @return A character string representing the Mermaid diagram.
#'
#' @examples
#' \dontrun{
#'   # Create a sample graphNEL object (using the graph package)
#'   library(graph)
#'   nodes <- c("A", "B", "C")
#'   edgeL <- list(
#'     A = list(edges = c(2, 3)),
#'     B = list(edges = c(3)),
#'     C = list(edges = c(1))
#'   )
#'   g <- new("graphNEL", nodes = nodes, edgeL = edgeL, edgemode = "directed")
#'   mermaid_code <- graphNEL_to_mermaid(g)
#'   cat(mermaid_code)
#' }
graphNEL_to_mermaid <- function(graph) {
  # Retrieve the vector of node names from the graphNEL object.
  nodes_vec <- nodes(graph)
  # Initialize the Mermaid diagram string with left-to-right orientation.
  mermaid_str <- "graph LR"  # Set LR orientation.
  n <- length(nodes_vec)

  # Loop over all unordered pairs of nodes to identify directed or undirected edges.
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      node_i <- nodes_vec[i]
      node_j <- nodes_vec[j]
      # Retrieve the list of neighbor indices for each node.
      edges_i <- graph@edgeL[[node_i]]$edges
      edges_j <- graph@edgeL[[node_j]]$edges

      # Check for the presence of edges in both directions.
      directed_i_j <- j %in% edges_i
      directed_j_i <- i %in% edges_j

      if (directed_i_j && directed_j_i) {
        # If both directions exist, print an undirected edge.
        mermaid_str <- paste0(mermaid_str, "\n    ", node_i, " --- ", node_j)
      } else if (directed_i_j) {
        # Only edge from node_i to node_j exists.
        mermaid_str <- paste0(mermaid_str, "\n    ", node_i, " --> ", node_j)
      } else if (directed_j_i) {
        # Only edge from node_j to node_i exists.
        mermaid_str <- paste0(mermaid_str, "\n    ", node_j, " --> ", node_i)
      }
    }
  }

  mermaid_str
}

#' Convert a CPDAG object to an adjacency matrix
#'
#' This function converts a completed partially directed acyclic graph (CPDAG)
#' into an adjacency matrix representation. The matrix is transposed before returning,
#' since the orientation used by the \code{pcalg} package is opposite to that of the
#' \code{graphNEL} class.
#'
#' @param cpdag A CPDAG object, typically represented as a \code{graphNEL} object.
#'
#' @return A square matrix of 0s and 1s where a 1 indicates the presence of an edge from
#'   the row node to the column node.
#'
#' @examples
#' \dontrun{
#'   # Assuming cpdag is a CPDAG object from a causal discovery package:
#'   adj_matrix <- cpdag_to_adj(cpdag)
#'   print(adj_matrix)
#' }
cpdag_to_adj <- function(cpdag) {
  # Retrieve the node names from the CPDAG object.
  nodes_vec <- nodes(cpdag)
  n <- length(nodes_vec)
  # Initialize an n x n adjacency matrix filled with zeros.
  adj <- matrix(0, n, n)
  rownames(adj) <- colnames(adj) <- nodes_vec

  # Loop over each node to update the adjacency matrix using the edge list.
  for (node in nodes_vec) {
    targets <- cpdag@edgeL[[node]]$edges
    if (!is.null(targets) && length(targets) > 0) {
      for (target_index in targets) {
        target <- nodes_vec[target_index]
        # Set the entry to 1 to indicate an edge from 'node' to 'target'.
        adj[node, target] <- 1
      }
    }
  }

  # Transpose the matrix to adjust for the orientation used by the pcalg package.
  t(adj)
}

#' Convert multiple DAGs to Mermaid diagram code
#'
#' This function converts a collection of Directed Acyclic Graphs (DAGs), each represented
#' as a flattened adjacency matrix, into Mermaid diagram code strings. Each DAG is reshaped
#' into an adjacency matrix and then converted into a Mermaid string using left-to-right orientation.
#'
#' @param allDags_res A list containing:
#'   \item{nodeNms}{A character vector of node names.}
#'   \item{dags}{A matrix where each row is a flattened adjacency matrix of a DAG, assumed to be
#'   in column-major order.}
#'
#' @return A list of character strings, each representing a Mermaid diagram for a DAG.
#'
#' @examples
#' \dontrun{
#'   # Example structure of allDags_res:
#'   # allDags_res$nodeNms = c("A", "B", "C")
#'   # allDags_res$dags is a matrix with each row representing a DAG
#'   mermaid_codes <- allDags_to_mermaid(allDags_res)
#'   cat(mermaid_codes[[1]])
#' }
allDags_to_mermaid <- function(allDags_res) {
  # Retrieve node names and determine the number of nodes.
  nodes <- allDags_res$nodeNms
  n <- length(nodes)

  # Determine the number of DAGs from the rows of the 'dags' matrix.
  num_dags <- nrow(allDags_res$dags)
  dag_mermaid <- vector("list", num_dags)

  # Process each DAG by reshaping its flattened adjacency matrix into a proper matrix.
  for (i in 1:num_dags) {
    dag_vec <- allDags_res$dags[i, ]
    # Reshape the vector into an n x n matrix assuming column-major order.
    dag_mat <- matrix(dag_vec, nrow = n, ncol = n, byrow = FALSE)
    rownames(dag_mat) <- nodes
    colnames(dag_mat) <- nodes

    # Start building the Mermaid diagram with left-to-right orientation.
    mermaid_str <- "graph LR"
    # Loop over the adjacency matrix to create directed edges.
    for (r in 1:n) {
      for (c in 1:n) {
        if (dag_mat[r, c] == 1) {
          mermaid_str <- paste0(mermaid_str, "\n    ", nodes[r], " --> ", nodes[c])
        }
      }
    }
    dag_mermaid[[i]] <- mermaid_str
  }

  dag_mermaid
}
