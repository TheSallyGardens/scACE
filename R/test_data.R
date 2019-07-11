#' test_data for pkgGCGFunction package
#'
#' A dataset containing single- cell chromatin accessibility data and gene expression data
#' The variables are as follows:
#'
#' \itemize{
#'   \item data_acc: chromatin accessibility data matrix, cells by regions
#'   \item data_exp: gene expression data matrix, cells by genes
#'   \item overlap_seq_acc: the linked features in chromatin accessibility data
#'   \item overlap_seq_exp: the linked features in gene expression data
#'   \item acc_true_cluster: the true cell labels for chromatin accessibility data, 1: K562 cell line, 2: HL60 cell line.
#'   \item exp_true_cluster: the true cell labels for gene expression data, 1: K562 cell line, 2: HL60 cell line.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name LZX_data
#' @usage pkgGCGFunction::LZX_data
#' @format large list (6 elements, 2 Mb)
"LZX_data"
