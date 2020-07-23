#' Original sequence position from aligned sequence
#'
#' Internal use.
#'
#' @param seq An aligned sequence in vector form. Gaps as "-".
#'
#' @return A tibble
seq_pos_from_aln <- function(seq){
  positions <- rep(NA, length(seq))
  ii <- seq != "-"
  positions[ii] <- cumsum(ii)[ii]

  tibble::tibble(char = as.character(seq),
                 aln_pos = 1:length(seq),
                 seq_pos = positions)
}

#' Match external features to alignment
#'
#' For internal use.
#'
#' @param aln Alignment in tibble form. Output from seq_pos_from_aln.
#' Must have columns 'seq_id', 'seq_pos'.
#' @param feats Feature table. Must have columns 'seq_id', 'feat_id',
#' 'start', 'end'. It should not have overlapping features in the same
#' sequence.
#' @importFrom magrittr %>%
#'
#' @return
feats_to_aln <- function(aln, feats){

  # Prepare column for features
  aln <- aln %>%
    dplyr::mutate(feat_id = NA)

  # Map features to aln
  for(f in unique(feats$feat_id)){
    # Select feature
    d <- feats %>%
      dplyr::filter(feat_id == f)

    # Map feature to aln
    aln <- aln %>%
      dplyr::left_join(d, by = "seq_id") %>%
      dplyr::mutate(feat_id.y = replace(feat_id.y,
                                        is.na(seq_pos) | seq_pos < start | seq_pos > end,
                                        NA))

    # Check overlapping features
    n_dups <- aln %>%
      dplyr::filter(!is.na(feat_id.x) & !is.na(feat_id.y)) %>%
      nrow()
    if(n_dups > 0)
      stop("ERROR: Some of your features seem to overlap within the same sequence.", call. = TRUE)

    aln <- aln %>%
      dplyr::mutate(feat_id = replace(feat_id.x,
                                      is.na(feat_id.x) & !is.na(feat_id.y),
                                      feat_id.y[ is.na(feat_id.x) & !is.na(feat_id.y) ])) %>%
      # dplyr::filter(!is.na(feat_id.x) | !is.na(feat_id.y))
      dplyr::select(-feat_id.y, -feat_id.x, -start, -end)
  }

  return(aln)
}

#' Plot alignment with features
#'
#' Takes an alignment in fasta format as defined in
#' seqinr, and a tibble of features. Produces plot with
#' annotated features.
#'
#' @param aln Alignment in fasta format from seqinr. See
#' examples.
#' @param feats Feature table. Must have columns 'seq_id', 'feat_id',
#' 'start', 'end'. It should not have overlapping features in the same
#' @param region Start and end positions of range to plot. If NULL, plot
#' all the alignment.
#'
#' @return A ggplot2 plot
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' # Creating alignment as defined by seqinr read.fasta
#' aln <- list(seq1 = c(rep('A', 10), rep('B', 10), rep('C', 10)),
#' seq2 = c(rep('A', 9), rep('B', 10), rep('C', 11)),
#' seq3 = c(rep('A', 11), rep('B', 12), rep('C', 7)))
#' aln
#'
#' # Table of features
#' feats <- tibble::tibble(seq_id = paste0('seq', 1:3),
#'                         feat_id = 'feat1',
#'                         start = c(11, 10, 12),
#'                         end = c(20, 19, 23)) %>%
#'   dplyr::bind_rows(tibble::tibble(seq_id = paste0('seq', c(1,3)),
#'                                   feat_id = 'feat2',
#'                                   start = c(2, 4),
#'                                   end = c(6, 7)))
#' feats
#'
#' # Plot
#' plot_aln_annots(aln = aln, feats = feats)
plot_aln_annots <- function(aln, feats, region = NULL){

  # Convert seqinr aln to tibble and optain original sequence positions
  aln <- aln %>%
    purrr::map_dfr(seq_pos_from_aln,
                   .id = "seq_id")

  # Match features to alignment
  dat <- feats_to_aln(aln = aln, feats = feats)

  # Determine region to plot
  if(is.null(region)){
    region <- range(dat$aln_pos)
  }

  # Plot
  dat %>%
    dplyr::filter(aln_pos >= region[1] & aln_pos <= region[2]) %>%

    ggplot2::ggplot(ggplot2::aes(x = aln_pos, y = seq_id)) +
    ggplot2::geom_tile(ggplot2::aes(fill = feat_id), col = NA) +
    # ggplot2::geom_text(ggplot2::aes(label = char)) +
    ggfittext::geom_fit_text(ggplot2::aes(label = char)) +
    # scale_fill_manual(values = c("pink", NA)) +
    ggplot2::scale_fill_brewer() +
    ggplot2::ylab(label = "Sequence") +
    ggplot2::xlab(label = "Position") +
    ggplot2::theme(plot.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
}
