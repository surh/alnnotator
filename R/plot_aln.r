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
#'
#' @return
feats_to_aln <- function(aln, feats){

  res <- aln %>%
    dplyr::left_join(feats, by = "seq_id") %>%
    dplyr::mutate(feat_id = replace(feat_id,
                             is.na(seq_pos) | seq_pos < start | seq_pos > end,
                             NA)) %>%
    dplyr::select(-start, -end)

  if(nrow(res) > nrow(aln))
    stop("ERROR: Some of your features seem to overlap within the same sequence.", call. = TRUE)

  return(res)
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
#'
#' @examples
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
    ggplot2::theme(plot.background = element_blank(),
                   panel.background = element_blank(),
                   panel.grid = element_blank())



}
