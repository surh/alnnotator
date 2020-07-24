#' Read InterPro 5 TSV output file
#'
#' Variable row lengths cause warnings that can be ignored.
#'
#' @param file TSV file from InterPro 5
#' @param date_format format for parsing date
#'
#' @return A tibble
#' @export
read_ipro_tsv <- function(file,
                          date_format = "%d-%m-%Y"){
  
  dat <- readr::read_lines(file, n_max = 20) %>%
    purrr::map(stringr::str_split, pattern = "\t", simplify = TRUE) %>%
    purrr::map(function(x){ x[ x == "-"] <- NA; x[ x == ""] <- NA; x}) %>%
    purrr::map_dfr(tibble::as_tibble)
  col_names <- c("prot_id",
                 "seq_md5",
                 "seq_length",
                 "analysis",
                 "sig_accession",
                 "sig_description",
                 "start",
                 "end",
                 "score",
                 "status",
                 "date",
                 "ipro_accession",
                 "ipro_description",
                 "GO",
                 "pathway")
  colnames(dat) <-col_names[1:ncol(dat)]
  
  dat$seq_length <- readr::parse_number(dat$seq_length)
  dat$start <- readr::parse_number(dat$start)
  dat$end <- readr::parse_number(dat$end)
  dat$score <- readr::parse_number(dat$score)
  dat$status <- readr::parse_logical(dat$status)
  dat$date <- readr::parse_date(dat$date, format = date_format)

  dat
}
