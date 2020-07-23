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
  # file <- "interpro_examples/iprscan5-R20200722-215734-0731-39996086-p2m.txt"
  # file <- "interpro_examples/iprscan5-R20200722-215836-0101-63677525-p2m.txt"
  # date_format <- "%d-%m-%Y"
  oldwarn <- getOption("warn")
  options(warn = -1)
  dat <- readr::read_tsv(file, col_names = FALSE,
                         na = c("", "NA", "-"),
                         col_types = readr::cols(X1 = readr::col_character(),
                                                 X2 = readr::col_character(),
                                                 X3 = readr::col_number(),
                                                 X4 = readr::col_character(),
                                                 X5 = readr::col_character(),
                                                 X6 = readr::col_character(),
                                                 X7 = readr::col_number(),
                                                 X8 = readr::col_number(),
                                                 X9 = readr::col_number(),
                                                 X10 = readr::col_logical(),
                                                 X11 = readr::col_date(format = date_format),
                                                 .default = readr::col_character()))
  options(warn = oldwarn)

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

  dat
}
