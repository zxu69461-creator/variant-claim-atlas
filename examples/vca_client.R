#!/usr/bin/env Rscript
# Client for the Variant–Trait Claim Atlas (https://zxu69461-creator.github.io/variant-claim-atlas/).
# Reads the static data files of the web tool and reproduces its verdict. Needs the package jsonlite.
#
#   Rscript vca_client.R check rs7903146 "type 2 diabetes"     # verdict for one variant-trait claim (names or abbreviations: T2D)
#   Rscript vca_client.R claims rs1800795 claims.csv           # every claim on a variant
#   Rscript vca_client.R batch my_list.txt results.csv         # rsIDs, or "rsID, trait" pairs, one per line
#
# Set VCA_BASE to a local copy of the repository to read it instead of the web site.
# As functions: source("vca_client.R"); vca_check("rs7903146", "T2D")

suppressPackageStartupMessages(library(jsonlite))
VCA_BASE <- Sys.getenv("VCA_BASE", "https://zxu69461-creator.github.io/variant-claim-atlas/")
.vca <- new.env()

vca_get <- function(path) {
  if (is.null(.vca[[path]])) {
    src <- if (grepl("^http", VCA_BASE)) paste0(sub("/?$", "/", VCA_BASE), path) else file.path(VCA_BASE, path)
    t <- fromJSON(src, simplifyVector = FALSE)
    .vca[[path]] <- t
  }
  .vca[[path]]
}
vca_rs <- function(x) as.integer(gsub("\\D", "", x))

# Claim records of one variant as a data frame.
vca_claims <- function(rs) {
  rs <- vca_rs(rs); meta <- vca_get("data/meta.json")
  t <- vca_get(sprintf("data/claims/%02d.json", rs %% meta$shards))
  rows <- Filter(function(r) r[[2]] == rs, t$rows)
  if (!length(rows)) return(data.frame())
  col <- function(i) sapply(rows, function(r) if (is.null(r[[i]])) NA else r[[i]])
  d <- as.data.frame(setNames(lapply(seq_along(t$cols), col), unlist(t$cols)), stringsAsFactors = FALSE)
  d
}

# Traits named by a query: names that contain it, and traits for which abstracts used it as a phrase (as on the web site).
vca_traits_for <- function(q) {
  q <- tolower(trimws(q))
  if (is.null(.vca$traits)) {
    .vca$traits <- sapply(vca_get("data/traits.json")$rows, `[[`, 1)
    lex <- vca_get("data/lexicon.json")$rows
    .vca$phr <- data.frame(trait = rep(sapply(lex, `[[`, 1), sapply(lex, function(r) length(r[[3]]))),
                           phrase = unlist(lapply(lex, `[[`, 3)), stringsAsFactors = FALSE)
  }
  # names match at the start of a word; short queries (abbreviations) only as a whole word
  rx <- paste0("(?<![a-z0-9])\\Q", q, "\\E", if (nchar(q) <= 4) "(?![a-z0-9])" else "")
  unique(c(.vca$traits[grepl(rx, tolower(.vca$traits), perl = TRUE)], .vca$phr$trait[.vca$phr$phrase == q]))
}

# Verdict for the claims on variant rs that concern the trait (all claims when trait is NULL).
vca_check <- function(rs, trait = NULL) {
  d <- vca_claims(rs)
  if (nrow(d) && !is.null(trait) && nzchar(trait)) d <- d[d$trait %in% vca_traits_for(trait), , drop = FALSE]
  rep <- d[d$status == "reported", , drop = FALSE]; tst <- rep[rep$testable, , drop = FALSE]; sup <- tst[tst$supported, , drop = FALSE]
  v <- if (!nrow(d)) "no claim on this trait" else if (!nrow(rep)) "no reported association" else if (!nrow(tst)) "not testable" else
    if (nrow(sup) / nrow(tst) >= 0.5) "supported" else if (nrow(sup)) "partly supported" else "not supported"
  yrs <- sup$first_gw_year[!is.na(sup$first_gw_year)]
  data.frame(rsid = paste0("rs", vca_rs(rs)), trait_query = if (is.null(trait)) "" else trait,
             matched_traits = paste(sort(unique(d$trait)), collapse = "; "), verdict = v, reported = nrow(rep),
             absent = nrow(d) - nrow(rep), testable = nrow(tst), supported = nrow(sup),
             first_genome_wide_year = if (length(yrs)) min(yrs) else NA, stringsAsFactors = FALSE)
}

vca_batch <- function(lines) {
  out <- list()
  for (l in trimws(lines)) {
    m <- regmatches(l, regexec("^(rs\\s?[0-9]+)\\s*[,;:\t]\\s*(\\S.*)$", l, ignore.case = TRUE))[[1]]
    if (length(m)) out[[length(out) + 1]] <- vca_check(m[2], gsub("^[\"' ]+|[\"' ]+$", "", m[3]))
    else for (r in regmatches(l, gregexpr("rs\\s?[0-9]+", l, ignore.case = TRUE))[[1]]) out[[length(out) + 1]] <- vca_check(r)
  }
  do.call(rbind, out)
}

if (sys.nframe() == 0) {
  a <- commandArgs(trailingOnly = TRUE)
  if (!length(a)) { cat("usage: Rscript vca_client.R check|claims|batch ...\n"); quit(status = 1) }
  res <- switch(a[1],
    check = vca_check(a[2], if (length(a) > 2) a[3] else NULL),
    claims = vca_claims(a[2]),
    batch = vca_batch(readLines(a[2], warn = FALSE)),
    stop("unknown command ", a[1]))
  dest <- if (a[1] %in% c("claims", "batch") && length(a) > 2) a[3] else ""
  if (nzchar(dest)) { write.csv(res, dest, row.names = FALSE); message(nrow(res), " rows written to ", dest) } else print(res, row.names = FALSE)
}
