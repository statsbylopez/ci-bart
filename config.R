# local configuration hack

if (Sys.info()["user"] == "liangyuanhu") {
  # For Liangyuan
  root <- "/Users/liangyuanhu/GoogleDrive/ci_bart/"
}

if (Sys.info()["user"] == "mlopez1") {
  # For Mike
  root <- "~/Dropbox/ci-bart/"
}

if (Sys.info()["user"] == "guchenyang") {
  # For Chenyang
  root <- ""
}

library(knitr)
library(tidyverse)
library(splines)
library(reshape2)