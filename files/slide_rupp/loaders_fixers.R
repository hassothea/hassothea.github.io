

pacman::p_load(tidyverse)
pacman::p_load(readr)
pacman::p_load(latex2exp)
pacman::p_load(kableExtra)
pacman::p_load(gt)
# old_theme <-theme_set(theme_xaringan())
#old_theme <-theme_set(theme_minimal(base_size=14, base_family = "Helvetica"))

pacman::p_load("knitr")

opts_chunk$set(warning = FALSE, 
               message = FALSE, 
               cache = TRUE, 
               autodep = TRUE, 
               tidy = FALSE, 
               fig.retina = 4)

require(data.table, quietly = TRUE, warn.conflicts = FALSE)

pacman::p_load(readr)
