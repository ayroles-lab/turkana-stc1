
library(tidyverse)

basic_theme <-
    theme_minimal() +
    theme(
	text=element_text(size=14),
	panel.grid.minor=element_blank()
    )