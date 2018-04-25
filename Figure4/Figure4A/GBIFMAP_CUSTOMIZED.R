library(rgbif)
library(ggplot2)

check_for_a_pkg <- function(x) {
  if (!requireNamespace(x, quietly = TRUE)) {
    stop("Please install ", x, call. = FALSE)
  } else {
    invisible(TRUE)
  }
}
gbif_capwords <- function(s, strict = FALSE, onlyfirst = FALSE) {
  cap <- function(s) paste(toupper(substring(s,1,1)),
                           {s <- substring(s,2); if(strict) tolower(s) else s}, sep = "", collapse = " " )
  if(!onlyfirst){
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
  } else {
    sapply(s, function(x)
      paste(toupper(substring(x,1,1)),
            tolower(substring(x,2)),
            sep="", collapse=" "), USE.NAMES=F)
  }
}

gbifmap_custom <- function(input = NULL, mapdatabase = "world", region = ".",
                           geom = geom_point, jitter = NULL, customize = NULL) {
  
  check_for_a_pkg("maps")
  if (!inherits(input, "data.frame")) stop("'input' must be a data.frame")
  tomap <- input[stats::complete.cases(input$decimalLatitude,
                                       input$decimalLatitude), ]
  tomap <- tomap[!tomap$decimalLongitude == 0 & !tomap$decimalLatitude == 0, ]
  tomap <- tomap[abs(tomap$decimalLatitude) <= 90 &
                   abs(tomap$decimalLongitude) <= 180, ]
  if (NROW(tomap) == 0) stop("after cleaning data, no records remaining")
  tomap$name <- as.factor(gbif_capwords(tomap$name, onlyfirst = TRUE))
  
  if (is.null(jitter)) {
    jitter <- position_jitter()
  }
  
  if (length(unique(tomap$name)) == 1) {
    theme2 <- theme(legend.position = "none")
  } else {
    theme2 <- NULL
  }
  
  world <- map_data(map = mapdatabase, region = region)
  message(paste("Rendering map...plotting ", nrow(tomap), " points", sep = ""))
  
  ggplot(world, aes(long, lat)) +
    geom_polygon(aes(group = group), fill = "white", color = "gray40",
                 size = 0.2) +
    geom(data = tomap, aes(decimalLongitude, decimalLatitude, colour = name),
         alpha = 0.6, size = 0.8, shape=16, position = jitter) +
    scale_color_brewer("", type = "qual", palette = 6) +
    labs(x = "", y = "") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", legend.key = element_blank()) +
    guides(col = guide_legend(nrow = 2)) +
    coord_fixed(ratio = 1) +
    blanktheme() +
    theme2 +
    customize
}
