
# This replace "SingleImagePlot" in seurat to make faster ImagePlots
message("Replace SingleImagePlot in Seurat with scattermore to make faster ImagePlots")
require(scattermore)

# Replace SingleImagePlot
# Date: 2024-05-01
# From: https://github.com/satijalab/seurat/blob/c54e57d3423b3f711ccd463e14965cc8de86c31b/R/visualization.R#L2362
SingleImagePlot <- function(
  data,
  col.by = NA,
  col.factor = TRUE,
  cols = NULL,
  shuffle.cols = FALSE,
  size = 0.1,
  molecules = NULL,
  mols.size = 0.1,
  mols.cols = NULL,
  mols.alpha = 1.0,
  alpha = molecules %iff% 0.3 %||% 0.6,
  border.color = 'white',
  border.size = NULL,
  na.value = 'grey50',
  dark.background = TRUE,
  ...
) {
	message("Replaced geom_point with geom_scattermore")
  # Check input data
  if (!rlang::is_na(x = data)) {
    if (!all(c('x', 'y', 'cell', 'boundary') %in% colnames(x = data))) {
      stop("Invalid data coordinates")
    }
    if (!rlang::is_na(x = col.by)) {
      if (!col.by  %in% colnames(x = data)) {
        warning(
          "Cannot find 'col.by' ('",
          col.by,
          "') in data coordinates",
          immediate. = TRUE
        )
        col.by <- NA
      } else if (isTRUE(x = col.factor) && !is.factor(x = data[[col.by]])) {
        data[[col.by]] <- factor(
          x = data[[col.by]],
          levels = unique(x = data[[col.by]])
        )
      } else if (isFALSE(x = col.factor) && is.factor(x = data[[col.by]])) {
        data[[col.by]] <- as.vector(x = data[[col.by]])
      }
    }
    if (rlang::is_na(x = col.by) && !is.null(x = cols)) {
      col.by <- RandomName(length = 7L)
      data[[col.by]] <- TRUE
    }
    if (!is.factor(x = data$boundary)) {
      data$boundary <- factor(
        x = data$boundary,
        levels = unique(x = data$boundary)
      )
    }
    # Determine alphas
    if (is.na(x = alpha)) {
      alpha <- 1L
    }
    alpha.min <- ifelse(
      test = alpha < 1L,
      yes = 1 * (10 ^ .FindE(x = alpha)),
      no = 0.1
    )
    if (alpha.min == alpha) {
      alpha.min <- 1 * (10 ^ (.FindE(x = alpha) - 1))
    }
    alphas <- .Cut(
      min = alpha.min,
      max = alpha,
      n = length(x = levels(x = data$boundary))
    )
  }
  # Assemble plot
  plot <- ggplot(
    data = data %NA% NULL,
    mapping = aes_string(
      x = 'y',
      y = 'x',
      alpha = 'boundary',
      fill = col.by %NA% NULL
    )
  )
  if (!rlang::is_na(x = data)) {
    plot <- plot +
      scale_alpha_manual(values = alphas) +
      if (anyDuplicated(x = data$cell)) {
        if (is.null(border.size)) {
          border.size <- 0.3
        }
        geom_polygon(
          mapping = aes_string(group = 'cell'),
          color = border.color,
          size = border.size
        )
      } else {
        # Default to no borders when plotting centroids
        if (is.null(border.size)) {
          border.size <- 0.0
        }
		geom_scattermore(
			shape = 21,
			color = border.color,
			stroke = border.size,
			size = size
		)
      }
    if (!is.null(x = cols)) {
      plot <- plot + if (is.numeric(x = cols) || cols[1L] %in% rownames(x = brewer.pal.info)) {
        palette <- brewer.pal(n = length(x = levels(x = data[[col.by]])), cols)
        if (length(palette) < length(x = levels(x = data[[col.by]]))) {
          num.blank <- length(x = levels(x = data[[col.by]])) - length(palette)
          palette <- c(palette, rep(na.value, num.blank))
        }
        if (isTRUE(shuffle.cols)) {
          palette <- sample(palette)
        }
        scale_fill_manual(values = palette, na.value = na.value)
      } else if (cols[1] %in% DiscretePalette(n = NULL)) {
        scale_fill_manual(
          values = DiscretePalette(
            n = length(x = levels(x = data[[col.by]])),
            palette = cols,
            shuffle = shuffle.cols
          ),
          na.value = na.value
        )
      } else {
        if (isTRUE(shuffle.cols)) {
          cols <- sample(cols)
        }
        scale_fill_manual(values = cols, na.value = na.value)
      }
    }
    if (length(x = levels(x = data$boundary)) == 1L) {
      plot <- plot + guides(alpha = 'none')
    }
    # Adjust guides
    if (isTRUE(x = col.factor) && length(x = levels(x = data[[col.by]])) <= 1L) {
      plot <- plot + guides(fill = 'none')
    }
  }
  # Add molecule sets
  if (is.data.frame(x = molecules)) {
    if (all(c('x', 'y', 'molecule') %in% colnames(x = molecules))) {
      if (!is.factor(x = molecules$molecule)) {
        molecules$molecule <- factor(
          x = molecules$molecule,
          levels = unique(x = molecules$molecule)
        )
      }
      plot <- plot + geom_scattermore(
        mapping = aes_string(fill = NULL, alpha = NULL, color = "molecule"),
        data = molecules,
        size = mols.size,
        alpha = mols.alpha,
        show.legend = c(color = TRUE, fill = FALSE, alpha = FALSE)
      ) +
        guides(color = guide_legend(override.aes = list(size = 3L))) +
        scale_color_manual(
          name = 'Molecules',
          values = mols.cols %||% brewer.pal(
            n = length(x = levels(x = molecules$molecule)),
            name = "Set1"
          ),
          guide = guide_legend()
        )
    } else {
      warning("Invalid molecule coordinates", immediate. = TRUE)
    }
  }

  if (isTRUE(dark.background)) {
    plot <- plot + DarkTheme()
  }
  return(plot)
}

# Replace
assignInNamespace(x = "SingleImagePlot", value = SingleImagePlot, ns = "Seurat")