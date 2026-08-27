# =========================================================================
# Figure 1A. Coastal source-soil sampling sites
# =========================================================================

library(dplyr)
library(ggplot2)
library(ggrepel)
library(sf)
library(rnaturalearth)

meta <- read.delim(
  "2-01_sample_metadata.tsv",
  stringsAsFactors = FALSE
)

map_df <- meta %>%
  filter(
    !is.na(source_soil_latitude),
    !is.na(source_soil_longitude)
  ) %>%
  transmute(
    sample_id,
    source_soil_site,
    latitude = source_soil_latitude,
    longitude = source_soil_longitude
  )

stopifnot(nrow(map_df) == 10)

# High-resolution Natural Earth land boundary (1:10m)
land <- ne_download(
  scale = 10,
  type = "land",
  category = "physical",
  returnclass = "sf"
)

p1a <- ggplot() +
  geom_sf(
    data = land,
    fill = "grey96",
    colour = "grey25",
    linewidth = 0.35
  ) +

  annotate(
    "text",
    x = 127.6,
    y = 36.0,
    label = "South Korea",
    size = 6,
    fontface = "bold",
    colour = "grey65"
  ) +

  geom_point(
    data = map_df,
    aes(
      x = longitude,
      y = latitude
    ),
    shape = 21,
    size = 3,
    fill = "white",
    colour = "black",
    stroke = 0.8
  ) +

  geom_text_repel(
    data = map_df,
    aes(
      x = longitude,
      y = latitude,
      label = sample_id
    ),
    size = 3.2,
    fontface = "bold",
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +

  coord_sf(
    xlim = c(125.7, 130.0),
    ylim = c(34.3, 36.5),
    expand = FALSE
  ) +

  labs(
    x = NULL,
    y = NULL,
    title = "Coastal source-soil sampling sites\n"
  ) +

  theme_void() +

  theme(
    plot.title = element_text(
      size = 11,
      face = "bold",
      hjust = 0.5
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.4
    )
  )

p1a


# =========================================================================
# Figure 1B. Experimental design and rhizosphere sampling
# =========================================================================

library(ggplot2)
library(grid)

# Global appearance settings
fig1b_font_size <- 3
fig1b_fill <- "#DAFBFB"

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

add_box <- function(
    p,
    xmin, xmax,
    ymin, ymax,
    label,
    fill = fig1b_fill,
    border = "grey30",
    size = fig1b_font_size,
    fontface = "plain",
    lineheight = 1.10
) {

  p +
    annotate(
      "rect",
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill,
      color = border,
      linewidth = 0.6
    ) +
    annotate(
      "text",
      x = (xmin + xmax) / 2,
      y = (ymin + ymax) / 2,
      label = label,
      size = size,
      fontface = fontface,
      lineheight = lineheight
    )
}

add_arrow <- function(
    p,
    x, xend,
    y, yend,
    linewidth = 0.6
) {

  p +
    annotate(
      "segment",
      x = x,
      xend = xend,
      y = y,
      yend = yend,
      linewidth = linewidth,
      lineend = "round",
      arrow = arrow(
        length = unit(0.17, "cm"),
        type = "closed"
      )
    )
}

add_line <- function(
    p,
    x, xend,
    y, yend,
    linewidth = 0.6
) {

  p +
    annotate(
      "segment",
      x = x,
      xend = xend,
      y = y,
      yend = yend,
      linewidth = linewidth,
      lineend = "round"
    )
}

# -------------------------------------------------------------------------
# Base canvas
# -------------------------------------------------------------------------

p1b <- ggplot() +
  coord_cartesian(
    xlim = c(0, 12),
    ylim = c(0, 18),
    clip = "off"
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 10,
      l = 15
    )
  )

# -------------------------------------------------------------------------
# Source-soil treatments
# -------------------------------------------------------------------------

p1b <- add_box(
  p1b,
  xmin = 0.7,
  xmax = 6.5,
  ymin = 14.8,
  ymax = 17.3,
  label = paste0(
    "10 coastal source-soil treatments\n",
    "Source soil + sterilized horticultural substrate\n",
    "1:9 (w/w)"
  )
)

# -------------------------------------------------------------------------
# Control
# -------------------------------------------------------------------------

p1b <- add_box(
  p1b,
  xmin = 7.5,
  xmax = 11.3,
  ymin = 14.8,
  ymax = 17.3,
  label = paste0(
    "Control\n",
    "Sterilized horticultural\n",
    "substrate only"
  )
)

# -------------------------------------------------------------------------
# Merge source-soil treatments and control
# -------------------------------------------------------------------------

p1b <- add_line(
  p1b,
  x = 3.6,
  xend = 3.6,
  y = 14.8,
  yend = 13.8
)

p1b <- add_line(
  p1b,
  x = 9.4,
  xend = 9.4,
  y = 14.8,
  yend = 13.8
)

p1b <- add_line(
  p1b,
  x = 3.6,
  xend = 9.4,
  y = 13.8,
  yend = 13.8
)

p1b <- add_arrow(
  p1b,
  x = 6.5,
  xend = 6.5,
  y = 13.8,
  yend = 13.0
)

# -------------------------------------------------------------------------
# Eleven experimental treatments
# -------------------------------------------------------------------------

p1b <- add_box(
  p1b,
  xmin = 4.1,
  xmax = 8.9,
  ymin = 11.2,
  ymax = 12.9,
  label = "11 experimental treatments",
  size = fig1b_font_size * 1.05,
  fontface = "bold"
)

# -------------------------------------------------------------------------
# Cucumber cultivation
# -------------------------------------------------------------------------

p1b <- add_arrow(
  p1b,
  x = 6.5,
  xend = 6.5,
  y = 11.2,
  yend = 10.4
)

p1b <- add_box(
  p1b,
  xmin = 3,
  xmax = 10,
  ymin = 7.8,
  ymax = 10.3,
  label = paste0(
    "Cucumber cv. Baekdadagi\n",
    "30 plants per treatment | 20 days\n",
    "Sown 14 May 2024 → rhizosphere collected 3 June 2024"
  )
)

# -------------------------------------------------------------------------
# Rhizosphere recovery
# -------------------------------------------------------------------------

p1b <- add_arrow(
  p1b,
  x = 6.5,
  xend = 6.5,
  y = 7.8,
  yend = 7.0
)

p1b <- add_box(
  p1b,
  xmin = 2.5,
  xmax = 10.5,
  ymin = 3.7,
  ymax = 6.9,
  label = paste0(
    "Rhizosphere recovery\n",
    "14 perimeter plants excluded | 8 roots × 2 batches\n",
    "500 mL PBS per batch | 20 min sonication in an ice-water bath\n",
    "Batches combined → centrifugation at 3,000 × g for 20 min"
  )
)

# -------------------------------------------------------------------------
# Final composite samples
# -------------------------------------------------------------------------

p1b <- add_arrow(
  p1b,
  x = 6.5,
  xend = 6.5,
  y = 3.7,
  yend = 2.9
)

p1b <- add_box(
  p1b,
  xmin = 2.5,
  xmax = 10.5,
  ymin = 0.5,
  ymax = 2.8,
  label = paste0(
    "One pooled composite rhizosphere sample per treatment\n",
    "11 composite samples total"
  ),
  size = fig1b_font_size * 1.05,
  fontface = "bold"
)

p1b <- p1b +
  ggtitle("Experimental design and rhizosphere sampling") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0.5,
      margin = margin(b = 8)
    )
  )

p1b


# =========================================================================
# Figure 1C. Genome-resolved metagenomic workflow
# =========================================================================

library(magick)
library(grid)
library(ggpubr)

fig1c_png <- png::readPNG(
  "fig1_asset/Fig1C_workflow.png"
)

p1c <- ggpubr::as_ggplot(
  grid::rasterGrob(
    fig1c_png,
    interpolate = TRUE
  )
) + ggtitle("Genome-resolved metagenomic workflow") +
  theme(
    plot.title = element_text(
      size = 11,
      face = "bold",
      hjust = 0.5
    )
  )

p1c


# =========================================================================
# Assemble Figure 1
# =========================================================================

fig1_left <- ggarrange(
  p1a,
  p1b,
  ncol = 1,
  nrow = 2,
  heights = c(1,1.2),
  labels = c(NA, "b"),
  font.label = list(
    size = 14,
    face = "bold"
  )
)

fig1 <- ggarrange(
  fig1_left,
  p1c,
  ncol = 2,
  widths = c(1, 1.5),
  labels = c("a", "c"),
  font.label = list(
    size = 14,
    face = "bold"
  ),
  align = "hv"
)

fig1


# =========================================================================
# Export Figure 1
# =========================================================================

ggsave(
  filename = "figure_and_table/Figure1.pdf",
  plot = fig1,
  width = 16.2,
  height = 8.2,
  units = "in",
  device = cairo_pdf
)

