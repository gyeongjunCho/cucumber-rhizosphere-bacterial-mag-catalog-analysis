# =========================================================================
# Figure 2. ANI95 catalogue structure and relationships
# =========================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexUpset)
library(ggpubr)
library(grid)


# -------------------------------------------------------------------------
# Input and shared settings
# -------------------------------------------------------------------------

mag <- read.delim(
  "2-02_mag_metadata.tsv",
  stringsAsFactors = FALSE
)

ani95_rep <- read.delim(
  "4-03_ani95_representative_metadata.tsv",
  stringsAsFactors = FALSE
)

ani95_cluster <- read.delim(
  "4-02_ani95_cluster_summary.tsv",
  stringsAsFactors = FALSE
)

ani95_membership <- read.delim(
  "4-01_ani95_cluster_membership.tsv",
  stringsAsFactors = FALSE
)

sample_order <- c(
  "Con", "CW1", "CW2", "GY1", "MP1", "MP2",
  "NH2", "NH3", "PH1", "PH2", "US1"
)

# Catalogue consistency checks
stopifnot(
  nrow(mag) == 6505,
  nrow(ani95_rep) == 4032,
  nrow(ani95_cluster) == 4032,

  !anyDuplicated(mag$mag_id),
  !anyDuplicated(ani95_rep$mag_id),
  !anyDuplicated(ani95_cluster$ani95_cluster_id),

  all(ani95_rep$mag_id %in% mag$mag_id),
  sum(ani95_cluster$cluster_size) == 6505
)


# =========================================================================
# Figure 2A. MAG catalogue subset relationships
# =========================================================================

# Define catalogue subsets
fig2a_dat <- mag %>%
  mutate(
    Primary = TRUE,

    Operational_HQ =
      completeness > 90 &
      contamination < 10,

    MIMAG_HQ =
      mimag_high_quality == 1,

    GenBank =
      !is.na(ncbi_genome_accession) &
      trimws(ncbi_genome_accession) != "",

    ANI95_rep =
      mag_id %in% ani95_rep$mag_id
  )

# Confirm expected set sizes
set_summary <- fig2a_dat %>%
  summarise(
    Primary = sum(Primary),
    Operational_HQ = sum(Operational_HQ),
    MIMAG_HQ = sum(MIMAG_HQ),
    GenBank = sum(GenBank),
    ANI95_rep = sum(ANI95_rep)
  )

print(set_summary)

stopifnot(
  set_summary$Primary == 6505,
  set_summary$Operational_HQ == 4049,
  set_summary$MIMAG_HQ == 3659,
  set_summary$GenBank == 4031,
  set_summary$ANI95_rep == 4032
)

# UpSet plot
p2a <- upset(
  fig2a_dat,

  intersect = c(
    "Primary",
    "Operational_HQ",
    "MIMAG_HQ",
    "GenBank",
    "ANI95_rep"
  ),

  wrap = TRUE,

  labeller = function(x) {
    dplyr::recode(
      x,
      Primary = "Primary catalogue",
      Operational_HQ = "Operational HQ",
      MIMAG_HQ = "MIMAG HQ",
      GenBank = "GenBank accessioned",
      ANI95_rep = "ANI95 representative"
    )
  },

  sort_intersections_by = c(
    "degree",
    "cardinality"
  ),
  sort_intersections = "descending",

  base_annotations = list(
    "Intersection size" =
      intersection_size(
        counts = TRUE,
        text = list(size = 3.5)
      )
  ),

  set_sizes = (
    upset_set_size(
      geom = geom_bar(width = 0.65),
      position = "right"
    ) +
      geom_text(
        aes(label = after_stat(count)),
        stat = "count",
        color = "white",
        hjust = 1.15,
        size = 3.5
      ) +
      labs(
        x = NULL,
        y = "Set size"
      )
  ),

  width_ratio = 0.25,
  height_ratio = 0.55,
  name = "Catalogue subset"
) +
  ggtitle(
    "MAG catalogue subset relationships"
  ) +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust=0.5
    )
  )

p2a


# =========================================================================
# Figure 2B. ANI95 cluster-size composition
# =========================================================================

# Catalogue-level ANI95 cluster statistics
cluster_stats <- ani95_cluster %>%
  summarise(
    n_clusters = n(),
    n_primary_MAGs = sum(cluster_size),
    stand_alone_clusters = sum(cluster_size == 1),
    stand_alone_percent = 100 * mean(cluster_size == 1),
    multi_member_clusters = sum(cluster_size > 1),
    median_cluster_size = median(cluster_size),
    mean_cluster_size = mean(cluster_size),
    max_cluster_size = max(cluster_size)
  )

print(cluster_stats)

# Group clusters according to the number of member MAGs
cluster_size_class <- ani95_cluster %>%
  mutate(
    size_class = case_when(
      cluster_size == 1  ~ "Stand-alone",
      cluster_size == 2  ~ "2-member",
      cluster_size <= 5  ~ "3–5-member",
      cluster_size <= 10 ~ "6–10-member",
      TRUE               ~ ">10-member"
    ),

    size_class = factor(
      size_class,
      levels = c(
        "Stand-alone",
        "2-member",
        "3–5-member",
        "6–10-member",
        ">10-member"
      )
    )
  ) %>%
  count(
    size_class,
    name = "n_clusters"
  ) %>%
  mutate(
    percent =
      100 * n_clusters / sum(n_clusters),

    label =
      sprintf(
        "%s (%.1f%%)",
        n_clusters,
        percent
      )
  )

print(cluster_size_class)

p2b <- ggplot(
  cluster_size_class,
  aes(
    x = n_clusters,
    y = size_class
  )
) +

  geom_col(
    width = 0.7
  ) +

  geom_text(
    aes(label = label),
    hjust = -0.1,
    size = 3.5
  ) +

  scale_x_continuous(
    expand = expansion(
      mult = c(0, 0.20)
    ),
    limits = c(0, 3200)
  ) +

  labs(
    x = "Number of ANI95 clusters",
    y = "MAGs per ANI95 cluster",
    title = "ANI95 cluster-size composition"
  ) +

  theme_classic(
    base_size = 10
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0.5
    )
  )

p2b


# =========================================================================
# Figure 2C. Sample occupancy of ANI95 clusters
# =========================================================================

# Number of ANI95 clusters observed across 1–11 samples
cluster_occupancy <- ani95_cluster %>%
  count(
    n_samples,
    name = "n_clusters"
  ) %>%
  complete(
    n_samples = 1:11,
    fill = list(
      n_clusters = 0
    )
  ) %>%
  arrange(n_samples) %>%
  mutate(
    percent =
      100 * n_clusters / sum(n_clusters),

    label =
      sprintf(
        "%s\n(%.1f%%)",
        n_clusters,
        percent
      )
  )

print(cluster_occupancy)

stopifnot(
  sum(cluster_occupancy$n_clusters) == 4032,
  all(cluster_occupancy$n_samples %in% 1:11)
)

p2c <- ggplot(
  cluster_occupancy,
  aes(
    x = n_samples,
    y = n_clusters
  )
) +

  geom_col(
    width = 0.7
  ) +

  geom_text(
    aes(label = label),
    vjust = -0.25,
    size = 3.2,
    lineheight = 0.9
  ) +

  scale_x_continuous(
    breaks = 1:11
  ) +

  scale_y_continuous(
    limits = c(0,3500),
    expand = expansion(
      mult = c(0, 0.12)
    )
  ) +

  labs(
    x = "Number of samples per ANI95 cluster",
    y = "Number of ANI95 clusters",
    title = "Sample occupancy of ANI95 clusters"
  ) +

  theme_classic(
    base_size = 10
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0.5
    )
  )

p2c


# =========================================================================
# Figure 2D. ANI95 cluster sharing among samples
# =========================================================================

# One record per ANI95 cluster and sample
cluster_presence <- ani95_membership %>%
  distinct(
    ani95_cluster_id,
    sample_id
  )

# Pairwise cluster sharing.
# Diagonal values represent the number of distinct ANI95 clusters
# recovered from each individual sample.
pairwise_sharing <- cluster_presence %>%
  inner_join(
    cluster_presence,
    by = "ani95_cluster_id",
    suffix = c("_x", "_y")
  ) %>%

  count(
    sample_id_x,
    sample_id_y,
    name = "n_shared"
  ) %>%

  complete(
    sample_id_x = sample_order,
    sample_id_y = sample_order,
    fill = list(
      n_shared = 0
    )
  ) %>%

  mutate(
    x_id = match(
      sample_id_x,
      sample_order
    ),
    y_id = match(
      sample_id_y,
      sample_order
    )
  ) %>%

  # Retain one half of the symmetric matrix
  filter(
    y_id >= x_id
  ) %>%

  mutate(
    sample_id_x = factor(
      sample_id_x,
      levels = sample_order
    ),

    sample_id_y = factor(
      sample_id_y,
      levels = rev(sample_order)
    )
  )

stopifnot(
  all(
    sample_order %in%
      cluster_presence$sample_id
  ),
  n_distinct(
    cluster_presence$ani95_cluster_id
  ) == 4032
)

p2d <- ggplot(
  pairwise_sharing,
  aes(
    x = sample_id_x,
    y = sample_id_y,
    fill = n_shared
  )
) +

  geom_tile(
    linewidth = 0.4,
    color = "white"
  ) +

  geom_text(
    aes(label = n_shared),
    size = 3
  ) +

  scale_fill_viridis_c(
    option = "C",
    direction = -1,
    name = "Shared ANI95\nclusters"
  ) +

  labs(
    x = NULL,
    y = NULL,
    title = "ANI95 cluster sharing among samples"
  ) +

  coord_fixed() +

  theme_classic(
    base_size = 10
  ) +

  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),

    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0.5
    ),

    legend.title = element_text(
      size = 9
    )
  )

p2d


# =========================================================================
# Assemble Figure 2
# =========================================================================

# ComplexUpset returns a wrapped patchwork object.
# Capture it as a grob before combining with ggpubr.
p2a_grob <- grid::grid.grabExpr(
  print(p2a)
)

p2a_gg <- ggpubr::as_ggplot(
  p2a_grob
)

# Lower-left panels
fig2_left <- ggarrange(
  p2b,
  p2c,
  ncol = 1,
  nrow = 2,
  labels = c("b", "c"),
  font.label = list(
    size = 20,
    face = "bold"
  ),
  align = "hv"
)

# Lower section
fig2_bottom <- ggarrange(
  fig2_left,
  p2d,
  ncol = 2,
  widths = c(1, 1.1),
  labels = c("", "d"),
  font.label = list(
    size = 20,
    face = "bold"
  )
)

# Final Figure 2
fig2 <- ggarrange(
  p2a_gg,
  fig2_bottom,
  ncol = 1,
  nrow = 2,
  labels = c("a", ""),
  font.label = list(
    size = 20,
    face = "bold"
  ),
  heights = c(1, 1.25)
)

fig2


# =========================================================================
# Export Figure 2
# =========================================================================

ggsave(
  filename = "figure_and_table/Figure2.pdf",
  plot = fig2,
  width = 12,
  height = 10,
  units = "in",
  device = cairo_pdf
)
