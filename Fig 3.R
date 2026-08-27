library(showtext)
library(sysfonts)

# -------------------------------------------------------------------------
# Font setup
# -------------------------------------------------------------------------

dejavu_regular <- system2(
  "fc-match",
  args = c(
    "-f",
    shQuote("%{file}"),
    shQuote("dejavu")
  ),
  stdout = TRUE
) |>
  trimws()

dejavu_bold <- system2(
  "fc-match",
  args = c(
    "-f",
    shQuote("%{file}"),
    shQuote("DejaVu Sans:style=Bold")
  ),
  stdout = TRUE
) |>
  trimws()

# Check detected font files
print(dejavu_regular)
print(dejavu_bold)

stopifnot(
  length(dejavu_regular) == 1,
  length(dejavu_bold) == 1,
  file.exists(dejavu_regular),
  file.exists(dejavu_bold)
)

sysfonts::font_add(
  family = "dejavu",
  regular = dejavu_regular,
  bold = dejavu_bold
)

showtext::showtext_auto()

# =========================================================================
# Figure 3. MAG quality and genome characteristics
# =========================================================================

library(dplyr)
library(ggplot2)
library(ggnewscale)
library(ggpubr)
library(VennDiagram)
library(grid)


# -------------------------------------------------------------------------
# Input and shared data preparation
# -------------------------------------------------------------------------

mag <- read.delim(
  "2-02_mag_metadata.tsv",
  stringsAsFactors = FALSE
)

sample_order <- c(
  "Con", "CW1", "CW2", "GY1", "MP1", "MP2",
  "NH2", "NH3", "PH1", "PH2", "US1"
)

mag <- mag %>%
  mutate(
    sample_id = factor(
      sample_id,
      levels = sample_order
    ),
    Operational_HQ =
      completeness > 90 &
      contamination < 10,
    GenBank =
      !is.na(ncbi_genome_accession) &
      trimws(ncbi_genome_accession) != "",
    MIMAG_HQ =
      mimag_high_quality == 1
  )

# Catalogue consistency checks
stopifnot(
  nrow(mag) == 6505,
  !anyDuplicated(mag$mag_id),
  sum(mag$Operational_HQ) == 4049,
  sum(mag$GenBank) == 4031,
  sum(mag$MIMAG_HQ) == 3659
)


# =========================================================================
# Figure 3A. MAG completeness and contamination by sample
# =========================================================================

hq_range <- data.frame(
  xmin = 90,
  xmax = 100,
  ymin = 0,
  ymax = 10,
  label = "Operational HQ range"
)

p3a <- ggplot(
  mag,
  aes(
    x = completeness,
    y = contamination
  )
) +

  # Operational-HQ range
  geom_rect(
    data = hq_range,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = label
    ),
    inherit.aes = FALSE,
    alpha = 0.15
  ) +

  scale_fill_manual(
    values = c(
      "Operational HQ range" = "gold"
    ),
    labels = c(
      "Operational HQ range" =
        "Operational HQ range\n(>90% completeness, <10% contamination)"
    ),
    name = NULL
  ) +

  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(alpha = 0.25)
    )
  ) +

  ggnewscale::new_scale_fill() +

  # Primary MAGs outside the operational-HQ range
  geom_point(
    data = filter(mag, !Operational_HQ),
    color = "grey70",
    shape = 16,
    size = 0.65,
    alpha = 0.40
  ) +

  # Operational-HQ MAGs
  geom_point(
    data = filter(mag, Operational_HQ),
    aes(
      shape = GenBank,
      fill = MIMAG_HQ
    ),
    color = "black",
    size = 1.9,
    stroke = 0.70,
    alpha = 0.90
  ) +

  # MIMAG genome-quality thresholds
  geom_vline(
    xintercept = 90,
    linetype = "dashed",
    linewidth = 0.4
  ) +

  geom_hline(
    yintercept = 5,
    linetype = "dashed",
    linewidth = 0.4
  ) +

  facet_wrap(
    ~ sample_id,
    ncol = 6
  ) +

  scale_shape_manual(
    values = c(
      `FALSE` = 21,
      `TRUE` = 24
    ),
    breaks = "TRUE",
    labels = "GenBank accessioned",
    name = NULL
  ) +

  scale_fill_manual(
    values = c(
      `FALSE` = "white",
      `TRUE` = "tomato"
    ),
    breaks = "TRUE",
    labels = paste0(
      "MIMAG HQ\n",
      "(>90% completeness, <5% contamination;\n",
      "5S/16S/23S rRNA, \u226518 tRNA aa types)"
    ),
    name = NULL
  ) +

  guides(
    shape = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(
        shape = 24,
        fill = "white",
        color = "black",
        size = 3,
        stroke = 0.9,
        alpha = 1
      )
    ),
    fill = guide_legend(
      order = 3,
      nrow = 1,
      override.aes = list(
        shape = 24,
        fill = "tomato",
        color = NA,
        size = 3,
        stroke = 0,
        alpha = 1
      )
    )
  ) +

  scale_x_continuous(
    limits = c(50, 100),
    breaks = c(50, 60, 70, 80, 90, 100)
  ) +

  scale_y_continuous(
    limits = c(0, 10),
    breaks = seq(0, 10, 2)
  ) +

  labs(
    x = "Completeness (%)",
    y = "Contamination (%)",
    title = "MAG completeness and contamination by sample"
  ) +

  theme_classic(base_size = 10) +

  theme(
    text = element_text(family = "dejavu"),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    panel.border = element_rect(
      fill = NA,
      linewidth = 0.4
    ),
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0.5
    ),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.justification = "center",
    legend.text = element_text(size = 8),
    legend.key = element_blank(),
    legend.spacing.x = unit(0.15, "cm")
  )


# =========================================================================
# Figure 3B. MAG contiguity by sample
# =========================================================================

contig_comp <- mag %>%
  mutate(
    contig_class = case_when(
      contig_count == 1  ~ "1",
      contig_count <= 5  ~ "2–5",
      contig_count <= 10 ~ "6–10",
      contig_count <= 25 ~ "11–25",
      TRUE               ~ "26–50"
    ),
    contig_class = factor(
      contig_class,
      levels = c(
        "1",
        "2–5",
        "6–10",
        "11–25",
        "26–50"
      )
    )
  ) %>%
  count(
    sample_id,
    contig_class,
    name = "n_MAGs"
  ) %>%
  group_by(sample_id) %>%
  mutate(
    percent = 100 * n_MAGs / sum(n_MAGs)
  ) %>%
  ungroup()

# Optional summary for console inspection
contig_summary <- mag %>%
  summarise(
    n_MAGs = n(),
    single_contig = sum(contig_count == 1),
    single_contig_percent = 100 * mean(contig_count == 1),
    median_contigs = median(contig_count),
    mean_contigs = mean(contig_count),
    max_contigs = max(contig_count)
  )

print(contig_summary)

p3b <- ggplot(
  contig_comp,
  aes(
    x = sample_id,
    y = percent,
    fill = contig_class
  )
) +

  geom_col(
    width = 0.75,
    color = "white",
    linewidth = 0.2
  ) +

  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0))
  ) +

  scale_fill_brewer(
    palette = "YlGnBu",
    direction = -1,
    name = "Contigs per MAG"
  ) +

  labs(
    x = NULL,
    y = "Primary MAGs (%)",
    title = "MAG contiguity by sample"
  ) +

  theme_classic(base_size = 10) +

  theme(
    text = element_text(family = "dejavu"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0.5
    ),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8)
  )


# =========================================================================
# Figure 3C. MIMAG high-quality criteria
# =========================================================================

mimag <- read.delim(
  "3-03_mimag_quality.tsv",
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(mimag) == 6505,
  !anyDuplicated(mimag$mag_id)
)

# Define MIMAG high-quality criterion sets
mimag_sets <- mimag %>%
  transmute(
    mag_id,
    genome_quality =
      pass_completeness_gt90 == 1 &
      pass_contamination_lt5 == 1,
    complete_rrna =
      pass_all_5S_16S_23S_rRNA == 1,
    trna_diversity =
      pass_trna_at_least_18_of_20_aa_types == 1,
    mimag_hq =
      mimag_high_quality == 1
  )

# Set sizes
n_genome <- sum(mimag_sets$genome_quality)
n_rrna   <- sum(mimag_sets$complete_rrna)
n_trna   <- sum(mimag_sets$trna_diversity)

# Pairwise intersections
n_genome_rrna <- sum(
  mimag_sets$genome_quality &
    mimag_sets$complete_rrna
)

n_genome_trna <- sum(
  mimag_sets$genome_quality &
    mimag_sets$trna_diversity
)

n_rrna_trna <- sum(
  mimag_sets$complete_rrna &
    mimag_sets$trna_diversity
)

# Three-way intersection
n_all <- sum(
  mimag_sets$genome_quality &
    mimag_sets$complete_rrna &
    mimag_sets$trna_diversity
)

stopifnot(
  n_all == 3659,
  all(
    mimag_sets$mimag_hq ==
      (
        mimag_sets$genome_quality &
          mimag_sets$complete_rrna &
          mimag_sets$trna_diversity
      )
  )
)

# Construct three-set Venn diagram
venn_glist <- VennDiagram::draw.triple.venn(
  area1 = n_genome,
  area2 = n_rrna,
  area3 = n_trna,

  n12 = n_genome_rrna,
  n13 = n_genome_trna,
  n23 = n_rrna_trna,
  n123 = n_all,

  category = c(
    "Genome quality\n>90% completeness\n<5% contamination",
    "Complete rRNA\n5S + 16S + 23S",
    "tRNA diversity\n≥18/20 aa types"
  ),

  fill = c(
    "gold",
    "skyblue",
    "tomato"
  ),
  alpha = rep(0.30, 3),

  col = rep("grey35", 3),
  lwd = 1,

  label.col = "black",
  cex = 0.9,
  fontface = "plain",
  fontfamily = "dejavu",

  cat.col = "black",
  cat.cex = 0.8,
  cat.fontfamily = "dejavu",
  cat.dist = c(
    0.25,  # Genome quality
    0.15,  # Complete rRNA
    0.12   # tRNA diversity
  ),
  cat.pos = c(
    -60,
    30,
    180
  ),

  rotation.degree = 20,
  rotation.centre = c(0.5, 0.5),

  euler.d = FALSE,
  scaled = FALSE,
  ind = FALSE
)

venn_grob <- gTree(
  children = venn_glist
)

p3c <- ggplot() +

  annotation_custom(
    grob = venn_grob,
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1
  ) +

  coord_fixed(
    ratio = 1,
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE,
    clip = "off"
  ) +

  labs(
    title = "MIMAG high-quality criteria"
  ) +

  theme_void(base_size = 8) +

  theme(
    text = element_text(family = "dejavu"),
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0.5,
      margin = margin(b = 8)
    ),
    plot.margin = margin(
      t = 10,
      r = 20,
      b = 10,
      l = 20
    )
  )


# =========================================================================
# Assemble Figure 3
# =========================================================================

fig3_bottom <- ggarrange(
  p3b,
  p3c,
  ncol = 2,
  widths = c(1.2, 1),
  labels = c("b", "c"),
  font.label = list(
    size = 13,
    face = "bold"
  )
)

fig3 <- ggarrange(
  p3a,
  fig3_bottom,
  ncol = 1,
  heights = c(1.3, 1),
  labels = c("a", ""),
  font.label = list(
    size = 13,
    face = "bold"
  )
)

fig3

# =========================================================================
# Export Figure 3
# =========================================================================

ggsave(
  filename = "figure_and_table/Figure3.pdf",
  plot = fig3,
  width = 11,
  height = 7,
  units = "in",
  device = cairo_pdf
)
