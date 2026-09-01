library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# -------------------------------------------------------------------------
# Input
# -------------------------------------------------------------------------

seqsum <- read.delim(
  "1-01_sequencing_assembly_summary.tsv",
  stringsAsFactors = FALSE
)

sample_order <- c(
  "Con", "CW1", "CW2", "GY1", "MP1", "MP2",
  "NH2", "NH3", "PH1", "PH2", "US1"
)

seqsum <- seqsum %>%
  mutate(
    sample_id = factor(sample_id, levels = sample_order),
    assembly_contig_N50_kb = assembly_contig_N50_bp / 1000
  )

stopifnot(
  nrow(seqsum) == 11,
  sum(seqsum$primary_MAG_count) == 6505
)

fig3_fill <- "grey55"

theme_fig3 <- theme_classic(base_size = 10) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(
      size = 11,
      face = "bold",
      hjust = 0.5
    )
  )

# =========================================================================
# Figure 3A. HiFi sequencing yield
# =========================================================================

p3a <- ggplot(
  seqsum,
  aes(x = sample_id, y = hifi_yield_Gb)
) +
  geom_col(
    width = 0.72,
    fill = fig3_fill
  ) +
  geom_text(
    aes(label = sprintf("%.1f", hifi_yield_Gb)),
    vjust = -0.35,
    size = 3
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    x = NULL,
    y = "HiFi sequencing yield (Gb)",
    title = "HiFi sequencing yield"
  ) +
  theme_fig3

# =========================================================================
# Figure 3B. Assembly size and contiguity
# =========================================================================
seqsum <- seqsum %>%
  mutate(
    label_b = paste0(
      "bold('", as.character(sample_id), "')~",
      "'(",
      format(
        assembly_contig_count,
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
      ),
      " contigs)'"
    )
  )

p3b <- ggplot(
  seqsum,
  aes(
    x = assembly_size_Gb,
    y = assembly_contig_N50_kb
  )
) +
  geom_point(
    size = 2.8,
    shape = 21,
    fill = fig3_fill,
    colour = "black"
  ) +
  geom_text_repel(
    aes(label = label_b),
    parse = TRUE,
    size = 3,
    lineheight = 0.9,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  labs(
    x = "Total assembly size (Gb)",
    y = "Assembly contig N50 (kb)",
    title = "Metagenome assembly"
  ) +
  theme_fig3 +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

# =========================================================================
# Figure 3C. Primary MAG count
# =========================================================================

p3c <- ggplot(
  seqsum,
  aes(x = sample_id, y = primary_MAG_count)
) +
  geom_col(
    width = 0.72,
    fill = fig3_fill
  ) +
  geom_text(
    aes(label = primary_MAG_count),
    vjust = -0.35,
    size = 3
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    x = NULL,
    y = "Number of primary MAGs",
    title = "Primary MAG recovery"
  ) +
  theme_fig3

# =========================================================================
# Figure 3D. MAG-associated HiFi bases
# =========================================================================

p3d <- ggplot(
  seqsum,
  aes(
    x = sample_id,
    y = MAG_associated_HiFi_bases_percent
  )
) +
  geom_col(
    width = 0.72,
    fill = fig3_fill
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f",
        MAG_associated_HiFi_bases_percent
      )
    ),
    vjust = -0.35,
    size = 3
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    x = NULL,
    y = "MAG-associated HiFi bases (%)",
    title = "MAG-associated HiFi signal"
  ) +
  theme_fig3

# =========================================================================
# Assemble Figure 3
# =========================================================================

fig3 <- ggarrange(
  p3a, p3b,
  p3c, p3d,
  ncol = 2,
  nrow = 2,
  labels = c("a", "b", "c", "d"),
  font.label = list(
    size = 18,
    face = "bold"
  ),
  align = "hv"
)

fig3

# =========================================================================
# Export Figure 3
# =========================================================================

ggsave(
  filename = "figure_and_table/Figure3.pdf",
  plot = fig3,
  width = 8,
  height = 6,
  units = "in",
  device = cairo_pdf
)
