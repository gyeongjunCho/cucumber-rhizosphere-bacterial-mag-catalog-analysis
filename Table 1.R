library(dplyr)
library(flextable)
library(officer)

# -------------------------------------------------------------------------
# Input
# -------------------------------------------------------------------------

meta <- read.delim(
  "2-01_sample_metadata.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_order <- c(
  "Con", "CW1", "CW2", "GY1", "MP1", "MP2",
  "NH2", "NH3", "PH1", "PH2", "US1"
)

# -------------------------------------------------------------------------
# Prepare Table 1
# -------------------------------------------------------------------------

tab1 <- meta %>%
  mutate(
    sample_id = factor(sample_id, levels = sample_order),
    `Source-soil site` = ifelse(
      is.na(source_soil_site), "—", source_soil_site
    ),
    `Collection date` = ifelse(
      is.na(source_soil_collection_date),
      "—",
      source_soil_collection_date
    ),
    `Source vegetation` = ifelse(
      is.na(source_habitat_vegetation),
      "—",
      source_habitat_vegetation
    ),
    Latitude = ifelse(
      is.na(source_soil_latitude),
      "—",
      sprintf("%.3f", source_soil_latitude)
    ),
    Longitude = ifelse(
      is.na(source_soil_longitude),
      "—",
      sprintf("%.3f", source_soil_longitude)
    ),
    `Field soil (%)` = sprintf("%.0f", field_soil_fraction * 100)
  ) %>%
  arrange(sample_id) %>%
  transmute(
    `Sample ID` = as.character(sample_id),
    `Source-soil site`,
    `Collection date`,
    `Source vegetation`,
    Latitude,
    Longitude,
    `Field soil (%)`
  )

# -------------------------------------------------------------------------
# Flextable
# -------------------------------------------------------------------------

ft <- flextable(tab1) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(
    j = c("Source-soil site", "Source vegetation"),
    align = "left",
    part = "body"
  ) %>%
  compose(
    j = "Source vegetation",
    value = as_paragraph(as_i(`Source vegetation`)),
    part = "body"
  ) %>%
  valign(valign = "center", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  padding(padding = 3, part = "all") %>%
  autofit()

# -------------------------------------------------------------------------
# Footnote
# -------------------------------------------------------------------------

note <- paste0(
  "Source-soil coordinates refer to the original coastal plant-habitat ",
  "sampling locations. All cucumber experiments were conducted at ",
  "35.83° N, 127.04° E using Cucumis sativus cv. Baekdadagi. ",
  "Rhizosphere samples were collected on 3 June 2024, 20 days after sowing. ",
  "Thirty plants were pooled per treatment to generate one composite ",
  "rhizosphere sample. Each sample was sequenced on the PacBio Revio ",
  "platform using two SMRT Cells."
)

# -------------------------------------------------------------------------
# Export to Word
# -------------------------------------------------------------------------

doc <- read_docx() %>%
  body_add_par(
    "Table 1. Sample and source-soil metadata.",
    style = "heading 1"
  ) %>%
  body_add_flextable(ft) %>%
  body_add_par("") %>%
  body_add_par(note)

print(
  doc,
  target = "figure_and_table/Table1.docx"
)
