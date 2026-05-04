library(dplyr)
library(ggplot2)
library(ggrepel)
library(maps)
library(readr)
library(stringr)

hospital_geo_path <- "/Users/saborpete/Desktop/Peter/Postdoc/CLIF-HROHCA/reference/clif_hospital_geography.csv"
output_dir <- "output/clif_hospital_map"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

hospitals <- read_csv(hospital_geo_path, show_col_types = FALSE) %>%
  mutate(
    hospital_county_fips = as.integer(hospital_county_fips),
    site_name = str_trim(site_name),
    hospital_id_name = str_trim(hospital_id_name)
  )

states <- map_data("state")
counties <- map_data("county") %>%
  mutate(polyname = paste(region, subregion, sep = ","))

county_fips <- maps::county.fips %>%
  as_tibble() %>%
  transmute(polyname = polyname, hospital_county_fips = as.integer(fips))

county_centroids <- counties %>%
  inner_join(county_fips, by = "polyname") %>%
  group_by(hospital_county_fips) %>%
  summarise(
    lon = mean(range(long), na.rm = TRUE),
    lat = mean(range(lat), na.rm = TRUE),
    .groups = "drop"
  )

hospital_counties <- hospitals %>%
  count(site_name, hospital_county_fips, hospital_county_name, name = "hospital_count") %>%
  left_join(county_centroids, by = "hospital_county_fips")

system_labels <- hospitals %>%
  left_join(county_centroids, by = "hospital_county_fips") %>%
  group_by(site_name) %>%
  summarise(
    lon = weighted.mean(lon, w = rep(1, n()), na.rm = TRUE),
    lat = weighted.mean(lat, w = rep(1, n()), na.rm = TRUE),
    hospital_count = n(),
    county_count = n_distinct(hospital_county_fips),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(site_name, " (", hospital_count, ")")
  )

system_palette <- c(
  "Emory" = "#2A9D8F",
  "JHU" = "#4E79A7",
  "Michigan" = "#EDC948",
  "NU" = "#B07AA1",
  "OHSU" = "#59A14F",
  "Penn" = "#E15759",
  "RUMC" = "#9C755F",
  "UCMC" = "#F28E2B",
  "UCSF" = "#76B7B2",
  "UMN" = "#5B5F97"
)

map_plot <- ggplot() +
  geom_polygon(
    data = states,
    aes(x = long, y = lat, group = group),
    fill = "#F7F4EE",
    color = "#D8D1C6",
    linewidth = 0.25
  ) +
  geom_point(
    data = hospital_counties,
    aes(x = lon, y = lat, size = hospital_count, fill = site_name),
    shape = 21,
    color = "white",
    stroke = 0.7,
    alpha = 0.95
  ) +
  geom_label_repel(
    data = system_labels,
    aes(x = lon, y = lat, label = label, fill = site_name),
    color = "white",
    family = "Helvetica",
    fontface = "bold",
    size = 3.8,
    label.size = 0,
    label.r = unit(0.1, "lines"),
    box.padding = 0.45,
    point.padding = 0.45,
    min.segment.length = 0,
    segment.color = "#6D6A65",
    segment.size = 0.35,
    seed = 20260427,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = system_palette) +
  scale_size_area(
    max_size = 12,
    breaks = c(1, 2, 3, 6, 9),
    name = "Hospitals in county"
  ) +
  coord_quickmap(
    xlim = c(-125, -66),
    ylim = c(24, 50)
  ) +
  labs(
    title = "CLIF Hospitals and Health Systems",
    subtitle = paste0(
      nrow(hospitals), " hospitals across ",
      n_distinct(hospitals$site_name), " health systems and ",
      n_distinct(hospitals$hospital_county_fips), " counties"
    ),
    caption = "Point size represents the number of CLIF hospitals in the same county. Labels show health system and hospital count."
  ) +
  theme_void(base_family = "Helvetica") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(
      color = "#252525",
      face = "bold",
      size = 24,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      color = "#4A4A4A",
      size = 12,
      margin = margin(b = 14)
    ),
    plot.caption = element_text(
      color = "#666666",
      size = 9,
      hjust = 0,
      margin = margin(t = 12)
    ),
    plot.margin = margin(22, 24, 18, 24),
    legend.position = "none"
  ) +
  guides(fill = "none", size = "none")

png_path <- file.path(output_dir, "clif_hospitals_health_systems_map.png")
pdf_path <- file.path(output_dir, "clif_hospitals_health_systems_map.pdf")

ggsave(png_path, map_plot, width = 13.333, height = 7.5, dpi = 300, bg = "white")
ggsave(pdf_path, map_plot, width = 13.333, height = 7.5, bg = "white")

message("Wrote: ", normalizePath(png_path))
message("Wrote: ", normalizePath(pdf_path))
