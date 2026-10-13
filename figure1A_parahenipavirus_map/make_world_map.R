# ============================================================
# Parahenipavirus World Distribution Map
# Slightly enhanced saturation
# ============================================================
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
# ------------------------------------------------------------
# 1. Load world map
# ------------------------------------------------------------
world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)
# ------------------------------------------------------------
# 2. Define country colors
# ------------------------------------------------------------
country_colors <- c(
  "MDG" = "#D07B42",  # Madagascar - orange
  "USA" = "#416FB5",  # USA - blue
  "AUS" = "#B99A72",  # Australia - warm taupe
  "KOR" = "#7D68B8",  # Republic of Korea - purple
  "CHN" = "#CC454A",  # China - red
"BEL" = "#F2C94C",  # Belgium - bright muted yellow
  "GHA" = "#83C363",  # Ghana - green
  "DEU" = "#627D82",  # Germany - slate blue-gray
  "GIN" = "#5EB8D6",  # Guinea - light blue
  "MYS" = "#CF6AB7",  # Malaysia - pink
  "TWN" = "#CC454A"   # Taiwan - same red as China
)
# ------------------------------------------------------------
# 3. Assign colors
# ------------------------------------------------------------
world$map_color <- country_colors[world$adm0_a3]
# All other countries
world$map_color[is.na(world$map_color)] <- "#E6E6E6"
# ------------------------------------------------------------
# 4. Countries with emphasized outlines
# ------------------------------------------------------------
outline_codes <- c(
  "MDG",
  "BEL",
  "GHA",
  "GIN",
  "DEU",
  "KOR",
  "MYS"
)
outline_countries <- world %>%
  filter(adm0_a3 %in% outline_codes)
# ------------------------------------------------------------
# 5. Draw map
# ------------------------------------------------------------
p <- ggplot() +
  # Base world map
  geom_sf(
    data = world,
    aes(fill = map_color),
    color = "#A6A6A6",
    linewidth = 0.12
  ) +
  # Stronger outlines for selected countries
  geom_sf(
    data = outline_countries,
    fill = NA,
    color = "#777777",
    linewidth = 0.22
  ) +
  scale_fill_identity() +
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    plot.margin = margin(5, 5, 5, 5)
  )
# ------------------------------------------------------------
# 6. Display
# ------------------------------------------------------------
p
# ------------------------------------------------------------
# 7. Export PNG
# ------------------------------------------------------------
ggsave(
  filename = "Parahenipavirus_world_map_saturated.png",
  plot = p,
  width = 12,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)
# ------------------------------------------------------------
# 8. Export PDF
# ------------------------------------------------------------
ggsave(
  filename = "Parahenipavirus_world_map_saturated.pdf",
  plot = p,
  width = 12,
  height = 6,
  units = "in",
  bg = "white"
)
