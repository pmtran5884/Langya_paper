# Figure 1A — parahenipavirus world distribution map

`make_world_map.R` is the figure code exactly as used for the paper. It colours
the countries with a reported parahenipavirus by a fixed ISO-A3 → hex lookup,
greys everything else, and writes the panel at 12 × 6 in, 600 dpi.

## Run

```bash
cd output && Rscript ../make_world_map.R
```

Writes `Parahenipavirus_world_map_saturated.png` and `.pdf` into the working
directory, which is why the command runs from `output/`. (Printing the ggplot
object also drops an incidental `Rplots.pdf`; it is gitignored and not part of
the figure.)

## Requirements

R ≥ 4.0 with `ggplot2`, `sf`, `dplyr`, `rnaturalearth`, `rnaturalearthdata`:

```r
install.packages(c("ggplot2", "sf", "dplyr", "rnaturalearth", "rnaturalearthdata"))
```

Country geometry comes from Natural Earth via `rnaturalearth`
(`ne_countries(scale = "medium")`); no other input data is needed, and the
version shipped in `rnaturalearthdata` is what fixes the borders, so the panel
is reproducible from the package versions alone.

## Committed output

`output/Parahenipavirus_world_map_saturated.{png,pdf}` — the rendered panel, so
a reviewer can compare without installing the R spatial stack.

Country–colour assignments (one entry per country in the figure) live in the
`country_colors` vector at the top of the script; `outline_codes` lists the
small countries that get a darker border so they stay visible at print size.
