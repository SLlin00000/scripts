variants <- read.delim(
  "haplotype/LCORL_missense_plot_positions.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

domains <- read.delim(
  "protein/interproscan/LCORL_interpro_domain_summary.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

protein_length <- 1904
domain_start <- domains$Start[1]
domain_end <- domains$End[1]

variants$Within_domain <-
  variants$Protein_position >= domain_start &
  variants$Protein_position <= domain_end

variant_colours <- ifelse(
  variants$Within_domain,
  "#D55E00",
  "#3F566B"
)

stem_heights <- c(1.05, 1.42, 1.05, 1.42, 1.05, 1.42)

draw_protein_map <- function() {
  par(
    mar = c(4.5, 2, 5.8, 1.5),
    xpd = NA,
    family = "sans"
  )

  plot(
    NA,
    xlim = c(0, protein_length + 30),
    ylim = c(-0.25, 1.85),
    axes = FALSE,
    xlab = "",
    ylab = "",
    bty = "n"
  )

  rect(
    1, 0.18,
    protein_length, 0.42,
    col = "#D9DEE5",
    border = "#566573",
    lwd = 1.4
  )

  rect(
    domain_start, 0.13,
    domain_end, 0.47,
    col = "#56B4E9",
    border = "#0072B2",
    lwd = 1.5
  )

  text(
    mean(c(domain_start, domain_end)),
    0.30,
    paste0("DUF4553\n", domain_start, "–", domain_end),
    cex = 0.82,
    font = 2
  )

  for (i in seq_len(nrow(variants))) {
    position <- variants$Protein_position[i]

    segments(
      position, 0.47,
      position, stem_heights[i],
      col = variant_colours[i],
      lwd = 1.6
    )

    points(
      position,
      stem_heights[i],
      pch = 21,
      bg = variant_colours[i],
      col = "white",
      cex = 1.55,
      lwd = 1.2
    )

    text(
      position,
      stem_heights[i] + 0.14,
      variants$Variant[i],
      col = variant_colours[i],
      font = 2,
      cex = 0.9
    )
  }

  axis(
    1,
    at = c(1, 500, 1000, 1500, 1904),
    labels = c("1", "500", "1000", "1500", "1904"),
    las = 1,
    col = "#566573",
    col.axis = "#263238"
  )

  mtext(
    "Amino-acid position",
    side = 1,
    line = 2.6,
    cex = 0.95
  )

  title(
    main = "LCORL missense variants and annotated protein domain",
    cex.main = 1.25,
    line = 3.2
  )

  mtext(
   paste(
      "Bos taurus | ENSBTAT00000066859.1 | ENSBTAP00000060198 |",
      protein_length,
      "aa"
    ),
    side = 3,
    line = 1.4,
    cex = 0.9,
    col = "#4D4D4D"
  )

  legend(
    "topleft",
    inset = c(0.01, 0.01),
    legend = c(
      "Outside annotated domain",
      "Within DUF4553",
      "DUF4553 (IPR028104; PF15090)"
    ),
    pch = c(21, 21, 22),
    pt.bg = c("#3F566B", "#D55E00", "#56B4E9"),
    col = c("#3F566B", "#D55E00", "#0072B2"),
    pt.cex = 1.2,
    bty = "n",
    cex = 0.82
  )
}

png(
  "figures/LCORL_protein_variant_map.png",
  width = 2700,
  height = 1350,
  res = 300,
  type = "cairo-png",
  bg = "white"
)
draw_protein_map()
dev.off()

pdf(
  "figures/LCORL_protein_variant_map.pdf",
  width = 10,
  height = 5,
  useDingbats = FALSE
)
draw_protein_map()
dev.off()

message("Created PNG and PDF in figures/")

