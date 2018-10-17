library("ggplot2")
library("grid")
library("gridExtra")

setwd("~/gitrepo/sgsdb");

GENES = c('PR', 'RT', 'IN')

PROFILE_DATA = read.table("data/report.csv", header=TRUE, sep=",", comment.char="");

NMUT = list(
  PR = as.double(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=PR, Category=All",]$value),
  RT = as.double(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=RT, Category=All",]$value),
  IN = as.double(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=IN, Category=All",]$value)
)
NMUT_SEQ = list(
  PR = as.double(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=PR, Category=All, NumSequences>1",]$value),
  RT = as.double(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=RT, Category=All, NumSequences>1",]$value),
  IN = as.double(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=IN, Category=All, NumSequences>1",]$value)
)
# NMUT_PT = list(
#   PR = as.integer(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=PR, Category=All, NumPatients>1",]$value),
#   RT = as.integer(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=RT, Category=All, NumPatients>1",]$value),
#   IN = as.integer(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=IN, Category=All, NumPatients>1",]$value)
# )

NUUM = list(
  PR = as.double(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=PR, Category=All, IsUnusual, NotAPOBEC",]$value),
  RT = as.double(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=RT, Category=All, IsUnusual, NotAPOBEC",]$value),
  IN = as.double(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=IN, Category=All, IsUnusual, NotAPOBEC",]$value)
)
NUUM_SEQ = list(
  PR = as.double(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=PR, Category=All, NumSequences>1, IsUnusual, NotAPOBEC",]$value),
  RT = as.double(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=RT, Category=All, NumSequences>1, IsUnusual, NotAPOBEC",]$value),
  IN = as.double(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=IN, Category=All, NumSequences>1, IsUnusual, NotAPOBEC",]$value)
)
# NUUM_PT = list(
#   PR = as.integer(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=PR, Category=All, NumPatients>1, IsUnusual, NotAPOBEC",]$value),
#   RT = as.integer(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=RT, Category=All, NumPatients>1, IsUnusual, NotAPOBEC",]$value),
#   IN = as.integer(PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=IN, Category=All, NumPatients>1, IsUnusual, NotAPOBEC",]$value)
# )

PUUM = list(
  PR = as.double(sub("%", "", PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=PR, Category=All, IsUnusual, NotAPOBEC",]$percent))/100,
  RT = as.double(sub("%", "", PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=RT, Category=All, IsUnusual, NotAPOBEC",]$percent))/100,
  IN = as.double(sub("%", "", PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=IN, Category=All, IsUnusual, NotAPOBEC",]$percent))/100
)
PUUM_SEQ = list(
  PR = as.double(sub("%", "", PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=PR, Category=All, NumSequences>1, IsUnusual, NotAPOBEC",]$percent))/100,
  RT = as.double(sub("%", "", PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=RT, Category=All, NumSequences>1, IsUnusual, NotAPOBEC",]$percent))/100,
  IN = as.double(sub("%", "", PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=IN, Category=All, NumSequences>1, IsUnusual, NotAPOBEC",]$percent))/100
)
PUUM_PT = list(
  PR = as.double(sub("%", "", PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=PR, Category=All, NumPatients>1, IsUnusual, NotAPOBEC",]$percent))/100,
  RT = as.double(sub("%", "", PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=RT, Category=All, NumPatients>1, IsUnusual, NotAPOBEC",]$percent))/100,
  IN = as.double(sub("%", "", PROFILE_DATA[PROFILE_DATA$name == "# Mutations per Patient" & PROFILE_DATA$subset == "Gene=IN, Category=All, NumPatients>1, IsUnusual, NotAPOBEC",]$percent))/100
)

APM = list(PR = 0.0554, RT = 0.0307, IN = 0.0621)
APM_SEQ = list(PR = 0.0359, RT = 0.0380, IN = 0.0513)
APM_PT = list(PR = 0.0415, RT = 0.0439, IN = 0.0562)

plot_nmut = function(gene) {
  plot = ggplot(data) +
    geom_histogram(aes(x = X..Mutations / X..Patients), bins=60) +
    geom_vline(aes(xintercept = NMUT[[gene]]), color = "red", linetype = "dashed") +
    geom_vline(aes(xintercept = NMUT_SEQ[[gene]]), color = "blue", linetype = "dashed") +
    # geom_vline(aes(xintercept = NMUT_PT[[gene]]), color = "green", linetype = "dashed") +
    scale_y_continuous(
      name = NULL
      # limits = c(0, 75)
    ) +
    scale_x_continuous(
      name = "# Mutations"
      # labels = if (gene == "IN") waiver() else NULL,
      # limits = c(150, 850)
    );
  plot;
}

plot_nuum = function(gene) {
  plot = ggplot(data) +
    geom_histogram(aes(x = X..Unusual.Mutations / X..Patients), bins=60) +
    geom_vline(aes(xintercept = NUUM[[gene]]), color = "red", linetype = "dashed") +
    geom_vline(aes(xintercept = NUUM_SEQ[[gene]]), color = "blue", linetype = "dashed") +
    # geom_vline(aes(xintercept = NUUM_PT[[gene]]), color = "green", linetype = "dashed") +
    scale_y_continuous(
      name = NULL
      # labels = NULL,
      # limits = c(0, 75)
    ) +
    scale_x_continuous(
      name = "# Unusual Mutations"
      # labels = if (gene == "IN") waiver() else NULL,
      # limits = c(0, 250)
    );
  plot;
}

plot_puum = function(gene) {
  plot = ggplot(data) +
    geom_histogram(aes(x = X..Unusual.Mutations.1), bins=60) +
    geom_vline(aes(xintercept = PUUM[[gene]]), color = "red", linetype = "dashed") +
    geom_vline(aes(xintercept = PUUM_SEQ[[gene]]), color = "blue", linetype = "dashed") +
    # geom_vline(aes(xintercept = PUUM_PT[[gene]]), color = "green", linetype = "dashed") +
    scale_y_continuous(
      name = NULL
      # labels = NULL,
      # limits = c(0, 75)
    ) +
    scale_x_continuous(
      name = "% Unusual Mutations",
      labels = scales::percent
      # labels = if (gene == "IN") scales::percent else NULL,
      # limits = c(0, 0.3)
    );
  plot;
}

for (gene in GENES) {
  filename = sprintf("local/permut.new/permut.%s.1000.o.txt", gene);
  data = read.table(filename, header=TRUE, sep="\t", comment.char="");
  plots = list(plot_nmut(gene), plot_nuum(gene), plot_puum(gene));
  pdf(sprintf('data/PermutationTestGraphOrdinaryMuts2x2%s.pdf', gene), width=8/3*2, height=4);
  print(grid.arrange(grobs = plots, ncol = 2, nrow = 2));
  dev.off();
}

plot_apm = function(gene) {
  plot = ggplot(data) +
    geom_histogram(aes(x = X..APOBEC.Mutations.1), binwidth=0.0002) +
    geom_vline(aes(xintercept = APM[[gene]]), color = "red", linetype = "dashed") +
    geom_vline(aes(xintercept = APM_SEQ[[gene]]), color = "blue", linetype = "dashed") +
    geom_vline(aes(xintercept = APM_PT[[gene]]), color = "green", linetype = "dashed") +
    scale_y_continuous(
      name = NULL,
      labels = NULL,
      limits = c(0, 75)
    ) +
    scale_x_continuous(
      name = if (gene == "IN") "APOBEC Mutations" else NULL,
      labels = if (gene == "IN") scales::percent else NULL,
      limits = c(0, 0.075)
    );
  plot;
}

