## Pour faire des graphes de la distribution a posteriori

library(ggplot2)
library(pscl)

Posterior <- as.data.frame(read.table("Posterior_anis_geom.txt"))

gg <- ggplot() +
    stat_density2d(data=Posterior, aes(x=V1, y=V2, fill=..level..), geom="polygon") +
    scale_fill_continuous(low="white",high="black") +
    theme_bw() +
    coord_fixed(ratio=1)
