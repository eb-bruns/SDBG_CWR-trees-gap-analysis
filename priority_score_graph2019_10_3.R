# H. Achicanoy
# CIAT, 2019

library(tidyverse)
library(grid)

fps_cajanus <- readxl::read_excel("D:/OneDrive - CGIAR/Documents/Publications/CWR_Capsicum/code/Capsicum_summary_2019_10_4.xlsx", sheet = 1)
#fps_cajanus <- readxl::read_excel("D:/ToBackup/others/characterization_spooner/_data/Lactuca_species_summary_2019-5_8_cleaned.xlsx", sheet = 1)
#fps_cajanus <- readxl::read_excel("D:/ToBackup/others/characterization_spooner/_data/Cucurbita_species_summary_2019-06-05.xlsx", sheet = 1)
#fps_cajanus

p <- ggplot2::ggplot(data=fps_cajanus, aes(x=`FCSc mean`, y=reorder(Taxon,-`FCSc mean`), fill = "final"), colour="black") + ggplot2::theme_bw()

p <- p + ggplot2::xlim(0,10) + ggplot2::geom_blank()

p <- p + ggplot2::annotate("rect", xmin=0, xmax=25, ymin=0, ymax=nrow(fps_cajanus)+1, alpha=.25, fill="red")
p <- p + ggplot2::annotate("rect", xmin=25, xmax=50, ymin=0, ymax=nrow(fps_cajanus)+1, alpha=.5, fill="orange")
p <- p + ggplot2::annotate("rect", xmin=50, xmax=75, ymin=0, ymax=nrow(fps_cajanus)+1, alpha=.4, fill="yellow")
p <- p + ggplot2::annotate("rect", xmin=75, xmax=100, ymin=0, ymax=nrow(fps_cajanus)+1, alpha=.3, fill="forestgreen")

p <- p + ggplot2::xlab(label="Conservation score") + ggplot2::ylab("")
p <- p + ggplot2::scale_x_continuous(breaks = seq(0, 100, 10))

# Ex situ scores
p <- p + ggplot2::geom_point(data=fps_cajanus, aes(x=SRSex, y=reorder(Taxon,`FCSc mean`), colour = "srs1"), size = 3.5)
p <- p + ggplot2::geom_point(data=fps_cajanus, aes(x=GRSex, y=reorder(Taxon,`FCSc mean`), colour = "grs1"), size = 3.5)
p <- p + ggplot2::geom_point(data=fps_cajanus, aes(x=ERSex, y=reorder(Taxon,`FCSc mean`), colour = "ers1"), size = 3.5)
p <- p + ggplot2::geom_point(data=fps_cajanus, aes(x=FCSex, y=reorder(Taxon,`FCSc mean`), colour = "fcs1"), size = 4)

# In situ scores
p <- p + ggplot2::geom_point(data=fps_cajanus, aes(x=SRSin, y=reorder(Taxon,`FCSc mean`), shape = "srs2"), colour = "dodgerblue", size = 3.5)
p <- p + ggplot2::geom_point(data=fps_cajanus, aes(x=GRSin, y=reorder(Taxon,`FCSc mean`), shape = "grs2"), colour = "slateblue", size = 3.5) # shape = 17
p <- p + ggplot2::geom_point(data=fps_cajanus, aes(x=ERSin, y=reorder(Taxon,`FCSc mean`), shape = "ers2"), colour = "forestgreen", size = 3.5)
p <- p + ggplot2::geom_point(data=fps_cajanus, aes(x=FCSin, y=reorder(Taxon,`FCSc mean`), shape = 'fcs2'), colour = "black", size = 4)

p <- p + ggplot2::geom_point(stat = "identity", shape = 18, colour = "red", size = 8)

# Organize elements
p <- p + ggplot2::theme(panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        panel.grid.major.y = element_line(colour="grey60", linetype="dashed"))
p <- p + ggplot2::theme(axis.text.x  = element_text(size=15),
                        axis.text.y  = element_text(face="italic",size=15),
                        axis.title.x = element_text(face="bold",size=15),
                        axis.title.y = element_text(face="bold",size=15))
p <- p + ggplot2::guides(size = FALSE)

# Legends
p <- p + ggplot2::scale_shape_manual(name   = "In situ score",
                                     values = c(srs2 = 17, grs2 = 17, ers2 = 17, fcs2 = 17),
                                     labels = c("SRSin","GRSin","ERSin","FCSin"),
                                     guide  = ggplot2::guide_legend(override.aes = list(size = c(rep(3.5,3),4),
                                                                                        colour = c("dodgerblue","slateblue","forestgreen","black"))))
p <- p + ggplot2::scale_colour_manual(name   = "Ex situ score",
                                      values = c(srs1 = "dodgerblue", grs1 = "slateblue", ers1 = "forestgreen", fcs1 = "black"),
                                      labels = c("SRSex","GRSex","ERSex","FCSex"),
                                      guide  = ggplot2::guide_legend(override.aes = list(size = c(rep(3.5,3), 4),
                                                                                         colour = c("dodgerblue","slateblue","forestgreen","black"))))
p <- p + ggplot2::scale_fill_manual(name   = "Final score",
                                    values = "red",
                                    labels = "FCSc mean",
                                    guide  = ggplot2::guide_legend(override.aes = list(size = 8,
                                                                                       shape = 18,
                                                                                       colour = "red")))
p <- p + ggplot2::annotate("text",
                           x      = c(12.5, 37.5, 62.5, 87.5),
                           y      = nrow(fps_cajanus) + 0.7, # 13.7
                           label  = c("HP","MP","LP","SC"),
                           colour = "black",
                           size   = 4)
p
ggsave("D:/Capsicum_Fig2a.png", plot = p, device = "png", units = "in", width = 10, height = 12, dpi = 320)
