library(hexSticker)
library(tidyverse)
library(RGCxGC)
library(colorRamps)
library(magick)
library(ggsci) # color palettes
library(showtext)

colfunc<- colorRampPalette(rev(c("red","yellow","springgreen",
                                 "blue", "white")))

GB_01 <- read_chrom(name = 'data/01_GB.cdf', 5L,
                    x_cut = c(26.56, 29), y_cut = c(2, 2.5))
GB_01_bsl <- baseline_corr(GB_01, lambda = 1e3)
GB_01_bsl@chromatogram[GB_01_bsl@chromatogram < 8e3] <- 0
plot(GB_01_bsl, col = colfunc(50), nlevels = 50, main = "TIC", type = "c")

png(filename = "./sticker/chrom.png",width = 4, height = 3, units = "in",
     res = 300, bg = "transparent", pointsize = 2)
par(bg = NA, mar = c(.1, .1, 1, .1))
plot(GB_01_bsl, color.palette = colfunc, nlevels = 50, main = "TIC")
dev.off()

font_add_google("Oswald", "os")
font_add_google("Lora", "os")
showtext_auto()


chrom <- image_read("sticker/chrom3.png")
sticker(chrom, package = "RGCxGC", s_width = 2.5,s_height = 1.25, 
        s_x = 1, s_y = 0.9, h_color = "#00A1AE", h_fill = "#EAEAFF",
        p_color = "black", p_size = 15, spotlight = T, p_y = 1.6,
        url = "git DanielQuiroz97/RGCxGC", u_size = 3.5, p_family = "os")

