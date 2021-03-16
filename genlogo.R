#install.packages("hexSticker")
library(hexSticker)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("emojifont")
library(emojifont)
#install.packages("fontawesome")
library(fontawesome)
aaa=fa(name = "fas fa-asterisk", fill = "cyan")

print_fontawesome<- function (x, view = interactive(), ...) 
  {
    dots <- list(...)
    html <- paste(x, collapse = "\n")
    c("<!DOCTYPE html>", "<html>", "<head>", "<meta charset=\"utf-8\">", 
      "</head>", "<body>", html, "</body>", "</html>") %>% 
      paste(collapse = "\n") %>% htmltools::HTML() %>% htmltools::html_print()
    return(html)
  }

library(rsvg)

cat(print_fontawesome(aaa, view=FALSE),file="man/figures/fa-asterisk.svg")
SVGres <- rsvg_svg("man/figures/fa-asterisk.svg","man/figures/Cairo-fa-asterisk.svg")
SVGres <- rsvg_svg("man/figures/fa-asterisk.svg")
PDFres <- rsvg_pdf("man/figures/fa-asterisk.svg","man/figures/Cairo-fa-asterisk.pdf")
PDFres <- rsvg_pdf("man/figures/fa-asterisk.svg")

library(png)
library(grid)
library(cowplot)
library(magick)
theme_set(theme_cowplot())
r = ggdraw() +
  draw_image("man/figures/Cairo-fa-asterisk.pdf", scale = .55)

sticker(
  r,
  package = "c060",
  p_size = 8,
  s_x = 0.95,
  s_y = 0.68,
  s_width = 1.7,
  s_height = 1.3,
  p_x = 1,
  p_y = 1.3,
  url = "https://cran.r-project.org/package=c060",
  u_color = "white",
  u_size = 1.1,
  h_fill = "black",
  h_color = "grey",
  filename = "man/figures/logo.png"
)
