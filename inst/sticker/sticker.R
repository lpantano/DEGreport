p <- ggplot(data.frame(x=1,y=1,image="inst/sticker/flag_hex.png"), aes(x,y)) +
    geom_image(aes(image=image), size=1.06) + theme_void()

hexSticker::sticker(p,
                    package="DEGreport", 
                    s_x = 1.03,
                    s_y = 0.98, 
                    s_width = 2.05,
                    s_heigh = 2.05,
                    p_x = 1,
                    p_y = 1.17,
                    h_color = "white",
                    h_fill = "transparent",
                    h_size = 1,
                    p_color = "black",
                    p_size = 24,
                    filename="inst/sticker/degreport.png")
