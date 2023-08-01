library(ggplot2)
world_map <- map_data("world")
world_map$subregion <- ifelse(world_map$region=="Madagascar","Madagascar",ifelse(
  world_map$region=="Tanzania","Tanzania",ifelse(
    world_map$region=="China","China",ifelse(
      world_map$region=="Taiwan","China",ifelse(
        world_map$region=="India","India",ifelse(
          world_map$region=="Denmark","Denmark",ifelse(
            world_map$region=="Estonia","Estonia",ifelse(
              world_map$region=="Finland","Finland",ifelse(
                world_map$region=="France","France",ifelse(
                  world_map$region=="Germany","Germany",ifelse(
                    world_map$region=="Iceland","Iceland",ifelse(
                      world_map$region=="Ireland","Ireland",ifelse(
                        world_map$region=="Italy","Italy",ifelse(
                          world_map$region=="Luxembourg","Luxembourg",ifelse(
                            world_map$region=="Netherlands","Netherlands",ifelse(
                              world_map$region=="Russia","Russia",ifelse(
                                world_map$region=="Sweden","Sweden",ifelse(
                                  world_map$region=="UK","UK",ifelse(
                                    world_map$region=="Canada","Canada",ifelse(
                                      world_map$region=="USA","USA","Others"))))))))))))))))))))
world_map$subregion <- factor(world_map$subregion,levels = c("Madagascar","Tanzania","China","India","Denmark","Estonia","Finland","France","Germany","Iceland","Ireland","Italy","Luxembourg","Netherlands","Russia","Sweden","UK","Canada","USA","Others"))
ggplot(world_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill= subregion), colour = "grey50") +
  scale_x_continuous(breaks = seq(-180, 210, 45), labels = function(x){paste0(x, "°")}) +
  scale_y_continuous(breaks = seq(-60, 100, 30), labels = function(x){paste0(x, "°")}) +
  scale_fill_manual(values = c("#4abbbe","#ee989e","#80b1d3","#ec9675","#ec5d69","#ddb743","#b9d1c5","#979a6f","#c9b27b","#db5b38","#8eb853","#eadfc3","#eec16c","#379160","#8dd3c7","#8684b9","#c3503f","#bebada","#fdb462","#f2f2f2"))+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank())
ggsave("map_new.pdf",width = 8.87,height = 4.2)