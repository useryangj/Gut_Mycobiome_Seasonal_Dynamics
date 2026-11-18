library(reshape2)
library(FSA)
library(ggplot2)
library(ggpubr)
library(cowplot)
macro <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/macro_fungi.txt", stringsAsFactors = FALSE)
genus <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/chouping/genus_table1.txt", stringsAsFactors = FALSE)
out   <- genus[ genus$Genus %in% macro$Genus , ]
write.table(out,"C:/Users/xuzy0/Desktop/xzy/split/beiyong/macro_table.txt",row.names = T,col.names = T,sep = "\t",quote = F)
data <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/beiyong/macro_table.txt",row.names = 1)
data <- sweep(data,2,colSums(data),"/")
data <- data[order(rowSums(data),decreasing = T),]
data <- data[c(1:10),]
x <- rownames(data)
data <- cbind(rownames(data),data)
data <- melt(data)
data$period <- ""
data[grep("^D",data$variable),4] <- "Spring"
data[grep("^A",data$variable),4] <- "Summer"
data[grep("^B",data$variable),4] <- "Autumn"
data[grep("^C",data$variable),4] <- "Winter"
table(data$period)
colnames(data) <- c("Species","sample","abundance","period")
data$Species <- factor(data$Species,levels = x)
data$period <- factor(data$period,levels = c("Spring", "Summer","Autumn","Winter"))
ggplot(data,aes(Species,abundance,fill=period)) +
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=1,color="black",size=0.5)+
  scale_fill_manual(values = c("#80cba4","#ebd49c","#ba3e45","#4965b0"))+
  scale_y_continuous(limits = c(0, NA),
                     expand = c(0, 0)) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
library(networkD3)
library(dplyr)
library(jsonlite)
library(htmlwidgets)
data$period <- factor(data$period, levels = c("Spring", "Summer", "Autumn", "Winter"))
links <- data %>%
  group_by(Species, period) %>%
  summarise(total_abundance = sum(abundance), .groups = 'drop') %>%
  filter(total_abundance > 0)
colnames(links) <- c("source", "target", "value")
species_nodes <- sort(unique(as.character(links$source)))
period_order <- c("Spring", "Summer", "Autumn", "Winter")
actual_periods_ordered <- period_order[period_order %in% links$target]
nodes <- data.frame(name = c(species_nodes, actual_periods_ordered))
links$source_id <- match(links$source, nodes$name) - 1
links$target_id <- match(links$target, nodes$name) - 1
species_colors <- c(
  "Russula" = "#B6b3d6", "Cortinarius" = "#ffd47d", "Pholiota" = "#f8b2a2",
  "Phaeotremella" = "#7bc0cd", "Auricularia" = "#F7B2C7", "Hypholoma" = "#51999f",
  "Stereum" = "#F0F8FF", "Trametes" = "#afcdb1", "Leucophleps" = "#bcbd22",
  "Lactarius" = "#c49892"
)
period_colors <- c(
  "Spring" = "#80cba4", "Summer" = "#fbda83", "Autumn" = "#e69a9b", "Winter" = "#c0dbe5"
)
color_map <- c(species_colors, period_colors)
domain_js <- toJSON(names(color_map), auto_unbox = TRUE)
range_js <- toJSON(unname(color_map), auto_unbox = TRUE)
color_scale_js <- paste0('d3.scaleOrdinal().domain(', domain_js, ').range(', range_js, ')')
sankey_plot <- sankeyNetwork(
  Links = links, Nodes = nodes, Source = "source_id", Target = "target_id",
  Value = "value", NodeID = "name", NodeGroup = "name",
  LinkGroup = "target", colourScale = color_scale_js,
  units = "Abundance", fontSize = 15, nodeWidth = 30,
  iterations = 0
)
sankey_plot_custom <- onRender(
  sankey_plot,
  '
    function(el, x) {
      var nodes = d3.select(el).selectAll(".node");
      var links = d3.select(el).selectAll(".link");
      var seasons = ["Spring", "Summer", "Autumn", "Winter"];
      nodes.select("text")
        .each(function(d) {
          if (seasons.indexOf(d.name) === -1) {
            d3.select(this).style("font-style", "italic");
          }
        });
      var sourceNodes = nodes.filter(function(d) { return seasons.indexOf(d.name) === -1; });
      var leftY0 = d3.min(sourceNodes.data(), function(d) { return d.y; });
      var leftY1 = d3.max(sourceNodes.data(), function(d) { return d.y + d.dy; });
      var totalLeftHeight = leftY1 - leftY0;
      var targetNodes = nodes.filter(function(d) { return seasons.indexOf(d.name) > -1; })
        .sort(function(a, b) { return seasons.indexOf(a.name) - seasons.indexOf(b.name); });
      var n = targetNodes.size();
      var totalTargetDy = 0;
      targetNodes.each(function(d) { totalTargetDy += d.dy; });
      var newGap = (totalLeftHeight - totalTargetDy) / (n - 1);
      var currentY = leftY0;
      targetNodes.each(function(d) {
        d.y = currentY;
        d3.select(this).attr("transform", "translate(" + d.x + "," + d.y + ")");
        currentY += d.dy + newGap;
      });
      links.attr("d", d3.sankey().link());
    }
  '
)
sankey_plot_custom
library(webshot2)
library(htmlwidgets)
saveWidget(sankey_plot_custom, "temp_sankey.html", selfcontained = TRUE)
webshot2::install_chrome()
webshot(
  url = "temp_sankey.html",
  file = "sankey_diagram.pdf",
  vwidth = 1000,
  vheight = 700
)
macro_fungi <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/macro_fungi.txt", stringsAsFactors = FALSE)
macro_fungi$Species <- paste(macro_fungi$Genus, macro_fungi$species, sep = "_")
species_table <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/chouping/species_table1.txt", stringsAsFactors = FALSE)
macro_species <- species_table[species_table$Species %in% macro_fungi$Species, ]
print("Original species_table dimensions:")
print(dim(species_table))
print("Filtered macro_species dimensions:")
print(dim(macro_species))
print("Matching species count:")
print(length(unique(macro_fungi$Species)))
print("First few matching species:")
print(head(macro_fungi$Species))
print("First few species_table species:")
print(head(species_table$Species))
write.table(macro_species, "C:/Users/xuzy0/Desktop/xzy/split/beiyong/macro_species.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
macro_species_data <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/beiyong/macro_species.txt", stringsAsFactors = FALSE)
macro_species_group <- data.frame(row.names = macro_species_data$Species,
                                  Jul = rowSums(macro_species_data[,2:11], na.rm = TRUE),
                                  Aug = rowSums(macro_species_data[,12:21], na.rm = TRUE),
                                  Sep = rowSums(macro_species_data[,22:26], na.rm = TRUE),
                                  Oct = rowSums(macro_species_data[,27:40], na.rm = TRUE),
                                  Nov = rowSums(macro_species_data[,41:51], na.rm = TRUE),
                                  Dec = rowSums(macro_species_data[,52:56], na.rm = TRUE),
                                  Jan = rowSums(macro_species_data[,57:66], na.rm = TRUE),
                                  Feb = rowSums(macro_species_data[,67:71], na.rm = TRUE),
                                  Mar = rowSums(macro_species_data[,72:76], na.rm = TRUE),
                                  Apr = rowSums(macro_species_data[,77:86], na.rm = TRUE),
                                  May = rowSums(macro_species_data[,87:96], na.rm = TRUE),
                                  Jun = rowSums(macro_species_data[,97:101], na.rm = TRUE))
x <- rownames(macro_species_group)
macro_species_group <- sweep(macro_species_group, 2, colSums(macro_species_group), "/")
rownames(macro_species_group) <- x
macro_species_group <- macro_species_group[order(rowSums(macro_species_group), decreasing = T), ]
macro_species_group <- macro_species_group[c(1:min(20, nrow(macro_species_group))), ]
if (nrow(macro_species_group) >= 20) {
  macro_species_group["Other", ] <- (1 - colSums(macro_species_group))
}
x <- rownames(macro_species_group)
macro_species_group <- cbind(rownames(macro_species_group), macro_species_group)
macro_species_group <- melt(macro_species_group)
colnames(macro_species_group) <- c("species", "months", "abundance")
macro_species_group$species <- factor(macro_species_group$species, levels = rev(x))
macro_species_group$months <- factor(macro_species_group$months, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                                                            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
p_macro <- ggplot(macro_species_group, aes(months, abundance, fill = species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(  "#d1e6c1", "#Face37",  "#E15759", "#76B7B2", "#666666", "#F28E2B",
                                 "#B07AA1", "#FF9DA7", "#8c564b",  "#BAB0AC", "#66A61E","#B6992D",
                                 "#D37295",  "#B2DF8A", "#9467bd", "#5ab1ef", "#F2F3C7","#4E79A7",
                                 "#C4A38a", "#86BCB6", "#f5e8f1"  )) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))+
guides(fill = guide_legend(ncol = 2))
p_macro
library(ggplot2)
library(reshape2)
library(ggalluvial)
species_table <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/beiyong/macro_species.txt", stringsAsFactors = FALSE)
print("First few rows of species_table:")
print(head(species_table))
print("Column names:")
print(colnames(species_table))
print("First column values:")
print(head(species_table$Species))
species_mat <- as.matrix(species_table[, -1])
rownames(species_mat) <- species_table$Species
month_groups <- sapply(strsplit(colnames(species_mat), "\\."), `[`, 2)
species_group <- t(rowsum(t(species_mat), group = month_groups))
month_order <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
species_group <- as.data.frame(species_group[, intersect(month_order, colnames(species_group))])
species_group <- sweep(species_group,2,colSums(species_group),"/")
species_group <- species_group[order(rowSums(species_group),decreasing = T),]
species_group <- species_group[c(1:20),]
species_group["Other",] <- (1-colSums(species_group))
x <- rownames(species_group)
print("Species names (x):")
print(x)
print("Column names of species_group:")
print(colnames(species_group))
species_group_df <- as.data.frame(species_group)
species_group_df$species <- rownames(species_group_df)
species_melted <- melt(species_group_df,
                       id.vars = "species",
                       measure.vars = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
                       variable.name = "Months",
                       value.name = "abundance")
species_melted$species <- factor(species_melted$species, levels = rev(x))
species_melted$seasons <- factor(species_melted$Months, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
p9 <- ggplot(species_melted, aes(seasons, abundance, alluvium = species)) +
  geom_alluvium(aes(fill = species), alpha = 1, linewidth = 0) +
  scale_fill_manual(
    values = setNames(c( "#d1e6c1", "#Face37",  "#E15759", "#0f4f94", "#666666", "#F28E2B",
                         "#B07AA1", "#FF9DA7", "#8c564b",  "#BAB0AC", "#66A61E","#B6992D",
                         "#D37295",  "#B2DF8A", "#9467bd", "#5ab1ef", "#F2F3C7","#4E79A7",
                         "#C4A38a", "#86BCB6", "#f5e8f1"), rev(x)),
    breaks = rev(x),
    labels = rev(x)
  ) +
  theme(
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14))+
  guides(fill = guide_legend(ncol = 2))
p9
ap_data <- read.delim("C:/Users/xuzy0/Desktop/xzy/split/beiyong/ap.txt", stringsAsFactors = FALSE)
ap_group <- data.frame(row.names = ap_data$Species,
                       Jul = rowSums(ap_data[,2:11], na.rm = TRUE),
                       Aug = rowSums(ap_data[,12:21], na.rm = TRUE),
                       Sep = rowSums(ap_data[,22:26], na.rm = TRUE),
                       Oct = rowSums(ap_data[,27:40], na.rm = TRUE),
                       Nov = rowSums(ap_data[,41:51], na.rm = TRUE),
                       Dec = rowSums(ap_data[,52:56], na.rm = TRUE),
                       Jan = rowSums(ap_data[,57:66], na.rm = TRUE),
                       Feb = rowSums(ap_data[,67:71], na.rm = TRUE),
                       Mar = rowSums(ap_data[,72:76], na.rm = TRUE),
                       Apr = rowSums(ap_data[,77:86], na.rm = TRUE),
                       May = rowSums(ap_data[,87:96], na.rm = TRUE),
                       Jun = rowSums(ap_data[,97:101], na.rm = TRUE))
x_ap <- rownames(ap_group)
ap_group <- ap_group[order(rowSums(ap_group), decreasing = T), ]
x_ap <- rownames(ap_group)
ap_group <- cbind(rownames(ap_group), ap_group)
ap_group <- melt(ap_group)
colnames(ap_group) <- c("species", "months", "abundance")
ap_group$species <- factor(ap_group$species, levels = rev(x_ap))
ap_group$months <- factor(ap_group$months, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                                      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
p_ap <- ggplot(ap_group, aes(months, abundance, fill = species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#d1e6c1", "#Face37", "#E15759", "#76B7B2", "#666666", "#F28E2B",
                               "#B07AA1", "#FF9DA7", "#8c564b", "#BAB0AC", "#66A61E", "#B6992D",
                               "#D37295", "#B2DF8A", "#9467bd", "#5ab1ef", "#F2F3C7", "#4E79A7")) +
  labs(x = "Month", y = "Abundance", fill = "Species") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  guides(fill = guide_legend(ncol = 2))
p_ap
