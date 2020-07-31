
library(circlize)
library(ComplexHeatmap)
library(grid)
library(gridBase)
library(RColorBrewer)
#Prepare the structure of circos plot
circos_plot_stucture = data.frame(group = 0, section = 0, x_value = 0, section_1 = 0)
for(i in 1:(31*5)){
  circos_plot_stucture = rbind(circos_plot_stucture, 0)
}
circos_plot_stucture = circos_plot_stucture[-1,]

circos_plot_stucture$x_value = rep(1:5, 31)
j = 1
for (i in 1:nrow(circos_plot_stucture)) {
  if((i-1) %% 5 == 0 & i != 1){
    j = j+1
  }
  circos_plot_stucture$section[i] = j
}

for (i in 1:nrow(circos_plot_stucture)) {
  if (circos_plot_stucture$section[i] %in% c(1:8)){
    circos_plot_stucture$group[i] = "Coxevac Unchallenge"
  }
  if (circos_plot_stucture$section[i] %in% c(9:16)){
    circos_plot_stucture$group[i] = "Naive Unchallenge"
  }
  if (circos_plot_stucture$section[i] %in% c(17:23)){
    circos_plot_stucture$group[i] = "Naive Challenge"
  }
  if (circos_plot_stucture$section[i] %in% c(24:31)){
    circos_plot_stucture$group[i] = "Coxevac Challenge"
  }
}

mice_8 = 0
mice_7 = 0

for(i in 1:8){
  x = rep(i,5)
  mice_8 = c(mice_8,x)
}
mice_8 = mice_8[-1]

for (i in 1:7) {
  x = rep(i,5)
  mice_7 = c(mice_7,x)
}
mice_7 = mice_7[-1]

circos_plot_stucture$section_1 = c(mice_8, mice_8, mice_7, mice_8)

circos_label =read.table("labeling.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

write.table(circos_plot_stucture, "structure.txt", row.names = F, sep = "\t")
#Prepare the structure of circos plot

library(circlize)

# #circos_plot_stucture = read.table("structure.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# #heatmap_data = read.table("heatmap_3.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# raw_heatmap_data = read.table("Heat_Day51 Grouped T_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# #raw_heatmap_data = read.table("Heat_Day51 Grouped Live_4clusters.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# heatmap_data = data.frame(tmpname = 0)
# for(i in 1:(ncol(raw_heatmap_data)-1)){
#   heatmap_data = cbind(heatmap_data, 0)
# }
# colnames(heatmap_data) = colnames(raw_heatmap_data)
# heatmap_data_structure = heatmap_data
# 
# #i = 1
# for (i in 1:nrow(raw_heatmap_data)) {
#   tmp_heatmap_data = heatmap_data_structure
#   for(j in 1:4){
#     if(j == 1){
#       tmp_heatmap_data[j,] = raw_heatmap_data[i,]
#     }else{
#       tmp_heatmap_data[j,] = raw_heatmap_data[i,]
#       tmp_heatmap_data$start[j] = tmp_heatmap_data$end[j-1]
#       tmp_heatmap_data$end[j] = tmp_heatmap_data$start[j] + 1
#     }
#     
#     if(j < 4){
#       tmp_heatmap_data = rbind(tmp_heatmap_data, 0)
#     }
#   }
#   heatmap_data = rbind(heatmap_data, tmp_heatmap_data)
# }
# heatmap_data = heatmap_data[-1,]
# 
# tmp_rescale <- function(x) (x-min(x))/(max(x) - min(x))
# heatmap_data_1 = heatmap_data
# for (i in 4:ncol(heatmap_data_1)) {
#   heatmap_data_1[,i] = tmp_rescale(heatmap_data_1[,i])
# }



heatmap_data_dealing = function(heatmap_data_file){
  #heatmap_data_file = "heatmap_D10.txt"
  
  heatmap_data = read.table(heatmap_data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  tmp_rescale <- function(x) (x-min(x))/(max(x) - min(x))
  
  for(j in 4:ncol(heatmap_data)){
    heatmap_data[,j] = tmp_rescale(heatmap_data[,j])
  }
  
  heatmap_data_1 <- heatmap_data
  heatmap_data_1 <- heatmap_data_1[which(rownames(heatmap_data_1) == 1),]
  heatmap_data_1[1,] <- NA
  
  chr <- unique(heatmap_data$chr)
  for(i in chr){
    #i = 1
    tmp.heatmap.data <- heatmap_data[heatmap_data$chr == i, ]
    
    for(j in 2:4){
      #j = 2
      tmp.heatmap.data = rbind(tmp.heatmap.data, NA)
      tmp.heatmap.data[j, ] = tmp.heatmap.data[j - 1, ]
      tmp.heatmap.data$start[j] = tmp.heatmap.data$start[j-1] + 1
      tmp.heatmap.data$end[j] = tmp.heatmap.data$end[j-1] + 1
    }
          
    heatmap_data_1 <- rbind(heatmap_data_1, tmp.heatmap.data)
  }
  
  heatmap_data_1 <- na.omit(heatmap_data_1)
  heatmap_data_1 <- heatmap_data_1[order(heatmap_data_1$chr, heatmap_data_1$start, decreasing = F), ]

  return(heatmap_data_1)
  # heatmap_data_1 = data.frame(tmp_name = 0)
  # for(i in 1: (ncol(heatmap_data)-1)){
  #   heatmap_data_1 = cbind(heatmap_data_1, 0)
  # }
  # colnames(heatmap_data_1) = colnames(heatmap_data)
  # for(i in unique(heatmap_data$chr)){
  #   #i = "Coxevac Unchallenge"
  #   tmp.data = heatmap_data[heatmap_data$chr == i, ]
  #   for(j in 4:ncol(tmp.data)){
  #     tmp.data[,j] = tmp_rescale(tmp.data[,j])
  #   }
  #   heatmap_data_1 = rbind(heatmap_data_1, tmp.data)
  # }
  # heatmap_data_1 = heatmap_data_1[-1,]
  
  
  # for (i in 4:ncol(heatmap_data_1)) {
  #   heatmap_data_1[,i] = tmp_rescale(heatmap_data_1[,i])
  # }
  
  
}

D10_data = heatmap_data_dealing(heatmap_data_file = "pos_cor.txt")
D35_data = heatmap_data_dealing(heatmap_data_file = "neg_cor.txt")
#D44_data = heatmap_data_dealing(heatmap_data_file = "heatmap_D44.txt")
#D51_data = heatmap_data_dealing(heatmap_data_file = "heatmap_D51.txt")

antibod_day10 = read.table("antibody_day10_2.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
antibod_day24 = read.table("antibody_day24_2.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
antibod_day35 = read.table("antibody_day35_2.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
antibody_data = list(antibod_day10, antibod_day24, antibod_day35)

spleen = read.table("spleen_2.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
heart = read.table("heart_2.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
liver = read.table("liver_2.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
lung = read.table("lung_2.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
pathology_score = list(heart, spleen, lung, liver)

body_weight = read.table("body_weight_3.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)

circos_plot = function(){
  # General parameters
  circos.par(#"track.height" = 0.1, 
    "cell.padding" = c(0, 0, 0, 0))
  
  # Initialise the chart giving factor and x-axis.
  fa = circos_plot_stucture$section
  f1 = factor(fa, levels = unique(fa))
  circos.initialize( factors= f1, x = circos_plot_stucture$x_value)
  
  # #Heatmap_1
  # f = colorRamp2(breaks = c(0, 0.5, 1), colors = c("#663366", "#33ffcc", "#990000"))
  # circos.genomicTrackPlotRegion(heatmap_data_1, stack = TRUE, 
  #                               panel.fun = function(region, value, ...) {
  #                                 #circos.axis(major.at = c(1:9))
  #                                 circos.genomicRect(region, value, col = f(value[[1]]), border = "black", lwd = 0.01, posTransform = posTransform.default, ...)
  #                               }, 
  #                               bg.border = NA, track.height = 0.5)
  
  
  # for (i in 1: nrow(circos_label)){
  #   highlight.sector(sector.index = circos_label$section[i], track.index = 1,
  #                    text = circos_label$label[i], text.vjust = "5mm", niceFacing = TRUE, col = "transparent", cex = 1
  #   )
  # }
  # 
  # highlight.sector(c("1":"8"), track.index = 1, text = "Coxevac Unchallenge", col = "transparent",
  #                  facing = "bending.inside", niceFacing = TRUE, text.vjust = "9mm", cex = 1)
  # highlight.sector(c("9":"16"), track.index = 1, text = "Coxevac Challenge", col = "transparent",
  #                  facing = "bending.inside", niceFacing = TRUE, text.vjust = "9mm", cex = 1)
  # highlight.sector(c("17":"23"), track.index = 1, text = "Naive Challenge", col = "transparent",
  #                  facing = "bending.inside", niceFacing = TRUE, text.vjust = "9mm", cex = 1)
  # highlight.sector(c("24":"31"), track.index = 1, text = "Naive Unchallenge", col = "transparent",
  #                  facing = "bending.inside", niceFacing = TRUE, text.vjust = "9mm", cex = 1)
  #?brewer.pal
  #Heatmap D10
  #display.brewer.pal(9, "YlOrRd")
  f = colorRamp2(breaks = c(0, 0.25, 0.5, 0.75, 1), colors = brewer.pal(5,"Oranges"))
  circos.genomicTrackPlotRegion(D10_data, stack = TRUE, 
                                panel.fun = function(region, value, ...) {
                                  #circos.axis(major.at = c(1:8))
                                  circos.genomicRect(region, value, col = f(value[[1]]), border = f(value[[1]]), lwd = 0.001, posTransform = posTransform.default, ...)
                                }, 
                                bg.border = "black", track.height = 0.1)
  #?circos.genomicTrackPlotRegion
  
  #labeling
  for (i in 1: nrow(circos_label)){
    highlight.sector(sector.index = circos_label$section[i], track.index = 1,
                     text = circos_label$label[i], text.vjust = "5mm", niceFacing = TRUE, col = "transparent", cex = 1
    )
  }
  
  highlight.sector(c("1":"8"), track.index = 1, text = "Coxevac Unchallenge", col = "transparent",
                   facing = "bending.inside", niceFacing = TRUE, text.vjust = "9mm", cex = 1)
  highlight.sector(c("9":"16"), track.index = 1, text = "Naive Unchallenge", col = "transparent",
                   facing = "bending.inside", niceFacing = TRUE, text.vjust = "9mm", cex = 1)
  highlight.sector(c("17":"23"), track.index = 1, text = "Naive Challenge", col = "transparent",
                   facing = "bending.inside", niceFacing = TRUE, text.vjust = "9mm", cex = 1)
  highlight.sector(c("24":"31"), track.index = 1, text = "Coxevac Challenge", col = "transparent",
                   facing = "bending.inside", niceFacing = TRUE, text.vjust = "9mm", cex = 1)
  
  
  #Heatmap D35
  f = colorRamp2(breaks = c(0, 0.25, 0.5, 0.75, 1), colors = brewer.pal(5,"Blues"))
  circos.genomicTrackPlotRegion(D35_data, stack = TRUE, 
                                panel.fun = function(region, value, ...) {
                                  #circos.axis(major.at = c(1:8))
                                  circos.genomicRect(region, value, col = f(value[[1]]), border = f(value[[1]]), lwd = 0.01, posTransform = posTransform.default, ...)
                                }, 
                                bg.border = "black", track.height = 0.1)
  
  # #Heatmap D44
  # f = colorRamp2(breaks = c(0, 0.25, 0.5, 0.75, 1), colors = brewer.pal(5,"Blues"))
  # circos.genomicTrackPlotRegion(D44_data, stack = TRUE, 
  #                               panel.fun = function(region, value, ...) {
  #                                 #circos.axis(major.at = c(1:8))
  #                                 circos.genomicRect(region, value, col = f(value[[1]]), border = f(value[[1]]), lwd = 0.01, posTransform = posTransform.default, ...)
  #                               }, 
  #                               bg.border = "black", track.height = 0.1)
  # 
  # #Heatmap D51
  # f = colorRamp2(breaks = c(0, 0.25, 0.5, 0.75, 1), colors = brewer.pal(5,"Reds"))
  # circos.genomicTrackPlotRegion(D51_data, stack = TRUE, 
  #                               panel.fun = function(region, value, ...) {
  #                                 #circos.axis(major.at = c(1:8))
  #                                 circos.genomicRect(region, value, col = f(value[[1]]), border = f(value[[1]]), lwd = 0.01, posTransform = posTransform.default, ...)
  #                               }, 
  #                               bg.border = "black", track.height = 0.1)
  
  # 
  # 
  # 
  # #Heatmap_3
  # f = colorRamp2(breaks = c(0, 0.5, 1), colors = c("#4daf4a", "#377eb8", "#e41a1c"))
  # circos.genomicTrackPlotRegion(heatmap_data_1, stack = TRUE, 
  #                               panel.fun = function(region, value, ...) {
  #                                 #circos.axis(major.at = c(1:8))
  #                                 circos.genomicRect(region, value, col = f(value[[1]]), border = "black", lwd = 0.2, posTransform = posTransform.default, ...)
  #                               }, 
  #                               bg.border = NA, track.height = 0.3)
  # #Heatmap_4
  # f = colorRamp2(breaks = c(0, 0.5, 1), colors = c("#7570b3", "#d95f02", "#1b9e77"))
  # circos.genomicTrackPlotRegion(heatmap_data_1, stack = TRUE, 
  #                               panel.fun = function(region, value, ...) {
  #                                 #circos.axis(major.at = c(1:8))
  #                                 circos.genomicRect(region, value, col = f(value[[1]]), border = "black", lwd = 0.2, posTransform = posTransform.default, ...)
  #                               }, 
  #                               bg.border = NA, track.height = 0.1)
  
  #Antibody
  col = c("#99cc33", "#ff6666", "#336699")
  circos.genomicTrackPlotRegion(antibody_data, ylim = c(0, 1.5), panel.fun = function(region, value, ...) {
    i = getI(...)
    circos.genomicPoints(region, value, col = col[i], cex = 0.7, pch = 16)
  }, track.height = 0.1)
  
  #pathology_score
  col = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
  circos.genomicTrack(pathology_score, ylim = c(0, 2.3), #stack = TRUE,
                      panel.fun = function(region, value, ...) {
                        i = getI(...)
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
                                           col = col[i], ...)
                      }, track.height = 0.1)
  
  #Body weitght
  min.bacterial.burden = min(log10(body_weight$bacteria_burden_1))
  max.bacterial.burden = max(log10(body_weight$bacteria_burden_1))
  bg.colour = colorRamp2(breaks = c(0, 1.5, 3, 4.5, 6.5), 
                         colors = c("#def0db", "#e6f49e", "#fddf92", "#faae69", "#d14252"))
  circos.trackPlotRegion(ylim = c(0, ceiling(max(body_weight$Body_weight))),
                         # bg.col = c("1" = "#fb8072"
                         #            #"Naive Challenge" = "#bebada", 
                         #            #"Naive Unchallenge" = "#fb8072", 
                         #            #"Coxevac Unchallenge" =  "#8dd3c7", 
                         #            #"Coxevac Challenge" = "#ffffb3"
                         #            ),
                         # bg.col = paste("#", body_weight$bg_colour, sep = ""),
                         bg.col = bg.colour(log10(body_weight$bacteria_burden_1)),
                         #bg.col = "#bebada",
                         factors = body_weight$chr, x=body_weight$start, 
                         y = body_weight$Body_weight, panel.fun = function(x, y) {
                           #circos.axis(major.at = c(1:8))
                           #bg.col = "red"
                         }, track.height = 0.1)
  circos.trackPoints(body_weight$chr, body_weight$end, body_weight$Body_weight, col = "#1a9641", pch = 16, cex = 0.8)
}


####legends#####
heatmap.T.cell.col.fun = colorRamp2(breaks = c(0, 0.25, 0.5, 0.75, 1), colors = brewer.pal(5,"Oranges"))
heatmap.T.cell <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), 
                         col_fun = heatmap.T.cell.col.fun, 
                         title_position = "topleft", 
                         title = "Track1-Positive Correlation")

heatmap.B.cell.col.fun = colorRamp2(breaks = c(0, 0.25, 0.5, 0.75, 1), colors = brewer.pal(5,"Blues"))
heatmap.B.cell <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), 
                         col_fun = heatmap.B.cell.col.fun, 
                         title_position = "topleft", 
                         title = "Track2-Negative Correlation")

# heatmap.innate.cell.col.fun = colorRamp2(breaks = c(0, 0.25, 0.5, 0.75, 1), colors = brewer.pal(5,"Blues"))
# heatmap.innate.cell <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), 
#                          col_fun = heatmap.innate.cell.col.fun, 
#                          title_position = "topleft", 
#                          title = "Track3-Day 44")
# 
# heatmap.D51.col.fun = colorRamp2(breaks = c(0, 0.25, 0.5, 0.75, 1), colors = brewer.pal(5,"Reds"))
# heatmap.D51.cell <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), 
#                            col_fun = heatmap.D51.col.fun, 
#                            title_position = "topleft", 
#                            title = "Track4-Day 51")

antibodys = Legend(at = c("Day10", "Day24", "Day35"), 
                   type = "points", 
                   legend_gp = gpar(col = c("#99cc33", "#ff6666", "#336699")), 
                   title_position = "topleft", 
                   title = "Track3-Antibody")

pathology.score = Legend(at = c("Heart", "Spleen", "Lung", "Liver"), 
                         type = "grid", 
                         legend_gp = gpar(fill = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")), 
                         title_position = "topleft", 
                         title = "Track4-Pathology")

body.weight = Legend(at = c("Body weight"), 
                     type = "points", 
                     legend_gp = gpar(col = c("#1a9641")), 
                     title_position = "topleft", 
                     title = "Track5-Body weight")

#min.bacterial.burden = min(log10(body_weight$bacteria_burden_1))
#max.bacterial.burden = max(log10(body_weight$bacteria_burden_1))
bg.colour = colorRamp2(breaks = c(0, 1.5, 3, 4.5, 6), 
                       colors = c("#def0db", "#e6f49e", "#fddf92", "#faae69", "#d14252"))
bacterial.burden = Legend(at = c(0, 1.5, 3, 4.5, 6), 
                          col_fun = bg.colour, 
                          title_position = "topleft", 
                          title = "Track5-Burden(log10.value)")

lgd.list <- packLegend(heatmap.T.cell, 
                       heatmap.B.cell, 
                       # heatmap.innate.cell, 
                       # heatmap.D51.cell,
                       antibodys, 
                       pathology.score, 
                       body.weight, 
                       bacterial.burden)


grid.newpage()
par(new=TRUE)
plot.new()

circle_size = unit(1, "snpc") # snpc unit gives you a square region

pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circos_plot()
upViewport()

draw(lgd.list, x = circle_size+unit(110, "mm"), just = "left")

