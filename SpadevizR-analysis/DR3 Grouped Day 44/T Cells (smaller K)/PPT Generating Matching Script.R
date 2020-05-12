#install.packages('officer') # Install
library('officer') # Load
library(magrittr)

#install.packages("magick")
library("magick")

#install.packages("DT")
library("DT")



dir <- getwd()
doc <- read_pptx() 
matching_data = read.csv("DR3 D44 T Cells Pearsons Coefficient.csv", header = T, sep = ",", stringsAsFactors = F)
matching_data = read.table("DR3 D44 T Cells smaller K Sorted Pearsons Data.txt", header = T, sep = "\t", stringsAsFactors = F)
rownumber = 0

for(i in 1:826){
 #   i = 1
  if(matching_data$experiment1_count[i] < 2000 | matching_data$experiment2_count[i] < 2000){
    next()
  } else{rownumber[i] = i}
  
  
  Grp1_cluster_number = substr(matching_data$experiment1_cluster[i], 14, 18)
  Grp2_cluster_number = substr(matching_data$experiment2_cluster[i], 14, 17)

  Grp1_cluster_file = paste(dir, "/Group1_ClusterImages/Cluster ", Grp1_cluster_number, ".jpeg", sep = "")
  Grp2_cluster_file = paste(dir, "/Group2_ClusterImages/Cluster ", Grp2_cluster_number, ".jpeg", sep = "")
  
  Grp1_scatter_file = paste(dir, "/Group1_Scatterplots/Cluster_", Grp1_cluster_number, ".jpeg", sep = "")
  Grp2_scatter_file = paste(dir, "/Group2_Scatterplots/Cluster_", Grp2_cluster_number, ".jpeg", sep = "")
  
  Grp1_cluster_file_magick = image_read(Grp1_cluster_file)
  #image_info(Grp1_cluster_file_magick)
  #Grp1_cluster_file_magick = image_crop(Grp1_cluster_file_magick, "2000x900")
  image_write(Grp1_cluster_file_magick, path = "Grp1_cluster_file.jpeg", format = "jpeg")
  
  Grp1_scatter_file_magick = image_read(Grp1_scatter_file)
  #image_info(Grp1_cluster_file_magick)
  image_write(Grp1_scatter_file_magick, path = "Grp1_scatter_file.jpeg", format = "jpeg")
  
  
  Grp2_cluster_file_magick = image_read(Grp2_cluster_file)
  #image_info(Grp1_cluster_file_magick)
  #Grp2_cluster_file_magick = image_crop(Grp2_cluster_file_magick, "2000x900")
  image_write(Grp2_cluster_file_magick, path = "Grp2_cluster_file.jpeg", format = "jpeg")
  
  Grp2_scatter_file_magick = image_read(Grp2_scatter_file)
  #image_info(Grp1_cluster_file_magick)
  image_write(Grp2_scatter_file_magick, path = "Grp2_scatter_file.jpeg", format = "jpeg")
  
  Grp1_cluster_file = paste(dir, "/Grp1_cluster_file.jpeg", sep = "")
  Grp2_cluster_file = paste(dir, "/Grp2_cluster_file.jpeg", sep = "")
  
  Grp1_scatter_file = paste(dir, "/Grp1_scatter_file.jpeg", sep = "")
  Grp2_scatter_file = paste(dir, "/Grp2_scatter_file.jpeg", sep = "")
  
  slide_title = paste(matching_data$experiment1_cluster[i], " VS ", matching_data$experiment2_cluster[i], sep = "")
  
  doc <- doc %>%
    add_slide(layout = "Two Content", master = "Office Theme") %>%
    #ph_with_text(type = "title", str = slide_title) %>%
    #ph_with_img(type = "body", str = "body (index 1) is text", index = 1) %>% 
    # ph_with_img(type = "body", index = 1, src = Grp1_cluster_file, height = 3.5, width = 4.58 ) %>%
    # ph_with_img(type = "body", index = 2, src = Grp2_cluster_file, height = 3.5, width = 4.58 )
  ph_with_img_at(src = Grp1_cluster_file, height = 3.4, width = 5.1, left = 0.1, top = 0) %>%
  ph_with_img_at(src = Grp2_cluster_file, height = 3.4, width = 5.1, left = 0.1, top = 3.4)%>%
  
  ph_with_img_at(src = Grp1_scatter_file, height = 2, width = 2, left = 5.32, top = 0.55) %>%
  ph_with_img_at(src = Grp2_scatter_file, height = 2, width = 2, left = 5.32, top = 3.75)
  
  file.remove(Grp1_cluster_file, Grp2_cluster_file, Grp1_scatter_file, Grp2_scatter_file)
  print(i)
}


print(doc, target = "DR3 D44 T Cells 99000per Matching Clusters SORTED all slides.pptx")


