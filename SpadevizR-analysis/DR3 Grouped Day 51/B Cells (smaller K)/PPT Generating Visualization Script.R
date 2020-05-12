

install.packages('officer') # Install
library('officer') # Load

install.packages('magittr')
library(magrittr)

install.packages("magick")
library("magick")

install.packages("DT")
library("DT")

#datatable(matching_data)

dir <- getwd()
doc <- read_pptx() 
#matching_data = read.csv("DR3 D51 T Cells Pearsons Coefficient.csv", header = T, sep = ",", stringsAsFactors = F)
matching_file_data = read.table("Matching Clusters for Visualization.txt", header = T, sep = "\t", stringsAsFactors = F)


screened_cluster_number = 0



for(i in 1:31){
# i=1
  #SETS OBJECT AS CLUSTER OF INTEREST
  ## make sure to change character count to match cluster number digits!!
  Grp1_cluster_number = substr(matching_file_data$Experiment1_Cluster[i], 1, 5)
  Grp2_cluster_number = substr(matching_file_data$Experiment2_Cluster[i], 1, 5)

  
  if(Grp1_cluster_number %in% screened_cluster_number & i > 1){
    next()
  }
  
  
  Grouped_cluster_number <- matching_file_data$new_name[match(as.numeric(Grp1_cluster_number), matching_file_data$Experiment1_Cluster)]
  
  
  if(is.na(Grouped_cluster_number)){
    
    print(paste("The cluster ", Grp1_cluster_number, " doesn't exist in matching file.", sep = ""))
    
    next()
    
  }
  
  screened_cluster_number[i] = Grp1_cluster_number
  
  
  #PASTES PARALLEL COORDINATE PLOT OF CLUSTERS OF INTEREST
  Grp1_cluster_file = paste(dir, "/Group1_ClusterImages/Cluster ", Grp1_cluster_number, ".jpeg", sep = "")
  Grp2_cluster_file = paste(dir, "/Group2_ClusterImages/Cluster ", Grp2_cluster_number, ".jpeg", sep = "")
  
  Grouped_cluster_file = paste(dir, "/Grouped_ClusterImages/Cluster_ ", Grouped_cluster_number, " .jpeg", sep = "")
  
    
  Volcano_Plot_CvNC = paste(dir, "/SpadevizR Analysis and Volcano Plots/Control vs Naive_chal", ".tiff", sep = "")
  Volcano_Plot_CvVU = paste(dir, "/SpadevizR Analysis and Volcano Plots/Control vs Vax_Unchal", ".tiff", sep = "")
  Volcano_Plot_CvVC = paste(dir, "/SpadevizR Analysis and Volcano Plots/Control vs Vax_Chal", ".tiff", sep = "")
  Volcano_Plot_NCvVU = paste(dir, "/SpadevizR Analysis and Volcano Plots/Naive_Chal vs Vax_Unchal", ".tiff", sep = "")
  Volcano_Plot_NCvVC = paste(dir, "/SpadevizR Analysis and Volcano Plots/Naive_Chal vs Vax_Chal", ".tiff", sep = "")
  Volcano_Plot_VUvVC = paste(dir, "/SpadevizR Analysis and Volcano Plots/Vax_Unchal vs Vax_Chal", ".tiff", sep = "")
  
  #PASTES LEGEND OF PARALLEL COORDINATE PLOTS OF CLUSTERS OF INTEREST
  #Legend_file = paste(dir, "/Group1_PhenoViewer/Cluster ", Grp1_cluster_number, ".jpeg", sep = "")
  
  #######   #######
  #########################
  ## NEED HELP HERE !!! ###
  #########################
  #######  want to set legend file = text or image of markers from PCP (from either group 1 or 2) #######
  
  #PASTES SCATTERPLOTS OF CLUSTERS OF INTEREST
  Grp1_scatter_file = paste(dir, "/Group1_Scatterplots/Cluster_", Grp1_cluster_number, ".jpeg", sep = "")
  Grp2_scatter_file = paste(dir, "/Group2_Scatterplots/Cluster_", Grp2_cluster_number, ".jpeg", sep = "")
  Grouped_scatter_file = paste(dir, "/Grouped_Scatterplots/V", Grouped_cluster_number, ".jpeg", sep = "")

  
  Volcano_Plot_magick_CvNC = image_read(Volcano_Plot_CvNC)
  Volcano_Plot_magick_CvVU = image_read(Volcano_Plot_CvVU)
  Volcano_Plot_magick_CvVC = image_read(Volcano_Plot_CvVC)
  Volcano_Plot_magick_NCvVU = image_read(Volcano_Plot_NCvVU)
  Volcano_Plot_magick_VUvVC = image_read(Volcano_Plot_VUvVC)
  Volcano_Plot_magick_NCvVC = image_read(Volcano_Plot_NCvVC)
  
  image_write(Volcano_Plot_magick_CvNC, path = "Control vs Naive_chal.jpeg", format = "jpeg")
  image_write(Volcano_Plot_magick_CvVU, path = "Control vs Vax_Unchal.jpeg", format = "jpeg")
  image_write(Volcano_Plot_magick_CvVC, path = "Control vs Vax_Chal.jpeg", format = "jpeg")
  image_write(Volcano_Plot_magick_NCvVU, path = "Naive_Chal vs Vax_Unchal.jpeg", format = "jpeg")
  image_write(Volcano_Plot_magick_VUvVC, path = "Vax_Unchal vs Vax_Chal.jpeg", format = "jpeg")
  image_write(Volcano_Plot_magick_NCvVC, path = "Naive_Chal vs Vax_Chal.jpeg", format = "jpeg")
  
  
  #######   #######
  #########################
  ## NEED HELP HERE !!! ###
  #########################
  #######  How to get PCP legend ?  #######
  
  #CROPS PCP FOR LEGEND
  # Legend_file_magick = image_read(Legend_file)
  # Legend_file_magick = image_crop(Legend_file_magick, "2000x500")
  # image_write(Legend_file_magick, path = "Legend_file.jpeg", format = "jpeg")


  #################
  
  #SET RESOLUTION OF PARALLEL COORDINATE PLOTS
  Grp1_cluster_file_magick = image_read(Grp1_cluster_file)
  Grp1_cluster_file_magick = image_crop(Grp1_cluster_file_magick, "2000x900")
  image_write(Grp1_cluster_file_magick, path = "Grp1_cluster_file.jpeg", format = "jpeg")
  
  Grp1_scatter_file_magick = image_read(Grp1_scatter_file)
  image_write(Grp1_scatter_file_magick, path = "Grp1_scatter_file.jpeg", format = "jpeg")
  
  Grp2_cluster_file_magick = image_read(Grp2_cluster_file)
  Grp2_cluster_file_magick = image_crop(Grp2_cluster_file_magick, "2000x900")
  image_write(Grp2_cluster_file_magick, path = "Grp2_cluster_file.jpeg", format = "jpeg")
  
  Grp2_scatter_file_magick = image_read(Grp2_scatter_file)
  image_write(Grp2_scatter_file_magick, path = "Grp2_scatter_file.jpeg", format = "jpeg")
  
  Grouped_cluster_file_magick = image_read(Grouped_cluster_file)
  Grouped_cluster_file_magick = image_crop(Grouped_cluster_file_magick, "2000x900")
  image_write(Grouped_cluster_file_magick, path = "Grouped_cluster_file.jpeg", format = "jpeg")
  
  Grouped_scatter_file_magick = image_read(Grouped_scatter_file)
  image_write(Grouped_scatter_file_magick, path = "Grouped_scatter_file.jpeg", format = "jpeg")

  #Legend_file = paste(dir, "/Grp1_cluster_file.jpeg", sep = "")
  
  Volcano_Plot_CvNC = paste(dir, "/Control vs Naive_chal.jpeg", sep = "")
  Volcano_Plot_CvVU = paste(dir, "/Control vs Vax_Unchal.jpeg", sep = "")
  Volcano_Plot_CvVC = paste(dir, "/Control vs Vax_Chal.jpeg", sep = "")
  Volcano_Plot_NCvVU = paste(dir, "/Naive_Chal vs Vax_Unchal.jpeg", sep = "")
  Volcano_Plot_VUvVC = paste(dir, "/Vax_Unchal vs Vax_Chal.jpeg", sep = "")
  Volcano_Plot_NCvVC = paste(dir, "/Naive_Chal vs Vax_Chal.jpeg", sep = "")
  
  
  Grp1_cluster_file = paste(dir, "/Grp1_cluster_file.jpeg", sep = "")
  Grp2_cluster_file = paste(dir, "/Grp2_cluster_file.jpeg", sep = "")
  
  Grp1_scatter_file = paste(dir, "/Grp1_scatter_file.jpeg", sep = "")
  Grp2_scatter_file = paste(dir, "/Grp2_scatter_file.jpeg", sep = "")
  
  Grouped_cluster_file = paste(dir, "/Grouped_cluster_file.jpeg", sep = "")
  Grouped_scatter_file = paste(dir, "/Grouped_scatter_file.jpeg", sep = "")
  
  #SETS SLIDE TITLE
  slide_title = paste(matching_file_data$Experiment1_Cluster[i], " VS ", matching_file_data$Experiment2_Cluster[i], sep = "")

  
  #SET POWERPOINT THEME
  doc <- doc %>%
    add_slide(layout = "Two Content", master = "Office Theme") %>%
    #ph_with_text(type = "title", str = slide_title) %>%
    #ph_with_img(type = "body", str = "body (index 1) is text", index = 1) %>% 
    
  #LEGEND SIZE AND LOCATION
  #ph_with_img_at(src = Legend_file, height = 2.21, width = 4.9, left = 0.1, top = 5.11)%>%

  #PARALLEL COORDINATE PLOTS SIZE AND LOCATION
  ph_with_img_at(src = Grouped_cluster_file, height = 2.04, width = 4.35, left = 0, top = 0) %>%
  ph_with_img_at(src = Grp1_cluster_file, height = 2.04, width = 4.35, left = 0, top = 2.73) %>%
  ph_with_img_at(src = Grp2_cluster_file, height = 2.04, width = 4.35, left = 0, top = 5)%>%
   
    
  #SCATTERPLOTS SIZE AND LOCATION
  ph_with_img_at(src = Grouped_scatter_file, height = 1.7, width = 1.7, left = 4.35, top = 0.5) %>%
  ph_with_img_at(src = Grp1_scatter_file, height = 1.7, width = 1.7, left = 4.35, top = 3.12) %>%
  ph_with_img_at(src = Grp2_scatter_file, height = 1.7, width = 1.7, left = 4.35, top = 5.45) %>% 
  
  
  #VOLCANO PLOTS SIZE AND LOCATION
  ph_with_img_at(src = Volcano_Plot_CvNC, height = 2.1, width = 1.22, left = 7.03, top = 0.0)%>%
  ph_with_img_at(src = Volcano_Plot_CvVU, height = 2.1, width = 1.22, left = 7.99, top = 0.0)%>%
    ph_with_img_at(src = Volcano_Plot_CvVC, height = 2.1, width = 1.22, left = 8.95, top = 0.0)%>%
    ph_with_img_at(src = Volcano_Plot_NCvVU, height = 2.1, width = 1.22, left = 7.03, top = 2.2)%>%
  ph_with_img_at(src = Volcano_Plot_NCvVC, height = 2.1, width = 1.22, left = 7.99, top = 2.2)%>%
  ph_with_img_at(src = Volcano_Plot_VUvVC, height = 2.1, width = 1.22, left = 8.95, top = 2.2)
  
  
  file.remove(Volcano_Plot_NCvVC, Volcano_Plot_VUvVC, Volcano_Plot_CvVC, Volcano_Plot_NCvVU, Volcano_Plot_CvNC, Volcano_Plot_CvVU, Grp1_cluster_file, Grp2_cluster_file, Grp1_scatter_file, Grp2_scatter_file, Grouped_cluster_file, Grouped_scatter_file)
  print(i)
}

#GENERATES POWERPOINT 
print(doc, target = "DR3 D51 B Cell 99000per 37p genbyR.pptx")


