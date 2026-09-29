library(readxl)
library(tidyverse)
library(Biostrings)
library(ggplot2)
library(zoo)
library(seqinr)
library(ggnewscale)
library(ggtree)
library("rjson")
library(ggmap)
library(ape)
library(lubridate)
library(geojsonsf)
library(ggpubr)
library(treeio)
library(circlize)
library("readxl")
library(ggh4x)
library(seqinr)
library(dplyr)
library("ggborderline")
# 
# setwd("/Users/shu/vjlab Dropbox/hu shu/flub-evo/analysis")
# setwd("/Users/ruopengxie/vjlab Dropbox/Ruopeng Xie/flub-evo/analysis")
# setwd("/Users/ruopeng/Library/CloudStorage/Dropbox-vjlab/Ruopeng Xie/flub-evo/analysis")


tree <- read.beast("../Beast/antigenic_NA/only_NI/antigenic/IBV_NA_antigenic.mean.mcc.tree")

(p <- ggtree(tree, mrsd = "2021-06-04") + 
    #geom_range(range='height_0.95_HPD', size=1, color="#BABABA",alpha=0.3)+
    theme_tree2()+
    #geom_rootedge(rootedge = 0.01, color = "#cc5252")+            #root length, set root (Goldfields-Esperance) color
    #scale_color_manual(values=location_color,breaks=node_location)+
    #geom_text2(aes(subset = !isTip, label=label),size = 1.5)+  #support value
    #geom_text(aes(label=round(as.numeric(location.prob), 2)))+
    #scale_x_continuous(limits = c(2001.7,2020.2),breaks = c(2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020))+
    scale_y_continuous(expand = c(0.01,0))
  #theme(axis.text.x = element_text(angle=25, hjust=1))
  #geom_hilight(node=1253, fill="steelblue", alpha=0.5) +
  #geom_hilight(node=1215, fill="#8852b8", alpha=0.5)
)

#antigenic

c.data.parent.node <- data.frame("node" = p$data$parent)
c.data <- data.frame("node" = p$data$node, 
                     "end_1"= p$data$antigenic1,
                     "end_2"= p$data$antigenic2,
                     "parent_node" = p$data$parent, 
                     "isTip" = p$data$isTip,
                     "strain" = p$data$label,
                     "height" = p$data$height,
                     row.names = NULL)
c.data.parent <- merge(c.data.parent.node,c.data,by = "node", all.x = TRUE) %>% 
  subset(,c("node","end_1","end_2","height")) %>% 
  dplyr::rename(parent_node="node",start_1="end_1",start_2="end_2",ancestor_height="height") %>% 
  unique()

diffusion.data <- merge(c.data,c.data.parent,by = "parent_node", all.x = TRUE)
diffusion.data$time <- 2021.423-diffusion.data$height
diffusion.data$ancestor_time <- 2021.423-diffusion.data$ancestor_height


lineage <- read_tsv("../Nextstrain/NA/data/IBV_NA_meta.tsv")

#antigenic dimension 1 over time
g <- ggplot()
location.list <- c()
for (i in 1:nrow(diffusion.data)) {
  location_1 <- c(diffusion.data[i,]$ancestor_time, diffusion.data[i,]$time)
  location_2 <- c(diffusion.data[i,]$start_1, diffusion.data[i,]$end_1)
  g <- g+geom_line(data= data.frame(location_1,location_2),aes(x= location_1,y= location_2))
  if(diffusion.data[i,]$isTip == TRUE){
    
    # remove the same location tips
    if(paste(diffusion.data[i,]$time,diffusion.data[i,]$end_1) %in% location.list){
      next
    } else{
      location.list <- rbind(location.list, paste(diffusion.data[i,]$time,diffusion.data[i,]$end_1))
    }
    
    if(diffusion.data[i,]$strain %in% lineage[which(lineage$lineage == "Y"),]$strain){
      g <- g+geom_point(data = diffusion.data[i,],aes(x=time,y=end_1),shape=16,stroke=0,size=3,color="#e41a1c",alpha=(50-diffusion.data[i,]$height)/50+0.05)
    } else if(diffusion.data[i,]$strain %in% lineage[which(lineage$lineage == "V"),]$strain){
      g <- g+geom_point(data = diffusion.data[i,],aes(x=time,y=end_1),shape=16,stroke=0,size=3,color="#377eb8",alpha=(50-diffusion.data[i,]$height)/50+0.05)
    } else{
      g <- g+geom_point(data = diffusion.data[i,],aes(x=time,y=end_1),shape=16,stroke=0,size=3,color="#636363",alpha=(max(diffusion.data$height)-diffusion.data[i,]$height)/max(diffusion.data$height)+0.1)
    }
  }
}

(g <- g + 
  xlab("year") +
  ylab("Antigenic dimension 1")+
  scale_y_continuous(breaks = seq(1,11,3))+
  scale_x_continuous(breaks = seq(1930,2024,10),guide = "axis_minor")+
  theme_classic())

ggsave("../figures/NA_antigenic_1_over_time.pdf",g,width=20,height=9,units="cm")


#antigenic dimension 2 over time
g <- ggplot()
location.list <- c()
for (i in 1:nrow(diffusion.data)) {
  location_1 <- c(diffusion.data[i,]$ancestor_time, diffusion.data[i,]$time)
  location_2 <- c(diffusion.data[i,]$start_2, diffusion.data[i,]$end_2)
  g <- g+geom_line(data= data.frame(location_1,location_2),aes(x= location_1,y= location_2))
  if(diffusion.data[i,]$isTip == TRUE){
    # remove the same location tips
    if(paste(diffusion.data[i,]$time,diffusion.data[i,]$end_2) %in% location.list){
      next
    } else{
      location.list <- rbind(location.list, paste(diffusion.data[i,]$time,diffusion.data[i,]$end_2))
    }
    
    if(diffusion.data[i,]$strain %in% lineage[which(lineage$lineage == "Y"),]$strain){
      g <- g+geom_point(data = diffusion.data[i,],aes(x=time,y=end_2),shape=16,stroke=0,size=3,color="#e41a1c",alpha=(50-diffusion.data[i,]$height)/50+0.05)
    } else if(diffusion.data[i,]$strain %in% lineage[which(lineage$lineage == "V"),]$strain){
      g <- g+geom_point(data = diffusion.data[i,],aes(x=time,y=end_2),shape=16,stroke=0,size=3,color="#377eb8",alpha=(50-diffusion.data[i,]$height)/50+0.05)
    } else{
      g <- g+geom_point(data = diffusion.data[i,],aes(x=time,y=end_2),shape=16,stroke=0,size=3,color="#636363",alpha=(max(diffusion.data$height)-diffusion.data[i,]$height)/max(diffusion.data$height)+0.1)
    }
  }
}
(g <- g + 
    xlab("year") +
    ylab("Antigenic dimension 2")+
    scale_x_continuous(breaks = seq(1930,2024,10),guide = "axis_minor")+
    theme_classic())

ggsave("../figures/NA_antigenic_2_over_time.pdf",g,width=20,height=9,units="cm")


#antigenic dimension 1 and 2

#df.cluster <- read_excel("../data/HA/MK_IBV_clusterIDs_v2.xlsx")

colors <- c("Anc40s" = "grey8",
            "Anc50-70" = "grey",
            "HK72" = "royalblue1",
            "Mal04" = "mediumblue",
            "Bris08 - V1A" = "deepskyblue",
            "7" = "dodgerblue4",
            "Col17 - V1A.1" = "slateblue",
            "Wash19 - V1A.3" = "darkviolet",
            "Austria21 - V1A.3a.2" = "purple4",
            "Yam90s" = "firebrick4",
            "Yam00s" = "orangered1",
            "Mass12 - Y2" = "goldenrod1",
            "Phu13 - Y3" = "red3")

g <- ggplot()
location.list <- c()
point.info <- data.frame()

for (i in 1:nrow(diffusion.data)) {
  location_1 <- c(diffusion.data[i,]$start_1, diffusion.data[i,]$end_1)
  location_2 <- c(diffusion.data[i,]$start_2, diffusion.data[i,]$end_2)
  g <- g+geom_line(data= data.frame(location_1,location_2),aes(x= location_1,y= location_2))
  #if(diffusion.data[i,]$isTip == TRUE){
  if(diffusion.data[i,]$isTip == TRUE){
    # if(diffusion.data[i,]$height < 0.65){
    #   print(diffusion.data[i,])
    #   g <- g+geom_point(data = diffusion.data[i,],aes(x=end_1,y=end_2),color="yellow",size=5,alpha=1)
    #   next
    # }
    # remove the same location tips
    # if(paste(diffusion.data[i,]$end_1,diffusion.data[i,]$end_2) %in% location.list){
    #   next
    # } else{
    #   location.list <- rbind(location.list, paste(diffusion.data[i,]$end_1,diffusion.data[i,]$end_2))
    # }
    # 
    # if (diffusion.data[i,]$strain %in% df.cluster$Isolate){
    #   g <- g+geom_point(data = diffusion.data[i,],aes(x=end_1,y=end_2),shape=16,stroke=0,size=3,
    #                     color=colors[df.cluster[which(df.cluster$Isolate == diffusion.data[i,]$strain),]$Cluster])
    # } else{
    #   g <- g+geom_point(data = diffusion.data[i,],aes(x=end_1,y=end_2),shape=16,stroke=0,size=3,
    #                     color="#BABABA",alpha=0.2)
    # }
    
    
    if(diffusion.data[i,]$strain %in% lineage[which(lineage$lineage == "Y"),]$strain){
      #g <- g+geom_point(data = diffusion.data[i,],aes(x=end_1,y=end_2),shape=16,stroke=0,size=3,color="#e41a1c",alpha=diffusion.data[i,]$height/50+0.2)
      g <- g+geom_point(data = diffusion.data[i,],aes(x=end_1,y=end_2),shape=16,stroke=0,size=3,color="#e41a1c",alpha=(50-diffusion.data[i,]$height)/50+0.05)
      
      point.alpha <- (50-diffusion.data[i,]$height)/50+0.05
      
    } else if(diffusion.data[i,]$strain %in% lineage[which(lineage$lineage == "V"),]$strain){
      #g <- g+geom_point(data = diffusion.data[i,],aes(x=end_1,y=end_2),shape=16,stroke=0,size=3,color="#377eb8",alpha=diffusion.data[i,]$height/50+0.2)
      g <- g+geom_point(data = diffusion.data[i,],aes(x=end_1,y=end_2),shape=16,stroke=0,size=3,color="#377eb8",alpha=(50-diffusion.data[i,]$height)/50+0.05)
      
      point.alpha <- (50-diffusion.data[i,]$height)/50+0.05
      
    } else{
      g <- g+geom_point(data = diffusion.data[i,],aes(x=end_1,y=end_2),shape=16,stroke=0,size=3,color="#636363",alpha=(max(diffusion.data$height)-diffusion.data[i,]$height)/max(diffusion.data$height)+0.1)
      
      point.alpha <- (max(diffusion.data$height)-diffusion.data[i,]$height)/max(diffusion.data$height)+0.1
    }
    
    # record information
    
    point.info <- rbind(
      point.info,
      data.frame(
        strain = diffusion.data[i,]$strain,
        alpha = point.alpha
      ))
  }
}
(g <- g + 
    xlab("Antigenic dimension 1") +
    ylab("Antigenic dimension 2")+
    theme_classic())

ggsave("../figures/NA_antigenic_1_2.pdf",g,width=20,height=9,units="cm")

write.csv(point.info, "../table/NA_antigenic_points_alpha.csv", row.names=FALSE)

#antigenic dimension 1 and 2 colored by year
library(scales)
g <- ggplot()
location.list <- c()
# create color function from 1940 to 2021.4
col_fun <- col_numeric(
  palette = colorRampPalette(c("#4daf4a", "#984ea3"))(220),
  domain = c(1940, 2021.43)
)

for (i in 1:nrow(diffusion.data)) {
  location_1 <- c(diffusion.data[i,]$start_1, diffusion.data[i,]$end_1)
  location_2 <- c(diffusion.data[i,]$start_2, diffusion.data[i,]$end_2)
  g <- g+geom_line(data= data.frame(location_1,location_2),aes(x= location_1,y= location_2))
  #if(diffusion.data[i,]$isTip == TRUE){
  if(diffusion.data[i,]$isTip == TRUE){
    g <- g+geom_point(data = diffusion.data[i,],aes(x=end_1,y=end_2),shape=16,stroke=0,size=3,color=col_fun(diffusion.data[i,]$time),alpha=0.8)
  }
}

(g <- g + 
    xlab("Antigenic dimension 1") +
    ylab("Antigenic dimension 2")+
    theme_classic())

ggsave("../figures/NA_antigenic_1_2_color_by_year.pdf",g,width=20,height=9,units="cm")


#antigenic dimension 1 and 2 colored by clusters

g <- ggplot()
location.list <- c()
for (i in 1:nrow(diffusion.data)) {
  location_1 <- c(diffusion.data[i,]$start_1, diffusion.data[i,]$end_1)
  location_2 <- c(diffusion.data[i,]$start_2, diffusion.data[i,]$end_2)
  g <- g+geom_line(data= data.frame(location_1,location_2),aes(x= location_1,y= location_2))

  if(diffusion.data[i,]$isTip == TRUE){
    
    if(diffusion.data[i,]$end_1 < -1 | diffusion.data[i,]$strain == "B/Lee/1940|1940"){  #cluster 1
      g <- g+geom_point(data = diffusion.data[i,],aes(x=end_1,y=end_2),shape=16,stroke=0,size=3,color="#59595A",alpha=0.8)
    } else if(diffusion.data[i,]$end_1 < 1){
      g <- g+geom_point(data = diffusion.data[i,],aes(x=end_1,y=end_2),shape=16,stroke=0,size=3,color="#A02B93",alpha=0.8)
    } else if(diffusion.data[i,]$end_1 < 3.01){
      g <- g+geom_point(data = diffusion.data[i,],aes(x=end_1,y=end_2),shape=16,stroke=0,size=3,color="#E87132",alpha=0.8)
    } else if(diffusion.data[i,]$end_1 < 7){
      g <- g+geom_point(data = diffusion.data[i,],aes(x=end_1,y=end_2),shape=16,stroke=0,size=3,color="#0070C0",alpha=0.8)
    } else{
      g <- g+geom_point(data = diffusion.data[i,],aes(x=end_1,y=end_2),shape=16,stroke=0,size=3,color="#1A6B23",alpha=0.8)
    }
  }
}
(g <- g + 
    xlab("Antigenic dimension 1") +
    ylab("Antigenic dimension 2")+
    theme_classic())

ggsave("../figures/NA_antigenic_1_2_cluster.pdf",g,width=20,height=9,units="cm")
