library(ggplot2)
library(dplyr)
library(readxl)
library(readr)

ha.map <- read.csv("../table/HI_only_tips_coordinates.csv") %>% dplyr::rename(ha_1 = x_coordinate, ha_2 = y_coordinate)
ha.cluster <- read_excel("../../data/HA/MK_IBV_clusterIDs_v2.xlsx") %>% 
  dplyr::rename(strain = Isolate) %>% 
  mutate(Cluster = gsub(" ","",Cluster)) %>% 
  dplyr::select(-c("Year"))

ha.map <- merge(x=ha.map, y=ha.cluster, by="strain",all.x=TRUE)

na.map <- read.csv("../table/NA_tips_coordinates.csv") %>% dplyr::rename(na_1 = x_coordinate, na_2 = y_coordinate)

ha.map$strain <- gsub('_A_E','',ha.map$strain)
ha.map$strain <- gsub('_V_E','',ha.map$strain)
ha.map$strain <- gsub('_V_C','',ha.map$strain)
ha.map$strain <- gsub('_Y_E','',ha.map$strain)
ha.map$strain <- gsub('_Y_C','',ha.map$strain)
ha.map$strain <- gsub('B/Ann_Arbor/1/86\\|1986','B/Ann_Arbor/1/1986|1986',ha.map$strain)
ha.map$strain <- gsub('B/BANGKOK/62/2002\\|2002-02-18','B/Bangkok/62/2002|2002-02-18',ha.map$strain)
ha.map$strain <- gsub('B/Bonn/43\\|1943','B/Bonn/1943|1943',ha.map$strain)
ha.map$strain <- gsub('B/BRISBANE/60/2008\\|2008-04-08','B/Brisbane/60/2008|2008-04-08',ha.map$strain)
ha.map$strain <- gsub('B/Colorado/06/2017\\|2017-02-25','B/Colorado/6/2017|2017-02-25',ha.map$strain)
ha.map$strain <- gsub('B/Crawley/46\\|1946','B/Crawley/1946|1946',ha.map$strain)
ha.map$strain <- gsub('B/Czechoslovakia/1/49\\|1949','B/Czechoslovakia/1/1949|1949',ha.map$strain)
ha.map$strain <- gsub('B/FLORIDA/4/2006\\|2006','B/Florida/4/2006|2006',ha.map$strain)
ha.map$strain <- gsub('B/Hong_Kong/05/1972\\|1972','B/Hong_Kong/5/1972|1972',ha.map$strain)
ha.map$strain <- gsub('B/HUBEI_WUJIAGANG/158/2009\\|2009','B/Hubei-Wujiagang/158/2009|2009',ha.map$strain)
ha.map$strain <- gsub('B/Johannesburg/33/1958\\|1958','B/Jonannesburg/33/1958|1958',ha.map$strain)
ha.map$strain <- gsub('B/MALAYSIA/2506/2004\\|2004','B/Malaysia/2506/2004|2004',ha.map$strain)
ha.map$strain <- gsub('B/Sydney/3/2009\\|2009-05-05','B/SYDNEY/3/2009|2009-05-05',ha.map$strain)
ha.map$strain <- gsub('B/SHANGDONG/7/1997\\|1997','B/Shangdong/7/1997|1997',ha.map$strain)
ha.map$strain <- gsub('B/SICHUAN/379/1999\\|1999','B/Sichuan/379/1999|1999',ha.map$strain)
ha.map$strain <- gsub('B/Singapore/222/79\\|1979','B/Singapore/222/1979|1979',ha.map$strain)
ha.map$strain <- gsub('B/Singapore/KK0148/2018\\|2018-01-29','B/Singapore/KK0148/2018.6|2018-01-29',ha.map$strain)
ha.map$strain <- gsub('B/SYDNEY/508/2010\\|2010-10-11','B/Sydney/508/2010|2010-10-11',ha.map$strain)


#ha.na.map[duplicated(ha.na.map$strain) | duplicated(ha.na.map$strain, fromLast = TRUE), ]
ha.na.map<- merge(x=na.map, y=ha.map, by="strain",all.x=TRUE) %>% 
  group_by(strain) %>% 
  slice_sample(n = 1) %>%   # Randomly keep one row per group
  ungroup() %>% 
  mutate(na_cluster_info = case_when(strain == "B/Lee/1940|1940" ~ "A",
                                     na_1 < -1 ~ "A",
                                     na_1 < 1 & na_1 >= -1 ~ "B",
                                     na_1 < 3.01 & na_1 >= 1 ~ "C",
                                     na_1 < 7  & na_1 >= 3.01 ~ "D",
                                     na_1 >= 7 ~ "E")) %>% 
  dplyr::rename(ha_cluster_info = Cluster)

write.csv(ha.na.map, file = "../table/HA_only_HI_NA_tips_coordinates.csv", row.names = FALSE,quote = F)

# map colored by lineage
lineage <- read_tsv("../Nextstrain/NA/data/IBV_NA_meta.tsv")

mycolors <- c("#636363", "#E41A1C","#377EB8")

ha.na.map.lineage <- merge(x=ha.na.map, y=lineage, by="strain",all.x=TRUE)

ha.na.map.lineage$lineage <- factor(ha.na.map.lineage$lineage,levels=c("A", "Y","V"))


# Read alpha point information
antigenic_points_alpha <- read.csv("../table/NA_antigenic_points_alpha.csv")

# Merge alpha values by strain

ha.na.map.lineage <- merge(
  ha.na.map.lineage,
  antigenic_points_alpha[, c("strain", "alpha")],
  by = "strain",
  all.x = TRUE
)

(g<- ggplot(ha.na.map.lineage, aes(x=ha_1, y=na_1, fill = lineage))+
  geom_point(aes(alpha = alpha), shape=21,stroke=0,size=3)+
  scale_fill_manual(values=mycolors)+
  theme_classic()+
    theme(legend.position=c(0.15,0.6)))

ggsave("../figures/HA_NA_antigenic_color_lineage_alpha.pdf",g,width=20,height=9,units="cm")


# map colored by year
ha.na.map.lineage$year <- sub("-.*", "", ha.na.map.lineage$date)
ha.na.map.lineage$year <- as.numeric(as.character(ha.na.map.lineage$year))


(g <- ggplot(ha.na.map.lineage, 
            aes(x = ha_1, y = na_1, fill = year)) +
  geom_point(shape = 21, stroke = 0, size = 3) +
  scale_fill_gradient(low = "#4daf4a", high = "#984ea3") +
  theme_classic() +
  theme(legend.position = c(0.15, 0.6)))

ggsave("../figures/HA_NA_antigenic_color_by_year.pdf",g,width=20,height=9,units="cm")


# map colored by NA cluster

mycolors <- c("grey8", "grey","deepskyblue","orangered1","royalblue1")

ha.na.map.lineage$na_cluster_info <- factor(ha.na.map.lineage$na_cluster_info,levels=c("A","B","C","D","E"))

(g<- ggplot(ha.na.map.lineage, aes(x=ha_1, y=na_1, fill = na_cluster_info, shape = lineage))+
    geom_point(stroke=0,size=3)+
    scale_fill_manual(values=mycolors)+
    # scale_shape_manual(values=c(21,22,23))+
    theme_classic()+
    theme(legend.position=c(0.15,0.6)))

ggsave("../figures/HA_NA_antigenic_color_na_cluster.pdf",g,width=20,height=9,units="cm")

# map colored by HA cluster

mycolors <- c("Anc40s" = "grey8",
            "Anc50-70" = "grey",
            "HK72" = "royalblue1",
            "Mal04" = "mediumblue",
            "Bris08-V1A" = "deepskyblue",
            "7" = "dodgerblue4",
            "Col17-V1A.1" = "slateblue",
            "Wash19-V1A.3" = "darkviolet",
            "Austria21-V1A.3a.2" = "purple4",
            "Yam90s" = "firebrick4",
            "Yam00s" = "orangered1",
            "Mass12-Y2" = "goldenrod1",
            "Phu13-Y3" = "red3")

ha.na.map.lineage$ha_cluster_info <- factor(ha.na.map.lineage$ha_cluster_info, levels = c("Anc40s","Anc50-70","HK72","Yam90s","Yam00s","Mal04","7",
                                                    "Bris08-V1A","Phu13-Y3","Mass12-Y2","Col17-V1A.1","Wash19-V1A.3",
                                                    "Austria21-V1A.3a.2"))

(g<- ggplot(ha.na.map.lineage, aes(x=ha_1, y=na_1, fill = ha_cluster_info))+
    geom_point(shape=21,stroke=0,size=3)+
    scale_fill_manual(values=mycolors)+
    theme_classic()+
    theme(legend.position=c(0.15,0.6)))

ggsave("../figures/HA_NA_antigenic_color_ha_cluster.pdf",g,width=20,height=9,units="cm")



