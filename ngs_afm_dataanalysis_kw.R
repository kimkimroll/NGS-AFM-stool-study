
setwd("//cdc.gov/project/CCID_NCIRD_DVD_PPLB/_PMDDL/Kim/R/terry")

library(ggplot2)

dat <-read.csv("201124_Polio83_AFM_FinalAnalysis_CC_V3_edit.csv", stringsAsFactors = FALSE, header = TRUE)

arm_col <-c("Arm1" = "#f46036",
           "Arm2" = "#3bacb2",
           "Arm2_nested" = "#2e294e")
tnames <-c("Arm1" = "Arm 1",
           "Arm2" = "Arm 2",
           "Arm2_nested" = "Arm 2 Nested")



#regrouping PV1
library(dplyr)
dat1 <-  dat %>%
  mutate(Type=replace(Type, Type=="PV1_SABIN", "PV1 Sabin")) %>%
  mutate(Type=replace(Type, Type=="PV1_NSL", "PV1 NSL")) %>%
  mutate(itd5_1=replace(itd5_1, itd5_1=="SL1", "PV1 Sabin")) %>%
  mutate(itd5_1=replace(itd5_1, itd5_1=="NSL1", "PV1 NSL")) %>%
  mutate(itd5_1=replace(itd5_1, itd5_1=="PV2 (sabin2ct)", "PV2")) %>%
  mutate(itd5_2=replace(itd5_2, itd5_2=="SL1", "PV1 Sabin")) %>%
  mutate(itd5_2=replace(itd5_2, itd5_2=="NSL1", "PV1 NSL")) %>%
  mutate(itd5_1=replace(itd5_1, itd5_1=="SL3", "PV3")) %>%
  mutate(itd5_2=replace(itd5_2, itd5_2=="SL3", "PV3"))
  #replace(PV1_NSL, PV1)
dat1$Type <- factor(dat1$Type, levels = c("Cox A13", "Noro GII.4", "PV1 Sabin", "PV1 NSL", "PV2", "PV3"))
dat1$itd5_1 <- factor(dat1$itd5_1, levels = c("PV1 Sabin", "PV1 NSL", "PV2", "PV3"))
dat1$itd5_2 <- factor(dat1$itd5_2, levels = c( "PV1 Sabin", "PV1 NSL", "PV2", "PV3"))
dat1$dash<-as.character(dat1$dash)

#--------------------------------------------------------editing up to here
 
#virus order
dat$Type <- factor(dat$Type, levels = c("Cox A13", "Noro GII.4",  "PV1_SABIN",  "PV1","PV1_NSL", "PV2", "PV3"))
dat$itd5_1<-as.character(dat$itd5_1)
dat$itd5_1 <- factor(dat$itd5_1, levels = c("SL1", "NSL1", "PV2 (sabin2ct)", "SL3"))
dat$itd5_2<-as.character(dat$itd5_2)
dat$itd5_2<- factor(dat$itd5_2, levels = c("SL1", "NSL1", "PV2 (sabin2ct)", "SL3"))
dat$dash<-as.character(dat$dash)

#arm order
dat$arm <- factor(dat$arm, levels = c("Arm1", "Arm2", "Arm2_nested"))
#sample id
dat$dash<-as.factor(dat$dash)

dat$coverage[dat$coverage == ""] <- NA
dat$ct_1[dat$ct_1 == ""] <- NA
dat$coverage<-as.numeric(dat$coverage)
dat$ct_1<-as.numeric(dat$ct_1)



#plot1

kt <-ggplot(dat, #aes(x = coverage, 
                   #y = sample_x, 
                     #group = arm, 
                     #color = arm
                   #),
            na.rm = TRUE) +
  geom_bar(aes(x = coverage, 
               y = dash, 
               fill = coverage, group = Type),
           
           size = 2,
           position="dodge", 
           stat="identity", 
           alpha = 0.6,
           na.rm = TRUE)+
  geom_text(aes(x = coverage, 
                y = dash,
                group = itd5_1,
                label = itd5_1), 
           # y = sample_x, 
            size = 3, position=position_dodge(width=0.8), vjust=.3, hjust = -.2)+
  geom_point(aes(x = ct_1*50, y = dash, 
                 group = itd5_1
                 ), shape = 10, size = 3)+
  geom_point(aes(x = ct_2*50, y = dash,
                 group = itd5_1), shape = 10, size = 3)+
 # scale_fill_manual(values = tcolor, labels = tnames, name = "Method")+
  #scale_color_manual(values = tcolor, labels = tnames, name = "Arm")+
  #scale_y_discrete(na.translate = FALSE)+
  scale_x_continuous(name="Coverage",  limits = c(0, 5000), sec.axis = sec_axis(~./50, name="Ct"))+
  theme_minimal()+
  theme(strip.text = element_text(size = 12))+
  #facet_grid(Type~.)+
  #lims(0, 5000)+
  labs(y = "Sample",
       x = "Coverage")

kt

ggsave("kim_tplot4.png", dpi = 500, units = "cm", height = 24, width = 34)


#------------------------------------------- CORRELATION

library(wesanderson)
col = wes_palette("Rushmore1", 19, type = "continuous")

t_c<-ggplot(dat)+
    geom_line(aes(x = Type, y = itd5_1, color = Type), 
             size = 1,
             #color = "black",
             #alpha = 0.4,
             na.rm = TRUE)+
 # geom_line(aes(x = Type, y = itd5_2), 
  #          size = 1,
  #          #alpha = 0.4,
  #          na.rm = TRUE)+
  geom_point(aes(x = Type, y = itd5_1, color = dash),
            size = 1,
            alpha = 0.4,
             na.rm = TRUE)+
  geom_point(aes(x = Type, y = itd5_2, color = dash), 
             size = 1,
            alpha = 0.4,
            na.rm = TRUE)+
  scale_y_discrete(limits = c("SL1", "NSL1", "PV2 (sabin2ct)", "SL3"),
                   #labels = c("SL1", "PV2 (sabin2ct)", "SL3", "NSL1")
                   )+
  scale_color_manual(values = col)+
  #scale_fill_viridis(option = "D", discrete = FALSE)+
  #theme_minimal()+
  theme(
    #legend.position = "none"
  )+
  labs(x = "Virus Type from Tree",
       y = "VIrus Type from ITD")

t_c

#-----------------------------------------------------this works

t_col <- c("2008713728" = "#a6cee3",
"2008716075" = "#1f78b4",
"2012708714" = "#b2df8a",
"2012709971" = "#33a02c",
"2012710075" = "#fb9a99",
"2012803357" = "#e31a1c",
"2013717836" = "#fdbf6f",
"2013725830" = "#ff7f00",
"2013726201" = "#cab2d6",
"2014706709" = "#6a3d9a",
"2014706710" = "#ffff99",
"2014707869" = "#b15928")



t_c<-ggplot(dat1)+
  geom_tile(aes(x = Type, y = itd5_1), fill = "#d6e2e9",
            na.rm = TRUE)+
  geom_tile(aes(x = Type, y = itd5_2), fill = "#d6e2e9",
            na.rm = TRUE)+
  geom_jitter(aes(x = Type, y = itd5_1, color = arm),
              shape = 1,
              stroke = 1,
             size = 4,
             width = 0.4, height = 0.4,
             na.rm = TRUE)+
  geom_jitter(aes(x = Type, y = itd5_2, color = arm), 
              shape = 1,
             size = 4,
             width = 0.4, height = 0.4,
                        stroke = 1,
             na.rm = TRUE)+
  geom_point(aes(x = Type, y = itd5_1, fill = dash),
             color = "darkgrey",
             shape = 3,
             size = 3,
              alpha = 1,
             na.rm = TRUE)+
  scale_y_discrete(limits = c("PV1 Sabin", "PV1 NSL", "PV2", "PV3")
  )+
  scale_color_manual(values = arm_col, name = "Method")+
  scale_fill_manual(name = "n = 12 samples", values = t_col)+
  theme_minimal()+
theme(
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 14)
  #legend.position = "none"
)+
  labs(x = "Virus Type from Tree",
       y = "Virus Type from ITD")
#  facet_grid(.~arm)
 

  
t_c
ggsave("kim_tplot_correlation2.png", dpi = 500, units = "cm", height = 18, width = 32)




#ct plotting----------------------------------------------------------------------


t_plot<-ggplot(dat1)+
 geom_bar(aes(x = coverage, y = dash, fill = arm),
          width = 1,
          stat = "identity",
          position = "dodge",
           #fill = "#013a63", 
          na.rm = TRUE)+
#  geom_point(aes(x = ct_1*125, y = dash),
#             size = 3)+
#  geom_point(aes(x = ct_2*125, y = dash), 
#             shape = 1,
#             size = 3)+
  scale_fill_manual(values = arm_col)+
#  scale_fill_gradientn(name = "Coverage", 
#                     #low = "white", high = "#2a6f97", 
#                      values = c(0, 500, 1000, 2500, 3500, 4500),
#                      colours = c("#eff3ff", "#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"))+
 # scale_y_continuous(name="Ct",  limits = c(35,15), trans = "reverse")+
  #scale_alpha_continuous(limits = c(0, 5000), breaks = c(0, 1000, 2500, 5000), labels = c(0, 1000, 2500, 5000))+
  #scale_alpha(range = c(5000, 0), limits = c(5000, 0))+
  scale_x_continuous(name="Coverage",  limits = c(0, 5000), 
                     #sec.axis = sec_axis(~./125, name="Ct")
                     )+
  theme_bw()+
  theme(
    axis.title = element_text(size = 20),
    legend.position = "bottom"
  )+
#  labs(x = "Virus Type from Tree",
#       y = "Virus Type from ITD")
  facet_grid(rows = vars(Type), 
             #space = "free", 
             scales = "free")



t_plot

ggsave("kim_tplot_ct1.png", dpi = 500, units = "cm", height = 12, width = 24)




#trying something----------------------------------------------------------------
arm_shape<-c("Arm1" = '49',
             "Arm2" = '50',
             "Arm2_nested" = '51')
arm_shape<-as.numeric(arm_shape)



t_plot<-ggplot()+
  geom_line(dat1,(aes(x = coverage, y = dash, color= ct_1)),
            #alpha = 0.2,
            size = 2,
           )+ 
  geom_line(dat1,(aes(x = coverage, y = dash, color= ct_2)),
            #alpha = 0.2,
            size = 2,
            )+
#  geom_point(aes(x = coverage, y = dash, fill = arm),
#             #shape = 1,
#             size = 3)+
#  scale_fill_manual(name = "Arm", values = arm_col)+

  scale_color_gradient(low = "white", high = "steelblue")

 # theme_bw()+
 # theme(
#    axis.title = element_text(size = 20),
#    legend.position = "bottom"
#  )#+
 
  #  labs(x = "Virus Type from Tree",
  #       y = "Virus Type from ITD")
  #facet_grid(rows = vars(Type), 
  #           #space = "free", 
  #           scales = "free")


 t_plot
 
 #----------------------------------------------------------
 #comparing ct values
 
 
 t_ct<-ggplot(dat1)+
   geom_point(aes(x = ct_1, y = dash, color = arm),
              size = 3)+
   geom_point(aes(x = ct_2, y = dash, color = arm),
              size = 3)+
   scale_color_manual(values = arm_col, name = "Arm")+
   scale_x_continuous(name = "Ct", 
                      limits = c(35,15), 
                      trans = "reverse"
                      )+
   theme_bw()+
 facet_grid(rows = vars(Type), 
            #space = "free", 
            scales = "free")
 
 
 
 t_ct
 
 ggsave("kim_tplot_ct2.png", dpi = 500, units = "cm", height = 12, width = 24)
 
 
 
#---------------------------------------------------------------------------------
 t_ct<-ggplot(dat1)+
#   geom_tile(aes(x = coverage, y = ct_1, fill = coverage),
#             )+
#   geom_tile(aes(x = coverage, y = ct_2, fill = coverage))+
#   geom_line(aes(x = coverage, y = ct_1, style = Type))+
   geom_point(aes(x = ct_1, y = coverage, color = dash),
              shape = 1,
              size = 3)+
   geom_point(aes(x = ct_2, y = coverage, color = dash),
              shape = 1,
              size = 3)+
   geom_bar(aes(x = ct_1, y = coverage, fill = dash),
            width = 1,
            stat = "identity",
            position = "dodge")+
   geom_bar(aes(x = ct_2, y = coverage, fill = dash),
            width = 1, 
            stat = "identity",
            position = "dodge")+
   scale_color_manual(values = t_col, name = "Sample")+
   scale_fill_manual(values = t_col, name = "Sample")+
   scale_x_continuous(name = "Ct", 
                      limits = c(35,15), 
                      trans = "reverse"
   )+
   theme_bw()+
   facet_grid(cols = vars(Type), 
              space = "free", 
              scales = "free")
 
 
 
 t_ct
 
 ggsave("kim_tplot_ct3samp.png", dpi = 500, units = "cm", height = 16, width = 24)
 
 
 t_ct<-ggplot(dat1)+
   #   geom_tile(aes(x = coverage, y = ct_1, fill = coverage),
   #             )+
   #   geom_tile(aes(x = coverage, y = ct_2, fill = coverage))+
   #   geom_line(aes(x = coverage, y = ct_1, style = Type))+
   geom_bar(aes(x = ct_1, y = coverage, fill = arm),
            #alpha = 0.8,
            width = 1,
            stat = "identity",
            #position = position_dodge(0.7)
            position = "identity"
            )+
   geom_bar(aes(x = ct_2, y = coverage, fill = arm),
            width = 1, 
            #alpha = 0.8,
            stat = "identity",
            #position = position_dodge(0.7)
            position = "identity"
            )+
   geom_point(aes(x = ct_1, y = coverage),
              color = "white",
              #color = 1,
              stroke = 1,
              size = 4.5)+
   geom_point(aes(x = ct_2, y = coverage),
              color = "white",
              #shape = 1,
              stroke = 1,
              size = 4.5)+
   geom_point(aes(x = ct_1, y = coverage, color = dash),
              #fill = "white",
              #shape = 1,
              #stroke = 2,
              size = 4)+
   geom_point(aes(x = ct_2, y = coverage, color = dash),
              #fill = "white",
              #shape = 1,
              #stroke = 2,
              size = 4)+
   scale_color_manual(values = t_col, name = "n = 12 samples")+
   scale_fill_manual(values = arm_col, name = "Arm")+
   scale_x_continuous(name = "Ct", 
                      limits = c(35,15), 
                      trans = "reverse"
   )+
   #scale_y_continuous(breaks = c(0,1000,2500,5000,))
   theme_bw()+
   facet_grid(cols = vars(Type), 
              space = "free", 
              scales = "free")
 
 
 
 t_ct
 
 ggsave("kim_tplot_ct3arm1.png", dpi = 500, units = "cm", height = 16, width = 36)
 
 #-------------------------------------------------------------------------------------------ahh
 

 library(dplyr)
 dat2 <-  dat %>%
   mutate(Type=replace(Type, Type=="PV1_SABIN", "PV1 Sabin")) %>%
   mutate(Type=replace(Type, Type=="PV1_NSL", "PV1 NSL")) %>%
   mutate(itd5_1=replace(itd5_1, itd5_1=="SL1", "PV1 Sabin")) %>%
   mutate(itd5_1=replace(itd5_1, itd5_1=="NSL1", "PV1 NSL")) %>%
   mutate(itd5_1=replace(itd5_1, itd5_1=="PV2 (sabin2ct)", "PV2")) %>%
   mutate(itd5_2=replace(itd5_2, itd5_2=="SL1", "PV1")) %>%
   mutate(itd5_2=replace(itd5_2, itd5_2=="NSL1", "PV1 NSL")) %>%
   mutate(itd5_1=replace(itd5_1, itd5_1=="SL3", "PV3")) %>%
   mutate(itd5_2=replace(itd5_2, itd5_2=="SL3", "PV3"))
 #replace(PV1_NSL, PV1)
 dat2$Type <- factor(dat2$Type, levels = c("Cox A13", "Noro GII.4", "PV1 Sabin", "PV1 NSL", "PV2", "PV3"))
 dat2$itd5_1 <- factor(dat2$itd5_1, levels = c("PV1 Sabin", "PV1 NSL", "PV2", "PV3"))
 dat2$itd5_2 <- factor(dat2$itd5_2, levels = c( "PV1 Sabin", "PV1 NSL", "PV2", "PV3"))
 
 library(dplyr)
 dat2 <- dat2 %>% filter(mixture == "y") 
 dat2$dash<-as.character(dat2$dash)
 
 
 library(stringr)
 dat2$dash1 <- str_sub(dat2$dash,-3)
 
 library(ggrepel)
 
 t_ct<-ggplot(dat2)+
   #   geom_tile(aes(x = coverage, y = ct_1, fill = coverage),
   #             )+
   #   geom_tile(aes(x = coverage, y = ct_2, fill = coverage))+
   #   geom_line(aes(x = coverage, y = ct_1, style = Type))+

   geom_point(aes(x = Type, y = ct_1, size = coverage, color = arm, group = arm),
              alpha = 0.5,
              stroke = 2,
              position = position_dodge(width=.8),
              #shape = 3,
              #color = "darkgrey",
              #size = 1
              na.rm = TRUE)+
   geom_point(aes(x = Type, y = ct_2, size = coverage, color = arm, group = arm),
              alpha = 0.5,
              stroke = 2,
              position = position_dodge(width=.8),
              #shape = 3,
              #color = "darkgrey",
              #size = 1
              na.rm = TRUE)+
   geom_point(aes(x = Type, y = ct_1, fill = dash),
              shape = 3,
              color = "black",
              alpha = 0.25,
              size = 2,
              na.rm = TRUE
   )+
   geom_point(aes(x = Type, y = ct_2, fill = dash),
              shape = 3,
              color = "black",
              alpha = 0.25,
              size = 2,
              na.rm = TRUE
   )+
   geom_text(aes(x = Type, y = ct_1, label = dash1),
              #shape = 3,
              color = "black",
              alpha = 0.5,
              size = 5,
             vjust = -0.5,
              na.rm = TRUE
   )+
   geom_text(aes(x = Type, y = ct_2, label = dash1),
             #shape = 3,
             color = "black",
             alpha = 0.5,
             size = 5,
             vjust = -0.5,
             na.rm = TRUE
   )+

   scale_y_continuous(name = "Ct", 
                      limits = c(35,15), 
                      trans = "reverse"
   )+ 
   scale_color_manual(values = arm_col, name = "Method")+
   scale_size(name = "Coverage", range = c(5,35), breaks = c(0,100,1000,2500,4000))+
   #scale_shape(name = "n=12 samples")+
   scale_x_discrete()+
   scale_fill_discrete(name = "n = 3 samples")+
   #scale_y_continuous(breaks = c(0,1000,2500,5000,))
   #guide(shape = "n=12")+
   theme_bw()+
   theme(
     #text = element_text(family = "URWGothic"),
     axis.title = element_text(size = 18),
     axis.text = element_text(size = 24),
     legend.text = element_text(size = 16),
     legend.title = element_text(size = 18),
     #legend.position = "bottom"
   )
 
 t_ct
 
 ggsave("kim_tplot_ct3arm2_edit_mixtures1.png", dpi = 500, units = "cm", height = 24, width = 42)
 

 
 dat3 <- dat1 %>% filter(mixture == "n") 
 dat3$dash1 <- str_sub(dat3$dash,-3)
 
 t_ct<-ggplot(dat3)+
   #   geom_tile(aes(x = coverage, y = ct_1, fill = coverage),
   #             )+
   #   geom_tile(aes(x = coverage, y = ct_2, fill = coverage))+
   #   geom_line(aes(x = coverage, y = ct_1, style = Type))+
   
   geom_point(aes(x = Type, y = ct_1, size = coverage, color = arm, group = arm),
              alpha = 0.5,
              stroke = 2,
              position = position_dodge(width=.8),
              #shape = 3,
              #color = "darkgrey",
              #size = 1
              na.rm = TRUE)+
   geom_point(aes(x = Type, y = ct_2, size = coverage, color = arm, group = arm),
              alpha = 0.5,
              stroke = 2,
              position = position_dodge(width=.8),
              #shape = 3,
              #color = "darkgrey",
              #size = 1
              na.rm = TRUE)+
   geom_point(aes(x = Type, y = ct_1, fill = dash),
              shape = 3,
              color = "black",
              alpha = 0.25,
              size = 2,
              na.rm = TRUE
   )+
   geom_point(aes(x = Type, y = ct_2, fill = dash),
              shape = 3,
              color = "black",
              alpha = 0.25,
              size = 2,
              na.rm = TRUE
   )+
   geom_text(aes(x = Type, y = ct_1, label = dash1),
             #shape = 3,
             color = "black",
             alpha = 0.5,
             size = 5,
             vjust = -0.5,
             na.rm = TRUE
   )+
   geom_text(aes(x = Type, y = ct_2, label = dash1),
             #shape = 3,
             color = "black",
             alpha = 0.5,
             size = 5,
             vjust = -0.5,
             na.rm = TRUE
   )+
   
   scale_y_continuous(name = "Ct", 
                      limits = c(35,15), 
                      trans = "reverse"
   )+ 
   scale_color_manual(values = arm_col, name = "Method")+
   scale_size(name = "Coverage", range = c(5,35), breaks = c(0,100,1000,2500,4000))+
   #scale_shape(name = "n=12 samples")+
   scale_x_discrete()+
   scale_fill_discrete(name = "n = 9 samples")+
   #scale_y_continuous(breaks = c(0,1000,2500,5000,))
   #guide(shape = "n=12")+
   theme_bw()+
   theme(
     #text = element_text(family = "URWGothic"),
     axis.title = element_text(size = 18),
     axis.text = element_text(size = 24),
     legend.text = element_text(size = 16),
     legend.title = element_text(size = 18),
     #legend.position = "bottom"
   )
 
 t_ct
 
 ggsave("kim_tplot_ct3arm2_edit_mixtures2.png", dpi = 500, units = "cm", height = 24, width = 42)
 

 
 
 library(dplyr)
 dat2 <-  dat %>%
   mutate(Type=replace(Type, Type=="PV1_SABIN", "PV1 Sabin")) %>%
   mutate(Type=replace(Type, Type=="PV1_NSL", "PV1 NSL")) %>%
   mutate(itd5_1=replace(itd5_1, itd5_1=="SL1", "PV1 Sabin")) %>%
   mutate(itd5_1=replace(itd5_1, itd5_1=="NSL1", "PV1 NSL")) %>%
   mutate(itd5_1=replace(itd5_1, itd5_1=="PV2 (sabin2ct)", "PV2")) %>%
   mutate(itd5_2=replace(itd5_2, itd5_2=="SL1", "PV1")) %>%
   mutate(itd5_2=replace(itd5_2, itd5_2=="NSL1", "PV1 NSL")) %>%
   mutate(itd5_1=replace(itd5_1, itd5_1=="SL3", "PV3")) %>%
   mutate(itd5_2=replace(itd5_2, itd5_2=="SL3", "PV3"))
 #replace(PV1_NSL, PV1)
 dat2$Type <- factor(dat2$Type, levels = c("Cox A13", "Noro GII.4", "PV1 Sabin", "PV1 NSL", "PV2", "PV3"))
 dat2$itd5_1 <- factor(dat2$itd5_1, levels = c("PV1 Sabin", "PV1 NSL", "PV2", "PV3"))
 dat2$itd5_2 <- factor(dat2$itd5_2, levels = c( "PV1 Sabin", "PV1 NSL", "PV2", "PV3"))
 
 library(dplyr)
 #dat2 <- dat2 %>% filter(mixture == "y") 
 dat2$dash<-as.character(dat2$dash)
 
 
 library(stringr)
 dat2$dash1 <- str_sub(dat2$dash,-3)
 
 t_ct<-ggplot(dat2)+
   #   geom_tile(aes(x = coverage, y = ct_1, fill = coverage),
   #             )+
   #   geom_tile(aes(x = coverage, y = ct_2, fill = coverage))+
   #   geom_line(aes(x = coverage, y = ct_1, style = Type))+
   
   geom_point(aes(x = Type, y = ct_1, size = coverage, color = arm, group = arm),
              alpha = 0.5,
              stroke = 2,
              position = position_dodge(width=.8),
              #shape = 3,
              #color = "darkgrey",
              #size = 1
              na.rm = TRUE)+
   geom_point(aes(x = Type, y = ct_2, size = coverage, color = arm, group = arm),
              alpha = 0.5,
              stroke = 2,
              position = position_dodge(width=.8),
              #shape = 3,
              #color = "darkgrey",
              #size = 1
              na.rm = TRUE)+
   geom_point(aes(x = Type, y = ct_1, fill = dash),
              shape = 3,
              color = "black",
              alpha = 0.25,
              size = 2,
              na.rm = TRUE
   )+
   geom_point(aes(x = Type, y = ct_2, fill = dash),
              shape = 3,
              color = "black",
              alpha = 0.25,
              size = 2,
              na.rm = TRUE
   )+
   geom_text(aes(x = Type, y = ct_1, label = dash1),
             #shape = 3,
             color = "black",
             alpha = 0.5,
             size = 5,
             vjust = -0.5,
             na.rm = TRUE
   )+
   geom_text(aes(x = Type, y = ct_2, label = dash1),
             #shape = 3,
             color = "black",
             alpha = 0.5,
             size = 5,
             vjust = -0.5,
             na.rm = TRUE
   )+
   
   scale_y_continuous(name = "Ct", 
                      limits = c(35,15), 
                      trans = "reverse"
   )+ 
   scale_color_manual(values = arm_col, name = "Method")+
   scale_size(name = "Coverage", range = c(5,35), breaks = c(0,100,1000,2500,4000))+
   #scale_shape(name = "n=12 samples")+
   scale_x_discrete()+
   scale_fill_discrete(name = "n = 12 samples")+
   #scale_y_continuous(breaks = c(0,1000,2500,5000,))
   #guide(shape = "n=12")+
   theme_bw()+
   theme(
     #text = element_text(family = "URWGothic"),
     axis.title = element_text(size = 18),
     axis.text = element_text(size = 24),
     legend.text = element_text(size = 16),
     legend.title = element_text(size = 18),
     #legend.position = "bottom"
   )
 
 t_ct
 
 ggsave("kim_tplot_ct3arm2_edit_mall.png", dpi = 500, units = "cm", height = 26, width =50)
 