## Sensitivity analysis##
#########################
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyverse)
######################


## day 60

df <- '
Parameter	Best_h	Best_l	Worst_l	Worst_h
mu	-8.002	-4.965	9.131	20.836
cf	-4.136	-1.907	2.586	4.866
tl	-8.41	-5.606	10.824	25.358
pH	-4.74	-2.316	3.054	5.652
' %>% read_table2()


order.parameters <- df %>% arrange(Worst_h) %>%
  mutate(Parameter=factor(x=Parameter, levels=Parameter)) %>%
  select(Parameter) %>% unlist() %>% levels()

colnames(df) = c("Parameter", "Best, high","Best, low", "Worst, low","Worst, high")
width = 0.95

df.2 <- df %>% 
  # gather columns Lower_Bound and Upper_Bound into a single column using gather
  gather(key='type', value='output.value', "Best, high":'Worst, high') %>%
  # just reordering columns
  select(Parameter, type, output.value) %>%
  # create the columns for geom_rect
  mutate(Parameter=factor(Parameter, levels=order.parameters),
         ymin=pmin(output.value, 0),
         ymax=pmax(output.value, 0),
         xmin=as.numeric(Parameter)-width/2,
         xmax=as.numeric(Parameter)+width/2)

df.2$ymax[1:4] = df.2$ymin[5:8]
df.2$ymin[13:16] = df.2$ymax[9:12]

##png(width = 960, height = 540)
ggplot() + 
  geom_rect(data = df.2, 
            aes(ymax=ymax, ymin=ymin, xmin=xmin, xmax=xmax, fill=type), alpha=0.7) +
  scale_fill_manual(values=c("#000000", "#e6e6e6", "#666666","#b3b3b3"))+
  theme_classic() + 
  theme(axis.title.y=element_blank(), legend.position = 'bottom',
        legend.title = element_blank()) + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        title=element_text(size=24),
        legend.text = element_text(size=16))+
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = c(1:length(order.parameters)), 
                     labels = order.parameters) +
  ylab("Change in proportions of late blown cheese at day 60 (%)")+
  coord_flip()
##dev.off()





## day120
df <- '
Parameter	Best_h	Best_l	Worst_l	Worst_h
mu	-50.849	-26.219	16.37	24.449
cf	-10.335	-4.614	3.113	5.888
tl	-30.938	-15.007	11.154	18.796
ph	-26.084	-11.927	8.124	13.641
' %>% read_table2()


order.parameters <- df %>% arrange(Worst_h) %>%
  mutate(Parameter=factor(x=Parameter, levels=Parameter)) %>%
  select(Parameter) %>% unlist() %>% levels()

colnames(df) = c("Parameter", "Best, high","Best, low", "Worst, low","Worst, high")
width = 0.95

df.2 <- df %>% 
  # gather columns Lower_Bound and Upper_Bound into a single column using gather
  gather(key='type', value='output.value', "Best, high":'Worst, high') %>%
  # just reordering columns
  select(Parameter, type, output.value) %>%
  # create the columns for geom_rect
  mutate(Parameter=factor(Parameter, levels=order.parameters),
         ymin=pmin(output.value, 0),
         ymax=pmax(output.value, 0),
         xmin=as.numeric(Parameter)-width/2,
         xmax=as.numeric(Parameter)+width/2)

df.2$ymax[1:4] = df.2$ymin[5:8]
df.2$ymin[13:16] = df.2$ymax[9:12]

##png(width = 960, height = 540)
ggplot() + 
  geom_rect(data = df.2, 
            aes(ymax=ymax, ymin=ymin, xmin=xmin, xmax=xmax, fill=type), alpha=0.7) +
  scale_fill_manual(values=c("#000000", "#e6e6e6", "#666666","#b3b3b3"))+
  theme_classic() + 
  theme(axis.title.y=element_blank(), legend.position = 'bottom',
        legend.title = element_blank()) + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        title=element_text(size=24),
        legend.text = element_text(size=16))+
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = c(1:length(order.parameters)), 
                     labels = order.parameters) +
  ylab("Change in proportions of late blown cheese at day 120 (%)")+
  coord_flip()
##dev.off()

