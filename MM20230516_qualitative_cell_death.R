# A) General basics ----
# 1. clean environment 
rm(list=ls())

# 2. open libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(readxl)
library(plotly)

# 3. set working directory
setwd("C:/Users/Mende012/Documents/Bioinformatics/Cell_death/Analysis_disease_symptoms_visual_scoring")

# 4. load file
df <- read_xlsx("MM20230516_summary_visual_score_cell_death.xlsx")

# B) get data(frame) into desired shape ----
# 5. rename columns and select only relevant columns
df <- data.frame("date_infiltration" = df$`Date infiltartion`,
                 "Effector_line" = df$`Effector line`, 
                 "nr_score_0" = df$disease_score_0,
                 "nr_score_1" = df$disease_score_1,
                 "nr_score_2" = df$disease_score_2,
                 "nr_score_3" = df$disease_score_3,
                 "nr_score_4" = df$disease_score_4,
                 "nr_infiltration_spots" = df$`n (infiltration spots)`,
                 "avg_score" = df$`average disease score`)


# drop first row that only contains numbers of disease scores
df <- df[!(rownames(df) == '1'),]

# C) calculate the percentages of disease score incidences ----
# 6. create input vectors

# vector of row numbers  
rows <- as.vector(1:length(df$Effector_line)) 

# Vector of colnames containing scores
n_scores <- grep("nr_score_", colnames(df), value = TRUE) 


# 7. Apply the p_calc function to each row for each nr_score
percentages <- lapply(n_scores, function(score) {
  lapply(rows, function(row_nr) {
    nr_sum <- df$nr_infiltration_spots[row_nr]
    nr_score <- df[[score]][row_nr]
    perc <- (nr_score / nr_sum) * 100
    return(perc)
  })
})

# transform that into a usable df

#----------------------------------------
# reshape df
# https://r-charts.com/part-whole/stacked-bar-chart-ggplot2/
# x axis - effector line
# fill - group (disease intensity 0-4 %)
# 
#-------------------------

# Basic stagged plot
library(ggplot2)

ggplot(df, aes(x = x, y = y, fill = group)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#DADAEB", "#9E9AC8", "#6A51A3")) 
#---------------

# Change from wide to long format
reshape(df, direction = "long". varyin  )




perc <- lapply(unique(df$assay_id), function(x) as.data.frame(DunnettTest(area_under_curve ~ treatment_id, 
                                                                                data = df[df$assay_id==x,], 
                                                                                control = "D36E_Bacterial PAMPs", 
                                                                                conf.level = 0.95)$`D36E_Bacterial PAMPs`)) %>% bind_rows()  

# pool control samples that are present in every experiment

# get a list of all the unique effector lines
e_lines <- unique(df$`Effector line`)

# write function to sum up all the +/- infiltration spots (to account for repeating control samples)

sum_d <-  lapply(e_lines, function(x){ # sum up of cell-death POSITIVE
          a <- df[(df$`Effector line` == x),] # subset datset by effector line
          b <- sum(a$`Positive = visual signs of cell death`) # sum infiltration spots WITH cell death
          return(b)
}) %>% unlist()

sum_n <- lapply(e_lines, function(x){ # sum up of cell-death NEGATIVE
  a <- df[(df$`Effector line` == x),] # subset datset by effector line
  b <- sum(a$`Negative = no visual signs of cell death`) # sum infiltration spots NO cell death
  return(b)
}) %>% unlist()


# store in new dataframe
simple_df_d <- as.data.frame(cbind(unlist(e_lines), as.numeric(sum_d))) # df cell death
simple_df_n <- as.data.frame(cbind(unlist(e_lines), as.numeric(sum_n))) # df no cell death


# rename columns of simple_df in a sensible way
colnames(simple_df_d)[c(1, 2)] = c("Effector", "count_death")
colnames(simple_df_n)[c(1, 2)] = c("Effector", "count_death")

# add column with "condition" which states cell death positive or negative 
# this is necessary to get data into the right shape for plotting
simple_df_d["condition"] <- as.vector(rep("death", times = length(e_lines)))
simple_df_n["condition"] <- as.vector(rep("no death", times = length(e_lines)))

# join the data frame for plotting
simple_df <- full_join(simple_df_d, simple_df_n)


##################
#################
################# probably add statistics here???? https://statsandr.com/blog/fisher-s-exact-test-in-r-independence-test-for-a-small-sample/

# function to sum up 
nr_spots <- lapply(e_lines, function(x){
  line <- simple_df[simple_df$Effector == x,]
  s <- sum(as.numeric(line$count_death))
  return(s)
}) 

# add nr all infiltration spots to df
simple_df["nr_spots"] <- unlist(nr_spots) 

# calculate percentag
percentage <- (as.numeric(simple_df$count_death)/as.numeric(simple_df$nr_spots))*100
simple_df["Percentage"] <- percentage

#simple_df <- simple_df[!(simple_df$condition== "no death"),]

#############################

# draw percent stacked barchart (https://r-graph-gallery.com/48-grouped-barplot-with-ggplot2.html)
  
g1 <- ggplot(simple_df, aes(fill = condition, x = Effector, y = Percentage)) +
         geom_col(position = position_fill(reverse = TRUE)) +
        scale_fill_manual(values = c("#E69F00", "#faebcc")) +
  
  # define axis limits if needed
  coord_cartesian(ylim = c(0.1, 1))  +

        # define the theme of the boxplot
  theme_bw() +  # make the bg white
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), # remove background, frame
        axis.line = element_line(colour = "black")) +

        # label the axises 
  xlab("Effector") +                
  ylab("% cell death positive infiltration spots)") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))

g1


ggsave("MM20230111_cell_death.svg",width = 7, height = 4)
g1
dev.off()

ggplotly(g1)

##############################################################################
#########################  export for summary heat map #######################
##############################################################################

# list all unique strains
e_lines

# only keep % of cell death positive incidence
d2 <- simple_df[!(simple_df$condition == "no death"),]

# drop all the columns not needed
d2 <- d2[c("Effector","Percentage")]

# Rename D36E, DC3000 to appear in the same order outputs in other analysis
d2[d2 == "AAAD36E"] <- "D36E"
d2[d2 == "ZZZDC3000"] <- "DC3000"

# export 
write.csv(d2, "./2_Cell_death_positive_percentage.csv")
