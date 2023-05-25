# A) General basics ----
# 1. clean environment 
rm(list=ls())

# generate prefix string
current_date <- Sys.Date() # Get the current date
print(MM, current_date)


# 2. open libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(readxl)
library(plotly)
library(tidyr)

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

# C) calculate the percentages of disease score incidences and store in df ----
# 6. create input vectors

# vector of row numbers  
rows <- as.vector(1:length(df$Effector_line)) 

# Vector of colnames containing scores
n_scores <- grep("nr_score_", colnames(df), value = TRUE) 

## 7. Design a function that can calculate the percentage of infiltration spots
# that was scored with a specific disease severity. Run the function to obtain a 
# list of values. 

result_perc <- do.call(cbind, lapply(n_scores, function(score) {
  t(sapply(rows, function(row_nr) {
    nr_sum <- df$nr_infiltration_spots[row_nr]
    nr_score <- df[[score]][row_nr]
    perc <- (nr_score / nr_sum) * 100
    return(perc)
  }))
})) %>% t()

# 8. Convert the original df to a long format so that the % results can be added
# 8.1 Convert the dataframe from wide to long format
df_long <- tidyr::pivot_longer(df, 
                               cols = starts_with("nr_score_"), 
                               names_to = "disease_score", 
                               values_to = "count")


# Sort the dataframe by first "disease_score", "disease_score", "Effector_line" columns in ascending order
df_sorted <- arrange(df_long, disease_score, v, Effector_line)

# 8.2 add the % (output of 7.) to the df_sorted, make EXTRA sure they are in the 
#     correct order
df_sorted$res_perc <- result_perc


# D) create a stagged bar graph from the df_sorted ----

# reshape df
# https://r-charts.com/part-whole/stacked-bar-chart-ggplot2/
# x axis - effector line
# fill - group (disease intensity 0-4 %)
# 
#-------------------------

ggplot(df_sorted, aes(x = Effector_line, y = res_perc, fill = disease_score)) + 
  geom_bar(stat = "identity")

# --> output super odd, did something go wrong before?

# E) Plot simple, stagged bar graph for all effector lines
simple_df <- data.frame(df_sorted$Effector_line, df_sorted$disease_score, df_sorted$count)
n_samples <- ###!!!

ggp <- ggplot(simple_df, aes(x = df_sorted.Effector_line, y = df_sorted.count,
                            fill = df_sorted.disease_score)) +
       geom_bar(position = "fill", stat = "identity") +
       # asign colours
       scale_fill_manual(values = c("#fffbf1", "#ffeabc", "#ffbe2d", "#E69F00", "#D55E00")) +
       # define axis limits if needed
       coord_cartesian(ylim = c(0.1, 1))  +
       # define the theme of the boxplot
       theme_bw() +  # make the bg white
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), # remove background, frame
        axis.line = element_line(colour = "black"))+
       # label the axises 
       xlab("Effector") +                
       ylab("% disease score") +
       theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
       # label individual bars
       geom_text(aes(label = n_samples), vjust = -0.5)




ggp                                 # Draw ggplot2 plot scaled to 1.00


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
