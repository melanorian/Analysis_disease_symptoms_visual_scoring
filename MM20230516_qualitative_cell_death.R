# A) General basics ----
# 1. clean environment 
rm(list=ls())

# Generate suffix for output names
SUF <- "1"

# generate prefix string including date/ Initials
current_date <- gsub("-", "", as.character(Sys.Date()))# Get the current date as character
pre <- paste("MM", current_date, sep = "")


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
df <- read_xlsx("MM20230601_summary_visual_score.xlsx")



# B) get data(frame) into desired shape ----
# 5. rename columns and select only relevant columns
df1 <- data.frame("date_infiltration" = df$`Date infiltartion`,
                 "Effector_line" = df$`Effector line`, 
                 "nr_score_0" = df$disease_score_0,
                 "nr_score_1" = df$disease_score_1,
                 "nr_score_2" = df$disease_score_2,
                 "nr_score_3" = df$disease_score_3,
                 "nr_score_4" = df$disease_score_4,
                 "nr_infiltration_spots" = df$`n (infiltration spots)`,
                 "avg_score" = df$`average disease score`)


# drop first row that only contains numbers of disease scores
df1 <- df1[!(rownames(df1) == '1'),]

# C) calculate the percentages of disease score incidences and store in df ----
# 6. create input vectors

# vector of row numbers  
rows <- as.vector(1:length(df1$Effector_line)) 

# Vector of colnames containing scores
n_scores <- grep("nr_score_", colnames(df1), value = TRUE) 

## 7. Design a function that can calculate the percentage of infiltration spots
# that was scored with a specific disease severity. Run the function to obtain a 
# list of values. 

result_perc <- do.call(cbind, lapply(n_scores, function(score) {
  t(sapply(rows, function(row_nr) {
    nr_sum <- df1$nr_infiltration_spots[row_nr]
    nr_score <- df1[[score]][row_nr]
    perc <- (nr_score / nr_sum) * 100
    return(perc)
  }))
})) %>% t()

# 8.!! Convert the original df to a long format so that the % results can be added
# 8.1 Convert the dataframe from wide to long format
df_long <- tidyr::pivot_longer(df1, 
                               cols = starts_with("nr_score_"), 
                               names_to = "disease_score", 
                               values_to = "count")


#!!! Sort the dataframe by first "disease_score", "disease_score", "Effector_line" columns in ascending order
df_sorted <- arrange(df_long, disease_score, Effector_line)

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

# E) Plot simple, stagged bar graph for all effector lines
simple_df <- data.frame(df_sorted$Effector_line, df_sorted$disease_score, df_sorted$count)

# get a list of all the unique effector lines
e_lines <- unique(df1$Effector_line)

# write function to sum up all the +/- infiltration spots (to account for repeating control samples)
sum_d <-  lapply(e_lines, function(x){ # sum up of cell-death POSITIVE
  a <- df1[(df1$Effector_line == x),] # subset datset by effector line
  b <- sum(a$nr_infiltration_spots) # sum infiltration spots WITH cell death
  return(b)
}) %>% unlist()

simple_df_d <- as.data.frame(cbind(unlist(e_lines), as.numeric(sum_d))) # df cell death'
simple_df_d <- arrange(simple_df_d, V1)
labels <- simple_df_d$V2

# add the labels to the simple_df
merged_df <- merge(simple_df, simple_df_d)
merged_df <- arrange(merged_df, merged_df$df_sorted.Effector_line)

# draw stagged bar graph
ggp <- ggplot(simple_df, aes(x = df_sorted.Effector_line, y = df_sorted.count,
                            fill = df_sorted.disease_score)) +
       geom_bar(position = "fill", stat = "identity") +
       # asign colours
       scale_fill_manual(values = c("#fffbf1", "#ffeabc", "#ffbe2d", "#E69F00", "#D55E00")) +
       # define axis limits if needed
       coord_cartesian(ylim = c(0.05, 1.0))  +
       # define the theme of the boxplot
       theme_bw() +  # make the bg white
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), # remove background, frame
        axis.line = element_line(colour = "black"))+
       # label the axises 
       xlab("Effector") +                
       ylab("% disease score") +
       theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
    #!!!   # label individual bars
       #geom_text(label = labels ,size = 3, hjust = 1)





ggp                                 # Draw ggplot2 plot scaled to 1.00


ggsave(paste(pre, sep = "","_cell_death.svg"),width = 7, height = 4)
ggp
dev.off()

# for exploring the plot ineractively 
# ggplotly(ggp)

############################################################################# 
############################# Statistics ####################################
#############################################################################

# E) Prepare data for Stats (partly redundand due to recycling of script) ----

# df with relevant information. Go back to original df 
### calculate the proportions in % of cell-death positive/negative infiltration spots
df2 <- data.frame("date_infiltration" = df$`Date infiltartion`,
                  "Effector_line" = df$`Effector line`, 
                  "infiltration_spots" = df$`n (infiltration spots)`,
                  "disease_pos" = df$`Positive = visual signs of cell death`,
                  "disease_neg" = df$`Negative = no visual signs of cell death`)
                  

# drop first row that only contains numbers of disease scores
df2 <- df2[!(rownames(df2) == '1'),]                 

# write function to sum up all the +/- infiltration spots (to account for repeating control samples)

sum_d <-  lapply(e_lines, function(x){ # sum up of cell-death POSITIVE
  a <- df2[(df2$Effector_line == x),] # subset datset by effector line
  b <- sum(a$disease_pos) # sum infiltration spots WITH cell death
  return(b)
}) %>% unlist()

sum_n <- lapply(e_lines, function(x){ # sum up of cell-death NEGATIVE
  a <- df2[(df2$Effector_line == x),] # subset datset by effector line
  b <- sum(a$disease_neg) # sum infiltration spots NO cell death
  return(b)
}) %>% unlist()


# store in new dataframe
simple_df_d <- as.data.frame(cbind(unlist(e_lines), as.numeric(sum_d))) # df cell death
simple_df_n <- as.data.frame(cbind(unlist(e_lines), as.numeric(sum_n))) # df no cell death


# rename columns of simple_df in a sensible way
colnames(simple_df_d)[c(1, 2)] = c("Effector", "count_death")
colnames(simple_df_n)[c(1, 2)] = c("Effector", "count_death")

# add column with "condition" which symptom positive or negative 
# this is necessary to get data into the right shape for plotting
simple_df_d["condition"] <- as.vector(rep("death", times = length(e_lines)))
simple_df_n["condition"] <- as.vector(rep("no death", times = length(e_lines)))

# join the data frame for plotting
simple_df_2 <- full_join(simple_df_d, simple_df_n)


# function to sum up 
nr_spots <- lapply(e_lines, function(x){
  line <- simple_df_2[simple_df_2$Effector == x,]
  s <- sum(as.numeric(line$count_death))
  return(s)
}) 

# add nr all infiltration spots to df
simple_df_2["nr_spots"] <- unlist(nr_spots) 

# calculate percentag
percentage <- (as.numeric(simple_df_2$count_death)/as.numeric(simple_df_2$nr_spots))*100
simple_df_2["Percentage"] <- percentage


# F) # Statistical analysis of results. Chi-square test cannot be performed due to too few
# samples (<5) in some groups --> Fisher's exact test ----

# bring datafram in a shape that's suitable for Fisher's exact test 

# re-name df so that info on cell-death-pos/neg can be merged based on effectors


simple_df_d$condition <- NULL
colnames(simple_df_d) <- c("Effector", "cell-death-positive")


simple_df_n$condition <- NULL
colnames(simple_df_n) <- c("Effector", "cell-death-negative")


# construct df for stat analysis
df_stat <- merge(simple_df_d, simple_df_n)

# Effectors as colname
rownames(df_stat) <- df_stat$Effector
df_stat$Effector <- NULL

# factors as numeric for applying statistics (oddly complicated but well...)

i <- c(1, 2)
df_stat[ , i] <- apply(df_stat[ , i], 2,          
                       function(x) as.numeric(as.character(x))) # first convert to character, important to retain the values


effect <- unique(rownames(df_stat))

# create a frequency table
#count(df_stat, df_stat$`cell-death-positive`)

#typeof(effect)
#df_stat <- as.factor(colnames(df_stat))

# expected frequencies per effector
#chisq.test(df_stat)$expected # --> results in error due to low expected frequencies

F_test <- sapply(effect, function(x){
  d <- df_stat[x,]                                #subset df to extract single eff row
  c <- df_stat[(rownames(df_stat) == "AAAD36E"),] #subset df to extract control row
  m <- rbind(c ,d)                                #join both rows to have a small df
  f <- fisher.test(m)                             #apply test
  p <- f$p.value
  return(p)}) %>% cbind() %>% data.frame()

colnames(F_test) <- "p-value"

# creat new column containing values of row names so that we have an identical
# column in both data sets that we can merge on 
F_test$Effectors <- row.names(F_test)
df_stat$Effectors <- row.names(df_stat)

# merge datafram
m <- merge(df_stat, F_test)

# safe results in a txt file
write.table(m, file = "./MM20230125_celldeath_assay_stats.txt")


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
