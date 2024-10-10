# !/usr/bin/env Rscript

# 1. clean environment 
rm(list=ls())

# 2. open libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
#library(readxl)
library(plotly)
library(tidyr)
library("cowplot")

# generate prefix string including date/ Initials
current_date <- gsub("-", "", as.character(Sys.Date()))# Get the current date as character
pre <- paste("MM", current_date, sep = "")

# 3. set working directory
setwd("/home/melanie/D36E_assays/disease_data")
dir_out <- "/home/melanie/D36E_assays/disease_output/"

col <- c("#fffbf1", "#ffeabc", "#ffbe2d", "#E69F00", "#D55E00")

# 4. load file
df <- read_csv("MM01102024_SEMIQUANT_visual_disease_score.csv")

# B) get data(frame) into desired shape ----
# 5. rename columns and select only relevant columns
df1 <- data.frame("date_scoring" = df$`Date scoring`,
                 "Effector_line" = df$`Effector line`, 
                 "OD" = df$Dilution,
                 "nr_score_0" = df$disease_score_0,
                 "nr_score_1" = df$disease_score_1,
                 "nr_score_2" = df$disease_score_2,
                 "nr_score_3" = df$disease_score_3,
                 "nr_score_4" = df$disease_score_4,
                 "nr_infiltration_spots" = df$`n (infiltration spots)`,
                 "avg_score" = df$`average disease score`
                 )


# drop first row that only contains numbers of disease scores
df1 <- df1[!(rownames(df1) == '1'),]
rm(df)

# C) calculate the percentages of disease score incidences and store in df ----
# 6. create input vectors

# vector of row numbers  
rows <- as.vector(1:length(df1$Effector_line)) 

# Vector of colnames containing scores
n_scores <- grep("nr_score_", colnames(df1), value = TRUE) 

# transform to long format so that each observation has its own row
df_long <- tidyr::pivot_longer(df1, 
                               cols = starts_with("nr_score_"), 
                               names_to = "disease_score", 
                               values_to = "count")

# add a column with the % of infiltration spots showing the symptom
df_long$result_perc <- (df_long$count / df_long$nr_infiltration_spots) * 100

head(df_long)

# D) create a stagged bar graph from the df_sorted ----

# Replace dates with dpi for easier interpretation
df_plot <- df_long %>%
  mutate(date_scoring = recode(date_scoring, 
                               "25/09/2024" = "Odpi", 
                               "27/09/2024" = "2dpi", 
                               "1/10/2024" = "6dpi")
         )


# 4.2 Define the function to generate and save plots for each subset
# Function for individual plots
f_dpi_plot <- function(df, titles) {

g <- ggplot(df, 
           aes(x = Effector_line, 
               y = count,
               fill = disease_score
           )
           ) +
 geom_bar(position = "fill", stat = "identity") +
 scale_fill_manual(values = col) +

coord_cartesian(ylim = c(0.05, 1.0))  +

#theme_bw() +
theme(
  plot.background = element_blank(),
  panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), # remove background, frame
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size=14),
      axis.title = element_text(size=14),
      legend.title = element_text(size=14), 
      legend.text = element_text(size=12)
      )+
ggtitle(paste0("timepoint: ", titles)) +
xlab("Effector line") +                
ylab("% disease score") +
theme(axis.text.x = element_text(angle = 45, hjust = 1)
      ) 

return(g)
}

# Function for grid of plots 
f_grid <- function(plot_list) {
#plot_list <- strain_plots
# Extract the legend from the first plot using cowplot::get_legend
legend <- cowplot::get_legend(plot_list[[1]])

# Remove the legend from all plots (except the one used to extract the legend)
plots_without_legend <- lapply(plot_list, function(p) {
p + theme(legend.position = "none")
})

# Use ggarrange to combine the plots into a grid and add the extracted legend
combined_plot <- ggarrange(
plotlist = plots_without_legend,
ncol = length(plot_list),           
common.legend = FALSE,
legend = NULL
)

# Add the extracted legend as a separate grob using cowplot::plot_grid
final_plot <- plot_grid(combined_plot, legend, 
                      ncol = 2,        # Arrange plots and legend side by side
                      rel_widths = c(3, 0.5))  # Adjust the space for plots vs legend

return(final_plot)
}

# Function to safe plots by pre/suffix as grids
f_save_grid_prefix <- function(plot_list, prefix){
  
  plot_names <- names(plot_list)
  plots_for_grid <- plot_list[str_extract(titles_dpi_OD, "^[^_]+") == prefix]
  title <- names(plots_for_grid)
  
  grid <- f_grid(plots_for_grid)
  
  ggsave(filename = paste0(dir_out, pre, "_", prefix, "_grid", ".svg"),
         plot = grid, 
         width = 9, 
         height = 4, 
         units = "in"
  )
  
}

# Function to draw grid and safed for sub-collections of plots
f_save_grid_suffix <- function(plot_list, suffix){
  
  plot_names <- names(plot_list)
  plots_for_grid <- plot_list[str_extract(plot_names, "(?<=_)[0-9.]+$") == suffix]
  title <- names(plots_for_grid)
  
  grid <- f_grid(plots_for_grid)
  
  ggsave(filename = paste0(dir_out, pre, "_", suffix, "_grid", ".svg"),
         plot = grid, 
         width = 9, 
         height = 4, 
         units = "in"
  )
  
}

# Plot + safe
grouped_dpi_OD <- df_plot %>%
  group_by(date_scoring, OD) %>%
  group_split() 

titles_dpi_OD <- map_chr(
  grouped_dpi_OD, ~ paste(unique(.x$date_scoring), 
                        "OD", unique(.x$OD),
                        sep = "_"
                        )
)

plots_dpi_OD <- map2(.x = grouped_dpi_OD,
                 .y = titles_dpi_OD,
                 .f = ~ f_dpi_plot(
                   df = .x, 
                   title = .y
                   )
                 )

# Name plots based on title
names(plots_dpi_OD) <- titles_dpi_OD

# extract unique factors by which plots will be safed
unique_dpi <- unique(df_plot$date_scoring)
unique_OD <- unique(df_plot$OD)

lapply(unique_dpi,
       f_save_grid_prefix,
       plot_list = plots_dpi_OD
       )

lapply(unique_OD,
       f_save_grid_suffix,
       plot_list = plots_dpi_OD
)



# # Save the plots using purrr::map2 and a lambda function
# map2(.x = plots_dpi_OD, 
# .y = titles_dpi_OD, 
# .f = ~ggsave(filename = paste0(dir_out, pre, "_boxplot_", SUF, .y, ".svg"),
#       plot = .x, 
#       width = 3.5, 
#       height = 4, 
#       units = "in", 
#       dpi = 300)
# )

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

# add column with % of disease-symptom positive infiltration spots
perc_pos <- subset(simple_df_2, !grepl("no death", condition))
perc_pos <- arrange(perc_pos, Effector )

m$perc_pos <- perc_pos$Percentage

# safe results in a txt file
write.table(m, file = paste(pre, sep = "","_stats_cell_death_2_3_4_positive.txt"))

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
