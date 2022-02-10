#text2SDG experiment

library(text2sdg)
library(tidyverse)
library(officer)
library(koRpus) # To calculate readability indices

txt <- readLines("effortcreeptxt.txt") # read txt file
# read_docx("")# read docx


hits <- detect_sdg(txt) # run sdg analysis

# Calculate sustainable development goals performance score

total <- sum(hits$hit)

# extract sustainable development goal 14 score (this is the relevant goal for us).

hitsdf <- as.data.frame(hits)

# Need to use tryCatch here, maybe not

df1 <- dim(hitsdf %>% select(sdg,hit) %>% filter(sdg=='SDG-01') %>% select(hit))[1]

if (df1 > 0) {
    df1 < - hitsdf %>% select(sdg,hit) %>% filter(sdg=='SDG-01') %>% select(hit)
}


tryCatch(
  {df2 <- hitsdf %>% select(sdg,hit) %>% filter(sdg=='SDG-02') %>% select(hit)},
  error={message('error')})
tryCatch(
  {df14 <- hitsdf %>% select(sdg,hit) %>% filter(sdg=='SDG-14') %>% select(hit)},
  error = {message('error')})


df14 <- hitsdf %>% select(sdg,hit) %>% filter(sdg=='SDG-14') %>% select(hit)

df14score <- sum(df14)







