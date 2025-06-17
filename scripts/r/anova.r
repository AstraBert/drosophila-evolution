library(dplyr)
library(purrr)

df <- read.csv("data/f4_stats_all/data_for_anova.csv")
# do anova by chromatin (het, nohet)
anova1 <- aov(F4 ~ Origin, data = df)
# view anova summary
summary(anova1)

# do post-hoc
TukeyHSD(anova1)

# do anova by origin (western, eastern...)
anova_by_origin <- df %>%
  group_by(Origin) %>%
  group_split() %>%
  map(~ aov(F4 ~ Group, data = .x)) %>%
  map(summary)

# view summary
anova_by_origin

# only map anova result (no summary)
anova_by_origin1 <- df %>%
  group_by(Origin) %>%
  group_split() %>%
  map(~ aov(F4 ~ Group, data = .x))

# do post-hoc
tukey_results <- map(anova_by_origin1, ~ {
  tryCatch(
    TukeyHSD(.x),
    error = function(e) NA
  )
})
# view post-hoc summary
tukey_results
