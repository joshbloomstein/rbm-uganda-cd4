library(tibble)
library(data.table)

# Load pathogen 2 x 2 table
# Pathogens present in fewer than 5 patients excluded due to lack of statistical power

df1 <- read.csv('07-29-21 All pathogens mortality at least 5 patients.csv', row.names = 1)

# Perform Fisher's exact test

values <- 1:dim(df1)[1]
values <- values[c(TRUE, FALSE)]

odd <- vector()
fisher.p <- vector()

for (i in values)    
{
  odd <- c(odd, fisher.test(df1[i:(i+1), 1:2])$estimate)
  fisher.p <- c(fisher.p, fisher.test(df1[i:(i+1), 1:2])$p.value)
}

df2 <- df1[seq(1, nrow(df1), 2), ]

df2$odds.ratio <- odd
df2$fisher.p.value <- fisher.p

# Perform Benjamini-Hochberg p value adjustment

df2$fisher.BH <- p.adjust(df2$fisher.p.value, method = 'BH')

write.csv(df2, file = 'Mortality Fisher exact testat at least 5 patients.csv')
