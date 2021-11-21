library(tibble)
library(data.table)

# Load patient and pathogen table

df1 <- read.csv('Allestablishedandpossiblepathogens_110420.csv', row.names = 1)

df1[is.na(df1)] <- 0

# Stratify by CD4 counts

low <- df1[df1$cd4rslt < 200, ]
high <- df1[df1$cd4rslt > 199, ]

low[nrow(low) + 1, ] <- colSums(low)
rownames(low)[rownames(low) == '167'] <- 'CD4.low'

high[nrow(high) + 1, ] <- colSums(high)
rownames(high)[rownames(high) == '52'] <- 'CD4.high'

low.counts <- low[167, 3:71]
high.counts <- high[52, 3:71]

all.counts <- rbind(low.counts, high.counts)

# Create 2 x 2 tables

low.no.path <- 167 - all.counts[1, 1:69]
high.no.path <- 52 - all.counts[2, 1:69]

no.path <- rbind(low.no.path, high.no.path)
colnames(no.path) <- paste('No', colnames(no.path), sep = '.')

all.counts.tables <- cbind(all.counts, no.path)[order(c(seq_along(all.counts), seq_along(no.path)))]

# Perform Fisher's exact test

values <- 1:dim(all.counts.tables)[2]
values <- values[c(TRUE, FALSE)]

odd <- vector()
fishp <- vector()

for (i in values)   
{
  odd <- c(odd, fisher.test(all.counts.tables[1:2, i:(i+1)])$estimate)
  fishp <- c(fishp, fisher.test(all.counts.tables[1:2, i:(i+1)])$p.value)
}

blank <- rep(NA, 69)

odd.combo <- rbind(matrix(odd, ncol = length(blank)), blank)
odd2 <- c(odd.combo)

all.counts.tables[nrow(all.counts.tables) + 1, ] <- odd2
rownames(all.counts.tables)[rownames(all.counts.tables) == '3'] <- 'Odds ratio'

fishp.combo <- rbind(matrix(fishp, ncol = length(blank)), blank)
fishp2 <- c(fishp.combo)

all.counts.tables[nrow(all.counts.tables) + 1, ] <- fishp2
rownames(all.counts.tables)[rownames(all.counts.tables) == '4'] <- 'Fisher.pvalue'

# Perform Benjamini-Hochberg p value adjustment

fishpadjust <- p.adjust(fishp, method = 'BH')
fishpadjust.combo <- rbind(matrix(fishpadjust, ncol = length(blank)), blank)
fishpadjust2 <- c(fishpadjust.combo)

all.counts.tables[nrow(all.counts.tables) + 1, ] <- fishpadjust2
rownames(all.counts.tables)[rownames(all.counts.tables) == '5'] <- 'Fisher.pvalue.BH'

all.counts.tables.transposed <- t(all.counts.tables)
all.counts.tables.transposed2 <- all.counts.tables.transposed[seq(1, nrow(all.counts.tables.transposed), 2), ]

write.csv(all.counts.tables.transposed2, file = 'CD4 Fisher exact test.csv')
