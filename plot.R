#!/usr/bin/env Rscript
suppressPackageStartupMessages(require("hotspotr"))
suppressPackageStartupMessages(require("ggplot2"))

get_colorscale <- function(low_colors = c("blue4", "blue3"), high_colors = c("red3", "red4"), inner_range = c("turquoise", "orange")) {
	cr <- colorRampPalette(inner_range)
	inner_colors <- cr(90)
	all_colors <- append(low_colors[1], append(rep(low_colors[2],4), inner_colors))

	all_colors <- append(append(all_colors, rep(high_colors[1],4)), high_colors[2])

	return(all_colors)
}

get_score_percentile <- function(scoredist, scores) {
	d <- ecdf(scoredist)
	p <- round(d(scores),2)
	p[p == 0] <- 0.01
	return(p)
}

make_hsmap <- function(data_df, scores, in_df) {
		
		cs <- get_colorscale()

		pval <- get_score_percentile(scores, data_df$score)

		data_df$pval <- pval
		data_df$score <- as.factor(pval)

		included_colors <- cs[sort(unique(pval*100))]

		return(hsmap(data_df, included_colors, in_df))
}

map_d <- read.csv("out.csv")
score_d <- read.csv("minmax.csv")
point_d <- read.csv("data.csv")

scores <- c(score_d$min, score_d$max)
hm <- make_hsmap(map_d, scores, point_d)

X11(type="cairo")
plot(hm)
message("Press Return To Continue")
invisible(readLines("stdin", n=1))