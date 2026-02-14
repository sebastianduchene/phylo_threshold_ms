setwd("~/Dropbox/projects_WORKING/temporal_signal_comment/temporal_signal_comment/temporal_sampling_bias_settings_3")
#setwd("~/Dropbox/projects_WORKING/temporal_signal_comment/temporal_signal_comment/temporal_sampling_bias_settings_4")

library(RColorBrewer)

pal <- brewer.pal(9, "Blues")


library(NELSI)

complete_trees <- read.tree("complete_sampling_trees.trees")
balanced_trees <- read.tree("balanced_sampling_trees.trees")
unbalanced_trees <- read.tree("unbalanced_sampling_trees.trees")

table(round(allnode.times(complete_trees[[17]], tipsonly = T, reverse = F), 2))
table(round(allnode.times(balanced_trees[[17]], tipsonly = T, reverse = F), 2))
table(round(allnode.times(unbalanced_trees[[17]], tipsonly = T, reverse = F), 2))


pdf("~/Dropbox/projects_WORKING/temporal_signal_comment/phylo_threshold_ms/sampling_bias_summary_trees_rates.pdf", useDingbats = FALSE, width = 11, height = 14)
par(mfrow = c(2, 1), mar = c(5, 5, 4, 1))
index = 17
plot(0, 0, xlim = c(0, 850), ylim = log10(c(1, 8000)), type = "n", bty = "n", 
     xaxt = "n", yaxt = "n", xlab = "Sampling scheme", ylab = "Time (years)", 
     main = expression("(a) Trees with temporal sampling"), 
     cex.lab = 1.3)
plot.tree.lines(complete_trees[[index]], plot.new = F, rotation.angle = pi*3/2, line.type = "l",  # nolint: line_length_linter.
                tips.colour = "grey", log.scale = T, lines.colour = pal[9])                
th <- log10(unique(round(allnode.times(complete_trees[[index]], tipsonly = T, reverse = F)), 3)) # nolint
for(i in 2:length(th)){
  lines(c(-20, 850), c(th[i], th[i]), lty = 3)
}
text(255, 0.05, label = "n=100")
text(255, th[2]+0.05, label = "n=100")
text(330, th[3]+0.05, label = "n=100")
text(420, th[4]+0.05, label = "n=100")
text(520, th[5]+0.05, label = "n=100")
#
lines(c(-20, 850), rep(0, 2), lty = 3)
times <- c(th, log10(max(allnode.times(complete_trees[[index]]))))
times[which.min(times)] <- 0
time_labels <- round(10^times)
time_labels[which.min(time_labels)] <- 0
axis(2, at = times, labels = time_labels, las = 2)
axis(1, at = c(240, 550, 730), labels = c("Complete data\n(n=500)", 
                                     "Time-uniform\n(n=100)", "Time-biased\n(n=100)"), tick = FALSE)
plot.tree.lines(balanced_trees[[index]], plot.new = F, rotation.angle = pi*3/2, line.type = "l", 
                x.offset = 550, tips.colour = "grey", log.scale = T, lines.colour = pal[6])
text(625, 0.05, label = "n=20")
text(625, th[2]+0.05, label = "n=20")
text(630, th[3]+0.05, label = "n=20")
text(655, th[4]+0.05, label = "n=20")
text(665, th[5]+0.05, label = "n=20")
plot.tree.lines(unbalanced_trees[[index]], plot.new = F, rotation.angle = pi*3/2, line.type = "l",
                x.offset = 700, tips.colour = "grey", log.scale = T, lines.colour = pal[3])
text(825, 0.05, label = "n=90")
text(825, th[2]+0.05, label = "n=5")
text(825, th[3]+0.05, label = "n=3")
text(825, th[4]+0.05, label = "n=1")
text(825, th[5]+0.05, label = "n=1")
############################# Adding the labels for number of tips at each epoch
##############################

th <- unique(round(allnode.times(balanced_trees[[index]], tipsonly = T, reverse = F)), 3)

#############################
############################# Now plot summary of estimates

complete_trees <- read.tree("~/Dropbox/projects_WORKING/temporal_signal_comment/temporal_signal_comment/temporal_sampling_bias_settings_3/complete_sampling_trees.trees")
balanced_trees <- read.tree("~/Dropbox/projects_WORKING/temporal_signal_comment/temporal_signal_comment/temporal_sampling_bias_settings_3/balanced_sampling_trees.trees")
unbalanced_trees <- read.tree("~/Dropbox/projects_WORKING/temporal_signal_comment/temporal_signal_comment/temporal_sampling_bias_settings_3/unbalanced_sampling_trees.trees")

# Just checking trees and alignments"
fasta_complete <- dir("~/Dropbox/projects_WORKING/temporal_signal_comment/temporal_signal_comment/temporal_sampling_bias_settings_3/", 
                         pattern = "complete.+fasta")
for(i in seq_along(complete_trees)){
     print(i)
     print(all(balanced_trees[[i]]$tip.labels %in% complete_trees[[i]]$tip.labels))
     print(all(unbalanced_trees[[i]]$tip.labels %in% complete_trees[[i]]$tip.labels))
     aln <- read.dna(fasta_complete[i], format = "fasta")
     print(all(rownames(aln[i]) %in% complete_trees[[i]]$tip.labels))
}

get_stats_ci <- function(log_file_name, tree){
  posterior <- tail(read.table(log_file_name, header = TRUE), 5000)
  clock_rate <- posterior$rate.mean
  true_tree_height <- max(allnode.times(tree))
  tree_height <- (true_tree_height - posterior$tree.height) / true_tree_height
  clock_rate_ci <- quantile(clock_rate, c(0.0125, 0.9875))
  tree_height_ci <- quantile(tree_height, c(0.0125, 0.9875))
  return(c(mean(clock_rate), clock_rate_ci, mean(tree_height), tree_height_ci))
}

log_files_dir <- "~/Dropbox/projects_WORKING/temporal_signal_comment/temporal_signal_comment/temporal_sampling_bias_settings_3/"
complete_log_files <- dir(pattern = "complete.+log")
balanced_log_files <- dir(pattern = "balanced.+log")
unbalanced_log_files <- dir(pattern = "unbalanced.+log")

results_table <- matrix(NA, nrow = length(complete_log_files), ncol = 18)
colnames(results_table) <- c("complete_clock_rate_mean", "complete_clock_rate_lower", 
                              "complete_clock_rate_upper", 
                              "complete_tree_height_mean", "complete_tree_height_lower", 
                              "complete_tree_height_upper", 
                              "balanced_mean_rate_mean", "balanced_mean_rate_lower", 
                              "balanced_mean_rate_upper", 
                              "balanced_tree_height_mean", "balanced_tree_height_lower", 
                              "balanced_tree_height_upper", 
                              "unbalanced_mean_rate_mean", "unbalanced_mean_rate_lower", 
                              "unbalanced_mean_rate_upper",
                              "unbalanced_tree_height_mean", "unbalanced_tree_height_lower", 
                              "unbalanced_tree_height_upper"
                              )
for(i in seq_along(complete_log_files)){
     complete_stats <- get_stats_ci(paste0(log_files_dir, complete_log_files[i]), complete_trees[[i]])
     balanced_stats <- get_stats_ci(paste0(log_files_dir, balanced_log_files[i]), balanced_trees[[i]])
     unbalanced_stats <- get_stats_ci(paste0(log_files_dir, unbalanced_log_files[i]), unbalanced_trees[[i]])
     results_table[i, ] <- c(complete_stats, balanced_stats, unbalanced_stats)
     print(results_table[i, ])
}
########################################################################################
results_table <- results_table[order(results_table[, 1]), ]
plot(0, 0, xlim = c(0, 101), ylim = log10(c(7e-6, 5e-5)), type = "n", xaxt = "n", 
     yaxt = "n", bty = "n", ylab = expression("Evol. rate"~10^{-5}~"subs/site/year"), 
     cex.lab = 1.3, xlab = "Simulation replicate", 
     main = expression("(b) Posterior of evol. rate for sampling schemes"))
axis(2, at = log10(c(7e-6, 1.5e-5, 3e-5, 5e-5)), labels = c(0.7, 1.5, 3, 5), las = 2)
axis(1, at = c(0, 101), labels = c("", ""))
lines(c(-5, 101), rep(log10(1.5e-5), 2), lty = 2, lwd = 2)
par(mfrow = c(2, 1), mar = c(2, 4.3, 4, 1))
results_table <- results_table[order(results_table[, 1]), ]
for(i in 1:(dim(results_table)[1])){
     points(i-0.2, log10(results_table[i, 13]), pch = 20, col = pal[3])
     lines(rep(i-0.2, 2), log10(results_table[i, 14:15]), col = pal[3], lwd = 2)
     points(i, log10(results_table[i, 7]), pch = 20, col = pal[6])
     lines(rep(i, 2), log10(results_table[i, 8:9]), col = pal[6], lwd = 2)
     points(i+0.2, log10(results_table[i, 1]), pch = 20, col = pal[9])
     lines(rep(i+0.2, 2), log10(results_table[i, 2:3]), col = pal[9], lwd = 2)
}
dev.off()

## Plot uncertainties
width_rate_complete <- (results_table[, "complete_clock_rate_upper"] - results_table[, "complete_clock_rate_lower"]) / results_table[, "complete_clock_rate_mean"]
width_rate_balanced <- (results_table[, "balanced_mean_rate_upper"] - results_table[, "balanced_mean_rate_lower"]) / results_table[, "balanced_mean_rate_mean"]
width_rate_unbalanced <- (results_table[, "unbalanced_mean_rate_upper"] - results_table[, "unbalanced_mean_rate_lower"]) / results_table[, "unbalanced_mean_rate_mean"]

width_tree_height_complete <- results_table[, 6] - results_table[, 5]
width_tree_height_balanced <- results_table[, 12] - results_table[, 11]
width_tree_height_unbalanced <- results_table[, 18] - results_table[, 17]


pdf("~/Dropbox/projects_WORKING/temporal_signal_comment/phylo_threshold_ms/sampling_bias_summary_rates.pdf", useDingbats = FALSE, width = 12, height = 7)
par(mfrow = c(1, 2), mar = c(4, 5, 5, 1))
plot(results_table[, "complete_clock_rate_mean"], results_table[, "balanced_mean_rate_mean"], 
     xlim = c(1.2e-5, 2.4e-5), ylim = c(1.2e-5, 2.4e-5), 
     col = pal[6], pch = 20, bty = "n", xlab = "Complete data (n=500)", 
     ylab = "Sampling (n=100)", cex.lab = 1.3, 
     main = expression("(a) Mean posterior evol. rate"~10^{-5}~"subs/site/year"), 
     cex.main = 1.3, cex = 1.5, xaxt= "n", yaxt = "n")
axis(1, at = c(1.2e-5, 1.5e-5, 1.8e-5, 2.1e-5, 2.4e-5), 
     labels = c(1.2, 1.5, 1.8, 2.1, 2.4))
axis(2, at = c(1.2e-5, 1.5e-5, 1.8e-5, 2.1e-5, 2.4e-5), 
     labels = c(1.2, 1.5, 1.8, 2.1, 2.4), las = 2)
lines(c(1.2e-5, 2.4e-5), c(1.2e-5, 2.4e-5), lty = 2)
points(results_table[, "complete_clock_rate_mean"], results_table[, "unbalanced_mean_rate_mean"], 
     pch = 20, col = pal[3], cex = 1.5)
legend("bottomright", legend = c("Time-uniform", "Time-biased"), 
       col = c(pal[6], pal[3]), pch = 20, bty = "n", cex = 1)
plot(width_rate_complete, width_rate_balanced,
     xlim = c(0.2, 0.6), ylim = c(0.2, 0.6), 
     col = pal[6], pch = 20, bty = "n", xlab = "Complete data (n=500)", 
     ylab = "", cex.lab = 1.3, xaxt = "n", yaxt = "n",
     main = expression("(b) Relative width of the 95% credible interval\nof the evol. rate"), 
     cex.main = 1.3, cex = 1.5)
axis(1, at = c(0.2, 0.3, 0.4, 0.5, 0.6),
     labels = c(0.2, 0.3, 0.4, 0.5, 0.6))
axis(2, at = c(0.2, 0.3, 0.4, 0.5, 0.6),
     labels = c(0.2, 0.3, 0.4, 0.5, 0.6), las = 2)
lines(c(0.2, 0.6), c(0.2, 0.6), lty = 2)
points(width_rate_complete, width_rate_unbalanced, 
     pch = 20, col = pal[3], cex = 1.5)
dev.off()

## Check bias
head(results_table)

mean( (1.5e-5 > results_table[, "complete_clock_rate_lower"]) & (1.5e-5 < results_table[, "complete_clock_rate_upper"]) )

mean( (1.5e-5 > results_table[, "balanced_mean_rate_lower"]) & (1.5e-5 < results_table[, "balanced_mean_rate_upper"]) )

mean( (1.5e-5 > results_table[, "unbalanced_mean_rate_lower"]) & (1.5e-5 < results_table[, "unbalanced_mean_rate_upper"]) )

## Check 
mean(  (results_table[, "balanced_mean_rate_mean"] - results_table[, "complete_clock_rate_mean"]) / results_table[, "balanced_mean_rate_mean"])

mean( results_table[, "unbalanced_mean_rate_mean"] - results_table[, "complete_clock_rate_mean"]  / results_table[, "complete_clock_rate_mean"] )


mean(results_table[, "balanced_mean_rate_mean"] - results_table[, "complete_clock_rate_mean"])

mean( results_table[, "unbalanced_mean_rate_mean"] - results_table[, "complete_clock_rate_mean"])


mean(  (results_table[, "balanced_mean_rate_mean"] > results_table[, "complete_clock_rate_mean"]))

mean(  (results_table[, "unbalanced_mean_rate_mean"] > results_table[, "complete_clock_rate_mean"]))


## Calculate coverage
coverage_complete <- sum( (1.5e-5 > results_table[, "complete_clock_rate_lower"]) & (1.5e-5 < results_table[, "complete_clock_rate_upper"]) )
coverage_balanced <- sum( (1.5e-5 > results_table[, "balanced_mean_rate_lower"]) & (1.5e-5 < results_table[, "balanced_mean_rate_upper"]) )
coverage_unbalanced <- sum( (1.5e-5 > results_table[, "unbalanced_mean_rate_lower"]) & (1.5e-5 < results_table[, "unbalanced_mean_rate_upper"]) )
coverage_complete
coverage_balanced
coverage_unbalanced

bias_complete <- mean((results_table[, "complete_clock_rate_mean"] - 1.5e-5) / 1.5e-5)
bias_balanced <- mean((results_table[, "balanced_mean_rate_mean"] - 1.5e-5) / 1.5e-5)
bias_unbalanced <- mean((results_table[, "unbalanced_mean_rate_mean"] - 1.5e-5) / 1.5e-5)
bias_complete
bias_balanced
bias_unbalanced

uncertainty_complete <- mean((results_table[, "complete_clock_rate_upper"] - results_table[, "complete_clock_rate_lower"]) / results_table[, "complete_clock_rate_mean"])
uncertainty_balanced <- mean((results_table[, "balanced_mean_rate_upper"] - results_table[, "balanced_mean_rate_lower"]) / results_table[, "balanced_mean_rate_mean"])
uncertainty_unbalanced <- mean((results_table[, "unbalanced_mean_rate_upper"] - results_table[, "unbalanced_mean_rate_lower"]) / results_table[, "unbalanced_mean_rate_mean"])
uncertainty_complete
uncertainty_balanced
uncertainty_unbalanced

# Make a table for each, coverage, bias, and uncertainty in each row
summary_table <- matrix(c(coverage_complete, bias_complete, uncertainty_complete,
                          coverage_balanced, bias_balanced, uncertainty_balanced,
                          coverage_unbalanced, bias_unbalanced, uncertainty_unbalanced), 
                        nrow = 3, byrow = T)
colnames(summary_table) <- c("Coverage", "Bias", "Uncertainty")
rownames(summary_table) <- c("Complete data (n=500)", "Time-uniform sampling (n=100)", "Time-biased sampling (n=100)")
summary_table <- t(summary_table)

#print table in latex style
library(xtable)
print(xtable(summary_table, digits = 3), include.rownames = T, include.colnames = T)

