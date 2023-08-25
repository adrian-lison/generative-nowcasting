suppressMessages({
symptom_onset_legend_with_missing <- cowplot::get_legend(
  ggplot(data = data.frame(class = forcats::fct_inorder(c("Simulated truth","Reported until nowcast (known onset date)","Reported until nowcast (missing onset date)"), ordered = T), x = 1:3, y = rnorm(3)), aes(x = x, y = y)) +
    geom_col(aes(fill = class)) + geom_line(aes(color = class)) +
    theme_bw() +
    theme(legend.position = "left") + #legend.title = element_blank()
    scale_color_manual(values = c("black", "darkgrey", "lightgrey"), name = "Symptom onsets") +
    scale_fill_manual(values = c("white", "darkgrey", "lightgrey"), name = "Symptom onsets") #+
  #guides(fill = guide_legend(override.aes = list(alpha = 0.5, colour = c(NA, NA, NA))))
)

symptom_onset_legend <- cowplot::get_legend(
  ggplot(data = data.frame(class = forcats::fct_inorder(c("Simulated truth","Reported until nowcast",""), ordered = T), x = 1:3, y = rnorm(3)), aes(x = x, y = y)) +
    geom_col(aes(fill = class)) + geom_line(aes(color = class)) +
    theme_bw() +
    theme(legend.position = "left") +
    scale_color_manual(values = c("black", "darkgrey", "white"), name = "Symptom onsets") +
    scale_fill_manual(values = c("white", "darkgrey", "white"), name = "Symptom onsets")
)

R_legend <- cowplot::get_legend(
  ggplot(data = data.frame(class = forcats::fct_inorder(c("Simulated truth", " ", ""), ordered = T), x = 1:3, y = rnorm(3)), aes(x = x, y = y)) +
    geom_col(aes(fill = class)) + geom_line(aes(color = class), linetype = "solid") +
    theme_bw() +
    theme(legend.position = "left") +
    scale_color_manual(values = c("black", "white", "white"), name = "Reproduction number") +
    scale_fill_manual(values = c("white", "white", "white"), name = "Reproduction number")
)

pattern_legend <- cowplot::get_legend(
  ggplot(data = data.frame(class = forcats::fct_inorder(c("Overprediction", "Dispersion", "Underprediction"), ordered = T), x = 3, y = rnorm(3)), aes(x=x, y = y)) +
    geom_col_pattern(aes(pattern = class), fill = "white", colour = "black", position = "stack") +
    theme_bw() +
    scale_pattern_discrete(choices = c("stripe", "circle", "crosshatch"), name = "Nowcast performance") +
    theme(legend.position = "left") + guides(pattern=guide_legend(ncol=1,reverse = F))
)
})