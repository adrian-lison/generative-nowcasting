suppressMessages({
symptom_onset_legend_with_missing <- cowplot::get_legend(
  ggplot(data = data.frame(class = forcats::fct_inorder(c("Consolidated estimate","Reported until nowcast (known onset date)",""), ordered = T), x = 1:3, y = rnorm(3)), aes(x = x, y = y)) +
    geom_col(aes(fill = class)) + geom_line(aes(color = class), linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "left") +
    scale_color_manual(values = c("black", "darkgrey", "white"), name = "Symptom onsets") +
    scale_fill_manual(values = c("white", "darkgrey", "white"), name = "Symptom onsets")
)

symptom_onset_legend <- cowplot::get_legend(
  ggplot(data = data.frame(class = forcats::fct_inorder(c("Consolidated estimate","Reported until nowcast",""), ordered = T), x = 1:3, y = rnorm(3)), aes(x = x, y = y)) +
    geom_col(aes(fill = class)) + geom_line(aes(color = class), linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "left") +
    scale_color_manual(values = c("black", "darkgrey", "white"), name = "Symptom onsets") +
    scale_fill_manual(values = c("white", "darkgrey", "white"), name = "Symptom onsets")
)

R_legend <- cowplot::get_legend(
  ggplot(data = data.frame(class = forcats::fct_inorder(c("Consolidated estimate", " ", ""), ordered = T), x = 1:3, y = rnorm(3)), aes(x = x, y = y)) +
    geom_col(aes(fill = class)) + geom_line(aes(color = class), linetype = "dashed") +
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