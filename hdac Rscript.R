library(limma)
library(EnhancedVolcano)


rownames(data) <- data$ID
expr_data <- data[, c("KO1", "KO2", "KO3", "WT1", "WT2", "WT3")]


expr_log2 <- log2(expr_data + 1)


group <- factor(c("KO", "KO", "KO", "WT", "WT", "WT"))
design <- model.matrix(~ group)

fit <- lmFit(expr_log2, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "groupWT", number = Inf)


png("HDAC1_volcano_fixed.png", width = 2400, height = 2400, res = 300)

EnhancedVolcano(
  toptable = res,
  lab = rownames(res),
  x = 'logFC',
  y = 'adj.P.Val',
  pCutoff = 0.05,
  FCcutoff = 1.0,
  title = "HDAC1 Knockout vs Wild Type (GSE5583)",
  subtitle = "Log2 Transformed Differential Expression",
  labSize = 4.0,
  drawConnectors = TRUE
)

dev.off()