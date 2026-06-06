ZITM: Zero-Inflated Telegraph Model
Single-cell RNA-seq differential expression analysis based on Zero-Inflated Telegraph Model

📦 Installation

install.packages("devtools")
devtools::install_github("Always-Stude/ZITM")

library(ZITM)

🚀 Quick Start Example

Simulate single-cell count matrix

mat <- matrix(
sample(0:50, size = 100*100, replace = TRUE),
nrow = 100, ncol = 100,
dimnames = list(paste0("gene",1:100), paste0("cell",1:100))
)

Experimental grouping

group <- factor(c(rep(1,50), rep(2,50)))
Core differential expression analysis

results <- Differential_analysis(counts = mat, group = group)

head(results)

📖 Documentation

?Differential_analysis

example(Differential_analysis)

📄 License

MIT License