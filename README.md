## Benchmarking miRNA Target Prediction Tools

mirtargetbenchmark is an R package aimed towards benchmarking miRNA target prediction tools.

### Features

- Filer gene/miRNA expression data based on variance and frequency.
- Build regression models on the expression data.
- Convert predictions from tools into matrix format.
- Convert gene/miRNA ids into desired format for expression/tool prediction matrices.
- Extract common gene/miRNA interactions from the tools and regression coefficient matrix.
- Extract the miRNAs with more than certain vaidated gene targets from [miRTarBase](https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/php/index.php).
- Benchmark the tools using binding sites.
- Benchmark the tools using confidence scores.
- Ensemble method for target prediction.

### Installation

Install the package using:

`devtools::install_github("biomedbigdata/mirtargetbenchmark")`
