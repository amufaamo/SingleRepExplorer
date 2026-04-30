# SingleRepExplorer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18190984.svg)](https://doi.org/10.5281/zenodo.18190984)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An interactive visualization tool for single-cell analysis.

## ⚡ TL;DR (Quick Start)

Run this command in your terminal:
```bash
docker run --rm -p 3838:3838 amufaamo/singlerepexplorer:latest
```
Then open http://localhost:3838 and click **[Demo Data]** to try it out!



## 📖 Manual

For detailed documentation, please visit:

**[SingleRepExplorer Manual](https://amufaamo.github.io/SingleRepExplorer/manual.html)**

## 🚀 Overview

SingleRepExplorer is a web-based interactive tool for **integrated single-cell RNA-seq and immune repertoire analysis**. It combines transcriptomic analysis with TCR/BCR repertoire profiling in a single platform. You can try the application immediately using the built-in **Demo Data** without needing to upload your own files.

### Transcriptome Analysis
- UMAP/t-SNE visualization (with **Harmony** Batch Effect Correction)
- Gene expression plots (Violin, Feature, Dot, Heatmap)
- Differential expression analysis (Supports **Multi-group Kruskal-Wallis** and Repertoire-specific grouping)
- Powerful Enrichment Pathways (**EnrichR** & **GSEA** via MSigDB)
- Automated cell type annotation (**scType** Marker-based & **SingleR** Reference-based)
- **Trajectory Inference / Pseudotime** Analysis with clone tracking (Slingshot)
- **Phenotype-Repertoire Correlation Engine** to find clone size-correlated genes
- **Cell-Cell Communication Analysis (CellChat)** to infer signaling networks

### Repertoire Analysis (TCR/BCR)
- Clonotype frequency and abundance
- Diversity analysis
- Clonal overlap between samples
- Public clonotype identification
- Clonotype tracking across conditions

## 🐳 Docker Usage

The easiest way to run SingleRepExplorer is using Docker Hub:
```bash
docker run --rm -p 3838:3838 amufaamo/singlerepexplorer:v1.0.0
```

Then open your browser and go to: **http://localhost:3838**

### Building from source (optional)

```bash
# Build the image
docker build -t singlerepexplorer .

# Run the container
docker run -p 3838:3838 singlerepexplorer
```

### Analysis Start (Demo Data)
If you want to try the application without uploading files, click the blue **[Demo Data]** button.
This will automatically load the provided example dataset and run the analysis.
> Note: The [Demo Data] button is located below the [Run] button.

## 📂 Project Structure

```
SingleRepExplorer/
├── app/              # Shiny application
├── docs/             # GitHub Pages documentation
├── Dockerfile        # Docker image configuration
└── manual.Rmd        # Manual source
```

## 📝 License

This project is released under the [MIT License](LICENSE).

## 📚 Citation

If you use SingleRepExplorer in your research, please cite it via the Zenodo concept DOI, which always resolves to the latest version:

> Ishikawa, M. (2026). *SingleRepExplorer: An Interactive Web Application for Integrated Single-Cell Transcriptomic and Immune Repertoire Analysis* (Version 1.0.1) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.18190984

BibTeX:

```bibtex
@software{ishikawa_singlerepexplorer_2026,
  author       = {Ishikawa, Masakazu},
  title        = {{SingleRepExplorer: An Interactive Web Application
                   for Integrated Single-Cell Transcriptomic and
                   Immune Repertoire Analysis}},
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0.1},
  doi          = {10.5281/zenodo.18190984},
  url          = {https://doi.org/10.5281/zenodo.18190984}
}
```

Full machine-readable metadata is in [`CITATION.cff`](CITATION.cff); GitHub renders a "Cite this repository" button on the repository page based on it.

## 🔗 Links

- [GitHub Pages](https://amufaamo.github.io/SingleRepExplorer/)
- [Manual](https://amufaamo.github.io/SingleRepExplorer/manual.html)
