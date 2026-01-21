# SingleRepExplorer

An interactive visualization tool for single-cell analysis.

### Analysis Start (Demo Data)
If you want to try the application without uploading files, click the blue **[Demo Data]** button.
This will automatically load the provided example dataset and run the analysis.
> Note: The [Demo Data] button is located below the [Run] button.

## 📖 Manual

For detailed documentation, please visit:

**[SingleRepExplorer Manual](https://amufaamo.github.io/SingleRepExplorer/manual.html)**

## 🚀 Overview

SingleRepExplorer is a web-based interactive tool for **integrated single-cell RNA-seq and immune repertoire analysis**. It combines transcriptomic analysis with TCR/BCR repertoire profiling in a single platform. You can try the application immediately using the built-in **Demo Data** without needing to upload your own files.

### Transcriptome Analysis
- UMAP/t-SNE visualization
- Gene expression plots (Violin, Feature, Dot, Heatmap)
- Differential expression analysis
- Automated cell type annotation

### Repertoire Analysis (TCR/BCR)
- Clonotype frequency and abundance
- Diversity analysis
- Clonal overlap between samples
- Public clonotype identification
- Clonotype tracking across conditions

## 🐳 Docker Usage

The easiest way to run SingleRepExplorer is using Docker Hub:

```bash
docker run --rm -p 3838:3838 amufaamo/singlerepexplorer:latest
```

Then open your browser and go to: **http://localhost:3838**

### Building from source (optional)

```bash
# Build the image
docker build -t singlerepexplorer .

# Run the container
docker run -p 3838:3838 singlerepexplorer
```

## 📂 Project Structure

```
SingleRepExplorer/
├── app/              # Shiny application
├── docs/             # GitHub Pages documentation
├── Dockerfile        # Docker image configuration
└── manual.Rmd        # Manual source
```

## 📝 License

Please refer to the manual for details.

## 🔗 Links

- [GitHub Pages](https://amufaamo.github.io/SingleRepExplorer/)
- [Manual](https://amufaamo.github.io/SingleRepExplorer/manual.html)
