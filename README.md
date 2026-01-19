# SingleRepExplorer

An interactive visualization tool for single-cell analysis.

## 📖 Manual

For detailed documentation, please visit:

**[SingleRepExplorer Manual](https://amufaamo.github.io/SingleRepExplorer/manual.html)**

## 🚀 Overview

SingleRepExplorer is a web-based interactive single-cell analysis tool built with Shiny.
Using Seurat objects, it provides the following analyses and visualizations:

- UMAP plots
- Violin plots
- Feature plots
- Dot plots
- Cluster analysis

## 🐳 Docker Usage

This project can be run as a Docker container.

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
