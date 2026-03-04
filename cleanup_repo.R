# RStudio APIを使用せず、直接パスを指定して安全に実行するバージョン
setwd("C:/Users/amufa/OneDrive/SingleRepExplorer")

message("Cleaning up the repository structure...")

# 1. Rename .Rproj file correctly and move it to root
if (file.exists("app/SiCR.Rproj")) {
  file.rename("app/SiCR.Rproj", "SingleRepExplorer.Rproj")
}

# 2. Move Docker files from app/ to root
docker_files <- c("Dockerfile.arm64", "docker-compose.yml", "build_and_run.sh")
for (f in docker_files) {
  if (file.exists(file.path("app", f))) {
    file.rename(file.path("app", f), f)
  }
}

# 3. Move images to docs/images/ if they exist in app/
images <- list.files("app", pattern = "^Figure_.*\\.(png|jpg)$", full.names = TRUE)
if (length(images) > 0) {
  if (!dir.exists("docs/images")) dir.create("docs/images", recursive = TRUE)
  for (img in images) {
    file.rename(img, file.path("docs/images", basename(img)))
  }
}

# 4. Remove redundant files in app/
if (file.exists("app/README.md")) unlink("app/README.md")
if (dir.exists("app/docs")) unlink("app/docs", recursive = TRUE)

# 5. Remove garbage and massive cache files
garbage_files <- c(
  ".DS_Store",
  ".RData 2",
  ".RDataTmp",
  "SiCR - ショートカット.lnk",
  "build_log.txt",
  "app/analysis_cache.rds",
  "app/filename.RData",
  "app/myReactives_data.RData",
  "app/seurat_object.rds",
  "app/df.csv",
  "app/bcr_df.csv",
  "app/tcr_df.csv",
  "app/.DS_Store"
)

for (gf in garbage_files) {
  if (file.exists(gf)) {
    unlink(gf)
  }
}

message("Cleanup finished successfully!")
