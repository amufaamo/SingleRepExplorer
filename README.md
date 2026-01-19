# SingleRepExplorer

SingleRepExplorerは、シングルセル解析のためのインタラクティブな可視化ツールです。

## 📖 マニュアル

詳細なマニュアルは以下のリンクからご覧いただけます：

**[SingleRepExplorer マニュアル](https://amufaamo.github.io/SingleRepExplorer/manual.html)**

## 🚀 概要

SingleRepExplorerは、Shinyを使用したWebベースのインタラクティブなシングルセル解析ツールです。
Seuratオブジェクトを使用して、以下のような解析と可視化が可能です：

- UMAPプロット
- バイオリンプロット
- フィーチャープロット
- ドットプロット
- クラスター解析

## 🐳 Dockerでの使用

このプロジェクトはDockerコンテナとして実行できます。

```bash
# イメージのビルド
docker build -t singlerepexplorer .

# コンテナの実行
docker run -p 3838:3838 singlerepexplorer
```

## 📂 プロジェクト構造

```
SingleRepExplorer/
├── app/              # Shinyアプリケーション
├── docs/             # GitHub Pagesドキュメント
├── Dockerfile        # Dockerイメージ設定
└── manual.Rmd        # マニュアルソース
```

## 📝 ライセンス

詳細はマニュアルをご確認ください。

## 🔗 リンク

- [GitHub Pages](https://amufaamo.github.io/SingleRepExplorer/)
- [マニュアル](https://amufaamo.github.io/SingleRepExplorer/manual.html)
