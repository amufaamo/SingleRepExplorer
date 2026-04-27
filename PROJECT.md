# SingleRepExplorer — Project Overview

## プロジェクト概要 / Project Summary

**SingleRepExplorer** は、プログラミングに精通していない生物学の研究者が、シングルセルRNA-seq解析（特に**シングルセルパスウェイ解析**や**シングルセル免疫学解析**）を自力で実施し、**論文を出版するまで**のすべての工程をノーコードで完結できることを目指したウェブアプリケーションです。

R Shiny + Docker を基盤とし、ブラウザ上で直感的に操作できるGUIを提供します。

---

## 目的 / Objectives

### 主目的
プログラミングスキルが限られている生物学研究者が、シングルセル免疫学解析およびパスウェイ解析を直感的に行い、トップレベルの学術誌に論文を投稿・出版できるようにサポートする。

### 対象ユーザー
- **プログラミングに精通していない生物学・医学研究者**
- シングルセル解析に取り組む免疫学研究者
- 臨床研究者・大学院生

### 最終目標
**Frontiers in Immunology** (Methods or Brief Research Report section) への投稿・出版。
バイオインフォマティクス専門誌ではなく、免疫学の最前線で戦う研究者が「ツールとして使いやすい」ことを証明し、アクセプトを勝ち取る。

---

### 解析機能 / Analysis Features

#### トランスクリプトーム・パスウェイ解析
- **シングルセルパスウェイ解析（Pathway Analysis）**: AUCell, UCell, または GSEA を用いた経路活性化スコアリング
- UMAP / t-SNE 次元削減（Harmony バッチ効果補正対応）
- 遺伝子発現プロット（Violin, Feature, Dot, Heatmap）
- 差次的発現解析（Multi-group Kruskal-Wallis、レパトアグループ対応）
- エンリッチメント解析（EnrichR & GSEA / MSigDB）
- 細胞型アノテーション（scType マーカーベース・SingleR リファレンスベース）
- 軌道推論 / 擬似時間解析（Slingshot）＋クローントラッキング
- 細胞間コミュニケーション解析（CellChat）

#### 免疫レパトア解析（Single-cell Immunology）
- クローンサイズ相関遺伝子探索（Phenotype-Repertoire Correlation Engine）
- クロノタイプ頻度・存在量
- 多様性指標解析
- サンプル間クローン重複
- パブリッククローン同定
- 条件間クロノタイプトラッキング

---

## 技術スタック / Tech Stack

| 要素 | 内容 |
|------|------|
| フロントエンド | R Shiny（ノーコードGUI） |
| バックエンド | R（Seurat, scRepertoire, Slingshot, CellChat 等） |
| インフラ | Docker / Docker Hub (`amufaamo/singlerepexplorer`) |
| ドキュメント | GitHub Pages（`amufaamo.github.io/SingleRepExplorer`） |
| データ形式 | 10x Genomics CellRanger出力（.h5）、VDJ CSV |

---

## 投稿状況 / Submission History

| 2025年 | BMC Bioinformatics（Software article） | **Reject** |
| 2026年 | **Frontiers in Immunology** | 再投稿準備中 |

### 再投稿に向けた戦略（Frontiers in Immunology ターゲット）
- **ターゲット部門**: *Technology and Code* を第一候補とする。免疫学者が「プログラミングを学ばなくても、自分のデータを解析して図を作れる」という実用性を実証した新しいソフトウェア・コードとして投稿する。
- **リジェクト理由の克服**: 前回の査読指摘を反映し、UIの安定性、解析手法の透明性（使用パッケージの明示）、および出力データの完全性を向上させる。
- **ターゲット候補**: **Frontiers in Immunology (第一候補)**, iScience, Bioinformatics.

---

## 開発状況 / Development Status

- 現行バージョン：`v1.0.0`（論文投稿用の最初の正式版。Docker Hubで `amufaamo/singlerepexplorer:v1.0.0` および `:latest` として公開予定）
- デモデータ：HCWワクチン接種コホート（COVID-19）データ内蔵
- 前バージョン：SiCR（SiCR2）→ SingleRepExplorer にリブランド

---

## リポジトリ構成 / Repository Structure

```
ruft_SingleRepExplorer/
├── app/              # Shiny アプリ本体
│   ├── app.R         # メインアプリ
│   ├── utils.R       # ユーティリティ関数
│   └── module/       # 各解析モジュール
├── docs/             # GitHub Pages ドキュメント
├── data/             # サンプルデータ・データベース
├── example/          # デモデータ
├── Dockerfile        # Docker イメージ設定
├── PROJECT.md        # このファイル（プロジェクト目標・概要）
└── README.md         # クイックスタートガイド
```

---

## 起動方法 / How to Launch (Windows / Docker)

Windows 環境（WSL推奨）でアプリを起動するには、以下の手順に従ってください。

### 方法 A: 自動スクリプトを使用する（推奨）
プロジェクトのルートにある `build_and_run.sh` を使用すると、既存コンテナの掃除からビルド、起動までを一括で行えます。

1.  **ターミナル（WSL等）を開く**
2.  **スクリプトを実行する**
    ```bash
    # 基本：ソースからビルドして起動
    ./build_and_run.sh

    # Docker Hub から最新イメージを pull して起動する場合
    ./build_and_run.sh --pull
    ```

### 方法 B: Docker Compose を直接使用する
1.  **Docker コンテナを起動する**
    ```bash
    docker-compose up -d
    ```

### 起動後の確認
- 起動したら、ブラウザで以下の URL にアクセスしてください：
  - [http://localhost:3838/](http://localhost:3838/)
- ログを確認する場合： `docker logs -f singlerepexplorer`

---

## 推奨される開発・リリースの流れ / Recommended Workflow

1.  **手元でテスト (Local Test)**
    ```bash
    ./build_and_run.sh
    ```
    - これにより、最新のソースコードからイメージがビルドされ、ローカルでコンテナが起動します。
    - ブラウザで `http://localhost:3838/` にアクセスして動作を確認します。
2.  **公開 (Public Release)**
    手元での動作が完璧であれば、Docker Hub へ公開します：
    ```bash
    # 例: v1.0.0 としてリリース
    ./build_and_push.sh v1.0.0
    ```
    - これにより、指定したバージョンタグと `latest` タグが Docker Hub にプッシュされます。
    - 以降、世界中のユーザーが `docker pull` や `./build_and_run.sh --pull` で最新版を利用できるようになります。

> [!TIP]
> **なぜ push 機能が分離されているのか？**
> 誤ってバグを含んだ状態で公開してしまうのを防ぎ、かつ手元での素早い動作確認を可能にするためです。「手元で確認（Test）してから公開（Push）」の 2 ステップを基本とします。

---

## 開発・リリース用スクリプト / Development & Release Scripts

管理・リリースのために以下のスクリプトが用意されています。これらは `app/app.R` 内のバージョン番号を自動的に取得してタグ付けを行います。

- **`build_and_run.sh`**: ローカル環境でのテスト・実行用。アーキテクチャ(x86_64/arm64)を自動判別します。
- **`build_and_push.sh`**: 新しいバージョンをビルドし、Docker Hub (`amufaamo/singlerepexplorer`) にプッシュします。
- **`build_and_push_auto.sh`**: 上記の簡易自動化版です。

---

---

---
## アプリ修正・改善状況 / Code Review & Fixes Status

> 2026-03-11のコードレビューに基づき、以下の修正を完了しました。

### ✅ 完了済みの修正 (Done)

1. **モジュール統合と依存関係の整理** (2026-03-11)
...
4. **リポジトリの軽量化** (2026-03-11)

5. **LLM Annotation のプロンプト最適化 & 機能拡張** (2026-03-12)
   - OpenAIへのプロンプトに「Species (種)」「Tissue (組織)」「Granularity (粒度)」のコンテキストを追加し、アノテーション精度を向上させました。
   - UIにこれらの値を指定できる入力項目を追加しました。
   - `app.R` に LLM Annotation モジュールが正しく `source()` および UI 統合されていなかった問題を修正し、正式にメニューに追加しました。

6. **CellChat の環境構築とUI表示の最適化** (2026-03-12)
   - Dockerfile に `CellChat` v2 がデフォルトで含まれていることを確認し、アプリ側の警告メッセージを「公式イメージ v2.1.0+ には内蔵済み」である旨に更新しました。

7. **ドキュメント (manual.Rmd) の大幅更新** (2026-03-12)
   - `manual.Rmd` に、これまで詳細な説明が不足していた以下の新機能の解説を追加しました：
     - Trajectory & Pseudotime (Slingshot)
     - Phenotype-Repertoire Correlation
     - Cell-Cell Communication (CellChat)
     - Multimodal Integration (CoNGA)
     - Automated Cell Type Annotation (LLM)
     - CDR3 Physicochemical Profiling
     - Sequence Similarity
     - BCR Lineage Tree
     - BCR SHM & Isotype

8. **Editable PowerPoint (.pptx) Plot Download 実装** (2026-03-16)
   - 全てのモジュールのプロットダウンロードを PDF から PowerPoint (.pptx) に変更しました。
   - `officer` + `rvg` (DrawingML) を採用し、PowerPoint 上で「グループ化解除」することで、フォント、色、線の太さなどを個別に編集可能です。
   - `utils.R` に共通ヘルパー関数 `save_plot_as_pptx` および `save_baseplot_as_pptx` を実装しました。
   - 全モジュールのダウンロードボタンをプロットの上に統一配置し、UI の一貫性を向上させました。

9. **Dockerコンテナのアップデート & バージョン整合 (v2.2.0)** (2026-03-17)
   - `app/app.R` 内のバージョン表記を v2.2.0 に更新。
   - Docker Hub へのイメージプッシュとコンテナ再起動を実施。
   - **[Tips] WSL環境でGoogle Driveマウントが外れた場合のコンテナ起動手順:**
     WSL上で「No such device」と出た場合は、以下の手順で手動で再マウントし、既存コンテナを削除してから再起動する。
     ```bash
     sudo umount /mnt/g
     sudo mount -t drvfs G: /mnt/g
     cd "/mnt/g/マイドライブ/ruft_SingleRepExplorer"
     docker rm -f singlerepexplorer
     docker-compose up -d
     ```

10. **Frontiers in Immunologyへの論文投稿準備 (2026-03-17 ~ 2026-03-23)**
    - `Frontiers in Immunology` の最新規定（Article Types）を調査し、投稿フォーマットの要件を `Frontiers_in_Immunology_Guidelines.md` にまとめた（投稿前の最終チェック時にはこのファイルを参照すること）。
    - 前回のBMC BioinformaticsのReject体験をふまえ、ツールの単なる紹介ではなく、新機能（LLM細胞アノテーション、CoNGA、CellChatなど）を活用して「シングルセルマルチオミクスデータからいかに免疫学的洞察を得られるか」を強調した論文草稿を `paper_FiI.md` に作成。
    - **方針の決定**: Article Type としては、新しいソフトウェアプラットフォームの発表に最も適した **「Technology and Code」** で投稿することを決定した。
    - それに伴い、草稿（`paper_FiI.md`）のアブストラクトの単一段落化、セクション番号の削除、Data Availability Statementの追加など、同カテゴリーのガイドラインに完全準拠するよう修正を完了。
11. **次期課題の追加対応 (2026-03-17)**
    - CellChat 実行時に `object 'CellChat' not found` が出る環境向けに、実行前に `CellChat` を明示ロードするように修正。
    - CoNGA Integration のダウンロードボタンを各プロットの上部に移動。
    - CDR3 Physicochemical Properties の計算時エラーを、欠損列や入力不足時に明確なエラーメッセージを表示するよう改善。
    - Repertoire Similarity の距離計算エラー（Hamming長さ不一致など）にフォールバック対応を追加。
    - Track Clonotype のプロットサイズ入力を `numericInput` に統一し、表示メッセージを英語化。
    - BCR Lineage Tree でクローン選択ができない問題を、全クローン選択可能にして解消。
    - BCR SHM & Isotype の表示が空になるケースで、原因が分かるバリデーションメッセージを追加。

12. **WSL Google Driveマウント不具合の修正とSkill化 (2026-03-19)**
    - WSLターミナルから `/mnt/g` 以下へのアクセスが `failed 19 (ENODEV)` で失敗する問題を、`sudo mount -t drvfs G: /mnt/g` による再マウントで解決。
    - この知見をプロジェクト内（`.agent/skills`）およびグローバル（`dotfiles` リポジトリの `CLAUDE.md`, `GEMINI.md`, `skills`）に反映。
    - WSLの `~/.bashrc` に起動時のマウントチェック関数とエイリアス `mount_g` を追加し、再発時の自動・半自動復旧を可能にした。
    - Google Drive内での `git clone` がパーミッションエラーで失敗することを踏まえ、Windows側のPowerShellを用いたクローンの手順を確立。

## 📝 次期アップデート・修正計画 (from `next.md`)

`next.md` に挙げられた20項目の改善要望を以下の3つのフェーズに分類し、順次対応する。

### フェーズ 1: エラー修正・バグ対応 (Bug Fixes)
- [x] **12.** `CellChat`モジュールのエラー対応（`Cell labels cannot contain '0'!`）→ クラスター名にプレフィックスを付けるなどで対処。
- [x] **15.** `Diversity`のテーブルがTab表示になっている＆エラーが表示される問題の修正。
- [x] **16.** `CDR3 Physicochemical Properties`のエラー（何も表示されない）の修正。
- [x] **17.** `Repertoire Similarity`のエラー（何も表示されない）の修正。
- [x] **19.** `BCR Lineage Tree`がボタンを押しても何も出てこない問題の修正。
- [x] **20.** `BCR SHM & Isotype`が何も出てこない問題の修正。
- [x] **5.** `Feature Plot` 初期表示の間延び問題（描画タイミングまたはPlotHeightの初期値の設定）の修正。
- [x] **9.** `Trajectory`モジュールにおけるColor Paletteが反映されない問題の修正。

### フェーズ 2: UI/UXの統一と見栄えの改善 (UI/UX Consistency)
- [x] **1.** Stepではないタイトル（Data Processing & UMAPなど）をStepタイトルと色分けして分かりやすくする。
- [x] **10.** UI設定パネルのタイトルを「Visualization Settings」や「Plot option」から、Plotに対する設定を明示する名前（例: Plot Options）に統一する。
- [x] **11.** `Base Font Size`, `Facet Label Size` をすべてのPlotに対して設定できるようにする。
- [x] **13.** `Bar Border Color` オプションの削除（不要なため）。
- [x] **14.** Repertoire側（Unique clonotype, Frequency analysisなど）のPlot色をSeuratのCluster色と統一する。
- [x] **18.** `Track clonotype` モジュールのDownload Plotボタンを、他モジュールに合わせてプロットの上部に配置。
- [x] **8.** `Trajectory Inference` で機能していないPseudotime（擬似時間）の記述やオプションを整理・削除。
- [x] **6b.** `Gene Expression` モジュールのFeature Plot/Violin Plot等について、ドロップダウンからタブ表示へと変更し、個別にDownloadボタンを配置。

### フェーズ 3: 機能拡張とインタラクティビティ (Feature Enhancements)
- [x] **6.** Clustering実行後も、QC Plotsをコントロール（フィルタリング再設定など）したまま、再度Runできるようにする。
- [x] **7.** `DEG`モジュールの Top Markers Expression で、DotPlotとHeatmapを選べるようにし、他のPlotと同様のオプション（Plot Height等）を付加する。
- [x] **2.** Dimensional Plotで、Lasso select等でSubsetting/Highlightした際、Lasso用のUMAP上にもHighlightが反映されるようにする。
- [x] **3.** Highlightして特定の細胞だけ表示し、それを基にSubsetting & Clusteringを行えるようにする。
- [x] **4.** Highlightした細胞をさらにLassoで囲み、名前変更やSubsetting & Clusteringを行えるようにする。

### フェーズ 4: 追加バグ修正とUI改善 (Additional Fixes - next.md)
- [x] **1.** `CellChat`モジュールの実行時に `Error: object 'CellChat' not found` が出る問題の修正。
- [x] **2.** `CoNGA-style Integration` のPlotのダウンロード等のボタンをプロットの上部に移動する。
- [x] **3.** `CDR3 Physicochemical Properties` で "Compute Properties" ボタンを押すとエラーになる問題の修正。
- [x] **4.** `Repertoire Similarity` で出るエラーを修正する。
- [x] **5.** `Track clonotype` モジュールの日本語のラベルを英語に修正し、Plot Heightを `numericInput` にする。
- [x] **6.** `BCR Lineage Tree` でCloneが選択できない問題の修正。
- [x] **7.** `BCR SHM & Isotype` が何も表示されない問題の修正。

---

---

## 📄 論文執筆状況 / Paper Writing Status (paper_FiI.md)

### 投稿先
**Frontiers in Immunology — Technology and Code** （第一候補）

### 本文の構成と完成度

| セクション | 状態 | 備考 |
|-----------|------|------|
| Abstract | ✅ 完成 | CoNGA主張は穏当な記述に留めている |
| Introduction | ✅ 完成 | |
| Methods (全セクション) | ✅ 完成 | CoNGA説明をPython→純R実装に修正済み |
| Results: Data Processing (Fig 2) | ✅ 文章完成 | 図は未作成 |
| Results: TCR Repertoire (Fig 3) | ✅ 文章完成 | 図は未作成 |
| Results: Clonotype + DEG (Fig 4, 5) | ✅ 文章完成 | 図は未作成 |
| Results: Advanced Multimodal (Fig 6, 7) | ✅ 文章完成 | CoNGA+CellChatを1セクションに統合済み（2026-03-30） |
| Results: B-Cell Lineage (Fig 8) | ✅ 文章完成 | 図は未作成 |
| Results: CDR3 Analysis (Fig 9) | ✅ 文章完成 | 図は未作成 |
| Discussion | ✅ 完成 | CoNGA主張を「全PBMC解析の限界」として正直に修正済み |
| Limitations | ✅ 完成 | CoNGA limitation更新済み |
| References | ✅ 完成 | |

### 本日（2026-03-30）実施した主な変更

1. **CoNGA解析のメモリ問題を修正**: `dist()`→`RANN::nn2()`に変更し、50k細胞でのOOM killを解消
2. **CoNGA可視化を2パネルに改善**: 左=スコアgradient / 右=閾値以上の細胞をcluster色でhighlight
3. **論文のCoNGAセクションを現実に合わせて修正**:
   - 「Cluster 7への集中」→削除（データが支持しないため）
   - CoNGA独立セクション→CellChatと統合して「Advanced Multimodal Analyses」セクションに
   - Methods: Pythonバックエンド記述→純R実装に修正
   - Discussion: 「receptor-constrained reprogrammingの明確な証拠」→「全PBMC解析の感度限界」に修正
   - Limitations: CoNGA Python依存→T細胞サブセット推奨に更新

### ⏭️ 次にやること（論文完成に向けた残タスク）

#### 優先度 高
1. **全Figureの実際の図を作成・保存** （最重要）
   - Fig 2: UMAP (clusters, TCR+ cells, feature plots, conditions) → Dimensional Reductionモジュール
   - Fig 3: Diversity, Frequency, Alluvial → Repertoireモジュール
   - Fig 4: Clonotype1 highlight UMAP + Lasso screenshot → Dimensional Reductionモジュール
   - Fig 5: Volcano plot + Dot plot → DEGモジュール
   - Fig 6: CoNGA 2パネルUMAP + CellChat circle/pathway/bubble → **ビルド後に取得（2026-04-01 コード完成）**
   - Fig 8: BCR lineage tree + SHM violin + isotype → BCR Lineageモジュール
   - Fig 9: CDR3 sequence similarity network + heatmap + hydrophobicity → CDR3モジュール

#### ✅ CellChatモジュール 2条件比較機能 実装完了（2026-04-01）

**実装内容**: `app/module/6_cellCommunication.R` を全面改修し、2条件比較モードを追加。

**新機能**:
- **Analysis Mode**: Single Condition / Compare Two Conditions 切り替え
- **比較モード**: Condition Column・Condition A・Condition B をUI上で選択
- **比較用プロットタイプ**:
  - `Pathway Strength Comparison` → `rankNet(mode="comparison", stacked=TRUE)` でggplot出力
  - `Differential Interaction Network` → `netVisual_diffInteraction(comparison=c(1,2))`
  - `Circle Plot (per condition)` → 条件別に`netVisual_aggregate`
  - `Bubble Plot (per condition)` → 条件別に`netVisual_bubble`（sources/targets選択可）
- **内部ヘルパー関数** `.run_single_cellchat()` を切り出し、単条件・比較モード両方で再利用
- **liftCellChat + mergeCellChat** による正確な条件間比較（クラスター空間の統一）
- ダウンサンプリング: 単条件10000細胞 / 比較モード各条件7500細胞

**コンテナ内テスト結果（2026-04-01 確認済み）**:
- 条件1(sample=1): 19 pathways, 条件2(sample=2): 21 pathways
- liftCellChat → mergeCellChat → rankNet → netVisual_diffInteraction 全て成功
- テストスクリプト: `/tmp/cellchat_compare_test.R`（コンテナ内）

**✅ ビルド済み・アプリで動作確認済み（2026-04-01）**

**サンプルマッピング（確定）**:
- `sample=1` → HCW_baseline（バーコードサフィックス -1）
- `sample=2` → HCW_convalescent（バーコードサフィックス -2）

**クラスターラベルについての重要な制約**:
- CellChatは `"0"` というラベルを仕様上拒否する（`Cell labels cannot contain '0'!`）
- CellChatは `"0"` というラベルを仕様上拒否する（`Cell labels cannot contain '0'!`）
- そのためseurat_clustersの数字に `Cluster ` プレフィックスを付与: Cluster 0, Cluster 1, Cluster 7... 
- ~~旧仕様: C0, C7 → 2026-04-08 に「Cluster 7」形式に変更済み~~

---

### 📋 2026-04-08 セッション記録

#### ✅ 完了した修正

1. **CellChat クラスター名の変更** (`6_cellCommunication.R`)
   - `C0`, `C7` → `Cluster 0`, `Cluster 7` に変更し、可読性を向上
   - サイドバーの `※C7 = Seurat cluster 7` 説明文を削除

2. **CellChat ダウンロードボタンの配置移動** (`6_cellCommunication.R`)
   - サイドバー最下部 → メインパネル（図の上）に移動

3. **BCR Lineage Tree: dowser依存の完全除去** (`4_11_phylogeneticTree.R`)
   - **根本原因**: `dowser::formatClones` は Immcantation/Change-O パイプラインで事前処理された IMGT-gapped aligned sequences を要求する。10x Genomics の Cell Ranger 出力にはこのフォーマットが存在しないため、ダミー列を追加しても解決不可能だった。
   - **解決策**: `dowser::formatClones` + `dowser::getTrees` を完全にバイパスし、以下のパッケージで直接系統樹を構築:
     - `msa::msa()` (ClustalOmega) → マルチプルアラインメント
     - `ape::dist.dna()` → 距離行列
     - `ape::nj()` → Neighbor-Joining 法
     - `phangorn::midpoint()` → Midpoint rooting
     - `ggtree` → 描画

4. **ClustalOmega 書き込み権限エラーの修正** (`4_11_phylogeneticTree.R`)
   - **エラー**: `FATAL: Sorry, I do not have permission to write to file 'tempClustalOmega.aln'.`
   - **原因**: ClustalOmega が一時ファイルをカレントディレクトリ（`/srv/shiny-server/`）に書こうとしたが、`shiny` ユーザーに書き込み権限なし
   - **修正**: `setwd(tempdir())` で書き込み可能ディレクトリに移動してから実行。失敗時は Muscle にフォールバック。

5. **S4 クラス `$` エラーの修正** (`4_11_phylogeneticTree.R`)
   - treedata (S4) オブジェクトに `$tip.label` でアクセスしようとして失敗するフォールバック処理を `@phylo$tip.label` に修正

6. **メタデータ結合方式の変更** (`4_11_phylogeneticTree.R`)
   - `treeio::full_join(treedata, ...)` → `ggtree::`%<+%`(tree, meta_df)` に変更（ggtree 標準方式）
   - `phylogenetic_tree()` は `list(tree = phylo, meta = data.frame)` を返すように変更

#### ⏳ 未完了・次回やること

1. **BCR Lineage Tree の `%<+%` 関数形式修正の動作確認**
   - `ggtree::`%<+%`()` の関数形式呼び出しが正しく動作するか、コンテナ再ビルド後に確認が必要
   - **手順**: `docker-compose up --build -d` → ブラウザリロード → BCR Lineage Tree → Build Tree
   - **期待結果**: Color Tips By (sample) と Label Tips By (BCR_pair_CTaa) が正しく反映される

2. **Figure 8 の作成** (BCR Lineage Tree が動作確認できたら)
   - (A) 系統樹を sample で色分けして .pptx ダウンロード
   - (B) SHM & Isotype モジュールで IgG vs IgM の mutation rate を描画
   - (C) Isotype distribution を描画

3. **Figure 6 の作成** (CellChat — 「Cluster 7」名前変更反映の確認含む)
   - CellChat Compare モード → Pathway Strength Comparison, Circle Plot, Bubble Plot

4. **残りの Figure 作成** (Fig 2, 3, 4, 5, 9)

5. **paper_FiI.md の微修正**
   - CellChat の記述で「C7」→「Cluster 7」に更新が必要な箇所がないか確認

**CXCLパスウェイが検出されない問題（2026-04-01 調査完了）**:
- CXCL13はCluster 7で4/1205細胞（0.3%）しか発現なし → CellChatのthreshold未達
- CXCLはこのデータセットでは検出不可能（downsampling関係なく発現量が不足）
- 代わりに **IFN-II（IFNG）** がC7→全クラスターへ明確に検出されている
- **paper_FiI.md のCXCL記述を全てIFN-IIに修正済み（2026-04-01）**

**paper_FiI.md の変更点（2026-04-01 更新済み）**:
- CellChatの記述を「2条件比較モード（liftCellChat/mergeCellChat使用）」として正確に記述
- Figure 6(B): CXCL circle plot → **IFN-II circle plot** に変更
- Figure 6(C): "CXCL and IFN-II" → **"IFN-II and MIF"** に変更
- Discussion の CXCL/IFN-II 記述 → IFN-II and MIF に変更

**次回やること（Figure 6 キャプチャ）**:
1. アプリを起動（コンテナ名: `singlerepexplorer`、`localhost:3838`）
2. セッションファイルをアップロード: `singlerepexplorer_session_20260326_051713.qs2`
3. Cell-Cell Communication タブ → Compare Two Conditions を選択
4. 設定:
   - Cell Annotation: `seurat_clusters`
   - Condition Column: `sample` / Condition A: `1` / Condition B: `2`
   - Signaling Database: Secreted Signaling
   - **Run CellChat Inference ボタンを押す**（約5〜10分かかる）
5. **Figure 6(C)** をキャプチャ: Plot Type → `Pathway Strength Comparison` → Download (.pptx)
6. **Figure 6(B)** をキャプチャ: Plot Type → `Circle Plot (per condition)` / Show Condition: Condition B / Pathway: **IFN-II** → Download
7. **Figure 6(D)** をキャプチャ: Plot Type → `Bubble Plot (per condition)` / Show Condition: Condition B / Source Clusters: `C7` / Target: 空欄（全クラスター） → Download

2. **Abstractの最終確認**
   - CoNGAに関する記述が現在の結果と矛盾していないか確認

3. **Word count確認**
   - Frontiers Technology and Code: 通常 ≤ 12,000 words

#### 優先度 中
4. **Figure 1（ワークフロー図）の作成**
   - 現在はプレースホルダー。PowerPointかFigmaで概念図を作成

5. **各Figure captionと本文の対応確認**
   - 特にFig 6（CoNGA）のキャプションが実際の図と一致するか確認

#### 優先度 低
6. **最終校正（英語表現）**
   - ネイティブチェックまたはAI校正

---

### ⏳ 残りの課題 (Remaining Tasks)

0. **✅ CoNGA Multimodal Integration — 2パネルUMAP可視化 (2026-03-30 完了)**
   - **背景**: 現在のCoNGAモジュールはスコアのgradient表示のみで、どのクラスターに高スコア細胞が集中しているかが一目でわからず、論文Figure 6として不十分。
   - **実装内容**:
     - 左パネル: `FeaturePlot` CoNGAスコア gradient（既存）
     - 右パネル: CoNGAスコア ≥ 閾値N の細胞のみクラスター色で強調、それ以外はgrey（"Receptor-driven cells"）
     - サイドバーに閾値スライダー（`sliderInput`）を追加
     - `patchwork` で2パネル結合してダウンロードもpatchwork対応
   - **目的**: 論文 `paper_FiI.md` の Figure 6 として使用できる publication-quality な図を生成する

1. **チュートリアル動画の作成**
   - 新機能を含めた操作ガイド動画を作成し、GitHub Pages にリンクを貼る。

2. **多言語対応 (Optional)**
   - UIの日本語化オプションの検討（現在は英語のみ）。

3. **全Plotのカスタマイズ機能統一** → [詳細仕様は下記セクション参照](#plot-customization-spec)

4. **TCR/BCR ファイルアップロードの拡張子フィルター (.csv 限定)**
   - 現在、.h5 アップロードは `.h5` ファイルのみ表示されるが、TCR/BCR のファイル選択ダイアログではすべてのファイルが表示されてしまう。
   - `fileInput` の `accept` パラメータを `.csv` に限定し、ユーザーが誤って別形式のファイルを選択するのを防ぐ。

5. **解析データの qs2 ダウンロード機能**
   - アプリの最上部（常に表示される位置）に、現在の解析セッションデータ（Seurat オブジェクト等）を `.qs2` 形式でダウンロードできるボタンを設置する。
   - ダウンロードには時間がかかるため、プログレスバーまたはスピナー等でユーザーに進捗を明示する。

6. **qs2 データの再アップロード機能**
   - 上記でダウンロードした `.qs2` ファイルを再度アップロードして解析を再開できる機能を、Data Upload セクションに追加する。
   - アップロード後、保存された Seurat オブジェクトやメタデータを復元し、各解析モジュールが即座に使用可能な状態にする。

---

## Plot カスタマイズ機能 統一実装仕様 / Plot Customization Spec {#plot-customization-spec}

> **目的**: 非プログラマーの研究者が論文にそのまま使える図を作れるよう、全プロットでサイズ・スタイル設定を統一する。
> **方針**: 既存のコントロールは保持し、不足しているモジュールに追加実装する。
> **実装順**: 番号順に実装予定（✅ = 完了、⏳ = 未実装）

---

### 【共通コントロール】全モジュールに揃えるべき標準セット

| コントロール | UI部品 | デフォルト値 |
|---|---|---|
| Plot Width (px) | `numericInput` | 700 |
| Plot Height (px) | `numericInput` | 500 |
| Legend Position | `selectInput` (right/left/top/bottom/none) | right |
| Base Font Size | `numericInput` | 12 |

---

### 各モジュール別 追加実装仕様

---

#### 1. `1_uploadCellranger.R` — QC Violin Plots ✅（実装済み）

現状: qc_plot_width/qc_plot_height/qc_pt_size 実装済み（QC VlnPlot・ScatterPlot両方に適用）

---

#### 2. `2_reductionPlot.R` — UMAP / t-SNE / PCA ✅（基本実装済み）

現状: width/height/point_size/legend 実装済み

| 追加コントロール | 理由 |
|---|---|
| Label Size | cluster ラベルのフォントサイズ |
| Label Toggle (Show/Hide) | ラベル表示/非表示の切り替え |
| Axis Title Toggle | 軸タイトルの表示/非表示 |

---

#### 3. `3_1_geneExpression.R` — Feature / Violin / Dot / Heatmap ✅（基本実装済み）

現状: width/height/point_size/legend 実装済み。Plot typeごとに追加コントロールが必要。

| Plot Type | 追加コントロール | 理由 |
|---|---|---|
| **Feature Plot** | Color Scale (viridis/RdBu/etc.) | 発現量の色調整 |
| **Violin Plot** | Add jitter points (on/off)、jitter point size | 個々の細胞の表示 |
| **Dot Plot** | Dot Scale factor、X-axis label angle (0/45/90°)、Color palette | Dotサイズの最大値制御・ラベル可読性 |
| **Heatmap** | Gene label font size、Cell label show/hide | 遺伝子名の視認性 |

---

#### 4. `3_2_differentialGeneExpression.R` — Volcano Plot / DEG Heatmap ✅（実装済み）

現状: width/height/点サイズ/legend_pos/axis_font_size/gene_label_size 実装済み

---

#### 5. `3_3_phenotypeRepertoire.R` — Scatter Plots (Clonesize vs Expression) ✅（実装済み）

現状: width/height/point_size/base_font_size/facet_font_size/legend_pos 実装済み

---

#### 6. `3_4_congaIntegration.R` — Feature Plot / Violin Plot ✅（実装済み）

現状: commonPlotOptions(width/height/legend)/point_size 実装済み。downloadハンドラーも対応。

---

#### 7. `4_3_uniqueClones.R` — Bar Chart ✅（実装済み）

現状: width/height/legend/bar_border_color/x_axis_angle 実装済み

---

#### 8. `4_4_frequencyPlot.R` — Bar / Bubble / Heatmap ✅（実装済み）

現状: width/height/legend/x_axis_angle/base_font_size 実装済み

---

#### 9. `4_5_clonalLength.R` — Boxplot / Histogram ✅（実装済み）

現状: width/height/base_font_size/legend_pos/x_axis_angle 実装済み

---

#### 10. `4_6_diversity.R` — Bar / Box Plot ✅（実装済み）

現状: width/height/legend/x_axis_angle/jitter(show_jitter+jitter_size) 実装済み

---

#### 11. `4_7_trackClonotype.R` — Bar / Alluvial Plot ✅（実装済み）

現状: width/height/legend 実装済み

追加不要。

---

#### 12. `4_8_clonalOverlap.R` — Heatmap / Alluvial ✅（実装済み）

現状: width/height/legend/heatmap_palette/heatmap_text_size 実装済み

---

#### 13. `4_9_clonalProportion.R` — Stacked Bar Chart ✅（実装済み）

現状: width/height/legend_position/x_axis_angle 実装済み

---

#### 14. `4_11_phylogeneticTree.R` — Phylogenetic Tree (ggtree) ✅（実装済み）

現状: width/height/tip_size/tree_line_width/legend_pos/label_size/label_offset 実装済み

---

#### 15. `4_12_clonalAbundance.R` — Line / Rank Plot ✅（実装済み）

現状: width/height/legend/line_width/point_size 実装済み

---

#### 16. `4_13_publicClonotype.R` — UpSet / Venn Diagram ✅（実装済み）

現状: width/height/set_name_size/text_size 実装済み

---

#### 17. `4_14_shmIsotype.R` — Boxplot / Violin / Bar (SHM & Isotype) ✅（実装済み）

現状: width/height/base_font_size/legend_pos/x_axis_angle/jitter 実装済み

---

#### 18. `4_15_sequenceSimilarity.R` — Network Graph / Distance Heatmap ✅（実装済み）

現状: width/height/node_size_min・max/node_label_size/edge_alpha/heatmap_text_size/legend_pos 実装済み

---

#### 19. `4_16_cdr3Properties.R` — Violin / Box Plot ✅（実装済み）

現状: width/height/legend/x_axis_angle/jitter(show_jitter+jitter_size) 実装済み

---

#### 20. `5_trajectoryInference.R` — Trajectory / Pseudotime Plot ✅（実装済み）

現状: width/height/point_size/curveWidth/legend_pos/color_palette 実装済み

---

#### 21. `6_cellCommunication.R` — Network Circle / Bubble Plot ✅（実装済み）

現状: width/height/vertex_size_max/edge_width_max/label_font_size 実装済み

---

### 実装優先順位

> すべてのHigh/Mediumタスクは完了。Lowタスクのみ残存。

| 優先度 | モジュール群 | 理由 |
|---|---|---|
| ~~**High**~~ | ~~`3_4`, `4_11`, `1_upload`~~ | ✅ 完了 |
| ~~**Medium**~~ | ~~`3_2`, `3_1`(Dot/Heatmap), `5_trajectory`~~ | ✅ 完了 |
| ~~**Low**~~ | ~~`4_3`, `4_4`, `4_6`, `4_8`, `4_9`, `4_12`, `4_16`~~ | ✅ 完了 |

---

## 主要リンク / Links

- GitHub Pages: https://amufaamo.github.io/SingleRepExplorer/
- Docker Hub: `amufaamo/singlerepexplorer:latest` (v1.0.0)
- マニュアル: https://amufaamo.github.io/SingleRepExplorer/manual.html

---

## PowerPoint (Editable) Plot Download 実装計画 {#pptx-download-plan}

> **目的**: 全てのプロットダウンロードを PDF → PowerPoint (.pptx, Editable/Ungroupable) に変更する。
> **方針**: `officer` + `rvg` パッケージを使い、DrawingML形式で出力。PowerPointで「Ungroup」すると各要素を個別に編集可能にする。
> **ボタン統一**: 全モジュールのダウンロードボタンをプロットの**上**に統一配置。ボタンラベルは `.pptx` 表記に統一。

---

### 実装方針

1. **`utils.R` にヘルパー関数 `save_plot_as_pptx()` を追加**
   - ggplot オブジェクトを受け取り、officer + rvg で `.pptx` に保存する共通関数
   - CellChat 等の base R plot 用に `save_baseplot_as_pptx()` も追加
2. **`Dockerfile` に `officer`, `rvg` を追加**
3. **`app.R` に `library(officer); library(rvg)` を追加**
4. **各モジュールの `downloadButton` / `downloadHandler` を順次変更**
   - `.pdf` → `.pptx` に変更
   - `ggsave()` / `pdf()` → `save_plot_as_pptx()` / `save_baseplot_as_pptx()` に変更

---

### 各モジュール実装チェックリスト

| # | ファイル | ダウンロード対象 | 状態 |
|---|---|---|---|
| 1 | `utils.R` | ヘルパー関数追加 | ✅ |
| 2 | `Dockerfile` | `officer`, `rvg` 追加 | ✅ |
| 3 | `app.R` | `library()` 追加 | ✅ |
| 4 | `2_reductionPlot.R` | DimPlot | ✅ |
| 5 | `3_1_geneExpression.R` | Feature/Violin/Dot/Heatmap | ✅ |
| 6 | `3_2_differentialGeneExpression.R` | Volcano/DotPlot/Enrichment | ✅ |
| 7 | `3_3_phenotypeRepertoire.R` | Scatter Plot | ✅ |
| 8 | `3_4_congaIntegration.R` | UMAP + Violin | ✅ |
| 9 | `4_3_uniqueClones.R` | Bar Chart | ✅ |
| 10 | `4_4_frequencyPlot.R` | Bar/Heatmap (2箇所) | ✅ |
| 11 | `4_5_clonalLength.R` | Box/Histogram | ✅ |
| 12 | `4_6_diversity.R` | Bar/Box | ✅ |
| 13 | `4_7_trackClonotype.R` | Bar/Alluvial | ✅ |
| 14 | `4_8_clonalOverlap.R` | Heatmap + Alluvial | ✅ |
| 15 | `4_9_clonalProportion.R` | Stacked Bar | ✅ |
| 16 | `4_11_phylogeneticTree.R` | Phylogenetic Tree | ✅ |
| 17 | `4_12_clonalAbundance.R` | Rank Plot | ✅ |
| 18 | `4_13_publicClonotype.R` | UpSet/Venn | ✅ |
| 19 | `4_14_shmIsotype.R` | SHM + Isotype Plot | ✅ |
| 20 | `4_15_sequenceSimilarity.R` | Network + Heatmap | ✅ |
| 21 | `4_16_cdr3Properties.R` | Violin/Box | ✅ |
| 22 | `5_trajectoryInference.R` | Trajectory Plot | ✅ |
| 23 | `6_cellCommunication.R` | Network/Bubble | ✅ |

---

## Excel (.xlsx) Dataframe Download Implementation Plan {#xlsx-download-plan}

> **目的**: データフレームのダウンロードを CSV から Excel (.xlsx) に変更する。
> **方針**: `openxlsx` パッケージを使用。
> **ボタン統一**: ダウンロードボタンをテーブルの**上**に統一配置。

### 各モジュール実装チェックリスト

| # | ファイル | ダウンロード対象 | 状態 |
|---|---|---|---|
| 1 | `3.2_differentialGeneExpression.R` | marker/DEG tables | ✅ |
| 2 | `3.3_phenotypeRepertoire.R` | correlation table | ✅ |
| 3 | `4.1_clonotypeInformation.R` | info table | ✅ |
| 4 | `4.3_uniqueClones.R` | unique clones table | ✅ |
| 5 | `4.4_frequencyPlot.R` | frequency table | ✅ |
| 6 | `4.5_clonalLength.R` | length table | ✅ |
| 7 | `4.6_diversity.R` | diversity table (Added) | ✅ |
| 8 | `4.7_trackClonotype.R` | tracking table | ✅ |
| 9 | `4.8_clonalOverlap.R` | overlap matrix (Added) | ✅ |
| 10 | `4.9_clonalProportion.R` | proportion table (Added) | ✅ |
| 11 | `4.10_antigenPrediction.R` | prediction table | ✅ |
| 12 | `4.12_clonalAbundance.R` | abundance table | ✅ |
| 13 | `5_trajectoryInference.R` | pseudotime table | ✅ |
| 14 | `7_llmAnnotation.R` | annotation table (Added) | ✅ |
