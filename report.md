# **シングルセルRNA-seqと免疫レパトアの統合解析における網羅的アプローチと先端手法**

## **1\. 導入**

適応免疫系の驚異的な特異性と多様性は、T細胞受容体（TCR）およびB細胞受容体（BCR）の遺伝子再構成プロセスであるV(D)J組換え、ならびに抗原曝露に伴うその後のクローン選択メカニズムによってもたらされる 1。歴史的に、免疫レパトアの解析はバルクRNAシーケンスやバルクDNAシーケンスに依存してきた。これらのバルク解析手法は、組織や血液サンプル全体のレパトア多様性を巨視的に捉えることには成功していたものの、個々の細胞内における重鎖・軽鎖（あるいはα鎖・β鎖）の正確なペアリング情報を取得することや、受容体配列の特性と細胞の転写状態（分化段階、活性化状態、疲弊化など）との直接的な連関を単一細胞レベルで解明することは技術的に不可能であった 2。

しかしながら、近年のシングルセル5'末端RNAシーケンス（scRNA-seq）技術、特に10x Genomics社などのプラットフォームの飛躍的な発展により、単一細胞レベルでの完全長に近いTCR/BCR配列（scVDJ-seq）の取得と、全転写産物の包括的なプロファイリングを同時に実行することが可能となった 1。この技術的ブレイクスルーは、腫瘍浸潤リンパ球（TIL）の微小環境における動態解析、自己免疫疾患における病原性クローンの特定、そして感染症（COVID-19など）に対する中和抗体の探索など、幅広い医学生物学分野にパラダイムシフトをもたらしている 5。例えば、Perezらはシングルセルシーケンスを大規模コホートに適用するためのコスト効率の高い高スループット手法「mux-seq」を開発し、全身性エリテマトーデス（SLE）患者と健常対照者から得られた120万個の末梢血単核球（PBMC）を解析することで、疾患メカニズムに関する貴重な洞察を得ている 4。

本報告書では、シングルセルRNA-seqと免疫レパトア（TCR/BCR）データの統合解析において実施される各種解析手法を可能な限り網羅的に体系化する。シーケンスデータの前処理から、多様性指標の数学的評価、体細胞高頻度突然変異（SHM）に基づく系統樹推論、細胞軌道へのクローン情報のマッピング、多層オミクス（マルチモーダル）の統合、そして最新の大規模言語モデル（LLM）や深層学習アプローチに至るまで、その理論的背景、生物学的意義、および実践的適用について詳述する。

## **2\. 解析パイプラインの基盤とデータ前処理**

シングルセルレパトア解析の第一段階は、次世代シーケンサーから出力された生の配列データを生物学的に意味のある構造化データへと変換するプロセスである。このプロセスは、リードのアセンブリ、アノテーション、そして厳密な品質管理（QC）から構成される。

### **2.1 配列のアセンブリとアノテーション**

取得されたFASTQファイルは、まずリファレンスゲノム（ヒトまたはマウスなど）に対してアライメントされ、シーケンスエラーの補正とコンティグ（Contig）の構築が行われる 8。この工程は一般的にCell Ranger V(D)Jパイプラインなどが担うが、組み立てられたコンティグデータは、その後さらに専門的なツール群を用いて詳細なアノテーションが付与される。ImmcantationスイートやMiXCRといったツール群は、V、D、Jの各遺伝子セグメントの対立遺伝子レベルでの割り当てを行い、同時に相補性決定領域（CDR1、CDR2、CDR3）およびフレームワーク領域（FRW）の境界を正確に決定する 9。

解析データの標準化という観点から、これらの前処理パイプラインから出力される中間データおよび最終データは、Adaptive Immune Receptor Repertoire（AIRR）フォーマットに準拠することが強く推奨されており、ScirpyやscRepertoireといった後続の解析パッケージはこのAIRRフォーマットをシームレスに読み込む機能を備えている 9。Scirpyの入出力モジュール（scirpy.io）は、10x Genomicsの出力だけでなく、TraCeR、BraCer、BD Rhapsodyなど、多様なプラットフォームからのデータを統合的にインポートできる柔軟性を持つ 11。

### **2.2 レパトア特有の品質管理とチェインのペアリング**

シングルセルデータ特有のノイズを排除するため、レパトアデータにはscRNA-seq側の指標（ミトコンドリア遺伝子含有率や検出遺伝子数など）とは独立した、免疫受容体特有の厳格な品質管理（QC）が適用される 12。

最初の関門は、プロダクティブ（機能的）な配列の選別である。V(D)J組換えの過程で生じるフレームシフト変異や、配列の途中に未成熟な終止コドンを含むコンティグは、機能的な受容体タンパク質として細胞表面に発現しないため、解析から除外される 9。Scirpyのchain\_qcモジュールなどは、受容体の生産性とペアリングの妥当性を体系的に評価する 11。

次に、ダブレット（細胞の二重取り）の特定と排除が行われる。シングルセル解析プラットフォームでは、1つの油滴（ドロップレット）またはウェルに2つ以上の細胞が封入される確率が一定割合で存在する。T細胞において2つ以上の異なるα鎖またはβ鎖、あるいはB細胞において複数の重鎖が高いリードカウントで検出された場合、それはダブレットである可能性が極めて高い。scRepertoireやScirpyは、このような異常なチェイン構成を持つ細胞を特定し、下流の解析からフィルタリングする強力な機能を備えている 10。

品質管理を通過したデータは、シングルセルID（10xバーコードなど）をキーとして、重鎖（IGH）と軽鎖（IGK/IGL）、またはTCRα鎖（TRA）とTCRβ鎖（TRB）が正確にペアリングされ、1つの統合された「クローン」として再構築される 3。この単一細胞解像度でのペアリングこそが、後述する抗原特異性の推論やクローン動態解析の基礎となる。

## **3\. 免疫レパトアの多様性とクローナリティの定量評価**

前処理が完了したレパトアデータを用いて、細胞集団におけるクローンの拡大状況やレパトア全体の多様性が数学的に定量化される。これは、個体の免疫状態（感染症への応答、腫瘍免疫の活性化、自己免疫の暴走など）を客観的に評価するための中心的な解析である 15。

### **3.1 クロノタイプの定義付け**

多様性を計算する前提として、「何をもって同一のクローン（クロノタイプ）とみなすか」という定義を明確にする必要がある。scRepertoireやDandelionなどのパッケージでは、研究の目的に応じてこの定義を柔軟に変更できる 16。 最も厳密な定義（Strict）では、V遺伝子、J遺伝子、およびCDR3領域の塩基配列がすべて完全に一致する細胞群のみを同一クロノタイプと見なす 16。一方で、B細胞のように体細胞高頻度突然変異（SHM）を蓄積する系譜を追跡する場合、塩基配列の完全一致を求めると同一系統のクローンが細分化されてしまう。そのため、同一のVおよびJ対立遺伝子を共有し、CDR3アミノ酸配列のハミング距離またはレーベンシュタイン距離が一定の閾値（例えば80%の類似性）以内であるものを、同一のクローンファミリーとしてグループ化するアプローチが採用される 9。

### **3.2 多様性と均等性の数学的指標**

免疫レパトアの多様性は、生態学における種多様性の概念を応用し、「豊富さ（Richness：存在する異なるクローンの総数）」と「均等性（Evenness：各クローンの出現頻度の偏りのなさ）」の2つの軸から評価される。これらの指標は、細胞のサンプリング深度（アンダーサンプリング）によるバイアスを受けやすいため、ダウンサンプリングやブートストラップ法を併用して平均値を算出するアプローチが推奨される 16。

以下の表は、シングルセルレパトア解析で適用される主要な多様性指標の特性をまとめたものである。

| 指標カテゴリ | 指標名 | 数学的定義および特徴 | 生物学的解釈と適用場面 | 参照 |
| :---- | :---- | :---- | :---- | :---- |
| **アルファ多様性** | Shannon Entropy | ![][image1] (![][image2]: クローン![][image3]の頻度)。情報理論に基づく指標。 | 豊富さと均等性の両方を考慮する。少数の巨大クローンによる独占があると値が低下するが、アンダーサンプリングの影響を強く受けるため少数の細胞しか得られないサンプルには不向きである。 | 15 |
| **アルファ多様性** | Gini-Simpson Index | ![][image4] | レパトアからランダムに2つの細胞を抽出した際、それらが異なるクローンに属する確率。均等性に関する情報を強く反映し、シミュレーションデータにおいてもダウンサンプリングに対して高い頑健性（ロバスト性）を示す。 | 15 |
| **アルファ多様性** | Chao1 / ACE | 観測されたクローン数に、頻度1および2のクローン数に基づく補正を加える。 | 未発見の希少なクローン数を推定し、レパトア全体の「真の豊富さ」を算出する。多様性の下限値を推定するのに優れている。 | 19 |
| **アルファ多様性** | D50 / Gini Index | D50はレパトア全体の50%を占める上位クローンの割合。Gini Indexは不平等の尺度。 | クローン拡大（Clonal expansion）の程度を直感的に評価する。値が0に近いほど、少数のクローンがレパトアを独占している状態（強い免疫応答またはクローン性増殖）を示す。 | 21 |
| **アルファ多様性** | Pielou's Evenness | ![][image5] (![][image6]: 観測されたクローン総数)。 | 均等性に特化した指標。値が1に近いほど、すべてのクローンが均等なサイズで存在することを示し、特定抗原への強い選択圧がかかっていない状態を示唆する。 | 19 |
| **ベータ多様性** | Morisita-Horn Index | 2つのサンプル間で共有されるクローンの相対頻度を基に計算される類似度指数。 | 血液と組織間、あるいは治療前後におけるレパトアの重複を評価する。Jaccard指数が「クローンの有無」に依存するのに対し、本指標は存在量の多いクローンの影響を強く受けるため、主要な免疫応答のトラッキングに適している。 | 18 |

### **3.4 遺伝子使用頻度とCDR3配列の構造的特性**

抗原に対する応答は、特定のV遺伝子やJ遺伝子の優先的な使用（V-J Skewing）として現れることが多い。PlatypusパッケージのVDJ\_VJ\_usage\_circos関数などは、V遺伝子とJ遺伝子のペアリング頻度をCircosプロットとして視覚化し、特定の疾患状態や組織特異的な遺伝子使用の偏りを明らかにする 25。

さらに、scRepertoireの最新バージョンでは、CDR3アミノ酸配列に沿った物理的特性（電荷、疎水性、かさ高さなど）のマッピングや、各アミノ酸残基ごとの位置的エントロピー（Positional Entropy）を定量化するpositionalProperty()およびpositionalEntropy()関数が実装されている。これにより、特定の抗原エピトープの認識に関与する可能性が高い可変モチーフや、構造的安定性に寄与する保存領域をピンポイントで特定することが可能となり、単なるクローンの集計を超えた機能的なレパトア解釈が実現している 26。

## **4\. B細胞特異的解析：親和性成熟、SHM、および系統樹推論**

T細胞のレパトア解析が主に「クローンサイズの動態」に焦点を当てるのに対し、B細胞受容体（BCR）の解析は、抗原感作後に胚中心（Germinal Center; GC）で生じる「抗体の進化」と「系統の分岐」の解明に重点が置かれる。これには、体細胞高頻度突然変異（SHM）とクラススイッチ再構成（CSR）の精密なトラッキングが不可欠である 3。

### **4.1 SHMとアイソタイプ・クラススイッチの評価**

抗原に遭遇し胚中心反応に移行したB細胞は、濾胞性ヘルパーT細胞（Tfh）からのシグナル（CD40L、IL-21など）を受け、活性化誘導シチジンデアミナーゼ（AID）の作用によりCDR領域を中心に高頻度で点突然変異を蓄積する 9。同時に、IgMやIgDからIgG、IgA、IgEといった機能的なアイソタイプへのクラススイッチ（CSR）を経験する 3。 シングルセル解析では、このSHM率とアイソタイプの情報が各細胞のバーコードと正確に紐付けられる。例えば、急性骨髄性白血病（AML）の再発・難治性モデルの解析では、再発時のBCRレパトアにおいて、特定の巨大クローンへの極端な偏り（多様性の低下）とともに、IgGおよびIgAへのクラススイッチの優位性、そして顕著に高いSHM頻度が観察されている。これは、再発に至る過程でB細胞が抗原駆動型の強い選択圧を受け、親和性成熟を遂げたことを示唆しており、疾患の進行とB細胞の動態が深くリンクしていることを裏付けている 13。

### **4.2 系統樹推論と祖先配列の復元**

同一のクローンファミリーに属するBCR配列群がどのように進化してきたかを可視化するため、系統樹（Phylogenetic tree）の構築が行われる。PlatypusのVDJ\_tree関数や、より高度に専門化されたSeQuoIA（Selection Quantification in Integrative AIRR data）パイプラインは、この系統樹推論を自動化する強力なツールである 9。

SeQuoIAでは、重鎖と軽鎖のペア情報を用いて最大節約法（Maximum Parsimony）に基づく系統樹（Parsimony Forest）を推論する 9。この過程で、配列変異の起点となる「祖先配列」の復元が行われる。ナイーブB細胞に由来する初期のクローンファミリーに対しては、リファレンスから推論された生殖系列（Germline）配列を祖先として設定する。一方、すでに親和性成熟が高度に進行しているメモリーB細胞などのファミリーに対しては、ファミリー内で最も変異の少ない配列群の厳密なコンセンサスから「最近隣共通祖先（Nearest Common Ancestor: NCA）」を再構築し、これを系統樹のルートに設定する 9。

### **4.3 生物学的制約に基づく系統樹の最適化と選択圧の推論**

数学的に構築された複数の系統樹候補の中には、生物学的に起こり得ない進化経路が含まれることがある。SeQuoIAは、B細胞の基礎生物学に基づいたフィルタリングを適用して最適な系統樹を選択する。具体的には、B細胞が生存を維持するためには受容体からの持続的なシグナル（Tonic signaling）が必須であるという事実に基づき、終止コドンやフレームシフト変異を含む「非機能的な内部ノード」を持つ系統樹を棄却する 9。さらに、独立したクラススイッチ（CSR）イベントが系統樹上で無数に発生するような非効率なツリーや、一度変異した塩基が元に戻る逆変異（Reversion）が多発するツリーもペナルティを与えられ除外される 9。

最適化された系統樹が構築されると、SeQuoIAはShazamパッケージなどを利用してAIDによるSHMの5-mer標的モデル（Mutability model）を構築し、系統樹上での実際の変異パターンと比較する。これにより、体細胞突然変異を「インビボの分子トレーサー」として活用し、胚中心の明領域（Light Zone）や暗領域（Dark Zone）におけるBCR駆動型の選択圧を定量化することが可能となる。これは、ワクチン接種後の効果判定や、B細胞リンパ腫の発生メカニズムを解明する上で極めて強力なアプローチである 9。

## **5\. T細胞受容体（TCR）の特異性、距離学習、およびパブリッククローン**

T細胞の解析においては、そのTCRがどのようなペプチド-MHC（pMHC）複合体を認識するのかという「抗原特異性」の推論が最大の課題となる。シングルセル解析によりα鎖とβ鎖のペア情報が得られることで、配列類似性に基づいた高度な距離学習モデルの適用が可能となった。

### **5.1 TCRdist3による距離学習と配列類似性ネットワーク**

同一の抗原エピトープを認識するTCRは、全く同じ塩基配列やアミノ酸配列を持っているわけではない。しかし、抗原との接触面に位置する重要なアミノ酸残基（モチーフ）を共有していることが多い 27。この潜在的な類似性を数学的に評価し、抗原特異的な細胞群を抽出するための代表的なアルゴリズムが**tcrdist3**である 27。

Dashらによって開発されたTCRdistアルゴリズムは、α鎖とβ鎖の両方を考慮した多重CDR距離（Multi-CDR distance）メトリクスを採用している。アミノ酸の置換ペナルティにはBLOSUM62マトリックス（0から-4のスコア）を使用し、ギャップには-4のペナルティを与える。特筆すべきは、抗原結合に対して最も直接的かつ支配的な影響を与えるCDR3領域でのアミノ酸置換に対して、他の領域の3倍の重み付け（3-fold weighting）を行う点である 28。 この距離マトリックスを用いることで、PlatypusやScirpyといった解析パッケージは、細胞間の配列類似性ネットワーク（Sequence Similarity Network）を構築する。一定の閾値以下の距離（例えば数個のアミノ酸変異のみの違い）にあるクローン同士をエッジで結びつけることで、バルクレパトアの中から抗原選択圧を受けているポリクローナルな受容体群を同定し、背景ノイズを補正したCDR3配列ロゴ（Logo plot）として可視化することが可能となる 11。

### **5.2 収斂組換えとパブリッククローンの動態**

免疫レパトア解析において、遺伝的に無関係な複数の個体間で共通して検出される同一のクロノタイプは「パブリッククローン（Public clones）」と呼ばれる 29。これらが共有される主要な生物学的メカニズムの一つが\*\*収斂組換え（Convergent Recombination）\*\*である。V(D)J組換えプロセスには、生殖系列遺伝子の使用頻度における固有のバイアスが存在し、異なるヌクレオチド配列が翻訳の結果として同一のアミノ酸配列（特にCDR3）をコードする現象が頻発する 30。

パブリッククローンと収斂組換えの解析は、特定の疾患や発達段階の免疫状態を理解する上で重要な鍵となる。例えば、超未熟児（妊娠23〜27週）のナイーブCD8+ T細胞レパトアの解析では、正期産児や成人と比較してPielou's EvennessやiChao1などの多様性指標が著しく低く、無関係な個人間で共有されるパブリッククローンの頻度が異常に高いことが示されている。これは、未熟児における胸腺でのT細胞選択が成人に比べて厳密でなく、特定の収斂した配列を持つT細胞が末梢に流出していることを示唆しており、彼らの感染症に対する高い脆弱性を説明する一因となっている 20。 一方で、特定のウイルス感染症に対する防御免疫の観点からもパブリッククローンは重要である。COVID-19回復期患者のB細胞レパトア解析では、SARS-CoV-2のスパイクタンパク質受容体結合ドメイン（RBD）に対して特異性を持つ、高度に収斂したメモリーB細胞クローン（IgG）が健常者モデルの予測をはるかに超える頻度で同定されている。これは、ウイルス抗原による強い「収斂選択（Convergent selection）」の結果であり、ワクチン設計や普遍的な中和抗体療法の開発に向けた強力なターゲットとなる 3。

### **5.3 既知抗原データベースとの照合**

シングルセルデータから得られたTCR/BCR配列を、VDJDBやIEDB（Immune Epitope Database）などのキュレーションされた抗原特異的受容体データベースに対してクエリ検索する手法も、Scirpyなどのパイプラインに標準実装されている。これにより、サンプル内に存在するCMV、インフルエンザ、SARS-CoV-2などのウイルス特異的クローンや、既知の腫瘍抗原特異的クローンの割合を定量化し、細胞の表現型と関連付けることができる 11。

## **6\. 転写産物データ（scRNA-seq）との統合解析手法**

シングルセルRNA-seqとTCR/BCR-seqを統合する最大の科学的価値は、受容体の特異性（クロノタイプ）と、その細胞が現在実行している転写プログラム（細胞の分化状態、機能、疲弊など）を細胞単位で直接的に結びつけることができる点にある 1。

### **6.1 クローンサイズの拡大と細胞表現型・疲弊化の相関**

免疫応答の最前線では、抗原を認識した特定のクローンが急速に増殖する（Clonal expansion）。この拡大したクローン（Expanded clones）と、単一細胞でのみ検出されるクローン（Singletons）の間で、遺伝子発現プロファイルを比較することで、免疫細胞の機能的二面性が明らかになる。

例えば、腎細胞がん（RCC）や肝細胞がん（HCC）における腫瘍浸潤リンパ球（TIL）の統合解析では、クローンサイズが拡大しているCD8+ T細胞集団において、グランザイムB（GZMB）、パーフォリン（PRF1）、NKG7、インターフェロンガンマ（IFNG）といった強い細胞傷害性マーカーの発現が認められると同時に、HAVCR2（TIM-3）、CD101、PD-1、ITGAE（CD103）といったT細胞疲弊（Exhaustion）および組織常在性マーカーの発現が著しく上昇していることが確認されている 7。対照的に、クローンサイズが小さいT細胞やシングルトンは、CCR7やSELLといったナイーブ／メモリー様マーカーを発現し、腫瘍内の異なる空間的密度領域に局在する傾向がある 35。 さらにHCCの解析では、抗PD-1抗体と抗VEGF抗体の併用療法に著効した患者において、拡大したCD8+ T細胞クローンやマクロファージ、樹状細胞に至るまで、広範な免疫細胞系譜においてインターフェロン（IFN）シグナル伝達経路への収束が観察されている。これは、レパトアの拡大情報と転写プログラムの統合が、免疫療法（ICI）の治療応答性を予測する信頼性の高いバイオマーカーの発見に直結することを示している 34。PlatypusやSeuratといったパッケージでは、これらのクローン情報をメタデータとして統合し、UMAP（Uniform Manifold Approximation and Projection）次元削減プロット上で拡大クローンの空間的分布を可視化したり、特定のクローンに属する細胞群を対象とした発現変動遺伝子（DEG）解析を自動化する機能が提供されている 12。

### **6.2 遺伝子発現空間とTCR/BCR空間のグラフベースの統合：CoNGA**

トランスクリプトームに基づくクラスタリングを先に行ってからTCRをマッピングする従来のアプローチでは、遺伝子発現の差異が微小である場合、TCR配列に依存した微妙な表現型の違いを見逃すリスクがある。この課題を解決するため、\*\*CoNGA（Clonotype Neighbor Graph Analysis）\*\*という画期的なアルゴリズムが開発された 38。

CoNGAの「Graph-vs-Graph」解析アプローチは、遺伝子発現プロファイルに基づく類似性グラフ（K近傍グラフ：KNNグラフ）と、TCR配列の類似性（tcrdistベース）に基づくグラフを独立して構築し、両者のグラフ間での統計的に有意な「細胞の重なり（Neighborhood overlap）」を体系的に探索する 38。これにより、事前のクラスタリングによるバイアスを排除した状態で、TCR配列と遺伝子発現が強く相関する細胞サブセットを教師なしで発見することが可能となる。実際、CoNGAを用いた解析により、これまで定義されていなかったHOBIT+/HELIOS+のヒト循環自然リンパ球様CD8+ T細胞集団や、胸腺細胞の分化状態を決定づけるTCR配列の決定要因など、バルク解析や単一モダリティ解析では決して見出せなかった新たな生物学的発見がもたらされている 32。

### **6.3 レパトア情報を組み込んだ細胞軌道と擬似時間（Pseudotime）解析**

細胞の分化や成熟の過程を連続的な軌道として推定する軌道推論（Trajectory inference）は、scRNA-seq解析の中核技術の一つである。しかし、転写産物のみに基づく軌道は、細胞が「同一のクローンから派生したかどうか」という系統的関係性（Clonal lineage）を全く反映していないという根本的な欠点を持つ 41。

この問題を解決するため、\*\*LRT（Lineage and RNA Trajectory）\*\*フレームワークや、**Dandelion**などのツールが開発され、scRNA-seqデータから構築された擬似時間に対してscTCR/BCR-seqのクローン情報を統合する機能を提供している 17。 LRTは、まず転写産物情報を用いて細胞の全体的な分化軌道を構築し、その後、TCR配列情報と細胞の表現型情報を利用して、特定の分化経路に対して強いバイアスを持つクロノタイプのクラスターを特定する 41。例えば、急性リンパ球性脈絡髄膜炎ウイルス（LCMV）感染モデルのCD4+ T細胞の解析において、LRTはscRNA-seq単独では不明瞭であった、Th1細胞と濾胞性ヘルパーT細胞（Tfh）への分岐軌道をクローン特異性の観点から明確に描き出した 41。これは、T細胞の運命決定（Fate decision）が単なる微小環境のシグナルだけでなく、TCRの抗原認識に基づく固有の特異性によって強く方向付けられていることを証明するものであり、軌道推論におけるレパトア統合の重要性を決定づけている。

## **7\. マルチモーダル統合と多層オミクス解析の最前線**

シングルセル解析技術は、RNAとV(D)J配列の同時取得にとどまらず、クロマチンアクセシビリティや細胞表面タンパク質など、さらに多くのモダリティ（多層オミクス）を同時にプロファイリングする方向へ進化している。

### **7.1 CITE-seqおよびscATAC-seqとの統合的アプローチ**

細胞表面タンパク質の絶対量をオリゴヌクレオチド標識抗体を用いて測定するCITE-seq（Cellular Indexing of Transcriptomes and Epitopes by Sequencing）や、オープンクロマチン領域をプロファイリングするscATAC-seqは、RNAの転写後制御やエピジェネティックな状態を解明するための強力なツールである 29。 LCMV感染モデルのCD8+ T細胞疲弊に関する研究では、scRNA-seq/scTCR-seqに加えてscATAC-seqが適用され、疲弊の進行（TexprogからTexint、そしてTextermへ）に伴う遺伝子座のアクセシビリティの変化と、特定のTCRクローンの拡大動態がエピジェネティックなレベルで関連付けられている 44。

データ解析の側面では、Seuratパッケージが提供する\*\*WNN（Weighted Nearest Neighbors）\*\*解析がこれらのマルチモーダルデータの統合において標準的な手法となっている 45。WNNは、異なるモダリティ間の情報量の重みを細胞ごとに最適化しながら単一の統合された隣接グラフを構築する。最新の応用例では、遺伝子発現データとCITE-seqのタンパク質データから第一のWNNグラフを構築し、さらにそのグラフをCoNGAやtcrdistに基づくTCR類似性グラフと結合する「Tri-modal embedding（3モダリティ統合）」というアプローチが報告されている。この極めて解像度の高い統合解析により、ワクチン投与によって誘導された真の「抗原特異的増殖CD8+ T細胞」と、サイトカインなどによって非特異的に活性化された「傍観者（Bystander）細胞」の細胞状態を、明確に区別し特徴づけることに成功している 45。

さらに、SAILERXのような最新のアルゴリズムは、scRNA-seqとscATAC-seqのマルチオームデータに対して、より頑健な遺伝子発現情報をリファレンスモダリティとして利用し、クロマチンアクセシビリティ空間の学習プロセスを正則化することで、モダリティ間の異質性を克服し、過学習を防ぐ統合メカニズムを提供している 47。

## **8\. 深層学習・機械学習と大規模言語モデル（LLM）の応用**

免疫レパトア解析の最大の手法的な壁は、TCRやBCRのアミノ酸配列（特にCDR3）が「可変長」である点にある。長さが異なる文字列データを、そのまま主成分分析（PCA）やクラスタリングなどの標準的な機械学習アルゴリズムに投入することは数学的に困難である。この問題と、細胞型アノテーションのボトルネックを解消するため、人工知能（AI）技術の導入が急速に進んでいる。

### **8.1 オートエンコーダによる受容体配列の潜在表現学習（Representation Learning）**

この可変長配列の問題を解決するために、scRepertoireのエコシステムに統合されている**Trex**（T細胞用）や**Ibex**（B細胞用）、さらにはそれらを支える**ImmApex**といったツール群が開発された 10。 これらは畳み込みニューラルネットワーク（CNN）ベースのオートエンコーダアーキテクチャを採用しており、可変長のCDR3アミノ酸配列を入力として受け取り、その物理化学的特性やモチーフのパターンを学習した上で、固定次元の「潜在ベクトル（Latent dimensions）」へと圧縮・変換する 26。このプロセスにより抽出された意味のある潜在次元は、SeuratやSingleCellExperimentなどの単一細胞解析ワークフローに数値行列として直接組み込むことが可能となる。この表現学習（Representation Learning）により、レパトア配列情報をscRNA-seqの連続的な数値データと同等に扱うことができるようになり、抗原特異性の予測モデルの構築や、マルチモーダルな細胞クラスタリングの精度が飛躍的に向上する 14。

### **8.2 大規模言語モデル（LLM）を用いた細胞型アノテーションの自動化**

シングルセル解析におけるもう一つの巨大な課題は、膨大な数の細胞クラスターに対して正確な「細胞型（Cell type）」の意味付け（アノテーション）を手作業で行うことの非効率性と主観性である。従来は、各クラスターで発現が変動しているマーカー遺伝子（DEG）のリストを文献と照らし合わせて手動で特定するか、特定の参照データセット（Reference dataset）に依存した自動ツールを使用するしかなかった。しかし、参照データセットに存在しない希少な細胞サブセットや、組織特異的な状態を持つ細胞の判定において、これらの手法は限界を迎えていた 50。

この課題に対する画期的な解決策として、**LICT（Large Language Model-based Identifier for Cell Types）やmLLMCelltype**といった、GPT-4、Claude 3、Geminiなどの大規模言語モデル（LLM）を活用した自動アノテーションフレームワークが実用化されている 50。 これらのツールは、静的な参照データに依存するのではなく、LLMのパラメータ内にエンコードされた膨大な生物医学的知識ベースを活用する。具体的には、各クラスターのDEGリストや細胞のコンテキスト情報をプロンプトとして複数のLLMに同時に入力し、出力されたアノテーション結果に対してモデル間でのクロスチェック（Multi-LLM consensus）を行う。mLLMCelltypeのようなフレームワークでは、モデル間の意見の不一致（不確実性）をシャノンエントロピーなどの指標を用いて定量化し、専門家によるマニュアルレビューが必要な曖昧なクラスターを明示する 51。 さらに、LLMの最大の強みである「Talk-to-machine」アプローチにより、「なぜその遺伝子発現パターンからその細胞型であると推論したのか」という論理的な推論過程（Reasoning chain）が透明性をもって出力される 50。このLLMベースの精密なアノテーションをTCR/BCRレパトア解析と組み合わせることで、「どの抗原特異的クローンが、どの精緻に定義された表現型サブセット（例：特定の疲弊段階にあるCD8+ T細胞や、特定のケモカイン受容体を発現するメモリーB細胞）に局在しているか」を、極めて客観的かつスケーラブルに評価することが可能となっている。

## **9\. 結論**

シングルセルRNA-seqと免疫レパトア解析（scTCR/BCR-seq）の統合的アプローチは、単なる細胞集団のカタログ化やクローン数のカウントという段階をとうに過ぎ、細胞の「状態（State）」と受容体の「特異性（Specificity）」の間の複雑な因果関係をシステムレベルで解明する次元へと到達している。

シーケンスデータの堅牢な前処理から始まり、クロノタイプの柔軟かつ厳密な定義、多様性・均等性の数理的評価を基盤とすることで、免疫応答の全体像を正確に把握することができる。さらに、B細胞における親和性成熟の軌跡を追う系統樹推論（SeQuoIA）や、T細胞における距離学習に基づく配列類似性ネットワーク（TCRdist3）の構築は、レパトアデータに隠された抗原選択圧のメカニズムを浮き彫りにする。そして、CoNGAやLRTといった先駆的なアルゴリズムを用いた転写産物とのグラフベース統合・軌道解析や、WNNを用いた多層オミクス（CITE-seq/scATAC-seq）のシームレスな統合は、従来のモダリティ単独の解析では決して到達し得なかった、高解像度での細胞状態の定義を可能にしている。

加えて、LLMによる細胞型アノテーションの自動化（LICT、mLLMCelltype）や、オートエンコーダを利用した受容体配列の潜在表現学習（Trex、Ibex）といった最先端の人工知能技術がパイプラインの中核に組み込まれたことで、解析のスループット、客観性、および予測モデリングの精度はかつてない水準に達している。これらの網羅的かつ統合的な解析手法の適用は、がん免疫療法（ICI）における治療応答性バイオマーカーの確立、自己免疫疾患における病原性メカニズムの解明、そして次世代のワクチンや抗体医薬の合理的設計において、今後ますます不可欠な基盤技術として機能し続けるであろう。

#### **引用文献**

1. Comprehensive Analysis of TCR and BCR Repertoires: Insights into Methodologies, Challenges, and Applications \- PMC, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC11853700/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11853700/)  
2. TCRscape: a single-cell multi-omic TCR profiling toolkit \- PMC \- NIH, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC12446293/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12446293/)  
3. BCR vs. TCR Sequencing: Key Differences and Applications \- CD Genomics, 3月 2, 2026にアクセス、 [https://www.cd-genomics.com/biomedical-ngs/resource/bcr-vs-tcr-sequencing-immune-repertoire.html](https://www.cd-genomics.com/biomedical-ngs/resource/bcr-vs-tcr-sequencing-immune-repertoire.html)  
4. Clinical application of single‐cell RNA sequencing in disease and therapy \- PMC, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC12576031/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12576031/)  
5. Analysis of single-cell TCR repertoires and gene expression from multi-modal scRNA-seq data \- DOI, 3月 2, 2026にアクセス、 [https://doi.org/10.1016/bs.mcb.2025.03.007](https://doi.org/10.1016/bs.mcb.2025.03.007)  
6. BCR Sequencing Applications \- Adaptive Biotechnologies, 3月 2, 2026にアクセス、 [https://www.adaptivebiotech.com/bcr-sequencing-applications/](https://www.adaptivebiotech.com/bcr-sequencing-applications/)  
7. Integrated TCR repertoire analysis and single-cell transcriptomic profiling of tumor-infiltrating T cells in renal cell carcinoma identifies shared and tumor-restricted expanded clones with unique phenotypes \- PMC, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC9515957/](https://pmc.ncbi.nlm.nih.gov/articles/PMC9515957/)  
8. Introduction to scRNA-Seq with R (Seurat) \- Bioinformatics \- National Cancer Institute, 3月 2, 2026にアクセス、 [https://bioinformatics.ccr.cancer.gov/docs/getting-started-with-scrna-seq/IntroToR\_Seurat/](https://bioinformatics.ccr.cancer.gov/docs/getting-started-with-scrna-seq/IntroToR_Seurat/)  
9. SeQuoIA: a single-cell BCR sequencing analysis pipeline for ..., 3月 2, 2026にアクセス、 [https://www.biorxiv.org/content/10.1101/2025.02.13.638132v2.full-text](https://www.biorxiv.org/content/10.1101/2025.02.13.638132v2.full-text)  
10. scRepertoire 2: Enhanced and efficient toolkit for single-cell immune profiling, 3月 2, 2026にアクセス、 [https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012760](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012760)  
11. Scirpy: single-cell immune receptor analysis in Python — scirpy, 3月 2, 2026にアクセス、 [https://scirpy.scverse.org/](https://scirpy.scverse.org/)  
12. Platypus: an open-access software for integrating lymphocyte single-cell immune repertoires with transcriptomes | bioRxiv, 3月 2, 2026にアクセス、 [https://www.biorxiv.org/content/10.1101/2020.11.09.374280v1.full-text](https://www.biorxiv.org/content/10.1101/2020.11.09.374280v1.full-text)  
13. Paired single-B-cell transcriptomics and receptor sequencing reveal activation states and clonal signatures that characterize B cells in acute myeloid leukemia \- PMC, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC10910691/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10910691/)  
14. scRepertoire 2: Enhanced and efficient toolkit for single-cell immune profiling \- PMC, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC12204475/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12204475/)  
15. Quantifying B-cell Clonal Diversity In Repertoire Data \- bioRxiv.org, 3月 2, 2026にアクセス、 [https://www.biorxiv.org/content/10.1101/2022.12.12.520133v2.full-text](https://www.biorxiv.org/content/10.1101/2022.12.12.520133v2.full-text)  
16. Calculate Clonal Diversity — clonalDiversity • scRepertoire \- Nick Borcherding, 3月 2, 2026にアクセス、 [https://www.borch.dev/uploads/screpertoire/reference/clonaldiversity](https://www.borch.dev/uploads/screpertoire/reference/clonaldiversity)  
17. Overview — dandelion documentation, 3月 2, 2026にアクセス、 [https://sc-dandelion.readthedocs.io/](https://sc-dandelion.readthedocs.io/)  
18. T Cell Clonal Analysis Using Single-cell RNA Sequencing and Reference Maps \- PMC, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC10450729/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10450729/)  
19. A comprehensive evaluation of diversity measures for TCR repertoire profiling \- PMC \- NIH, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC12080070/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12080070/)  
20. Public Clonotypes and Convergent Recombination Characterize the Naïve CD8+ T-Cell Receptor Repertoire of Extremely Preterm Neonates \- PubMed, 3月 2, 2026にアクセス、 [https://pubmed.ncbi.nlm.nih.gov/29312340/](https://pubmed.ncbi.nlm.nih.gov/29312340/)  
21. BCR/TCR Sequencing Repertoire Profiling | Human Immunology Core | Perelman School of Medicine at the University of Pennsylvania, 3月 2, 2026にアクセス、 [https://www.med.upenn.edu/humanimmunologycore/bcr-tcr-sequencing.html](https://www.med.upenn.edu/humanimmunologycore/bcr-tcr-sequencing.html)  
22. Diversity and Clonality | Platforma docs, 3月 2, 2026にアクセス、 [https://docs.platforma.bio/guides/vdj-analysis/diversity-analysis/](https://docs.platforma.bio/guides/vdj-analysis/diversity-analysis/)  
23. Repertoire-based mapping and time-tracking of T helper cell subsets in scRNA-Seq \- PMC, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC12006041/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12006041/)  
24. T-Cell Receptor Repertoire Analysis with Computational Tools—An Immunologist's Perspective \- PMC, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC8700004/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8700004/)  
25. Platypus: an open-access software for integrating lymphocyte single ..., 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC8046018/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8046018/)  
26. scRepertoire 2: Enhanced and Efficient Toolkit for Single-Cell Immune Profiling | bioRxiv, 3月 2, 2026にアクセス、 [https://www.biorxiv.org/content/10.1101/2024.12.31.630854v1.full-text](https://www.biorxiv.org/content/10.1101/2024.12.31.630854v1.full-text)  
27. Flexible Distance-Based TCR Analysis in Python with tcrdist3 \- ResearchGate, 3月 2, 2026にアクセス、 [https://www.researchgate.net/publication/363446286\_Flexible\_Distance-Based\_TCR\_Analysis\_in\_Python\_with\_tcrdist3](https://www.researchgate.net/publication/363446286_Flexible_Distance-Based_TCR_Analysis_in_Python_with_tcrdist3)  
28. Flexible Distance-Based TCR Analysis in Python with tcrdist3 \- PMC \- NIH, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC9719034/](https://pmc.ncbi.nlm.nih.gov/articles/PMC9719034/)  
29. Single Cell and Immunoanalysis Core \- The Human Skin Disease Resource Center, 3月 2, 2026にアクセス、 [https://humanskin.bwh.harvard.edu/human-analytic-techniques/singlecellimmunoanalysis/](https://humanskin.bwh.harvard.edu/human-analytic-techniques/singlecellimmunoanalysis/)  
30. Public Clonotypes and Convergent Recombination Characterize the Naïve CD8+ T-Cell Receptor Repertoire of Extremely Preterm Neonates \- Frontiers, 3月 2, 2026にアクセス、 [https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2017.01859/full](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2017.01859/full)  
31. Modeling and predicting the overlap of B- and T-cell receptor repertoires in healthy and SARS-CoV-2 infected individuals \- PMC, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC10075420/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10075420/)  
32. Multimodal T cell analysis with CoNGA \- PMC, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC9667881/](https://pmc.ncbi.nlm.nih.gov/articles/PMC9667881/)  
33. Integrating T-cell receptor and transcriptome for large-scale single-cell immune profiling analysis \- bioRxiv, 3月 2, 2026にアクセス、 [https://www.biorxiv.org/content/10.1101/2021.06.24.449733v2.full.pdf](https://www.biorxiv.org/content/10.1101/2021.06.24.449733v2.full.pdf)  
34. Integrated single-cell transcriptome and TCR profiles of hepatocellular carcinoma highlight the convergence on interferon signaling during immunotherapy, 3月 2, 2026にアクセス、 [https://jitc.bmj.com/content/12/11/e010534](https://jitc.bmj.com/content/12/11/e010534)  
35. Single-cell, Spatially-Resolved TCR Profiling Links T Cell Phenotype and Clonality in Human Tumors | bioRxiv, 3月 2, 2026にアクセス、 [https://www.biorxiv.org/content/10.1101/2025.08.06.668925v1.full-text](https://www.biorxiv.org/content/10.1101/2025.08.06.668925v1.full-text)  
36. Integrated TCR repertoire analysis and single-cell transcriptomic profiling of tumor-infiltrating T cells in renal cell carcinoma identifies shared and tumor-restricted expanded clones with unique phenotypes \- Frontiers, 3月 2, 2026にアクセス、 [https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.952252/full](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.952252/full)  
37. Integrative analysis in Seurat v5 \- Satija Lab, 3月 2, 2026にアクセス、 [https://satijalab.org/seurat/articles/seurat5\_integration](https://satijalab.org/seurat/articles/seurat5_integration)  
38. Integrating T cell receptor sequences and transcriptional profiles by clonotype neighbor graph analysis (CoNGA) \- PMC, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC8832949/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8832949/)  
39. phbradley/conga: Clonotype Neighbor Graph Analysis \- GitHub, 3月 2, 2026にアクセス、 [https://github.com/phbradley/conga](https://github.com/phbradley/conga)  
40. Linking T cell receptor sequence to transcriptional profiles with clonotype neighbor graph analysis (CoNGA) \- bioRxiv.org, 3月 2, 2026にアクセス、 [https://www.biorxiv.org/content/10.1101/2020.06.04.134536v1.full-text](https://www.biorxiv.org/content/10.1101/2020.06.04.134536v1.full-text)  
41. LRT: T Cell Trajectory Inference by Integrative Analysis of Single-Cell TCR-seq and RNA-seq data | bioRxiv, 3月 2, 2026にアクセス、 [https://www.biorxiv.org/content/10.1101/2022.04.14.488320v2.full-text](https://www.biorxiv.org/content/10.1101/2022.04.14.488320v2.full-text)  
42. LRT: Integrative analysis of scRNA-seq and scTCR-seq data to investigate clonal differentiation heterogeneity \- PubMed, 3月 2, 2026にアクセス、 [https://pubmed.ncbi.nlm.nih.gov/37428794/](https://pubmed.ncbi.nlm.nih.gov/37428794/)  
43. Ingesting multimodal (CITE-Seq \+ VDJ) datasets from cellranger outputs, 3月 2, 2026にアクセス、 [https://panpipes-tutorials.readthedocs.io/en/latest/ingesting\_multimodal\_data/ingesting\_multimodal\_data.html](https://panpipes-tutorials.readthedocs.io/en/latest/ingesting_multimodal_data/ingesting_multimodal_data.html)  
44. Divergent clonal differentiation trajectories of T cell exhaustion \- PMC, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC11225711/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11225711/)  
45. Single-cell immune repertoire analysis 1 2 Sergio E ... \- UQ eSpace, 3月 2, 2026にアクセス、 [https://espace.library.uq.edu.au/view/UQ:5a08d41/scVDJseq\_review\_accepted\_version.pdf](https://espace.library.uq.edu.au/view/UQ:5a08d41/scVDJseq_review_accepted_version.pdf)  
46. Using Seurat with multimodal data \- Satija Lab, 3月 2, 2026にアクセス、 [https://satijalab.org/seurat/articles/multimodal\_vignette.html](https://satijalab.org/seurat/articles/multimodal_vignette.html)  
47. Integrated analysis of multimodal single-cell data with structural similarity \- PMC, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC9757079/](https://pmc.ncbi.nlm.nih.gov/articles/PMC9757079/)  
48. A toolkit for single-cell immune receptor profiling • scRepertoire \- borch.dev, 3月 2, 2026にアクセス、 [https://www.borch.dev/uploads/screpertoire/](https://www.borch.dev/uploads/screpertoire/)  
49. Combining Deep Learning and TCRs with Trex • scRepertoire \- Nick Borcherding, 3月 2, 2026にアクセス、 [https://www.borch.dev/uploads/screpertoire/articles/trex](https://www.borch.dev/uploads/screpertoire/articles/trex)  
50. Evaluation of cell type annotation reliability using a large language model-based identifier, 3月 2, 2026にアクセス、 [https://pmc.ncbi.nlm.nih.gov/articles/PMC12462508/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12462508/)  
51. Large Language Model Consensus Substantially Improves the Cell Type Annotation Accuracy for scRNA-seq Data | bioRxiv, 3月 2, 2026にアクセス、 [https://www.biorxiv.org/content/10.1101/2025.04.10.647852v1.full-text](https://www.biorxiv.org/content/10.1101/2025.04.10.647852v1.full-text)  
52. Which is the Best method for cell type annotation for single cell characterization ? | ResearchGate, 3月 2, 2026にアクセス、 [https://www.researchgate.net/post/Which\_is\_the\_Best\_method\_for\_cell\_type\_annotation\_for\_single\_cell\_characterization](https://www.researchgate.net/post/Which_is_the_Best_method_for_cell_type_annotation_for_single_cell_characterization)

[image1]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAHkAAAAWCAYAAADkWDPGAAABE0lEQVR4Xu3YMWpCQRAG4BFTBLSUQKqAhZ2V2FkKdl4iRU6SC+QO6XIEC8E+XS5gY5UyjYI6/5tdngzLMwgiI/8HP/j8WZtx1/cUISIiIqKb+NAcCuml3r+PLFJHQXQ1z5pfsQEO03WG17PUrdI11lBAeZeWvGt2mokvKI4HsQFvfCG2a5eaH6mPcApoLDZk7FhvLta9+YJiwQD3mqkvxAaPIeOLQEHhCMZRvNWsC8m/1a284AzsfP8ZTfnWDKqVdDW4mcJNVemohqYbMgriS2yIfV+ojliHnU6BYYAYJAbqYfDoPn1BsTQdx5c8Hz+K/WHy3zxp2tVKugo+H9+xF6l38GleNSPNX6ErPV4RERERUYMjV0tP4AsTA/EAAAAASUVORK5CYII=>

[image2]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAA0AAAAXCAYAAADQpsWBAAAA+klEQVR4XmNgGAWjgGqAE4hdgZgVyucAYgcgdoMpwAZ2AvE+IL4IxHJA/AmItwHxHSCOA2JGhFIIEAFiSyBmAeL/QLwHVRosthxNjMEFiAUZIJpBCqJRpcFih5H4hUhshnQGiCQvkhjIMHTbHWAMkJuXAvEkuBQEGDNANLWiiYOBKRB/A2JNJDFQiO4A4ukMEP9OAOJaIO6FKQA5DWQisqZgIP4LxPJQfg0Q8wHxSpgCkNNAmmYxQOLKCojfA/EVmAIg0GeAhHIETOA3EL9lgGgQZ4BELjoAhe5VBqQ4A9lyGi6NHXgA8S8gVgdxpBmwxw868APiJ0BcAAASuSrAkxKx5AAAAABJRU5ErkJggg==>

[image3]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAYAAAAVCAYAAABljp99AAAAeElEQVR4XmNgGAqAE4gjgJgRXWI6EP8HYk10CUkgdkYXxAlkoRgF7AFiISD+CsRXYYJiQGwLZYMshktwQGklIP4NxDYwCRjIAeITQMyPLGgJxD+BWAfK54FJtDJAzAcZCxL0g0kcZoC4COSy7TBBGACZLYwuSDwAANlwEA4sQUchAAAAAElFTkSuQmCC>

[image4]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGYAAAAWCAYAAAAy/emjAAABU0lEQVR4Xu2YMUoDQRSGn4hgMI1aiGATcgAFr2BhIx5AsPYAkoB9Wg9gk1NIJFfIAawVvYBgY6H+P282bB67y06zQzLvg6/IvF3YnT87O29FHMdxHMfJmgn8a/AHTpdHp+ce9u3gJrILj+GnaBDD8LtwHsZv4VY4JwVHcATfRK8rGzj5v3Yw0BOtv9hCBxzAa3gIx5JZMHwSOPEftlCC9W872DHZBXMqOvG88TqKd05KsgvmRnQZu7CFEh5Mx+zDBXwW3QhUwWPaLGUn8D3SGLIKhsvYF3ywhRLnosEwwJRkFQx7FE46/+1VcEfG3Rg3BgNT65qYYHg/j/CppSlbgUr4FDCYumXsUvT9c2cLCYgJZu1p6l/YQ7A+s4UatmW1OW1jDNkE09S/nMFX0aVuz9RSsfHB8MaK7W+V/D52BXeKExLCa2UY9hppU9/lOI7jOE7gH38qZEdTiUxlAAAAAElFTkSuQmCC>

[image5]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAFkAAAAWCAYAAACrBTAWAAAA2klEQVR4Xu3XMQ4BQRiG4V9CQqMhEaF2Ax2NU+iU3ECtVytEJREncBTRKMQhhEjw/dkZ2Z2GtTOa/Z7kLcZMNVkzuyJEREQkdbRHT6dVfBFlU0A11EYL9EAjVI4vIj/sE70TbnAwXXRBU3eC/NlIdBa33Anyx15+PCoC0g3WS+8bVXROWe7pG4Zu8sn5nTzSS083mZdeQGN0Rz13gn5zk+ip7ZhxAx3QFhXtog/0eGmmLFfsZ3PFjCdm3H+voMzm6IqGaI2OaJBYQZnpX12PiiWaoVJymoiI/uAFXp0tzou0nckAAAAASUVORK5CYII=>

[image6]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAsAAAAXCAYAAADduLXGAAAA3klEQVR4Xu2RvQ4BQRSFryAhOp3CE0h0GpFso6DQ8QxKHU+gUBKFaLYSiTeQKHRKlTfQiEJIKIifc3N3ZndupdH5ki+TPedOZneH6M93NOAMjmBVdQ4ePMMhnMM7bDoToACPcKUL8NbBGL5gTRfgooMrfMCKLsBWBzeS4xa6AFkd7EmG2SVsubVLGg7gk8JN7Ck6pEnBCTxQuCFhyiQsmwdFn2Q4Z4IS3NjapUcynDGBD3e2dlmTXJSFB/kyNDGSvBsN+TI4rMN4kPHagVOSb7IUgzUP2yQD/I/te/6WD6d0LXoXDa+zAAAAAElFTkSuQmCC>