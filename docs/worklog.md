# 作業日誌

---

## 2026-05-25

### 作業内容

**リポジトリの立ち上げと課題理解**

1. **課題の把握**
   - 2次元強磁性イジング模型の転移温度をモンテカルロ法で求める課題の内容を整理
   - 求められる成果物を確認：MC実装・有限サイズスケーリング・$O(L^2)$確認・ユニットテスト

2. **overview.md の作成・改訂**
   - 課題の概要ノートを作成（モデル定義・手法・成果物チェックリスト）
   - オンサーガーの厳密解の説明を追加
   - モンテカルロ法（MCMC・重要サンプリング・詳細釣り合い・メトロポリス法）の説明を追加
   - 説明が冗長との指摘を受け、箇条書き中心の簡潔な形式に改訂

3. **background.md の作成・拡充**
   - 歴史・理論的背景を overview.md から分離して別ファイルに移動
   - オンサーガー解（転送行列法）・MCMC理論・有限サイズスケーリングの理論的背景を収録
   - 有限サイズスケーリング節を大幅拡充：
     - スケーリング仮説と各物理量のスケーリング形
     - 2Dイジングの臨界指数一覧表（$\nu=1,\ \beta=1/8,\ \gamma=7/4,\ \alpha=0$）
     - Binder キュムラントの交差原理の導出
     - $T_c$ 推定の実践手順（磁化率ピーク外挿法・Binder キュムラント交差法）

4. **数式レンダリングの修正**
   - ブロック数式を `$$数式$$`（1行）から `$$\n数式\n$$`（3行）形式に統一
   - `.vscode/settings.json` を作成し VS Code のマークダウン数式表示を有効化

5. **README.md・CLAUDE.md の整備**
   - README.md：ファイル構成・ドキュメントの読み方・数式表示方法を追記
   - CLAUDE.md：新規作成。ファイル構成・コードスタイル・ドキュメント方針・Git ルールを記述

6. **GitHub へのアップロード**
   - public リポジトリ `Sayuri-Miyake/AI_kadai1` を作成してプッシュ
   - URL：https://github.com/Sayuri-Miyake/AI_kadai1

### 現在のファイル構成

```
AI_kadai1/
├── CLAUDE.md
├── README.md
├── overview.md
├── background.md
├── docs/
│   └── worklog.md    ← このファイル
└── .vscode/
    └── settings.json
```

### 次回以降の作業予定

- [ ] シミュレーションコードの実装（Python 推奨）
- [ ] ユニットテストの作成
- [ ] 数値実験の実行と結果の整理
- [ ] 有限サイズスケーリング解析

---

## 2026-06-01

### 作業内容

**Rust によるモンテカルロシミュレーション実装と可視化**

1. **`simulation/` ディレクトリの新規作成（Rust プロジェクト）**
   - `Cargo.toml`：クレート設定（依存：`rand`・`csv`）
   - `src/ising.rs`：`IsingModel` 構造体の実装
     - 指数テーブルの事前計算（$\exp(-\Delta E / k_B T)$）でホットループを高速化
     - 周期境界条件の実装
     - メトロポリス法による1スピンフリップ更新
   - `src/observables.rs`：`Accumulator` 構造体の実装
     - 磁化 $m$・磁化率 $\chi$・比熱 $C$・Binder キュムラント $U_4$ の逐次集計
   - `src/main.rs`：温度スキャンのメインループ
     - $L = 8, 16, 32, 64$ の4サイズを対象
     - 各温度点でサーモ化ステップ後に物理量を測定し CSV 出力

2. **ユニットテストの作成（11件）**
   - エネルギー計算・$\Delta E$ 計算の正確性テスト
   - 低温・高温極限における磁化の振る舞いテスト
   - 全テストが `cargo test` でパスすることを確認

3. **`simulation/results.html`：インタラクティブ可視化の作成**
   - Plotly.js を用いた4パネル表示（磁化・磁化率・比熱・Binder キュムラント）
   - ブラウザ単体で動作するスタンドアロン HTML

4. **`background.md` の MCMC 節を拡充**
   - マルコフ連鎖理論の記述追加
   - Glauber ダイナミクス vs. メトロポリス法の比較
   - 自己相関時間と臨界スローダウン（$\tau \sim L^z$）の説明追加

5. **`.gitignore` の追加**
   - Rust ビルド成果物（`target/`）と生成 CSV ファイルを除外

### 現在のファイル構成

```
AI_kadai1/
├── CLAUDE.md
├── README.md
├── overview.md
├── background.md
├── .gitignore
├── docs/
│   └── worklog.md    ← このファイル
├── .vscode/
│   └── settings.json
└── simulation/       ← 新規追加
    ├── Cargo.toml
    ├── Cargo.lock
    ├── results.html
    └── src/
        ├── main.rs
        ├── ising.rs
        └── observables.rs
```

### 次回以降の作業予定

- [ ] シミュレーションの実行と CSV データの取得
- [ ] 有限サイズスケーリング解析（$T_c$ 推定・臨界指数の確認）
- [ ] ユニットテストのカバレッジ拡充
- [ ] Python / Jupyter での解析スクリプト作成（必要に応じて）
