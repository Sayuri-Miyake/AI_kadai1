# 2D Ferromagnetic Ising Model — Monte Carlo Simulation

2次元強磁性イジング模型の転移温度を古典モンテカルロ法（メトロポリス法）で求める課題です。

## 課題内容

- $L \times L$ 正方格子上のイジング模型をモンテカルロシミュレーションで解析
- 有限サイズスケーリングにより熱力学極限の転移温度 $T_c$ を推定
- オンサーガーの厳密解 $T_c/J \approx 2.2692$ と比較して精度を検証
- 計算時間が $O(L^2)$ でスケールすることを確認
- 小さな $L$ に対するユニットテストを実装

## ファイル構成

```
AI_kadai1/
├── README.md         # このファイル
├── CLAUDE.md         # Claude Code 向けプロジェクト設定
├── overview.md       # 課題の要点（モデル・手法・成果物チェックリスト）
├── background.md     # 背景・歴史的解説（オンサーガー解・MCMC・有限サイズスケーリング理論）
├── docs/
│   └── worklog.md    # 作業日誌
└── .vscode/
    └── settings.json # VS Code マークダウン数式レンダリング設定
```

## ドキュメントの読み方

| ファイル | 用途 |
|----------|------|
| [overview.md](overview.md) | まずここ。課題に必要な式・手順・チェックリスト |
| [background.md](background.md) | 余裕があれば。理論的背景・歴史・導出 |

## 数式の表示について

マークダウン中の数式（`$...$` / `$$...$$`）を表示するには以下のいずれかが必要です：

- **GitHub**：このリポジトリをブラウザで開けば自動で表示
- **VS Code**：バージョン 1.88 以降。プロジェクト内の `.vscode/settings.json` で設定済み

## 参考文献

- L. Onsager, Phys. Rev. **65**, 117 (1944)：2Dイジング模型の厳密解
- N. Metropolis et al., J. Chem. Phys. **21**, 1087 (1953)：メトロポリスアルゴリズム
- K. Binder, Z. Phys. B **43**, 119 (1981)：Binder キュムラント法
