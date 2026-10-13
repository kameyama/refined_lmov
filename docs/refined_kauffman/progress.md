# 進捗記録（第 3 の研究：refined Kauffman の large N duality）

## 現在の状況（再開用、2026-10-14）
- 開始。ユーザーの指示に従い、(1) unrefined の Kauffman の large N duality の文献調査、(2) 2018 年の試行錯誤（Dropbox の Mathematica ノートブック 42 本）のミスの調査、の順に進める。
- ブランチ：claude/refined-kauffman-large-n-00d58e（worktree）。

## 2026-10-14
- ディレクトリ docs/refined_kauffman/ を作成（README.md、progress.md）。
- 過去のノートブック 42 本（2018-07-10〜10-02）を scripts/nb2txt.py でテキスト化して読解を開始。同梱の PDF は Cherednik 1111.6195（DAHA-Jones、B/C 型の hyperpolynomial）と Cherednik–Elliot 1503.01441（refined composite invariants of torus knots via DAHA）（arXiv API で題名を確認）。
- 文献の取得：justfile の papers に 0904.1088（Mariño）、hep-th/0405083・hep-th/0411227（Bouchard–Florea–Mariño）、hep-th/0505008（Borhade–Ramadevi）、1003.5282（Paul–Borhade–Ramadevi）、1003.2861（Stevan）、hep-th/0512298（Gukov–Walcher）、1202.4456（Aganagic–Schaeffer）、1111.6195、1503.01441、1310.2240（Nawata–Ramadevi–Zodinmawia）、1007.1656（Chen–Chen）、hep-th/0108064（Mariño–Vafa）を追加して `just fetch-refs` を実行（題名・内容は原典で確認してから literature.md に書く）。
