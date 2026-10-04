set shell := ["bash", "-euo", "pipefail", "-c"]

# サンドボックス外で実行するのは、このファイルのレシピだけ。
# 作業ディレクトリはこの justfile のあるディレクトリ（refined_cs/）。

papers := "1105.5117 1202.2489 1703.05408v3 1106.4305 1201.3339 1304.3328 1408.4348 1509.08351 2410.14703 1410.0859 1605.00978 1211.4596"

# レシピ一覧
default:
    just --list

# Python venv を作り、依存パッケージを入れる（PyPI へ接続）
setup:
    python3 -m venv .venv
    .venv/bin/pip install -q --upgrade pip
    .venv/bin/pip install -q -r requirements.txt

# arXiv から論文の PDF と TeX ソース（ancillary ファイルを含む）を refs/ に取得する
fetch-refs:
    mkdir -p refined_cables/refs
    for id in {{papers}}; do \
      base="${id%v*}"; dir="refined_cables/refs/$base"; mkdir -p "$dir"; \
      [ -s "$dir/$base.pdf" ] || curl -sSfL -A "refined-cs-research" -o "$dir/$base.pdf" "https://arxiv.org/pdf/$id"; \
      [ -s "$dir/src.tar" ] || curl -sSfL -A "refined-cs-research" -o "$dir/src.tar" "https://arxiv.org/e-print/$id"; \
      sleep 3; \
    done
