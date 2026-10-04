set shell := ["bash", "-euo", "pipefail", "-c"]

# サンドボックス外で実行するのは、このファイルのレシピだけ。
# 作業ディレクトリはこの justfile のあるディレクトリ（リポジトリのルート）。

papers := "1105.5117 1202.2489 1703.05408v3 1106.4305 1201.3339 1304.3328 1408.4348 1509.08351 2410.14703 1410.0859 1605.00978 1211.4596 math/0505662"

# レシピ一覧
default:
    just --list

# .venv を作り、rcsknots を editable インストールする（PyPI へ接続）
setup:
    python3 -m venv .venv
    .venv/bin/pip install -q --upgrade pip
    .venv/bin/pip install -q -e ".[dev]"

# arXiv から論文の PDF と TeX ソース（ancillary ファイルを含む）を references/<id>/ に取得して展開する
fetch-refs:
    for id in {{papers}}; do \
      base="${id%v*}"; base="${base//\//_}"; dir="references/$base"; mkdir -p "$dir"; \
      [ -s "$dir/$base.pdf" ] || { curl -sSfL -A "refined-cs-research" -o "$dir/$base.pdf" "https://arxiv.org/pdf/$id"; sleep 3; }; \
      [ -s "$dir/src.tar" ] || { curl -sSfL -A "refined-cs-research" -o "$dir/src.tar" "https://arxiv.org/e-print/$id"; sleep 3; }; \
      if [ ! -d "$dir/src" ]; then \
        mkdir -p "$dir/src"; \
        case "$(file -b "$dir/src.tar")" in \
          gzip*) tar xzf "$dir/src.tar" -C "$dir/src" 2>/dev/null || gunzip -c "$dir/src.tar" > "$dir/src/main.tex" ;; \
          *tar*) tar xf "$dir/src.tar" -C "$dir/src" ;; \
        esac; \
      fi; \
    done

# Wolfram Language スクリプトを実行する（ライセンス読み込みのためサンドボックス外で動かす）
wolfram file:
    wolframscript -file {{file}}
