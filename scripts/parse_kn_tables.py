"""1703.05408 付録 C の BPS 表（\\hat N_{ρ,g,β,J_r}）を TeX ソースから読み取り、data/reference/kn_tables.json に保存する。

各表：caption の (ρ, 結び目) と、{(g, β, 2J_r): N}（0 の項は省く）。
"""
import json
import re
from pathlib import Path

SRC = Path("references/1703.05408/src/refined-integrality.tex")
OUT = Path("data/reference/kn_tables.json")


def cells(row):
    return [c.strip() for c in row.split("&")]


def parse_tabular(body):
    rows = [r.strip() for r in re.split(r"\\\\", body)]
    rows = [r.replace("\\hline", "").strip() for r in rows]
    rows = [r for r in rows if r]
    beta_row = next(r for r in rows if "\\beta" in cells(r)[0])
    jr_row = next(r for r in rows if "2J_r" in r)
    betas = []
    for c in cells(beta_row)[1:]:
        m = re.match(r"\\multicolumn\{(\d+)\}\{[^}]*\}\{(-?\d+(?:/\d+)?)\}", c)
        n, b = (int(m.group(1)), m.group(2)) if m else (1, c)
        betas += [b] * n
    jrs = cells(jr_row)[1:]
    assert len(jrs) == len(betas), (len(jrs), len(betas))
    data = {}
    for r in rows[rows.index(jr_row) + 1:]:
        cs = cells(r)
        g = int(cs[0])
        for b, j, v in zip(betas, jrs, cs[1:]):
            v = int(v)
            if v:
                data[f"{g},{b},{j}"] = v
    return data


def main():
    text = SRC.read_text()
    tables = []
    for m in re.finditer(r"\\begin\{tabular\}\{[^}]*\}(.*?)\\end\{tabular\}(.*?)\\caption\{([^\n]*)", text, re.S):
        body, caption = m.group(1), m.group(3)
        cm = re.search(r"\\widehat\{\{?N\}?\}_\{([^}]*)\}\((.*?)\)", caption)
        if not cm:
            continue
        sub, knot = cm.group(1), cm.group(2)
        rho = re.findall(r"\[([\d,]+)\]", sub)
        knot = knot.replace("\\mathbf{", "").replace("}", "").replace("T_{", "T").replace(",", "_")
        tables.append({"rho": rho, "knot": knot, "caption": caption, "data": parse_tabular(body)})
    OUT.write_text(json.dumps(tables, indent=1))
    for tb in tables:
        print(tb["knot"], tb["rho"], len(tb["data"]))


if __name__ == "__main__":
    main()
