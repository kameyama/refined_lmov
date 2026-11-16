"""検出スクリプトと不変量スクリプトの機能テスト。

実行：python3 -m pytest .claude/skills/hep-th-paper/tests
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent
SCRIPTS = HERE.parent / "scripts"


def run(script, *args):
    return subprocess.run([sys.executable, str(SCRIPTS / script), *map(str, args)], capture_output=True, text=True, cwd=HERE)


def test_tells_on_ai_paragraph():
    r = run("check_tells.py", "fixture_ai.tex", "--bib", "fixture.bib")
    assert r.returncode == 1
    for expected in [
        "[hard 誠実性] \\cite のキーが .bib にない: TODO-missing",
        "DOI も eprint もない（実在を原典で確認）: Macdonald95",
        "出典のない「よく知られている」",
        "内容のない導入: In recent years",
        "まとめ文: Overall",
        "not only ... but also",
        "ヘッジの重なり: may potentially",
        "意義づけで終わる文: underscoring",
        "[soft 語彙] 同じ段落に 12 語",  # コメントの中の 2 語は数えず、\% のあとの delve は数える
    ]:
        assert expected in r.stdout, expected
    assert "KN17" not in r.stdout  # .bib にあり、DOI もあるキーは指摘しない


def test_tells_after_fix():
    r = run("check_tells.py", "fixture_fixed.tex", "--bib", "fixture.bib")
    assert "語彙 0 回" in r.stdout
    for gone in ["内容のない導入", "まとめ文", "not only", "ヘッジ", "意義づけ", "[soft 語彙]"]:
        assert gone not in r.stdout, gone
    assert "出典のない「よく知られている」" in r.stdout  # 誠実性の項目は直さずに残る


def test_voice_removes_the_authors_usual_words(tmp_path):
    (tmp_path / "old.tex").write_text("We delve into a novel and robust construction.\n" * 3, encoding="utf-8")
    r = run("check_tells.py", "fixture_ai.tex", "--bib", "fixture.bib", "--voice", tmp_path)
    assert "同じ段落に 8 語" in r.stdout  # 12 語から delve 2、novel 1、robust 1 が外れる


def test_abstract_checks_and_proper_names(tmp_path):
    long_sentence = "We compute the invariants of torus knots. " * 50
    (tmp_path / "a.tex").write_text(
        "\\begin{abstract}\nWe refine the duality of \\cite{KN17}. " + long_sentence + "\n\\end{abstract}\n"
        "The result is reconstructed by Garner's algorithm, a pivotal step.\n", encoding="utf-8")
    r = run("check_tells.py", tmp_path / "a.tex", "--all")
    assert "[hard 構造] abstract の中の \\cite（著者名だけを書く）: KN17" in r.stdout
    assert "arXiv の abstract 欄は 1,920 文字まで" in r.stdout
    assert "単発: pivotal" in r.stdout and "garner" not in r.stdout.lower()  # 人名の Garner は語彙に数えない


def test_invariants_pass_for_a_prose_only_fix():
    r = run("check_invariants.py", "--before", "fixture_ai.tex", "fixture_fixed.tex")
    assert r.returncode == 0, r.stdout


def test_invariants_catch_changes(tmp_path):
    changed = (HERE / "fixture_fixed.tex").read_text(encoding="utf-8")
    changed = changed.replace("q^{-1/2}", "q^{1/2}").replace("12 knots", "13 knots").replace("in \\cite{KN17}", "before")
    (tmp_path / "changed.tex").write_text(changed, encoding="utf-8")
    r = run("check_invariants.py", "--before", "fixture_ai.tex", tmp_path / "changed.tex")
    assert r.returncode == 1
    for name in ["数式・表: 不一致", "引用: 不一致", "本文の数値: 不一致"]:
        assert name in r.stdout
    assert "参照: 一致" in r.stdout
