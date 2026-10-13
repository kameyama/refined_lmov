#!/bin/zsh
cd /Users/masaya/projects/ai_research/refined_cs/.claude/worktrees/rename-literature-survey-doc-9cc5fc
export PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov
PY=/Users/masaya/projects/ai_research/refined_cs/.venv/bin/python
echo "=== $(date) K10 の 3 箱の開始（高速化の後、commit 646422a 以降）"
/usr/bin/time -l $PY scripts/framed_knot/v0_iterated.py --workers 3 3 21 111 && \
/usr/bin/time -l $PY scripts/framed_knot/v0_K10_lmov.py --workers 3 && \
echo "=== $(date) K10 の 3 箱：完了"
