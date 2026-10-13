(* 1703.05408 の ancillary ノートブックから、トーラス結び目と 4_1 の refined CS 不変量（reduced）を書き出す。
   出力：data/framed_knot/reference/notebook_torus.txt
   各行：knot<TAB>color<TAB>式（a = A^2, q = Q^2, t = T^2 と置いた Laurent 多項式、InputForm）
   reduced 値 = ノートブックの unreduced 値 / 同じ彩色の unknot（torusknot*[.., 0, ..]）。
   実行：just wolfram scripts/framed_knot/export_notebook_torus.wls *)

nb = FileNameJoin[{Directory[], "references", "1703.05408", "src", "anc", "refined-LMOV-non-torus.nb"}];
cells = NotebookImport[nb, "Input" -> "HeldExpression"];
ReleaseHold[cells[[3]]];   (* データの定義セル（torusknotsym, torusknotanti, torusknot21, tripletcT34 など） *)
ReleaseHold[cells[[4]]];   (* 非トーラス結び目のデータ（4_1：singboxfig, horizontal2fig, ...） *)

toLaurent[expr_] := Module[{e},
  e = FunctionExpand[expr] /. {a -> A^2, q -> Q^2, t -> T^2};
  e = PowerExpand[e];
  Expand[Factor[Together[e]]]];

unknotSym[n_] := torusknotsym[n, 0, a, q, t];
unknotAnti[n_] := torusknotanti[n, 0, a, q, t];
unknot21 := torusknot21[0, a, q, t];

rows = {};
add[knot_, color_, value_] := AppendTo[rows, knot <> "\t" <> color <> "\t" <> ToString[toLaurent[value], InputForm]];

Do[
  k = "T2_" <> ToString[2 p + 1];
  add[k, "1", torusknotsym[1, p, a, q, t]/unknotSym[1]];
  add[k, "2", torusknotsym[2, p, a, q, t]/unknotSym[2]];
  add[k, "11", torusknotanti[2, p, a, q, t]/unknotAnti[2]];
  add[k, "3", torusknotsym[3, p, a, q, t]/unknotSym[3]];
  add[k, "111", torusknotanti[3, p, a, q, t]/unknotAnti[3]];
  add[k, "21", torusknot21[p, a, q, t]/unknot21],
  {p, 1, 4}];

add["T3_4", "1", tripletcT34[1, a, q, t]/unknotSym[1]];
add["T3_4", "2", tripletcT34[2, a, q, t]/unknotSym[2]];
add["T3_4", "11", tripletcT34anti2[a, q, t]/unknotAnti[2]];

(* 4_1（figure-eight）：§6 の入力。[2,1] は KN が構成した候補（式 (fig8-21)） *)
add["4_1", "1", singboxfig[a, q, t]/unknotSym[1]];
add["4_1", "2", horizontal2fig[a, q, t]/unknotSym[2]];
add["4_1", "11", vertical2fig[a, q, t]/unknotAnti[2]];
add["4_1", "3", horizontal3fig[a, q, t]/unknotSym[3]];
add["4_1", "111", vertical3fig[a, q, t]/unknotAnti[3]];
add["4_1", "21", figeight21fig[a, q, t]/unknot21];

out = FileNameJoin[{Directory[], "data", "reference", "notebook_torus.txt"}];
Quiet[CreateDirectory[DirectoryName[out]]];
Export[out, StringRiffle[rows, "\n"], "Text"];
Print["wrote ", Length[rows], " rows to ", out];
