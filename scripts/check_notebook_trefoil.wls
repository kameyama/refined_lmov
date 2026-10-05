(* 1703.05408 のトーラス用ノートブックで、三葉結び目（|μ|=3 は a^(1/2) を掛けて a の冪を整数にしてから数える） T(2,3) の f^q, f^{\bar t}, \hat f を knotparam で太字にした式の係数の符号を数える。
   ノートブックの定義セル（関数定義と代入）だけを評価し、Results の三葉結び目のセルと同じ式を評価する。
   実行：just wolfram scripts/check_notebook_trefoil.wls *)

nb = FileNameJoin[{Directory[], "references", "1703.05408", "src", "anc", "refined-LMOV-torus.nb"}];
cells = NotebookImport[nb, "Input" -> "HeldExpression"];
(* definition cells: input cells 2-10 (Mq2, Mq3, torus-knot data, fq/ft, RFq2/RFq3, knotparam, gfactor) *)
Scan[ReleaseHold, cells[[2 ;; 10]]];

count[expr_] := Module[{num, rules, neg, pos},
  num = Expand[Numerator[Together[PowerExpand[expr]]]];
  rules = CoefficientRules[num, {A, Q, T}];
  neg = Count[rules[[All, 2]], c_ /; c < 0];
  pos = Count[rules[[All, 2]], c_ /; c > 0];
  {Length[rules], neg, pos}];

show[label_, list_, names_] := Do[
  Module[{c = count[list[[i]]]},
   Print[label, " ", names[[i]], ": terms ", c[[1]], ", negative ", c[[2]], ", positive ", c[[3]],
    If[c[[2]] > 0 && c[[3]] > 0, "   <-- mixed signs", ""]]], {i, Length[list]}];

prep[x_] := Simplify[x // FunctionExpand // PowerExpand // Together // Factor, a > 0 && q > 0 && t > 0] // knotparam // PowerExpand;

Print["T(2,3) in the torus notebook (Results > Trefoil (T23)); knotparam = {a -> -T A^2, q -> Q^2 T^2, t -> Q^2}"];
show["input cell 23: f^q   (RFq2[1])", prep /@ RFq2[1], {"[2]", "[1,1]"}];
show["input cell 26: f^t   (RFt2[1])", prep /@ RFt2[1], {"[2]", "[1,1]"}];
show["input cell 24: f^q   (a^(1/2) RFq3[1])", prep /@ (a^(1/2) RFq3[1]), {"[3]", "[2,1]", "[1,1,1]"}];
show["input cell 28: f^t   (a^(1/2) RFt3[1])", prep /@ (a^(1/2) RFt3[1]), {"[3]", "[2,1]", "[1,1,1]"}];
show["input cell 29+: hat f (IMq2.RFq2[1])", prep /@ (IMq2 . RFq2[1]), {"[2]", "[1,1]"}];
show["hat f (IMq3.RFq3[1]/a^(1/2))", prep /@ (IMq3 . RFq3[1]/a^(1/2)), {"[3]", "[2,1]", "[1,1,1]"}];
