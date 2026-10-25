(* 別研究 [2,1] の X1：KNTZ 1902.02275 の補助ファイル HOMFLYPTinvariants.m をそのまま実行し、照合用の値を書き出す。 *)
(* 実行：just wolfram scripts/framed_knot/x1_kntz_package.wls（サンドボックス外）。出力：data/framed_knot/reference/kntz_package/values.txt *)
(* 行：種類<TAB>引数<TAB>InputForm の値。種類 F = RefinedFFactor[λ, m, s]、P = PoincarePolynomial`DoubleBraidPolynomial[r, s, m, n]、 *)
(*     H = HOMFLYPT`DoubleBraidPolynomial[r, s, m, n]（非 refined） *)
dir = DirectoryName[$InputFileName] <> "../";
Get[dir <> "references/1902.02275/src/anc/HOMFLYPTinvariants.m"];
out = OpenWrite[dir <> "data/framed_knot/reference/kntz_package/values.txt"];
emit[kind_, args_, val_] := (WriteString[out, kind <> "\t" <> ToString[args] <> "\t" <> ToString[InputForm[val]] <> "\n"]; Print[kind, " ", args]);
Do[emit["F", {lam, m, Length[lam]}, Factor[Together[HOMFLYPTinvariants`TwistFactor`RefinedFFactor[lam, m, Length[lam]]]]],
  {lam, {{1}, {2}, {1, 1}, {2, 1}, {2, 2}}}, {m, {2, 3}}];
Do[emit["P", {rs[[1]], rs[[2]], mn[[1]], mn[[2]]},
    Expand[Together[HOMFLYPTinvariants`PoincarePolynomial`DoubleBraidPolynomial @@ Join[rs, mn]]]],
  {rs, {{1, 1}, {2, 1}, {1, 2}, {2, 2}}}, {mn, {{1, 1}, {-1, 1}, {2, 1}, {-2, 1}, {2, 2}, {-2, 2}}}];
Do[emit["H", {rs[[1]], rs[[2]], mn[[1]], mn[[2]]},
    Expand[Together[HOMFLYPTinvariants`HOMFLYPT`DoubleBraidPolynomial @@ Join[rs, mn]]]],
  {rs, {{2, 1}, {2, 2}}}, {mn, {{-1, 1}, {2, 1}}}];
Close[out];
Print["kntz package done"];
