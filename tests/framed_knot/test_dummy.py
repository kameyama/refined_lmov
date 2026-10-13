"""環境確認用のダミーテスト：パッケージと依存ライブラリが import できること。"""


def test_imports():
    import flint
    import mpmath
    import numpy
    import sympy

    import rcsknots

    assert rcsknots.__doc__
    assert flint.nmod(3, 7) * 5 == 1
    assert numpy and mpmath and sympy
