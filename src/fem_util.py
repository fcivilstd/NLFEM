import numpy as np
import math

class Util:
    """
    計算ツール
    """

    def __init__(self) -> None:
        """
        コンストラクタ
        """
        pass

    @staticmethod
    def t2dott2(s: np.ndarray, t: np.ndarray) -> float:
        """
        2階のテンソルの内積

        Args:
            s (numpy.ndarray): 列ベクトル表示した2階のテンソル
            t (numpy.ndarray): 列ベクトル表示した2階のテンソル

        Returns:
            float: 計算した内積
        """
        tmp: float = s.T[0][3] * t[3][0] + s.T[0][4] * t[4][0] + s.T[0][5] * t[5][0]
        return (s.T @ t)[0][0] + tmp

    @staticmethod
    def t2xt2(s: np.ndarray, t: np.ndarray) -> np.ndarray:
        """
        2階のテンソルのテンソル積

        Args:
            s (numpy.ndarray): 列ベクトル表示した2階のテンソル
            t (numpy.ndarray): 列ベクトル表示した2階のテンソル

        Returns:
            numpy.ndarray: 計算したテンソル積
        """
        ret: np.ndarray = np.zeros((6, 6))

        for i in range(6):
            for j in range(6):
                ret[i][j] = s[i][0] * t[j][0]

        return ret

    @staticmethod
    def norm2(t: np.ndarray) -> float:
        """
        2階のテンソルのノルム

        Args:
            t (numpy.ndarray): 列ベクトル表示した2階のテンソル

        Returns:
            float: 計算したノルム
        """
        return math.sqrt(Util.t2dott2(t, t))

    @staticmethod
    def trace2(t: np.ndarray) -> float:
        """
        2階のテンソルのトレース

        Args:
            t (numpy.ndarray): 列ベクトル表示した2階のテンソル

        Returns:
            float: 計算したトレース
        """
        return t[0][0] + t[1][0] + t[2][0]

    @staticmethod
    def i2() -> np.ndarray:
        """
        2階の恒等テンソル

        Returns:
            numpy.ndarray: 列ベクトル表示した恒等テンソル
        """
        ret: np.ndarray = np.zeros((6, 1))

        ret[0][0] = 1.
        ret[1][0] = 1.
        ret[2][0] = 1.

        return ret

    @staticmethod
    def i4s() -> np.ndarray:
        """
        4階の対称恒等テンソル

        Returns:
            numpy.ndarray: マトリクス表示した対称恒等テンソル
        """
        ret: np.ndarray = np.zeros((6, 6))

        for i in range(6):
            if i < 3:
                ret[i][i] = 1.
            else:
                ret[i][i] = 0.5

        return ret

    @staticmethod
    def i4d() -> np.ndarray:
        """
        4階の偏差射影テンソル

        Returns:
            numpy.ndarray: マトリクス表示した偏差射影テンソル
        """
        return Util.i4s() - Util.t2xt2(Util.i2(), Util.i2()) / 3.

    @staticmethod
    def trueStrain(strain: np.ndarray) -> np.ndarray:
        """
        工学ひずみの2倍されている成分を2で割る

        Args:
            t (numpy.ndarray): 工学ひずみ

        Returns:
            numpy.ndarray: 列ベクトル表示した真のひずみ
        """
        ret: np.ndarray = np.zeros((6, 1))

        for i in range(6):
            if i < 3:
                ret[i][0] = strain[i][0]
            else:
                ret[i][0] = strain[i][0] / 2.

        return ret