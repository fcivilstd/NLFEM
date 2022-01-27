from fem_util import *

import numpy as np
import math
import sys

class Material:
    """
    材料に関するクラス
    """

    def __init__(self, young: float, poisson: float) -> None:
        """
        コンストラクタ

        Args:
            young   (float): ヤング率(E)
            poisson (float): ポアソン比(ν)
        """

        self._isYield: bool = False                 # 降伏しているかのフラグ
        self._young  : float = young                # ヤング率(E)
        self._poisson: float = poisson              # ポアソン比(ν)
        self._dMatrix: np.ndarray = self.dEMatrix() # コンシステント接線(D)

    def young(self) -> float:
        """
        ヤング率(E)

        Returns:
            float: メンバ変数に格納されているヤング率
        """

        return self._young

    def poisson(self) -> float:
        """
        ポアソン比(ν)

        Returns:
            float: メンバ変数に格納されているポアソン比
        """

        return self._poisson

    def g(self) -> float:
        """
        せん断弾性係数(G)

        Returns:
            float: 計算したせん断弾性係数
        """

        return self._young / 2. / (1. + self._poisson)

    def k(self) -> float:
        """
        体積弾性係数(K)

        Returns:
            float: 計算した体積弾性係数
        """

        return self._young / 3. / (1. - 2. * self._poisson)

    def dMatrix(self) -> np.ndarray:
        """
        コンシステント接線(D)

        Returns:
            numpy.ndarray: メンバ変数に格納されているコンシステント接線
        """
        
        return self._dMatrix

    def dEMatrix(self) -> np.ndarray:
        """
        弾性コンシステント接線(De)

        Returns:
            numpy.ndarray: 計算した弾性コンシステント接線
        """

        return 2. * self.g() * Util.i4s() + (self.k() - 2. / 3. * self.g()) * Util.t2xt2(Util.i2(), Util.i2())

class VonMises(Material):
    """
    Von Misesに関する変数や関数を定義したクラス
    """

    def __init__(self, young: float, poisson: float) -> None:
        """
        コンストラクタ

        Args:
            young   (float): ヤング率(E)
            poisson (float): ポアソン比(ν)
        """

        super().__init__(young, poisson)

        self._nItr: int = 10            # 反復回数の上限
        self._tol: float = 0.01         # 収束比の閾値
        self._yieldStress: float = 295. # 初期降伏応力

    def yieldStress(self, cumPStrain: float) -> float:
        """
        降伏応力

        Args:
            cumPStrain (float): 累積塑性ひずみ

        Returns:
            float: 計算した降伏応力
        """

        # ヤング率の100分の1で硬化するモデル
        return self._yieldStress + self._young / 100. * cumPStrain

    def hardening(self, cumPStrain) -> float:
        """
        硬化係数

        Args:
            cumPStrain (float): 累積塑性ひずみ

        Returns:
            float: 計算した硬化係数
        """

        # ヤング率の100分の1で硬化するモデル
        return self._young / 100.

    def stateUpdate(
            self,
            dltStrain: np.ndarray,
            prevEStrain: np.ndarray,
            prevPStrain: np.ndarray,
            prevCumPStrain: float,
        ) -> tuple:
        """
        状態変数の更新

        Args:
            dltStrain      (numpy.ndarray): ひずみ増分
            prevEStrain    (numpy.ndarray): 直前の弾性ひずみ
            prevPStrain    (numpy.ndarray): 直前の塑性ひずみ
            prevCumPStrain (float)        : 直前の累積塑性ひずみ

        Returns:
            numpy.ndarray: 現在の応力
            numpy.ndarray: 現在の弾性ひずみ
            numpy.ndarray: 現在の塑性ひずみ
            float        : 現在の累積塑性ひずみ
        """

        y: float = self.yieldStress(prevCumPStrain)                   # 降伏応力

        eStrainT: np.ndarray = prevEStrain + dltStrain                # 試行弾性ひずみ
        evStrainT: float = Util.trace2(eStrainT)                      # 試行弾性体積ひずみ
        edStrainT: np.ndarray = eStrainT - evStrainT / 3. * Util.i2() # 試行弾性偏差ひずみ
        pT: float = self.k() * evStrainT                              # 試行静水圧
        sT: np.ndarray = 2. * self.g() * Util.trueStrain(edStrainT)   # 試行偏差応力
        qT: float = math.sqrt(3. / 2. * Util.t2dott2(sT, sT))         # 試行弾性相当応力

        eStrain: np.ndarray = eStrainT                                # 弾性ひずみ
        p: float = pT                                                 # 静水圧
        s: np.ndarray = sT                                            # 偏差応力
        pStrain: np.ndarray = prevPStrain                             # 塑性ひずみ
        cumPStrain: float = prevCumPStrain                            # 累積塑性ひずみ
        stress: np.ndarray = s + p * Util.i2()                        # 応力

        # 降伏判定
        if qT - y > 0.:
            self._isYield = True
            dltG: float = 0.                                          # 塑性乗数の増分
            phi: float = qT - y                                       # 降伏関数

            # ニュートン・ラプソン法
            for itr in range(self._nItr + 1):

                if itr == self._nItr:
                    sys.exit("not converged")

                h: float = self.hardening(prevCumPStrain + dltG)      # 硬化係数
                d: float = -3. * self.g() - h                         # d(phi)/d(dltG)
                dltG -= phi / d

                phi = qT - 3. * self.g() * dltG - self.yieldStress(prevCumPStrain + dltG)

                # 収束判定
                if abs(phi) <= self._tol:
                    s = (1. - dltG * 3. * self.g() / qT) * sT
                    eStrain = s / 2. / self.g() + evStrainT / 3. * Util.i2()
                    pStrain = prevPStrain + dltG * math.sqrt(3. / 2.) * s / Util.norm2(s)
                    cumPStrain = prevCumPStrain + dltG
                    stress = s + p * Util.i2()
                    h = self.hardening(prevCumPStrain + dltG)
                    n: np.ndarray = sT / Util.norm2(sT)               # 単位流れベクトル

                    # 弾塑性コンシステント接線を更新
                    self._dMatrix = 2. * self.g() * (1. - dltG * 3. * self.g() / qT) * Util.i4d() \
                                    + 6. * self.g() ** 2 * (dltG / qT - 1. / (3. * self.g() + h)) * Util.t2xt2(n, n) \
                                    + self.k() * Util.t2xt2(Util.i2(), Util.i2())

                    break
        else:
            # 塑性から弾性に変わった場合弾性コンシステント接線に更新
            if self._isYield == True:
                self._isYield = False
                self._dMatrix = self.dEMatrix()

        return stress, eStrain, pStrain, cumPStrain
