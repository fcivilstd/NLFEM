from fem_geometry import *
from fem_material import *

import numpy as np

class Node:
    """
    節点に関するクラス
    """

    def __init__(self, id: int, x: float, y: float, z: float) -> None:
        """
        コンストラクタ

        Args:
            id (int)  : 節点ID
            x  (float): x座標
            y  (float): y座標
            z  (float): z座標
        """

        self._id: int = id  # 節点ID
        self._x : float = x # x座標
        self._y : float = y # y座標
        self._z : float = z # z座標

    def id(self) -> int:
        """
        節点ID

        Returns:
            int: 節点ID
        """

        return self._id

    def x(self) -> float:
        """
        x座標

        Returns:
            float: x座標
        """

        return self._x
    
    def y(self) -> float:
        """
        y座標

        Returns:
            float: y座標
        """

        return self._y

    def z(self) -> float:
        """
        z座標

        Returns:
            float: z座標
        """

        return self._z

class Element:
    """
    要素に関するクラス
    """

    def __init__(self, id: int, nodesIds: list, geometry: Geometry, material: Material):
        """
        コンストラクタ

        Args:
            id       (int)     : 要素ID
            nodesIds (list)    : 要素の構成節点
            geometry (Geometry): 幾何情報
            material (Material): 材料情報
        """

        self._id        : int = id                                                    # 要素ID
        self._nodesIds  : list = nodesIds                                             # 要素の節点構成
        self._geometry  : Geometry = geometry                                         # 幾何情報
        self._material  : Material = material                                         # 材料情報
        self._nNodes    : int = len(self._nodesIds)                                   # 節点数
        self._kMatrix   : np.ndarray = np.zeros((3 * self._nNodes, 3 * self._nNodes)) # 要素剛性マトリクス(K)
        self._dltDisp   : np.ndarray = np.zeros((3 * self._nNodes, 1))                # 変位増分
        self._dltStrain : list = [np.zeros((6, 1))] * self._geometry.gauss().nGauss() # ひずみ増分
        self._stress    : list = [np.zeros((6, 1))] * self._geometry.gauss().nGauss() # 弾性ひずみ
        self._eStrain   : list = [np.zeros((6, 1))] * self._geometry.gauss().nGauss() # 塑性ひずみ
        self._pStrain   : list = [np.zeros((6, 1))] * self._geometry.gauss().nGauss() # 塑性ひずみ
        self._cumPStrain: list = [0.] * self._geometry.gauss().nGauss()               # 累積塑性ひずみ
        self._intF      : np.ndarray = np.zeros((3 * self._nNodes, 1))                # 内力ベクトル

    def id(self) -> int:
        """
        要素ID

        Returns:
            int: 要素ID
        """

        return self._id

    def kMatrix(self) -> np.ndarray:
        """
        要素剛性マトリクス(K)

        Returns:
            numpy.ndarray: 要素剛性マトリクス(K)
        """

        return self._kMatrix

    def calcKMatrix(self) -> None:
        """
        要素剛性マトリクス(K)の計算
        """

        self._kMatrix = np.zeros((3 * self._nNodes, 3 * self._nNodes))
        gauss: Gauss = self._geometry.gauss()
        for i in range(gauss.nGauss()):
            self._kMatrix += gauss.weight(i) * self._geometry.detJac(i) \
                            * (self._geometry.bMatrix(i).T @ self._material.dMatrix() @ self._geometry.bMatrix(i))       

    def calcStrain(self) -> None:
        """
        ひずみの計算
        """

        for i in range(self._geometry.gauss().nGauss()):
            self._dltStrain[i] = self._geometry.bMatrix(i) @ self._dltDisp

    def stateUpdate(self) -> None:
        """
        状態変数の更新
        """

        for i in range(self._geometry.gauss().nGauss()):
            self._stress[i], self._eStrain[i], self._pStrain[i], self._cumPStrain[i] \
                 = self._material.stateUpdate(self._dltStrain[i], self._eStrain[i], self._pStrain[i], self._cumPStrain[i])

    def calcInternalForces(self) -> None:
        """
        内力ベクトルの計算
        """

        self._intF = np.zeros((3 * self._nNodes, 1))
        gauss: Gauss = self._geometry.gauss()
        for i in range(gauss.nGauss()):
            self._intF += gauss.weight(i) * self._geometry.detJac(i) \
                        * (self._geometry.bMatrix(i).T @ self._stress[i])

    def setDltDisp(self, dltDisp: np.ndarray) -> None:
        """
        変位増分の設定

        Args:
            dltDisp (numpy.ndarray): 変位増分
        """

        self._dltDisp = np.zeros((3 * self._nNodes, 1))
        for i in range(self._nNodes):
            for j in range(3):
                self._dltDisp[3 * i + j][0] = dltDisp[3 * self._nodesIds[i] + j][0]

    def getKMatrix(self, kMat: np.ndarray) -> None:
        """
        要素剛性マトリクスから全体剛性マトリクスの組み立て

        Args:
            kMat (numpy.ndarray): 全体剛性マトリクス
        """

        for i in range(self._nNodes):
            for j in range(self._nNodes):
                for k in range(3):
                    for l in range(3):
                        kMat[3 * self._nodesIds[i] + k][3 * self._nodesIds[j] + l] += self._kMatrix[3 * i + k][3 * j + l]

    def getStrain(self, strain: np.ndarray) -> None:
        """
        ひずみの取得

        Args:
            strain (numpy.ndarray): ひずみ
        """

        gauss: Gauss = self._geometry.gauss()
        for i in range(gauss.nGauss()):
            for j in range(6):
                strain[6 * self._id + j][0] += self._eStrain[i][j][0] + self._pStrain[i][j][0]

        for j in range(6):
            strain[6 * self._id + j][0] /= gauss.nGauss()

    def getStress(self, stress: np.ndarray) -> None:
        """
        応力の取得

        Args:
            stress (numpy.ndarray): 応力
        """

        gauss: Gauss = self._geometry.gauss()
        for i in range(gauss.nGauss()):
            for j in range(6):
                stress[6 * self._id + j][0] += self._stress[i][j][0]
        
        for j in range(6):
            stress[6 * self._id + j][0] /= gauss.nGauss()

    def getInternalForces(self, intF: np.ndarray) -> None:
        """
        内力ベクトルの取得

        Args:
            intF (numpy.ndarray): 内力ベクトル
        """

        for i in range(self._nNodes):
            for j in range(3):
                intF[3 * self._nodesIds[i] + j][0] += self._intF[3 * i + j][0]
