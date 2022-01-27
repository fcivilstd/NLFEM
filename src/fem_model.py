from fem_common     import *
from fem_geometry   import *
from fem_material   import *
from fem_visualizer import *

import numpy as np
import pickle

class Model:
    """
    解析モデルに関するクラス
    """

    def __init__(self) -> None:
        """
        コンストラクタ
        """

        self._filename   : str = ""                   # 入力ファイル名

        self._nNodes     : int = 0                    # 節点数
        self._nElms      : int = 0                    # 要素数
        self._nPLoad     : int = 0                    # 節点荷重数
        self._nFixedNodes: int = 0                    # 拘束節点数
        self._nodes      : list = []                  # 節点情報のリスト
        self._elms       : list = []                  # 要素情報のリスト
        self._nodesIds   : list = []                  # 要素の構成節点
        self._exF        : np.ndarray = np.zeros((1)) # 外力ベクトル
        self._fixedNodes : list = []                  # 拘束節点のリスト

class NLModel(Model):
    """
    非線形モデル
    """

    def __init__(self) -> None:
        """
        コンストラクタ
        """

        super().__init__()

        self._nInc     : int = 100      # インクリメント数
        self._nItr     : int = 10       # 最大イテレーション数
        self._tol      : float = 0.1    # 収束判定の閾値
        self._lastInc  : int = 0        # 解が得られた最後のインクリメント

        self._exFInc   : list = []      # インクリメント毎の外力ベクトル
        self._intFInc  : list = []      # インクリメント毎の内力ベクトル
        self._dispInc  : list = []      # インクリメント毎の変位
        self._strainInc: list = []      # インクリメント毎のひずみ
        self._stressInc: list = []      # インクリメント毎の応力

    def load(self, filename: str) -> None:
        """
        入力ファイルの読み込み

        Args:
            filename (str): 入力ファイルの名前
        """

        young  : float = 200000
        poisson: float = 0.3

        self._filename = filename

        with open("../in/" + filename + ".txt") as f:
            self._nNodes, self._nElms = map(int, f.readline().split())

            for nodeId in range(self._nNodes):
                x: float
                y: float
                z: float
                x, y, z = map(float, f.readline().split())
                self._nodes.append(Node(nodeId, x, y, z))

            for elmId in range(self._nElms):
                nodesIds: list = list(map(lambda i: int(i) - 1, f.readline().split()))
                nodes: list = []
                x: list = []
                y: list = []
                z: list = []
                for nodeId in nodesIds:
                    nodes.append(self._nodes[nodeId])
                    x.append(self._nodes[nodeId].x())
                    y.append(self._nodes[nodeId].y())
                    z.append(self._nodes[nodeId].z())
                self._elms.append(
                    Element(
                        elmId,
                        nodesIds,
                        Tetrahedron(x, y, z),
                        VonMises(young, poisson)
                    )
                )
                self._nodesIds.append(nodesIds)

            self._nPLoad = int(f.readline())

            self._exF = np.zeros((3 * self._nNodes, 1))
            for i in range(self._nPLoad):
                line: list = f.readline().split()

                nodeId: int = int(line[0]) - 1

                fx: float
                fy: float
                fz: float
                fx, fy, fz = map(float, line[1:4])
                self._exF[3 * nodeId][0] = fx
                self._exF[3 * nodeId + 1][0] = fy
                self._exF[3 * nodeId + 2][0] = fz

            self._nFixedNodes = int(f.readline())
            for i in range(self._nFixedNodes):
                self._fixedNodes.append(int(f.readline()) - 1)

    def setBoundaryCondition(self, kMat: np.ndarray, r: np.ndarray) -> tuple:
        """
        境界条件の設定
        
        Args:
            kMat (np.ndarray): 全体剛性マトリクス(K)
            r    (np.ndarray): 残差力ベクトル(r)

        Returns:
            numpy.ndarray: 境界条件を設定した全体剛性マトリクス(K)
            numpy.ndarray: 境界条件を設定した残差力ベクトル(r)
        """

        _kMatBC: np.ndarray = np.copy(kMat)
        _rBC   : np.ndarray = np.copy(r)

        for i in self._fixedNodes:
            for j in range(3):
                _rBC[3 * i + j] = 0.

                _kMatBC[:, 3 * i + j] = 0.
                _kMatBC[3 * i + j, :] = 0.
                _kMatBC[3 * i + j, 3 * i + j] = 1.

        return _kMatBC, _rBC

    def solve(self) -> None:
        """
        インクリメント毎の解を求める
        """

        self._exFInc = [np.zeros((3 * self._nNodes, 1))]
        self._intFInc = [np.zeros((3 * self._nNodes, 1))] * (self._nInc + 1)
        self._dispInc = [np.zeros((3 * self._nNodes, 1))] * (self._nInc + 1)
        self._strainInc = [np.zeros((6 * self._nElms, 1))] * (self._nInc + 1)
        self._stressInc = [np.zeros((6 * self._nElms, 1))] * (self._nInc + 1)

        # 各インクリメントにおける外力ベクトルを計算
        for i in range(self._nInc):
            self._exFInc.append(self._exF * float(i + 1) / float(self._nInc))

        for i in range(self._nInc):
            inc: int = i + 1

            print("increment: " + str(inc))

            r       : np.ndarray = self._exFInc[inc] - self._intFInc[inc - 1]
            prevIntF: np.ndarray = np.copy(self._intFInc[inc - 1])
            disp    : np.ndarray = np.copy(self._dispInc[inc - 1])

            kMat: np.ndarray = np.zeros((3 * self._nNodes, 3 * self._nNodes))
            for elm in self._elms:
                elm.calcKMatrix()
                elm.getKMatrix(kMat)

            kMatBC: np.ndarray
            rBC   : np.ndarray
            kMatBC, rBC = self.setBoundaryCondition(kMat, r)

            # ニュートン・ラプソン法
            for itr in range(self._nItr + 1):

                if itr == self._nItr:
                    print("not converged")
                    return

                print("iteration: " + str(itr))

                dltDisp: np.ndarray = np.linalg.solve(kMatBC, rBC)
                disp += dltDisp

                intF: np.ndarray = np.zeros((3 * self._nNodes, 1))
                for elm in self._elms:
                    elm.setDltDisp(dltDisp)
                    elm.calcStrain()
                    elm.stateUpdate()
                    elm.calcInternalForces()
                    elm.getInternalForces(intF)

                r -= intF - prevIntF
                prevIntF = intF

                kMat = np.zeros((3 * self._nNodes, 3 * self._nNodes))
                for elm in self._elms:
                    elm.calcKMatrix()
                    elm.getKMatrix(kMat)

                kMatBC, rBC = self.setBoundaryCondition(kMat, r)

                print("convergence ratio: " + str(np.linalg.norm(dltDisp, ord = np.inf) / np.linalg.norm(disp - self._dispInc[inc - 1], ord = np.inf)))

                # 収束判定
                if np.linalg.norm(dltDisp, ord = np.inf) / np.linalg.norm(disp - self._dispInc[inc - 1], ord = np.inf) <= self._tol:
                    strain: np.ndarray = np.zeros((6 * self._nElms, 1))
                    stress: np.ndarray = np.zeros((6 * self._nElms, 1))
                    for elm in self._elms:
                        elm.getStrain(strain)
                        elm.getStress(stress)

                    self._intFInc[inc] = intF
                    self._dispInc[inc] = disp
                    self._strainInc[inc] = strain
                    self._stressInc[inc] = stress

                    self._lastInc = inc
                    break

    def save(self) -> None:
        """
        解析結果をバイナリデータで保存する
        """

        with open("../out/" + self._filename + ".pickle", mode = "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def restore(filename: str) -> "NLModel":
        """
        解析結果をバイナリデータから復元する

        Args:
            filename (str): 解析結果が保存してあるファイル名
        """

        with open("../out/" + filename + ".pickle", mode = "rb") as f:
            return pickle.load(f)

    def principalStrain(self, inc: int = -1) -> list:
        """
        最大主ひずみを計算して取得する
        
        Args:
            inc (int): 最大主ひずみを計算するインクリメント

        Returns:
            list: 最大主ひずみの要素スカラー値
        """

        if inc == -1:
            inc = self._lastInc

        ps: list = []
        for elmId in range(self._nElms):
            a: np.ndarray = np.zeros((3, 3))
            a[0][0] = self._strainInc[inc][6 * elmId + 0][0]
            a[1][1] = self._strainInc[inc][6 * elmId + 1][0]
            a[2][2] = self._strainInc[inc][6 * elmId + 2][0]
            a[0][1] = self._strainInc[inc][6 * elmId + 3][0]
            a[1][0] = self._strainInc[inc][6 * elmId + 3][0]
            a[1][2] = self._strainInc[inc][6 * elmId + 4][0]
            a[2][1] = self._strainInc[inc][6 * elmId + 4][0]
            a[0][2] = self._strainInc[inc][6 * elmId + 5][0]
            a[2][0] = self._strainInc[inc][6 * elmId + 5][0]

            eig: np.array
            v  : np.array
            eig, v = np.linalg.eig(a)
            ps.append(max(eig))

        return ps

    def principalStress(self, inc: int = -1) -> list:
        """
        最大主応力を計算して取得する
        
        Args:
            inc (int): 最大主応力を計算するインクリメント

        Returns:
            list: 最大主応力の要素スカラー値
        """

        if inc == -1:
            inc = self._lastInc

        ps: list = []
        for elmId in range(self._nElms):
            a: np.ndarray = np.zeros((3, 3))
            a[0][0] = self._stressInc[inc][6 * elmId + 0][0]
            a[1][1] = self._stressInc[inc][6 * elmId + 1][0]
            a[2][2] = self._stressInc[inc][6 * elmId + 2][0]
            a[0][1] = self._stressInc[inc][6 * elmId + 3][0]
            a[1][0] = self._stressInc[inc][6 * elmId + 3][0]
            a[1][2] = self._stressInc[inc][6 * elmId + 4][0]
            a[2][1] = self._stressInc[inc][6 * elmId + 4][0]
            a[0][2] = self._stressInc[inc][6 * elmId + 5][0]
            a[2][0] = self._stressInc[inc][6 * elmId + 5][0]

            eig: np.array
            v  : np.array
            eig, v = np.linalg.eig(a)
            ps.append(max(eig))

        return ps

    def meaningScalarsToNode(self, elmScalars: list) -> list:
        """
        要素スカラー値を節点で平均化

        Args:
            elmScalars (list): 要素スカラー値

        Returns:
            list: 節点で平均化した要素スカラー値
        """

        nodeScalars: list = [[] for i in range(self._nNodes)]
        mNodeScalars: list = []
        for i in range(len(elmScalars)):
            for j in range(len(self._nodesIds[i])):
                nodeScalars[self._nodesIds[i][j]].append(elmScalars[i])

        for nodeScalar in nodeScalars:
            m = 0
            for scalar in nodeScalar:
                m += scalar
            mNodeScalars.append(m / len(nodeScalar))

        return mNodeScalars

    def visualize(self, scalars: list = [], inc: int = -1) -> None:
        """
        解析結果の描画
        
        Args:
            inc (int): 表示するインクリメント
        """

        if inc == -1:
            inc = self._lastInc

        if scalars == []:
            scalars = self.meaningScalarsToNode(self.principalStress())

        x: list = []
        y: list = []
        z: list = []
        dispX: list = [0] * self._nNodes
        dispY: list = [0] * self._nNodes
        dispZ: list = [0] * self._nNodes
        for node in self._nodes:
            nodeId = node.id()
            x.append(node.x())
            y.append(node.y())
            z.append(node.z())
            dispX[nodeId] = self._dispInc[inc][3 * nodeId][0]
            dispY[nodeId] = self._dispInc[inc][3 * nodeId + 1][0]
            dispZ[nodeId] = self._dispInc[inc][3 * nodeId + 2][0]

        vis: Visualizer = Visualizer(
            self._nNodes, self._nElms,
            x, y, z,
            self._nodesIds,
            dispX, dispY, dispZ
        )

        print("\nvisualized increment: " + str(inc))

        vis.visualize(scalars)