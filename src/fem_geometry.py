import numpy as np

class Gauss:
    """
    ガウス点に関するクラス
    """

    def __init__(self) -> None:
        """
        コンストラクタ
        """

        pass

class ShapeFunction:
    """
    形状関数に関するクラス
    """

    def __init__(self) -> None:
        """
        コンストラクタ
        """

        pass
    
class Geometry:
    """
    幾何に関するクラス
    """

    def __init__(self, x: list, y: list, z: list) -> None:
        """
        コンストラクタ

        Args:
            x (list): x座標のリスト
            y (list): y座標のリスト
            z (list): z座標のリスト
        """

        self._nNodes: int = len(x)
        self._x     : list = x
        self._y     : list = y
        self._z     : list = z

    def nNodes(self) -> int:
        """
        節点数

        Returns:
            int: 節点数
        """

        return self._nNodes

    def nodesIds(self) -> list:
        """
        要素の構成節点

        Returns:
            list: 要素の構成節点
        """

        return self._nodesIds

    def x(self) -> list:
        """
        x座標のリスト

        Returns:
            list: x座標のリスト
        """

        return self._x

    def y(self) -> list:
        """
        y座標のリスト

        Returns:
            list: y座標のリスト
        """

        return self._y

    def z(self) -> list:
        """
        z座標のリスト

        Returns:
            list: z座標のリスト
        """

        return self._z

class GaussTetrahedron(Gauss):
    """
    四面体のガウス点
    """

    def __init__(self) -> None:
        """
        コンストラクタ
        """

        super().__init__()

        self._nGauss: int = 1          # ガウス点数

        self._weight: list = [1. / 6.] # ガウス点毎の重み

        self._l1: list = [1. / 4.]     # ガウス点の座標(l1)
        self._l2: list = [self._l1[0]] # ガウス点の座標(l2)
        self._l3: list = [self._l1[0]] # ガウス点の座標(l3)
        self._l4: list = [self._l1[0]] # ガウス点の座標(l4)

    def nGauss(self) -> int:
        """
        ガウス点数

        Returns:
            int: ガウス点数
        """

        return self._nGauss

    def weight(self, i: int) -> float:
        """
        ガウス点の重み

        Args:
            i (int): ガウス点ID 

        Returns:
            float: ガウス点の重み
        """

        return self._weight[i]

    def g(self, i: int) -> float:
        """
        ガウス点の座標(ξ)

        Args:
            i (int): ガウス点ID 

        Returns:
            float: ガウス点の座標(ξ)
        """

        return self._l1[i]

    def e(self, i: int) -> float:
        """
        ガウス点の座標(η)

        Args:
            i (int): ガウス点ID 

        Returns:
            float: ガウス点の座標(η)
        """

        return self._l2[i]

    def ze(self, i: int) -> float:
        """
        ガウス点の座標(ζ)

        Args:
            i (int): ガウス点ID 

        Returns:
            float: ガウス点の座標(ζ)
        """

        return self._l3[i]

class ShapeFunctionTetrahedron(ShapeFunction):
    """
    四面体の形状関数
    """

    def __init__(self, gauss: GaussTetrahedron) -> None:
        """
        コンストラクタ

        Args:
            gauss (GaussTetrahedron): 四面体のガウス点
        """

        super().__init__()

        self._gauss: GaussTetrahedron = gauss                         # 四面体のガウス点

        self._dNdG: list = [np.zeros((4, 1))] * self._gauss.nGauss()  # ガウス点毎のd(N)/d(ξ)
        self._dNdE: list = [np.zeros((4, 1))] * self._gauss.nGauss()  # ガウス点毎のd(N)/d(η)
        self._dNdZe: list = [np.zeros((4, 1))] * self._gauss.nGauss() # ガウス点毎のd(N)/d(ζ)

        self._dNdL1: list = [0.] * 4                                   # d(N)/d(L1)
        self._dNdL2: list = [0.] * 4                                   # d(N)/d(L2)
        self._dNdL3: list = [0.] * 4                                   # d(N)/d(L3)
        self._dNdL4: list = [0.] * 4                                   # d(N)/d(L4)

        self._dNdL1[0] = 1.
        self._dNdL2[1] = 1.
        self._dNdL3[2] = 1.
        self._dNdL4[3] = 1.

        for i in range(gauss.nGauss()):
            for j in range(4):
                self._dNdG[i][j][0] = self._dNdL1[j] - self._dNdL4[j]
                self._dNdE[i][j][0] = self._dNdL2[j] - self._dNdL4[j]
                self._dNdZe[i][j][0] = self._dNdL3[j] - self._dNdL4[j]

    def dNdG(self, i: int) -> np.ndarray:
        """
        ガウス点のd(N)/d(ξ)

        Args:
            i (int): ガウス点ID

        Returns:
            numpy.ndarray: d(N)/d(ξ)
        """

        return self._dNdG[i]

    def dNdE(self, i: int) -> np.ndarray:
        """
        ガウス点のd(N)/d(η)

        Args:
            i (int): ガウス点ID

        Returns:
            numpy.ndarray: d(N)/d(η)
        """

        return self._dNdE[i]
        
    def dNdZe(self, i: int) -> np.ndarray:
        """
        ガウス点のd(N)/d(ζ)

        Args:
            i (int): ガウス点ID

        Returns:
            numpy.ndarray: d(N)/d(ζ)
        """

        return self._dNdZe[i]

class Tetrahedron(Geometry):
    """
    四面体に関するクラス
    """

    def __init__(self, x: list, y: list, z: list) -> None:
        """
        コンストラクタ

        Args:
            x (list): x座標のリスト
            y (list): y座標のリスト
            z (list): z座標のリスト
        """

        super().__init__(x, y, z)

        self._gauss: GaussTetrahedron = GaussTetrahedron()                                    # 四面体のガウス点
        self._shapeFunction: ShapeFunctionTetrahedron = ShapeFunctionTetrahedron(self._gauss) # 四面体の形状関数
        self._dNdX: list = [np.zeros((self._nNodes, 1))] * self._gauss.nGauss()               # ガウス点毎のd(N)/d(x)
        self._dNdY: list = [np.zeros((self._nNodes, 1))] * self._gauss.nGauss()               # ガウス点毎のd(N)/d(y)
        self._dNdZ: list = [np.zeros((self._nNodes, 1))] * self._gauss.nGauss()               # ガウス点毎のd(N)/d(z)
        self._detJac: list = []                                                               # ガウス点毎のヤコビアンの行列式
        self._bMatrix: list = []                                                              # ガウス点毎のひずみ-変位関係マトリクス(B)

        for i in range(self._gauss.nGauss()):
            jac: np.ndarray = self.jacobian(i)
            invJac: np.ndarray = np.linalg.inv(jac)

            self._detJac.append(np.linalg.det(jac))

            for j in range(self._nNodes):
                tmp1: np.ndarray = np.zeros((3, 1))
                tmp1[0][0] = self._shapeFunction.dNdG(i)[j][0]
                tmp1[1][0] = self._shapeFunction.dNdE(i)[j][0]
                tmp1[2][0] = self._shapeFunction.dNdZe(i)[j][0]

                tmp2: np.ndarray = invJac @ tmp1
                self._dNdX[i][j][0] = tmp2[0][0]
                self._dNdY[i][j][0] = tmp2[1][0]
                self._dNdZ[i][j][0] = tmp2[2][0]

        for i in range(self._gauss.nGauss()):
            bMat: np.ndarray = np.zeros((6, 3 * self._nNodes))

            for j in range(self._nNodes):
                bMat[0][3 * j] = self._dNdX[i][j][0]
                bMat[1][3 * j + 1] = self._dNdY[i][j][0]
                bMat[2][3 * j + 2] = self._dNdZ[i][j][0]
                bMat[3][3 * j] = self._dNdY[i][j][0]
                bMat[3][3 * j + 1] = self._dNdX[i][j][0]
                bMat[4][3 * j + 1] = self._dNdZ[i][j][0]
                bMat[4][3 * j + 2] = self._dNdY[i][j][0]
                bMat[5][3 * j] = self._dNdZ[i][j][0]
                bMat[5][3 * j + 2] = self._dNdX[i][j][0]

            self._bMatrix.append(bMat)

    def gauss(self) -> GaussTetrahedron:
        """
        四面体のガウス点

        Returns:
            GaussTetrahedron: 四面体のガウス点
        """

        return self._gauss

    def jacobian(self, i: int) -> np.ndarray:
        """
        ガウス点のヤコビアン

        Args:
            i (int): ガウス点ID

        Returns:
            numpy.ndarray: ガウス点のヤコビアン
        """

        jac = np.zeros((3, 3))

        for j in range(self._nNodes):
            jac[0][0] += self._shapeFunction.dNdG(i)[j][0] * self._x[j]
            jac[0][1] += self._shapeFunction.dNdG(i)[j][0] * self._y[j]
            jac[0][2] += self._shapeFunction.dNdG(i)[j][0] * self._z[j]

            jac[1][0] += self._shapeFunction.dNdE(i)[j][0] * self._x[j]
            jac[1][1] += self._shapeFunction.dNdE(i)[j][0] * self._y[j]
            jac[1][2] += self._shapeFunction.dNdE(i)[j][0] * self._z[j]

            jac[2][0] += self._shapeFunction.dNdZe(i)[j][0] * self._x[j]
            jac[2][1] += self._shapeFunction.dNdZe(i)[j][0] * self._y[j]
            jac[2][2] += self._shapeFunction.dNdZe(i)[j][0] * self._z[j]

        return jac

    def detJac(self, i: int) -> float:
        """
        ガウス点のヤコビアンの行列式

        Args:
            i (int): ガウス点ID

        Returns:
            numpy.ndarray: ガウス点のヤコビアンの行列式
        """

        return self._detJac[i]

    def bMatrix(self, i: int) -> np.ndarray:
        """
        ガウス点のひずみ-変位関係マトリクス(B)

        Args:
            i (int): ガウス点ID

        Returns:
            numpy.ndarray: ガウス点のひずみ-変位関係マトリクス(B)
        """

        return self._bMatrix[i]