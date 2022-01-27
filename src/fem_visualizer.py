import open3d as o3d
import colorsys

class Visualizer:
    """
    ビジュアライザー
    """

    def __init__(
        self,
        nNodes: int, nElms: int,
        x: list, y: list, z: list,
        nodesIds: list,
        dispX: list, dispY: list, dispZ: list
    ) -> None:
        """
        コンストラクタ

        Args:
            nNodes   (int) : 節点数
            nElms    (int) : 要素数
            x        (list): 初期のx座標
            y        (list): 初期のy座標
            z        (list): 初期のz座標
            nodesIds (list): 各要素の構成節点ID(0-indexed)
            dispX    (list): 節点変位x
            dispY    (list): 節点変位y
            dispZ    (list): 節点変位z
        """

        self._nNodes: int = nNodes          # 節点数
        self._nElms: int = nElms            # 要素数
        self._x: float = x                  # 初期のx座標
        self._y: float = y                  # 初期のy座標
        self._z: float = z                  # 初期のz座標
        self._nodesIds: list = nodesIds     # 各要素の構成節点ID(0-indexed)
        self._dispX: list = dispX           # 節点変位x
        self._dispY: list = dispY           # 節点変位y
        self._dispZ: list = dispZ           # 節点変位z
        self._nLayers: int = 0              # ヒートマップの層数
        self._heatmapColors: list = []      # ヒートマップで使用する色
        self._mapedHeatmapColors: list = [] # 節点の色
        self._originalMesh = None           # 初期のメッシュ(Open3D)
        self._lineSet = None                # 初期のメッシュ(ワイヤーフレーム, Open3D)
        self._axis = None                   # 座標軸のメッシュ(Open3D)

    def getTriangleForMesh(self) -> list:
        """
        三角形頂点の組み合わせを取得

        Returns:
            list: 三角頂点の組み合わせ
        """

        triangles: list = []
        for ids in self._nodesIds:
            triangles.append([ids[0], ids[1], ids[2]])
            triangles.append([ids[0], ids[3], ids[1]])
            triangles.append([ids[0], ids[2], ids[3]])
            triangles.append([ids[1], ids[3], ids[2]])

        return triangles

    def getOriginalVerticesForMesh(self) -> list:
        """
        変形前の頂点座標の取得

        Returns:
            list: 変形前の頂点座標
        """

        originalVertices: list = []
        for i in range(self._nNodes):
            originalVertices.append([self._x[i], self._y[i], self._z[i]])

        return originalVertices

    def getDeformedVerticesForMesh(self) -> list:
        """
        変形後の頂点座標の取得

        Returns:
            list: 変形後の頂点座標
        """

        deformedVertices: list = []
        for i in range(self._nNodes):
            deformedVertices.append([
                                    self._x[i] + self._dispX[i],
                                    self._y[i] + self._dispY[i],
                                    self._z[i] + self._dispZ[i]
                                ])

        return deformedVertices

    def genRGBHeatmap(self, s: float, t: float) -> list:
        """
        ヒートマップ用のRGB配列を生成

        Args:
            s (float): hsv色空間のhにおける始点(0~360)
            t (float): hsv色空間のhにおける終点(0~360)

        Returns:
            list: 各層のRGB値(0~1)
        """

        rgb: list = []
        h  : list = []
        for i in range(self._nLayers + 1):
            tmp = s + abs(t - s) / self._nLayers * i
            if tmp > 360.:
                tmp -= 360.

            h.append(tmp / 360.)

        for _h in h:
            rgb.append(colorsys.hsv_to_rgb(_h, 1, 1))

        return rgb

    def colorMapping(self, heatmapColors: list, scalars: list) -> list:
        """
        節点にヒートマップで使用する色を割り当てる

        Args:
            heatmapColors (list): 使用する色のリスト
            scalars       (list): 節点スカラー値

        Returns:
            list: 割り当てられた節点の色
        """

        minScalar: float = min(scalars)
        maxScalar: float = max(scalars)

        interval  : float = (maxScalar - minScalar) / self._nLayers
        valueLayer: list = []

        for i in range(self._nLayers + 1):
            valueLayer.append(minScalar + interval * i)

        print("\n-----Heatmap Value-----\n")
        for i in range(len(valueLayer)):
            if i == 0:
                print("\033[33m" + str(valueLayer[len(valueLayer) - 1]) + "\033[0m")
            elif i == len(valueLayer) - 1:
                print("\033[34m" + str(valueLayer[0]) + "\033[0m")
            else:
                print(str(valueLayer[len(valueLayer) - i - 1]))

        mapedHeatmapColors: list = [[0, 0, 0] for i in range(len(scalars))]

        for i in range(len(scalars)):
            scalar = scalars[i]
            for j in range(self._nLayers):
                if scalar >= valueLayer[j] and scalar <= valueLayer[j + 1]:
                    mapedHeatmapColors[i] = heatmapColors[j]

        return mapedHeatmapColors

    def visualize(self, scalar: list, nLayers: int = 10) -> None:
        """
        節点にヒートマップで使用する色を割り当てる

        Args:
            scalars (list): 節点スカラー値
            nLayers (int) : ヒートマップの層数
        """

        self._nLayers = nLayers

        self._heatmapColors = self.genRGBHeatmap(240, 60)

        # ヒートマップ表示したいスカラー値の設定
        self._mapedHeatmapColors = self.colorMapping(self._heatmapColors, scalar)
    
        # 変形前のメッシュ
        self._originalMesh = o3d.geometry.TriangleMesh()
        self._originalMesh.vertices = o3d.utility.Vector3dVector(self.getOriginalVerticesForMesh())
        self._originalMesh.triangles = o3d.utility.Vector3iVector(self.getTriangleForMesh())
        self._originalMesh.compute_vertex_normals()

        # 変形前のメッシュのエッジ
        self._lineSet = o3d.geometry.LineSet.create_from_triangle_mesh(self._originalMesh)
        self._lineSet.colors = o3d.utility.Vector3dVector(
                                [[1, 0, 0] for i in range(len(self._lineSet.lines))]
                            )

        # 変形後のメッシュ
        self._deformedMesh = o3d.geometry.TriangleMesh()
        self._deformedMesh.vertices = o3d.utility.Vector3dVector(self.getDeformedVerticesForMesh())
        self._deformedMesh.triangles = o3d.utility.Vector3iVector(self.getTriangleForMesh())
        self._deformedMesh.vertex_colors = o3d.utility.Vector3dVector(self._mapedHeatmapColors)
        self._deformedMesh.compute_vertex_normals()

        # 座標軸のメッシュ
        self._axis = o3d.geometry.TriangleMesh.create_coordinate_frame(
                                                                size = 10,
                                                                origin = [-10, -10, -10]
                                                            )

        # 描画
        o3d.visualization.draw_geometries(
                                        [self._lineSet, self._deformedMesh, self._axis],
                                        "Result",
                                        1920,
                                        1080,
                                        0,
                                        0
                                    )