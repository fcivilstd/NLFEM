from fem_model import *

def solve():
    case = "case1"

    # 解析モデルの構築
    model = NLModel()
    model.load(case)

    # 求解
    model.solve()

    # 解析結果の出力
    model.save()

    # 描画
    model.visualize()

def main():
    solve()

if __name__ == '__main__':
    main()