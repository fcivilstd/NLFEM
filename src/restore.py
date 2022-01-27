from fem_model import *

def restore():
    case = "case1"

    # 解析モデルの復元
    model = NLModel.restore(case)

    # 描画
    model.visualize()

def main():
    restore()

if __name__ == '__main__':
    main()