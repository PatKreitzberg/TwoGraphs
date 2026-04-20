from TwoGraph import *
from MatrixMath import *



if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python two_graph.py <file>")
        sys.exit(1)

    g = TwoGraph(sys.argv[1])
