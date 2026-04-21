from TwoGraph import *
from ClaudeHomologyCalculator import *
from GeminiHomologyCalculator import *

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python two_graph.py <file>")
        sys.exit(1)

    g = TwoGraph(sys.argv[1])

    print("Edge order", [e.label for e in g.edges])

    print("D matrices")

    print("d1", g.d_1.matrix)
    print("d2", g.d_2.matrix)


    result = claude_calculate_homology(g.d_1.matrix, g.d_2.matrix)
    print("H1:", result["H1"]["group"])
    print("H2:", result["H2"]["group"])

    exit()
    print(gemini_integer_homology(Matrix(g.d_1.matrix), Matrix(g.d_2.matrix)))
