"""
TwoGraph: parse and compute homology of a 2-dimensional directed graph (2-complex).

File format
-----------
Lines of the form:

    <label> <source> <target>        — defines a directed edge (1-cell)
    <a> <b> = <c> <d>                — commuting square: a·b = c·d (2-cell)

where every token is a string label.

A "commuting square"  a·b = c·d  means the composite of edges a then b equals
the composite of edges c then d.  This gives a 2-cell whose boundary is

    ∂₂(square) = a + b − c − d   (in C₁, with integer coefficients)

Homology
--------
The chain complex is

    C₂  --∂₂-->  C₁  --∂₁-->  C₀

    H₁ = ker ∂₁ / im ∂₂
    H₂ = ker ∂₂            (there are no 3-cells)

Groups are displayed as  Z{<generators>} / Z{<generators>}  where generators
are named linear combinations of edges (for H₁) or squares (for H₂).

  • The kernel generators come from sympy's integer nullspace (DomainMatrix over ZZ).
  • The image generators are exactly the columns of ∂₂ — each is the boundary
    of one commuting square, written as a signed sum of edge labels.
  • The rank / torsion structure is read from sympy's smith_normal_form.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Union

from sympy import ZZ, Matrix
from sympy.matrices.normalforms import smith_normal_form as _sympy_snf
from sympy.polys.matrices import DomainMatrix


# ---------------------------------------------------------------------------
# Integer linear-algebra helpers (all backed by sympy)
# ---------------------------------------------------------------------------

def _int_nullspace(mat: list[list[int]]) -> list[list[int]]:
    """
    Return a ℤ-basis for the (right) kernel of an integer matrix,
    as a list of integer column vectors (each a Python list of ints).

    Uses sympy's DomainMatrix.nullspace() over ZZ, which produces a minimal
    basis without any denominators.
    """
    if not mat or not mat[0]:
        # Degenerate: zero rows → kernel is all of ℤ^ncols
        ncols = len(mat[0]) if mat else 0
        return [[(1 if i == j else 0) for i in range(ncols)] for j in range(ncols)]
    M = Matrix(mat)
    nrows, ncols = M.shape
    if ncols == 0:
        return []
    DM = DomainMatrix.from_Matrix(M).convert_to(ZZ)
    K = DM.nullspace().to_Matrix()   # rows of K are kernel basis vectors
    if K.rows == 0:
        return []
    return [list(map(int, K.row(i))) for i in range(K.rows)]


def _snf_diagonal(mat: list[list[int]]) -> list[int]:
    """
    Return the diagonal of the Smith normal form of an integer matrix,
    as a list of non-negative integers (including trailing zeros so that
    the caller can read off the rank).

    Uses sympy.matrices.normalforms.smith_normal_form.
    """
    if not mat or not mat[0]:
        return []
    M = Matrix(mat)
    D = _sympy_snf(M, domain=ZZ)
    nrows, ncols = D.shape
    return [int(D[i, i]) for i in range(min(nrows, ncols))]


def _vec_to_label(vec: list[int], names: list[str]) -> str:
    """
    Format an integer coefficient vector as a signed linear combination
    of generator names.

        [1, -2, 0, 3], ['a','b','c','d']  →  "a - 2b + 3d"
    """
    terms = []
    for coeff, name in zip(vec, names):
        if coeff == 0:
            continue
        if not terms:
            if coeff == 1:
                terms.append(name)
            elif coeff == -1:
                terms.append(f"-{name}")
            else:
                terms.append(f"{coeff}{name}")
        else:
            if coeff == 1:
                terms.append(f"+ {name}")
            elif coeff == -1:
                terms.append(f"- {name}")
            elif coeff > 0:
                terms.append(f"+ {coeff}{name}")
            else:
                terms.append(f"- {abs(coeff)}{name}")
    return " ".join(terms) if terms else "0"


def _format_group(
    ker_vecs: list[list[int]],
    ker_names: list[str],
    im_vecs: list[list[int]],
    im_names: list[str],
) -> str:
    """
    Format  Z{ker generators} / Z{im generators}.

    Special cases:
      ker trivial  →  "0"
      im trivial   →  "Z{ … }"   (no denominator)
    """
    if not ker_vecs:
        return "0"
    ker_strs = [_vec_to_label(v, ker_names) for v in ker_vecs]
    num = "Z{ " + ", ".join(ker_strs) + " }"
    if not im_vecs:
        return num
    im_strs = [_vec_to_label(v, im_names) for v in im_vecs]
    den = "Z{ " + ", ".join(im_strs) + " }"
    return f"{num} / {den}"


# ---------------------------------------------------------------------------
# TwoGraph
# ---------------------------------------------------------------------------

class TwoGraph:
    """
    A directed 2-graph (2-dimensional cell complex) built from:
      - 0-cells  (vertices)
      - 1-cells  (directed labelled edges)
      - 2-cells  (commuting squares, each asserting a·b = c·d)

    Attributes
    ----------
    vertices  : list[str]
    edges     : list[dict]   — keys: 'label', 'source', 'target'
    squares   : list[dict]   — keys: 'label', 'left', 'right'
                               (left, right are 2-tuples of edge labels)
    """

    def __init__(self) -> None:
        self.vertices: list[str] = []
        self.edges: list[dict] = []
        self.squares: list[dict] = []
        self._vertex_index: dict[str, int] = {}
        self._edge_index: dict[str, int] = {}
        self._square_index: dict[str, int] = {}

    # ------------------------------------------------------------------
    # Parsing
    # ------------------------------------------------------------------

    @classmethod
    def from_file(cls, path: Union[str, Path]) -> "TwoGraph":
        """Parse a text file and return a TwoGraph."""
        g = cls()
        with open(path) as f:
            for lineno, raw in enumerate(f, 1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                g._parse_line(line, lineno)
        return g

    @classmethod
    def from_string(cls, text: str) -> "TwoGraph":
        """Parse a multi-line string and return a TwoGraph."""
        g = cls()
        for lineno, raw in enumerate(text.splitlines(), 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            g._parse_line(line, lineno)
        return g

    def _parse_line(self, line: str, lineno: int = 0) -> None:
        if "=" in line:
            self._parse_square(line, lineno)
        else:
            self._parse_edge(line, lineno)

    def _parse_edge(self, line: str, lineno: int) -> None:
        parts = line.split()
        if len(parts) != 3:
            raise ValueError(
                f"Line {lineno}: edge line must be '<label> <source> <target>', "
                f"got: {line!r}"
            )
        label, source, target = parts
        self._add_vertex(source)
        self._add_vertex(target)
        self._add_edge(label, source, target, lineno)

    def _parse_square(self, line: str, lineno: int) -> None:
        m = re.fullmatch(r"(\S+)\s+(\S+)\s*=\s*(\S+)\s+(\S+)", line)
        if not m:
            raise ValueError(
                f"Line {lineno}: square line must be '<a> <b> = <c> <d>', "
                f"got: {line!r}"
            )
        a, b, c, d = m.group(1), m.group(2), m.group(3), m.group(4)
        label = f"({a}·{b}={c}·{d})"
        self._add_square(label, (a, b), (c, d), lineno)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _add_vertex(self, v: str) -> None:
        if v not in self._vertex_index:
            self._vertex_index[v] = len(self.vertices)
            self.vertices.append(v)

    def _add_edge(self, label: str, source: str, target: str, lineno: int = 0) -> None:
        if label in self._edge_index:
            existing = self.edges[self._edge_index[label]]
            if existing["source"] != source or existing["target"] != target:
                raise ValueError(
                    f"Line {lineno}: edge {label!r} redefined with different endpoints."
                )
            return
        self._edge_index[label] = len(self.edges)
        self.edges.append({"label": label, "source": source, "target": target})

    def _add_square(
        self,
        label: str,
        left: tuple[str, str],
        right: tuple[str, str],
        lineno: int = 0,
    ) -> None:
        for el in left + right:
            if el not in self._edge_index:
                raise ValueError(
                    f"Line {lineno}: square references unknown edge {el!r}. "
                    f"Declare edges before squares."
                )
        self._square_index[label] = len(self.squares)
        self.squares.append({"label": label, "left": left, "right": right})

    # ------------------------------------------------------------------
    # Boundary maps
    # ------------------------------------------------------------------

    def _boundary1_matrix(self) -> list[list[int]]:
        """
        ∂₁ : C₁ → C₀.  Shape (|V|, |E|).
        Column for edge e: −1 at row source(e), +1 at row target(e).
        """
        nv, ne = len(self.vertices), len(self.edges)
        mat = [[0] * ne for _ in range(nv)]
        for j, edge in enumerate(self.edges):
            mat[self._vertex_index[edge["source"]]][j] -= 1
            mat[self._vertex_index[edge["target"]]][j] += 1
        return mat

    def _boundary2_matrix(self) -> list[list[int]]:
        """
        ∂₂ : C₂ → C₁.  Shape (|E|, |squares|).
        Column k for square a·b = c·d: +a + b − c − d in edge coordinates.
        """
        ne, ns = len(self.edges), len(self.squares)
        mat = [[0] * ns for _ in range(ne)]
        for k, sq in enumerate(self.squares):
            a, b = sq["left"]
            c, d = sq["right"]
            mat[self._edge_index[a]][k] += 1
            mat[self._edge_index[b]][k] += 1
            mat[self._edge_index[c]][k] -= 1
            mat[self._edge_index[d]][k] -= 1
        return mat

    # ------------------------------------------------------------------
    # Homology
    # ------------------------------------------------------------------

    def first_homology(self) -> str:
        """
        H₁ = ker ∂₁ / im ∂₂.

        Numerator generators: integer nullspace of ∂₁, each expressed as a
        signed sum of edge labels (these are the independent 1-cycles).

        Denominator generators: columns of ∂₂, each expressed as a signed
        sum of edge labels (these are the boundary of each commuting square).

        Example output:
            Z{ a - b, c } / Z{ a + b - c - d }
        """
        edge_names = [e["label"] for e in self.edges]
        if not self.edges:
            return "0"

        ker1 = _int_nullspace(self._boundary1_matrix())
        if not ker1:
            return "0"

        if not self.squares:
            return _format_group(ker1, edge_names, [], edge_names)

        # Columns of ∂₂ are the image generators (boundaries of squares).
        # Skip zero columns — a square with ∂ = 0 lives in ker ∂₂ and
        # contributes to H₂, not to a relation in H₁.
        d2 = self._boundary2_matrix()
        im2 = [
            [d2[r][k] for r in range(len(self.edges))]
            for k in range(len(self.squares))
            if any(d2[r][k] != 0 for r in range(len(self.edges)))
        ]

        return _format_group(ker1, edge_names, im2, edge_names)

    def second_homology(self) -> str:
        """
        H₂ = ker ∂₂.

        Generators: integer nullspace of ∂₂, each expressed as a signed sum
        of square labels.  (Always free abelian — no 3-cells to impose torsion.)

        Example output:
            Z{ (a·b=b·a) }
        """
        square_names = [s["label"] for s in self.squares]
        if not self.squares:
            return "0"

        ker2 = _int_nullspace(self._boundary2_matrix())
        if not ker2:
            return "0"

        ker_strs = [_vec_to_label(v, square_names) for v in ker2]
        return "Z{ " + ", ".join(ker_strs) + " }"

    def homology_ranks(self) -> dict:
        """
        Return a dict with Betti numbers and torsion coefficients for H₁ and H₂,
        computed via Smith normal form.

        Keys: 'beta_1', 'torsion_1', 'beta_2'
        """
        edge_names = [e["label"] for e in self.edges]

        # β₁ = dim ker ∂₁ − rank ∂₂
        d1_diag = _snf_diagonal(self._boundary1_matrix()) if self.edges else []
        rank_d1 = sum(1 for d in d1_diag if d != 0)
        dim_ker_d1 = len(self.edges) - rank_d1

        d2_diag = _snf_diagonal(self._boundary2_matrix()) if self.squares else []
        rank_d2 = sum(1 for d in d2_diag if d != 0)
        torsion_1 = [d for d in d2_diag if d > 1]

        beta_1 = max(dim_ker_d1 - rank_d2, 0)
        beta_2 = max(len(self.squares) - rank_d2, 0)

        return {"beta_1": beta_1, "torsion_1": torsion_1, "beta_2": beta_2}

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def _square_boundary_str(self, s: dict) -> str:
        edge_names = [e["label"] for e in self.edges]
        col = [0] * len(self.edges)
        a, b = s["left"]
        c, d = s["right"]
        col[self._edge_index[a]] += 1
        col[self._edge_index[b]] += 1
        col[self._edge_index[c]] -= 1
        col[self._edge_index[d]] -= 1
        return _vec_to_label(col, edge_names)

    def summary(self) -> str:
        sep = "=" * 64
        mid = "-" * 64
        lines = [sep, "TwoGraph", mid]

        lines.append(f"  Vertices ({len(self.vertices)}):  {', '.join(self.vertices)}")

        lines.append(f"  Edges ({len(self.edges)}):")
        for e in self.edges:
            lines.append(f"    {e['label']:22s}  {e['source']} → {e['target']}")

        lines.append(f"  Squares ({len(self.squares)}):")
        if self.squares:
            for s in self.squares:
                bnd = self._square_boundary_str(s)
                lines.append(f"    {s['label']:38s}  ∂ = {bnd}")
        else:
            lines.append("    (none)")

        ranks = self.homology_ranks()
        torsion_note = (
            f"  [torsion coefficients: {ranks['torsion_1']}]"
            if ranks["torsion_1"] else ""
        )

        lines += [
            mid,
            f"  H₁ = {self.first_homology()}",
            *(([torsion_note]) if torsion_note else []),
            f"  H₂ = {self.second_homology()}",
            sep,
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"TwoGraph(vertices={len(self.vertices)}, "
            f"edges={len(self.edges)}, squares={len(self.squares)})"
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python two_graph.py <file>")
        sys.exit(1)

    g = TwoGraph.from_file(sys.argv[1])
    print(g.summary())
