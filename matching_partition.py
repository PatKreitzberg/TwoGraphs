from collections import defaultdict


def find_matching_and_partition(vertices):
    """
    Parameters
    ----------
    vertices : list of (e, f, u)
        e, f : objects with a .degree attribute that is 1 or 2
        u    : int

    Each vertex is matched to exactly one other vertex (g, h, w) satisfying:
        - e.degree != g.degree
        - u == w
    Matched first-elements (e, g) are required to land in the same partition set.

    Returns
    -------
    matching : list of (i, j) int pairs
        Indices into `vertices`. Each pair is a valid match.
    E1, E2 : sets
        Partition of S = {v[0] for v in vertices} into two non-empty disjoint sets.
        For every matched pair (i, j), vertices[i][0] and vertices[j][0] are
        in the same set.
    """
    n = len(vertices)
    if n == 0:
        raise ValueError("No vertices provided")
    if n % 2 != 0:
        raise ValueError("Odd number of vertices — perfect matching impossible")

    # ── Step 1: perfect matching ──────────────────────────────────────────────
    # Group vertex indices by u value.
    # Within each group, pair every degree-1 vertex with a degree-2 vertex.
    by_u = defaultdict(list)
    for i, (e, f, u) in enumerate(vertices):
        by_u[u].append(i)

    matching = []
    for u, indices in by_u.items():
        deg1 = [i for i in indices if vertices[i][0].degree == 1]
        deg2 = [i for i in indices if vertices[i][0].degree == 2]
        if len(deg1) != len(deg2):
            raise ValueError(
                f"u={u}: {len(deg1)} degree-1 vs {len(deg2)} degree-2 vertices — "
                "perfect matching impossible"
            )
        for i, j in zip(deg1, deg2):
            matching.append((i, j))

    # ── Step 2: union-find on elements of S ───────────────────────────────────
    # Matched first-elements must share an E set, so union them together.
    all_elems = list({v[0] for v in vertices})
    idx_of = {e: k for k, e in enumerate(all_elems)}
    parent = list(range(len(all_elems)))

    def find(x):
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:           # path compression
            parent[x], x = root, parent[x]
        return root

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for i, j in matching:
        union(idx_of[vertices[i][0]], idx_of[vertices[j][0]])

    # ── Step 3: collect components and partition ──────────────────────────────
    components = defaultdict(set)
    for elem in all_elems:
        components[find(idx_of[elem])].add(elem)
    comp_list = list(components.values())

    if len(comp_list) < 2:
        raise ValueError(
            "All elements of S form one connected component — "
            "cannot split into two non-empty sets"
        )

    # Any split across component boundaries is valid; put first component in E1.
    E1 = comp_list[0]
    E2 = set().union(*comp_list[1:])

    return matching, E1, E2


# ── Example usage ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    class Elem:
        def __init__(self, name, degree):
            self.name = name
            self.degree = degree
        def __repr__(self):
            return f"{self.name}(deg={self.degree})"

    a = Elem("a", 1)
    b = Elem("b", 2)
    c = Elem("c", 1)
    d = Elem("d", 2)
    e = Elem("e", 1)
    f = Elem("f", 2)

    vertices = [
        (a, b, 0),   # u=0, degree-1
        (b, a, 0),   # u=0, degree-2
        (c, d, 1),   # u=1, degree-1
        (d, c, 1),   # u=1, degree-2
        (e, f, 2),   # u=2, degree-1
        (f, e, 2),   # u=2, degree-2
    ]

    matching, E1, E2 = find_matching_and_partition(vertices)

    S = {v[0] for v in vertices}

    print("S  =", S)
    print("E1 =", E1)
    print("E2 =", E2)
    print("Matching (index pairs):", matching)
    print("Valid partition:", (E1 | E2 == S) and (len(E1 & E2) == 0) and bool(E1) and bool(E2))
