import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import numpy as np


def draw_multidigraph(g):
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("#1a1a2e")
    ax.set_facecolor("#1a1a2e")

    # --- Layout ---
    pos = nx.spring_layout(g, seed=42)

    # Spread nodes a bit more if only a few
    if len(pos) <= 3:
        keys = list(pos.keys())
        if len(keys) == 1:
            pos[keys[0]] = np.array([0.0, 0.0])
        elif len(keys) == 2:
            pos[keys[0]] = np.array([-0.5, 0.0])
            pos[keys[1]] = np.array([0.5, 0.0])

    # --- Draw edges ---
    # Group edges by (u, v) pair so we can fan them out
    edge_groups = {}
    for u, v, data in g.edges(data=True):
        key = (u, v)
        edge_groups.setdefault(key, []).append(data)

    label_positions = []  # [(x, y, label, color)]

    for (u, v), edges in edge_groups.items():
        n = len(edges)
        is_loop = u == v

        for i, data in enumerate(edges):
            color = data.get("color", "white")
            label = data.get("label", "")

            if is_loop:
                # Draw a self-loop as a circular arc above the node
                x, y = pos[u]
                angle_offset = (i - (n - 1) / 2) * 0.35  # fan loops sideways
                loop_radius = 0.12 + i * 0.04

                theta = np.linspace(0 + angle_offset, 2 * np.pi + angle_offset, 200)
                lx = x + loop_radius * np.cos(theta)
                ly = y + loop_radius * 1.5 * np.sin(theta) + loop_radius * 1.2

                ax.plot(lx, ly, color=color, linewidth=1.8, zorder=2)

                # Arrow at the end of the loop
                dx = lx[-1] - lx[-2]
                dy = ly[-1] - ly[-2]
                ax.annotate(
                    "",
                    xy=(lx[-1], ly[-1]),
                    xytext=(lx[-1] - dx * 3, ly[-1] - dy * 3),
                    arrowprops=dict(
                        arrowstyle="-|>",
                        color=color,
                        lw=1.5,
                        mutation_scale=14,
                    ),
                    zorder=3,
                )

                # Label at the top of the loop
                mid_idx = len(lx) // 2
                label_positions.append((lx[mid_idx], ly[mid_idx] + 0.03, label, color))

            else:
                # Fan out multiple edges between same pair using arc curvature
                # Spread: center edge is straight-ish, others curve left/right
                spread = 0.25
                if n == 1:
                    rad = 0.15  # slight curve even for single edge
                else:
                    rad = -spread + i * (2 * spread / (n - 1)) if n > 1 else 0.0

                src = pos[u]
                dst = pos[v]

                # Draw curved arrow
                arrow = FancyArrowPatch(
                    posA=src,
                    posB=dst,
                    connectionstyle=f"arc3,rad={rad:.3f}",
                    arrowstyle="-|>",
                    color=color,
                    linewidth=1.8,
                    mutation_scale=16,
                    zorder=2,
                    shrinkA=12,
                    shrinkB=12,
                )
                ax.add_patch(arrow)

                # Label at midpoint of the arc
                mx = (src[0] + dst[0]) / 2
                my = (src[1] + dst[1]) / 2
                # Offset perpendicular to the edge
                dx = dst[0] - src[0]
                dy = dst[1] - src[1]
                length = np.sqrt(dx**2 + dy**2) or 1
                perp = np.array([-dy, dx]) / length
                offset = perp * rad * 0.7
                lx = mx + offset[0]
                ly = my + offset[1]
                label_positions.append((lx, ly, label, color))

    # --- Draw nodes ---
    node_radius = 0.07
    for node, (x, y) in pos.items():
        circle = plt.Circle(
            (x, y),
            node_radius,
            color="#e0e0ff",
            zorder=4,
            linewidth=2,
            ec="#aaaacc",
        )
        ax.add_patch(circle)
        ax.text(
            x,
            y,
            str(node),
            ha="center",
            va="center",
            fontsize=14,
            fontweight="bold",
            color="#1a1a2e",
            zorder=5,
            fontfamily="monospace",
        )

    # --- Draw edge labels ---
    for lx, ly, label, color in label_positions:
        ax.text(
            lx,
            ly,
            label,
            ha="center",
            va="center",
            fontsize=9,
            color=color,
            fontfamily="monospace",
            bbox=dict(
                boxstyle="round,pad=0.2",
                fc="#1a1a2e",
                ec=color,
                alpha=0.85,
                linewidth=1,
            ),
            zorder=6,
        )

    # --- Legend for edge colors ---
    seen_colors = {}
    for _, _, data in g.edges(data=True):
        c = data.get("color", "white")
        if c not in seen_colors:
            seen_colors[c] = c
    legend_patches = [
        mpatches.Patch(color=c, label=c) for c in seen_colors
    ]
    ax.legend(
        handles=legend_patches,
        loc="lower right",
        framealpha=0.3,
        facecolor="#1a1a2e",
        edgecolor="#555577",
        labelcolor="white",
        fontsize=10,
    )

    ax.set_title(
        "MultiDiGraph Visualization",
        color="#ccccee",
        fontsize=14,
        pad=12,
        fontfamily="monospace",
    )

    plt.tight_layout()
    plt.savefig("multidigraph.png", dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    print("Saved to multidigraph.png")
    plt.show()


# ── Example ──────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    g = nx.MultiDiGraph()
    edge_list = [
        ("v", "u", {"color": "dodgerblue",  "label": "b_1"}),
        ("v", "u", {"color": "dodgerblue",  "label": "b_2"}),
        ("v", "u", {"color": "dodgerblue",  "label": "b_3"}),
        ("v", "v", {"color": "tomato",       "label": "r_1"}),
        ("v", "v", {"color": "tomato",       "label": "r_2"}),
    ]
    g.add_edges_from(edge_list)
    draw_multidigraph(g)
