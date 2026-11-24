#!/usr/bin/env python3
from __future__ import annotations
import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple
import matplotlib.pyplot as plt


@dataclass
class Node:
    name: str = ""
    length: float = 0.0
    children: List["Node"] = field(default_factory=list)
    x: float = 0.0
    y: float = 0.0

    def is_leaf(self) -> bool:
        return not self.children


def split_top_level(text: str, sep: str = ",") -> List[str]:
    parts: List[str] = []
    depth = 0
    in_quote = False
    start = 0
    for i, ch in enumerate(text):
        if ch == "'":
            in_quote = not in_quote
        elif not in_quote:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            elif ch == sep and depth == 0:
                parts.append(text[start:i].strip())
                start = i + 1
    parts.append(text[start:].strip())
    return [p for p in parts if p]


def parse_label_and_length(text: str) -> Tuple[str, float]:
    s = text.strip()
    if not s:
        return "", 0.0
    in_quote = False
    idx = -1
    for i in range(len(s) - 1, -1, -1):
        ch = s[i]
        if ch == "'":
            in_quote = not in_quote
        if ch == ":" and not in_quote:
            idx = i
            break
    if idx != -1:
        name = s[:idx].strip().strip("'")
        try:
            length = float(s[idx + 1 :].strip())
        except ValueError:
            length = 0.0
        return name, length
    return s.strip("'"), 0.0


def parse_newick(text: str) -> Node:
    s = text.strip()
    if s.endswith(";"):
        s = s[:-1]

    def parse_subtree(chunk: str) -> Node:
        chunk = chunk.strip()
        if chunk.startswith("("):
            depth = 0
            in_quote = False
            close_idx = -1
            for i, ch in enumerate(chunk):
                if ch == "'":
                    in_quote = not in_quote
                elif not in_quote:
                    if ch == "(":
                        depth += 1
                    elif ch == ")":
                        depth -= 1
                        if depth == 0:
                            close_idx = i
                            break
            inside = chunk[1:close_idx]
            suffix = chunk[close_idx + 1 :].strip()
            name, length = parse_label_and_length(suffix)
            node = Node(name=name, length=length)
            for part in split_top_level(inside, ","):
                node.children.append(parse_subtree(part))
            return node
        name, length = parse_label_and_length(chunk)
        return Node(name=name, length=length)

    return parse_subtree(s)


def assign_layout(root: Node, cladogram: bool = True) -> List[Node]:
    leaves: List[Node] = []

    def walk(node: Node, parent_x: float = 0.0) -> None:
        node.x = parent_x
        for child in node.children:
            step = 1.0 if cladogram else (child.length if child.length > 0 else 1.0)
            walk(child, node.x + step)

    walk(root, 0.0)

    y = 0

    def place(node: Node) -> None:
        nonlocal y
        if node.is_leaf():
            node.y = float(y)
            leaves.append(node)
            y += 1
            return
        for c in node.children:
            place(c)
        node.y = sum(c.y for c in node.children) / len(node.children)

    place(root)
    return leaves


def draw(root: Node, title: str, out_png: Path, out_svg: Path) -> None:
    leaves = assign_layout(root, cladogram=True)
    fig_h = max(6.0, 0.35 * len(leaves))
    fig, ax = plt.subplots(figsize=(16, fig_h))

    def edges(node: Node) -> None:
        for c in node.children:
            ax.plot([node.x, c.x], [c.y, c.y], color="#1a1a1a", lw=1.8)
            ax.plot([node.x, node.x], [node.y, c.y], color="#1a1a1a", lw=1.8)
            edges(c)

    edges(root)

    xmax = max(l.x for l in leaves) if leaves else 1.0
    pad = max(0.1, 0.03 * xmax)
    for l in leaves:
        ax.scatter([l.x], [l.y], s=25, color="black", zorder=3)
        ax.text(l.x + pad, l.y, l.name, va="center", ha="left", fontsize=11)

    ax.set_title(title, fontsize=15, pad=10)
    ax.set_xlabel("Cladogram depth (topology only)")
    ax.set_yticks([])
    ax.grid(axis="x", color="#d0d0d0", lw=0.7)
    ax.set_xlim(0, xmax + 6 * pad)
    ax.set_ylim(-1, len(leaves))
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tree", required=True)
    p.add_argument("--title", required=True)
    p.add_argument("--out-prefix", required=True)
    a = p.parse_args()

    tree = parse_newick(Path(a.tree).read_text())
    pref = Path(a.out_prefix)
    draw(tree, a.title, pref.with_suffix(".png"), pref.with_suffix(".svg"))


if __name__ == "__main__":
    main()
