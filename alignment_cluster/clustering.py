"""UPGMA hierarchical clustering over the raw pairwise distances.

Inter-cluster distance is the arithmetic mean of all original sequence-pair
distances across the two clusters, so clusters of different sizes are not
given equal weight.

Tie breaking and child ordering use each cluster's minimum leaf identifier:
candidate pairs are compared lexicographically as ordered identifier pairs.
Distances are never rounded before choosing a merge.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Cluster:
    """Mutable working cluster used while building the tree."""

    members: frozenset[str]
    min_identifier: str
    size: int


@dataclass
class TreeNode:
    """Serializable tree node."""

    id: str
    kind: str
    height: float
    branch_length: float
    children: list["TreeNode"] = field(default_factory=list)
    members: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        data: dict = {
            "id": self.id,
            "type": self.kind,
            "height": self.height,
            "branch_length": self.branch_length,
            "members": list(self.members),
        }
        if self.kind == "merged":
            data["children"] = [child.to_dict() for child in self.children]
        return data


def _cluster_distance(
    first: Cluster,
    second: Cluster,
    pair_distances: dict[tuple[str, str], float],
) -> float:
    total = 0.0
    for member_a in sorted(first.members):
        for member_b in sorted(second.members):
            total += pair_distances[(member_a, member_b)]
    return total / (len(first.members) * len(second.members))


def upgma(
    identifiers: list[str],
    pair_distances: dict[tuple[str, str], float],
) -> tuple[list[dict], TreeNode]:
    """Run UPGMA and return merge records plus the final tree.

    Args:
        identifiers: Identifiers in lexicographic order.
        pair_distances: Symmetric mapping of original pairwise distances;
            self distances must be zero.

    Returns:
        A tuple of merge records in occurrence order and the root node.
    """

    clusters: dict[int, Cluster] = {}
    leaf_nodes: dict[int, TreeNode] = {}
    for identifier in identifiers:
        node_id = identifier
        clusters[node_id] = Cluster(
            members=frozenset({identifier}),
            min_identifier=identifier,
            size=1,
        )
        leaf_nodes[node_id] = TreeNode(
            id=node_id,
            kind="leaf",
            height=0.0,
            branch_length=0.0,
            members=[identifier],
        )

    alive: set[int] = set(clusters)
    merge_records: list[dict] = []
    merged_nodes: dict[int, TreeNode] = dict(leaf_nodes)
    next_merge_number = 1

    def node_sort_key(node_id: int) -> str:
        return clusters[node_id].min_identifier

    while len(alive) > 1:
        alive_sorted = sorted(alive, key=node_sort_key)
        best: tuple[int, int] | None = None
        best_distance = 0.0
        best_key: tuple[str, str] | None = None

        for index, id_a in enumerate(alive_sorted):
            for id_b in alive_sorted[index + 1 :]:
                cluster_a = clusters[id_a]
                cluster_b = clusters[id_b]
                distance = _cluster_distance(
                    cluster_a, cluster_b, pair_distances
                )
                tie_key = (
                    cluster_a.min_identifier,
                    cluster_b.min_identifier,
                )
                if (
                    best is None
                    or distance < best_distance
                    or (distance == best_distance and tie_key < best_key)
                ):
                    best = (id_a, id_b)
                    best_distance = distance
                    best_key = tie_key

        assert best is not None and best_key is not None
        id_a, id_b = best
        child_a = merged_nodes[id_a]
        child_b = merged_nodes[id_b]
        if clusters[id_a].min_identifier > clusters[id_b].min_identifier:
            id_a, id_b = id_b, id_a
            child_a, child_b = child_b, child_a

        merged_cluster = Cluster(
            members=clusters[id_a].members | clusters[id_b].members,
            min_identifier=best_key[0],
            size=clusters[id_a].size + clusters[id_b].size,
        )
        merge_id = f"merge_{next_merge_number}"
        next_merge_number += 1
        height = best_distance / 2.0

        ordered_children = [child_a, child_b]
        for child in ordered_children:
            child.branch_length = height - child.height

        members = sorted(merged_cluster.members)
        merged_node = TreeNode(
            id=merge_id,
            kind="merged",
            height=height,
            branch_length=0.0,
            children=ordered_children,
            members=members,
        )

        merge_records.append(
            {
                "order": len(merge_records) + 1,
                "id": merge_id,
                "child_ids": [child_a.id, child_b.id],
                "child_min_identifiers": [
                    clusters[id_a].min_identifier,
                    clusters[id_b].min_identifier,
                ],
                "members": members,
                "inter_cluster_distance": best_distance,
                "height": height,
            }
        )

        alive.remove(id_a)
        alive.remove(id_b)
        clusters[merge_id] = merged_cluster
        merged_nodes[merge_id] = merged_node
        alive.add(merge_id)

    root_id = next(iter(alive))
    root = merged_nodes[root_id]
    root.branch_length = 0.0
    return merge_records, root
