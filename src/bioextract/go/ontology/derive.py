from collections import defaultdict, deque

from .constant import HIERARCHICAL_RELATION_TYPES, LARGE_DISTANCE
from .model import (
    AncestorColumnBuffer,
    DepthColumnBuffer,
    EdgeColumnBuffer,
    TermColumnBuffer,
)


# #region GraphDerivation
def derive_topo_order(
    indegree: dict[str, int],
    children: dict[str, list[str]],
) -> list[str]:
    root_nodes = sorted(
        go_id for go_id, indegree_value in indegree.items() if indegree_value == 0
    )
    queue_zero_indegree_nodes = deque(root_nodes)
    topo_order_root_to_leaf: list[str] = []
    map_indegree_remaining = indegree.copy()

    while queue_zero_indegree_nodes:
        current_node = queue_zero_indegree_nodes.popleft()
        topo_order_root_to_leaf.append(current_node)
        for child_node in sorted(children.get(current_node, [])):
            map_indegree_remaining[child_node] -= 1
            if map_indegree_remaining[child_node] == 0:
                queue_zero_indegree_nodes.append(child_node)

    if len(topo_order_root_to_leaf) != len(indegree):
        raise ValueError(
            "GO edge graph contains a cycle; ancestor and depth derivation require a DAG."
        )

    return topo_order_root_to_leaf


def calculate_depth_range_from_parents(
    parents: list[str],
    min_depth_from_root: dict[str, int],
    max_depth_from_root: dict[str, int],
) -> tuple[int, int]:
    parent_min_depths = [min_depth_from_root[parent_go_id] for parent_go_id in parents]
    parent_max_depths = [max_depth_from_root[parent_go_id] for parent_go_id in parents]
    return min(parent_min_depths) + 1, max(parent_max_depths) + 1


def derive_ancestor_min_distance_from_parents(
    parents: list[str],
    min_distance_to_ancestors: dict[str, dict[str, int]],
) -> dict[str, int]:
    map_min_distance_to_ancestors: dict[str, int] = {}
    for parent_node in parents:
        map_min_distance_to_ancestors[parent_node] = min(
            map_min_distance_to_ancestors.get(parent_node, LARGE_DISTANCE),
            1,
        )
        for ancestor_node, min_distance in min_distance_to_ancestors[
            parent_node
        ].items():
            distance_to_ancestor_via_parent = min_distance + 1
            map_min_distance_to_ancestors[ancestor_node] = min(
                map_min_distance_to_ancestors.get(ancestor_node, LARGE_DISTANCE),
                distance_to_ancestor_via_parent,
            )

    return map_min_distance_to_ancestors


def derive_graph_metrics(
    nodes: set[str],
    parents: dict[str, list[str]],
    topo_order: list[str],
) -> tuple[dict[str, dict[str, int]], dict[str, int], dict[str, int]]:
    map_min_distance_to_ancestors: dict[str, dict[str, int]] = {
        node: {} for node in nodes
    }
    map_min_depth_from_root: dict[str, int] = {}
    map_max_depth_from_root: dict[str, int] = {}

    for node in topo_order:
        direct_parent_nodes = parents.get(node, [])
        if not direct_parent_nodes:
            map_min_depth_from_root[node] = 0
            map_max_depth_from_root[node] = 0
            continue

        map_min_depth_from_root[node], map_max_depth_from_root[node] = (
            calculate_depth_range_from_parents(
                parents=direct_parent_nodes,
                min_depth_from_root=map_min_depth_from_root,
                max_depth_from_root=map_max_depth_from_root,
            )
        )
        map_min_distance_to_ancestors[node] = derive_ancestor_min_distance_from_parents(
            parents=direct_parent_nodes,
            min_distance_to_ancestors=map_min_distance_to_ancestors,
        )

    return (
        map_min_distance_to_ancestors,
        map_min_depth_from_root,
        map_max_depth_from_root,
    )


def derive_graph_tables(
    term_data: TermColumnBuffer,
    edge_data: EdgeColumnBuffer,
) -> tuple[AncestorColumnBuffer, DepthColumnBuffer]:
    map_namespace = dict(zip(term_data.go_id, term_data.namespace, strict=True))
    set_nodes = set(map_namespace)
    map_parents: dict[str, list[str]] = defaultdict(list)
    map_children: dict[str, list[str]] = defaultdict(list)
    map_indegree: dict[str, int] = dict.fromkeys(set_nodes, 0)

    for child_node, parent_node, relation_type in zip(
        edge_data.child_go_id,
        edge_data.parent_go_id,
        edge_data.relation_type,
        strict=True,
    ):
        if relation_type not in HIERARCHICAL_RELATION_TYPES:
            continue
        map_parents[child_node].append(parent_node)
        map_children[parent_node].append(child_node)
        map_indegree[child_node] = map_indegree.get(child_node, 0) + 1
        map_indegree.setdefault(parent_node, 0)

    topo_order_root_to_leaf = derive_topo_order(map_indegree, map_children)
    (
        map_min_distance_to_ancestors,
        map_min_depth_from_root,
        map_max_depth_from_root,
    ) = derive_graph_metrics(
        nodes=set_nodes,
        parents=map_parents,
        topo_order=topo_order_root_to_leaf,
    )

    ancestor_cols = AncestorColumnBuffer()
    depth_cols = DepthColumnBuffer()
    for node in sorted(set_nodes):
        for ancestor_node, min_distance in sorted(
            map_min_distance_to_ancestors[node].items()
        ):
            if node == ancestor_node:
                continue
            ancestor_cols.go_id.append(node)
            ancestor_cols.ancestor_go_id.append(ancestor_node)
            ancestor_cols.min_distance.append(min_distance)

        depth_cols.go_id.append(node)
        depth_cols.namespace.append(map_namespace[node])
        depth_cols.min_depth_from_root.append(map_min_depth_from_root[node])
        depth_cols.max_depth_from_root.append(map_max_depth_from_root[node])

    return ancestor_cols, depth_cols

    # #endregion

    return ancestor_cols, depth_cols


# #endregion
