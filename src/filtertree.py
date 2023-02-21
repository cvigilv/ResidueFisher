import os
import sys
import pprint

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

from tqdm import tqdm
from itertools import product
from collections import defaultdict
from ete3 import Tree, NodeStyle
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

if __name__ == "__main__":
    print("filtertree.py")

    # Read tree file into string
    treestring = open(sys.argv[1], "r").read().replace("\n", "")
    output_source = os.path.basename(sys.argv[1]).split('.')[1]
    output_path = os.path.dirname(sys.argv[1])

    # Convert tree string into Tree object
    t = Tree(treestring, format=1)
    nleafs = len(t.get_leaves())
    print(f"Number of leafs: {nleafs}")

    # Get distances between all leaves
    source = []
    target = []
    distance = []
    for n1, n2 in product(t.get_leaves(), repeat=2):
        source.append(n1.name)
        target.append(n2.name)
        distance.append(t.get_distance(n1, n2))

    dmat = pd.DataFrame(
        {"source": source, "target": target, "distance": distance}
    ).pivot_table(index="source", columns="target", values="distance")

    # Cluster leafs based in distance using k-means
    clustering = defaultdict(list)
    if nleafs > 100:
        print("Number of hits {nleafs} > 100. Testing 100 evenly spaced number of clusters.")
    clusters = np.unique(np.floor(np.linspace(2, nleafs, 100)).astype(int))
    for i in tqdm(clusters):
        kmeans = KMeans(n_clusters=i, random_state=0).fit(dmat.values)
        silhouette = silhouette_score(
            dmat.values, kmeans.labels_, metric="precomputed", random_state=0
        )

        clustering["nclusters"].append(i)
        clustering["silhouette"].append(silhouette)
        clustering["labels"].append(kmeans.labels_)

    cluster_results = pd.DataFrame(clustering)
    best_cluster = (
        cluster_results.sort_values(by="silhouette", ascending=False).iloc[0].to_dict()
    )

    print(f"Best cluster found:")
    pprint.pprint(best_cluster)

    # Create visualization
    cmap = plt.get_cmap("rainbow", max(best_cluster["labels"]) + 1)

    sns.clustermap(
        dmat, cmap="Spectral_r", row_colors=[cmap(l) for l in best_cluster["labels"]]
    )
    plt.savefig(f"{output_path}/tree.matrix_complete.{output_source}.pdf", dpi=300)

    leaf2group = dict(zip(dmat.columns, best_cluster["labels"]))
    for i, leaf in enumerate(t.get_leaves()):
        nstyle = NodeStyle()
        nstyle["bgcolor"] = mpl.colors.to_hex(cmap(leaf2group[leaf.name]))
        leaf.set_style(nstyle)
    t.render(f"{output_path}/tree.complete.{output_source}.pdf", dpi=300)
    print("Leaf group assignation:")
    pprint.pprint(leaf2group)

    group_representatives = {}
    for leaf, group in leaf2group.items():
        group_representatives[group] = leaf
    print("Clusters representative hits:")
    pprint.pprint(group_representatives)

    sns.clustermap(
        dmat.loc[group_representatives.values(), group_representatives.values()],
        cmap="Spectral_r",
        row_colors=[cmap(l) for l in group_representatives.keys()],
    )
    plt.savefig(f"{output_path}/tree.matrix_representatives.{output_source}.pdf", dpi=300)

    t2 = t.copy()
    for leaf in t2.get_leaves():
        if leaf.name not in group_representatives.values():
            leaf.delete()
    t2.render(f"{output_path}/tree.representatives.{output_source}.pdf", dpi=300)
    print(f"Number of leafs after filtering: {len(t2.get_leaves())}")

    # Write filtered tree to file
    with open(sys.argv[2], "w+") as io:
        io.write(t2.write(format=1))
