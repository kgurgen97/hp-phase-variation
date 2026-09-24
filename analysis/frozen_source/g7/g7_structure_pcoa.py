#!/usr/bin/env python3
"""G7 section 4: PCoA (classical MDS) and simple hierarchical clustering on
the SCREENING_STRUCTURE_ONLY k-mer sketch distance matrix
(results/g7/g7_distance_matrix.tsv). Pure Python (no numpy/scipy available
in this environment); eigendecomposition via the cyclic Jacobi algorithm,
appropriate for the small (32x32) symmetric matrix here.

Run: python3 scripts/g7/g7_structure_pcoa.py
"""
import csv
import os

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DIST_TSV = os.path.join(REPO, "results", "g7", "g7_distance_matrix.tsv")
MANIFEST_TSV = os.path.join(REPO, "metadata", "g7", "g7_population_structure_manifest.tsv")
OUT_COORDS = os.path.join(REPO, "results", "g7", "g7_structure_coordinates.tsv")
OUT_CLUSTERS = os.path.join(REPO, "results", "g7", "g7_structure_cluster_summary.tsv")


def load_matrix():
    with open(DIST_TSV, newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        ids = header[1:]
        mat = []
        for row in r:
            mat.append([float(x) for x in row[1:]])
    return ids, mat


def jacobi_eigen(a, tol=1e-10, max_sweeps=100):
    """Cyclic Jacobi eigenvalue algorithm for a symmetric matrix `a`
    (list of lists). Returns (eigenvalues, eigenvectors) with eigenvectors
    as columns of a matrix (list of rows)."""
    n = len(a)
    a = [row[:] for row in a]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(max_sweeps):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        if off < tol:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-14:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                sign = 1.0 if theta >= 0 else -1.0
                t = sign / (abs(theta) + (theta ** 2 + 1) ** 0.5)
                c = 1.0 / (t ** 2 + 1) ** 0.5
                s = t * c
                app, aqq, apq = a[p][p], a[q][q], a[p][q]
                a[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq
                a[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq
                a[p][q] = 0.0
                a[q][p] = 0.0
                for i in range(n):
                    if i != p and i != q:
                        aip, aiq = a[i][p], a[i][q]
                        a[i][p] = c * aip - s * aiq
                        a[p][i] = a[i][p]
                        a[i][q] = s * aip + c * aiq
                        a[q][i] = a[i][q]
                for i in range(n):
                    vip, viq = v[i][p], v[i][q]
                    v[i][p] = c * vip - s * viq
                    v[i][q] = s * vip + c * viq
    eigvals = [a[i][i] for i in range(n)]
    return eigvals, v


def classical_mds(dist, n_components=3):
    n = len(dist)
    d2 = [[dist[i][j] ** 2 for j in range(n)] for i in range(n)]
    row_mean = [sum(d2[i]) / n for i in range(n)]
    grand_mean = sum(row_mean) / n
    b = [[-0.5 * (d2[i][j] - row_mean[i] - row_mean[j] + grand_mean) for j in range(n)] for i in range(n)]
    eigvals, eigvecs = jacobi_eigen(b)
    order = sorted(range(n), key=lambda k: -eigvals[k])[:n_components]
    coords = [[eigvecs[i][k] * (max(eigvals[k], 0.0) ** 0.5) for k in order] for i in range(n)]
    top_eigvals = [eigvals[k] for k in order]
    total_var = sum(abs(e) for e in eigvals)
    pct = [abs(e) / total_var * 100 if total_var else 0.0 for e in top_eigvals]
    return coords, top_eigvals, pct


def hierarchical_average_linkage(ids, dist):
    n = len(ids)
    clusters = {i: [ids[i]] for i in range(n)}
    active = list(range(n))
    cur_dist = {(i, j): dist[i][j] for i in range(n) for j in range(n) if i != j}
    next_id = n
    merges = []
    while len(active) > 1:
        best = None
        for ii in range(len(active)):
            for jj in range(ii + 1, len(active)):
                i, j = active[ii], active[jj]
                d = cur_dist[(i, j)]
                if best is None or d < best[0]:
                    best = (d, i, j)
        d, i, j = best
        merges.append((clusters[i], clusters[j], d))
        new_members = clusters[i] + clusters[j]
        clusters[next_id] = new_members
        new_active = [x for x in active if x not in (i, j)]
        for k in new_active:
            di = cur_dist.get((i, k), cur_dist.get((k, i)))
            dj = cur_dist.get((j, k), cur_dist.get((k, j)))
            ni, nj = len(clusters[i]), len(clusters[j])
            avg = (di * ni + dj * nj) / (ni + nj)
            cur_dist[(next_id, k)] = avg
            cur_dist[(k, next_id)] = avg
        active = new_active + [next_id]
        next_id += 1
    return merges


def main():
    ids, dist = load_matrix()
    with open(MANIFEST_TSV, newline="") as fh:
        manifest = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    coords, eigvals, pct = classical_mds(dist, n_components=3)
    with open(OUT_COORDS, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "PCo1", "PCo2", "PCo3", "role", "disease_label", "study", "country"])
        for idx, sid in enumerate(ids):
            m = manifest.get(sid, {})
            w.writerow([sid] + [f"{c:.6f}" for c in coords[idx]] + [
                m.get("role", "NA"), m.get("disease_label", "NA"),
                m.get("study", "NA"), m.get("country", "NA")])

    merges = hierarchical_average_linkage(ids, dist)
    with open(OUT_CLUSTERS, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["merge_order", "cluster_size", "merged_distance", "members"])
        for order, (a, b, d) in enumerate(merges, start=1):
            members = a + b
            w.writerow([order, len(members), f"{d:.6f}", ";".join(members)])

    print("PCo eigenvalues (top 3):", [round(e, 4) for e in eigvals])
    print("PCo variance pct (top 3):", [round(p, 2) for p in pct])
    print("wrote", OUT_COORDS)
    print("wrote", OUT_CLUSTERS)


if __name__ == "__main__":
    main()
