#!/usr/bin/env python3
"""Simulate a small kairos-lr (e2e) results/ directory for trying out hgt_viz.

Everything is synthetic: taxa, contigs, mobile elements, CARD and RIP-DB hits
(the CARD accessions are made up). Writes

  OUT/results/detect/NAME.{kairos_region_metadata, kairos_region_overlaps,
                           assignRes2, potential_hgt_regions.taxonomy,
                           ffn_coords}.tsv
  OUT/results/detect/NAME.features/events.tsv
  OUT/results/detect/NAME-card.tsv, NAME-ripdb.tsv
  OUT/results/call/insitu_events.tsv
  OUT/aro_index.tsv

and, with --get-hgt path/to/kairos-lr/scripts/get_hgt-v2.py, runs it to fill
OUT/results/detect/NAME.hgt/.
"""
import argparse
import os
import random
import subprocess
import sys

import pandas as pd

RANKS = ["d", "p", "c", "o", "f", "g", "s"]
RANK_NAMES = ["domain", "phylum", "class", "order", "family", "genus", "species"]
GENE = 1000  # every simulated gene occupies one 1 kb slot

CARD_GENES = [  # (name, AMR gene family, drug class, mechanism)
    ("TEM-1", "TEM beta-lactamase", "cephalosporin;penam", "antibiotic inactivation"),
    ("CTX-M-15", "CTX-M beta-lactamase", "cephalosporin", "antibiotic inactivation"),
    ("OXA-48", "OXA beta-lactamase", "carbapenem;cephalosporin;penam", "antibiotic inactivation"),
    ("tet(M)", "tetracycline-resistant ribosomal protection protein", "tetracycline",
     "antibiotic target protection"),
    ("tet(A)", "major facilitator superfamily (MFS) antibiotic efflux pump", "tetracycline",
     "antibiotic efflux"),
    ("sul1", "sulfonamide resistant sul", "sulfonamide", "antibiotic target replacement"),
    ("sul2", "sulfonamide resistant sul", "sulfonamide", "antibiotic target replacement"),
    ("aph(3'')-Ib", "APH(3'')", "aminoglycoside", "antibiotic inactivation"),
    ("aadA2", "ANT(3'')", "aminoglycoside", "antibiotic inactivation"),
    ("ermB", "Erm 23S ribosomal RNA methyltransferase", "macrolide;lincosamide;streptogramin",
     "antibiotic target alteration"),
    ("qnrS1", "quinolone resistance protein (qnr)", "fluoroquinolone", "antibiotic target protection"),
    ("vanA", "glycopeptide resistance gene cluster", "glycopeptide", "antibiotic target alteration"),
    ("mcr-1", "MCR phosphoethanolamine transferase", "peptide", "antibiotic target alteration"),
    ("catA1", "chloramphenicol acetyltransferase (CAT)", "phenicol", "antibiotic inactivation"),
]
RIP_FAMILIES = ["rep3_1", "rep3_2", "IncFII_repA", "repB_RCR", "rep1", "pri",
                "rep_trans", "repB_HUH", "orc"]


def make_tree(rng, n_species):
    phyla = {"Bacteria": ["Pseudomonadota", "Bacillota", "Bacillota_A", "Bacteroidota",
                          "Actinomycetota", "Chloroflexota", "Patescibacteria"],
             "Archaea": ["Halobacteriota", "Thermoproteota"]}
    lineages = set()
    while len(lineages) < n_species:
        d = "Bacteria" if rng.random() < 0.9 else "Archaea"
        p = rng.choice(phyla[d])
        c = f"{p.split('_')[0][:6]}ia_c{rng.randint(1, 3)}"
        o = f"{c}_o{rng.randint(1, 2)}"
        f = f"{o}_f{rng.randint(1, 2)}"
        g = f"Gen{rng.randint(1, 60)}" + ("_A" if rng.random() < 0.2 else "")
        s = f"{g} sp{rng.randint(1, 50)}"
        lineages.add((d, p, c, o, f, g, s))
    return sorted(lineages)


def tax_row(ident, lin, depth, rng):
    if depth == 0:
        return [ident, 1, "no rank", "unclassified", 3, 0, 0, 0.0, ""]
    frag = rng.randint(2, 40)
    ret = rng.randint(1, frag)
    agr = rng.randint(1, ret)
    return [ident, 1000 + abs(hash(lin[:depth])) % 10**6, RANK_NAMES[depth - 1],
            lin[depth - 1], frag, ret, agr, round(100 * agr / ret, 2),
            ";".join(f"{r}_{t}" for r, t in zip(RANKS[:depth], lin[:depth]))]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", default="demo")
    ap.add_argument("--contigs", type=int, default=600)
    ap.add_argument("--elements", type=int, default=60)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--get-hgt", help="path to kairos-lr scripts/get_hgt-v2.py")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    det_dir = os.path.join(a.out, "results", "detect")
    os.makedirs(os.path.join(det_dir, f"{a.name}.features"), exist_ok=True)
    os.makedirs(os.path.join(a.out, "results", "call"), exist_ok=True)
    det = os.path.join(det_dir, a.name)

    lineages = make_tree(rng, 80)
    contigs = []
    for i in range(a.contigs):
        n_genes = int(min(1500, max(8, rng.lognormvariate(4.2, 0.9))))
        contigs.append({"contig": f"contig_{i+1}", "lin": rng.choice(lineages),
                        "n_genes": n_genes, "len": n_genes * GENE + rng.randint(0, 900),
                        "depth": rng.choices([7, 6, 5, 4, 0], [55, 25, 10, 5, 5])[0]})
    by_name = {c["contig"]: c for c in contigs}

    # genes: gene k (1-based) occupies [(k-1)*1000+51, k*1000-50]
    genes = []
    for c in contigs:
        for k in range(1, c["n_genes"] + 1):
            genes.append((f"{c['contig']}_{k}", (k - 1) * GENE + 51, k * GENE - 50,
                          rng.choice([1, -1]), f"ID={k};partial=00"))
    strand_of = {g[0]: g[3] for g in genes}
    coords = {g[0]: (g[1], g[2]) for g in genes}

    # mobile elements with cargo
    elements = []
    for _ in range(a.elements):
        n = rng.randint(2, 25)
        cargo = {}
        if rng.random() < 0.45:
            for j in rng.sample(range(n), min(n, rng.randint(1, 3))):
                cargo[j] = ("card", rng.choice(CARD_GENES)[0])
        if rng.random() < 0.35:
            j = rng.randrange(n)
            cargo[j] = ("ripdb", rng.choice(RIP_FAMILIES))
            if n > 2 and rng.random() < 0.5:
                cargo[(j + 1) % n] = ("ripdb", "parA")
        elements.append({"origin": rng.choice(lineages), "n": n, "cargo": cargo,
                         "carriers": rng.sample(contigs, rng.choices([2, 3, 5, 10, 25],
                                                                     [40, 25, 20, 10, 5])[0])})

    placements = {c["contig"]: [] for c in contigs}
    for ei, el in enumerate(elements):
        for c in el["carriers"]:
            placements[c["contig"]].append(("el", ei, el["n"]))
    for c in contigs:
        for _ in range(rng.choices([0, 1, 2], [60, 30, 10])[0]):
            placements[c["contig"]].append(("island", None, rng.randint(1, 20)))

    meta, region_tax, region_orfs = [], [], {}
    for cname, plc in placements.items():
        c = by_name[cname]
        plc = [p for p in plc if p[2] < c["n_genes"] * 0.8]
        rng.shuffle(plc)
        cursor, k = 1, 0
        for kind, ei, n in plc:
            r = rng.random()
            if r < 0.08:
                first = cursor
            elif r < 0.16:
                first = c["n_genes"] - n + 1
            else:
                first = rng.randint(cursor, max(cursor, c["n_genes"] - n + 1))
            if first < cursor or first + n - 1 > c["n_genes"]:
                continue
            last = first + n - 1
            cursor = last + 2
            k += 1
            ridx = f"{cname}_region_{k}"
            start, stop = coords[f"{cname}_{first}"][0] - 1, coords[f"{cname}_{last}"][1]
            if stop > c["len"] - 1000:
                rtype = "_end_edge_"
            elif start < 1001:
                rtype = "_start_edge_"
            else:
                rtype = "_internal_"
            orfs = [f"{cname}_{j}" for j in range(first, last + 1)]
            region_orfs[ridx] = (kind, ei, orfs)
            title = f"{ridx}{rtype}length={stop - start}"
            meta.append({"contig": cname, "start": start, "stop": stop, "orf_cts": n,
                         "orfs_encompassed": ",".join(orfs), "region_idx": ridx,
                         "region_len": stop - start, "contig_length": c["len"],
                         "end_edge_cutoff": c["len"] - 1000, "start_edge_cutoff": 1001,
                         "region_type": rtype, "region_title": title})
            src = (elements[ei]["origin"] if kind == "el"
                   else rng.choice(lineages) if rng.random() < 0.6 else c["lin"])
            depth = rng.choices([7, 6, 5, 4, 3, 0], [40, 25, 12, 8, 5, 10])[0]
            region_tax.append(tax_row(f"{title}::{cname}:{start}-{stop}", src, depth, rng))

    meta_df = pd.DataFrame(meta)
    meta_df.to_csv(f"{det}.kairos_region_metadata.tsv", sep="\t", index=False)
    pd.DataFrame(region_tax).to_csv(f"{det}.potential_hgt_regions.taxonomy.tsv",
                                    sep="\t", index=False, header=False)
    pd.DataFrame([tax_row(c["contig"], c["lin"], c["depth"], rng) for c in contigs]).to_csv(
        f"{det}.assignRes2.tsv", sep="\t", index=False, header=False)
    pd.DataFrame(genes).to_csv(f"{det}.ffn_coords.tsv", sep="\t", index=False, header=False)

    # overlaps (one direction per pair, as regions-v2.py writes) and events
    by_elem = {}
    for ridx, (kind, ei, orfs) in region_orfs.items():
        if kind == "el":
            by_elem.setdefault(ei, []).append((ridx, orfs))
    contig_of = dict(zip(meta_df["region_idx"], meta_df["contig"]))
    reg = meta_df.set_index("region_idx")
    ov, ev = [], []
    for ei, regs in by_elem.items():
        for i in range(len(regs)):
            for j in range(i + 1, len(regs)):
                (r1, o1), (r2, o2) = regs[i], regs[j]
                for g in sorted(rng.sample(range(len(o1)), rng.randint(1, len(o1)))):
                    ov.append({"contig1": contig_of[r1], "contig2": contig_of[r2],
                               "region_idx_contig1": r1, "region_idx_contig2": r2,
                               "orf_region1": o1[g], "orf_region2": o2[g]})
                ev.append({"event_id": f"EV{len(ev) + 1:07d}",
                           "region_idx_1": r1, "contig_1": contig_of[r1],
                           "shared_start_1": reg.at[r1, "start"] + rng.randint(-40, 40),
                           "shared_end_1": reg.at[r1, "stop"] + rng.randint(-40, 40),
                           "bp_left_1": "", "bp_right_1": "", "empty_site_1": "",
                           "region_idx_2": r2, "contig_2": contig_of[r2],
                           "shared_start_2": reg.at[r2, "start"] + rng.randint(-40, 40),
                           "shared_end_2": reg.at[r2, "stop"] + rng.randint(-40, 40),
                           "bp_left_2": "", "bp_right_2": "", "empty_site_2": "",
                           "strand": 1, "identity": round(rng.uniform(0.985, 1), 4),
                           "shared_len": int(reg.at[r1, "region_len"]), "refined": True,
                           "region_type_1": reg.at[r1, "region_type"],
                           "region_type_2": reg.at[r2, "region_type"]})
    pd.DataFrame(ov).to_csv(f"{det}.kairos_region_overlaps.tsv", sep="\t", index=False)
    pd.DataFrame(ev).to_csv(os.path.join(f"{det}.features", "events.tsv"), sep="\t", index=False)

    # diamond best hits: element cargo + a few chromosomal genes
    def qtitle(orf):
        s, e = coords[orf]
        return f"{orf} # {s} # {e} # {strand_of[orf]} # ID=x;partial=00"

    card_acc = {g[0]: f"ARO:39{i:05d}" for i, g in enumerate(CARD_GENES)}
    card, rip = [], []
    for ridx, (kind, ei, orfs) in region_orfs.items():
        if kind != "el":
            continue
        for j, (db, label) in elements[ei]["cargo"].items():
            pid = round(rng.uniform(92, 100), 1)
            if db == "card":
                card.append([qtitle(orfs[j]), f"gb|DEMO{j:04d}.1|{card_acc[label]}|{label} [demo]",
                             pid, round(rng.uniform(300, 900), 1), 1e-120])
            else:
                rip.append([qtitle(orfs[j]), f"DEMO_{ei}_{j}|{label}",
                            round(rng.uniform(25, 100), 1), round(rng.uniform(80, 600), 1), 1e-40])
    in_region = {o for (_, _, orfs) in region_orfs.values() for o in orfs}
    for c in rng.sample(contigs, len(contigs) // 15):
        orf = f"{c['contig']}_{rng.randint(1, c['n_genes'])}"
        if orf not in in_region:
            label = rng.choice(CARD_GENES)[0]
            card.append([qtitle(orf), f"gb|DEMOCHR.1|{card_acc[label]}|{label} [demo]",
                         round(rng.uniform(90, 100), 1), 500.0, 1e-100])
    pd.DataFrame(card).to_csv(f"{det}-card.tsv", sep="\t", index=False, header=False)
    pd.DataFrame(rip).to_csv(f"{det}-ripdb.tsv", sep="\t", index=False, header=False)
    pd.DataFrame([{"ARO Accession": card_acc[n], "CVTERM ID": "", "Model Sequence ID": "",
                   "Model ID": "", "Model Name": n, "ARO Name": n, "Protein Accession": "",
                   "DNA Accession": "", "AMR Gene Family": fam, "Drug Class": dc,
                   "Resistance Mechanism": mech, "CARD Short Name": n}
                  for n, fam, dc, mech in CARD_GENES]).to_csv(
        os.path.join(a.out, "aro_index.tsv"), sep="\t", index=False)

    # a handful of in situ calls on random events
    ins = []
    for e in rng.sample(ev, min(40, len(ev))):
        ins.append({"group": "G", "timepoint": rng.randint(2, 50), "event_id": e["event_id"],
                    "direction": "directed", "tier": rng.choice("AB"),
                    "recipient_contig": e["contig_1"], "donor_contig": e["contig_2"],
                    "persistence": rng.choice(["yes", "no", "partial"])})
    pd.DataFrame(ins).to_csv(os.path.join(a.out, "results", "call", "insitu_events.tsv"),
                             sep="\t", index=False)

    print(f"{len(contigs)} contigs, {len(meta_df)} regions, {len(ev)} events, "
          f"{len(card)} CARD hits, {len(rip)} RIP-DB hits")
    if a.get_hgt:
        subprocess.run([sys.executable, a.get_hgt, "--path", det_dir + "/", "--sample", a.name,
                        "--out-dir", f"{det}.hgt"], check=True, stdout=subprocess.DEVNULL)
        print(f"get_hgt-v2.py outputs in {det}.hgt/")


if __name__ == "__main__":
    main()
