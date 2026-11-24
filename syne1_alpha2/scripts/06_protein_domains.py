#!/usr/bin/env python3
"""Step 6 - protein domain map for the giant isoform of every view.

Human domains come from UniProt Q8NF91 (SYNE1_HUMAN; canonical sequence identical to
NP_892006.3): calponin-homology domains CH1/CH2, 74 spectrin repeats, the transmembrane
segment and the KASH domain, plus the LEWD kinesin-binding motif (found by sequence).
They are projected onto the giant protein of mouse, rat and zebrafish (and onto the GRCz12tu
models) through a pairwise FAMSA alignment with the human protein. The browser then projects
these onto any other isoform through shared codons (viewer/variant_engine.js).
"""
from __future__ import annotations

import sys
from pathlib import Path

from pyfamsa import Aligner, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, http, load_json, save_json  # noqa: E402

GIANT_REF = {
    "human_SYNE1": "NM_182961.4", "mouse_Syne1": "NM_001347732.2", "rat_Syne1": "NM_001419849.1",
    "zebrafish_syne1a": "ENSDART00000139350", "zebrafish_syne1b": "ENSDART00000148520",
    "zebrafish_syne1a_GRCz12tu": "XM_073933275.1", "zebrafish_syne1b_GRCz12tu": "XM_087227815.1",
}


def human_domains(uni: dict) -> list[dict]:
    out = []
    for f in uni["features"]:
        a, b = f["location"]["start"]["value"], f["location"]["end"]["value"]
        d = f.get("description", "")
        if f["type"] == "Domain" and d.startswith("Calponin"):
            out.append(dict(start=a, end=b, label=d.replace("Calponin-homology (CH)", "CH").strip(), type="CH"))
        elif f["type"] == "Repeat" and d.startswith("Spectrin"):
            out.append(dict(start=a, end=b, label="SR" + d.split()[-1], type="SR"))
        elif f["type"] == "Transmembrane":
            out.append(dict(start=a, end=b, label="TM", type="TM"))
        elif f["type"] == "Domain" and d == "KASH":
            out.append(dict(start=a, end=b, label="KASH", type="KASH"))
    seq = uni["sequence"]["value"]
    i = seq.find("LEWD")
    if i >= 0:
        out.append(dict(start=i + 1, end=i + 4, label="LEWD", type="motif"))
    order = {"CH": 0, "SR": 1, "TM": 2, "KASH": 3, "motif": 4}
    return sorted(out, key=lambda d: (order[d["type"]], d["start"]))


def project(human_seq: str, target_seq: str, doms: list[dict]) -> list[dict]:
    msa = Aligner(guide_tree="upgma").align([Sequence(b"h", human_seq.encode()), Sequence(b"t", target_seq.encode())])
    rows = {s.id.decode(): s.sequence.decode() for s in msa}
    h, t = rows["h"], rows["t"]
    hmap, hi, ti = {}, 0, 0
    for a, b in zip(h, t):
        if a != "-":
            hi += 1
        if b != "-":
            ti += 1
        if a != "-" and b != "-":
            hmap[hi] = ti
    out = []
    for d in doms:
        mapped = [hmap[r] for r in range(d["start"], d["end"] + 1) if r in hmap]
        if len(mapped) >= 0.5 * (d["end"] - d["start"] + 1):
            out.append(dict(d, start=min(mapped), end=max(mapped)))
    for d in out:  # re-anchor the LEWD motif on the target sequence
        if d["type"] == "motif":
            j = target_seq.find("LEWD", max(0, d["start"] - 40))
            if j >= 0 and abs(j + 1 - d["start"]) < 40:
                d["start"], d["end"] = j + 1, j + 4
    return out


def main() -> None:
    uni = http("https://rest.uniprot.org/uniprotkb/Q8NF91.json", kind="ensembl", as_json=True)
    built = load_json(DATA / "models_built.json")

    def prot(key):
        tid = GIANT_REF[key]
        for t in built[key]["transcripts"]:
            if t["id"].split(".")[0] == tid.split(".")[0]:
                return t["protein"]
        raise KeyError(tid)
    hs = prot("human_SYNE1")
    assert hs == uni["sequence"]["value"], "UniProt Q8NF91 canonical differs from NP_892006.3"
    hd = human_domains(uni)
    out = {"source": f"UniProt Q8NF91 (entry version {uni['entryAudit']['entryVersion']}, "
                     f"sequence version {uni['entryAudit']['sequenceVersion']}); projected by FAMSA alignment",
           "views": {"human_SYNE1": hd}}
    for key in GIANT_REF:
        if key == "human_SYNE1":
            continue
        out["views"][key] = project(hs, prot(key), hd)
        n = {t: sum(1 for d in out["views"][key] if d["type"] == t) for t in ("CH", "SR", "TM", "KASH", "motif")}
        print(f"{key}: {n}")
    print("human:", {t: sum(1 for d in hd if d["type"] == t) for t in ("CH", "SR", "TM", "KASH", "motif")})
    save_json(DATA / "domains.json", out)


if __name__ == "__main__":
    main()
