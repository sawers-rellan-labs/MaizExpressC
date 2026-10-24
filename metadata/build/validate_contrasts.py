#!/usr/bin/env python3
"""Check metadata/contrasts.tsv against metadata/samples.tsv (cold).

Checks: columns; class vocabulary; sample ids exist, are included and belong to the contrast's
study; numerator and denominator are disjoint and n counts match; each side is one condition group
(a panel side: one group per inbred); sign convention (treatment in the numerator, a control at
the control temperature in the denominator, the numerator colder unless it is a recovery group);
held-constant factors really are equal on both sides; times match unless tagged not_time_matched;
no two pooled contrasts of a study use the same samples. Also lists included samples that no
contrast uses, for review. Usage: python3 metadata/build/validate_contrasts.py
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_contrasts import COLS, GROUP_KEY, ROOT, is_control, num  # noqa: E402

CLASSES = {"chilling", "freezing", "recovery", "germination", "chronic", "field_cold"}


def main():
    samples = {s["sample_id"]: s for s in csv.DictReader(
        open(ROOT / "metadata" / "samples.tsv"), delimiter="\t")}
    rows = list(csv.reader(open(ROOT / "metadata" / "contrasts.tsv"), delimiter="\t"))
    problems = []
    if rows[0] != COLS:
        problems.append(f"header differs: {rows[0]}")
    contrasts = [dict(zip(rows[0], r)) for r in rows[1:]]
    used = set()
    seen_sets = defaultdict(dict)
    for c in contrasts:
        cid = c["contrast_id"]
        tags = set(filter(None, c["tags"].split(";")))
        if c["contrast_class"] not in CLASSES:
            problems.append(f"{cid}: class {c['contrast_class']}")
        ns = c["numerator_samples"].split(";")
        ds = c["denominator_samples"].split(";")
        if set(ns) & set(ds):
            problems.append(f"{cid}: numerator and denominator overlap")
        if len(ns) != int(c["n_numerator"]) or len(ds) != int(c["n_denominator"]):
            problems.append(f"{cid}: n counts do not match sample lists")
        if (c["replicated"] == "TRUE") != (len(ns) >= 2 and len(ds) >= 2):
            problems.append(f"{cid}: replicated flag does not match n")
        side = {}
        for name, ids in (("num", ns), ("den", ds)):
            for i in ids:
                s = samples.get(i)
                if s is None:
                    problems.append(f"{cid}: unknown sample {i}")
                    continue
                if s["include"] != "TRUE":
                    problems.append(f"{cid}: excluded sample {i}")
                if s["study_id"] != c["study_id"]:
                    problems.append(f"{cid}: sample {i} belongs to {s['study_id']}")
                used.add(i)
            side[name] = [samples[i] for i in ids if i in samples]
        if not side.get("num") or not side.get("den"):
            continue
        # each side is one condition group (a panel side: one group per inbred, same design otherwise)
        for name in ("num", "den"):
            keys = [k for k in GROUP_KEY if not ("panel" in tags and k == "genotype")]
            if len({tuple(s[k] for k in keys) for s in side[name]}) > 1:
                problems.append(f"{cid}: {name} samples are not one condition group")
        a, b = side["num"][0], side["den"][0]
        # held constant
        held = ["tissue_detail"] if "tissue_mismatch" not in tags else []
        if "panel" not in tags:
            held.append("genotype")
        if not tags & {"chilling_site", "combined_treatment"}:
            held.append("other_factor_level")
        for f in held:
            if a[f] != b[f]:
                problems.append(f"{cid}: {f} differs ({a[f]!r} vs {b[f]!r})")
        if "panel" in tags and sorted(s["genotype"] for s in side["num"]) != sorted(s["genotype"] for s in side["den"]):
            problems.append(f"{cid}: panel sides do not hold the same inbreds")
        if a["time_h"] != b["time_h"] and b["time_h"] != "" and "not_time_matched" not in tags:
            problems.append(f"{cid}: time differs ({a['time_h']} vs {b['time_h']}) without tag")
        # sign convention: treatment / control
        if is_control(a) or not is_control(b):
            problems.append(f"{cid}: numerator must be a treatment group and denominator a control")
        if a["phase"] != "recovery" and not num(a["temperature_C"]) < num(b["temperature_C"]):
            problems.append(f"{cid}: numerator {a['temperature_C']} °C not colder than denominator {b['temperature_C']} °C")
        if (c["contrast_class"] == "recovery") != (a["phase"] == "recovery"):
            problems.append(f"{cid}: class {c['contrast_class']} vs numerator phase {a['phase']}")
        if c["pool"] == "TRUE":
            key = (frozenset(ns), frozenset(ds))
            if key in seen_sets[c["study_id"]]:
                problems.append(f"{cid}: same samples as {seen_sets[c['study_id']][key]}")
            seen_sets[c["study_id"]][key] = cid

    unused = defaultdict(list)
    for i, s in samples.items():
        if s["include"] == "TRUE" and i not in used:
            unused[s["study_id"]].append(i)
    print(f"{len(contrasts)} contrasts; {len(problems)} problems")
    for p in problems:
        print("  -", p)
    print(f"included samples used by no contrast: {sum(map(len, unused.values()))}")
    for st, ids in sorted(unused.items()):
        print(f"  {st} ({len(ids)}): {', '.join(sorted(ids))[:200]}")
    return problems


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
