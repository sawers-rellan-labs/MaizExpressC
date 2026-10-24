#!/usr/bin/env python3
"""Build metadata/contrasts.tsv from metadata/samples.tsv (cold).

Each contrast compares two condition groups of one study (samples sharing every design field
except the replicate). Sign convention (HANDOFF §3): log2FC = log2(cold / control), positive =
higher under cold; recovery contrasts put the recovered group in the numerator.

A group is a control when it sits at the study's control temperature (phase control, or chronic
at the control temperature) and is not a recovery group; every other group is a treatment.
Generic rule: inside each stratum of genotype, tissue and other factor level, each treatment
group is paired with
  1. the control at the same time_h (time-matched), else
  2. the control without time_h (a single control, taken as matched), else
  3. the time-0 control (tag not_time_matched; kept in pooling, as in MaizExpressN).
Classes, first match: recovery (phase recovery); germination (seed or germinating tissue, or
germination from sowing); chronic (grown from sowing at the cool temperature); freezing (below
0 °C); chilling (everything else, including 0 °C root chilling in iced water).
Tags added here: not_time_matched, unreplicated / single_library_arm, mild_cold (cold >= 15 °C),
cool_control (control < 20 °C); study rules add the rest. pool = FALSE marks display-only designs.
The rules() table gives each study's notes, tags, blocking factor and any special design.
Usage: python3 metadata/build/build_contrasts.py
"""
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "metadata" / "contrasts.tsv"
COLS = ["contrast_id", "study_id", "contrast_class", "numerator_group", "denominator_group",
        "held_constant", "temp_numerator", "temp_denominator", "cold_C", "control_C", "time_h",
        "recovery_h", "replicated", "n_numerator", "n_denominator", "sign_note", "notes",
        "numerator_samples", "denominator_samples", "genotype", "tissue_group", "tissue_detail",
        "other_factor_level", "tags", "pool", "batch_term"]
# replicates = samples sharing every design field (same key as library_exclusions.py)
GROUP_KEY = ["study_id", "genotype", "tissue_group", "tissue_detail", "temp_label_paper", "temperature_C",
             "temperature_night_C", "control_temperature_C", "phase", "time_h", "recovery_h", "light_condition",
             "other_factor", "other_factor_level"]
BY_DEFAULT = ("genotype", "tissue_detail", "other_factor_level")
GERMINATION = {"Li2022d"}            # whole seedlings germinated from sowing at 10 vs 25 °C
MILD_C, COOL_CONTROL_C = 15, 20
# short tissue names for ids where a study has several tissue details
TISSUE_ID = {
    "Crow2020": [("S2 leaf, 5 cm of leaf base", "V3S2LeafBase"), ("S2 leaf, 5 cm of leaf tip", "V3S2LeafTip"),
                 ("5 cm of leaf base", "V3LeafBase"), ("V1 pooled leaf", "V1Leaf"), ("V3 leaf tissue", "V3Leaf"),
                 ("leaf sheath", "V3Sheath"), ("primary root (V1)", "V1Root"), ("primary root (V3)", "V3Root"),
                 ("stem and shoot apical", "SAM")],
    "Xue2021": [("top crown", "TopCrown"), ("bottom crown", "BottomCrown")],
    "Zhou2024": [("cortex", "Root-Cortex"), ("stele", "Root-Stele"), ("elongation", "Root-EZ"),
                 ("meristematic", "Root-MZ"), ("whole primary root", "Root-Whole")],
}
OTHER_ID = {"greenhouse": "GH", "whole": "wholePlant", "shoot": "shootOnly", "root": "rootOnly"}
UNPAIRED = []


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def load_groups():
    groups = defaultdict(list)
    for s in csv.DictReader(open(ROOT / "metadata" / "samples.tsv"), delimiter="\t"):
        if s["include"] == "TRUE":
            groups[tuple(s[k] for k in GROUP_KEY)].append(s)
    out = defaultdict(list)
    for key, samples in groups.items():
        g = dict(zip(GROUP_KEY, key))
        g["samples"] = sorted(samples, key=lambda s: int(s["replicate"]))
        out[g["study_id"]].append(g)
    return out


def is_control(g):
    return g["phase"] != "recovery" and num(g["temperature_C"]) == num(g["control_temperature_C"])


def klass(g):
    if g["phase"] == "recovery":
        return "recovery"
    if g["tissue_group"] == "Seed/germination" or g["study_id"] in GERMINATION:
        return "germination"
    if g["phase"] == "chronic":
        return "chronic"
    return "freezing" if num(g["temperature_C"]) < 0 else "chilling"


def temp_text(g):
    night = f"/{g['temperature_night_C']}" if g["temperature_night_C"] not in ("", "NR") else ""
    return f"{g['temperature_C']}{night} °C"


def fmt_h(x):
    return ("%g" % float(x)) if num(x) is not None else ""


def id_time(t):
    """Time label for ids: minutes below 1 h, else hours."""
    t = num(t)
    if t is None:
        return ""
    return f"{round(t * 60)}min" if t < 1 else f"{t:g}h"


def id_temp(t):
    return ("minus%gC" % -t) if t < 0 else "%gC" % t


def id_geno(g):
    return re.sub(r"[^A-Za-z0-9.-]", "", re.sub(r"\s*\(.*$", "", g))


def id_tissue(study, g):
    for key, short in TISSUE_ID.get(study, []):
        if key in g["tissue_detail"]:
            return short
    return g["tissue_group"].replace("Whole seedling", "Seedling").replace("Seed/germination", "Seed")


def id_other(x):
    if not x:
        return ""
    w = x.split()[0]
    if w.isupper():                    # HY, TC, NM (Zhou2022 experiments)
        return w
    return OTHER_ID.get(w.lower(), "".join(v[:1].upper() + v[1:] for v in re.findall(r"[A-Za-z0-9]+", x)[:2]))


def label(gs):
    g = gs[0]
    geno = g["genotype"] if len({x["genotype"] for x in gs}) == 1 else f"{len(gs)} inbreds"
    parts = [geno, g["tissue_group"], g["temp_label_paper"]]
    if g["time_h"] not in ("", "NR") and fmt_h(g["time_h"]) not in g["temp_label_paper"]:
        parts.append(f"{fmt_h(g['time_h'])} h")
    if g["other_factor_level"]:
        parts.append(g["other_factor_level"])
    return " | ".join(parts)


def one(vals):
    vals = list(dict.fromkeys(vals))
    return vals[0] if len(vals) == 1 else " / ".join(vals)


def make(study, num_g, den_g, *, cls=None, tags=(), notes="", pool=True, batch="", cold_C=None,
         sign=None, held=None, idp=None):
    """One contrast; num_g and den_g are a group or a list of groups (a panel side)."""
    ng = num_g if isinstance(num_g, list) else [num_g]
    dg = den_g if isinstance(den_g, list) else [den_g]
    a, b = ng[0], dg[0]
    ns = [s["sample_id"] for g in ng for s in g["samples"]]
    ds = [s["sample_id"] for g in dg for s in g["samples"]]
    cls = cls or klass(a)
    cold = cold_C if cold_C is not None else num(a["temperature_C"])
    ctrl = num(b["temperature_C"])
    tags = list(tags)
    if b["time_h"] == "0" and a["time_h"] not in ("0", ""):
        tags.append("not_time_matched")
    if len(ns) < 2 or len(ds) < 2:
        tags.append("unreplicated" if (len(ns) < 2 and len(ds) < 2) else "single_library_arm")
    if sign is None:
        what = "after recovery from cold" if cls == "recovery" else "under cold"
        sign = f"log2FC = log2({a['temp_label_paper']} / {b['temp_label_paper']}); positive = higher {what}"
    if held is None:
        held = {"genotype": one(g["genotype"] for g in ng + dg), "tissue": one(g["tissue_detail"] for g in ng + dg)}
        if a["other_factor_level"] and a["other_factor_level"] == b["other_factor_level"]:
            held["other"] = a["other_factor_level"]
        if a["time_h"] == b["time_h"] and a["time_h"] not in ("", "NR"):
            held["time_h"] = fmt_h(a["time_h"])
    t = a["time_h"] if a["time_h"] != "NR" else ""
    ip = {"geno": id_geno(a["genotype"]), "tissue": id_tissue(study, a), "other": id_other(a["other_factor_level"]),
          "temp": id_temp(cold), "suffix": "",
          "time": ("rec" + id_time(a["recovery_h"])) if cls == "recovery" and a["recovery_h"] else id_time(t)}
    ip.update(idp or {})
    return {
        "contrast_id": None, "_id": ip,
        "study_id": study, "contrast_class": cls,
        "numerator_group": label(ng), "denominator_group": label(dg),
        "held_constant": "; ".join(f"{k}={v}" for k, v in held.items() if v != ""),
        "temp_numerator": temp_text(a), "temp_denominator": temp_text(b),
        "cold_C": "%g" % cold, "control_C": "%g" % ctrl,
        "time_h": fmt_h(t), "recovery_h": a["recovery_h"],
        "replicated": "TRUE" if len(ns) >= 2 and len(ds) >= 2 else "FALSE",
        "n_numerator": len(ns), "n_denominator": len(ds),
        "sign_note": sign, "notes": notes,
        "numerator_samples": ";".join(ns), "denominator_samples": ";".join(ds),
        "genotype": one(g["genotype"] for g in ng + dg) if len(ng) == 1 else f"panel of {len(ng)} inbreds",
        "tissue_group": a["tissue_group"],
        "tissue_detail": one(g["tissue_detail"] for g in ng + dg),
        "other_factor_level": one(g["other_factor_level"] for g in ng + dg),
        "tags": ";".join(dict.fromkeys(tags)), "pool": "TRUE" if pool else "FALSE",
        "batch_term": batch,
    }


def pair_control(t, ctrls):
    """Control for treatment group t: time-matched, else untimed, else time 0."""
    for test in (lambda c: t["time_h"] != "" and c["time_h"] == t["time_h"],
                 lambda c: c["time_h"] == "",
                 lambda c: c["time_h"] == "0"):
        hit = [c for c in ctrls if test(c)]
        if hit:
            if len(hit) > 1:
                sys.exit(f"{t['study_id']}: several controls for {label([t])}")
            return hit[0]
    return None


def generic(study, groups, by=BY_DEFAULT, keep=lambda g: True, **kw):
    strata = defaultdict(list)
    for g in groups:
        if keep(g):
            strata[tuple(g[b] for b in by)].append(g)
    out = []
    for gs in strata.values():
        ctrls = [g for g in gs if is_control(g)]
        for t in gs:
            if is_control(t):
                continue
            c = pair_control(t, ctrls)
            if c is None:
                UNPAIRED.append(f"{study}: {label([t])}")
                continue
            out.append(make(study, t, c, **kw))
    return out


def pick(groups, **kv):
    return [g for g in groups if all(re.fullmatch(v, g[k]) for k, v in kv.items())]


def adjust(cs, when, **changes):
    """Change pool, cold_C, batch, tags (added) or notes (appended) of the contrasts where when(c) holds."""
    for c in cs:
        if not when(c):
            continue
        if "pool" in changes:
            c["pool"] = "TRUE" if changes["pool"] else "FALSE"
        if "cold_C" in changes:
            c["cold_C"] = "%g" % changes["cold_C"]
            c["_id"]["temp"] = id_temp(changes["cold_C"])
        if "batch" in changes:
            c["batch_term"] = changes["batch"]
        if "tags" in changes:
            c["tags"] = ";".join(dict.fromkeys(filter(None, c["tags"].split(";") + list(changes["tags"]))))
        if "drop_tags" in changes:
            c["tags"] = ";".join(x for x in c["tags"].split(";") if x not in changes["drop_tags"])
        if "notes" in changes:
            c["notes"] = changes["notes"]
    return cs


def is_recovery(c):
    return c["contrast_class"] == "recovery"


# ---------------------------------------------------------------- per-study rules
def rules():
    R = {}

    def simple(study, **kw):
        R[study] = lambda gs: generic(study, gs, **kw)

    # acute chilling with a matched or single control
    for st, note in [
        ("Jiang2022", "4 °C for 12 h vs untreated (25 °C), per genotype."),
        ("Liao2026", "4 °C for 12 h vs 25 °C, per genotype; which leaf NR."),
        ("Li2017", "4 °C for 2 h vs the 28/22 °C control (n = 2 vs 2). ABA-treated samples are not a cold contrast."),
        ("Li2020a", "Three chilling temperatures (16, 10 and 4 °C) for 2 h, each vs the same 25 °C control."),
        ("Lu2017", "Teosinte (Zea mays ssp. mexicana); roots, stems and leaves pooled; 4 °C for 12 h vs 25 °C."),
        ("Zhao2021", "Whole roots, 5 °C for 3 d vs 22 °C."),
        ("Sommer2025", "Root hairs of the primary root; 24 h at mild (16/12 °C) or severe (10/6 °C) cold vs 27/22 °C."),
        ("Frey2020", "Primary root tips, 8 h at 12 °C vs 18 °C after 4 d at 26/18 °C; one contrast per inbred "
                     "(tolerant and susceptible lines)."),
        ("Friero2023", "Root tips of germinated seeds, 20 °C vs 30 °C at 6 and 12 h after the shift (time-matched)."),
        ("Yang2025", "10 °C for 7 d vs 25 °C, water-treated plants only; the ABA, BR and ABA + BR groups have no "
                     "matched control and are not used."),
        ("Jia2020", "10 °C for 3 d (seedling10) and 10 °C for 3 d then 2 °C for 1 d (seedling2, cumulative) vs "
                    "28 °C, per genotype."),
        ("Zhou2024", "Primary root zones of PE0075, 16/12 °C vs 26/18 °C; T1 after the 16 h day phase, T2 after the "
                     "following 8 h night (time-matched controls). Laser-capture cell types not used."),
        # time-zero controls (not_time_matched, kept in pooling)
        ("Wang2023", "4 °C for 2 h, 24 h and 7 d vs the 0 h sample at 25 °C (genotype NR; not age-matched at 7 d)."),
        ("Yu2021", "5 °C for 12 and 24 h vs the 0 h control, per inbred."),
        ("Li2025", "4 °C for 24 and 48 h vs the 0 h control, per genotype."),
        ("Zeng2021", "4 °C for 12 h vs the 0 h control. Its LH244 wild-type libraries are also Zeng2024's "
                     "(counted once, here)."),
        ("Zeng2024", "ZmTIP4-OE only: 4 °C for 12 h vs the 0 h control; the wild type is Zeng2021's."),
        ("Li2016", "Freezing: -1 °C for 1 h vs the 0 h control at 25/20 °C. KR701 has one cold library."),
        # germination
        ("Jiang2023", "Germinating seeds after 3 d at 15 °C vs 25 °C (16 h light); dry seeds not used."),
        ("Zhang2020", "Embryos of imbibed seeds, 2 h at 10 °C vs 25 °C after 6 h imbibition at 25 °C."),
        ("Zhang2024", "Radicles of germinated seeds, 5 °C for 48 h vs the 15 °C germination temperature "
                      "(the paper's control)."),
    ]:
        simple(st, notes=note)
    simple("Gao2024", tags=["temp_assumed"],
           notes="Cold vs untreated for 12 h; 4 vs 25 °C assumed from the same lab's studies (Methods not checked).")
    simple("Makarevitch2015", batch="replicate (replicates sequenced in separate batches)", tags=["batch_blocked"],
           notes="Whole above-ground seedlings, 5 °C for 16 h vs 24 °C. B73 and B37 blocked on replicate; "
                 "Oh43 has one library per condition.")
    simple("Crow2020", tags=["stage_matched"],
           notes="Grown from sowing in cool (22/11 °C) vs warm (32/22 °C) chambers; sampled at the same V stage "
                 "(V1 or V3), so cool plants are older. Shallow libraries (2 M mapped-read floor).")
    simple("Li2022d", tags=["stage_matched"],
           notes="Whole seedlings germinated at 10 °C (15 d) vs 25 °C (3 d), matched on developmental stage.")

    R["Avila2018"] = lambda gs: adjust(
        generic("Avila2018", gs, batch="replicate (independent repeats with matched controls)", tags=["batch_blocked"],
                notes="Second leaf blade, 24 h into 3 d at 14/2 °C vs a stage-matched 24/14 °C control."),
        is_recovery, cold_C=14, notes="Second leaf blade, 1 d at 24/14 °C after 3 d at 14/2 °C vs the day-4 control.")

    def li2022c(gs):
        out = generic("Li2022c", gs, notes="4 °C for 12 h; wild type vs the 25 °C 12 h control (time-matched), "
                                           "bzip68 mutant vs its 0 h control.")
        wt = pick(gs, genotype="LH244")
        out.append(make("Li2022c", pick(wt, phase="cold")[0], pick(wt, time_h="0")[0], pool=False,
                        tags=["secondary_control"], idp={"suffix": "vs0h"},
                        notes="Secondary contrast: wild type 4 °C 12 h vs the 0 h control (display only; the "
                              "time-matched 25 °C contrast is the pooled one)."))
        return out
    R["Li2022c"] = li2022c

    R["He2023a"] = lambda gs: adjust(
        generic("He2023a", gs, notes="2 °C for 72 and 120 h vs time-matched 25 °C controls, per genotype."),
        is_recovery, pool=False, tags=["timing_unknown"], cold_C=2,
        notes="Recovery (24 h at 25 °C) vs the plain control; the length of the cold before recovery and the "
              "control's age are not reported, so display only until the authors confirm.")
    R["Meng2022a"] = lambda gs: adjust(
        generic("Meng2022a", gs, tags=["not_age_matched"]), lambda c: True, cold_C=13,
        notes="Germinating seeds after 96 h at 13 °C and 48 h recovery at 25 °C vs seeds kept at 25 °C (not "
              "age-matched; no sample at the end of the cold).")
    R["He2016"] = lambda gs: adjust(
        generic("He2016", gs, notes="Embryos of imbibing seeds, 2 h at 5 °C vs 25 °C; one library per condition."),
        is_recovery, cold_C=5, notes="Embryos of imbibing seeds, 1 h back at 25 °C after 2 h at 5 °C vs the 3 h "
                                     "25 °C control; one library per condition.")

    def meng2025(gs):
        out = []
        for geno, setid in (("JD558", "setA"), ("NR", "setB")):
            mine = pick(gs, genotype=geno)
            for t in [g for g in mine if g["phase"] in ("cold", "recovery")]:
                c = pick(mine, phase="control", time_h=t["time_h"])[0]
                mism = t["tissue_detail"].split()[0] != c["tissue_detail"].split()[0]
                out.append(make("Meng2025", t, c, pool=False, cold_C=4,
                                tags=["combined_treatment"] + (["tissue_mismatch"] if mism else []),
                                held={"genotype": geno, "time_h": fmt_h(t["time_h"])},
                                idp={"geno": setid, "other": ""},
                                notes="Combined cold (4 °C) and high soil moisture (25%) vs 20/15 °C at 15%, "
                                      "time-matched; display only (cold cannot be separated from moisture)."
                                      + (" Control libraries are embryos, stressed ones whole seeds." if mism else "")))
        return out
    R["Meng2025"] = meng2025

    R["Zhang2017"] = lambda gs: adjust(
        generic("Zhang2017", gs, batch="replicate (b1-b3: planting and sequencing batch)", tags=["batch_blocked"],
                notes="All above-ground tissue, 6 °C vs 29/23 °C at 0.5-24 h after cold onset at lights-off "
                      "(time-matched)."),
        lambda c: c["other_factor_level"] != "", pool=False, batch="", tags=["combined_treatment"],
        drop_tags=["batch_blocked"],
        notes="Greenhouse experiment: 1 d at 6 °C in a growth chamber vs plants left in the greenhouse "
              "(29/23 °C, more light), so cold is combined with the move; display only.")

    def xue2021(gs):
        note = "QuantSeq 3' tag libraries. "
        out = generic("Xue2021", gs, keep=lambda g: g["other_factor"] == "", tags=["tag_seq"],
                      notes=note + "Leaves, 24 h at 4 °C vs 26 °C (GSE167156, n = 2 vs 2), per inbred.")
        site = [g for g in gs if g["other_factor"]]
        where = {"4S4R": "shoot and root at 4 °C", "4S26R": "shoot at 4 °C, root at 26 °C",
                 "26S0R": "shoot at 26 °C, root in iced water at 0 °C"}
        for td in sorted({g["tissue_detail"] for g in site}):
            mine = [g for g in site if g["tissue_detail"] == td]
            ctrl = pick(mine, temp_label_paper="26S26R")[0]
            tg = mine[0]["tissue_group"]
            arms = [("4S4R", "whole", True, [])]
            if tg in ("Leaf", "Root"):
                local, distal = ("4S26R", "26S0R") if tg == "Leaf" else ("26S0R", "4S26R")
                arms += [(local, "local", False, ["local_only"]), (distal, "distal", False, ["distal"])]
            for lab, kind, pool, extra in arms:
                t = pick(mine, temp_label_paper=lab)[0]
                sign = None
                if kind == "distal":
                    sign = (f"log2FC = log2({lab} / 26S26R); positive = higher in the {tg.lower()} of plants chilled "
                            f"elsewhere (the sampled {tg.lower()} stayed at 26 °C)")
                what = {"whole": "whole-plant chilling vs 26 °C.",
                        "local": "chilled organ only (local response); display only, the whole-plant contrast is "
                                 "the pooled one.",
                        "distal": "organ not chilled (distal response); display only."}[kind]
                out.append(make("Xue2021", t, ctrl, pool=pool, tags=["tag_seq", "chilling_site"] + extra, sign=sign,
                                held={"genotype": "B73", "tissue": td, "time_h": "24"},
                                notes=note + f"Split-temperature experiment (GSE167881), 24 h, {where[lab]}; {what}"))
        return out
    R["Xue2021"] = xue2021

    def zhou2022(gs):
        exp = lambda g: g["other_factor_level"].split()[0]
        out = generic("Zhou2022", gs, notes="V2 leaves, 6/2 °C vs 30/20 °C with the time-matched control of the "
                                            "same experiment (HY replicated; TC time course and NM panel have one "
                                            "library per genotype and time).")
        hy = [g for g in gs if exp(g) == "HY"]
        for t in [g for g in hy if g["phase"] == "cold"]:
            c0 = [g for g in hy if g["genotype"] == t["genotype"] and g["phase"] == "control" and g["time_h"] == "0"][0]
            out.append(make("Zhou2022", t, c0, pool=False, tags=["secondary_control"], idp={"suffix": "vs0h"},
                            notes="Secondary contrast vs the 0 h control of the HY experiment (display only; the "
                                  "time-matched contrast is the pooled one)."))
        nm = [g for g in gs if exp(g) == "NM"]
        for h in ("1", "25"):
            cold = {g["genotype"]: g for g in nm if g["phase"] == "cold" and g["time_h"] == h}
            ctrl = {g["genotype"]: g for g in nm if g["phase"] == "control" and g["time_h"] == h}
            both = sorted(set(cold) & set(ctrl))
            out.append(make("Zhou2022", [cold[x] for x in both], [ctrl[x] for x in both], tags=["panel"],
                            batch="genotype (paired: one cold and one control library per inbred)",
                            held={"tissue": nm[0]["tissue_detail"], "other": nm[0]["other_factor_level"],
                                  "time_h": h},
                            idp={"geno": "NMpanel"},
                            notes=f"Diversity panel at {h} h: cold vs control across the {len(both)} inbreds with "
                                  "both libraries, paired by inbred (inbreds serve as replicates; model ~ genotype "
                                  "+ temperature within these libraries)."
                                  + (" KI11 has no 1 h cold library." if h == "1" else "")))
        return out
    R["Zhou2022"] = zhou2022
    return R


def finalize(contrasts):
    """Temperature tags from the final cold_C / control_C; no blocking factor on unreplicated contrasts."""
    for c in contrasts:
        tags = [x for x in c["tags"].split(";") if x and x not in ("mild_cold", "cool_control")]
        if float(c["cold_C"]) >= MILD_C:
            tags.append("mild_cold")
        if float(c["control_C"]) < COOL_CONTROL_C:
            tags.append("cool_control")
        if c["replicated"] == "FALSE" and c["batch_term"]:
            c["batch_term"] = ""
            tags = [x for x in tags if x != "batch_blocked"]
        c["tags"] = ";".join(tags)


def assign_ids(contrasts):
    """Ids from study, class, genotype, tissue, other factor, temperature and time; parts that do not vary
    within a study are left out (tissue is always kept, other factor too when it only renames the tissue)."""
    by_study = defaultdict(list)
    for c in contrasts:
        by_study[c["study_id"]].append(c)
    for study, cs in by_study.items():
        vary = {k for k in ("geno", "other", "temp", "time") if len({c["_id"][k] for c in cs}) > 1}
        pairs = {(c["_id"]["tissue"], c["_id"]["other"]) for c in cs}
        if len(pairs) == len({p[0] for p in pairs}) == len({p[1] for p in pairs}):
            vary.discard("other")
        for c in cs:
            ip = c.pop("_id")
            parts = [study, c["contrast_class"], ip["geno"] if "geno" in vary else "", ip["tissue"],
                     ip["other"] if "other" in vary else "", ip["temp"] if "temp" in vary else "",
                     ip["time"] if "time" in vary else "", ip["suffix"]]
            c["contrast_id"] = "_".join(p for p in parts if p)


def main():
    groups = load_groups()
    R = rules()
    contrasts, missing = [], []
    for study in sorted(groups):
        if study not in R:
            missing.append(study)
            continue
        contrasts += R[study](groups[study])
    if missing:
        sys.exit(f"no contrast rule for: {missing}")
    finalize(contrasts)
    assign_ids(contrasts)
    ids = [c["contrast_id"] for c in contrasts]
    dup = sorted({i for i in ids if ids.count(i) > 1})
    if dup:
        sys.exit(f"duplicate contrast ids: {dup}")
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(contrasts)
    print(f"{len(contrasts)} contrasts from {len({c['study_id'] for c in contrasts})} studies -> {OUT}",
          file=sys.stderr)
    for u in UNPAIRED:
        print(f"  no control for {u}", file=sys.stderr)


if __name__ == "__main__":
    main()
