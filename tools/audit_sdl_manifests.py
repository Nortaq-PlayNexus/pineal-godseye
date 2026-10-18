"""Does ScientificDiscoveryLab's hash manifest actually pin anything?

Pineal-godseye's failure mode was hand-transcribed numbers with no receipt. The
fix adopted there was per-cell provenance. ScientificDiscoveryLab appears to
already do something like this -- CANONICAL_HASHES_NON_DESTRUCTIVE.json and
results_manifest.json.

A recorded hash is not a receipt unless someone CHECKS it. This script tests
whether the recorded digests match the files they claim to pin. If they match,
the project has genuine integrity checking and the pineal-godseye pattern does
not recur. If they do not match, or the manifest is a hand-written list with no
verification code, the same weakness is present.

Read-only.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(r"C:\Users\natha\AI_RESEARCH\ScientificDiscoveryLab")


def sha256_file(p, chunk=1 << 20):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(chunk), b""):
            h.update(c)
    return h.hexdigest()


def load(p):
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception as e:
        return {"__error__": str(e)}


def main():
    print("ScientificDiscoveryLab -- do the manifests actually pin anything?\n")

    # 1. CANONICAL_HASHES
    ch = ROOT / "AUDIT" / "CANONICAL_HASHES_NON_DESTRUCTIVE.json"
    if ch.exists():
        d = load(ch)
        print(f"CANONICAL_HASHES_NON_DESTRUCTIVE.json  ({ch.stat().st_size:,} B)")
        print(f"  top-level keys: {list(d)[:8] if isinstance(d, dict) else type(d)}")
        entries = []
        if isinstance(d, dict):
            for k, v in d.items():
                if isinstance(v, dict):
                    for kk, vv in v.items():
                        if isinstance(vv, str) and len(vv) == 64:
                            entries.append((f"{k}/{kk}", vv))
                        elif isinstance(vv, dict):
                            dg = vv.get("sha256") or vv.get("digest") or \
                                 vv.get("hash")
                            fp = vv.get("path") or vv.get("file")
                            if dg and len(str(dg)) == 64:
                                entries.append((fp or f"{k}/{kk}", str(dg)))
                elif isinstance(v, str) and len(v) == 64:
                    entries.append((k, v))
        print(f"  digests found: {len(entries)}")
        checked = matched = mismatched = missing = 0
        bad = []
        for ref, dg in entries:
            cand = ROOT / ref
            if not cand.exists():
                missing += 1
                continue
            got = sha256_file(cand)
            checked += 1
            if got.lower() == dg.lower():
                matched += 1
            else:
                mismatched += 1
                bad.append((ref, dg[:12], got[:12]))
        print(f"  checked={checked} matched={matched} "
              f"MISMATCH={mismatched} path-not-found={missing}")
        for b in bad[:8]:
            print(f"    MISMATCH {b[0]}\n      recorded {b[1]}  actual {b[2]}")
        if checked and mismatched == 0:
            print("  -> Digests genuinely match the files. Real integrity check.")
        elif checked == 0:
            print("  -> No resolvable file references; cannot confirm the")
            print("     digests are checked against anything.")
        print()

    # 2. results_manifest
    rm = ROOT / "AUDIT" / "results_manifest.json"
    if rm.exists():
        d = load(rm)
        print(f"results_manifest.json  ({rm.stat().st_size:,} B)")
        print(f"  top-level: {list(d)[:8] if isinstance(d, dict) else type(d)}")
        print()

    # 3. Is verification code present, or only the manifest files?
    print("VERIFICATION CODE")
    hits = 0
    for py in ROOT.rglob("*.py"):
        try:
            t = py.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue
        if "CANONICAL_HASHES" in t or "canonical_hash" in t.lower():
            rel = py.relative_to(ROOT)
            lines = [ln.strip() for ln in t.splitlines()
                     if "CANONICAL_HASH" in ln or "canonical_hash" in ln.lower()]
            print(f"  {rel}")
            for ln in lines[:3]:
                print(f"      {ln[:100]}")
            hits += 1
        if hits >= 6:
            break
    if hits == 0:
        print("  NONE FOUND. The manifest exists but nothing reads it.")
        print("  -> A hash file nobody opens is documentation, not a check.")
    print()

    print("INTERPRETATION")
    print("-" * 60)
    print("ScientificDiscoveryLab records hashes and pins files. That is the")
    print("structural protection pineal-godseye lacked and later adopted.")
    print("Whether it is ENFORCED depends on whether code compares the recorded")
    print("digest against the file, not merely whether a manifest exists.")


if __name__ == "__main__":
    main()
