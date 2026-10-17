import datetime
import pathlib

PATH = pathlib.Path(r"C:\Users\natha\memory\DECISIONS.md")

ENTRY = """
## 2026-10-18 -- pineal-godseye: Sections 2-3 executed; gene-loss hypothesis SCOPED, walrus flag unverified

### Section 2 -- knowledge landscape written
`projects/pineal-godseye/KNOWLEDGE_LANDSCAPE.md` (4 parts, ~14.8 KB).
Definition (4 competing framings: endocrine gland / neuroendocrine transducer /
CVO / pineal complex), history (Galen -> Vesalius -> Descartes -> Lerner 1958 ->
Axelrod pathway -> SCN mapping -> Nichols 2018 -> gene-loss genomics 2019-21 ->
Yale organoids 2024), hierarchical map Levels 1-6, coverage map.
The calcite/piezoelectric/EMF claim is recorded as **unsupported by
peer-reviewed evidence** -- flagged, not repeated.

### Section 3 -- tested Huelsmann 2019's hypothesis against our own data
Huelsmann (Sci Adv 5:eaaw6671) proposed melatonin-gene loss "may have been a
precondition for the evolution of unihemispheric sleep." Tested on the 25-taxon
phylogeny via `tools/test_gene_loss_hypothesis.py` +
`tools/validate_usws_flag.py`.

**RESULT: the hypothesis is SCOPED, not falsified.** Unihemispheric sleep
occurs in taxa with a FULLY INTACT melatonin pathway (birds: chicken, house
sparrow -- textbook USWS with contralateral eye closure). So melatonin loss is
NOT a precondition for the trait across vertebrates; it is at most a
precondition for the mammalian aquatic forms. Two independent studies of
cetaceans agree gene loss and unihemispheric sleep co-occur THERE.

**Robust to the one bad flag:** excluding the unverified walrus row leaves
3 intact / 3 lost -- the association does not depend on it.

### Third dataset error found (audit pattern continues)
`unihemispheric_sleep=1` for **Walrus** is **not supported**. Primary walrus
sleep work (Pryaslova et al. 2022, PMC8740614) documents QS and REM, not
EEG-confirmed unihemispheric sleep. Lyber & Bovet 2002 state mammalian USWS is
"restricted to aquatic species (Cetacea, eared seals, manatees)"; the walrus
flag rests on eye-closure in ONE captive young animal. NOT corrected in the
data -- awaiting a decision, because the honest value may be 0/1/unknown rather
than binary, and the schema has no unknown.

Also: **the alligator flag is weaker than it looks.** Kelly et al. 2015 J Exp
Biol 218:3175 is titled "Unihemispheric sleep in crocodilians?" -- it
demonstrates unilateral eye CLOSURE, not EEG-confirmed hemispheric asymmetry.
Emerling 2021 confirms crocodylians retain all four melatonin genes, so this
row is intact-pathway + USWS-flagged.

### Standing lesson (4th time this project)
Every value in this dataset that I have checked against primary literature has
needed checking. 3 of 3 audit rounds found an error. The `unihemispheric_sleep`
column is a hand-assigned literature flag, NOT a measured trait -- it is the
weakest evidence in the whole file and should be treated as such. Do not run
statistics on it without reading the flag's provenance first.

### Tooling note
The `write` tool fails with a provider error above roughly 800 words of
content. Large documents must be written as multiple smaller files and
concatenated. Not a content problem -- a size ceiling.
"""


def main():
    raw = PATH.read_bytes()
    body = ENTRY if raw.endswith(b"\n") else "\n" + ENTRY
    PATH.write_bytes(raw + body.encode("utf-8"))
    after = PATH.read_bytes()
    assert after.count(b"\r\n") == raw.count(b"\r\n"), "CRLF count changed!"
    print("appended; CRLF count unchanged:", after.count(b"\r\n"))


if __name__ == "__main__":
    main()
