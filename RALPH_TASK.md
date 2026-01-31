---
task: Canonical Repeat Units (Catalog Loci)
test_command: "cargo test --features python && cd python && python -m pytest -W error"
---

# Task: Canonical Repeat Units (Catalog Loci)

Make **catalog loci/clusters** the canonical unit for repeats. Integrated-callset SV records that overlap repeats contribute to the repeat unit's presence/zygosity and proxy allele representation. TRGT contributes an alternate allele representation that can override proxy by QC policy. Do **not** emit `repeat_<trid>` as a separate default-tested unit; keep it as a representation (and optionally as a shadow unit only in comparison mode). SV units remain one-per-SV only for SVs **outside** repeat catalog regions.

---

## Design

1. **Canonical repeat unit = catalog locus**
   - The single canonical "repeat" test unit for a region is the **catalog locus** (from BED/JSON), not a separate SV unit or a separate TRGT locus unit.
   - One test unit per catalog locus that is either:
     - overlapped by ≥1 SV (so we can build proxy alleles), and/or
     - covered by TRGT (so we have true repeat alleles).
   - Identify catalog loci as usual (BED + optional JSON). Optionally support locus clusters if the catalog defines them.

2. **SVs overlapping repeats → contribute to the repeat unit only**
   - Integrated-callset SV records that overlap a catalog repeat locus **do not** get their own `sv_<id>` test unit.
   - They contribute only to the **repeat (catalog) unit** for that locus:
     - **Presence/zygosity:** SV presence/GT feeds into whether the repeat unit is "present" and zygosity for that sample.
     - **Proxy allele representation:** SV-based proxy (e.g. INS/DEL length) is the allele representation when TRGT is absent or overridden by QC.
   - SVs inside catalog repeat regions → no separate Sv unit; they feed the single repeat (catalog) unit.

3. **SVs outside repeat catalog regions**
   - SV records that do **not** overlap any catalog repeat locus remain **one-per-SV** test units (`sv_<id>`) as today.

4. **TRGT as alternate representation; override by QC**
   - TRGT provides an **alternate allele representation** for a catalog locus (true repeat lengths).
   - When both proxy (SV-derived) and TRGT data exist for the same catalog locus and sample, a **QC policy** decides which representation to use (e.g. prefer TRGT when pass, else fall back to proxy).
   - TRGT is stored as a **representation** on the canonical (catalog) repeat unit, not as a separate default-tested unit.

5. **Do not emit `repeat_<trid>` as a default-tested unit**
   - **Default:** Do **not** create a separate test unit `repeat_<trid>` for each TRGT locus. TRGT data is merged into the canonical catalog-locus unit (matched by locus id / position / catalog id as appropriate).
   - **Optional comparison mode:** Allow a "shadow" or secondary unit per TRGT locus **only** when a comparison mode is enabled (e.g. to compare catalog-unit results vs raw TRGT-unit results). This is not the default; it is for analysis/validation only.

---

## Success Criteria

- [x] Catalog loci are the only repeat test units in the default output (no separate `repeat_<trid>` by default).
- [x] SVs overlapping a catalog repeat locus do not get an `sv_<id>` unit; they only contribute to that repeat unit's presence/zygosity and proxy allele representation.
- [x] SVs outside all catalog repeat regions still get one `sv_<id>` unit each.
- [x] TRGT data is merged into the canonical catalog-locus unit (matching by catalog id / locus); QC policy can prefer TRGT over proxy when both exist.
- [x] Optional: comparison mode can emit shadow `repeat_<trid>` units for validation/comparison; default is off.
- [x] Tests and notebook updated to reflect canonical repeat units and SV-only-for-non-repeat regions.

---

## Context

### Current Behavior

- **Sv:** One test unit per SV record (`sv_<id>`), regardless of overlap with catalog or TRGT.
- **RepeatProxy:** One test unit per catalog locus that overlaps ≥1 SV (`proxy_<catalog_id>`); genotypes from SVs via Resolver.
- **TrueRepeat:** One test unit per TRGT locus (`repeat_<trid>`); genotypes from TRGT.
- The same genomic region can thus have three units (Sv + RepeatProxy + TrueRepeat). No merging; no "canonical" repeat unit.

### Desired Behavior

- **Repeat (canonical):** One test unit per catalog locus that has SV overlap and/or TRGT coverage. Genotypes combine proxy (from SVs) and TRGT; QC policy chooses which allele representation to use when both exist. No separate `repeat_<trid>` by default.
- **Sv:** One test unit per SV **only** for SVs that do **not** overlap any catalog repeat locus. SVs inside repeat regions contribute only to the repeat unit.
- **Optional:** Comparison mode can emit shadow `repeat_<trid>` units for validation.

### Key Files

| File | Change |
|------|--------|
| `src/lib.rs` | Build logic: catalog-first repeat units; SV units only for non-overlapping SVs; merge TRGT into catalog units; no default TrueRepeat units. |
| `src/resolver.rs` | Resolve genotypes for catalog units (proxy from SVs; merge or override with TRGT per QC policy). |
| `src/testunit.rs` | Unit types / IDs for canonical repeat vs shadow; metadata for representation source. |
| `src/cache/` | Schema for genotypes (allele representation source; QC); provenance for canonical vs shadow. |
| `notebooks/storm_demo.ipynb` | Stats and text for canonical repeat units; optional comparison-mode note. |
| Tests | Fixtures and expectations for canonical units, SV-only-outside-repeat, TRGT merged into catalog. |

---

## Implementation Plan

### Phase 1: Catalog-first repeat units

1. **Define canonical repeat units**
   - From catalog (BED + optional JSON), build the set of catalog loci (id, chrom, start, end, motif, etc.).
   - For each catalog locus, determine:
     - Which SVs overlap it (existing `map_svs_to_catalog` / overlap logic).
     - Which TRGT loci match it (by catalog id / LocusId or by position overlap; define matching rule).
   - Emit **one** test unit per catalog locus that has at least one overlapping SV **or** TRGT coverage (or both). Unit id e.g. `repeat_<catalog_id>` or keep `proxy_<catalog_id>` naming; document as the canonical repeat unit.

2. **Genotypes for canonical repeat units**
   - **Proxy (SV-derived):** For samples with overlapping SV calls, compute presence/zygosity and proxy allele (e.g. INS/DEL length) using existing Resolver logic.
   - **TRGT:** For samples with TRGT genotypes at the matching locus, attach true repeat allele representation.
   - **QC policy:** When both proxy and TRGT exist for a sample at a locus, apply policy (e.g. prefer TRGT if pass, else proxy). Store chosen representation and optionally store both in genotype row for downstream comparison.
   - Persist in cache: per (unit, sample) one chosen allele representation plus optional source flag (proxy vs TRGT vs merged).

### Phase 2: SV units only outside repeat regions

3. **Filter SVs for standalone Sv units**
   - After building canonical repeat units, identify SVs that overlap **any** catalog repeat locus.
   - Emit `sv_<id>` test units **only** for SVs that do **not** overlap any catalog locus. SVs that overlap repeats are not emitted as separate Sv units; they already contributed to the repeat unit.

### Phase 3: No default TrueRepeat units; optional shadow

4. **Stop emitting `repeat_<trid>` by default**
   - Remove (or gate) the loop that creates a test unit per TRGT locus. TRGT data is only used to fill the canonical catalog-locus units (matching TRGT locus to catalog locus by id/position).

5. **Optional comparison mode**
   - Add a build/cache option (e.g. `emit_trgt_shadow_units: bool` or `comparison_mode: bool`). When enabled, additionally emit shadow test units `repeat_<trid>` with TRGT-only genotypes so users can compare association results: canonical repeat unit vs raw TRGT unit. Default: off.

### Phase 4: Tests and docs

6. **Tests**
   - Update or add tests: cache build with catalog + SV + TRGT produces canonical repeat units only (no standalone `repeat_<trid>`), SV units only for non-overlapping SVs, and TRGT merged into repeat unit genotypes. Add test for comparison mode emitting shadow units if implemented.
   - Update fixtures if needed (catalog + SV + TRGT overlap scenarios).

7. **Notebook and docs**
   - Update combination stats / "test units by source" to describe canonical repeat units and SV-only-outside-repeat. Remove or relabel "TrueRepeat" for default output; add note on optional comparison (shadow) mode if implemented.
   - README or RALPH_TASK context: short description of canonical repeat units and QC policy.

---

## Acceptance Criteria

When complete:

- Canonical repeat units are the only repeat units in the default cache (one per catalog locus with SV and/or TRGT).
- SVs overlapping a catalog repeat do not have their own `sv_<id>` unit.
- SVs outside catalog repeat regions have one `sv_<id>` unit each.
- TRGT data is merged into canonical repeat units; QC policy selects proxy vs TRGT when both exist.
- No `repeat_<trid>` units in default output; optional shadow units only when comparison mode is on.
- `cargo test --features python` and `cd python && python -m pytest -W error` pass.
- Notebook and docs reflect the new design.

---

## Completion Condition

When:

- All success-criteria checkboxes above are checked.
- `cargo test --features python` passes.
- `cd python && python -m pytest -W error` passes.
- Cache build produces canonical repeat units and SV-only-outside-repeat; TRGT merged into repeat units; no default TrueRepeat units.

Write **DONE** in `.ralph/progress.md` and stop.
