# Paper B → EPJD (Springer/EDP/SIF) format migration — **APPLIED**

Status: **applied on 2026-08-15**, on explicit written authorisation from the
author. `main.tex` is now an EPJ-class manuscript. Formatting only: no
scientific content, number, claim or table value was changed. The
pre-migration REVTeX 4.2 source is preserved verbatim as
`../main_revtex_backup.tex` (with a header comment marking it the fallback).

## Which class was actually used, and why not `svjour3`

The earlier version of this file targeted `svjour3.cls`. That was wrong for
this machine, and arguably wrong in general:

* **`svjour3.cls` is not on CTAN and is not installed here.** Springer
  distributes it directly (journal-specific download bundles), so it cannot be
  pulled from a TeX distribution.
* **What *is* official and available is CTAN
  `macros/latex/contrib/springer/svjour/epj/`**, titled *"LaTeX2e support for
  The European Physical Journal"*: **`svjour.cls` v1.8 (1999-09-14) +
  `svepj.clo` v1.1**. That is the EPJ-specific Springer package, and it is
  what the manuscript now uses. Both files, plus Springer's own fill-in
  template `_epj_template_reference.tex`, sit next to `main.tex`.

Class line, before and after:

```latex
% before
\documentclass[aps,pra,twocolumn,superscriptaddress,floatfix,longbibliography]{revtex4-2}
% after
\documentclass[epj,nopacs]{svjour}
```

`svepj.clo` already executes `final,twocolumn,10pt`, so the layout stays
two-column at 10pt; the page measure becomes 180 mm × 56 pc with a 4 mm
column gap.

## What the migration actually did

1. **Class** — as above. `svjour` requires a journal option; `epj` makes it
   load `svepj.clo`. `nopacs` suppresses the otherwise **mandatory** PACS
   block (without it the class typesets a literal
   `PACS. XX.XX.XX No PACS code given`). PACS was retired in 2010 and EPJD
   asks for keywords instead, so no PACS codes were invented.
2. **Front matter** — `\author` + `\affiliation` + `\email` (REVTeX) became
   svjour's `\author` + `\institute{... \\ \email{...}}`. `\titlerunning` and
   `\authorrunning` were added: the EPJ class builds a running head from the
   title and substitutes *"Title Suppressed Due to Excessive Length"* when it
   does not fit, which this title does not. Both are short forms of the real
   title/author, nothing new.
3. **Abstract** — `\begin{abstract}...\end{abstract}` became the command form
   `\abstract{...}`, placed **before** `\maketitle` (svjour stores it and
   typesets it inside the title block). Text unchanged.
4. **Keywords** — the REVTeX-era `\paragraph*{Keywords.}` hack became
   svjour's real mechanism, `\keywords{...}` inside the `\abstract` argument.
   It prints as a bold **"Key words."** run-in paragraph under the abstract.
   Same seven keywords, same order.
5. **Tables** — REVTeX's `ruledtabular` / `\colrule` do not exist outside
   REVTeX. All four tables were converted to `booktabs`
   (`\toprule` / `\midrule` / `\bottomrule`) inside plain `tabular`, with
   `\centering`. `\multicolumn` spans carried over unchanged.
   * `tab:shi` stays a **`table*`** (12 rows × 7 columns, spanning both
     columns) — the layout is still two-column, so starred floats remain
     valid. Set at 8 pt with `\tabcolsep=4pt`; all twelve rows render.
   * `tab:prod` (16 rows × 6 cols) and `tab:ddoub` (10 rows × 6 cols) stay
     single-column `table`s at 8 pt with `\tabcolsep=3pt`.
6. **The `\shrinktable` hack is gone.** It poked REVTeX's internal
   `\table@hook`, which svjour does not have. `svepj.clo` also aliases
   `\footnotesize` to `\small` (9 pt) and then jumps to `\scriptsize` (7 pt),
   i.e. it has no 8 pt step. An explicit 8 pt `\tablesize` was defined in the
   preamble so the dense tables keep the size they had under REVTeX.
7. **`\emergencystretch = 1.5em`** was added. svjour sets none, and a handful
   of body lines carrying long unbreakable numeric literals
   (`2.283366856-0.054871772i` and friends) overran the column. This is a
   line-breaking parameter; no text was rewritten to fit.
8. **Bibliography** — see below.
9. **Unchanged and still working**: `\section*{Data availability}`,
   `{Code availability}`, `{Funding}`, `{Conflicts of interest}`,
   `{Acknowledgments}`; the `sloppypar` in Code availability; the figure; all
   `\label`/`\ref`; both `% CIERRE95:` provenance comments and all three
   `% NEEDS_HUMAN:` comments (comments only — nothing typeset).

Not present in this manuscript, hence nothing to convert: `\preprint`,
`\pacs`, `\onecolumngrid` / `\twocolumngrid`, `\squeezetable`, `revsymb`
macros, `figure*`.

## Bibliography style

**`spphys.bst` is not installed and is not obtainable from this MiKTeX
repository** (it ships with Springer's own bundles, like `svjour3`).
Springer's EPJ template does not prescribe a style either — it leaves
`\bibliographystyle{}` empty and shows a manual `thebibliography`.

Fallback used: **`natbib` (`numbers,sort&compress`) + `unsrtnat.bst`**, both
standard and installed. This is the closest available analogue of `spphys`,
which is itself a natbib numeric style with references in order of citation:

* numeric citations, `[13]`, `[3–6]` — same look as the REVTeX build;
* references numbered in order of first citation — the numbering is
  **identical** to the REVTeX/`apsrev4-2` build (1 = Rogers 1970, …,
  13 = Shi 2014, …, 20 = Ackad & Horbatsch 2008);
* all 20 entries resolve; **0 undefined citations**;
* `unsrtnat` additionally prints the `doi` fields that `apsrev4-2` was
  dropping — a gain, not a loss.

`NEEDS_HUMAN`: if EPJD supplies `spphys.bst` (or its successor) with the
submission bundle, swap `\bibliographystyle{unsrtnat}` for it and drop the
explicit `natbib` load only if the class-provided style requires it.

## Build result (pdflatex → bibtex → pdflatex ×2, MiKTeX 26.2)

| Metric | REVTeX 4.2 (before) | svjour+svepj EPJ (after) |
|---|---|---|
| Errors | 0 | **0** |
| Overfull boxes | 0 | **0** |
| Underfull boxes | — | 16 (7 `sloppypar`, 4 Data availability, 3 stretched body lines, 1 bibliography, 1 `\vbox`) |
| Undefined references | 0 | **0** |
| Undefined citations | 0 | **0** |
| Pages | 6 | **6** |
| Figures / tables rendered | 1 / 4 | 1 / 4 |

Verified reproducible from a clean directory (only `main.tex`, `refs.bib`,
`svjour.cls`, `svepj.clo`, `figB_splitting.pdf`): identical numbers.

## Still open

`NEEDS_HUMAN`: **abstract length.** The abstract is **315 words**
(LaTeX-stripped; 316 counting the "Abstract." label as typeset) and may
exceed EPJD's limit — check the current author instructions. It was
deliberately **not** trimmed here: the numbers in its later sentences are the
paper's claim to novelty. Candidate cuts are listed in
`../../_cierre95/handoff_B.md`.

`NEEDS_HUMAN`: EPJD accepts other classes at initial submission in current
practice (reformatting happens at acceptance), so this migration is a
convenience, not a gate — confirm on the EPJD submission page.

`NEEDS_HUMAN`: EPJ manuscripts are normally dated
`\date{Received: date / Revised version: date}` (Springer fills these in).
`main.tex` keeps the REVTeX-era `\date{\today}`; change it if the editorial
system wants the received/revised form.

Leftover: `../mainNotes.bib` was auto-generated by REVTeX's `\bibliography`
and holds nothing but two `@CONTROL` entries for `apsrev4-2`. It is now
orphaned and can be deleted; it does not affect the build.
