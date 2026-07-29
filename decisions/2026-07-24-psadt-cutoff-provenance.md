# The 10-month doubling-time cutoff is data-derived, not a fixed convention

**Date:** 2026-07-24
**Author of this note:** Travis (corresponding), for Josephine (first author) to review
**Target journal:** European Urology
**Status:** implemented on branch `two-critique-prose`

This note records a reference-review finding that changed one argument in the
Discussion, and the reasoning behind the fix. It matters mostly as a defensive point:
the original wording was a soft target for a knowledgeable reviewer.

## What the earlier draft claimed

The Discussion listed cautions bounding what the literature supports. The second one
said the cutpoint-bias result (data-derived cutpoints inflate false positives) does not
apply to our endpoints "because neither the 50% decline nor the 10-month doubling time
is chosen from the data; both are fixed conventions."

The conclusion is right. The stated reason is wrong for the doubling-time half, and
wrong in a way that our own citation exposes.

## What the evidence actually shows

Pound 1999, which we cite, produced the 10-month value by recursive partitioning:
it searched its own cohort for the split that best predicted metastasis. That is the
"optimal cutpoint" procedure the dichotomization literature warns against. So the
10-month threshold has a data-derived origin, not a round-number one like PSA50's 50%.

It also never settled. The threshold moved into the castration-resistant drug trials
through the denosumab-147 program (Smith, JCO 2013: high-risk enrollment at PSA >=8
and/or PSADT <=10 months), and SPARTAN, PROSPER, and ARAMIS adopted <=10 months.
But EMBARK (Freedland, NEJM 2023), the trial that carried this drug class into high-risk
biochemical recurrence, uses <=9 months. Cohort studies keep reporting still other
"optimal" values (D'Amico's <3 months; Freedland 2005's <3 / 3-8.9 / 9-14.9 / >=15 bins,
built by combining categories with similar hazard ratios; a 2024 Japanese survival-tree
analysis landing at 2.85 and 4.65 months and arguing 10 is not optimal for its population).

The negative evidence on doubling time (Vickers 2009's review that it adds little beyond
absolute PSA) is in the untreated/pretreatment setting, which the paper already excludes.
So none of this weakens the paper. It sharpens it.

## The distinction that holds

What protects a registration trial from cutpoint-induced type-I-error inflation is not
that the threshold is arbitrary. It is that the threshold is pre-specified rather than
re-optimized against that trial's own outcome data. A fixed-in-advance cutpoint, whatever
its historical pedigree, does not dredge the current dataset. That is the honest version
of the argument, and it survives a reviewer who knows where 10 months came from.

The instability of the "optimal" value across populations and endpoints then becomes a
bonus: it is itself a reason to analyze the continuous kinetic rather than any single
threshold of it, which is what the paper recommends.

## What changed

- Rewrote the second caution in the Discussion to say the above, citing Pound 1999 for
  the recursive-partitioning origin and contrasting the <=10-month CRPC trials with
  EMBARK's <=9 months.
- Added `freedland_2023_embark` to `references.bib` (NEJM 2023;389(16):1453-1465,
  doi 10.1056/NEJMoa2303974), metadata verified against Crossref.
- Added citations that were missing from the sentence: Royston 2006 for the
  optimal-cutpoint result itself, and D'Amico 2003 (<3 months) and Freedland 2005
  (binned) for the "cohort studies report still other values" clause, which had been
  asserted without support.

## Still yours to decide

- Whether to name the denosumab-147 program and the Japanese survival-tree cutoff in the
  text, or leave them as the unnamed "cohort studies reporting still other optimal
  values." The current wording keeps the sentence short; naming them adds rigor at the
  cost of length.
- EMBARK is a hormone-sensitive biochemical-recurrence trial, not castration-resistant.
  If a reviewer reads the progression-endpoint section as CRPC-specific, it may be worth
  one clause noting the doubling-time threshold spans both settings.

## How this was checked

Reference metadata (all 21 original citations plus EMBARK) verified against Crossref;
cutoff history traced through PubMed and the trial primary publications; claims read
against source text, not recalled. The provenance point turns on Pound 1999's own methods
section, which we already hold in `resources/`.
