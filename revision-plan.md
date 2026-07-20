# PSA50 Paper — Revision Plan

**Purpose:** self-contained plan to execute the next round of work on this paper. Written to be picked up in a fresh Claude Code session. Combines (a) a statistical/code review of the current draft and simulations and (b) a verified literature foundation on why percentage change from baseline and responder dichotomization are statistically weak.

**Repo:** `github.com/jkaminags/PSA50Paper` (public). Any Issues/PRs go here, **not** the MSKCC enterprise host. As of 2026-07-20 the repo has no open Issues or PRs.

**Author note:** Josephine Kaminaga (first author), Travis Gerke (corresponding). Draft is in [`index.qmd`](index.qmd); simulations in [`simulations/`](simulations/).

---

## 1. The one decision that reframes the whole paper

The draft blends **two distinct statistical critiques** as if they were one. They have different mechanisms, different fixes, and different literature. Separating them is the highest-leverage edit.

1. **Dichotomization** — collapsing a continuous PSA trajectory into a binary "≥50% decline" flag (this *is* PSA50). Loses information and power; a median-type split costs roughly the power of discarding a third of the sample (Royston/Altman/Sauerbrei 2006; Altman & Royston 2006).
2. **Percentage change as a continuous metric** — analyzing `(baseline − follow-up)/baseline` directly. Fails to remove baseline imbalance (regression to the mean), is a ratio so it's skewed/asymmetric and non-normal, is the least efficient of the standard analysis choices, and its power collapses when baseline variance shifts (Vickers 2001; Kaiser 1989).

The current title ("PSA50 is a Statistically Incorrect Endpoint") and Background ([index.qmd:22–24](index.qmd)) slide between these. The cited percent-change literature (Vickers/Kaiser/Bland) is about the *continuous* metric, but PSA50 is the *dichotomized* one.

**Recommended thesis:** *PSA50 discards information twice over — it dichotomizes, and the thing it dichotomizes (percent change) is itself a poor scale — and there are statistically efficient, clinically interpretable alternatives (log-scale ANCOVA / raw change, with per-patient percent-change waterfalls for description).*

**Decision needed from Travis:** confirm this framing before prose is rewritten. It also determines which simulation "alternative hypothesis" design to keep (see §4).

---

## 2. Verified literature foundation

A deep-research pass (25 claims, adversarially verified 3 votes each, 0 refuted) established the following. Use these as the citation backbone. Access notes tell you which PDFs to grab.

### 2a. Core methodological references (the spine of the argument)

| Reference | What it establishes | Access |
|---|---|---|
| **Vickers 2001**, *BMC Med Res Methodol* 1:6 — "The use of percentage change from baseline… is statistically inefficient: a simulation study" | Percent change (FRACTION) has the **lowest power** of POST/CHANGE/FRACTION/ANCOVA and its power **collapses when baseline SD doubles**; ratios have "no analytic reason" to be normal. Recommends: don't analyze percent change; use ANCOVA and back-convert to a percentage using *mean* baseline scores. Power at r=0.2: ANCOVA 72.3% vs FRACTION 45.1%; at r=0.5: FRACTION 67.0% vs ANCOVA 82.3%. | **Open** (PMC34605). Already in [references.bib](references.bib) as `vickers_2001_the`. |
| **Bland & Altman 2011**, *Trials* 12:264 — "Comparisons against baseline within randomised groups… can be highly misleading" | Testing change vs baseline **within** each arm (instead of between arms) inflates Type I error to **~50% (two arms), ~75% (three arms)**. "Significant in one arm, not the other" is not a valid between-group test. Within-arm before/after significance can't establish an effect (time trends + RTM). | **Open** (PMC3286439). Already in bib as `bland_2011_comparisons`. |
| **Vickers & Altman 2001**, *BMJ* 323:1123 — "Statistics Notes: Analysing controlled trials with baseline and follow up measurements" | "ANCOVA is the preferred general approach." Change scores don't remove baseline imbalance (RTM). Efficiency at r=0.6: 85 (follow-up only) vs 68 (change) vs 54 (ANCOVA) patients. | **Open** (PMC1121605); bmj.com bot-blocks direct fetch (403). **Add to bib.** |
| **Kaiser 1989**, *Stat Med* 8(10):1183–1190 — "Adjusting for baseline: change or percentage change?" | **More nuanced than a blanket anti-percent-change paper — read it as a "when is percent change OK" reference.** Neither change nor percent change is automatically correct; pick the one that is *independent of baseline* by plotting each against baseline. The wrong choice costs power: in his simulations you need ~**40–45% more patients** (sample-size ratio ≈ 140–145%). His likelihood-ratio rule (Appendix I) shows percent change is the right scale when the response SD is proportional to baseline (multiplicative error) — the bridge to the log-scale argument. Positions ANCOVA as the general case when the response-on-baseline slope is near neither 0 nor 1. **Read (`resources/kaiser1989.pdf`).** In bib as `kaiser_1989_adjusting`. |
| **Royston, Altman & Sauerbrei 2006**, *Stat Med* 25(1):127–141 — "Dichotomizing continuous predictors in multiple regression: a bad idea" | Exact efficiency figures: median dichotomization of a **normal** predictor → 65% asymptotic efficiency ("effectively equivalent to losing a third of the data"); of an **exponential/right-skewed** predictor → only **48%** (worse). Since PSA is right-skewed, cite the ~48% end. Data-derived "optimal" cutpoints inflate Type I error to ~25–50% and bias effect sizes. | **Paywalled** at Wiley; **open author manuscript** at Oxford ORA (uuid:31fc8902…). **Read (`resources/royston2005.pdf`). Add to bib.** |
| **Altman & Royston 2006**, *BMJ* 332:1080 — "The cost of dichotomisation" | Median-split power loss ≈ discarding one-third of the data. "Deliberately discarding data is surely inadvisable." | **Open** (bmj.com 403s automated fetch). **Add to bib.** |
| **Harrell, BBR**, "Change from baseline" chapter — hbiostat.org/bbr/change | Textbook case against change scores / percent change; ANCOVA preferred. Sharper than the current "seven assumptions" gloss. | **Open.** Already in bib as `harrell_2025_biostatistics` (cite the change chapter specifically). |

> ⚠️ **Two scope caveats to state honestly in the paper:**
> 1. Royston/Altman/Sauerbrei and Altman & Royston concern dichotomizing *predictors/covariates*; PSA50 dichotomizes an *outcome*. The 48–65% efficiency figures transfer by analogy, not direct measurement — present them as illustrative, not exact for PSA50.
> 2. **Do not invoke the "optimal cutpoint" bias critique for PSA50.** Royston's Type I error inflation (25–50%) applies to *data-derived* cutpoints. PSA50's 50%-decline threshold is a *fixed, a priori* convention, so that bias argument doesn't apply. The applicable Royston critique is information/power loss only. Using the wrong argument here is the kind of thing a reviewer will catch.

### 2b. When percent / relative change IS defensible (the honest counterpoint)

This is the "maybe it's okay sometimes" branch Travis asked about. The answer: **yes, on the log scale.**

- **Cole 2000**, *Stat Med* 19:3109–3125 — "Sympercents: symmetric percentage differences on the 100·logₑ scale." Read (`resources/cole2000.pdf`). Two concrete defects of conventional percent change, with Cole's own examples: it is **asymmetric** (British adults: men are 8.4% taller than women, but women are 7.7% shorter than men — same two numbers, different magnitude) and **non-additive** (three successive 8% rises give 26%, not 24%). Logs fix both, because `log(x₂) − log(x₁) = log(x₂/x₁)`: differences on the log scale are symmetric and additive. Multiply by 100 and you get a **"sympercent" (s%)** — a symmetric percentage difference (equivalently, percent change with the mean of the two values as denominator). **This legitimizes the paper's log-ratio as the *correct interpretable relative-change scale*.** PDF mirror: gwern.net/doc/statistics/order/comparison/2000-cole.pdf. **Add to bib.**
- **Synthesis for the paper (resolves the draft's tension):** Cole justifies the log-ratio for **presentation**; Vickers/Kaiser justify baseline adjustment for **inference**. Combine them: run **ANCOVA on `log(PSA_followup)` with `log(PSA_baseline)` as covariate** (PSA is right-skewed / roughly log-normal, mCRPC values span single digits to thousands), then **report the group effect as a sympercent** (a symmetric percent change). The draft's mistake is using a log-ratio *t*-test for inference — that's still an unadjusted change score (Kaiser Model 1 vs 2; Vickers & Altman: change scores don't remove baseline imbalance). Right presentation scale, wrong inference engine.

### 2c. PSA / oncology-specific references (the applied layer)

These were surfaced during fetching but fell outside the final top-25 verification budget, so treat their specific numbers as needing a confirming read (grab the PDFs). They are directly on-point:

- **"Improving power in PSA response analyses of metastatic castration-resistant prostate cancer"**, *BMC Cancer* 2022 (s12885-022-09227-7). The single most on-point applied reference: directly critiques dichotomizing continuous PSA change into PSA50 ("a 31% reduction is treated the same as 90% but completely differently from 29%") and argues for continuous/higher-power analysis. **Open access.** **Add to bib; read closely — this paper may partly pre-empt or strongly support your thesis.**
- **PCWG3 (Prostate Cancer Working Group 3), Scher et al.**, *JCO* 2016 (PMC4872347). Recommends presenting PSA response as **per-patient percent change waterfall plots** and absolute change over time — i.e., the continuous distribution, not a single responder proportion. Useful: the field's own guideline already leans continuous. **Open access. Add to bib.**
- **JCO 2013** (PMC3805930), Armstrong/Halabi-type work: PSA declines ≥30%/≥50% evaluated as surrogates for overall survival in mCRPC. Context for why responder cutpoints persist (surrogacy claims). **Add to bib.**

### 2d. Open questions the research did **not** close (flag as limitations or resolve)

1. Direct power-loss quantification for dichotomizing the PSA *outcome* (vs predictor) — not found; your own simulation can supply this.
2. Canonical log-scale references beyond Cole (Törnqvist/Vartia on log ratios; Senn on log-scale analysis) — not verified in this pass; optional to chase.
3. Regulatory/clinical acceptability of a continuous replacement in the Phase II single-arm settings where PSA50 dominates — a Discussion point, not a stats point.

---

## 3. Draft (`index.qmd`) revisions

Ordered by importance.

1. **Retitle + reframe** around the two-critique structure (§1). Current title overclaims; "statistically incorrect" invites a fight. Consider "PSA50 is a statistically inefficient and information-losing endpoint."
2. **Background para 2 ([index.qmd:24](index.qmd)) — demote the ceiling-effect argument.** It's the weakest link: percent change is bounded above at 100% but PSA can rise unboundedly, so it isn't the clean floor/ceiling story the text implies. Lead instead with regression to the mean + inefficiency (Vickers, Bland & Altman).
3. **Fix the "seven assumptions" citation ([index.qmd:24](index.qmd)).** Harrell's actual argument is narrower and stronger (RTM + inefficiency → ANCOVA). Cite the specific BBR "change" chapter claim, not "seven assumptions," which will draw scrutiny.
4. **Resolve the log-ratio vs ANCOVA tension ([index.qmd:135, 141](index.qmd)).** The draft proposes a *t*-test on log ratios then dismisses ANCOVA for "linearity." That's backwards. A log-ratio *t*-test is still a change-score analysis (no baseline adjustment, still RTM-prone). Recommend **log-scale ANCOVA** (`log(FUP) ~ log(baseline) + group`) as the inference method, and **report the group effect as a sympercent** (Cole 2000) so you keep the interpretable relative-change number the draft wants. They coincide when the log-baseline slope = 1; ANCOVA estimates that slope instead of assuming it. Note the draft's "linearity" objection to ANCOVA actually argues *for* the log scale, not against ANCOVA — on the log scale the PSA-vs-baseline relationship is closer to linear. Back with Kaiser 1989 (Model 1 vs 2), Vickers & Altman 2001, Cole 2000.
5. **Keep the intuitive example ([index.qmd:121](index.qmd))** (both arms ~100% PSA50, Fisher n.s., but arm B far better) — it's a good *dichotomization/ceiling* illustration. But relabel it as demonstrating information loss from dichotomization, not "percent change is wrong."
6. **Abstract, figure captions, Results/Discussion** are stubs; write after the framing and simulation are settled.
7. Run the **`avoid-ai-writing`** skill over new prose before finalizing (Travis's standing preference).

---

## 4. Simulation fixes (`simulations/`)

There are **two parallel, partly-incompatible frameworks.** First decision: **consolidate to one.**

- [`simstemp.qmd`](simulations/simstemp.qmd) — older approach (Fisher vs log-ratio vs raw; simple mean-shift alternative). **Broken against current code:** calls `generate_PSA_sim()` with 8 args (current signature takes 7) and calls `raw_diff_test()`/`raw_after_test()` which don't exist in [`sim_functions.R`](simulations/sim_functions.R). Also an append bug at [simstemp.qmd:136](simulations/simstemp.qmd) (`raw_sim` where `after_sim` intended).
- [`sims.qmd`](simulations/sims.qmd) — newer approach (responder-curve mixtures; Fisher/%change/raw/ANCOVA/AUC). Incomplete (Months 6/9/12 empty, no results tabulation) and buggy.

**Recommendation:** keep `sims.qmd`'s richer generative model, delete or archive `simstemp.qmd`, and port the log-ratio and log-scale-ANCOVA tests into the `sims.qmd` framework.

### Confirmed bugs (fix these)

| # | Location | Bug | Effect |
|---|---|---|---|
| 1 | [sims.qmd:327](simulations/sims.qmd) | Medium-baseline **alternative** passes `props[[response[1]]]` for *both* A and B (should be `response[2]` for B) | Med-alt is secretly a null → medium-baseline power is wrong |
| 2 | [sims.qmd:191,233,285,326,368](simulations/sims.qmd) | Every loop uses `low_corrs[corr_num]` regardless of baseline; `med_corrs`/`high_corrs` ([sims.qmd:82–83](simulations/sims.qmd)) are never used | Baseline↔change correlation is not what's claimed for med/high |
| 3 | [sim_functions.R:391](simulations/sim_functions.R) | `halflife_calc` filters `id == 1` (should be `id == i`) | Every "average" half-life is just patient 1 |
| 4 | [sim_functions.R:264–282](simulations/sim_functions.R) | `log_percent_test` reads columns `psa_baseline`/`psa_after`/`patient_id` that don't exist in `df_psa_wide` (it has `baseline`/`FUP`/`id`); and `sims.qmd` never calls it | **The paper's own proposed method isn't actually simulated.** Wire it in and fix columns. |
| 5 | [simstemp.qmd:136](simulations/simstemp.qmd) | Appends `raw_sim` instead of `after_sim` | Copy-paste error (moot if file is archived) |

### Design issues (not bugs, but decide deliberately)

- **Two different estimands hide in the two frameworks.** `simstemp` alternative = a true mean-decline shift (tests raw power to detect a real PSA difference). `sims.qmd` alternative = different responder *mixtures* that can yield similar PSA50 rates but different kinetics (tests whether continuous metrics recover information PSA50 misses). **These support different papers.** Pick the one that matches the §1 thesis. The mixture design is the stronger match for "PSA50 loses information."
- **"Power" is not comparable across tests targeting different estimands** (Fisher on PSA50 vs *t*-test on AUC). Frame the comparison as *sensitivity to a clinically real difference*, not head-to-head power at a shared null. State the estimand explicitly.
- **Add log-scale ANCOVA** (`aov`/`lm` of `log(FUP) ~ log(baseline) + group`) as the recommended comparator, alongside the existing raw-scale `ancova_test` ([sim_functions.R:315](simulations/sim_functions.R)).
- Thin variance structure: within a response class at a given month, all patients get an identical deterministic curve value; the only noise is the `psa_FUP` term ([sim_functions.R:173](simulations/sim_functions.R)). Consider per-patient curve jitter so the null isn't artificially clean.
- Reproducibility: set seeds per chunk and record `renv`/session info before the final run.

---

## 5. references.bib additions

Add: Vickers & Altman 2001 (BMJ), Royston/Altman/Sauerbrei 2006 (Stat Med), Altman & Royston 2006 (BMJ), Cole 2000 (Stat Med), the BMC Cancer 2022 PSA-power paper, PCWG3/Scher 2016 (JCO), and the JCO 2013 PSA-surrogate paper. Verify the existing `harrell_2025_biostatistics` cite points at the BBR "change" chapter.

---

## 6. Work items (track from this doc — do NOT file as GitHub Issues)

Travis opted to work from this plan rather than file Issues. Treat the list below as the working checklist.

1. **Reframe thesis: separate dichotomization from percent-change critiques** (§1). Blocks prose rewrite.
2. **Fix medium-baseline alternative bug (`sims.qmd:327`)** — group B uses A's response mixture.
3. **Fix correlation-SD selection for med/high baseline** — `med_corrs`/`high_corrs` unused.
4. **`halflife_calc` uses `id == 1`** instead of `id == i` (`sim_functions.R:391`).
5. **Wire the log-ratio test into `sims.qmd` and fix its column names** (`log_percent_test`).
6. **Consolidate simulation frameworks** — archive `simstemp.qmd`; port needed tests into `sims.qmd`.
7. **Add log-scale ANCOVA as the recommended comparator.**
8. **Complete Results tabulation + Months 6/9/12**, or scope the paper to one timepoint.

---

## 7. Reference PDFs

**Obtained and read** (in gitignored `resources/`, local only — never commit/upload; copyright): Kaiser 1989, Royston/Altman/Sauerbrei 2006, Cole 2000. Their specifics are already folded into §2 above.

Everything else in §2 is open access. Still worth pulling and reading closely before writing the applied section: the **BMC Cancer 2022** PSA-power paper (§2c) — it's the closest prior work and may overlap the thesis.

---

## 8. Suggested kickoff order (fresh session)

1. Confirm the §1 thesis with Travis (blocks 3 & 4).
2. Fix the confirmed simulation bugs (§4) — mechanical, low-risk, unblocks trustworthy numbers.
3. Add log-scale ANCOVA + wire in log-ratio; decide on one framework/estimand.
4. Re-run simulations; build the Results table.
5. Rewrite Background + Methods prose against §2 citations; update references.bib.
6. Run `avoid-ai-writing` over new prose.
7. File Issues (§6) if using them for tracking.
