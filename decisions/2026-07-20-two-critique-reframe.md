# Reframing the paper around two separate critiques

**Date:** 2026-07-20
**Author of this note:** Travis (corresponding), for Josephine (first author) to review and iterate
**Target journal:** European Urology
**Status:** adopted, open to revision (see "Where you get to decide" at the end)

This note explains a framing decision for the paper and the simulation changes that
follow from it. It is written so you can see the reasoning, agree or push back, and
change direction if you land somewhere different. Nothing here is locked. If you read
the sources and disagree, that is a good outcome, not a problem.

## The short version

The current draft treats PSA50 as failing for one reason. It actually mixes two
different statistical problems that happen to point the same direction. They have
different mechanisms, different fixes, and different supporting literature. Pulling
them apart is the single change that does the most for the paper, because right now
the citations in the draft are about one problem while PSA50 is really about the other.

The two problems:

1. **Dichotomization.** PSA50 takes a continuous PSA trajectory and collapses it to
   a yes/no flag: did the patient hit a 50% decline or not. A 51% decline and a 95%
   decline become the same "responder"; a 49% decline and no response at all become
   the same "non-responder." Splitting a continuous measure into two bins throws away
   information and statistical power. This is what PSA50 *is*.

2. **Percent change as a continuous scale.** Analyzing `(baseline - follow-up) / baseline`
   as a number. This is a different issue. Even without any dichotomizing, percent
   change is a poor scale: it does not remove baseline imbalance (regression to the
   mean), it is a ratio so it is skewed and asymmetric, and it is the least efficient
   of the standard analysis choices.

Most of the methods literature in the draft (Vickers, Kaiser, Bland & Altman) is about
problem 2, the continuous percent-change metric. But PSA50 is problem 1, the
dichotomized flag. The draft slides between them as if they were one argument. A
sharp reviewer will notice that the cited evidence does not quite match the target.

## Why separate them instead of just picking one

Because the fixes are different, and the paper is stronger if it names both honestly.

- The fix for dichotomization is: don't dichotomize; analyze the continuous outcome.
- The fix for percent change is: don't analyze percent change on its raw scale. Adjust
  for baseline, and if you want a relative-change number, compute it on the log scale.

If the paper only argued "PSA50 dichotomizes," a reader could reasonably say "fine, then
just use raw percent change as a continuous variable," and walk straight into problem 2.
If it only argued "percent change is a bad scale," it would be citing the right papers
for the wrong target, because PSA50 is not raw percent change; it is a threshold on it.
Naming both closes the escape hatches and matches the literature to the claim.

**Proposed thesis:** PSA50 discards information twice. It dichotomizes, and the thing
it dichotomizes (percent change) is itself a poor scale. And there are efficient,
clinically interpretable alternatives to both problems.

**Proposed title direction:** something like "PSA50 is a statistically inefficient and
information-losing endpoint" rather than "statistically incorrect." "Incorrect" invites
a fight over a word; "inefficient and information-losing" is precise and defensible.

## The evidence, matched to each critique

These are the core references. The specifics below came from a verified literature pass;
the PDFs for the three starred ones are in the local `resources/` folder (gitignored;
do not commit them, they are copyrighted).

**For dichotomization (problem 1):**
- **Royston, Altman & Sauerbrei 2006**: dichotomizing a continuous variable at the
  median keeps about 65% of the information for a normal variable, and only about 48%
  for a right-skewed one. PSA is right-skewed, so the ~48% figure is the relevant one.
  Read it as "you are discarding roughly half your information."
- **Altman & Royston 2006**: the companion short piece. Median-split power loss is
  roughly like throwing away a third of your patients.

**For percent change as a scale (problem 2):**
- **Vickers 2001**: a simulation study. Percent change has the lowest power of the
  standard options, and its power collapses when baseline variability increases. His
  recommendation: don't analyze percent change directly; use ANCOVA.
- **Kaiser 1989**: the "when is percent change actually okay" reference. His answer:
  pick the scale (change vs percent change) that is independent of baseline, and check
  by plotting each against baseline. The wrong choice costs about 40 to 45% more
  patients for the same power.
- **Bland & Altman 2011**: testing before-vs-after within each arm separately (rather
  than between arms) inflates the false-positive rate badly. Relevant to how single-arm
  PSA50 studies get interpreted.

**For the fix (the honest counterpoint that percent change is fine on the log scale):**
- **Cole 2000** ("Sympercents"): conventional percent change has two concrete defects.
  It is asymmetric (a 10% rise and a 10% fall are not the same size), and it is
  non-additive (three successive 8% rises give 26%, not 24%). Logs fix both, because a
  difference of logs is a log of the ratio. Multiply a log difference by 100 and you get
  a "sympercent," a symmetric percentage difference. This is what legitimizes reporting
  our effect as a relative change: on the log scale, it behaves properly.

## The closest prior work: Grayling et al. 2022

There is one applied paper we have to position against carefully, because it is the
nearest neighbor to this thesis and it is good:

Grayling MJ, McMenamin M, Chandler R, Heer R, Wason JMS. "Improving power in PSA
response analyses of metastatic castration-resistant prostate cancer trials."
*BMC Cancer* 2022;22:111. Open access.

What they did: a literature review of 64 articles covering 78 mCRPC treatment arms with
waterfall plots, then applied the "augmented binary" method. The idea is to keep the
binary PSA50 responder endpoint but recover the continuous percent-change data that sits
underneath it (Box-Cox transformed to normality), which lets you estimate the responder
proportion far more precisely. Their headline number: a median efficiency gain of
**103.2% (IQR 89.8 to 190.9%)**, meaning the augmented analysis reaches the precision a
trial would otherwise need roughly double the sample size to achieve. A 100-patient trial
gets the precision of about 200 patients, with no extra data collection.

Why this matters for us, in three parts:

1. **It supports critique 1, and it does so in exactly our setting.** They have already
   quantified, for mCRPC PSA specifically, how much is lost by throwing away the
   continuous data behind the binary flag. So the paper should *not* present "dichotomizing
   PSA loses power" as a new finding. Cite Grayling as having shown it, and build on it.

2. **Their fix and ours are different, and that difference is our contribution.** They
   keep the binary PSA50 responder rate as the thing being estimated and make that
   estimate more efficient. This paper argues something they do not: that the binary
   estimand *itself* discards clinically meaningful PSA kinetics. Two patients with very
   different response curves (a fast deep responder and a slow shallow one) can share the
   same PSA50 status, and the survival literature says those curves are not equivalent.
   Augmentation cannot recover that, because it is still aimed at the responder proportion.
   The kinetics argument, and metrics like AUC that capture it, is where this paper adds
   something Grayling does not.

3. **They leave problem 2 untouched.** Their method does not address regression to the
   mean, baseline adjustment, or ANCOVA at all. So our log-scale ANCOVA recommendation is
   complementary to their work, not a competitor to it. Also worth noting: they Box-Cox
   transform toward normality rather than committing to the log scale; we argue for the
   log scale specifically, on interpretability grounds (Cole's sympercents).

Net position: Grayling does not pre-empt the paper. It is the strongest piece of support
for the dichotomization half, and the honest way to use it is to say "the efficiency cost
of the binary flag is established (Grayling); we go further and show the binary flag also
discards response-kinetics information that no amount of augmentation recovers, and we
give interpretable, baseline-adjusted alternatives." Their analysis code and a comparison
Shiny app are public, which matters for the simulation decision below.

## The method we're recommending, and why the draft's version is off

The draft proposes a *t*-test on log ratios, then dismisses ANCOVA for "linearity."
This is backwards. The distinction is subtle, but it decides which test you should run,
so it is spelled out here.

A *t*-test on log ratios is still a **change score**: it takes each patient's
`log(follow-up) - log(baseline)` and compares the two groups. A change score does not
adjust for baseline. It assumes the right way to remove baseline is to subtract it
(here, on the log scale). That assumption is exactly what regression to the mean
violates. So the log-ratio *t*-test inherits problem 2's baseline issue even though it
uses the better (log) scale.

The fix that keeps the good scale and removes the baseline problem is **ANCOVA on the
log scale**: `log(follow-up) ~ log(baseline) + group`. It adjusts for baseline instead
of assuming subtraction is correct, and it estimates the baseline relationship from the
data rather than fixing it at 1 (which is what a change score silently assumes). Then
you report the group effect back-transformed as a sympercent, so you still get the
interpretable relative-change number the draft wanted.

The draft's own "linearity" objection to ANCOVA actually argues *for* the log scale,
not against ANCOVA. On the log scale the PSA-versus-baseline relationship is closer to
linear. So keep ANCOVA and move it to the log scale, rather than dropping it.

Sources for this: Kaiser 1989 (his Model 1 vs Model 2), Vickers & Altman 2001 (ANCOVA
is the preferred general approach), Cole 2000 (the sympercent for reporting).

## Two things to be careful about (so a reviewer can't catch us)

1. **The efficiency numbers (48%, 65%) are about dichotomizing a *predictor*, not an
   *outcome*.** PSA50 dichotomizes an outcome. The mechanism (information loss) is the
   same, but the exact percentages transfer by analogy, not by direct measurement.
   Present them as illustrative, and let our own simulation (and Grayling's applied
   figure) supply the outcome-specific numbers. Do not state the 48% as if it were
   measured for PSA50.

2. **Do not use the "optimal cutpoint" bias argument.** Royston also shows that choosing
   a cutpoint from the data inflates false positives to 25 to 50%. That does *not* apply
   to PSA50, because 50% is a fixed, agreed-upon threshold, not one we searched for in the
   data. Using that argument here would be wrong, and it is exactly the kind of thing a
   reviewer will flag. The argument that *does* apply is information and power loss only.

## What changed in the simulation code

The simulation was updated so it matches this framing and produces trustworthy numbers:

- **Fixed three bugs** that were corrupting results: the medium-baseline alternative
  scenario was using group A's response mixture for group B (so it was secretly a null);
  every loop was using the low-baseline noise level regardless of baseline (the medium
  and high noise levels were defined but never used); and the half-life function was
  averaging patient 1 sixty times instead of averaging all sixty patients.
- **Wired in the log-ratio test and added a log-scale ANCOVA test.** The log-ratio test
  existed but read column names that don't exist and was never called; it now runs. The
  log-scale ANCOVA (`log(FUP) ~ log(baseline) + group`) is the method this note
  recommends, so the simulation now actually compares it against the others.
- **Floored simulated PSA at 0.1 ng/mL.** The generator could produce negative PSA
  values, which are nonphysical and break the log-scale tests. 0.1 ng/mL is the
  functional sensitivity of a standard clinical PSA assay, which is the assay class used
  in mCRPC trials (ultrasensitive assays that reach 0.01 ng/mL and below are a
  post-prostatectomy biochemical-recurrence tool, not a metastatic-disease one). Because
  mCRPC PSA values run high, this floor is essentially non-binding; it only removes the
  rare nonphysical negatives. This is settled, not an open question.
- **Archived `simstemp.qmd`** to `archive/` (older framework, broken against the current
  functions). The project now consolidates on `sims.qmd`, which has the richer generative
  model (response-curve mixtures). That mixture design is the better match for the
  dichotomization argument, because it produces scenarios where two groups can have
  similar PSA50 rates but different kinetics, which is precisely the information PSA50
  throws away.

## Where you get to decide

These are open, and your call as first author:

- **The framing itself.** If you read Kaiser and Cole and come away thinking the
  one-critique version is cleaner, say so. The two-critique version is a recommendation,
  not a verdict.
- **Whether to add Grayling's augmented binary method as a comparator.** It is the
  strongest competing approach, and including it in the simulation (their code and a
  Shiny app are public) would make the head-to-head much more convincing to a European
  Urology reviewer. It is also a real amount of work: their method needs a Box-Cox
  transform and delta-method variance, not a one-line test like the others. My lean is to
  include it, but it is a scope decision worth making deliberately.
- **Estimand and scope.** The mixture design tests "can continuous metrics recover
  information PSA50 misses." An alternative design tests raw power to detect a real mean
  difference. These support slightly different papers, so we should be explicit about
  which one we are making. Related: do we report all of months 3/6/9/12, or scope to one
  timepoint cleanly? Right now only month 3 is built out.

For context on the last big open question: the applied prostate-cancer literature is
worth leaning on, and Grayling is the anchor for it (see the section above). European
Urology is a clinical journal, so the paper should lead with the clinical message (what
PSA50 hides from a trialist and a patient) and keep the statistical machinery accessible,
with the kinetics/information argument doing the heavy lifting that sets us apart from the
pure-efficiency framing of Grayling.

If any of this doesn't sit right, that is the point of writing it down, so we can argue
about it with the reasoning in front of us.
