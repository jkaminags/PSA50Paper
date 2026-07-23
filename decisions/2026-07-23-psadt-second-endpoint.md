# Adding PSA doubling time as a second endpoint

**Date:** 2026-07-23
**Author of this note:** Travis (corresponding), for Josephine (first author) to review and iterate
**Target journal:** European Urology
**Status:** implemented on branch `two-critique-prose`, open to revision

This note explains why the paper now covers PSA doubling time alongside PSA50, what
that changed about the framing, and which parts are still yours to decide. It follows
`2026-07-20-two-critique-reframe.md`, which set up the two-critique structure this
builds on.

## The short version

The two-critique framing turned out not to be about PSA50. It is about how the field
summarizes a PSA trajectory, and PSA50 is one instance. PSA doubling time is another,
and it makes exactly the same two mistakes. Covering both makes the argument
substantially stronger, because a pattern that shows up twice in unrelated settings is
a property of the endpoints rather than a quirk of one simulation.

## Why doubling time belongs in this paper

PSA50 is a **response** measure. It reads how far PSA fell, in disease where treatment
is driving PSA down. PSA doubling time is a **progression** measure. It reads how fast
PSA is climbing, in biochemical recurrence and non-metastatic castration-resistant
disease. Different questions, different patients, opposite direction of travel.

That difference is worth stating plainly rather than smoothing over, because it is the
thing that makes the unification non-trivial. What a trial actually observes is a
patient's PSA over time, a curve that is close to a straight line on the log scale.
Every PSA endpoint is a summary of that one curve. Response endpoints read its drop.
Progression endpoints read its slope. Area under the curve reads all of it. Because
they are summaries of a common object, a design choice that damages one damages the
other, and doubling time is damaged the same two ways PSA50 is:

1. **A distorting scale.** PSADT is `log(2)` divided by the slope of log PSA on time.
   The slope is what the regression estimates and it is well behaved. The doubling time
   inverts it. Inverting a slope stretches the upper tail and fails outright when the
   slope is zero or negative, which is what happens whenever PSA is not rising. The
   PSA Working Group guideline that codifies the calculation
   ([Arlen 2008](https://doi.org/10.1016/j.juro.2008.01.099)) acknowledges this case
   and offers no remedy.

2. **A threshold.** Doubling time is habitually split into rapid and slow. The 10-month
   cutpoint governed eligibility for SPARTAN, PROSPER, and ARAMIS.

## What the simulation found

A new longitudinal engine (`simulations/lognormal_longitudinal_functions.R`, run by
`run_kinetics.R`) generates quarterly PSA over a year under a random-slope log-linear
model, then compares five analyses using the same size-adjusted power machinery the
response simulation uses.

To detect a slowing of median doubling time from 9 to 11 months, the rapid/slow flag
needs about 135 patients per arm against about 87 for a test on the log-scale growth
rate: a 1.56-fold difference.

The decomposition is the part worth your attention, because it reproduces the response
side almost exactly:

| | dichotomization costs | scale costs |
|---|---|---|
| Response (PSA50) | ~5 patients/arm | ~78 patients/arm |
| Progression (PSADT) | ~4 patients/arm | ~45 patients/arm |

Two simulations with no shared generative model, no shared clinical setting, and
opposite directions of PSA change produce the same asymmetry. The threshold is nearly
free in power terms; the scale is what costs the trial.

One finding is immediately actionable: simply taking the **logarithm of the doubling
time**, changing nothing else about the analysis, recovers about 40 of the 45 patients
per arm the reciprocal costs. A trial committed to reporting doubling times can get
most of the available efficiency without abandoning the quantity.

## Two things that were checked and are worth knowing

**The mixed model earns its place only under dropout.** Under a balanced visit
schedule the two-stage t-test on per-patient growth rates and the linear mixed model
tie exactly (0.798 vs 0.798 calibrated power), which is what theory predicts: with
equal visits every patient's slope has equal precision, so the simple average is
already optimal. Under dropout, where patients contribute two to five visits, they
separate (0.406 vs 0.585). The mixed model's advantage is correct weighting of
unequally informative patients, not a better scale. The paper says this rather than
implying the mixed model is always better, because a reviewer would find the balanced
case and the overclaim would cost more than the result is worth.

**The two scales fail for different reasons and the paper keeps them distinct.**
Percentage change is inefficient because it divides by a noisy baseline: that is
regression to the mean. Doubling time is inefficient because it inverts a slope: that
is a nonlinear transform with a pole at zero. Same pathology, different mechanism.
The RTM argument from the response side is *not* transplanted onto PSADT.

## Where you get to decide

- **The title.** It is now "PSA50 and PSA Doubling Time are Statistically Inefficient
  and Information-Losing Endpoints," which is the minimal extension of the previous
  one. A version that leads with the principle rather than the two endpoints is
  defensible and might read better for European Urology. Your call.
- **How much weight doubling time carries.** The current draft keeps PSA50 as the
  flagship and uses doubling time as the case that proves the problem belongs to the
  trajectory. It could be balanced more evenly if you think the progression setting is
  where the clinical appetite actually is.
- **The Vickers question.** Andrew Vickers has written skeptically about PSA velocity
  and doubling time, but only in the *pretreatment and screening* setting; he is
  explicit that the value of PSA kinetics in recurrent and advanced disease "has never
  been seriously questioned." Since this paper lives in the post-recurrence setting,
  citing him as a general critic would misrepresent him, and a European Urology
  reviewer who knows this literature would catch it. The draft cites him only for the
  setting-independent point that calculation methods proliferate and disagree. If you
  want to engage the pretreatment critique, it needs its own careful paragraph.
- **The novelty claim.** A literature search found no prior work framing the doubling
  time as a reciprocal reparametrization of the growth rate. That is good for the
  paper, but it means the argument rests on the derivation and simulation here rather
  than on established results. The Discussion says so. If you find prior art, that
  paragraph needs to change.
