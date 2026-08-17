### Supplementary Figure 5 | Negative control: pure isolates fed through the decomposition

**a**, Group confusion matrix for the uncorrected pipeline. Archaea acts as a
false-positive sink (Plantae-to-Archaea 23%, Bacteria-to-Archaea 17%). **b**,
Group confusion matrix for the corrected pipeline adopted in Fig. 5 (internal-standard
normalisation, response-efficiency correction with a 0.20 floor and an
uncalibrated fallback of 1.0, the ArchLips-restricted archaeal reference, and
distribution of enriched-feature weight across the phyla sharing a feature). The
Archaea sink is removed (non-archaeal-to-Archaea at most 0.8%) and archaeal
self-recovery rises from 68% to 93%; the dominant group is correct for 79% of
isolates (130/164), and for every group except Plantae and Animalia the
correction raises self-recovery (Bacteria 63 to 79%, Fungi 48 to 69%,
Protozoa 51 to 59%). **c**, Per-group self-recovery, uncorrected versus corrected.

All 164 isolates were decomposed leave-one-out, with the phylum-centroid
reference rebuilt after removing the test isolate, so no sample is matched
against itself. Both panels use the same 736-feature substrate and the same
isolates, so the difference between them is the correction stack alone. Analysis
units are the 16 phyla of the locked release; the six groups shown are display
summaries, as in Fig. 5.

Panels **a**-**c** use the fold-change-weighted similarity estimator, the
living-community estimator that produces the bars in Fig. 5. Under the
marker-panel estimator shown as diamonds in Fig. 5, overall dominant-group
accuracy is higher (83.5%, 137/164) and the plant and animal isolates that this
estimator recovers poorly are recovered well (Plantae 47 to 95%, Animalia 47 to
82%), while archaeal recovery falls to 58%. The two estimators fail on different
groups, consistent with their different estimands: the similarity family
estimates living community composition, the marker panel the provenance of
matched lipid signal.

Persistent limitation, unchanged by the correction: plant and animal isolates
are the least well recovered under the similarity estimator, mirroring the
reference-panel bias flagged for Animalia in Fig. 5.

Values are a documented reimplementation. The submitted producer
(`analysis-19/16_negative_control/build_figure.py`) is not in the repository, so
the previously published panel values are not reproducible and the two sets of
numbers are not interchangeable.
