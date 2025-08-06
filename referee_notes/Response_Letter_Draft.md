# Cover letter

## Major Points

**Major Point 1:** We agree with the referee's concern about unrealistic assumptions. We have **revised the paper to present results using realistic EB parameter diversity** sampled from the Kepler EB catalog, rather than fixed parameters. This approach provides **realistic performance expectations** for PLATO operations.

We emphasize that **our goal is not to provide realistic EB contamination rate estimates for PLATO** - such comprehensive analyses have been performed by other studies (e.g., Bray et al. 2023; Prša et al. 2011, 2016 for Kepler/TESS populations). Instead, we **focus on detection method performance** under controlled conditions. While all contaminants are assumed to be EBs (representing a conservative upper bound), **each EB now uses authentic transit depth and duration pairs** from observed Kepler eclipsing binaries, addressing the parameter realism concern.

The previous approach using fixed parameters (providing absolute worst-case bounds) is available in the first author's PhD thesis for reference.

**Major Point 2:** This concern is directly addressed by our revised approach using **realistic parameter diversity**. Rather than assuming uniform transit depths and durations, we now **sample (δ, td) pairs from the Kepler Eclipsing Binary Catalog** (Kirk et al. 2016), which provides authentic parameter distributions from observed systems.

This approach ensures that **each contaminant EB has realistic, varied parameters** rather than the artificial uniformity of fixed values. The parameter sampling preserves correlations between depth and duration found in real EB populations, providing **much more realistic performance assessments** for PLATO detection methods.

As noted in our response to Major Point 1, the uniform parameter approach (representing absolute worst-case bounds) remains available in the first author's PhD thesis for reference.

**Manuscript Changes:**

- Section 5.2: Describes realistic parameter sampling methodology
- Section 6: Updated calculations using variable parameters
- Section 7: Results show efficiency under realistic parameter diversity
- All efficiency plots and tables now reflect variable parameter scenario

**Major Point 3:** 
## Minor Points

**Minor Point 1**: We agree and the line in the manuscript saying that Kepler performed on-board photometry extraction was deleted in the introduction.

**Minor Point 2**: We agree with the referee. It has stated in the manuscript that L3 data products is the final planetary catalog from PLATO, and that catalog is not exclusively coming from the P1 sample.

**Minor Point 3:** We agree that this required clearer justification. We have added explicit explanation of why double-aperture photometry is cheaper than centroid measurements in Section 7:

_"Flux measurements require only a single scalar value per aperture, while centroid measurements must compute and transmit both x and y coordinates."_

This factor-of-two difference in data volume directly translates to the 50% reduction in CPU and telemetry costs. While specific resource allocation details are available in PLATO consortium technical documentation, the fundamental cost difference stems from this basic data volume requirement (1 value vs 2 values per measurement).

**Manuscript Changes:**

- Section 7: Added detailed justification for computational cost differences
- Clarified that cost reduction stems from 1 vs 2 data values per measurement

**Minor Point 4:** We have added a brief explanation in Section 2.4 clarifying that Equations (23) and (24) correct the placement of the (1 - SPR_tot) factor in the original Marchiori et al. formulations. **Detailed mathematical derivations are available upon request if required.**

**Minor Point 5:** We agree that the extended mask design choice required better justification. We have added a reference to previous analysis showing that **1-pixel ring extension provides a good performance**.

The choice of single-pixel extension is based on **systematic analysis of extended mask performance** (Gutiérrez-Canales 2025, Fig. 7.1), which demonstrates that **larger extensions (+1 pixels) lead to significant efficiency degradation** for extended flux measurements. This occurs because larger masks include more background noise while adding minimal contaminant signal.

**Manuscript Changes:**

- Section 3.1: Added justification for 1-pixel extended mask design with thesis reference
- Clarified that larger extensions decrease rather than improve efficiency

**Minor Point 6:** We have added clarification in Section 3.1 using inverse logic: **if a true planet were transiting the target star**, the extended light curve would show **equal or shallower depth** than the nominal curve due to dilution from additional flux in the larger aperture.

**The key insight is that deeper extended transits are physically incompatible with target-star planets**, making contamination the only viable explanation. This eliminates false positives while preserving all genuine planetary signals.

**Manuscript Changes:**

- Section 3.1: Added explanation of why deeper extended transits indicate contamination rather than alternative phenomena

**Minor Point 7:** SPRk values are computed theoretically from flux ratios and stellar positions, so there are no numerical precision errors or measurement uncertainties associated with individual SPRk calculations.

The probability that two contaminants in the same imagette would have identical SPRk values is negligible, given that SPRk depends on multiple parameters (stellar magnitude, position, PSF characteristics). In the extremely unlikely event of a tie, our algorithm would select the first contaminant encountered or, alternatively, could prioritize the brighter contaminant for mask centering.

This methodological detail does not affect our scientific conclusions, as ties have not been observed in our analysis of >7000 targets.

**Manuscript Changes:** None required.

**Minor Point 8:**  The factors of 10^-6 convert transit depths from parts per million (ppm) to dimensionless fractions required for the physical calculations. We have clarified this in Section 4.
**Manuscript Changes:**

- Section 4: Added clarification of ppm to dimensionless fraction conversion

**Minor Point 9:** Overlapping masks are technically possible in our methodology and do not create computational problems. However, mask overlap occurs primarily when target and contaminant stars are very close, making it difficult to disentangle their individual contributions.

This physical limitation is already captured in our efficiency measurements: secondary flux methods achieve ~90% efficiency, with the remaining ~10% likely representing cases where target-contaminant separation is insufficient for reliable secondary mask performance. When stars are too close, the secondary mask cannot effectively isolate the contaminant signal from the target signal.

Our methodology appropriately handles this by measuring **actual achieved efficiency** rather than assuming perfect performance, so these challenging configurations are already accounted for in our results.

**Manuscript Changes:** None required.

**Minor Point 10**: We agree that this required clarification. We have explicitly stated that Eq. (39) represents the uncertainty of averaged centroid measurements rather than single-point measurements. The clarifications emphasize that centroid shifts are averaged over one hour and multiple cameras to reduce statistical uncertainty.

**Manuscript Changes:**

- Before Eq. (39): Added explanation that centroid measurements are "averaged over a duration of one hour and a given number of cameras, NT, **in order to accordingly reduce the uncertainty**"
- Figure caption: Clarified that the plot shows "**Distribution of Eq. (39)**" (the averaged uncertainty) along target magnitude

**Minor Point 11:** The magnitude cutoff does not correspond to a specific PLATO sensitivity limit but rather reflects the completeness and reliability limits of our stellar catalog source.

**Manuscript Changes:**

- Section 5.1: Clarified that magnitude cutoff reflects Gaia DR3 catalog limits

**Minor Point 12:**  We agree that the efficiency computation conditions should be stated more explicitly. We have **added comprehensive specifications** for all detection thresholds and efficiency criteria:
**For flux measurements:**
 - Extended flux: η > 3.0 AND δ_extended > δ_nominal + 3σ_depth
 - Secondary flux: η > 3.0 AND δ_secondary > δ_nominal + 3σ_depth

**For centroid measurements:** 
- Nominal centroids: η > 3.0
- Extended centroids: η > 3.0
- Secondary centroids: η > 3.0

 Regarding the comment about "how quickly the threshold is considered to be passed" - Our statistical significances are theoretical computations (like SPR_k) with no numerical precision limitations, so we're uncertain how this applies to our threshold evaluations.
**Manuscript Changes:** 
- Section 6: Added explicit threshold criteria for all efficiency calculations

**Minor Point 13:** We have clarified that all P magnitude values refer to target star magnitudes, including for secondary mask cases. We have updated Section 7 text and figure captions to explicitly state "target star P magnitude." Additionally, we removed secondary masks from Figure B.2 since secondary mask sizes are determined by contaminant star properties rather than target star properties, which would create inconsistency in the magnitude reference framework. This ensures all remaining plots use target star magnitude consistently.

**Minor Point 14:** We agree that the objective of this analysis should be stated explicitly. We have clarified that this analysis verifies whether introducing extended and secondary masks stays within PLATO's on-board mask library constraints.

The critical concern is that adding new mask types could potentially exhaust the limited mask library budget (8,000 unique shapes total), which would make the approach operationally infeasible. We have added this explanation in Section 7:
"This is important because the introduction of extended and secondary masks could increase the total number of mask shapes such that the allowed, on-board limit of mask shapes is reached. In more detail, the on-board software has a mask library that consists of 8,000 different shapes that are shared among the three types of masks."
**Manuscript Changes:**
- Section 7: Added explicit explanation of why mask shape counting is operationally critical
