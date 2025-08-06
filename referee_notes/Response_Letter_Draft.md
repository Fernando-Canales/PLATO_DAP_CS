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

**Minor Point 6:**

**Minor Point 13:** We have clarified that all P magnitude values refer to target star magnitudes, including for secondary mask cases. We have updated Section 7 text and figure captions to explicitly state "target star P magnitude." **Additionally, we removed secondary masks from Figure B.2** since secondary mask sizes are determined by contaminant star properties rather than target star properties, which would create inconsistency in the magnitude reference framework. This ensures all remaining plots use target star magnitude consistently.