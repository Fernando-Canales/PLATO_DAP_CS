# Cover letter

## Major Points

**Major Point 1:** "The paper assumes that all background stars are eclipsing binaries undergoing transit-like events, and therefore potentially capable of generating a false positive signal. This is an unrealistic assumption, and I think an unnecessary one. I appreciate that, as stated in section 5.2, it is not possible to know the details of individual contamination sources. However, as also stated in section 5.2, there are studies about the number and occurrence rates of eclipsing binaries. Could these have been used to assign eclipsing binary status to suitable fractions of the stars in the sample, allowing for a more realistic proportion of possible eclipsing binary contaminants in each window? Such an approach should be considered."

**Response:** We emphasize that our goal is not to provide realistic EB contamination rate estimates for PLATO - such comprehensive analyses have been performed by other studies (e.g., Bray et al. 2023; Prša et al. 2011, 2016 for Kepler/TESS populations). Instead, we focus on detection method performance under controlled conditions. We have added at the end of the second paragraph that "Furthermore, even if it would be useful, it is not in the scope of this paper to realistically determine expected EB populations for PLATO." While all contaminants are assumed to be EBs (representing a conservative upper bound), each EB now uses authentic transit depth and duration pairs from observed Kepler eclipsing binaries, addressing the parameter realism concern.  We have added in Sect. 5.2 that "we don't aim to estimate realistic EB rates for PLATO, as we mentioned in the introduction."

**Manuscript changes**
- Introduction: Mentioned the realistic EB rate is out of the scope of this paper
- Section 5.2: Mentioned the papers that actually explored EB rates for PLATO and other space missions

**Major Point 2:** "The eclipse depth of the eclipsing binary appears in equation (10) and propagates forwards into the equations for statistical significance for a transit on a contaminant, and thus into the equations for efficiency of the extended mask and secondary mask flux methods. In section 2.4 the paper notes that the intrinsic eclipse depth is assumed to be the same for all contaminant eclipsing binaries. This seems like a very limiting assumption, and while I understand that it makes the study simpler to perform, I wonder whether it is limiting the study and strongly constraining the applicability of the results. Would the results change if a range of eclipse depths is present among the contaminating stars, or is only the mean eclipse depth relevant? This assumption that all eclipsing binaries have the same eclipse depth should be explained in more detail, including information about the justification for it, whether any other scenarios were considered, and how it limits the results that are presented. "

**Response:** We agree with the referee's concern about unrealistic assumptions. We have revised the paper to present results using realistic EB parameter diversity sampled from the Kepler EB catalog, rather than fixed parameters. This approach provides realistic performance expectations for PLATO operations. This concern is directly addressed by our revised approach using realistic parameter diversity. Rather than assuming uniform transit depths and durations, we now sample (δ, td) pairs from the Kepler Eclipsing Binary Catalog (Kirk et al. 2016), which provides authentic parameter distributions from observed systems.

This approach ensures that each contaminant EB has realistic pair of parameters.

**Manuscript Changes:**

- Section 5.2: Described the sampled methodology
- Section 6: Updated calculations using variable parameters
- Section 7: Results show efficiency under realistic parameter assumptions
- All efficiency plots and tables are changed to show the variable transit parameter scenario

**Major Point 3:** "The paper misrepresents the process by which PLATO will validate its exoplanet candidates. On page 2, in the left-hand column, starting at line 50, it refers to PLATO relying heavily on an initial validation process performed onboard. There is no validation performed on-board the spacecraft; the centroid data and/or double-aperture photometry data are collected onboard, but the validation is performed on-ground as part of the exoplanet data processing pipeline. The paper then states that the initial PLATO strategy for detecting false positives consisted only of computing centroid shift measurements. This is incorrect; other methods for detecting false positives were considered right from the start. It is more correct to say that centroid shift measurements were part of the initial false positive detection strategy, but the only additional data from the spacecraft that were considered (beyond the basic photometry). Similarly, section 6.2 states that nominal centroids are the strategy envisaged for 5% of the P5 sample, which is incorrect; other methods will be employed, with the centroids providing an additional discriminatory ability for those 5% of the targets. Furthermore, the paper also states that that current strategy to detect false positives is to use centroid shifts and/or double-aperture photometry. This is stated in the introduction (page 2, right-hand column) and in Section 4 (page 6, left-hand column). While they will be important, these methods will form only part of the false positive identification strategy for the mission, which will also use a range of other techniques. The description of the role that centroids and double-aperture photometry play in candidate validation should be updated to more accurately represent both the history here, and what will be done"

**Response:** We have made substantial revisions to correct these misconceptions:
- Data collection vs validation distinction: We have clarified that centroid and double-aperture data are collected on-board for the majority of PLATO targets, but that validation analysis occurs on-ground for all type of light-curves. However, for targets for which the flux is computed on board, validation based on PLATO is at this stage only possible whenever centroid is also computed on board

- Historical accuracy: We have corrected the claim that centroids were the initial strategy, acknowledging that multiple false positive detection methods were considered from the beginning.
- 5% allocation clarification: We cited the PLATO SciRD as the reference where it can be found that  5% of the P5 targets will have on-board nominal centroid measurements. We also mentioned that the focus of our work are the targets that will have their photometry obtained on-board, that is in fact the majority of PLATO targets.

**Manuscript Changes:**

- Page 2: Corrected on-board vs on-ground processing description
- Introduction: Re positioned centroid/double-aperture methods as part of broader validation strategy
- Section 4: Clarified role as contributing methods rather than primary strategy

**Major Point 4:** "Aside from the aforementioned assumptions regarding the eclipsing binary contaminants, the paper makes a number of other assumptions: a. The L1 pipeline has perfectly removed any systematics and background flux. This seems reasonable, given that the true performance of the pipeline in this respect cannot be evaluated until real data is available. b. There is no residual drift of stars across the CCD. This seems problematic to me; is ‘residual drift’ the drift remaining after correction? Will the drift vary throughout the field of view, and will the effectiveness of the correction be equally good across that field of view? c. When defining the equations for the centroid shifts, that there is no stellar variability. d. When assessing instrument performance, that each camera has the same set of PSFs. Won’t the optics, and thus the PSFs, of each camera vary slightly? The level of this variation will not be known until flight data is available, but on-ground testing data may provide some useful insights already and could be used to define or constrain some small random differences between simulated PSFs for each camera. These varied assumptions are not justified in the paper. Such justifications should be added."

**Response:** We agree that the mentioned assumptions required explicit justification. We have added such justifications as follow:

a) L1 assumption: Nothing to add

b) Residual Stellar Drift: We have clarified in Section 2.3 that "residual drift" refers to remaining positional drift after L1 correction**, acknowledging that drift and correction effectiveness vary across the FoV and that detailed assessment is ongoing consortium work.

c) Stellar Variability: We have added in Section 5.2 that "Another important assumption of this work is related to stellar variability, that is not explicitly modeled in our flux and centroid shift calculations. This simplification is justified because intrinsic stellar variations typically occur on timescales different from the transit detection averaging process described in Sections 2.3 and 4.2 "

d) PSF Uniformity: We have justified this in the following way

- Section 5.3: PSF set reflects realistic camera-to-camera variations from Monte Carlo simulations 
- Section 8.4: We refer to Gutierrez-Canales (2025) to show that the efficiency of the metrics does not change a lot under a series of different PSF configurations (temperature, diffusion kernel value, etc.)

**Manuscript Changes:**

- Section 2.3: Added residual drift clarification and ongoing work
- Section 5.2: Added explanation about stellar variability not relevant for our computations
- Section 5.3: Added PSF variation justification
- Section 7: Added real camera PSF validation reference

**Major Point 5:** "The results presented in Section 7 are not sufficiently described or explained. a. In Figure 4b, why is the secondary mask equally efficient for the 6 camera and 24 camera cases at P magnitude 11.0? This is a notable outlier from the behaviour evident at the other magnitudes considered. b. I’m surprised that in Figure 4a, the 24 camera and 6 camera cases are so similar, and that for the secondary mask centroid the 6 camera case is sometimes more efficient. Why is this the case? c. Is the y-axis of Figure 5 the noise of the extended and nominal mask metrics, or the uncertainty in the metric? Per equation (39), this is the error in the absolute centroid shift, averaged over 1 hour and Nt cameras, which could be taken to be either. This should be clarified."

**Response:** We appreciate these detailed observations. So far, in the paper (see response to Major Point 2) the adopted scenario is that every contaminant is an EB but now the transit depth and duration of each EB is randomly sampled from observed distributions. Therefore, the unusual behaviors noted by the referee in a) and b) (equal efficiencies at P=11.0, 6-camera outperforming 24-camera) no longer appear in our revised results using realistic parameter diversity. These anomalies were artifacts of the unrealistic fixed parameter assumption. The revised results show consistent expected behavior: 24 cameras outperform 6 cameras across all methods and magnitudes.

Furthermore, significant changes were performed to the efficiency computations. The main changes are the following:
1) There was a bug in the python code regarding the expressions for the significant transit depth conditions in Eq.(50) and Eq.(57). The bug was fixed and the overall efficiency of both extended and secondary fluxes increased by ~2%. 
2) We introduced a third condition for the efficiency computation of extended, nominal and secondary centroid shifts (namely, Eq.(62), Eq.(65) and Eq.(68)). This third condition has to do directly with Eq.(36) , that we derived based on the hypothesis that centroid error should be significantly smaller than centroid shifts. However, we noticed that Eq.(36) didn't hold for some of the measurements. In particular, we noticed that we included in our results centroid shifts computed with secondary masks of just 1 pixel in size.  Naturally, centroid shifts computed with 1-pixel secondary masks are obviously zero and therefore Eq.(36) no longer holds. Therefore, the condition that we added to Eq.(62), Eq.(65) and Eq.(68) is (Delta C > 10 sigma/sqrt(ntr td) where we used the corresponding values of Delta C and sigma for each mask) makes sure that centroid shifts are always bigger than the corresponding centroid errors.  When implementing this new condition to all the centroid metrics, the overall efficiency of the extended, nominal and secondary centroid shifts changed significantly. Now the extended centroids have an efficiency of ~87%, nominal centroids of ~83% and secondary centroids of ~75%. Due to this decrease in nominal, extended and secondary centroids efficiencies, secondary flux measurements are now the most efficient metric overall.

c) Figure Y-axis Clarification: We have clarified that the y-axis shows "uncertainty in the centroid shift" (not raw noise), representing σ_centroid from Eq. (39) - the error in absolute centroid shift averaged over 1 hour and N_T cameras.

Enhanced Results Description: We have added comprehensive explanations for each method's performance characteristics and physical interpretations for observed efficiency patterns.

**Manuscript Changes:**

- Figures: Updated to Scenario B results, eliminating previous anomalies
- Section 7: Added detailed explanations of efficiency patterns and physical reasoning
- Figure captions: Clarified y-axis as "uncertainty" rather than "noise"
- Added diagnostic analysis: Explaining magnitude-dependent camera sensitivity (P=10.5, P=11.5)

**Major Point 6:** "It is not clear to me how the magnitude of the target stars is incorporated into the calculations performed, as it does not seem to be mentioned in the derivations worked though in sections 2 to 6. Table 2 shows the efficiencies determined across a magnitude range. For these values, I guess that the efficiency is determined for all target stars in this magnitude range, and the mean/median and standard deviation determined to give the values in the table. Is this correct? Whether it is or not, the steps by which the values in Table 2 were determined should be described. In Figure 4, the paper displays efficiencies at specific magnitudes. How is the magnitude used here? Presumably not all of the target stars are exactly the magnitudes plotted, so were they rounded? Or are these magnitude bins (as suggested by footnote 7), and if so, how are the bins defined? The steps by which the values plotted in Figure 4 were determined should be described."

**Response:** We agree that the magnitude incorporation methodology required explicit description. We have added explanations of both the efficiency plot generation and table value derivation:

Efficiency Plot Methodology: We have described the complete binning approach:

- 7 magnitude bins of 0.5 magnitudes each (P = 10.0 to 13.0)
- 1000 randomly selected targets per bin from P ± 0.25 range
- Each data point represents efficiency for all 1000 targets in that bin

Table Value Derivation: We have explicitly linked table and figures, stating:  "Table 2 shows the averaged efficiency of the metrics presented in this work. The values presented are the magnitude-averaged efficiencies across the entire P5 range of the values in Figures 5a and 5b"

Target vs Contaminant Magnitude: We have explicitly stated that all magnitude references correspond to target star magnitudes, including for secondary mask analyses, maintaining methodological consistency.

**Manuscript Changes:**

- Before efficiency plots: Added a description of the binning methodology
- Before efficiency table: Clarified that values are magnitude-averaged from plots
- Section 7: Added target star magnitude clarification for all methods
- Figure captions: Updated to specify "target star P magnitude"

**Major Point 7:** "The discussion of the results in Section 8.1 is not sufficiently detailed and overstates the conclusions that can be drawn from Table 2. The paper states (section 8.1, page 13) that the secondary flux measurements are less efficient at detecting false positives than centroid measurements. However, the results in Table 2 show that this is not strictly correct. Secondary flux measurements are less efficient than nominal and extended centroid shift measurements, true, but their efficiency is equal to that of secondary centroid shift measurements, within uncertainties. One could conclude that they are even marginally more efficient than secondary centroid shift measurements. Similarly, the first sentence of section 8.2 states that centroid measurements are by far the most efficient method for detecting false positives. But this is only true of nominal mask centroid shifts and extended mask centroid shifts, which are equally efficient (within uncertainties) according to Table 2. Secondary mask centroid shifts, on the other hand, are no better than secondary flux measurements, as I noted above. The paper also states (Section 8.1, final paragraph, on page 13) that on average, the statistical significance of the secondary flux measurement is greater than the statistical significance of the centroid shift in the nominal mask (this is stated as an equality). No quantitative evidence is given to support this conclusion, as Table 2 presents efficiencies, and Figure 4 plots efficiencies. Support for the conclusion can be taken from Table D.1, but this presents the various significance values for only one star and its contaminants. This is not enough evidence to support the statement about the relationship of the average significances. Further discussion of the results is needed, with more quantitative support for the conclusions reached, and careful wording that more accurately describes the results."
**Response:**

**Response**: We agree that in Section 8.1 the descriptions were overstated from Table 2 and misleading. Furthermore, due to the changes in the code and the conditions for centroid metrics, the results referred to in this Major Point are changed. However, Section 8.1 now clearly states the efficiencies of the new results and discuss them in detail, alongside with the values for the new comparative table/matrices between metrics that are mentioned in the next Major Point. We recall that the paper has new efficiency results for the metrics and that these new results are explained in the new Sections 7 and 8. More details are given in the response to Major Point 5.

We also agree that on the manuscript we mention that, on average,  the statistical significane of the secondary flux is greater than the statistical significance of the nominal centroid shift but no quantitative evidence was given. We have added to Section 8.1 what is now Fig.7, a figure that contains two histograms, the one on the left refers to the counts of log_10(eta_sec/eta_nom_cob) and the one on the right refers to log_10(eta_sec/eta_ext_cob).  Both histograms support our claim that secondary flux significances are on average greater than the one of nominal centroids and also extended centroids.



**Major Point 8:** "Presentation and discussion of the percentage of false positives detected by each metric, and the comparisons between them, is incomplete. There are five metrics being compared in this paper, giving 10 possible comparisons. But the set of bullet points in the right-hand column of page 13 gives only four comparisons. The secondary flux measurements, which the paper is proposing are the most effective metric, are not mentioned at all, which seems a significant omission. How do they compare to the other metrics in this context? Perhaps the paper is using “extended” in this section in the on-board software sense, as described in footnote 5? If so, this is unnecessarily confusing given that the rest of the paper has used “extended” to mean specifically one of the two double aperture approaches, and it obfuscates the results. Furthermore, the language of these bullet points could be improved as it is somewhat confusing as currently written. The current format is "% of false positives detected only by metric a but not by metric b". If the false positives are detectable 'only' by metric a then mentioning metric b is unnecessary. On the other hand, if the percentages are for false positives detected by metric a but not metric b, do the other metrics detect them or not? Section 8.2 should be updated to better present the comparison of false positive detection percentages between metrics. A table or pairwise matrix may be a better way to present these results, rather than a list of bullet points, and further discussion may be warranted. In addition, Footnote 5 should be removed, as it introduces unnecessary ambiguity, and the text updated (if necessary) to clarify the language used to discuss the results. Finally, the summary of the comparison results in the text should be checked for consistency against the full set of comparison results (currently it reads ~27% vs 28%+/-1.4%). "

**Response:** We completely agree with this comprehensive review and have  restructured Section 8.2 to address all concerns:

1. Complete Comparison Coverage: We have replaced the incomplete 4-point bullet list with two pairwise matrix tables (Tables 3 and 4) showing the possible comparisons between the metrics. We have included the abbreviations EFX for extended flux, ECOB for extended centroid shifts, NCOB for nominal centroid shifts, SFX for secondary flux and SCOB for secondary centroid shifts. Table 3 shows the comparison between extended mask and nominal mask-based methods since for those metrics we consider the first 10 contaminants in terms  of SPRk in each window, while for secondary mask-based methods we only consider one contaminant per window, the one with the highest SPRk value. This is related to the next point.

2. Secondary Flux Inclusion: Table 4 now explicitly includes secondary flux measurements (SFX), addressing the significant omission noted by the referee. Since secondary-mask based method focus only on one contaminant in each window, we cannot directly compare them with extended mask or nominal mask-based methods. Therefore, we created Table 4, including only comparison between SFX and SCOB.
3. Language Clarity and precision: We have clarified the language - Tables 3 and 4 show "percentage of FPs detectable only by the method in the row but not by the method in the column," instead of the confusing "only by this method... but not by this other method" . Also, all the numbers in the Tables have the same, consistent number of decimal figures.

We mention as well that all the numbers and figures regarding this Major Point are now different in the new version of the paper. This is because we changed the efficiency conditions for the centroid shift efficiencies, as mentioned in more detail on the response of Major Point 5. 

**Manuscript Changes:**

- Replaced: 4-point bullet list with comprehensive matrix tables
- Added: All missing method comparisons, especially secondary flux
- Enhanced: Discussion based on complete pairwise analysis
- Clarified: Percentage interpretation language throughout
- Removed: Confusing footnote 5 and associated ambiguity

**Major Point 9:** "Appendix A discusses the detectability of planets, but as with sections 8.1 and 8.2, the discussion is not sufficiently detailed and is imprecise in its language. Three classes of planets are considered, but these classes are not defined. What radii (or range of radii) are used for this analysis? ‘super-Earth’, in particular, has a number of different definitions in the literature. The last sentence of the section states that the results show that using 6 cameras is sufficient to detect super-Earths, but no evidence is provided to support this statement. Figure A1 shows the case of super-Earths observed with 24 cameras, not 6 cameras. "

**Response:** We thank the referee for catching this inconsistency, which revealed **a caption error**. Figure A1 **does indeed show Super-Earths with 6 cameras** (orange stars), providing direct evidence for our conclusion about 6-camera sufficiency. The figure caption incorrectly stated "24 cameras" - this has been corrected.
The evidence supporting our claim is the following: Super-Earth significance values (orange stars) are well above the η_min detection threshold when using only 6 cameras, demonstrating adequate detectability.

We have also added explicit planet definitions following Borucki et al.(1996) classifications:

- Earth-like planets: δ_p = 84 ppm, t_d = 13 hr
- Super-Earths: δ_p = 522 ppm, t_d = 42 hr (interpolated between Earth and Neptune)
- Jovian planets: δ_p = 10,100 ppm, t_d = 29.6 hr

**Manuscript Changes:**

- Corrected figure caption: Super-Earths now correctly stated as "6 cameras"
- Added explicit planet parameter definitions with literature references
- Clarified that significance above η_min threshold demonstrates detectability

## Minor Points

**Minor Point 1**: "The introduction references Kepler performing on-board photometry. Kepler did not do this; it downloaded pixel data for defined postage stamps around its targets. This statement should be removed, or updated to better match / explain the intended meaning." 

**Response:** We agree and the line in the manuscript saying that Kepler performed on-board photometry extraction was deleted in the introduction.

**Minor Point 2**: "In section 2.3, on page 3, in the last paragraph of the right-hand column, the paper references PLATO’s Level 3 data products consisting of the catalogue of confirmed planetary systems from P1. This is not quite correct, as the catalogue is not confined a priori to the P1 sample. This statement should be corrected."

**Response:** We agree with the referee. It has stated in the manuscript that L3 data products is the final planetary catalog from PLATO, and that catalog is not exclusively coming from the P1 sample.

**Minor Point 3:** "On page 2 the paper references double aperture photometry requiring 50% less CPU and telemetry than centroid shifts. Similar statements about the reduced CPU and telemetry consumption of double aperture photometry are present on page 13. The justification for these statements is not provided, nor is a citation given to support these statements. Such a justification or supporting citation should be added."

**Response:** We agree that this required clearer justification. We have added explicit explanation of why double-aperture photometry is cheaper than centroid measurements in Section 7:

_"Flux measurements require only a single scalar value per aperture, while centroid measurements must compute and transmit both x and y coordinates."_

This factor-of-two difference in data volume directly translates to the 50% reduction in CPU and telemetry costs. While specific resource allocation details are available in PLATO consortium technical documentation, the fundamental cost difference stems from this basic data volume requirement (1 value vs 2 values per measurement).

**Manuscript Changes:**

- Section 7: Added detailed justification for computational cost differences

**Minor Point 4:** "In section 2.4, equations (23) and (24) are referred to as corrected versions of equations (28) and (23) in Marchiori et al. (2019). What has been corrected for this paper, and why did it need to be corrected? These statements should be clarified."

**Response:** We have added a brief explanation in Section 2.4 clarifying that Equations (23) and (24) correct the placement of the (1 - SPR_tot) factor in the original Marchiori et al.(2019) formulations. Detailed mathematical derivations are available upon request if required.

**Minor Point 5:** "In section 3.1, the paper defines the optimal extended mask as being an expansion of the nominal mask. The description of the size of this expansion is unclear, as it refers to surrounding the nominal mask by one pixel, but Figure 1b nicely illustrates the idea. This expansion by one pixel in all directions is described as the optimal definition, with reference made to minimizing the increase in noise, but surely the truly optimal way in that case would be add just a small number of pixels covering the contaminants, rather than a full ring as shown in Figure 1b? This would have the benefit of not including pixels with minimal flux, this limiting the additional background contribution. What I think is missing here is an explanation of the process for determining the optimal size of the extended mask. Section 3.1 mentions a process of increasing the size of the extended mask by one more pixel and checking the impact on efficiency, with a note that this has been done before. Where has this been done, and by whom? Is this a reference to section 8.1 on page 13, where efficiency values are shown for the case of the extended mask considering just the most significant contaminant? Please support this statement with further explanation, a cross-reference to section 8.1, or a suitable citation. "

**Response:** We agree that the extended mask design choice required better justification. We have added a reference to previous analysis showing that 1-pixel ring extension provides a good performance. Furthermore, the choice of single-pixel extension is based on systematic analysis of extended mask performance. For instance, Fig 7.1 in Gutiérrez-Canales 2025 (which is available here: https://theses.hal.science/tel-05165095) shows that even larger extensions of the nominal mask (+2 pixels) lead to significant efficiency degradation for extended flux measurements. This occurs because larger masks include more background noise while adding minimal contaminant signal.  However, we acknowledge that more optimal ways to build extended masks can be implemented. One of these is, as suggested by the referee, to add a small number of pixels covering contaminants. Similar ways to described optimal extended masks are on-going work within the consortium. This is mentioned at the end of Sect. 3.1 and also at Section 8.4 


**Manuscript Changes:**

- Section 3.1: Added justification for 1-pixel extended mask design with thesis reference

**Minor Point 6:** "In section 3.1, on page 5, right-hand column, the paper states that if there is a deeper transit in the light curve from the extended mask than in the light curve from the nominal mask, then this means that the signal in the target is due to one of the contaminants. Can this be concluded definitively? Or in this situation is it just very probably that the signal is due to one of the contaminants, but there could still be other explanations?"

**Response:** We have added clarification in Section 3.1 using inverse logic: if a true planet were transiting the target star, the extended light curve would show equal or shallower depth than the nominal curve due to dilution from additional flux in the larger aperture.

**The key insight is that deeper extended transits are physically incompatible with target-star planets**, making contamination the only viable explanation. This eliminates false positives while preserving all genuine planetary signals.

**Manuscript Changes:**

- Section 3.1: Added explanation of why deeper extended transits indicate contamination rather than alternative phenomena

**Minor Point 7:** "When calculating SPRk to identify the most significant contaminant for defining the secondary mask, to what precision is this calculated? If two contaminants in the imagettes have equal SPRk, which is used to define the secondary mask? "

**Response:** SPRk values are computed theoretically from flux ratios and stellar positions, so there are no numerical precision errors or measurement uncertainties associated with individual SPRk calculations.

The probability that two contaminants in the same imagette would have identical SPRk values is negligible, given that SPRk depends on multiple parameters (stellar magnitude, position, PSF characteristics). In the extremely unlikely event of a tie, our algorithm would select the first contaminant encountered or, alternatively, could prioritize the brighter contaminant for mask centering.

This methodological detail does not affect our scientific conclusions, as ties have not been observed in our analysis of >7000 targets.

**Manuscript Changes:** None required.

**Minor Point 8:** "In Equation (32), where do the factors of 10-6 come from?" 

**Response:** The factors of 10^-6 convert transit depths from parts per million (ppm) to dimensionless fractions required for the physical calculations. We have clarified this in Section 4.

**Manuscript Changes:**

- Section 4: Added clarification of ppm to dimensionless fraction conversion

**Minor Point 9:** "Can the secondary mask pixels overlap with the nominal mask pixels? This point is not clear from the paper, and the example in Figure 1 is just one possible configuration. For example, if the most significant contaminant in the Figure 1 field was the star located at (2.5, 1.5), would the secondary mask method still work, and would it still be an improvement over the extended mask? In such cases of overlapping nominal and secondary masks, would differential analysis have merit?"

**Response:** Overlapping masks are technically possible in our methodology and do not create computational problems. However, mask overlap occurs primarily when target and contaminant stars are very close, making it difficult to disentangle their individual contributions. Furthermore, we changed the overall look of Fig. 1, now the displayed secondary mask is not only one pixel in size.

This physical limitation is already captured in our efficiency measurements: secondary flux methods achieve ~90% efficiency, with the remaining ~10% likely representing cases where target-contaminant separation is insufficient for reliable secondary mask performance. When stars are too close, the secondary mask cannot effectively isolate the contaminant signal from the target signal.

Our methodology appropriately handles this by measuring **actual achieved efficiency** rather than assuming perfect performance, so these challenging configurations are already accounted for in our results.

**Manuscript Changes:** None required.

**Minor Point 10**: "Equation (39) is described as the average centroid measurement over a duration of one hour and a given number of cameras. But the equation itself seems to be for the uncertainty in the averaged centroid. This should be checked and either the text or the equation updated if necessary."

**Response:** We agree that this required clarification. We have explicitly stated that Eq. (39) represents the uncertainty of averaged centroid measurements rather than single-point measurements. The clarifications emphasize that centroid shifts are averaged over one hour and multiple cameras to reduce statistical uncertainty. We mention as well that we average over one hour because transits typically last only a few hours and also because the values of the transit duraition from the catalog we are using are in hours.

RS->  remind why this choice of time-scale (this is linked to the transit duration, which typical last few hours)

**Manuscript Changes:**

- Before Eq. (39): Added explanation that centroid measurements are "averaged over a duration of one hour and a given number of cameras, NT, **in order to accordingly reduce the uncertainty**"
- Figure caption: Clarified that the plot shows "**Distribution of Eq. (39)**" (the averaged uncertainty) along target magnitude

**Minor Point 11:** "Why is magnitude 21 the faint cutoff for the stellar sample used in this work? Does this correspond to the sensitivity limit of PLATO? Is it the magnitude at which stars stop being relevant as potential contaminants? "

**Response:** The magnitude cutoff does not correspond to a specific PLATO sensitivity limit but rather reflects the completeness and reliability limits of our stellar catalog source. We have added: in Section 5.1 "This magnitude range and cutoff for faint targets comes from Gaia DR3,  rather than PLATO sensitivity constraints."

**Manuscript Changes:**

- Section 5.1: Clarified that magnitude cutoff reflects Gaia DR3 catalog limits

**Minor Point 12:** "In section 6.1.1, the three conditions used in the computation of NextFP are not sufficiently well-defined in the text. At the start of the left-hand column of page 9. First, what is the threshold used for determining that the transit depth difference is statistically significant? Per equation (45) it seems to be 3-sigma, but this should be stated explicitly. Secondly, at what precision are the statistical significances and thresholds (eta terms) measured? This affects how quickly the threshold is considered to be passed. These conditions should be expanded on in the text. Similar comments apply to the conditions in sections 6.1.2, 6.2.1 and 6.2.2."

**Response:** We agree that the efficiency computation conditions should be stated more explicitly. We have specified the criteria for all detection thresholds and efficiency criteria:
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

**Minor Point 13:** "In section 7 and Appendix B, the P magnitude values are referenced as belonging to the target star. This is clear for the nominal mask and extended mask cases, but for the secondary mask case, is it still the magnitude of the target that is relevant, or the magnitude of the most significant contaminant (i.e. the one used to position the secondary mask)? This should be clarified in the text, and also in the appropriate figure captions and/or axis labels."

**Response:** We have clarified that all P magnitude values refer to target star magnitudes, including for secondary mask cases. We have updated Section 7 text and figure captions to explicitly state "target star P magnitude." Additionally, we removed secondary masks from Figure B.2 since secondary mask sizes are determined by contaminant star properties rather than target star properties, which would create inconsistency in the magnitude reference framework. This ensures all remaining plots use target star magnitude consistently.

**Minor Point 14:** "The relevance of the discussion of the number of unique mask shapes is not clear to me and should be stated in the paper. Is the point to show that the increase in possible mask configurations introduced by the secondary mask does not reach the limit allowed by the size of the on-board library?" 

**Response:** We agree that the objective of this analysis should be stated explicitly. We have clarified that this analysis verifies whether introducing extended and secondary masks stays within PLATO's on-board mask library constraints.

The critical concern is that adding new mask types could potentially exhaust the limited mask library budget (8,000 unique shapes total). We have added this explanation in Section 7:
"This is important because the introduction of extended and secondary masks could increase the total number of mask shapes such that the allowed, on-board limit of mask shapes is reached. In more detail, the on-board software has a mask library that consists of 8,000 different shapes that are shared among the three types of masks."

**Manuscript Changes:**
- Section 7: Added explicit explanation of why mask shape counting is operationally critical

## Editorials
We have attended all the editorial points from the referee's report.
