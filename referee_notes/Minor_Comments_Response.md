#### First minor point
>	I. The introduction references Kepler performing on-board photometry. Kepler did not do this; it downloaded pixel data for defined postage stamps around its targets. This statement should be removed, or updated to better match / explain the intended meaning. 

**Response**: We agree and the line has to be removed since Kepler did not use on-board photometry

#### Second minor point
>	II. In section 2.3, on page 3, in the last paragraph of the right-hand column, the paper references PLATO’s Level 3 data products consisting of the catalogue of confirmed planetary systems from P1. This is not quite correct, as the catalogue is not confined a priori to the P1 sample. This statement should be corrected.

**Response**: We agree. A quick look at the Science Requirement Document or at the last PLATO paper could fix this one.

#### Third minor point
> 	 III. On page 2 the paper references double aperture photometry requiring 50% less CPU and telemetry than centroid shifts. Similar statements about the reduced CPU and telemetry consumption of double aperture photometry are present on page 13. The justification for these statements is not provided, nor is a citation given to support these statements. Such a justification or supporting citation should be added.

**Response**: We agree. I have to provide a specific sentence saying that for centroid shift measurements we need two quantities (x and y centroid position) and that for double-aperture photometry we only need one quantity (flux measurements). We can't however refer to any ATBD or internal technical note because such documents are not available to the public. 

#### Fourth minor point
> 	IV. In section 2.4, equations (23) and (24) are referred to as corrected versions of equations (28) and (23) in Marchiori et al. (2019). What has been corrected for this paper, and why did it need to be corrected? These statements should be clarified. 

**Response**: We agree. There is no specific mention on what was changed. I could mention in a couple of sentences the terms that disappeared from Marchiori equations in order to get our equations. We have PDF files on the correction of such equations. So, if the sentences are not clear, worst-case-scenario we can provide the full correction as an Appendix.

#### Fifth minor point
> 	  V. In section 3.1, the paper defines the optimal extended mask as being an expansion of the nominal mask. The description of the size of this expansion is unclear, as it refers to surrounding the nominal mask by one pixel, but Figure 1b nicely illustrates the idea. This expansion by one pixel in all directions is described as the optimal definition, with reference made to minimising the increase in noise, but surely the truly optimal way in that case would be add just a small number of pixels covering the contaminants, rather than a full ring as shown in Figure 1b? This would have the benefit of not including pixels with minimal flux, this limiting the additional background contribution. What I think is missing here is an explanation of the process for determining the optimal size of the extended mask. Section 3.1 mentions a process of increasing the size of the extended mask by one more pixel and checking the impact on efficiency, with a note that this has been done before. Where has this been done, and by whom? Is this a reference to section 8.1 on page 13, where efficiency values are shown for the case of the extended mask considering just the most significant contaminant? Please support this statement with further explanation, a cross-reference to section 8.1, or a suitable citation. 

**Response**: We agree. This is a long paragraph, but the point of it is that I should mention that I already did an analysis for the efficiency of extended masks built with a 2-pixel ring around the corresponding nominal mask. I could even mention my PhD thesis as a reference. I can also include some of the discussion from my thesis to this comment about optimal extended masks.

#### Sixth minor point
>	VI. In section 3.1, on page 5, right-hand column, the paper states that if there is a deeper transit in the light curve from the extended mask than in the light curve from the nominal mask, then this means that the signal in the target is due to one of the contaminants. Can this be concluded definitively? Or in this situation is it just very probably that the signal is due to one of the contaminants, but there could still be other explanations?

**Response**: We agree. The question is interesting. 'Could there be other explanations?' Réza suggested to answer it in an 'negative/inverse way'. This means to write a sentence saying that if there is a true planet transiting the target star and we have a nominal mask and an extended mask in that imagette, then the extended light curve will be, at most, equal in depth to the nominal light curve. As a mater of fact, in such a case the extended light curve will be less deep because of the dilution. This shows that extended light curves couldn't trigger a false negative and also that is difficult for extended light curves that are deeper than nominal light curves to be originated by a transit in the target. 

#### Seventh minor point
> 	VII. When calculating SPRk to identify the most significant contaminant for defining the secondary mask, to what precision is this calculated? If two contaminants in the imagettes have equal SPRk, which is used to define the secondary mask? 

**Repsonse**: The SPRk computation is completely 'theoretical' so there are no numerical precision or error measurements associated to it. Also, the chance that two contaminants in a given imagette have the same SPRk value is very unlikely. It could happen that two contaminants have very similar SPRk values, but in that case we could take the first one or the brighter one. 

#### Eight minor point
> 	VIII. In Equation (32), where do the factors of 10-6 come from?

**Response**: I have to stay clearly in the text that the factor comes from the use of ppm in the transit depth.

#### Ninth minor point
> 	IX. Can the secondary mask pixels overlap with the nominal mask pixels? This point is not clear from the paper, and the example in Figure 1 is just one possible configuration. For example, if the most significant contaminant in the Figure 1 field was the star located at (2.5, 1.5), would the secondary mask method still work, and would it still be an improvement over the extended mask? In such cases of overlapping nominal and secondary masks, would differential analysis have merit?

**Response**: We agree. Interesting question. In principle there is no problem when it comes to overlapping masks in an imagette. However, if a target and a contaminant are very close then the efficiency of secondary masks are not that good since it becomes difficult to disentangle the contribution from the target and the contaminant. Currently, secondary fluxes have an efficiency of 90% so the remaining 10% could be related to cases where the target and the contaminant are very close to each other.

#### Tenth minor point
> 	X. Equation (39) is described as the average centroid measurement over a duration of one hour and a given number of cameras. But the equation itself seems to be for the uncertainty in the averaged centroid. This should be checked and either the text or the equation updated if necessary. 

**Response**: We agree. I have to say clearly that we are dealing with a centroid shift that is averaged over one hour and a given number of cameras to reduce the uncertainty accordingly. Therefore the quantity in Eq. (39) is the error or uncertainty associated to the averaged centroid.

#### Eleventh minor point
> 	XI. Why is magnitude 21 the faint cutoff for the stellar sample used in this work? Does this correspond to the sensitivity limit of PLATO? Is it the magnitude at which stars stop being relevant as potential contaminants? 

**Response**: The cutoff does not correspond to any PLATO sensitivity limit (not that I know of). The cutoff comes from GAIA DR3.

#### Twelfth minor point
> 	XII. In section 6.1.1, the three conditions used in the computation of NextFP are not sufficiently well-defined in the text. At the start of the left-hand column of page 9. First, what is the threshold used for determining that the transit depth difference is statistically significant? Per equation (45) it seems to be 3-sigma, but this should be stated explicitly. Secondly, at what precision are the statistical significances and thresholds (eta terms) measured? This affects how quickly the threshold is considered to be passed. These conditions should be expanded on in the text. Similar comments apply to the conditions in sections 6.1.2, 6.2.1 and 6.2.2.

**Response**: We agree. I have to stay more explicitly all the conditions for the efficiency computation of the metrics (extended and secondary fluxes as well as nominal, extended and secondary centroids). The statistical significances are also theoretical computations, like SPRk, so no 'precision' is involved. However, the remark 'this affects how quickly the threshold is considered to be passed' is a bit unclear and I don't really understand what does that mean in this context.

#### Thirteenth minor point
> 	XIII. In section 7 and Appendix B, the P magnitude values are referenced as belonging to the target star. This is clear for the nominal mask and extended mask cases, but for the secondary mask case, is it still the magnitude of the target that is relevant, or the magnitude of the most significant contaminant (i.e. the one used to position the secondary mask)? This should be clarified in the text, and also in the appropriate figure captions and/or axis labels. 

**Response**: Appendix B is more detailed in my PhD thesis. I will consider what I wrote there in order to be more clear and answer this comment.

#### Fourteenth minor point
> 	XIV. The relevance of the discussion of the number of unique mask shapes is not clear to me and should be stated in the paper. Is the point to show that the increase in possible mask configurations introduced by the secondary mask does not reach the limit allowed by the size of the on-board library? 

**Response**: We agree. I should state clearly in the paper that one of the objectives is to check that the number of unique mask shapes for the nominal, extended and secondary masks is within the limits of the current mask library budget. 